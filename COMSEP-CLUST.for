c=====================================================================
c This is the source code of the method "COMSEP-Clust" introduced in  
c "A Novel Optimization Approach Towards Improving Compactness 
c                                        and Separability of Clusters"
c Code is written by Adil Bagirov (Novemeber 2021) and has two phases: 
c Phase 1. Compactness of clustering
c Phase 2. Separability of clustering
c The name of input file should be: "datainput"
c The results will be written in output files:
c Clust_member: Cluster memberships
c Center: Center of clusters
c Results: Objective function value, CPU time, DB index, Purity, 
c number of misclassified points & value of total misclassification cost 
c=====================================================================
      PARAMETER(maxdim=4000, maxrec=435000, maxclust=15, maxnft=3,
     1 maxclass=1, maxsize=1000000)
     
      implicit double precision(a-h,o-z)
      double precision x(maxdim),a(maxrec,maxnft),dminim(maxrec)
     1 ,x2(maxsize),amed(maxnft),z(maxdim),x3(maxdim),xbest(maxdim)
     2 ,fval(maxrec),x5(maxsize),x6(maxdim),fval1(maxrec),xb(maxdim)
     3 ,plabel(maxclass),dcent(maxclust,maxclust)
      integer nk(maxclust,maxrec),nel(maxclust),nob(maxclass)
     1 ,lcand(maxrec)
      common /csize/m,/c22/a,/anclust/nclust,/cnft/nft,/cmf/mf,/ceps/eps
     1 ,/crecord/nrecord,/cnc/nc,/cnk/nk,nel,/adminim/dminim,/cnf/nf
     2 ,/cmultiple/multiple,/ccand/ncand,lcand,/ctnorm/tnorm,/cns/ns
     3 ,/cnob/plabel,nob,/cnclass/nclass,/cnpnout/npurity,/cdcent/dcent
     4 ,/cgamma/gamma1,gamma2,/ctoler/toler,/cnversion/nversion
     5 ,/ctau/tau1,tau2,/cnphase/nphase
      open(39,file='Clust_member.txt')      
      open(40,file='Clusters.txt')
      open(41,file='Center.txt')
      open(42,file='Results.txt')
      open(78,file='datainput.txt',status='old',
     1 form='formatted')
c========================================================================
      PRINT *,' '
      PRINT *,'Enter the number of features:'
      read *,nft        
      PRINT *,' '
      PRINT *,'Class outputs?'
      PRINT *,' '
      PRINT *,'      1 - if yes'
      PRINT *,'      2 - if no'
      read *,npurity
      PRINT *,' '
      if (npurity.eq.1) then
         PRINT *,'Number of classes:'
         read *,nclass
         PRINT *,' '
         PRINT *,'Output column:'
         read *,noutcom
         PRINT *,' '
      end if
      PRINT *,'Maximum number of clusters:'
      read *,nclust       
      PRINT *,' '
      PRINT *,' '
      PRINT *,'Version?'
      PRINT *,' '
      PRINT *,'      1 - only clustering'
      PRINT *,'      2 - new model'
      read *,nversion
      PRINT *,' '
c===================================================================
      tau1=9.0d-01
      tau2=1.0d+01
c===================================================================
      tlimit=7.2d+04
      do i=1,maxrec
       read(78,*,END=901,ERR=900) (a(i,k),k=1,nft)
       nrecord=i
      end do
  900 stop 'Error in input file'       
  901 WRITE(40,*) 'Input complete. Number of records: ',nrecord
      WRITE(40,*)
c===================================================================
      if(nrecord.le.200) then
       gamma1=1.0d+00
       gamma2=1.0d+00
       gamma3=1.0d+00
      end if

      if((nrecord.gt.200).and.(nrecord.le.1000)) then
       gamma1=8.0d-01
       gamma2=5.0d-01
       gamma3=4.0d-01
      end if
      
      if((nrecord.gt.1000).and.(nrecord.le.5000)) then
       gamma1=7.0d-01
       gamma2=4.0d-01
       gamma3=3.0d-01
      end if

      if((nrecord.gt.5000).and.(nrecord.le.15000)) then
       gamma1=2.0d-01
       gamma2=1.0d-01
       gamma3=1.0d-01
      end if
      if((nrecord.gt.15000).and.(nrecord.le.50000)) then
       gamma1=4.0d-01
       gamma2=2.0d-01
       gamma3=1.0d-01
      end if
      if(nrecord.gt.50000) then
       gamma1=3.0d-01
       gamma2=1.0d-01
       gamma3=1.0d-01
      end if
c===================================================================
      write(40,311) gamma1,gamma2,gamma3
  311 format('Par-ter:',' g1=',f7.4,' g2=',f7.4,' g3=',f7.4)
      write(40,*)
c==================================================================
      IF(npurity.eq.1) mf=nft-1
      IF(npurity.eq.2) mf=nft
      IF(npurity.eq.2) noutcom=0
      tnorm=0.0d+00
      call cpu_time(time1)

      if ((npurity.eq.1).AND.(noutcom.lt.nft)) then
       do i=1,nrecord
        j1=0
        do j=1,nft
         if (j.ne.noutcom) then
          j1=j1+1
          amed(j1)=a(i,j)
         end if
        end do
        amed(nft)=a(i,noutcom)
        do j=1,nft
         a(i,j)=amed(j)
        end do
       END do
      end if
      call cpu_time(time1)
c===========================================================
      IF(npurity.eq.1) THEN
       do i=1,nclass
        nob(i)=0
       end do
       a2=a(1,noutcom)
       plabel(1)=a2
       nob(1)=1
       n2=1
       do i=2,nrecord
        clabel=a(i,noutcom)
        do j=1,n2
         IF(clabel.eq.plabel(j)) GO TO 22
        end do
        n2=n2+1
        plabel(n2)=clabel
        j=n2
 22     nob(j)=nob(j)+1
       end do
       do i=1,nclass
         WRITE(40,41) i,plabel(i),nob(i)
       end do
      END if
 41   FORMAT('Class',i4,' label',f7.0,'  number of instances',i7)
c================================================================
       WRITE(42,*)
       WRITE(42,701)
 701   FORMAT('#Clust','                  fval ','        penalty',
     1 '   negatives','    purity','      DB')
       WRITE(42,*) 
c================================================================
      do nc=1,nclust  
       PRINT 42,nc
 42    FORMAT('Cluster No.:',i10)
       if(nc.eq.1) then
                   call step1(f,x)
                   toler=1.0d-02*f/dble(nrecord)
                   go to 1
       END if
       nphase=1
       call step2(nstart,x2)
       fbarmin=1.0d+28
       fbarmax=0.0d+00
       do j=1,nstart
        do k=1,mf
         z(k)=x2(k+(j-1)*mf)
        end do
        ns=1
        m=mf
        call optim(z,x6,barf)
        fval(j)=barf
        fbarmin=dmin1(fbarmin,barf)
        fbarmax=dmax1(fbarmax,barf)
        do k=1,mf
         x2(k+(j-1)*mf)=x6(k)
        end do
       end do
c============================================================
       fbarmin=fbarmin+gamma3*(fbarmax-fbarmin)
       nstart1=0
       do j=1,nstart
        if (fval(j).le.fbarmin) then
         nstart1=nstart1+1
         do k=1,mf
          x5(k+(nstart1-1)*mf)=x2(k+(j-1)*mf)
         end do
         fval(nstart1)=fval(j)
        end if
       end do

       nstart=nstart1
       do i=1,nstart
        do k=1,mf
         x2(k+(i-1)*mf)=x5(k+(i-1)*mf)
        end do
       end do
c===========================================================
       do k=1,mf
        x5(k)=x2(k)
       end do
       fval1(1)=fval(1)
       nstart2=1
       do j=2,nstart
        do j1=1,nstart2
         f31=0.0d+00
         do k=1,mf
          f31=f31+(x5(k+(j1-1)*mf)-x2(k+(j-1)*mf))**2
         end do
         IF(f31.LE.toler) THEN
          if(fval1(j1).ge.fval(j)) then
           fval1(j1)=fval(j)
           do k=1,mf
            x5(k+(j1-1)*mf)=x2(k+(j-1)*mf)
           end do
          end if
          GO TO 1200
         END IF 
        end do
        nstart2=nstart2+1
        do k=1,mf
         x5(k+(nstart2-1)*mf)=x2(k+(j-1)*mf)
        end do
        fval1(nstart2)=fval(j)
1200   end do
       do i=1,nstart2
        do k=1,mf
         x2(k+(i-1)*mf)=x5(k+(i-1)*mf)
        end do
       end do
       nstart=nstart2
c==============================================
       m=mf*nc
       fbest=1.0d+28
       do j=1,nstart
        do i=1,mf
         x(i+(nc-1)*mf)=x2(i+(j-1)*mf)
        END do
        do j1=1,m
         x3(j1)=x(j1)
        end do

        ns=2
        call optim(x3,x6,fcurrent)
        if (fcurrent.lt.fbest) then
         fbest=fcurrent
         do j1=1,m
          xbest(j1)=x6(j1)
          xb(j1)=x6(j1)
         end do
        end if
       end do

       call clusters2(xbest,f)
       call silhouette(nneg0,ratio,pen)
       if(ratio.gt.9.95d-01) go to 3
       
       if((nversion.eq.2).and.(nc.gt.1)) then
        nphase=2
        ns=2
        tau2=1.0d+00
        do it=1,3 
         do i=1,m
          x(i)=xb(i)
         END do
         call optim(x,x6,fcurrent)
         call clusters2(x6,f)
         dif=f-fbest-7.5d-02*fbest
         if(dif.gt.0.0d+00) go to 3
         call silhouette(nneg,ratio,pen)
         if(ratio.gt.9.99d-01) then
          do i=1,m
           xbest(i)=x6(i)
          end do
          go to 3
         end if
         if(nneg.lt.nneg0) then
          nneg0=nneg
          do i=1,m
           xbest(i)=x6(i)
          end do
          if(it.eq.1) tau2=4.0d+00
          if(it.eq.2) tau2=1.6d+01  
          if(it.eq.3) tau2=6.4d+01  
         end if        
         if(nneg.ge.nneg0) then
          if(it.eq.1) tau2=4.0d+00
          if(it.eq.2) tau2=1.6d+01       
          if(it.eq.3) tau2=6.4d+01
         end if        
        end do 
       end if
       
  3    do j1=1,m
        x(j1)=xbest(j1)
       end do
       
       call clusters2(x,f)
       call silhouette(nneg,ratio,pen)

       do i=1,nc
        dcent(i,i)=0.0d+00
       end do
       
       dist1=0.0d+00
       do i=1,nc
        dist2=0.0d+00
        do j=i+1,nc
         dcent(i,j)=0.0d+00
         do j1=1,mf
          dcent(i,j)=dcent(i,j)+(x(j1+(i-1)*mf)-x(j1+(j-1)*mf))**2
         end do
         dist2=dmax1(dist2,1.0d+00/dcent(i,j))
         dcent(j,i)=dcent(i,j)
        end do
        dist1=dist1+dist2
       end do
c================================================================
  1    call clusters(x,f)
c========================================
c  Printing cluster centers
c========================================
       if(nc.gt.1) then
        write(41,*)
        write(41,411) nc
        write(41,*)
        do i=1,nc
         write(41,412) i,(x(j+(i-1)*mf),j=1,mf)
        end do
       end if
411   format('The number of clusters:',i10)
412   format('Center No.:',i3,' = (',59f15.5,')')   
c========================================
       if(nc.gt.1) call finresult(x,f,purity,db)
       call cpu_time(time3)
       time4=time3-time1
       write(40,603) nc,f,tnorm,time4
 603   format(i6,f28.4,f18.0,f11.3)
c===================================================================
      if(nc.gt.1) then
       write(42,499) nc,f,pen,nneg,purity,db
      end if  
 499  format(i6,f28.4,f10.4,i10,f12.3,f10.4)     
       IF(time4.gt.tlimit) GO TO 2
      end do
   2  continue
      CLOSE(39)
      close(40)
      CLOSE(41)
      CLOSE(42)
      close(78)
      stop
      end

c====================================================================
c The first phase: Compactness step starts here
c====================================================================
      subroutine clusters(x,f)
      PARAMETER(maxdim=4000, maxrec=435000, maxclust=15, maxnft=3,
     1 maxclass=1)
      implicit double precision (a-h,o-z)
      double precision x(maxdim),a(maxrec,maxnft),dminim(maxrec)
     1 ,rad(maxclust),radmax(maxclust),dcent(maxclust,maxclust)
     2 ,ratio(maxclust)
      integer nk(maxclust,maxrec),nel(maxclust),lcand(maxrec)
     1 ,list1(maxrec),lcand1(maxrec)
      common /c22/a,/cmf/mf,/anclust/nclust,/crecord/nrecord,/cnc/nc
     1 ,/cnk/nk,nel,/cnclass/nclass,/ccand/ncand,lcand,/adminim/dminim
     3 ,/ctnorm/tnorm,/clist1/list1,/cdcent/dcent
     
      do j=1,nclust
       nel(j)=0
       rad(j)=0.0d+00
       radmax(j)=0.0d+00
      END do

      f=0.0d+00
      do k=1,nrecord
       f2=1.0d+28
       do j=1,nc
        if(j.gt.1) then
         dif1=dcent(j,jmin)-4.0d+00*f2
         if(dif1.ge.0.0d+00) go to 7
        end if  

        f1=0.0d+00
        do k1=1,mf
         f1=f1+(a(k,k1)-x(k1+(j-1)*mf))**2
        END do
        tnorm=tnorm+1.0d+00
        if(f2.gt.f1) then
         f2=f1
         jmin=j
        end if
  7    END do
       dminim(k)=f2
       f=f+f2
       nel(jmin)=nel(jmin)+1
       nk(jmin,nel(jmin))=k
       list1(k)=jmin
       rad(jmin)=rad(jmin)+f2
       radmax(jmin)=dmax1(radmax(jmin),f2)
      end do

      if(nc.gt.1) then
       write(39,*)
       write(39,301) nc
       write(39,*)
       do i=1,nrecord
        k1=list1(i)
        write(39,302) k1
       end do
      end if
301   format('The number of clusters:',i10)
302   format(i8)

      do k=1,nc
       if(nel(k).gt.0) then
          rad(k)=rad(k)/dble(nel(k))
        else 
          rad(k)=0.0d+00
       end if
      end do

      if(nrecord.lt.500) then
       do k=1,nc
        rad(k)=0.0d+00
       end do
      end if

      if((nc.gt.5).and.(nrecord.gt.500)) then
       ratmin=1.0d+28
       do k=1,nc
        if((rad(k).gt.0.0d+00).and.(radmax(k).gt.0.0d+00)) then 
         ratio(k)=radmax(k)/rad(k)
         ratmin=dmin1(ratmin,ratio(k))
        end if 
       end do      
       do k=1,nc
        step1=5.0d-01*ratmin/ratio(k)
        rad(k)=rad(k)+step1*(radmax(k)-rad(k))
       end do
      end if

      if(nc.lt.nclust) then
       ncand=0
       do k=1,nc
        if(nel(k).gt.2) then
         toler3=5.0d-01*rad(k)
         ncand1=0
         do i=1,nel(k)
          i1=nk(k,i)
          if(ncand1.eq.0) then
           if(dminim(i1).gt.rad(k)) then
            ncand=ncand+1
            ncand1=ncand1+1
            lcand(ncand)=i1
            lcand1(ncand1)=i1           
           end if
          end if
          if(ncand1.gt.0) then
           if(dminim(i1).gt.rad(k)) then
            do j=1,ncand1
             j1=lcand1(j)
             d4=0.0d+00
             do k1=1,mf
              d4=d4+(a(i1,k1)-a(j1,k1))**2
             end do
             tnorm=tnorm+1.0d+00
             if(d4.le.toler3) go to 11
            end do
            ncand=ncand+1
            lcand(ncand)=i1
            ncand1=ncand1+1
            lcand1(ncand1)=i1
           end if
          end if
  11     end do
        end if
       end do 
      end if 
c============================================================
      return
      end

c=====================================================================
      subroutine clusters2(x,f)
      PARAMETER(maxdim=4000, maxrec=435000, maxclust=15, maxnft=3,
     1 maxclass=1)
      implicit double precision (a-h,o-z)
      double precision x(maxdim),a(maxrec,maxnft)
      integer nk(maxclust,maxrec),nel(maxclust),list1(maxrec)
      common /c22/a,/cmf/mf,/anclust/nclust,/crecord/nrecord,/cnc/nc
     1 ,/cnk/nk,nel,/ctnorm/tnorm,/clist1/list1
     
      do j=1,nclust
       nel(j)=0
      END do

      f=0.0d+00
      do k=1,nrecord
       f2=1.0d+28
       do j=1,nc
        f1=0.0d+00
        do k1=1,mf
         f1=f1+(a(k,k1)-x(k1+(j-1)*mf))**2
        END do
        tnorm=tnorm+1.0d+00
        if(f2.gt.f1) then
         f2=f1
         jmin=j
        end if
       END do
       f=f+f2
       nel(jmin)=nel(jmin)+1
       nk(jmin,nel(jmin))=k
       list1(k)=jmin
      end do
c============================================================
      return
      end

c=========================================================================
c  Step1 computes the centroid and the value of clustering function at centroid
c=========================================================================
      subroutine step1(f,x)
      PARAMETER(maxdim=4000, maxrec=435000, maxnft=3)
      implicit double precision (a-h,o-z)
      double precision x(maxdim),a(maxrec,maxnft)
      common /c22/a,/cmf/mf,/crecord/nrecord,/ctnorm/tnorm
      do j=1,mf
       x(j)=0.0d+00
       do i=1,nrecord
        x(j)=x(j)+a(i,j)
       END do
       x(j)=x(j)/DBLE(nrecord)
      END do
      f=0.0d+00
      do i=1,nrecord
       f1=0.0d+00
       do j=1,mf
        f1=f1+(a(i,j)-x(j))**2
       END do
       f=f+f1
      END do
      tnorm=tnorm+dble(nrecord)
      return
      end

c=========================================================================
c  Step2 computes clusters for each data point
c=========================================================================
      subroutine step2(nstart,x2)
      PARAMETER(maxdim=4000, maxrec=435000, maxclust=15, maxnft=3
     1 ,maxsize=1000000)
      implicit double precision (a-h,o-z)
      double precision x2(maxsize),a(maxrec,maxnft),dminim(maxrec)
     1 ,fmin1(maxrec),x4(maxnft),fval(maxrec)
      integer lcand(maxrec),lcand1(maxrec),nclose(maxrec)
     1 ,lclose(1000,maxrec)
      common /c22/a,/cmf/mf,/adminim/dminim,/crecord/nrecord
     1 ,/ctnorm/tnorm,/ccand/ncand,lcand,/cgamma/gamma1,gamma2
     2 ,/ctoler/toler

      if(nrecord.le.200) then
       do i1=1,ncand
        i=lcand(i1)
        do k=1,mf
         x2(k+(i1-1)*mf)=a(i,k)
        end do
       end do
       nstart=ncand
       return
      end if

      nstart=0
      fmin=0.0d+00
      fmax=-1.0d+28
      do i1=1,ncand
       i=lcand(i1)
       nclose(i)=0
       d21=0.0d+00
       do l=1,nrecord
        d3=0.0d+00
        tnorm=tnorm+1.0d+00
        do k=1,mf
         d3=d3+(a(i,k)-a(l,k))**2
        end do
        if(d3.lt.dminim(l)) then
         d21=d21+d3-dminim(l)
         nclose(i)=nclose(i)+1
         lclose(i,nclose(i))=l
        end if
       end do
       fmin1(i)=d21
       fmin=dmin1(fmin,d21)
       fmax=dmax1(fmax,d21)
      end do
      
      fmin2=fmin+gamma1*(fmax-fmin)
      ncand1=0
      do i1=1,ncand
       i=lcand(i1)
       if (fmin1(i).le.fmin2) then
        ncand1=ncand1+1
        lcand1(ncand1)=i
       end if
      end do

      ncand=ncand1
      do i1=1,ncand
       lcand(i1)=lcand1(i1)
      end do

      do i1=1,ncand
       i=lcand(i1)
       IF(nclose(i).eq.0) GO TO 1
       do k=1,mf
        d3=0.0d+00
        do j=1,nclose(i)
         j1=lclose(i,j)
         d3=d3+a(j1,k)
        end do
        x4(k)=d3/DBLE(nclose(i))
       end do
       do j=1,nstart
        d4=0.0d+00
        tnorm=tnorm+1.0d+00
        do k=1,mf
         d4=d4+(x2(k+(j-1)*mf)-x4(k))**2
        end do
        IF(d4.le.toler) GO TO 1
       end do
       nstart=nstart+1
       do k=1,mf
        x2(k+(nstart-1)*mf)=x4(k)
       end do
   1  end do
      d2=0.0d+00
      d6=-1.0d+28
      do j=1,nstart
       d21=0.0d+00
       do l=1,nrecord
        d3=0.0d+00
        tnorm=tnorm+1.0d+00
        do k=1,mf
         d3=d3+(x2(k+(j-1)*mf)-a(l,k))**2
        end do
        d21=d21+dmin1(0.0d+00,d3-dminim(l))
       end do
       fval(j)=d21
       d2=dmin1(d2,d21)
       d6=dmax1(d6,d21)
      end do
      
      d2=d2+gamma2*(d6-d2)      

      nstart1=0
      do j=1,nstart
       if (fval(j).le.d2) then
        nstart1=nstart1+1
        do k=1,mf
         x2(k+(nstart1-1)*mf)=x2(k+(j-1)*mf)
        end do
       end if
      end do
      nstart=nstart1
      return
      end

      subroutine auxfunc(x,fval)
      implicit double precision (a-h,o-z)
      PARAMETER(maxdim=4000, maxrec=435000, maxclust=15, maxnft=3)
      double precision x(maxdim),a(maxrec,maxnft),dminim(maxrec)
      common /c22/a,/crecord/nrecord,/cmf/mf,/adminim/dminim
     1 ,/ctnorm/tnorm

      fval=0.0d+00
      do i=1,nrecord
       f1=dminim(i)
       tnorm=tnorm+1.0d+00
       f2=0.0d+00
       do j=1,mf
        f2=f2+(a(i,j)-x(j))**2
       end do
       f3=dmin1(f1,f2)
       fval=fval+f3
      end do
      return
      end

      subroutine clusterfunc(x,f)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=4000, maxrec=435000, maxclust=15, maxnft=3)
      double precision x(maxvar),a(maxrec,maxnft),dclust(maxclust)
     1 , dnorm(maxrec,maxclust)  
      integer l1(maxrec),l2(maxrec,maxclust),l3(maxrec),nel(maxclust)
     1 ,nk(maxclust,maxrec) 
      common /c22/a,/crecord/nrecord,/cnc/nc,/cmf/mf,/ctnorm/tnorm
     1 ,/cl1/l1,/cnversion/nversion,/cl3/l2,l3,/ctau/tau1,tau2
     2 ,/cnel1/nel,nk,/cdnorm/dnorm,/cnphase/nphase

      do k=1,nc
       dclust(k)=0.0d+00
      end do

      f1=0.0d+00
      do i=1,nrecord
       f4=1.0d+28
       do k=1,nc
        f3=0.0d+00
        do j=1,mf
         f3=f3+(a(i,j)-x(j+(k-1)*mf))**2
        end do
        dnorm(i,k)=f3
        if(f4.gt.f3) then
          f4=f3
          l1(i)=k                 
        end if
       end do
       f1=f1+f4
       k1=l1(i)
       dclust(k1)=dclust(k1)+f4
      end do
      f=f1
      tnorm=tnorm+dble(nrecord)*dble(nc)
c=====================================================
      if(nphase.eq.2) then
       do k=1,nc
        nel(k)=0
       end do
      
       do i=1,nrecord
        k1=l1(i)
        nel(k1)=nel(k1)+1
        nk(k1,nel(k1))=i
       end do

       do k=1,nrecord
        l3(k)=0
       end do

       do k=1,nrecord
        do j=1,nc
         l2(k,j)=0
        end do
       end do
c====================================================
       fhat=0.0d+00
       do i=1,nrecord
         i2=l1(i)
         dj=dnorm(i,i2)
         if(nel(i2).gt.1) then
            dj1=(dclust(i2)-dj)/dble(nel(i2)-1)
            dj2=dj-dj1
            if(dj2.ge.0.0d+00) l2(i,i2)=1
           else
            dj1=0.0d+00
            dj2=0.0d+00
         end if
         dj8=0.0d+00
         do j=1,nc
          if(j.ne.i2) then
           dj3=dnorm(i,j)
           dj4=dclust(j)/dble(nel(j))
           dj5=dj3-dj4
           if(dj5.ge.0.0d+00) l2(i,j)=1
           dj6=dmax1(0.0d+00,dj5)
           dj7=dj2-tau1*dj6
           if(dj7.ge.0.0d+00) then
            if(dj8.le.dj7) then
             dj8=dj7
             l3(i)=j
            end if   
           end if
          end if
         end do
         fhat=fhat+dj8
       end do
       f=f+tau2*fhat
      end if
c=========================================================
      return
      end

      subroutine fgrad1(x1,grad1)
       implicit double precision (a-h,o-z)
       PARAMETER(maxdim=4000, maxrec=435000, maxclust=15, maxnft=3)
       double precision a(maxrec,maxnft), grad1(maxdim), x1(maxdim)
       common /c22/a,/crecord/nrecord,/cnc/nc,/cmf/mf,/csize/m,/cns/ns

       do i=1,m
        grad1(i)=0.0d+00
       end do

       if (ns.eq.1) then
        do i=1,nrecord
         do j=1,m
          grad1(j)=grad1(j)+2.0d+00*(x1(j)-a(i,j))
         end do
        END do
       end if

       if (ns.eq.2) then
        do i=1,nrecord
         do j=1,nc
          do k=1,mf
           grad1(k+(j-1)*mf)=grad1(k+(j-1)*mf)
     1          +2.0d+00*(x1(k+(j-1)*mf)-a(i,k))
          end do
         end do
        END do
       end if
      end subroutine

      subroutine fgrad2(x,grad2)
       implicit double precision (a-h,o-z)
       PARAMETER(maxdim=4000, maxrec=435000, maxclust=15, maxnft=3)
       double precision x(maxdim), a(maxrec,maxnft), d2(maxclust)
     1  ,dminim(maxrec), grad2(maxdim)
      common /c22/a,/crecord/nrecord,/cnc/nc,/cmf/mf,/ctnorm/tnorm
     1 ,/csize/m,/cns/ns,/adminim/dminim

      do i=1,m
       grad2(i)=0.0d+00
      end do

      if (ns.eq.1) then
        do i=1,nrecord
         tnorm=tnorm+1.0d+00
         f3=0.0d+00
         do j=1,mf
          f3=f3+(a(i,j)-x(j))**2
         end do
         if (f3.ge.dminim(i)) then
          do j=1,m
           grad2(j)=grad2(j)+2.0d+00*(x(j)-a(i,j))
          end do
         end if
        END do
      end if

      if (ns.eq.2) then
        do i=1,nrecord
         do j=1,nc
          tnorm=tnorm+1.0d+00
          d2(j)=0.0d+00
          do k=1,mf
           d2(j)=d2(j)+(a(i,k)-x(k+(j-1)*mf))**2
          end do
         end do
         d3=0.0d+00
         do j=1,nc
          d4=0.0d+00
          do k=1,nc
           if (k.ne.j) then
            d4=d4+d2(k)
           end if
          end do
          if (d3.lt.d4) then
           d3=d4
           jindex=j
          end if
         end do
         do j=1,nc
          if (j.ne.jindex) then
           do k=1,mf
            grad2(k+(j-1)*mf)=grad2(k+(j-1)*mf)
     1       +2.0d+00*(x(k+(j-1)*mf)-a(i,k))
           end do
          end if
         end do
        end do
       end if
      end subroutine

c=====================================================================
      subroutine funcgrad(x,grad)
      PARAMETER(maxvar=4000, maxrec=435000, maxclust=15, maxnft=3)
      implicit double precision (a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),grad(maxvar)
     1 ,dnorm(maxrec,maxclust), da(maxvar) 
      integer l1(maxrec),l2(maxrec,maxclust),l3(maxrec),nel(maxclust)
     1 ,nk(maxclust,maxrec) 
      common /c22/a,/cmf/mf,/crecord/nrecord,/cnc/nc,/csize/m,/cl1/l1
     1 ,/cl3/l2,l3,/ctau/tau1,tau2,/cdnorm/dnorm,/cnel1/nel,nk

      call clusterfunc(x,f)

      do i=1,m
       grad(i)=0.0d+00
      end do
            
      do i=1,nrecord
       k=l1(i)
       do j=1,mf
        grad(j+(k-1)*mf)=grad(j+(k-1)*mf)+2.0d+00*(x(j+(k-1)*mf)-a(i,j))
       end do
      end do
c==================================================================      
      do k=1,m
       da(k)=0.0d+00
      end do
      do i=1,nrecord
       i2=l1(i)
       if(l3(i).gt.0) then
        if(l2(i,i2).eq.1) then
         do j=1,mf
          da(j+(i2-1)*mf)=da(j+(i2-1)*mf)+2.0d+00*(x(j+(i2-1)*mf)
     1     -a(i,j))
         end do          
         do j2=1,nel(i2)
          j1=nk(i2,j2)
          if(j1.ne.i) then
           do j=1,mf
            da(j+(i2-1)*mf)=da(j+(i2-1)*mf)-2.0d+00*(x(j+(i2-1)*mf)
     1       -a(j1,j))/dble(nel(i2)-1)
           end do
          end if          
         end do
        end if

        k1=l3(i)
        if(l2(i,k1).eq.1) then
         do j=1,mf
          da(j+(k1-1)*mf)=da(j+(k1-1)*mf)-2.0d+00*tau1*(x(j+(k1-1)*mf)
     1     -a(i,j))
         end do          
         do j2=1,nel(k1)
          j1=nk(k1,j2)
          do j=1,mf
           da(j+(k1-1)*mf)=da(j+(k1-1)*mf)+2.0d+00*tau1*(x(j+(k1-1)*mf)
     1       -a(j1,j))/dble(nel(k1))
          end do
         end do
        end if
       end if
      end do  
      do k=1,m
        grad(k)=grad(k)+tau2*da(k)
      end do 
c==================================================================      
      return
      end

c=====================================================================
      SUBROUTINE optim(x0,x,fvalue)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=4000, maxrec=435000, maxclust=15)
      double precision x0(maxvar),x(maxvar)
      COMMON /cmf/mf,/cns/ns,/cnc/nc
c=======================================================
c  n - number of variables
c=======================================================
      IF(ns.eq.1) n = mf
      IF(ns.eq.2) n = nc*mf
c=======================================================
      do j=1,n
       x(j)=x0(j)
      end do
      call qsm(n,x,fvalue)
c=======================================================
      RETURN
      END

c======================================================================
      subroutine func(x,objf)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=4000)
      double precision x(maxvar)
      COMMON /cns/ns
c===================================================================
      IF(ns.eq.1) call auxfunc(x,f)
      IF(ns.eq.2) call clusterfunc(x,f)
      objf=f
c===================================================================
      return
      end

c======================================================================
c Without scaling
c======================================================================
      subroutine qsm(nvar,x,f)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=4000)
      double precision x(maxvar)
      COMMON /csize/m,/citer/maxiter,niter,/cnf/nf
      m=nvar
      nf=0
      niter=0
      maxiter=5000
      call optimum(x)
      call fv(x,f)
c===========================================================
      return
      end

c============================================================
      subroutine optimum(x)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=4000, maxdg=1000, maxit=10000)
      double precision x(maxvar),x1(maxvar),g(maxvar),v(maxvar)
     1 ,w(maxdg,maxvar),prod(maxdg,maxdg),z(maxdg),fvalues(maxit)
     2 ,dg3(maxvar)
      INTEGER ij(maxdg)
      common /csize/m,/citer/maxiter,niter,/cij/ij,jvertex,/cz/z
     1 ,/ckmin/kmin,/cns/ns,/cdg3/dg3,/cnphase/nphase
c====================================================================
      dist1=1.0d-07
      step0=-2.0d-01
      div=1.0d-01
      eps0=1.0d-07
      slinit=1.0d+00
      slmin=1.0d-05
      sdif=1.0d-05
      mturn=4
      nbundle=MIN(17,m+3)
c====================================================================
      sl=slinit/div
      call fv(x,f2)
  1   sl=div*sl
      IF(sl.lt.slmin) return
      do i=1,m
       g(i)=1.0d+00/dsqrt(DBLE(m))
      end do
      nnew=0
c================================================================
   2  niter=niter+1
        IF(niter.gt.maxiter) RETURN
        nnew=nnew+1
        f1=f2
        fvalues(niter)=f1
c---------------------------------------------------------------
        if (nnew.gt.mturn) then
         mturn2=niter-mturn+1
         ratio1=(fvalues(mturn2)-f1)/(dabs(f1)+1.0d+00)
         IF(ratio1.LT.sdif) GO TO 1
        end if
        if (nnew.GE.(2*mturn)) then
         mturn2=niter-2*mturn+1
         ratio1=(fvalues(mturn2)-f1)/(dabs(f1)+1.0d+00)
         IF(ratio1.LT.(1.0d+01*sdif)) GO TO 1
        end if
c--------------------------------------------------------------
        do ndg=1,nbundle
            call dgrad(ndg,x,sl,g,v)
            dotprod=0.0d+00
            do i=1,m
             dotprod=dotprod+v(i)*v(i)
            end do
            r=dsqrt(dotprod)
            IF(r.lt.eps0) GO TO 1
            IF(ndg.eq.1) then
                         rmean=r
                         kmin=1
                         rmin=r
            END if
            IF(ndg.gt.1) then
                         rmin=dmin1(rmin,r)
                         IF(r.eq.rmin) kmin=ndg
                         rmean=((ndg-1)*rmean+r)/ndg
            END if
            toler=dmax1(eps0,dist1*rmean)
            do i=1,ndg-1
             prod(ndg,i)=0.0d+00
             do j=1,m
              prod(ndg,i)=prod(ndg,i)+w(i,j)*v(j)
             end do
             prod(i,ndg)=prod(ndg,i)
            end do
            prod(ndg,ndg)=dotprod
c====================================================================
            do i=1,m
             w(ndg,i)=v(i)
            end do
            call wolfe(ndg,prod)
c================================
            do i=1,m
             v(i)=0.0d+00
             do j=1,jvertex
              v(i)=v(i)+w(ij(j),i)*z(j)
             END do
            END do
c================================
            r=0.0d+00
            do i=1,m
             r=r+v(i)*v(i)
            end do
            r=dsqrt(r)
            if(r.lt.toler) GO TO 1
c===========================================================
             do i=1,m
              g(i)=-v(i)/r
              x1(i)=x(i)+sl*g(i)
             end do
c===========================================================
             call fv(x1,f4)
             f3=(f4-f1)/sl
             decreas=step0*r
             if(f3.lt.decreas) then
                        call armijo(x,g,f1,f5,f4,sl,step,r)
                        f2=f5
                        do i=1,m
                         x(i)=x(i)+step*g(i)
                        end do
                        sl=1.2d+00*sl
                        GO TO 2
             end if
         END do
c=====================================================
      go to 1
      return
      end

c=====================================================================
c Subroutine dgrad calculates subgradients or discrete gradients
c=====================================================================
      subroutine dgrad(ndg,x,sl,g,dg)
       implicit double precision (a-h,o-z)
       PARAMETER(maxvar=4000, maxdg=1000)
       double precision x1(maxvar),g(maxvar),x(maxvar),dg(maxvar)
     1  ,dg1(maxvar),dg2(maxvar),dg3(maxvar)
       common /csize/m,/cdg3/dg3,/cnphase/nphase

       do k=1,m
         x1(k)=x(k)+sl*g(k)
       end do
       
       if(nphase.eq.1) then
        call fgrad1(x1,dg1)
        if(ndg.eq.1) then
         call fgrad2(x,dg2)
         do i=1,m
          dg3(i)=dg2(i)
         end do
        end if
        do i=1,m
         dg(i)=dg1(i)-dg3(i)
        end do
       end if
       if(nphase.eq.2) then
        call funcgrad(x1,dg)
       end if
      end subroutine

c===========================================================
c Line search (Armijo-type), Step 5 Algorithm 2.
c===========================================================
      subroutine armijo(x,g,f1,f5,f4,sl,step,r)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=4000, maxdg=1000)
      common /csize/m
      double precision x(maxvar),g(maxvar),x1(maxvar)
      step=sl
      f5=f4
      s1=sl
      k=0
  1   k=k+1
      IF(k.gt.20) RETURN
      s1=2.0d+00*s1
      do i=1,m
       x1(i)=x(i)+s1*g(i)
      end do
      call fv(x1,f50)
      f30=f50-f1+5.0d-02*s1*r
      IF(f30.gt.0.0d+00) RETURN
      step=s1
      f5=f50
      GO TO 1
      end subroutine

      subroutine fv(x,f)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=4000)
      double precision x(maxvar)
      COMMON /cnf/nf
c================================================================
      call func(x,f)
      nf=nf+1
c=================================================================
      end subroutine

c====================================================================
c  Subroutines Wolfe and Equations solves quadratic
c  programming problem, to find
c  descent direction, Step 3, Algorithm 2.
c===============================================================

      subroutine wolfe(ndg,prod)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=4000, maxdg=1000)
      common /csize/m,/w01/a,/cij/ij,jvertex,/cz/z,/ckmin/kmin
      INTEGER ij(maxdg)
      double precision z(maxdg),z1(maxdg),a(maxdg,maxdg)
     1 ,prod(maxdg,maxdg)
      j9=0
      jmax=500*ndg
      jvertex=1
      ij(1)=kmin
      z(1)=1.d+00
c=======================================
c  To calculate X
c=======================================
 1    r=0.0d+00
      do i=1,jvertex
       do j=1,jvertex
        r=r+z(i)*z(j)*prod(ij(i),ij(j))
       end do
      end do
      IF(ndg.eq.1) RETURN
c========================================
c  To calculate <X,P_J> and J
c========================================
      t0=1.0d+28
      do i=1,ndg
        t1=0.0d+00
        do j=1,jvertex
          t1=t1+z(j)*prod(ij(j),i)
        end do
        if(t1.lt.t0) then
                     t0=t1
                     kmax=i
        end if
      end do
c========================================
c  First stopping criterion
c========================================
      rm=prod(kmax,kmax)
      do j=1,jvertex
       rm=dmax1(rm,prod(ij(j),ij(j)))
      end do
      r2=r-1.d-12*rm
      if(t0.gt.r2) RETURN
c========================================
c  Second stopping criterion
c========================================
      do i=1,jvertex
       if(kmax.eq.ij(i)) RETURN
      end do
c========================================
c Step 1(e) from Wolfe's algorithm
c========================================
      jvertex=jvertex+1
      ij(jvertex)=kmax
      z(jvertex)=0.0d+00
c========================================
 2    do i=1,jvertex
       do j=1,jvertex
        a(i,j)=1.0d+00+prod(ij(i),ij(j))
       end do
      end do
      j9=j9+1
      if(j9.gt.jmax) RETURN
      call equations(jvertex,z1)
      do i=1,jvertex
       if(z1(i).le.1.0d-10) go to 3
      end do
      do i=1,jvertex
       z(i)=z1(i)
      end do
      go to 1
  3   teta=1.0d+00
      do i=1,jvertex
       z5=z(i)-z1(i)
       if(z5.gt.1.0d-10) teta=dmin1(teta,z(i)/z5)
      end do
      do i=1,jvertex
       z(i)=(1.0d+00-teta)*z(i)+teta*z1(i)
       if(z(i).le.1.0d-10) then
                          z(i)=0.0d+00
                          kzero=i
       end if
      end do
      j2=0
      do i=1,jvertex
       IF(i.ne.kzero) then
                     j2=j2+1
                     ij(j2)=ij(i)
                     z(j2)=z(i)
       END if
      end do
      jvertex=j2
      go to 2
      return
      end

      subroutine equations(n,z1)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=4000, maxdg=1000)
      common /w01/a
      double precision a(maxdg,maxdg),z1(maxdg),b(maxdg,maxdg)
      do i=1,n
       do j=1,n
        b(i,j)=a(i,j)
       end do
       b(i,n+1)=1.0d+00
      end do
      do i=1,n
       r=b(i,i)
       do j=i,n+1
        b(i,j)=b(i,j)/r
       end do
       do j=i+1,n
        do k=i+1,n+1
         b(j,k)=b(j,k)-b(i,k)*b(j,i)
        end do
       end do
      end do
      z1(n)=b(n,n+1)
      do i=1,n-1
        k=n-i
        z1(k)=b(k,n+1)
        do j=k+1,n
         z1(k)=z1(k)-b(k,j)*z1(j)
        END do
      end do
      z2=0.d+00
      do i=1,n
       z2=z2+z1(i)
      end do
      do i=1,n
       z1(i)=z1(i)/z2
      end do
      return
      end

c==================================================================
      subroutine silhouette(nneg,ratio,pen)
      PARAMETER(maxrec=435000, maxclust=15, maxnft=3)
      implicit double precision(a-h,o-z)
      double precision a(maxrec,maxnft),sa(maxclust),sil(maxrec)
      integer list1(maxrec),nel(maxclust)
      common /c22/a,/cmf/mf,/crecord/nrecord,/clist1/list1,/cnc/nc
c====================================================================
      do k=1,nc
       nel(k)=0
      end do
      do i=1,nrecord
       k=list1(i)
       nel(k)=nel(k)+1
      end do
c====================================================================
      ns=0
      nneg=0
      pen=0.0d+00
      do i=1,nrecord
        k1=list1(i)
        do k=1,nc
         sa(k)=0.0d+00
        end do
        do j=1,nrecord
         if(j.ne.i) then
          d1=0.0d+00
          do j2=1,mf
           d1=d1+(a(i,j2)-a(j,j2))**2
          end do
         ! d1=dsqrt(d1) ! commented by Sona 31/08/21 *******************
          k2=list1(j)
          sa(k2)=sa(k2)+d1
         end if
        end do
        if(nel(k1).gt.1) then
           da1=sa(k1)/dble(nel(k1)-1)
          else
           da1=0.0d+00
        end if
        da2=1.0d+28
        do k=1,nc
         if(k.ne.k1) then
          if(nel(k).gt.0) then
            da3=sa(k)/dble(nel(k))
            da2=dmin1(da2,da3)
          end if  
         end if
        end do 
        sil(i)=(da2-da1)/dmax1(da1,da2)
        if(sil(i).gt.0.0d+00) ns=ns+1
        if(sil(i).lt.0.0d+00) then
         nneg=nneg+1
         pen=pen+sil(i)
        end if 
      end do
      ratio=dble(ns)/dble(nrecord)
c====================================================================
      return
      end 

c==================================================================
      subroutine finresult(x,f,purity,db)
      PARAMETER(maxvar=4000, maxrec=435000, maxclust=15, maxnft=3
     1 ,maxclass=1)
      implicit double precision(a-h,o-z)
      double precision a(maxrec,maxnft),rad(maxclust)
     1 ,dminim(maxrec),plabel(maxclass),fk(maxclust),x(maxvar)
      integer nel(maxclust),nk(maxclust,maxrec),list1(maxrec)
     1 ,nob(maxclass),nob1(maxclass)
      common /c22/a,/cnft/nft,/cmf/mf,/cnclass/nclass,/crecord/nrecord
     1 ,/clist1/list1,/cnob/plabel,nob,/cnpnout/npurity,/cnk/nk,nel
     2 ,/cnc/nc,/crad/rad,/adminim/dminim
c====================================================================
      f=0.0d+00
      do k=1,nrecord
       f=f+dminim(k)
      end do
c====================================================================
c Purity
c====================================================================
      if (npurity.eq.1) then
       n3=0
       do i=1,nc
        do j=1,nclass
         nob1(j)=0
        end do
        do j=1,nel(i)
         j1=nk(i,j)
         do k=1,nclass
          c1=plabel(k)
          IF(a(j1,nft).EQ.c1) nob1(k)=nob1(k)+1
         end do
        end do
        n2=nob1(1)
        do k=2,nclass
         if (n2.lt.nob1(k)) then
          n2=nob1(k)
         end if
        end do
        n3=n3+n2
       end do
       purity=1.0d+02*dble(n3)/dble(nrecord)
      end if
      if(npurity.eq.2) purity=0.0d+00
c====================================================================
c Davies-Bouldin (DB) validity index
c=====================================================
      do i=1,nc 
       fk(i)=0.0d+00
      end do

      do i=1,nrecord
       k=list1(i)
       fk(k)=fk(k)+sqrt(dminim(i))
      end do
      do i=1,nc
       IF(nel(i).gt.0) fk(i)=fk(i)/DBLE(nel(i))
      end do
      fdb=0.0d+00
      do k=1,nc
       fm=0.0d+00
       do i=1,nc
        if (i.ne.k) then
         fk2=fk(i)+fk(k)
         f12=0.0d+00
         do j=1,mf
          f12=f12+(x(j+(i-1)*mf)-x(j+(k-1)*mf))**2
         end do
         f12=dsqrt(f12)
         f2=fk2/f12
         fm=dmax1(fm,f2)
        end if
       end do
        fdb=fdb+fm
      end do
      db=fdb/DBLE(nc)
c====================================================================
      return
      end 
