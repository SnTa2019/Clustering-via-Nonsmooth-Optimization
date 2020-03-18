c=============================================================
c DG-CLUST: discrete gradient clustering algorithm 
c
c The name of input file should be: datainput.txt
c
c The name of output file is: DGMResults.txt
c
c=============================================================
c     main programm
      PARAMETER(maxvar=2000, maxrec=200000, maxclust=100, maxnft=100
     1 ,maxclass=100, maxsize=7000000)
      implicit double precision(a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),dminim(maxrec)
     1,x2(maxsize),amed(maxnft),plabel(maxclass),xbest(maxvar),z(maxvar)
     2 ,fval(maxrec),dcent(maxclust,maxclust),x5(maxvar),x3(maxvar)
      integer nk(maxclust,maxrec),nel(maxclust),nob(maxclass)
      common /csize/m,/c22/a,/anclust/nclust,/cnft/nft,/cmf/mf
     1 ,/crecord/nrecord,/cnc/nc,/cnk/nk,nel,/adminim/dminim,/cns/ns
     2 ,/cnob/plabel,nob,/ctnorm/tnorm,/cind/ind,/cgamma/gamma1,gamma2
     3 ,/cnpnout/npurity,/cnclass/nclass
      open(40,file='DGMResults.txt')
      open(78,file='datainput.txt',status='old',form='formatted')
c========================================================================
      PRINT *,' '
      PRINT *,'Number of features:'
      read *,nft
      PRINT *,' '
      PRINT *,'Class outputs?'
      PRINT *,' '
      PRINT *,'   1 - if yes'
      PRINT *,'   2 - if no'
      read *,npurity
      PRINT *,' '
      if (npurity.eq.1) then
         PRINT *,'Number of classes in your dataset:'
         read *,nclass
         PRINT *,' '
         PRINT *,'Output column:'
         read *,noutcom
         PRINT *,' '
      end if
      PRINT *,'Maximum number of clusters:'
      read *,nclust
      PRINT *,' '
      PRINT *,'Compute indices (DB, Dunn, Purity, Silhouettes)?'
      PRINT *,' '
      PRINT *,'   1 - yes'
      PRINT *,'   2 - no'
      read *,ind
      PRINT *,' ' 
c===========================================================
      tlimit=7.2d+04
      do i=1,maxrec
       read(78,*,END=901,ERR=900) (a(i,k),k=1,nft)
       nrecord=i
      end do
  900 stop 'Error in input file'       
  901 WRITE(40,*) 'Input complete. Number of records: ',nrecord
      WRITE(40,*)
      
      IF(npurity.eq.1) mf=nft-1
      IF(npurity.eq.2) mf=nft
      IF(npurity.eq.2) noutcom=0
      tnorm=0.0d+00
c==============================================================
      if ((npurity.eq.1).AND.(noutcom.gt.0)) then
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
       a2=a(1,nft)
       plabel(1)=a2
       nob(1)=1
       n2=1
       do i=2,nrecord
        clabel=a(i,nft)
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
c=================================================================
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
       gamma1=5.0d-01
       gamma2=3.0d-01
       gamma3=2.0d-01
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
      call cpu_time(time1)
c===================================================================
       WRITE(40,*)
       WRITE(40,*)
       WRITE(40,*) '___________________________________________________'
       WRITE(40,*)
       WRITE(40,701)
 701   FORMAT('#Clusters','       Function value','    DB index'
     1 ,'   Dunn index','   Purity','  #+Silhouet'
     2 ,'   #Dist_f_eval','      CPU')
       WRITE(40,*)
c===================================================================
      do nc=1,nclust
       PRINT 42,nc
 42    FORMAT('Cluster No.:',i10)
       if(nc.eq.1) then
                   call step1(f,x)
                   toler=1.0d-02*f/dble(nrecord)
                   go to 1
       END if
       call step2(toler,nstart,x2)
       m=mf
       fbarmin=1.0d+26
       fbarmax=0.0d+00
       do j=1,nstart
        do k=1,mf
         z(k)=x2(k+(j-1)*mf)
        end do
        ns=1
        call dgm(z,barf)
        fval(j)=barf
        fbarmin=dmin1(fbarmin,barf)
        fbarmax=dmax1(fbarmax,barf)
        do k=1,mf
         x2(k+(j-1)*mf)=z(k)
        end do
       end do
c==============================================
       fbarmin=fbarmin+gamma3*(fbarmax-fbarmin)
       nstart1=0
       do j=1,nstart
        if (fval(j).le.fbarmin) then
         nstart1=nstart1+1
         do k=1,mf
          x5(k+(nstart1-1)*mf)=x2(k+(j-1)*mf)
         end do
        end if
       end do

       nstart=nstart1
       do i=1,nstart
        do k=1,mf
         x2(k+(i-1)*mf)=x5(k+(i-1)*mf)
        end do
       end do
c==============================================
       do k=1,mf
        x5(k)=x2(k)
       end do
       nstart2=1
       do j=2,nstart
        do j1=1,nstart2
         f31=0.0d+00
         do k=1,mf
          f31=f31+(x5(k+(j1-1)*mf)-x2(k+(j-1)*mf))**2
         end do
         IF(f31.LE.toler) GO TO 1200
        end do
        nstart2=nstart2+1
        do k=1,mf
         x5(k+(nstart2-1)*mf)=x2(k+(j-1)*mf)
        end do
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
        call dgm(x3,fcurrent)
        if (fcurrent.lt.fbest) then
         fbest=fcurrent
         do j1=1,m
          xbest(j1)=x3(j1)
         end do
        end if
       end do
       f=fbest
       do j1=1,m
        x(j1)=xbest(j1)
       end do
       
       do i=1,nc
        dcent(i,i)=0.0d+00
       end do
       
       do i=1,nc
        do j=i+1,nc
         dcent(i,j)=0.0d+00
         do j1=1,mf
          dcent(i,j)=dcent(i,j)+(x(j1+(i-1)*mf)-x(j1+(j-1)*mf))**2
         end do
         dcent(j,i)=dcent(i,j)
        end do
       end do
c================================================================
  1    call clusters(x,f,db,dunn,purity,nsil)
       if(npurity.eq.2) purity=0.0d+00

       call cpu_time(time3)
       time4=time3-time1
       write(40,603) nc,f,db,dunn,purity,nsil,tnorm,time4
 603   format(i6,f24.4,2f12.6,f10.2,i10,f18.0,f11.3)

       IF(time4.gt.tlimit) GO TO 2
      end do
   2  continue
      close(40)
      close(78)
      stop
      end
c================================================================
      subroutine clusters(x,f,db,dunn,purity,ns)
      PARAMETER(maxdim=2000, maxrec=200000, maxclust=100, maxnft=100,
     1 maxclass=100)
      implicit double precision (a-h,o-z)
      double precision x(maxdim),a(maxrec,maxnft),dminim(maxrec)
     1 ,fk(maxclust),plabel(maxclass),rad(maxclust),radmax(maxclust)
     2 ,dcent(maxclust,maxclust),sa(maxclust),sil(maxrec)
     3 ,ratio(maxclust)
      integer nk(maxclust,maxrec),nel(maxclust),nob(maxclass)
     1 ,nob1(maxclass),lcand(maxrec),list1(maxrec),lcand1(maxrec)
      common /c22/a,/cmf/mf,/anclust/nclust,/crecord/nrecord
     1 ,/cnc/nc,/cnk/nk,nel,/cnob/plabel,nob,/cnclass/nclass
     2 ,/cnpnout/npurity,/ccand/ncand,lcand,/adminim/dminim
     3 ,/cdcent/dcent,/ctnorm/tnorm,/clist1/list1,/cind/ind
      do j=1,nclust
       nel(j)=0
       rad(j)=0.0d+00
       radmax(j)=0.0d+00
      END do

      f=0.0d+00
      do k=1,nrecord
       f2=1.0d+22
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
       ratmin=1.0d+26
       do k=1,nc
        ratio(k)=radmax(k)/rad(k)
        ratmin=dmin1(ratmin,ratio(k))
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
      if (npurity.eq.1) then
       n3=0
       do i=1,nc
        do j=1,nclass
         nob1(j)=0
        end do
        do j=1,nel(i)
         j1=nk(i,j)
         do k=1,nclass
          IF(a(j1,mf+1).EQ.plabel(k)) nob1(k)=nob1(k)+1
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
c=====================================================
      if(ind.eq.2) return
c=====================================================
c Calculation of Davies-Bouldin (DB) validity index
c=====================================================
      do i=1,nc
       fk(i)=0.0d+00
      end do

      fdb=0.0d+00
      do i=1,nrecord
       k=list1(i)
       fk(k)=fk(k)+dminim(i)
      end do
      do i=1,nc
       IF(nel(i).gt.0) fk(i)=fk(i)/DBLE(nel(i))
      end do
      do k=1,nc
       fm=0.0d+00
       do i=1,nc
        if (i.ne.k) then
         fk2=fk(i)+fk(k)
         f12=0.0d+00
         do j=1,mf
          f12=f12+(x(j+(i-1)*mf)-x(j+(k-1)*mf))**2
         end do
         f2=fk2/f12
         fm=dmax1(fm,f2)
        end if
       end do
       fdb=fdb+fm
      end do
      db=fdb/DBLE(nc)
c============================================================
c Calculation of Dunn validity index
c============================================================
      if(nc.eq.1) dunn=1.0d+00
      if(nc.gt.1) then
       dn=1.0d+26
       do i=1,nc
        dn2=1.0d+26
        do j=i+1,nc
         dm=1.0d+26
         do j1=1,nel(i)
          k1=nk(i,j1)
          do j2=1,nel(j)
           k2=nk(j,j2)
           dm1=0.0d+00
           do j3=1,mf
            dm1=dm1+(a(k1,j3)-a(k2,j3))**2
           end do
           dm=dmin1(dm,dm1)
          end do
         end do
         dn2=dmin1(dn2,dm)
        end do
        dn=dmin1(dn,dn2)
       end do

       dn3=0.0d+00
       do i=1,nc
        diam=0.0d+00
        do j=1,nel(i)
         j1=nk(i,j)
         do j2=j+1,nel(i)
          j3=nk(i,j2)
          d4=0.0d+00
          do j4=1,mf
           d4=d4+(a(j1,j4)-a(j3,j4))**2
          end do         
          diam=dmax1(diam,d4)
         end do
        end do
        dn3=dmax1(dn3,diam)
       end do
       dunn=dn/dn3
      end if
c============================================================
c Calculation of silhouette
c============================================================
      ns=0
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
         k2=list1(j)
         sa(k2)=sa(k2)+d1
        end if
       end do
       da1=sa(k1)/dble(nel(k1))
       da2=1.0d+26
       do k=1,nc
        if(k.ne.k1) then
         da3=sa(k)/dble(nel(k))
         da2=dmin1(da2,da3)
        end if
       end do 
       sil(i)=(da2-da1)/dmax1(da1,da2)
       if(sil(i).gt.0.0d+00) ns=ns+1
      end do
c======================================================================
      return
      end

c=====================================================
c Step1 calculates the center of the dataset
c=====================================================
      subroutine step1(f,x)
      PARAMETER(maxvar=2000, maxrec=200000, maxclust=100, maxnft=100)
      implicit double precision (a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft)
      common /c22/a,/cmf/mf,/crecord/nrecord
      do i=1,mf
       x(i)=0.0d+00
       do j=1,nrecord
        x(i)=x(i)+a(j,i)
       END do
       x(i)=x(i)/dble(nrecord)
      END do
      f=0.0d+00
      do i=1,nrecord
       f1=0.0d+00
       do j=1,mf
        f1=f1+(a(i,j)-x(j))*(a(i,j)-x(j))
       END do
       f=f+f1
      END do
      return
      end

c=========================================================================
c  Step2 computes clusters for each data point
c=========================================================================
      subroutine step2(toler,nstart,x2)
      PARAMETER(maxdim=2000, maxrec=200000, maxclust=100, maxnft=100
     1 ,maxsize=7000000)
      implicit double precision (a-h,o-z)
      double precision x2(maxsize),a(maxrec,maxnft),dminim(maxrec)
     1 ,fmin1(maxrec),x4(maxnft),fval(maxrec)
      integer l4(maxrec),lcand(maxrec),lcand1(maxrec)
      common /c22/a,/cmf/mf,/adminim/dminim,/crecord/nrecord
     1 ,/ctnorm/tnorm,/ccand/ncand,lcand,/cgamma/gamma1,gamma2

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
      fmax=-1.0d+26
      do i1=1,ncand
       i=lcand(i1)
       d21=0.0d+00
       do l=1,nrecord
        d3=0.0d+00
        tnorm=tnorm+1.0d+00
        do k=1,mf
         d3=d3+(a(i,k)-a(l,k))**2
        end do
        d21=d21+dmin1(0.0d+00,d3-dminim(l))
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
       nclose=0
       do j=1,nrecord
        d1=0.0d+00
        tnorm=tnorm+1.0d+00
        do k=1,mf
         d1=d1+(a(i,k)-a(j,k))**2
        end do
        IF(d1.LT.dminim(j)) then
         nclose=nclose+1
         l4(nclose)=j
        END if
       end do

       IF(nclose.eq.0) GO TO 1
       do k=1,mf
        d3=0.0d+00
        do j=1,nclose
         j1=l4(j)
         d3=d3+a(j1,k)
        end do
        x4(k)=d3/DBLE(nclose)
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
      d6=-1.0d+26
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
c================================================================
      subroutine auxfunc(f,x)
      PARAMETER(maxvar=2000, maxrec=200000, maxclust=100, maxnft=100)
      implicit double precision (a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),dminim(maxrec)
     1 ,fval2(maxrec)
      common /c22/a,/cmf/mf,/crecord/nrecord,/adminim/dminim
     1 ,/cval2/fval2,/ctnorm/tnorm
       f=0.0d+00
       do i=1,nrecord
        f4=dminim(i)
        f3=0.0d+00
        do j=1,mf
         f3=f3+(a(i,j)-x(j))**2
        END do
        tnorm=tnorm+1.0d+00
        fval2(i)=f3
        f4=dmin1(f4,f3)
        f=f+f4
      end do
      return
      end
c================================================================
      subroutine auxgrad(x,grad,h)
      PARAMETER(maxvar=2000, maxrec=200000, maxclust=100, maxnft=100)
      implicit double precision (a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),dminim(maxrec)
     1 ,fval2(maxrec),grad(maxvar)
      common /c22/a,/cmf/mf,/crecord/nrecord,/adminim/dminim
     1 ,/cval2/fval2,/cicor/icor
      f1=0.0d+00
      do i=1,nrecord
       f1=f1+dmin1(dminim(i),fval2(i))
      end do
      do k=1,mf
       if(k.ne.icor) then
        f2=0.0d+00
        do i=1,nrecord
         fval2(i)=fval2(i)+2.0d+00*h*(x(k)-a(i,k))+h**2
         f2=f2+dmin1(dminim(i),fval2(i))
        end do
        grad(k)=(f2-f1)/h
        f1=f2
       end if
      end do
      return
      end
c================================================================
      subroutine func(f,x)
      PARAMETER(maxvar=2000, maxrec=200000, maxclust=100, maxnft=100)
      implicit double precision (a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),fval1(maxrec,maxclust)
      common /c22/a,/cmf/mf,/crecord/nrecord,/cnc/nc,/cval1/fval1
     1 ,/ctnorm/tnorm 
       f=0.0d+00
       do i=1,nrecord
        f4=1.0d+26
        do k=1,nc
         f3=0.0d+00
         do j=1,mf
          f3=f3+(a(i,j)-x(j+(k-1)*mf))**2
         END do
         tnorm=tnorm+1.0d+00
         fval1(i,k)=f3
         f4=dmin1(f4,f3)                 
        end do
        f=f+f4
       end do
      return
      end
c================================================================
      subroutine funcgrad(x,grad,h)
      PARAMETER(maxvar=2000, maxrec=200000, maxclust=100, maxnft=100)
      implicit double precision (a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),fval1(maxrec,maxclust)
     1 ,grad(maxvar) 
      common /c22/a,/cmf/mf,/crecord/nrecord,/cval1/fval1,/cnc/nc
     1 ,/csize/m,/cicor/icor
      f1=0.0d+00
      do i=1,nrecord
       f3=1.0d+26
       do k=1,nc
        f3=dmin1(f3,fval1(i,k))
       end do
       f1=f1+f3
      end do
      do k=1,m
       if(k.ne.icor) then
        k1=(k-1)/mf+1
        k2=k-(k1-1)*mf
        do i=1,nrecord
         fval1(i,k1)=fval1(i,k1)+2.0d+00*h*(x(k)-a(i,k2))+h**2
        end do
        f2=0.0d+00
        do i=1,nrecord
         f3=1.0d+26
         do k3=1,nc
          f3=dmin1(f3,fval1(i,k3))
         end do
         f2=f2+f3
        end do
        grad(k)=(f2-f1)/h
        f1=f2
       end if
      end do
      return
      end

c=====================================================================
      subroutine dgm(x,f2)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=2000, maxdg=1000,maxit=100000)
      double precision x(maxvar),x1(maxvar),g(maxvar),v(maxvar)
     1 ,w(maxdg,maxvar),prod(maxdg,maxdg),z(maxdg),fvalues(maxit)
      INTEGER ij(maxdg)
      common /csize/m,/cij/ij,jvertex,/cz/z,/ckmin/kmin
c====================================================================
      dist1=1.0d-07
      step0=-2.0d-01
      div=1.0d-01
      eps0=1.0d-07
      slinit=1.0d+00
      slmin=1.d-05*slinit
      pwt=1.0d-06
      sdif=1.0d-05
      mturn=4
      maxiter=5000
      niter=0
      nbundle=min(m+3,10)
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
         ratio1=(fvalues(mturn2)-f1)/(dabs(f1)+1.d+00)
         IF(ratio1.LT.sdif) GO TO 1
        end if
        if (nnew.GE.(2*mturn)) then
         mturn2=niter-2*mturn+1
         ratio1=(fvalues(mturn2)-f1)/(dabs(f1)+1.d+00)
         IF(ratio1.LT.(1.d-01*sdif)) GO TO 1
        end if
c--------------------------------------------------------------
        do ndg=1,nbundle
            call dgrad(x,sl,g,v,f1,f4,ndg,pwt)
            dotprod=0.d+00
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
             prod(ndg,i)=0.d+00
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
             v(i)=0.d+00
             do j=1,jvertex
              v(i)=v(i)+w(ij(j),i)*z(j)
             END do
            END do
c================================
            r=0.d+00
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
                        IF(ndg.le.2) sl=1.1d+00*sl
                        GO TO 2
             end if
         END do
c=====================================================
      go to 1
      return
      end

c==============================================================
c  Subroutines Wolfe and Equations solves quadratic programming problem, 
c  to find descent direction, Step 3, Algorithm 2.
c===============================================================

      subroutine wolfe(ndg,prod)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=2000, maxdg=1000)
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
      t0=1.0d+12
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
c===============================================================
      subroutine equations(n,z1)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=2000, maxdg=1000)
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

c=====================================================================
c Subroutine dgrad calculates discrete gradients
c=====================================================================
      subroutine dgrad(x,sl,g,dg,f1,f4,ndg,pwt)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=2000, maxdg=1000)
      double precision x1(maxvar),g(maxvar),x(maxvar),dg(maxvar)
      common /csize/m,/cicor/icor,/cns/ns

      a1=0.0d+00
      do k=1,m
       d1=dabs(g(k))
       if (a1.lt.d1) then
         a1=d1
         icor=k
       end if
      end do

      do k=1,m
       x1(k)=x(k)+sl*g(k)
      end do

      IF(ndg.gt.1) r2=f4
      IF(ndg.eq.1) call fv(x1,r2)

      if(ns.eq.1) call auxgrad(x1,dg,pwt)
      if(ns.eq.2) call funcgrad(x1,dg,pwt)
      
      dsum=0.0d+00
      do k=1,m
       if (k.ne.icor) then
        dsum=dsum+dg(k)*g(k)
       end if
      end do
      dsum=sl*dsum
      dg(icor)=(r2-f1-dsum)/(sl*g(icor))
      return
      end

c===========================================================
c Line search (Armijo-type), Step 5 Algorithm 2.
c===========================================================
      subroutine armijo(x,g,f1,f5,f4,sl,step,r)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=2000, maxdg=1000)
      common /csize/m
      double precision x(maxvar),g(maxvar),x1(maxvar)
      step=sl
      f5=f4
  1   step=2.0d+00*step
      do i=1,m
       x1(i)=x(i)+step*g(i)
      end do
      call fv(x1,f6)
      f3=f6-f1+5.0d-02*step*r
      IF(f3.gt.0.0d+00) then
       step=step/2.0d+00
       return
      END IF
      f5=f6
      GO TO 1
      return
      end
c========================================
      subroutine fv(x,f)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=2000)
      double precision x(maxvar)
      common /cns/ns
c================================================================
      if(ns.eq.1) call auxfunc(f,x)
      if(ns.eq.2) call func(f,x)
c=================================================================
      return
      end
