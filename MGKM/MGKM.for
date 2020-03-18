c=============================================================
c MGKM: modified global k-means algorithm
c
c The name of input file should be: datainput
c
c The name of output file is: MGKMResuslts
c
c with gkmeans0: to calculate the next cluster center
c without gkmeans0: to calculate distance matrix
c
c====================================================================
      PARAMETER(maxdim=2000, maxrec=200000, maxclust=100, maxnft=100
     1 ,maxclass=100, maxsize=7000000)
      implicit double precision(a-h,o-z)
      double precision x(maxdim),a(maxrec,maxnft),dminim(maxrec)
     1 ,x2(maxsize),amed(maxnft),x3(maxdim),xbest(maxdim)
     2 ,x5(maxsize),plabel(maxclass),dcent(maxclust,maxclust)
      integer nk(maxclust,maxrec),nel(maxclust),nob(maxclass)
     1 ,lcand(maxrec)
      common /c22/a,/anclust/nclust,/cnft/nft,/cmf/mf,/cnclass/nclass
     1 ,/crecord/nrecord,/cnc/nc,/cnk/nk,nel,/adminim/dminim
     2 ,/ccand/ncand,lcand,/ctnorm/tnorm,/cnob/plabel,nob
     3 ,/cnpnout/npurity,/cgamma/gamma1,gamma2,/cind/ind
     4 ,/cdcent/dcent,/cbest/xbest
      open(40,file='MGKMResuslts.txt')
      open(78,file='datainput.txt',status='old',form='formatted')     
c====================================================================
      PRINT *,'Number of features:'
      read *,nft
      PRINT *,' '
      PRINT *,'Class outputs?'
      PRINT *,' '
      PRINT *,'              1 - if yes'
      PRINT *,'              2 - if no'
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
      PRINT *,'Compute indices (DB, Dunn, Purity, Silhouettes)?'
      PRINT *,' '
      PRINT *,'   1 - yes'
      PRINT *,'   2 - no'
      read *,ind
      PRINT *,' ' 
c======================================================================
      tlimit=7.2d+04
      do i=1,maxrec
       read(78,*,END=901,ERR=900) (a(i,k),k=1,nft)
       nrecord=i
      end do
  900 stop 'Error in input file'       
  901 WRITE(40,*) 'Input complete. Number of records: ',nrecord
      WRITE(40,*)
c=======================================================================
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
c=====================================================================
      IF(npurity.eq.1) mf=nft-1
      IF(npurity.eq.2) mf=nft
      IF(npurity.eq.2) noutcom=0
      tnorm=0.0d+00
c=====================================================================
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
       plabel(1)=a(1,nft)
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
       if(nc.eq.1) then
                   call step1(f,x)
                   do i=1,mf
                    xbest(i)=x(i)
                   end do
                   toler=1.0d-02*f/dble(nrecord)
                   go to 1
       END if
       call step2(toler,nstart,x2)
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
         tnorm=tnorm+1.0d+00
         IF(f31.LE.(1.0d-01*toler)) THEN
          GO TO 1200
         END IF 
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
       fbest=1.0d+26
       do j=1,nstart
        do i=1,mf
         x(i+(nc-1)*mf)=x2(i+(j-1)*mf)
        END do
        do j1=1,m
         x3(j1)=x(j1)
        end do
        call gkmean1(f,x3)
        if (f.lt.fbest) then
         fbest=f
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
         tnorm=tnorm+1.0d+00
         dcent(j,i)=dcent(i,j)
        end do
       end do
       
       PRINT 123,nc
 123   FORMAT(i20)
c================================================================
  1    call clusters(x,f,db,dunn,purity,ns)
       if(npurity.eq.2) purity=0.0d+00
       call cpu_time(time3)
       time4=time3-time1
       write(40,603) nc,f,db,dunn,purity,ns,tnorm,time4
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
         toler3=9.0d-01*rad(k)
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
c====================================================
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
c============================================================
      return
      end

c=======================================================================
c  Step1 computes centroid and value of clustering function at centroid
c=======================================================================
      subroutine step1(f,x)
      PARAMETER(maxdim=2000, maxrec=200000, maxclust=100, maxnft=100)
      implicit double precision (a-h,o-z)
      double precision x(maxdim),a(maxrec,maxnft)
      common /c22/a,/cmf/mf,/crecord/nrecord,/ctnorm/tnorm
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
       tnorm=tnorm+1.0d+00
       do j=1,mf
        f1=f1+(a(i,j)-x(j))*(a(i,j)-x(j))
       END do
       f=f+f1
      END do
      return
      end

c=======================================================================
c  Step2 computes clusters for each data point
c=======================================================================
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
      subroutine gkmean0(x,fval)
      implicit double precision (a-h,o-z)
      PARAMETER(maxdim=2000, maxrec=200000, maxclust=100, maxnft=100)
      double precision x(maxdim),a(maxrec,maxnft),dminim(maxrec)
     1 ,xbest(maxdim),dmin(maxclust)
      integer lmembers(maxrec),l5(maxrec),list1(maxrec)
      common /c22/a,/crecord/nrecord,/cmf/mf,/adminim/dminim
     1 ,/ctnorm/tnorm,/cbest/xbest,/cnc/nc,/clist1/list1

      nchange=nrecord/1000         
      nmembers=0
  2   l6=0
      n1=0

      do i=1,nc-1
       tnorm=tnorm+1.0d+00
       dmin(i)=0.0d+00
       do k=1,mf
        dmin(i)=dmin(i)+(xbest(k+(i-1)*mf)-x(k))**2
       end do
      end do

      do i=1,nrecord
       f1=dminim(i)
       f10=4.0d+00*f1
       i1=list1(i)
       if(dmin(i1).lt.f10) then
        tnorm=tnorm+1.0d+00
        f3=0.0d+00
        do j=1,mf
         f3=f3+(a(i,j)-x(j))**2
        end do
        IF(f3.lt.f1) THEN
         do j=1,nmembers
          if(i.eq.lmembers(j)) go to 1
         end do
         l6=l6+1
 1       n1=n1+1
         l5(n1)=i
        END if
       end if 
      end do

      nmembers=n1
      do j=1,n1
        lmembers(j)=l5(j)
      end do

      do j=1,mf
       f4=0.0d+00
       do i=1,nmembers
        i1=lmembers(i)
        f4=f4+a(i1,j)
       end do
       if(n1.gt.0) x(j)=f4/dble(n1)
      end do

      IF(l6.gt.nchange) GO TO 2
      fval=0.0d+00
      do i=1,nrecord
       f3=0.0d+00
       tnorm=tnorm+1.0d+00
       do j=1,mf
        f3=f3+(a(i,j)-x(j))**2
       end do
       fval=fval+dmin1(dminim(i),f3)
      end do
      return
      end
c================================================================
      subroutine gkmean1(f,x)
      implicit double precision (a-h,o-z)
      PARAMETER(maxdim=2000, maxrec=200000, maxclust=100, maxnft=100)
      double precision x(maxdim),a(maxrec,maxnft),f4(maxnft),
     1 dminim(maxrec),dcent1(maxclust,maxclust)
      integer lmembers(maxclust,maxrec),nmembers(maxclust)
     1 ,l5(maxclust,maxrec),n1(maxclust)
      common /c22/a,/crecord/nrecord,/cnc/nc,/adminim/dminim
     1 ,/cmf/mf,/almembers/lmembers,nmembers,/ctnorm/tnorm

      nchange=nrecord/1000

      do k=1,nc
       nmembers(k)=0
      end do

  2   l6=0
      do k=1,nc
       n1(k)=0
      end do

      do i=1,nc
       dcent1(i,i)=0.0d+00
      end do

      do i=1,nc
       do j=i+1,nc
        dcent1(i,j)=0.0d+00
        do k=1,mf
         dcent1(i,j)=dcent1(i,j)+(x(k+(i-1)*mf)-x(k+(j-1)*mf))**2
        end do
        tnorm=tnorm+1.0d+00
        dcent1(j,i)=dcent1(i,j)
       end do
      end do

      do i=1,nrecord
       f1=1.0d+26
       do k=1,nc
        if(k.gt.1) then
         d0=dcent1(k,kmin)
         d1=4.0d+00*f1
         if(d0.ge.d1) go to 3
        end if
        tnorm=tnorm+1.d+00
        f3=0.d+00
        do j=1,mf
         f3=f3+(a(i,j)-x(j+(k-1)*mf))**2
        end do
        IF(f1.ge.f3) THEN
            f1=f3
            kmin=k
        END if
 3     end do
       do j=1,nmembers(kmin)
        if(i.eq.lmembers(kmin,j)) go to 1
       end do
       l6=l6+1
 1     n1(kmin)=n1(kmin)+1
       l5(kmin,n1(kmin))=i
      end do

      do i=1,nc
       nmembers(i)=n1(i)
       do j=1,n1(i)
        lmembers(i,j)=l5(i,j)
       end do
      end do

      do k=1,nc
       do j=1,mf
        f4(j)=0.0d+00
        do i=1,nmembers(k)
         i1=lmembers(k,i)
         f4(j)=f4(j)+a(i1,j)
        end do
       end do
       do j=1,mf
        if(nmembers(k).gt.0) x(j+(k-1)*mf)=f4(j)/dble(nmembers(k))
       end do
      end do

      IF(l6.gt.nchange) GO TO 2
      f=0.0d+00
      do k=1,nc
       do i=1,nmembers(k)
        tnorm=tnorm+1.0d+00
        i1=lmembers(k,i)
        f1=0.0d+00
        do j=1,mf
         f1=f1+(a(i1,j)-x(j+(k-1)*mf))**2
        end do
        f=f+f1
       end do
      end do
c=====================================================
      return
      end
