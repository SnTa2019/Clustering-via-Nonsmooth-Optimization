c============================================================================================
c Code for new Absolute indices: determining compactness, separability and number of clusters
c============================================================================================
c     main programm
      PARAMETER(maxvar=7000, maxrec=200000, maxclust=100, maxnft=300
     1 ,maxclass=100, maxsize=7000000)
      implicit double precision(a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),dminim(maxrec)
     1,x2(maxsize),amed(maxnft),plabel(maxclass),xbest(maxvar),z(maxvar)
     2 ,fval(maxrec),x5(maxvar),x3(maxvar),xglob(maxvar)
     3 ,dminim1(maxrec),dk1(maxrec),dk2(maxrec),dk3(maxrec)
      integer nk(maxclust,maxrec),nel(maxclust),nob(maxclass)
      common /csize/m,/c22/a,/anclust/nclust,/cnft/nft,/cmf/mf,/cnc/nc
     1 ,/crecord/nrecord,/cnk/nk,nel,/adminim/dminim,/cns/ns
     2 ,/cnob/plabel,nob,/cgamma/gamma1,gamma2,/ctnorm/tnorm,/cnel2/nel2
     3 ,/cnpnout/npurity,/cnclass/nclass,/craver/raver,/cglob/xglob
     4 ,/cepsil/epsil1,/cmed/dmed0,/cmean/dmean0,/cstd/std0,/cnnn/nmin
      open(49,file='Results.txt')
      open(78,file='inputdata.txt',status='old',form='formatted')

c======================================================================
c users to provide 
c======================================================================
      PRINT *,' '
      PRINT *,'Number of features:'
      read *,nft
      PRINT *,' '
      PRINT *,'Maximum number of clusters:'
      read *,nclust
      PRINT *,' '
      PRINT *,'Minimum number of clusters:'
      read *,nmin
      PRINT *,' '
c===========================================================
      tlimit=1.8d+04
      do i=1,maxrec
       read(78,*,END=901,ERR=900) (a(i,k),k=1,nft)
       nrecord=i
      end do
  900 stop 'Error in input file'       
  901 WRITE(49,*) 'Input complete. Number of records: ',nrecord
      WRITE(49,*)
      npurity =2
      IF(npurity.eq.1) mf=nft-1
      IF(npurity.eq.2) mf=nft
      IF(npurity.eq.2) noutcom=0
      tnorm=0.0d+00
c==============================================================
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
 22      nob(j)=nob(j)+1
        end do
        do i=1,nclass
         WRITE(49,41) i,plabel(i),nob(i)
        end do
       END if
 41    FORMAT('Class',i4,' label',f7.0,'  number of instances',i7) 

c======================================================================
c First method - L2, rad - radius of a(.,.), 
c dk2 - list of different distance, dk3 - list of increasing distances
c======================================================================
      do i=1,mf
       x(i)=0.0d+00
       do j=1,nrecord
        x(i)=x(i)+a(j,i)
       end do
       x(i)=x(i)/dble(nrecord)
       xglob(i)=x(i)
      end do
      rad=0.0d+00
      idif=0
      do i=1,nrecord
       dist=0.0d+00
       do j=1,mf
        dist=dist+(x(j)-a(i,j))**2
       end do
       d1=sqrt(dist)
       do j=1,idif
        if(d1.eq.dk2(j)) go to 201
       end do
       idif=idif+1
       dk2(idif)=d1
 201   rad=dmax1(rad,d1)
      end do
      rad=rad/dble(idif)
      epsil1=rad
    
c---------------------------------------------------
      call sortinc(dk2,dk3,idif)
      do i=1,idif
       dk3(i+1)=dk3(i)
      end do
      dk3(1)=0.0d+00
      nel3=idif+1
c=====================================================================
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
       gamma1=2.0d-01
       gamma2=1.0d-01
       gamma3=5.0d-02
      end if
      if(nrecord.gt.50000) then
       gamma1=5.0d-02
       gamma2=1.0d-02
       gamma3=1.0d-03
      end if
      call cpu_time(time1)
      
c====================================================================
      do nc=1,nclust
       if(nc.eq.1) then
                   call clusters(x,f)
                   raver=0.0d+00
                   nel(nc)=nrecord
                   do i=1,nrecord
                    dminim1(i)=dsqrt(dminim(i))
                    raver=raver+dminim1(i)
                    dk1(i)=dminim1(i)
                    nk(nc,i)=i
                   end do
                   raver=raver/dble(nrecord)
                   call compmedian(dk1,nrecord,dmed)
                   dmed0=dmed
                   call compstdev(dk1,dmean,std,nrecord)
                   dmean0=dmean
                   std0=std
                   call finresult(x)
                   write(49,*)
                   write(49,149) f
 149               format('The value of function f1 =',f30.4)             
                   toler=1.0d-02*f/dble(nrecord)

c====================================================================
       WRITE(49,*)
       WRITE(49,*)
       WRITE(49,704)
 704   FORMAT('Notations:')
       WRITE(49,*)
       WRITE(49,702)
 702   FORMAT('#Clust    - the number of clusters')  
       WRITE(49,705)
 705   FORMAT('Func. val - the optimal value of the clustering funct.')
       WRITE(49,706)
 706   FORMAT('Comp-f    - compactness of cluster distribution')
       WRITE(49,726)
 726   FORMAT('Comp-s    - compactness of reduced clust. distr.')
       WRITE(49,707)
 707   FORMAT('S-ratio-f - ratio of separated of clusters')
       WRITE(49,708)
 708   FORMAT('Margin-f  - margin for cluster distrbution')
       WRITE(49,727)
 727   FORMAT('Sep-err-f - fraction of correctly classified points')
       WRITE(49,737)
 737   FORMAT('S-ratio-s - ratio of separated of reduced clust. distr ')
       WRITE(49,738)
 738   FORMAT('Margin-s  - margin for reduced cluster distrbution')
       WRITE(49,739)
 739   FORMAT('Sep-er-s  - fraction of correctly class. points(reduced')
       WRITE(49,740)
 740   FORMAT('S-ratio-m - ratio of separated of clust. using medians')
       WRITE(49,741)
 741   FORMAT('Margin-m  - margin for cluster distrbution using median')
       WRITE(49,742)
 742   FORMAT('Sep-er-m  - fraction of correctly class.  with medians')
       WRITE(49,709)
 709   FORMAT('Sil-proc  - procent of points with positive silhouettes')
       WRITE(49,710)
 710   FORMAT('Sil-aver  - average value of silhouette')
       WRITE(49,711)
 711   FORMAT('DB        - Davies-Bouldin index')
       WRITE(49,712)
 712   FORMAT('XB        - Xie-Beni index')
       WRITE(49,715)
 715   FORMAT('DI        - Dunn index')
       WRITE(49,716)
 716   FORMAT('CH        - Calinski-Harabasz index')
       WRITE(49,717)
 717   FORMAT('Gstr      - Strict G-index')
       WRITE(49,718)
 718   FORMAT('Grex      - relaxed G-index')
       WRITE(49,*)
       WRITE(49,*)
c===================================================================
       WRITE(49,*)
       WRITE(49,701)
 701   FORMAT('#Clust','           Func. val.',' |','  Comp-f.'
     1 ,' Comp-s.',' Comp-m.','|',' S-rat-f','  M-min-f','   M-sum-f'
     1 ,' |',' S-rat-s','  M-min-s','   M-sum-s',' |',' S-rat-m'
     2 ,'  M-min-m','   M-sum-m',' |',' Sil-pr',' Sil-av','     DB'
     3 ,'     XB','      DI','          CH','        Gstr','     Grex')
       WRITE(49,*)
c====================================================================
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
       fbest=1.0d+30
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
c================================================================
       call clusters(x,f)
       if(npurity.eq.2) purity=0.0d+00
       call cpu_time(time3)
       if(nc.gt.1) call finresult(x)
       time4=time3-time1
       IF(time4.gt.tlimit) GO TO 2
  1   end do
  2   continue

      CLOSE(49)
      close(78)
      stop
      end

      subroutine clusters(x,f)
      PARAMETER(maxdim=2000, maxrec=200000, maxclust=100, maxnft=300,
     1 maxclass=100)
      implicit double precision (a-h,o-z)
      allocatable lcand1(:)
      double precision x(maxdim),a(maxrec,maxnft),dminim(maxrec)
     1 ,rad(maxclust),radmax(maxclust),ratio(maxclust)
     2 ,dist(maxclust,maxclust),dc(maxrec,maxclust),rad2(maxclust)
      integer nk(maxclust,maxrec),nel(maxclust),list1(maxrec)
     1 ,lcand(maxrec)
      common /c22/a,/cmf/mf,/anclust/nclust,/crecord/nrecord,/crad/rad
     1 ,/cnc/nc,/cnk/nk,nel,/cnclass/nclass,/cnpnout/npurity
     2 ,/ccand/ncand,lcand,/adminim/dminim,/ctnorm/tnorm,/clist1/list1
     3 ,/crma/radmax,/cdist/dist,/cdc/dc
       allocate(lcand1(maxrec))

      do j=1,nclust
       nel(j)=0
       rad(j)=0.0d+00
       radmax(j)=0.0d+00
      END do

      f=0.0d+00
      do k=1,nrecord
       f2=1.0d+30
       do j=1,nc
        f1=0.0d+00
        do k1=1,mf
         f1=f1+(a(k,k1)-x(k1+(j-1)*mf))**2
        END do
        dc(k,j)=f1
        if(f2.gt.f1) then
         f2=f1
         jmin=j
        end if
       END do
       dminim(k)=f2
       f=f+f2
       nel(jmin)=nel(jmin)+1
       nk(jmin,nel(jmin))=k
       list1(k)=jmin
       rad(jmin)=rad(jmin)+f2
       radmax(jmin)=dmax1(radmax(jmin),f2)
      end do
      tnorm=tnorm+dble(nc*nrecord)

      do i=1,nc
       dist(i,i)=0.0d+00
       do j=i+1,nc
        dist(i,j)=0.0d+00
        do k=1,mf
         dist(i,j)=dist(i,j)+(x(k+(i-1)*mf)-x(k+(j-1)*mf))**2
        end do
        dist(j,i)=dist(i,j)
       end do
      end do
 
      do k=1,nc
       if(nel(k).gt.0) then
          rad(k)=rad(k)/dble(nel(k))
          rad2(k)=rad(k)
        else 
          rad(k)=0.0d+00
          rad2(k)=0.0d+00
       end if
      end do

      if(nrecord.lt.500) then
       do k=1,nc
        rad2(k)=0.0d+00
       end do
      end if

      if((nc.gt.5).and.(nrecord.gt.500)) then
       ratmin=1.0d+26
       do k=1,nc
        if(rad(k).gt.0.0d+00) then
            ratio(k)=radmax(k)/rad(k)
           else
            ratio(k)=1.0d+00 
        end if
        ratmin=dmin1(ratmin,ratio(k))
       end do      
       do k=1,nc
        step1=5.0d-01*ratmin/ratio(k)
        rad2(k)=rad2(k)+step1*(radmax(k)-rad2(k))
       end do
      end if 

      if(nc.lt.nclust) then
       ncand=0
       do k=1,nc
        if(nel(k).gt.2) then
         toler3=5.0d-01*rad2(k)
         ncand1=0
         do i=1,nel(k)
          i1=nk(k,i)
          if(ncand1.eq.0) then
           if(dminim(i1).gt.rad2(k)) then
            ncand=ncand+1
            ncand1=ncand1+1
            lcand(ncand)=i1
            lcand1(ncand1)=i1           
           end if
          end if
          if(ncand1.gt.0) then
           if(dminim(i1).gt.rad2(k)) then
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
c======================================================================
      return
      end

c=========================================================================
c  Step2 computes clusters for each data point
c=========================================================================
      subroutine step2(toler,nstart,x2)
      PARAMETER(maxdim=2000, maxrec=200000, maxclust=100, maxnft=300
     1 ,maxsize=7000000)
      implicit double precision (a-h,o-z)
      allocatable fmin1(:),fval(:),l4(:),lcand1(:)
      double precision x2(maxsize),a(maxrec,maxnft),dminim(maxrec)
     1 ,x4(maxnft)
      integer lcand(maxrec)
      common /c22/a,/cmf/mf,/adminim/dminim,/crecord/nrecord
     1 ,/ctnorm/tnorm,/ccand/ncand,lcand,/cgamma/gamma1,gamma2
      allocate(fmin1(maxrec))
      allocate(fval(maxrec))
      allocate(l4(maxrec))
      allocate(lcand1(maxrec))

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

      subroutine auxfunc(f,x)
      PARAMETER(maxvar=7000, maxrec=200000, maxclust=100, maxnft=300)
      implicit double precision (a-h,o-z)
      allocatable dist1(:)
      double precision x(maxvar),a(maxrec,maxnft),dminim(maxrec)
      integer iact(maxrec)
      common /c22/a,/cmf/mf,/crecord/nrecord,/adminim/dminim,/cic/iact
     1 ,/ctnorm/tnorm
      allocate(dist1(maxrec)) 

       do i=1,nrecord
        dist1(i)=0.0d+00
        do j=1,mf
         dist1(i)=dist1(i)+(x(j)-a(i,j))**2
        end do
       end do
       tnorm=tnorm+dble(nrecord)
       f=0.0d+00
       do i=1,nrecord
        f5=dist1(i)
        f4=dminim(i)
        if(f5.le.f4) then
           f=f+f5
           iact(i)=1
          else
           f=f+f4
           iact(i)=2 
        end if
       end do
      return
      end
c=====================================================================

      subroutine auxgrad(x,grad)
      PARAMETER(maxvar=7000, maxrec=200000, maxclust=100, maxnft=300)
      implicit double precision (a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),dminim(maxrec)
     1 ,grad(maxvar)
      integer iact(maxrec)  
      common /c22/a,/cmf/mf,/crecord/nrecord,/adminim/dminim,/cic/iact
     
      do j=1,mf
       grad(j)=0.0d+00
      end do
      do i=1,nrecord
       if(iact(i).eq.1) then
        do j=1,mf
         grad(j)=grad(j)+2.0d+00*(x(j)-a(i,j))
        end do
       end if 
      end do
      return
      end

c=====================================================================
      subroutine func(f,x)
      PARAMETER(maxvar=7000, maxrec=200000, maxclust=100, maxnft=300)
      implicit double precision (a-h,o-z)
      allocatable dist2(:,:)
      double precision x(maxvar),a(maxrec,maxnft)
      integer lact3(maxrec)
      common /c22/a,/cmf/mf,/crecord/nrecord,/cnc/nc,/ctnorm/tnorm
     2 ,/cic1/lact3
      allocate(dist2(maxrec,maxclust))
      
      do i=1,nrecord
       do j=1,nc
        dist2(i,j)=0.0d+00
        do k=1,mf
         dist2(i,j)=dist2(i,j)+(x(k+(j-1)*mf)-a(i,k))**2
        end do
       end do
      end do
      tnorm=tnorm+dble(nrecord*nc)
   
      f=0.0d+00
      do i=1,nrecord
       f4=1.0d+30
       do k=1,nc
        f5=dist2(i,k)
        if(f4.gt.f5) then
         f4=f5
         lact3(i)=k
        end if                 
       end do
       f=f+f4
      end do
      return
      end
c=====================================================================

      subroutine funcgrad(x,grad)
      PARAMETER(maxvar=7000, maxrec=200000, maxclust=100, maxnft=300)
      implicit double precision (a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),grad(maxvar)
      integer lact3(maxrec)
      common /c22/a,/cmf/mf,/crecord/nrecord,/cnc/nc,/csize/m
     1 ,/cic1/lact3

      do i=1,m
       grad(i)=0.0d+00
      end do     

      do i=1,nrecord
       j1=lact3(i)
       do j=1,mf
        k1=(j1-1)*mf+j
        grad(k1)=grad(k1)+2.0d+00*(x(k1)-a(i,j))
       end do
      end do
      return
      end

c=====================================================================
      subroutine dgm(x,f2)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=7000, maxdg=1000,maxit=100000)
      allocatable fvalues(:),prod(:,:),w(:,:)
      double precision x(maxvar),x1(maxvar),g(maxvar),v(maxvar),z(maxdg)
      INTEGER ij(maxdg)
      common /csize/m,/cij/ij,jvertex,/cz/z,/ckmin/kmin,/cnc/nc,/cns/ns
      allocate(fvalues(maxit))
      allocate(prod(maxdg,maxdg))
      allocate(w(maxdg,maxvar))
c====================================================================
      dist1=1.0d-07
      step0=-5.0d-02
      div=1.0d-01
      eps0=1.0d-07
      slinit=1.0d+00
      slmin=1.0d-05*slinit
      sdif=1.0d-05
      mturn=4
      maxiter=5000
      niter=0
      nbundle=min(m+3,40)
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
            call dgrad(x,sl,g,v)
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
c  Subroutines Wolfe and Equations solves quadratic
c  programming problem, to find
c  descent direction, Step 3, Algorithm 2.
c===============================================================

      subroutine wolfe(ndg,prod)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=7000, maxdg=1000)
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
      t0=1.0d+30
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
c=====================================================================

      subroutine equations(n,z1)
      implicit double precision (a-h,o-z)
      allocatable b(:,:)
      PARAMETER(maxvar=7000, maxdg=1000)
      common /w01/a
      double precision a(maxdg,maxdg),z1(maxdg)
      allocate(b(maxdg,maxdg))
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
      subroutine dgrad(x,sl,g,dg)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=7000)
      double precision x1(maxvar),g(maxvar),x(maxvar),dg(maxvar)
      common /csize/m,/cns/ns
       do k=1,m
        x1(k)=x(k)+sl*g(k)
       end do
       if(ns.eq.1) call auxgrad(x1,dg)
       if(ns.eq.2) call funcgrad(x1,dg)
      return
      end

c===========================================================
c Line search (Armijo-type), Step 5 Algorithm 2.
c===========================================================
      subroutine armijo(x,g,f1,f5,f4,sl,step,r)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=7000)
      common /csize/m
      double precision x(maxvar),g(maxvar),x1(maxvar)
      step=sl
      f5=f4
  1   step=2.0d+00*step
      do i=1,m
       x1(i)=x(i)+step*g(i)
      end do
      call fv(x1,f6)
      f3=f6-f1+1.0d-02*step*r
      IF(f3.gt.0.0d+00) then
       step=step/2.0d+00
       return
      END IF
      f5=f6
      GO TO 1
      return
      end

      subroutine fv(x,f)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=7000)
      double precision x(maxvar)
      common /cns/ns
c================================================================
      if(ns.eq.1) call auxfunc(f,x)
      if(ns.eq.2) call func(f,x)
c=================================================================
      return
      end

c=================================================================
c  x      - vector of all cluster centers
c  a      - data set
c  rad    - (squared) average radius
c  dminim - (squared) distance between a given point and its center
c  radmax - radius of a cluster
c  dist   - (squared) distance between cluster centers
c  dist1  - distance between points and cluster centers
c  fk     - cluster error value (with real distance) for a cluster (DB)
c  sa     - sum of distances for a given cluster (for Silhouettes)
c  tmar   - margins between clusters
c  y      - auxiliary variable used to keep current cluster center (no need?)
c  comp   - compactness number of a cluster
c  xglob  - center of a data set
c  dminim1- distance between a given point and its center
c  dk1    - distance (with square root) between a point and its cluster center
c  dk2    - auxiliary vector used for distances
c  dk3    - ordered distances
c  cor    - core radii of clusters
c  ext    - external radii of clusters (cor+2*dispersion)
c  oirex  - relaxed absolute index of a cluster
c  oistr  - strict absolute index of a cluster
c  list   - cluster membership list
c  nel    - number of points in a given cluster
c  nk     - list of points in a given cluster
c  NIQ    - (i,j) number of points from cluster j neihboring cluster i
c  nk1    - list of points in an annulus
c  ln     - list of neibhoring clusters for a given cluster
c  nn     - number of neigbhors for a given cluster
c  IQ     - (i,j) list of points from cluster j neihboring cluster i
c  w      - degree of membership of point in fuzzy clustering
c  ind    - indicator showing whether annulus is empty or nonempty
c==================================================================
      subroutine finresult(x)
      PARAMETER(maxvar=7000, maxrec=200000, maxclust=100, maxnft=300
     1 ,maxclass=100, maxdiv=200000)
      implicit double precision(a-h,o-z)
      allocatable:: dist1(:,:), w(:,:),IQ(:,:,:),dk1(:),dk3(:)
     1 ,dist2(:,:),tmar(:,:)
      allocatable dminim1(:),lind(:)
      double precision x(maxvar),a(maxrec,maxnft),dminim(maxrec)
     1 ,radmax(maxclust),dist(maxclust,maxclust),fk(maxclust)
     2 ,sa(maxclust),xglob(maxvar),compact(3,maxclust)
     3 ,comp(maxclust),cor(maxclust),ext(maxclust),oirex(maxclust)
     4 ,oistr(maxclust),y(maxvar),compn(maxclust),radrex(maxclust)
     5 ,radmed(maxclust),compm(maxclust),sepbar(3,maxclust)
     6 ,sephat(3,maxclust),silave(maxclust),dab(maxclust),xben(maxclust)
     7 ,cen(maxclust),parc(maxclust),dunn(maxclust),cal(maxclust)
     8 ,str(maxclust),rex(maxclust)
      integer list(maxrec),nel(maxclust),nk(maxclust,maxrec)
     2 ,NIQ(maxclust,maxclust),nk1(maxrec),ln(maxclust,maxclust)
     3 ,nn(maxclust)
      common /c22/a,/cnft/nft,/cmf/mf,/cnclass/nclass,/cnc/nc
     1 ,/crecord/nrecord,/clist1/list,/cepsil/epsil1,/adminim/dminim
     2 ,/cnk/nk,nel,/cdist/dist,/cnk1/nk1,nel1,/craver/raver
     3 ,/cglob/xglob,/cnel2/nel2,/cmed/dmed0,/cmean/dmean0,/cstd/std0
     4 ,/clustnumber/compact,sepbar,sephat,silave,dab,xben,cen,parc
     5 ,dunn,cal,str,rex,/anclust/nclust
      allocate(dk1(maxrec))
      allocate(dk3(maxrec))
      allocate(dminim1(maxrec))
      allocate(dist1(maxrec,maxclust))
      allocate(w(maxrec,maxclust))
      allocate(IQ(maxclust,maxclust,maxrec))
      allocate(lind(maxrec))
      allocate(dist2(maxclust,maxclust))
      allocate(tmar(maxclust,maxclust))
c==============================================================
c Function value
c==============================================================
      f=0.0d+00
      do i=1,nrecord
       f = f+dminim(i)
       dminim1(i)=sqrt(dminim(i))
      end do
      do k=1,nc
       k4=0
       radmax(k)=0.0d+00
       do j=1,nel(k)
        k2=nk(k,j)
        radmax(k)=dmax1(radmax(k),dminim1(k2))                          ! radius of cluster 
        k4=k4+1
        dk1(k4)=dminim1(k2)
       end do
       call compstdev(dk1,dmean,std,nel(k))
       radrex(k)=dmean+2.0d+00*std                                      ! radius mean+2*standard deviation
       radrex(k)=dmin1(radrex(k),radmax(k))                             
       call compmedian(dk1,nel(k),dmed1)
       radmed(k)=dmed1                                                  ! radius using median
      end do
      do i=1,nc
        dist(i,i)=0.0d+00
        do j=i+1,nc
         dist(i,j)=0.0d+00
         do k=1,mf
          dist(i,j)=dist(i,j)+(x(k+(i-1)*mf)-x(k+(j-1)*mf))**2          ! dist(i,j) - distance between clusters i and j
         end do
         dist(i,j)=dsqrt(dist(i,j))
         dist(j,i)=dist(i,j)
        end do
      end do

      do i=1,nrecord
        do j=1,nc
         dist1(i,j)=0.0d+00
         do k=1,mf
          dist1(i,j)=dist1(i,j)+(x(k+(j-1)*mf)-a(i,k))**2               ! dist1(i,j) - distance between point i and cluster j
         end do
         dist1(i,j)=dsqrt(dist1(i,j))
        end do
      end do
c=====================================================================
c  Calculating compactness of clusters (full version)
c=====================================================================
      eps=1.0d+00*epsil1
      if(eps.eq.0.0d+00) eps=1.0d-02
      
      compf=0.0d+00
      compmean=0.0d+00
      compmed=0.0d+00

      do k=1,nc
       do i=1,mf
        y(i)=x(i+(k-1)*mf)
       end do
       call intervals(k,nstep,dminim1,dk3)

c      full cluster
       comp2=0.0d+00
       comp3=0.0d+00
       do i=2,nstep
        dk4=dk3(i)-dk3(i-1)
        if(dk4.le.eps) then
         dk5=dk3(i)+eps
         dk6=dk3(i)-eps
         nel1=0
         do k1=1,nel(k)
          k2=nk(k,k1)
          if((dminim1(k2).gt.dk6).and.(dminim1(k2).le.dk5)) then
             nel1=nel1+1
             nk1(nel1)=k2
          end if
         end do
         if(nel1.gt.0) then
           call sphere(y,c1)
           comp2=comp2+(1.0d+00-c1)*dk4
         end if  
        end if
        if(dk4.gt.eps) then
          comp3=comp3+dk4-eps         
        end if
       end do
       if(radmax(k).gt.0.0d+00) then
            comp(k)=1.0d+00-(comp2+comp3)/radmax(k)                     ! compactness of a given cluster k
         else
            comp(k)=1.0d+00
       end if

       comp2=0.0d+00
       comp3=0.0d+00
       do i=2,nstep
        if(dk3(i).le.radrex(k)) then
         d8=dk3(i)
         dk4=dk3(i)-dk3(i-1)
         if(dk4.le.eps) then
          dk5=dk3(i)+eps
          dk6=dk3(i)-eps
          nel1=0
          do k1=1,nel(k)
           k2=nk(k,k1)
           if((dminim1(k2).gt.dk6).and.(dminim1(k2).le.dk5)) then
              nel1=nel1+1
              nk1(nel1)=k2
           end if
          end do
          if(nel1.gt.0) then
            call sphere(y,c1)
            comp2=comp2+(1.0d+00-c1)*dk4
          end if  
         end if
         if(dk4.gt.eps) then
           comp3=comp3+dk4-eps         
         end if
        end if
       end do
       if(d8.gt.0.0d+00) then
            compn(k)=1.0d+00-(comp2+comp3)/d8
         else
            compn(k)=1.0d+00
       end if
       compn(k)=dmin1(1.0d+00,compn(k))
       comp2=0.0d+00
       comp3=0.0d+00
       do i=2,nstep
        if(dk3(i).le.radmed(k)) then
         d8=dk3(i)
         dk4=dk3(i)-dk3(i-1)
         if(dk4.le.eps) then
          dk5=dk3(i)+eps
          dk6=dk3(i)-eps
          nel1=0
          do k1=1,nel(k)
           k2=nk(k,k1)
           if((dminim1(k2).gt.dk6).and.(dminim1(k2).le.dk5)) then
              nel1=nel1+1
              nk1(nel1)=k2
           end if
          end do
          if(nel1.gt.0) then
            call sphere(y,c1)
            comp2=comp2+(1.0d+00-c1)*dk4
          end if  
         end if
         if(dk4.gt.eps) then
           comp3=comp3+dmax1(0.0d+00,dk4-eps)         
         end if
        end if
       end do
       if(d8.gt.0.0d+00) then
            compm(k)=1.0d+00-(comp2+comp3)/d8
         else
            compm(k)=1.0d+00
       end if
       compm(k)=dmin1(1.0d+00,compm(k))
       
       compf=compf+dble(nel(k))*comp(k)/dble(nrecord)
       compmean=compmean+dble(nel(k))*compn(k)/dble(nrecord)            ! total for k-clustering
       compmed=compmed+dble(nel(k))*compm(k)/dble(nrecord)
      end do
      compact(1,nc)=compf
      compact(2,nc)=compmean
      compact(3,nc)=compmed
      
      if(nc.eq.1) then
       write(49,*)
       write(49,249) nc,f,compf,compmean,compmed
249    format(i8,f24.6,3f12.4)
       return       
      end if
c======================================================================= 
c Calculation of sets Q1 and Q2 - full version
c n(i,j) is the number of points from i-th cluster close to j-th cluster
c iq(i,j,k) is their list 
c=======================================================================
      if(nc.gt.1) then
       do i=1,nc
        do j=1,nc
         niq(i,j)=0
        end do
       end do

       do i=1,nc
        k1=nel(i)
        do j=i+1,nc
         k3=nel(j)
         do k=1,k1
          k2=nk(i,k)
          if(dist1(k2,j).le.dist(i,j)) then
           niq(i,j)=niq(i,j)+1
           iq(i,j,niq(i,j))=k2
          end if
         end do
         do k=1,k3
           k2=nk(j,k)
           if(dist1(k2,i).le.dist(j,i)) then
            niq(j,i)=niq(j,i)+1
            iq(j,i,niq(j,i))=k2
           end if
         end do
        end do
       end do
      end if
c=======================================================================
c Calculation of margins - full version-tmar(i,j)-margin between i and j
c=======================================================================
      if(nc.gt.1) then
       do i=1,nc
        do j=1,nc
         tmar(i,j) = 0.0d+00 
        end do
       end do

       do i=1,nc
        do j=i+1,nc
         k1=niq(i,j)
         if(k1.gt.0) then
          r1=0.0d+00
          do j1=1,k1
           k2=iq(i,j,j1)
           r2=dist1(k2,i)
           r1=dmax1(r1,r2)
          end do
          dist2(i,j)=r1
          k3=niq(j,i)
          r3=0.0d+00
          do j1=1,k3
           k4=iq(j,i,j1)
           r2=dist1(k4,j)
           r3=dmax1(r3,r2)
          end do
          dist2(j,i)=r3
          tmar(i,j)=(dist(i,j)-r1-r3)/dist(i,j)
         end if 
         if(k1.eq.0) then
          tmar(i,j)=0.0d+00
         end if
         tmar(j,i)=tmar(i,j)
        end do
       end do
      end if
c-------------------------------------------------------------------------
c identification of neighbors-full version: ln(i,j)-list of neighbors of i
c-------------------------------------------------------------------------
      if(nc.eq.2) then
       nn(1)=1
       nn(2)=1
       ln(1,nn(1))=2
       ln(2,nn(2))=1
      end if
      if(nc.gt.2) then
       do i=1,nc
        nn(i)=0
       end do
       do i=1,nc
        do j=1,nc
         if(j.ne.i) then
          k2=0
          do k=1,nc 
           if((k.ne.i).and.(k.ne.j)) then
            sp=0.0d+00
            do k1=1,mf
             sp=sp+(x(k1+(j-1)*mf)-x(k1+(i-1)*mf))*(x(k1+(k-1)*mf)
     1          -x(k1+(i-1)*mf))
            end do
            sp=sp/(dist(i,j)*dist(i,k))
            if((sp.gt.0.5).and.(dist(i,k).lt.dist(i,j))) k2=k2+1
           end if
          end do
          if(k2.eq.0) then
           nn(i)=nn(i)+1
           ln(i,nn(i))=j
          end if
         end if
        end do
       end do
      end if
c-----------------------------------------------------
c Calculation of total margins
c-----------------------------------------------------
      if(nc.gt.1) then
       skbar=0.0d+00
       skhat=0.0d+00
       n1=0
       n2=0
       do i=1,nc
        sk=1.0d+30
        sk1=0.0d+00
        do j2=1,nn(i)
         n1=n1+1
         j=ln(i,j2)
         if(j.ne.i) then
          if(tmar(i,j).gt.0.0d+00) n2=n2+1
          sk=dmin1(sk,tmar(i,j))
          sk1=sk1+tmar(i,j)
         end if
        end do
        skbar=skbar+sk
        skhat=skhat+sk1
       end do
       sep1=dble(n2)/dble(n1)
       skbar=skbar/dble(nc)
       skhat=skhat/dble(n1)
      end if
      if(nc.eq.1) then
       sepbar(1,1)=0.0d+00
       sephat(1,1)=0.0d+00
       skbar1=0.0d+00
       skhat1=0.0d+00
      end if
      if(nc.gt.1) then
       sepbar(1,nc)=skbar
       sephat(1,nc)=skhat
       skbar1=skbar
       skhat1=skhat
      end if
c----------------------------------------------------
c calculation of misclassified points - full version
c----------------------------------------------------
      if(nc.gt.1) then
       do i=1,nrecord
        lind(i)=0
       end do
       do i=1,nc
        do j=1,nc
         if(j.ne.i) then
          n1=niq(i,j)
          do k=1,n1
           k1=iq(i,j,k)
           if(dist1(k1,j).le.dist2(j,i)) then
            lind(k1)=1
           end if
          end do
         end if
        end do
       end do
       mis=0
       do i=1,nrecord
        mis=mis+lind(i)
       end do
       cmis1=1.0d+00-dble(mis)/dble(nrecord) 
      end if
c====================================================================      
c  Calculation of sets Q1 and Q2 - relaxed version
c============================================================
      if(nc.gt.1) then
       do i=1,nc
        do j=1,nc
         niq(i,j)=0
        end do
       end do

       do i=1,nc
        k1=nel(i)
        do j=i+1,nc
         k3=nel(j)
         do k=1,k1
          k2=nk(i,k)
          if(dist1(k2,i).le.radrex(i)) then
           if(dist1(k2,j).le.dist(i,j)) then
            niq(i,j)=niq(i,j)+1
            iq(i,j,niq(i,j))=k2
           end if
          end if
         end do
         do k=1,k3
           k2=nk(j,k)
           if(dist1(k2,j).le.radrex(j)) then
            if(dist1(k2,i).le.dist(j,i)) then
             niq(j,i)=niq(j,i)+1
             iq(j,i,niq(j,i))=k2
            end if
           end if
         end do
        end do
       end do
      end if
c=============================================================
c Calculation of margins - relaxed version
c=============================================================
      if(nc.gt.1) then
       do i=1,nc
        do j=1,nc
         tmar(i,j) = 0.0d+00 
        end do
       end do

       do i=1,nc
        do j=i+1,nc
         k1=niq(i,j)
         if(k1.gt.0) then
          r1=0.0d+00
          do j1=1,k1
           k2=iq(i,j,j1)
           r2=dist1(k2,i)
           r1=dmax1(r1,r2)
          end do
          dist2(i,j)=r1
          k3=niq(j,i)
          r3=0.0d+00
          do j1=1,k3
           k4=iq(j,i,j1)
           r2=dist1(k4,j)
           r3=dmax1(r3,r2)
          end do
          dist2(j,i)=r3
          tmar(i,j)=(dist(i,j)-r1-r3)/dist(i,j)
         end if 
         if(k1.eq.0) then
          tmar(i,j)=0.0d+00
         end if
         tmar(j,i)=tmar(i,j)
        end do
       end do
      end if
c----------------------------------------------------
      if(nc.gt.1) then
       skbar=0.0d+00
       skhat=0.0d+00
       n1=0
       n2=0
       do i=1,nc
        sk=1.0d+30
        sk1=0.0d+00
        do j2=1,nn(i)
         n1=n1+1
         j=ln(i,j2)
         if(j.ne.i) then
          if(tmar(i,j).gt.0.0d+00) n2=n2+1
          sk=dmin1(sk,tmar(i,j))
          sk1=sk1+tmar(i,j)
         end if
        end do
        skbar=skbar+sk
        skhat=skhat+sk1
       end do
       seps1=dble(n2)/dble(n1)
       skbar=skbar/dble(nc)
       skhat=skhat/dble(n1)
      end if

      if(nc.eq.1) then
       sepbar(2,1)=0.0d+00
       sephat(2,1)=0.0d+00
       skbar2=0.0d+00
       skhat2=0.0d+00
      end if
      if(nc.gt.1) then
       sepbar(2,nc)=skbar
       sephat(2,nc)=skhat
       skbar2=skbar
       skhat2=skhat
      end if
c--------------------------------------------------------
c calculation of misclassified points - relaxed version
c--------------------------------------------------------
      if(nc.gt.1) then
       do i=1,nrecord
        lind(i)=0
       end do
       do i=1,nc
        do j=1,nc
         if(j.ne.i) then
          n1=niq(i,j)
          do k=1,n1
           k1=iq(i,j,k)
           if(dist1(k1,j).le.dist2(j,i)) then
            lind(k1)=1
           end if
          end do
         end if
        end do
       end do
       mis=0
       do i=1,nrecord
        mis=mis+lind(i)
       end do
       cmis2=1.0d+00-dble(mis)/dble(nrecord) 
      end if
c====================================================================      
c  Calculation of sets Q1 and Q2 - core version
c============================================================
      if(nc.gt.1) then
       do i=1,nc
        do j=1,nc
         niq(i,j)=0
        end do
       end do

       do i=1,nc
        k1=nel(i)
        do j=i+1,nc
         k3=nel(j)
         do k=1,k1
          k2=nk(i,k)
          if(dist1(k2,i).le.radmed(i)) then
           if(dist1(k2,j).le.dist(i,j)) then
            niq(i,j)=niq(i,j)+1
            iq(i,j,niq(i,j))=k2
           end if
          end if
         end do
         do k=1,k3
           k2=nk(j,k)
           if(dist1(k2,j).le.radmed(j)) then
            if(dist1(k2,i).le.dist(j,i)) then
             niq(j,i)=niq(j,i)+1
             iq(j,i,niq(j,i))=k2
            end if
           end if
         end do
        end do
       end do
      end if
c=============================================================
c Calculation of margins - core version
c=============================================================
      if(nc.gt.1) then
       do i=1,nc
        do j=1,nc
         tmar(i,j) = 0.0d+00 
        end do
       end do

       do i=1,nc
        do j=i+1,nc
         k1=niq(i,j)
         if(k1.gt.0) then
          r1=0.0d+00
          do j1=1,k1
           k2=iq(i,j,j1)
           r2=dist1(k2,i)
           r1=dmax1(r1,r2)
          end do
          dist2(i,j)=r1
          k3=niq(j,i)
          r3=0.0d+00
          do j1=1,k3
           k4=iq(j,i,j1)
           r2=dist1(k4,j)
           r3=dmax1(r3,r2)
          end do
          dist2(j,i)=r3
          tmar(i,j)=(dist(i,j)-r1-r3)/dist(i,j)
         end if 
         if(k1.eq.0) then
          tmar(i,j)=0.0d+00
         end if
         tmar(j,i)=tmar(i,j)
        end do
       end do
      end if
c----------------------------------------------------
      if(nc.gt.1) then
       skbar=0.0d+00
       skhat=0.0d+00
       n1=0
       n2=0
       do i=1,nc
        sk=1.0d+30
        sk1=0.0d+00
        do j2=1,nn(i)
         j=ln(i,j2)
         n1=n1+1
         if(j.ne.i) then
          if(tmar(i,j).gt.0.0d+00) n2=n2+1
          sk=dmin1(sk,tmar(i,j))
          sk1=sk1+tmar(i,j)
         end if
        end do
        skbar=skbar+sk
        skhat=skhat+sk1
       end do
       sepm1=dble(n2)/dble(n1)
       skbar=skbar/dble(nc)
       skhat=skhat/dble(n1)
      end if

      if(nc.eq.1) then
       sepbar(3,1)=0.0d+00
       sephat(3,1)=0.0d+00
       skbar3=0.0d+00
       skhat3=0.0d+00
      end if
      if(nc.gt.1) then
       sepbar(3,nc)=skbar
       sephat(3,nc)=skhat
       skbar3=skbar
       skhat3=skhat
      end if
c--------------------------------------------------------
c calculation of misclassified points - core version
c--------------------------------------------------------
      if(nc.gt.1) then
       do i=1,nrecord
        lind(i)=0
       end do
       do i=1,nc
        do j=1,nc
         if(j.ne.i) then
          n1=niq(i,j)
          do k=1,n1
           k1=iq(i,j,k)
           if(dist1(k1,j).le.dist2(j,i)) then
            lind(k1)=1
           end if
          end do
         end if
        end do
       end do
       mis=0
       do i=1,nrecord
        mis=mis+lind(i)
       end do
       cmis3=1.0d+00-dble(mis)/dble(nrecord) 
      end if
c======================================================================
c   Silhouette coeffcients: ratio - #of points with positive sil to total
c    #of points, sav - average of all silhouttes
c======================================================================
      if(nc.gt.1) then
       ns=0
       sav=0.0d+00
       do i=1,nrecord
         k1=list(i)
         do k=1,nc
           sa(k)=0.0d+00
         end do
         do j=1,nrecord
           d10=0.0d+00
           do j2=1,mf
            d10=d10+(a(i,j2)-a(j,j2))**2
           end do
           d10=dsqrt(d10)
           k2=list(j)
           sa(k2)=sa(k2)+d10
         end do
         if(nel(k1).gt.1) then
            da1=sa(k1)/dble(nel(k1)-1)
           else
            da1=0.0d+00 
         end if
         da2=1.0d+26
         do k=1,nc
          if(k.ne.k1) then
           if(nel(k).gt.0) then
             da3=sa(k)/dble(nel(k))
            else
             da3=1.0d+25 
           end if
           da2=dmin1(da2,da3)
          end if
         end do 
         sil=(da2-da1)/dmax1(da1,da2)
         if(sil.gt.0.0d+00) ns=ns+1
         sav=sav+sil
       end do
       ratio=dble(ns)/dble(nrecord)
       sav=sav/dble(nrecord)
      end if
      if(nc.eq.1) silave(1)=0.0d+00
      if(nc.gt.1) silave(nc)=sav
c====================================================================
c Davies-Bouldin (DB) validity index
c====================================================================
      if(nc.gt.1) then
       do i=1,nc 
        fk(i)=0.0d+00
       end do

       do i=1,nrecord
        k=list(i)
        fk(k)=fk(k)+dminim1(i)
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
          if(dist(i,k).gt.0.0d+00) then
              f2=fk2/dist(i,k)
             else
              f2=0.0d+00 
          end if  
          fm=dmax1(fm,f2)
         end if
        end do
        fdb=fdb+fm
       end do
       db=fdb/DBLE(nc)
      end if
      if(nc.eq.1) dab(1)=1.0d+03
      if(nc.gt.1) dab(nc)=db      
c=====================================================================
c Fuzzy membership
c=====================================================================
      if(nc.gt.1) then
       mc1=2
       mc=ceiling(dble(2)/dble(mc1-1))
       do i=1,nrecord
        do k=1,nc
          if(dist1(i,k).le.1.0d-04) then
           w(i,k)=1.0d+00
           do j=1,nc
            if(j.ne.k) w(i,j)=0.0d+00
           end do
           go to 111
          end if
        end do
        do k=1,nc
         w1=0.0d+00
         do k1=1,nc
          w1=w1+(dist1(i,k)/dist1(i,k1))**mc
         end do
         w(i,k)=1.0d+00/w1
        end do
111    end do
      end if
      if(nc.eq.1) then
       do i=1,nrecord
        w(i,1)=1.0d+00
       end do
      end if
c=====================================================================
c XB - Xie and Beni index
c=====================================================================
      if(nc.gt.1) then
       xb1=0.0d+00
       do k=1,nc
        do i=1,nrecord
         xb1=xb1+(w(i,k)*dist1(i,k))**2
        end do
       end do
       xb2=1.0d+30
       do i=1,nc
        do k=i+1,nc
         xb2=dmin1(xb2,dist(i,k)**2)
        end do
       end do
       xb2=dble(nrecord)*xb2
       xb=xb1/xb2
      end if
      if(nc.eq.1) xb=1.0d+00
      if(nc.eq.1) xben(1)=1.0d+00
      if(nc.gt.1) xben(nc)=xb   
c=====================================================================
c CE - Classification entropy index
c=====================================================================
      if(nc.gt.1) then
       ce=0.0d+00
       do k=1,nc
        do i=1,nrecord
         if(w(i,k).gt.1.0d-06) then
           ce=ce+w(i,k)*log(w(i,k))
         end if
        end do
       end do
       ce=-ce/dble(nrecord)
      end if
      if(nc.eq.1) ce=1.0d+00
      if(nc.eq.1) cen(1)=1.0d+00
      if(nc.gt.1) cen(nc)=ce 
c=====================================================================
c PC - Partition coefficient
c=====================================================================
      if(nc.gt.1) then
       pc=0.0d+00
       do k=1,nc
        do i=1,nrecord
         pc=pc+w(i,k)**2
        end do
       end do
       pc=pc/dble(nrecord)
      end if
      if(nc.eq.1) pc=1.0d+00
      if(nc.eq.1) parc(1)=1.0d+00
      if(nc.gt.1) parc(nc)=pc   
c=====================================================================
c DI - Dunn index
c=====================================================================
      if(nc.gt.1) then
       dt=0.0d+00
       do i=1,nc
        dt1=0.0d+00
        k1=nel(i)
        do k2=1,k1
         k3=nk(i,k2)
         do k4=k2+1,k1
          k5=nk(i,k4)
          dt2=0.0d+00
          do j=1,mf
           dt2=dt2+(a(k3,j)-a(k5,j))**2
          end do
          dt2=dsqrt(dt2)
          dt1=dmax1(dt1,dt2)
         end do
        end do       
        dt=dmax1(dt,dt1)       
       end do
      
       dt3=1.0d+30
       do i=1,nc
        i1=nel(i)
        do i2=i+1,nc
         i3=nel(i2)
         do i4=1,i1
          i5=nk(i,i4)
          do i6=1,i3
           i7=nk(i2,i6)
           dt4=0.0d+00
           do j=1,mf
            dt4=dt4+(a(i5,j)-a(i7,j))**2
           end do
           dt4=dsqrt(dt4)
           dt3=dmin1(dt3,dt4)
          end do
         end do
        end do
       end do
       di=dt3/dt
      end if 
      if(nc.eq.1) di=1.0d+00
      if(nc.eq.1) dunn(1)=1.0d+00
      if(nc.gt.1) dunn(nc)=di   
c=====================================================================
c CH - Calinski- Harabasz index
c=====================================================================
      if(nc.gt.1) then
       ds1=0.0d+00
       do i=1,nc
        ds2=0.0d+00
        do j=1,mf
         ds2=ds2+(x(j+(i-1)*mf)-xglob(j))**2
        end do
        ds2=dsqrt(ds2)
        ds1=ds1+dble(nel(i))*ds2
       end do
       ds3=0.0d+00
       do i=1,nrecord
        ds3=ds3+dminim1(i)
       end do
       ds4=dble(nrecord-nc)/dble(nc-1)
       CH=ds4*ds1/ds3
      end if
      if(nc.eq.1) CH=0.0d+00
      if(nc.eq.1) cal(1)=0.0d+00
      if(nc.gt.1) cal(nc)=ch   
c======================================================================
c Calculation G-indices
c======================================================================
      if(nc.gt.1) then
       Cor0=dmed0
       Ext0=dmean0+2.0d+00*std0
       do i=1,nc
        i1=nel(i)
        do j=1,i1
         i2=nk(i,j)
         dk1(j)=dminim1(i2)
        end do
        call compmedian(dk1,i1,dmed)
        call compstdev(dk1,dmean,std,i1)
        cor(i)=dmed
        ext(i)=dmean+2.0d+00*std
       end do
       do i=1,nc
        oirex(i)=1.0d+20
        oistr(i)=1.0d+20
        do j=1,nc
         if(j.ne.i) then
          oirex(i)=dmin1(oirex(i),dist(i,j)-cor(i)-cor(j))
          oistr(i)=dmin1(oistr(i),dist(i,j)-ext(i)-ext(j))
         end if
        end do
       end do
       gstr1=0.0d+00
       grex1=0.0d+00
       gstr2=0.0d+00
       grex2=0.0d+00
       do k=1,nc
        gstr1=gstr1+oistr(k)*dble(nel(k))
        gstr2=gstr2+ext(k)*dble(nel(k))
        grex1=grex1+oirex(k)*dble(nel(k))
        grex2=grex2+cor(k)*dble(nel(k))
       end do
       gstr=gstr1/gstr2
       grex=grex1/grex2
      end if
      if(nc.eq.1) then
       gstr=1.0d+00
       grex=1.0d+00
       str(1)=1.0d+00
       rex(1)=1.0d+00
      end if
      if(nc.gt.1) then
       str(nc)=gstr
       rex(nc)=grex
      end if
      
      if(nc.gt.1) then
       write(49,248) nc,f,compf,compmean,compmed,sep1,skbar1,skhat1
     1  ,seps1,skbar2,skhat2,sepm1,skbar3,skhat3,ratio,sav,db,xb,di,ch
     2  ,gstr,grex
      end if
248   format(i5,f22.4,' |',3f8.4,' |',f6.3,2f11.3,'|',f6.3,2f11.3,'|'
     1 ,f6.3,2f11.3,'|',2f7.3,f8.3,f8.3,f10.6,f12.3,2f9.3)
c====================================================================
      if(nc.eq.nclust) call cnumber
      return
      end 

c==================================================================
      subroutine sphere(x,coef)
      PARAMETER(maxvar=7000, maxrec=200000, maxnft=300)
      implicit double precision(a-h,o-z)
      allocatable anorm(:)
      double precision x(maxnft),a(maxrec,maxnft),d(maxvar)
      integer nk1(maxrec)
      common /c22/a,/cmf/mf,/cnk1/nk1,nel1
      allocate(anorm(maxrec))
c=====================================================================
      do i=1,nel1
       k=nk1(i)
       d1=0.0d+00
       do j=1,mf
        d1=d1+(a(k,j)-x(j))**2       
       end do
       d1=sqrt(d1)
       anorm(i)=d1
      end do

      md=2*mf
      nd=0
      do i=1,md
       do j=1,mf
        d(j)=0.0d+00
       end do
       if(i.le.mf) then
        d(i)=1.0d+00
        k=i
       end if
       if(i.gt.mf) then
        k=i-mf
        d(k)=-1.0d+00
       end if
       do j1=1,nel1
        j=nk1(j1)
        pr=(a(j,k)-x(k))*d(k)
        pr=pr/anorm(j1)
        if(pr.ge.2.5d-01) then
          nd=nd+1
          go to 111
        end if  
       end do
 111  end do
      coef=dble(nd)/dble(md)
c============================================================  
      return
      end 

c=========================================================================
c Finding the median of a set
c=========================================================================
      subroutine compmedian(dk,nel,dmed)
       PARAMETER(maxrec=200000)
       implicit double precision (a-h,o-z)
       allocatable b(:), b1(:)
       double precision dk(maxrec)
       allocate(b(maxrec))
       allocate(b1(maxrec))
       do i=1,nel
        b(i)=dk(i)
       end do
       call sortinc(b,b1,nel)
       n2=nel/2
       n3=nel-2*n2
       if(n3.eq.0) then
          dmed=(b1(n2)+b1(n2+1))/2.0d+00
         else
          dmed=b1(n2+1)
       end if
      end subroutine

c=======================================================================
c Finding the mean value and standard deviation
c=======================================================================
      subroutine compstdev(dk,dmean,std,nel)
       PARAMETER(maxrec=200000)
       implicit double precision (a-h,o-z)
       double precision dk(maxrec)
       dmean=0.0d+00 
       do i=1,nel
        dmean=dmean+dk(i)
       end do 
       dmean=dmean/dble(nel)
       std=0.0d+00
       do i=1,nel
        std=std+(dk(i)-dmean)**2
       end do
       if(nel.gt.1) then
            std=std/dble(nel-1)
         else
            std1=0.0d+00   
       end if
       std=dsqrt(std)
      end subroutine

c=========================================================================
c Finding intervals of a cluster
c=========================================================================
      subroutine intervals(nc,nstep,dminim1,dk3)
       PARAMETER(maxrec=200000, maxclust=100)
       implicit double precision (a-h,o-z)
       allocatable dk2(:)
       double precision dminim1(maxrec),dk3(maxrec)
       integer nel(maxclust),nk(maxclust,maxrec)
       common /cnk/nk,nel
       allocate(dk2(maxrec))
        i1=0
        n1=nel(nc)
        do i=1,n1
         i2=nk(nc,i)
         d1=dminim1(i2)
         do j=1,i1
          if(d1.eq.dk2(j)) go to 201
         end do
         i1=i1+1
         dk2(i1)=d1
 201    end do
c---------------------------------------------------
        call sortinc(dk2,dk3,i1)
        do i=1,i1
         dk3(i+1)=dk2(i)
        end do
        nstep=i1+1
        dk3(1)=0.0d+00
      end subroutine

c=======================================================================
c Ordering numbers in increasing order
c=======================================================================
      subroutine sortinc(b,b1,nel)
       PARAMETER(maxrec=200000)
       implicit double precision (a-h,o-z)
       double precision b(maxrec),b1(maxrec)
        do i=1,nel
         bmin=b(i)
         bmin1=b(i)
         kmin=i
         do k=i+1,nel
          if(b(k).lt.bmin) then
           kmin=k
           bmin=b(k)
          end if
         end do
         b(i)=b(kmin)
         b(kmin)=bmin1       
        end do
        do i=1,nel
         b1(i)=b(i)
        end do 
      end subroutine

c=======================================================================
c Determining the number of clusters
c=======================================================================
      subroutine cnumber
       PARAMETER(maxclust=100)
       implicit double precision (a-h,o-z)
       double precision compact(3,maxclust),sepbar(3,maxclust)
     2 ,sephat(3,maxclust),silave(maxclust),dab(maxclust),xben(maxclust)
     3 ,cen(maxclust),parc(maxclust),dunn(maxclust),cal(maxclust)
     4 ,str(maxclust),rex(maxclust)
       integer list(maxclust)
       common /clustnumber/compact,sepbar,sephat,silave,dab,xben,cen
     1  ,parc,dunn,cal,str,rex,/anclust/nclust,/cnnn/n2
       d1=-1.0d+10
       d2=1.0d+05
       d3=1.0d+05
       d4=-1.0d+10
       d5=-1.0d+10
       d6=-1.0d+10
       d7=-1.0d+10
       do i=n2,nclust
        if(silave(i).gt.d1) then
         d1=silave(i)
         nsil=i
        end if
        if(dab(i).lt.d2) then
         d2=dab(i)
         ndab=i
        end if
        if(xben(i).lt.d3) then
         d3=xben(i)
         nben=i
        end if
        if(dunn(i).gt.d4) then
         d4=dunn(i)
         ndunn=i
        end if
        if(cal(i).gt.d5) then
         d5=cal(i)
         ncal=i
        end if
        if(str(i).gt.d6) then
         d6=str(i)
         nstr=i
        end if
        if(rex(i).gt.d7) then
         d7=rex(i)
         nrex=i
        end if
       end do

       c1=-1.0d+10
       c2=-1.0d+10
       c3=-1.0d+10
       c4=-1.0d+10
       c5=-1.0d+10
       c6=-1.0d+10
       c7=-1.0d+10
       c8=-1.0d+10
       c9=-1.0d+10
       do i=n2,nclust
        if(compact(1,i).gt.c1) then
         c1=compact(1,i)
         ncomf=i
        end if
        if(compact(2,i).gt.c2) then
         c2=compact(2,i)
         ncoms=i
        end if
        if(compact(3,i).gt.c3) then
         c3=compact(3,i)
         ncomm=i
        end if
        if(sepbar(1,i).gt.c4) then
         c4=sepbar(1,i)
         nbarf=i
        end if
        if(sepbar(2,i).gt.c5) then
         c5=sepbar(2,i)
         nbars=i
        end if
        if(sepbar(3,i).gt.c6) then
         c6=sepbar(3,i)
         nbarm=i
        end if
        if(sephat(1,i).gt.c7) then
         c7=sephat(1,i)
         nhatf=i
        end if
        if(sephat(2,i).gt.c8) then
         c8=sephat(2,i)
         nhats=i
        end if
        if(sephat(3,i).gt.c9) then
         c9=sephat(3,i)
         nhatm=i
        end if
       end do
c-------------------------------------------------
c  full version
c-------------------------------------------------
        eps=2.0d-02
        l1=0
        do i=n2,nclust
         e1=c1-compact(1,i)
         if(e1.le.eps) then
          l1=l1+1
          list(l1)=i
         end if
        end do
        
        b1=-1.0d+10
        b2=-1.0d+10        
        do i=1,l1
         i1=list(i)
         if(sepbar(1,i1).gt.b1) then
          b1=sepbar(1,i1)
          nebar1=i1
         end if
         if(sephat(1,i1).gt.b2) then
          b1=sephat(1,i1)
          nehat1=i1
         end if
        end do
c---------------------------------------------------
c  mean+2*std version
c-------------------------------------------------
        l1=0
        do i=n2,nclust
         e1=c2-compact(2,i)
         if(e1.le.eps) then
          l1=l1+1
          list(l1)=i
         end if
        end do
        
        b3=-1.0d+10
        b4=-1.0d+10        
        do i=1,l1
         i1=list(i)
         if(sepbar(2,i1).gt.b3) then
          b3=sepbar(2,i1)
          nebar2=i1
         end if
         if(sephat(2,i1).gt.b4) then
          b4=sephat(2,i1)
          nehat2=i1
         end if
        end do
c---------------------------------------------------
c  median (core) version
c-------------------------------------------------
        l1=0
        do i=n2,nclust
         e1=c3-compact(3,i)
         if(e1.le.eps) then
          l1=l1+1
          list(l1)=i
         end if
        end do
        
        b5=-1.0d+10
        b6=-1.0d+10        
        do i=1,l1
         i1=list(i)
         if(sepbar(3,i1).gt.b5) then
          b5=sepbar(3,i1)
          nebar3=i1
         end if
         if(sephat(3,i1).gt.b6) then
          b6=sephat(3,i1)
          nehat3=i1
         end if
        end do
     
c------------------------------------------------------------
      end subroutine
