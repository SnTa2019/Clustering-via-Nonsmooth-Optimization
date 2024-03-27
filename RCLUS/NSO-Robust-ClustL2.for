c=============================================================
c RCLUST: Robust clustering algorithm: the use of soft trimming approach 
c=============================================================

      PARAMETER(maxvar=5000, maxrec=200000, maxclust=100, maxnft=100
     1 ,maxclass=100, maxsize=7000000)
      implicit double precision(a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),x2(maxsize)
     1 ,amed(maxnft),plabel(maxclass),xbest(maxvar),x3(maxvar)
     2 ,z(maxvar),x5(maxvar),fval(maxrec)
      integer nob(maxclass),list1(maxrec)
      common /csize/m,/c22/a,/anclust/nclust,/cnft/nft,/cmf/mf,/cnc/nc
     1 ,/crecord/nrecord,/cnob/plabel,nob,/ctnorm/tnorm,/cns/ns
     2 ,/cgamma/gamma1,gamma2,/cnpnout/npurity,/cnclass/nclass,/ctau/tau
     3 ,/ceps/eps,/cweight/weight,/ceps2/eps2,/clist1/list1
      open(39,file='Centers.txt')   
      open(40,file='Results.txt')
      open(78,file='inputdata.txt',status='old',form='formatted')        

c========================================================================
      nft = 21          !Number of features
      nclust = 10    !'Maximum number of clusters
      npurity = 2
      
c===========================================================
      tlimit=1.8d+04
      do i=1,maxrec
       read(78,*,END=901,ERR=900) (a(i,k),k=1,nft)  ! read input data 
       nrecord=i
      end do
  900 stop 'Error in input file'       
  901 WRITE(40,*) 'Input complete. Number of records: ',nrecord  !number of points
      WRITE(40,*)
      
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
 22     nob(j)=nob(j)+1
       end do
      END if
c================================================================
       WRITE(40,*)
       WRITE(40,701)
 701   FORMAT('#Clust','                  fval ','          CPU')
       WRITE(40,*) 
c================================================================
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
c==========================================
      eps=5.0d-01   
      eps2=5.0d-02  
      tau=5.0d-01    
      
      do nc=1,nclust
       if(nc.le.5) then
           eps=1.0d+00
           tau=5.0d-01
         else
           eps=1.0d+00
           tau=5.0d-01
       end if    

       if(nc.eq.1) then
           call step1(f,x)
           toler=1.0d-02*f/dble(nrecord)
           go to 1
       END if
       weight=tau*f/dble(nrecord)   
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
c================================================================
  1    call distribution(x,f)
       if(nc.gt.1) then
        write(39,*)
        write(39,*)
        write(39,940) nc
        write(39,*)
        do i=1,nc
         write(39,941) (x(j+(i-1)*mf),j=1,mf)
        end do
       end if
 940   format('The number of clusters:',i10)       
 941   format(14f11.4)   
       call cpu_time(time3)
       time4=time3-time1
       write(39,*)
       write(39,942) time4 
       write(39,*)
 942   format(' CPU time = ',f11.4)  
       write(40,240) nc,f,time4
       IF(time4.gt.tlimit) GO TO 2
      end do
   2  continue
 240  format(i6,f28.4,f12.4)
      close(39)      
      close(40)
      close(78)
      stop
      end
c==============================================================

      subroutine distribution(x,f)
      PARAMETER(maxdim=5000, maxrec=200000, maxclust=100, maxnft=100,
     1 maxclass=100)
      implicit double precision (a-h,o-z)
      double precision x(maxdim),a(maxrec,maxnft),dminim(maxrec)
     2 ,rad(maxclust),rad2(maxclust)
      integer nk(maxclust,maxrec),nel(maxclust),list1(maxrec)
     1 ,lcand(maxrec),lcand1(maxrec),lact3(maxrec)
      common /c22/a,/cmf/mf,/anclust/nclust,/crecord/nrecord,/crad/rad
     1 ,/cnc/nc,/adminim/dminim,/clcand/ncand,lcand,/cnk/nk,nel
     2 ,/ctnorm/tnorm,/clist1/list1,/ctoler/toler,/crscale/rsc
     3 ,/clast3/lact3,/ceps/eps,/ctau/tau,/ceps3/eps3
c==============================================================
       do j=1,nclust
        nel(j)=0
        rad(j)=0.0d+00
        rad2(j)=0.0d+00
       END do
      
       f=0.0d+00
       do k=1,nrecord
        f20=1.0d+28
        do j=1,nc
         f1=0.0d+00
         do k1=1,mf
          f1=f1+(a(k,k1)-x(k1+(j-1)*mf))**2
         END do
         tnorm=tnorm+1.0d+00
         if(f20.gt.f1) then
          f20=f1
          jmin=j
         end if
        END do
        dminim(k)=f20
        f=f+f20
        nel(jmin)=nel(jmin)+1
        nk(jmin,nel(jmin))=k
        list1(k)=jmin
        rad(jmin)=rad(jmin)+f20
        rad2(jmin)=rad2(jmin)+dble(lact3(k))*f20
       end do

       if(nc.gt.1) then
        write(39,*)
        write(39,*)
        write(39,940) nc
        write(39,*)
        write(39,942) (nel(i),i=1,nc)        
        write(39,*)        
       end if
 940   format('The number of clusters:',i10)       
 942   format('Distribution of points over clusters:',10i8) 

       do j=1,nc
        if(nel(j).gt.0) then
          rad(j)=rad(j)/dble(nel(j))
          rad2(j)=rad2(j)/dble(nel(j))
         else
          rad(j)=0.0d+00
          rad2(j)=0.0d+00
        end if
       end do

       rsc1=0.0d+00
       do i=1,nc
        rsc1=rsc1+rad2(i)
       end do
       rsc1=rsc1/dble(nc)
       tau1=5.5d-01
       tau2=6.0d-01
       if(nc.le.3) rsc=tau1*rsc1
       if(nc.gt.3) rsc=tau2*rsc1
 
c==================================================================
      n2=0
      do i=1,nrecord
       k=list1(i)
       if(dminim(i).gt.rad(k)) n2=n2+1
      end do
      s1=dble(n2)/dble(nrecord)
c==================================================================
       if(nc.lt.nclust) then
        ncand=0
        do i=1,nc
         k1=nel(i)
         do j=1,k1
          k2=nk(i,j)
          if(dminim(k2).gt.rad(i)) then
            ncand=ncand+1
            lcand(ncand)=k2
          end if
         end do
        end do
        ncand1=0
        do i=1,ncand
         i1=lcand(i)
         do j=1,ncand1
          j1=lcand1(j)
          d1=0.0d+00
          do k=1,mf
           d1=d1+(a(i1,k)-a(j1,k))**2
          end do
          if(d1.le.toler) go to 1
         end do
         ncand1=ncand1+1
         lcand1(ncand1)=i1
   1    end do

        ncand=ncand1
        do i=1,ncand
         lcand(i)=lcand1(i)
        end do
       end if
c========================================================
      return
      end

c=====================================================
c Step1 calculates the center of the dataset
c=====================================================
      subroutine step1(f,x)
      PARAMETER(maxvar=5000, maxrec=200000, maxnft=100)
      implicit double precision (a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),dist2(maxrec)
     1 ,x2(maxnft),dist1(maxrec) 
      integer lact3(maxrec)
      common /c22/a,/cmf/mf,/crecord/nrecord,/ctnorm/tnorm,/ceps2/eps2
     1 ,/clast3/lact3,/ceps/eps,/ctau/tau,/crscale/rsc,/ceps3/eps3 
      do i=1,mf
       x(i)=0.0d+00
       do j=1,nrecord
        x(i)=x(i)+a(j,i)
       END do
       x(i)=x(i)/dble(nrecord)
      END do

      f=0.0d+00
      r2=0.0d+00
      do i=1,nrecord
       f1=0.0d+00
       do j=1,mf
        f1=f1+(a(i,j)-x(j))*(a(i,j)-x(j))
       END do
       f=f+f1
       dist2(i)=f1
       r2=dmax1(r2,dist2(i))
      END do
      tnorm=tnorm+dble(nrecord)
c==============================================================
      do i=1,mf
       x2(i)=0.0d+00
      end do
      r1=f/dble(nrecord)

      r3=r1+8.0d-01*(r2-r1)      
      n2=0
      do i=1,nrecord
       if(dist2(i).le.r3) then
        n2=n2+1
        do j=1,mf
         x2(j)=x2(j)+a(i,j)
        end do
       end if
      end do 
      do j=1,mf
       x2(j)=x2(j)/dble(n2)
      end do

      f2=0.0d+00
      do i=1,nrecord
       if(dist2(i).le.r1) then
        f3=0.0d+00
        do j=1,mf
         f3=f3+(x2(j)-a(i,j))**2
        end do
        f2=f2+f3
       end if
      end do 
      rsc=f2/dble(nrecord)

c====================================================================
      r3=0.0d+00
      r4=0.0d+00
      do i=1,nrecord
       d2=0.0d+00
       do j=1,nrecord
        if(j.ne.i) then
         d1=0.0d+00
         do k=1,mf
          d1=d1+(a(i,k)-a(j,k))**2
         end do
         d2=d2+d1
        end if
       end do
       dist1(i)=d2/dble(nrecord-1)
       r3=dmax1(r3,dist1(i))
       r4=r4+dist1(i)
      end do
      r4=r4/dble(nrecord)
      r5=r3/r4

      if(r5.le.3.0d+00) tau3=9.5d-01
      if((r5.gt.3.0d+00).and.(r5.le.5.0d+00)) tau3=9.0d-01
      if((r5.gt.5.0d+00).and.(r5.le.1.0d+01)) tau3=8.5d-01
      if((r5.gt.1.0d+01).and.(r5.le.2.5d+01)) tau3=7.5d-01
      if((r5.gt.2.5d+01).and.(r5.le.5.0d+01)) tau3=3.0d-01
      if((r5.gt.5.0d+01).and.(r5.le.1.0d+02)) tau3=6.0d-01
      if((r5.gt.1.0d+02).and.(r5.le.2.5d+02)) tau3=5.0d-01
      if((r5.gt.2.5d+02).and.(r5.le.5.0d+02)) tau3=8.0d-03
      if((r5.gt.5.0d+02).and.(r5.le.1.0d+03)) tau3=3.0d-02
      if(r5.gt.1.0d+03) tau3=2.0d-03
      r3=tau3*r3
      do i=1,nrecord
       lact3(i)=0
      end do
      
      do i=1,nrecord
       if(dist1(i).le.r3) then
        lact3(i)=1
       end if
      end do 
c=====================================================================
      n1=0
      do i=1,nrecord
       if(lact3(i).eq.0) n1=n1+1
      end do
      print *,n1,n2
           
c=====================================================================      
      return
      end

c=========================================================================
c  Step2 computes clusters for each data point
c=========================================================================
      subroutine step2(toler,nstart,x2)
      PARAMETER(maxdim=5000, maxrec=200000, maxclust=100, maxnft=100
     1 ,maxsize=7000000)
      implicit double precision (a-h,o-z)
      double precision x2(maxsize),a(maxrec,maxnft),dminim(maxrec)
     1 ,fmin1(maxrec),x4(maxnft),fval(maxrec)
      integer l4(maxrec),lcand(maxrec),lcand1(maxrec)
      common /c22/a,/cmf/mf,/adminim/dminim,/crecord/nrecord
     1 ,/ctnorm/tnorm,/clcand/ncand,lcand,/cgamma/gamma1,gamma2

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
      end subroutine

c=====================================================================
      subroutine auxfunc(f,x)
      PARAMETER(maxvar=5000, maxrec=200000, maxclust=100, maxnft=100)
      implicit double precision (a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),dminim(maxrec)
     1 ,dist1(maxrec) 
      integer iact(maxrec)
      common /c22/a,/cmf/mf,/crecord/nrecord,/adminim/dminim,/cic/iact
     1 ,/ctnorm/tnorm

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
        f4=dmin1(dminim(i),f5)
        if(f4.eq.dminim(i)) iact(i)=1
        if(f4.eq.f5) iact(i)=2
        f=f+f4
       end do
      return
      end
      
c==============================================================

      subroutine auxgrad(x,grad)
      PARAMETER(maxvar=5000, maxrec=200000, maxclust=100, maxnft=100
     1 ,maxcannot=10000, maxmust=10000)
      implicit double precision (a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),dminim(maxrec)
     1 ,grad(maxvar)
      integer iact(maxrec)  
      common /c22/a,/cmf/mf,/crecord/nrecord,/adminim/dminim
     1 ,/cic/iact
     
      do i=1,mf
       grad(i)=0.0d+00
      end do
      do i=1,nrecord
       if(iact(i).eq.2) then
        do j=1,mf
         grad(j)=grad(j)+2.0d+00*(x(j)-a(i,j))
        end do
       end if 
      end do
      return
      end

c=====================================================================
      subroutine func(f,x)
      PARAMETER(maxvar=5000, maxrec=200000, maxclust=100, maxnft=100)
      implicit double precision (a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),dist2(maxrec,maxclust)
      integer lact(maxrec),lact2(maxrec),lact3(maxrec)
      common /c22/a,/cmf/mf,/crecord/nrecord,/cnc/nc,/ctnorm/tnorm
     2 ,/cic1/lact,/ceps/eps,/cdist2/dist2,/cweight/weight,/crscale/rsc
     3 ,/cic2/lact2,/clast3/lact3

      do i=1,nrecord
       do j=1,nc
        dist2(i,j)=0.0d+00
        do k=1,mf
         dist2(i,j)=dist2(i,j)+(x(k+(j-1)*mf)-a(i,k))**2
        end do
       end do
      end do
      tnorm=tnorm+dble(nrecord)*dble(nc)
   
      f=0.0d+00
      do i=1,nrecord
       f4=1.0d+30
       do k=1,nc
        f5=dist2(i,k)
        if(f4.gt.f5) then
         f4=f5
         lact(i)=k
        end if                 
       end do
       f=f+dble(lact3(i))*f4
      end do

      do i=1,nrecord
       lact2(i)=0
      end do

      h=0.0d+00
      do i=1,nrecord
       if(lact3(i).eq.1) then
        k=lact(i)
        ta=dist2(i,k)/rsc
        w=1.0d+00+(1.7d+00-ta)/sqrt(6.0d-01+(1.5d+00-ta)**2)
        w=-w+eps
        h=h+dmax1(0.0d+00,w)
        if(w.ge.0.0d+00) lact2(i)=1
       end if
      end do
      f=f+weight*h
      return
      end
c==============================================================

      subroutine funcgrad(x,grad)
      PARAMETER(maxvar=5000, maxrec=200000, maxclust=100, maxnft=100)
      implicit double precision (a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),grad(maxvar)
     1 ,dist2(maxrec,maxclust) 
      integer lact(maxrec),lact2(maxrec),lact3(maxrec)
      common /c22/a,/cmf/mf,/crecord/nrecord,/cnc/nc,/csize/m,/cic1/lact
     1 ,/cdist2/dist2,/ceps/eps,/crscale/rsc,/cic2/lact2,/cweight/weight
     2 ,/clast3/lact3
      call func(f,x)
      do i=1,m
       grad(i)=0.0d+00
      end do     
      
      do i=1,nrecord
       j1=lact(i)
       do j=1,mf
        k1=(j1-1)*mf+j
        grad(k1)=grad(k1)+2.0d+00*dble(lact3(i))*(x(k1)-a(i,j))
       end do
      end do
      
      do i=1,nrecord
       jj=lact2(i)*lact3(i)
       if(jj.eq.1) then
        j1=lact(i)
        ta=dist2(i,j1)/rsc
        h1=6.0d-01+(ta-1.5d+00)**2
        h2=-h1+(ta-1.7d+00)*(ta-1.5d+00)
        ht=h2/(h1*sqrt(h1))
        do j=1,mf
          k1=(j1-1)*mf+j
          grad(k1)=grad(k1)-2.0d+00*ht*weight*(x(k1)-a(i,j))/rsc
        end do
       end if
      end do
     
      return
      end

c=====================================================================
      subroutine dgm(x,f2)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=5000, maxdg=1000,maxit=100000)
      double precision x(maxvar),x1(maxvar),g(maxvar),v(maxvar)
     1 ,w(maxdg,maxvar),prod(maxdg,maxdg),z(maxdg),fvalues(maxit)
      INTEGER ij(maxdg)
      common /csize/m,/cij/ij,jvertex,/cz/z,/ckmin/kmin
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
c  programming problem, to find descent direction
c===============================================================

      subroutine wolfe(ndg,prod)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=5000, maxdg=1000)
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
c==============================================================

      subroutine equations(n,z1)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=5000, maxdg=1000)
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
      subroutine dgrad(x,sl,g,dg)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=5000, maxdg=1000)
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
c Line search (Armijo-type)
c===========================================================
      subroutine armijo(x,g,f1,f5,f4,sl,step,r)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=5000, maxdg=1000)
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
c==============================================================

      subroutine fv(x,f)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=5000)
      double precision x(maxvar)
      common /cns/ns
      if(ns.eq.1) call auxfunc(f,x)
      if(ns.eq.2) call func(f,x)
      return
      end
