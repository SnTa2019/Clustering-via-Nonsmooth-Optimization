!*************************************************************************
!*                                                                       *
!*     Subroutines for LMBM-Clust                                        *
!*     (version 2.1, last modified 11.11.2017).                          *
!*                                                                       *
!*************************************************************************
!*
!*     Modules included:
!*
!*     clusteringmod       !
!*

MODULE clusteringmod       ! Subroutines for clustering software

    USE r_precision, ONLY : prec      ! Precision for reals.
    IMPLICIT NONE

    ! MODULE clusteringmod includes the following subroutines (S) and functions (F).

    PUBLIC :: &
        check, &           ! S  Checking the results and and calculation of
                           !    validity indices.
        step1, &           ! S  Computation of the centroid and the value of
                           !    clustering function at centroid.
        step2              ! S  Computation of clusters for each data point.
!    PRIVATE :: &

CONTAINS

    SUBROUTINE check(x,db,db2,dn,dn3,sep)

        USE param, ONLY : zero, one, two, large
        USE initclust, ONLY : &
            maxdim, &         ! Maximum number of variables in optimization,
                              !   maxdim >= nclust * nft (parameter).
            a, &              ! Data matrix.
            mf, &             ! Number of used features:
                              !   mf = nft, when no classes,
                              !   mf = nft - 1, when classes.
            nclust, &         ! Maximum number of clusters.
            tnorm, &          ! The total number of norms computed.
            nrecord, &        ! Number of instances in data.
            m1, &             ! Maximum number of initial solutions.
            nc, &             ! Current number of clusters, loops from 1 to nclust.
            nk, &             ! nk(nclust,nrecord)
            nel, &            ! nel(nclust), nel(i)=number of records in cluster i
            ncand, &          ! Number of canditate points.
            lcand, &          ! lcand(nrecord)
            list1, &          ! list1(nrecord), list1(i)=the cluster where point i belongs
            dminim, &         ! dminim(nrecord), the distance of a(i) and the nearest centroid
            dcent             ! dcent(nclust,nclust), Distance (affinity) matrix for cluster centers

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::      &
            db,db2,dn,sep,      & ! Davies Boulding and Dunn indices
            dc(nclust,nclust),  &
            dc1,dn1,dn2,dn3,radm, &
            x(maxdim)

        ! Local variables
        REAL(KIND=prec) ::      &
            f10,                &
            f,f1,f2,            &
            fdb,                &
            fk2,fk3,            &
            fm,fm3,             &
            fk(nclust),         &
            rad(nclust),        & ! Average distance of data points in cluster i to its centroid.
            radmax(nclust),     & ! Maximum distance of data points in cluster i to its centroid.
            r1,r2,r3,r4,r5

        ! Local integer variables
        INTEGER ::  &
            j,jmin, &
            i,      &
            k,      &
            k1,     &
            ncand1, &
            n1


        DO j=1,nclust
            nel(j)=0
            rad(j)=zero
            radmax(j)=zero
        END DO
        
        f=zero
        jmin=1
        outerloop: DO k=1,nrecord
            f2=large
            f10=large
            innerloop: DO j=1,nc
                IF(j > 1) THEN
                    IF(dcent(j,jmin) >= f10) CYCLE innerloop
                END IF
                f1=zero
                DO k1=1,mf
                    f1=f1+(a(k1,k)-x(k1+(j-1)*mf))**2

                END DO
                tnorm=tnorm+one
                IF(f1 < f2) THEN
                    f2=f1
                    f10=4.0_prec*f1
                    jmin=j
                END IF
            END DO innerloop
            dminim(k)=f2
            f=f+f2
            nel(jmin)=nel(jmin)+1
            nk(jmin,nel(jmin))=k
            list1(k)=jmin
            rad(jmin)=rad(jmin)+f2
            radmax(jmin)=MAX(radmax(jmin),f2)
        END DO outerloop

        DO k=1,nc
            IF(nel(k) > 0) THEN
                rad(k)=rad(k)/REAL(nel(k),prec)
            ELSE
                rad(k)=zero
            END IF
        END DO

        ncand=0
        loop1: DO k=1,nc
            IF(nel(k) == 1) CYCLE loop1
            n1=nel(k)*m1/nrecord

            n1=MAX(1,n1)
            ncand1=0
            r1=radmax(k)
            r2=rad(k)
            r5=1.0E-03_prec*(r1-r2)
            IF(ABS(r1-r2) < 1.0E-06_prec) r5=1.0E-06_prec
            r3=r1+r5
            DO
                r3=r3-r5
                r4=r3-r5
                IF(r3 < r2) CYCLE loop1
                DO j=1,nel(k)
                    i=nk(k,j)
                    IF((dminim(i) <= r3).and.(dminim(i) > r4)) THEN
                        ncand=ncand+1
                        ncand1=ncand1+1
                        lcand(ncand)=i
                        IF(ncand1 >= n1) CYCLE loop1
                    END IF
                END DO
            END DO
        END DO loop1

        !=====================================================
        ! Calculation of distances between cluster centers
        !=====================================================

        do i=1,nc
            dc(i,i) = zero
        end do

        do i=1,nc
            do j=i+1,nc
                dc1 = zero
                do k=1,mf
                    dc1=dc1+(x(k+(i-1)*mf)-x(k+(j-1)*mf))**2
                end do
                dc(i,j)=SQRT(dc1)
                dc(j,i)=dc(i,j)
            end do
        end do


        !=====================================================
        ! Calculation of Davies-Bouldin (DB) validity index
        ! (last changed 25 August 2017)
        ! db is the index using Euclidean distances while
        ! db2 uses squared Euglidean distances
        !=====================================================

        DO i=1,nc
            fk(i)=zero
        END DO

        fdb=zero
        db2=zero

        DO i=1,nc
            fk(i)=SQRT(rad(i))
        END DO


        DO k=1,nc
            fm=zero
            fm3=zero
            DO i=1,nc
                IF (i.ne.k) THEN
                    fk2=fk(i)+fk(k)
                    fk3=(rad(i)+rad(k))/(dc(i,k)*dc(i,k))
                    f2=fk2/dc(i,k)
                    fm=MAX(fm,f2)
                    fm3=MAX(fm3,fk3)
                END IF
            END DO
            fdb=fdb+fm
            db2=db2+fm3
        END DO
        db=fdb/REAL(nc,prec)
        db2=db2/REAL(nc,prec)


        !============================================================
        ! Calculation of Dunn validity index
        !============================================================
        radm= zero
        DO i=1,nc
            radm=MAX(radm,radmax(i)) ! maximum radius of all clusters
        END DO
        radm = SQRT(radm)


        dn3 = zero
        dn = large
        DO i=1,nc
            dn1 = large
            DO j=1,nc
                IF (j.NE.i) THEN
                    dn2=dc(i,j)/radm
                    dn1=MIN(dn2,dn1) ! distance of cluster i to the closest cluster j
                END IF
            END DO
            dn=MIN(dn,dn1)
            dn3=MAX(dn,dn1)
        END DO
        IF (nc == 1) dn=zero
        IF (nc == 1) dn3=zero

        !============================================================
        ! Calculation of the quality of separation
        !============================================================
        sep=zero
        !        do i=1,nc
        !            m2=nel(i)
        !            m3=0
        !            loop_sep: do j=1,m2
        !                k=nk(i,j)
        !                do k1=1,nc
        !                    if(k1.ne.i) then
        !                        d1=zero
        !                        do k2=1,mf
        !                            d1=d1+(a(k2,k)-x(k2+(k1-1)*mf))**2
        !                            !d1=d1+(a(k,k2)-x(k2+(k1-1)*mf))**2
        !                        end do
        !                        if(d1<radmax(k1)) then
        !                            m3=m3+1
        !                            CYCLE loop_sep
        !                        end if
        !                    end if
        !                end do
        !                end do loop_sep
        !                sep=sep+REAL(m3,prec)/REAL(nel(i),prec)
        !            end do
        !            sep=sep/REAL(nc,prec)
        !============================================================


        RETURN

    END SUBROUTINE check




    !===============================================================================
    !  Step1 computes the centroid and the value of clustering function at centroid
    !===============================================================================
    SUBROUTINE step1(f,x)

        USE param, ONLY : zero, one
        USE initclust, ONLY : &
            a, &           ! Data matrix
            mf, &          ! Number of used features:
                           !   mf = nft, when no classes,
                           !   mf = nft - 1, when classes.
            nrecord, &     ! Number of reports in data
            dminim, &      ! Distance of a(i) and the nearest centroid.
            tnorm          ! The total number of norms computed.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            f,              &
            x(mf)

        ! Local variables
        REAL(KIND=prec) ::  &
            f1

        INTEGER ::      &
            i,j


        DO i=1,mf
            x(i)=zero
            DO j=1,nrecord
                x(i)=x(i)+a(i,j) ! Note, here is inefficient order of array a(i,j)
            END DO
            x(i)=x(i)/REAL(nrecord,prec)
        END DO
        f=zero

        DO i=1,nrecord
            f1=zero
            tnorm=tnorm+one
            DO j=1,mf
                f1=f1+(a(j,i)-x(j))*(a(j,i)-x(j))
            END DO
            dminim(i)=f1
            f=f+f1
        END DO

        RETURN

    END SUBROUTINE step1


    !=========================================================================
    !  Step2 computes clusters for each data point
    !=========================================================================
    SUBROUTINE step2(toler,nstart,x2)

        USE param, ONLY : zero, one, two
        USE initclust, ONLY : &
            maxsize, &     !
            a, &           ! Data matrix
            mf, &          ! Number of used features:
                           !   mf = nft, when no classes,
                           !   mf = nft - 1, when classes.
            nrecord, &     ! Number of reports in data.
            dminim, &      ! Distance of a(i) and the nearest centroid.
            ncand, &       ! Number of canditate points.
            lcand, &       !
            ng1, &         !
            tnorm          ! The total number of norms computed.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            toler,          &
            dcand(ncand,ncand), &
            x2(maxsize)

        INTEGER ::  &
            nstart

        ! Local variables
        REAL(KIND=prec) ::  &
            d1,d3,d4,d21,   &
            fmin,fmin2,     &
            fmin1(nrecord), &
            f1,f10,f2,      &
            fmin0,          &
            gamma0,         &
            toler2

        INTEGER ::           &
            ncand1,          &
            nclose,          &
            nstart0,         &
            l4(nrecord),     &
            lcand1(nrecord), &
            icand,           &
            i1,              &
            i,j,             &
            j1,j2,           &
            k,l

        INTRINSIC :: MAX,MIN      ! mimimum

        nstart0=MAX(1,ng1*ncand/10)

        nstart=0
        fmin=zero

        DO i1=1,ncand
            i=lcand(i1)
            f2=dminim(i)
            d21=zero
            loop1: DO l=1,nrecord
                IF (f2<(4.0_prec*dminim(l))) THEN                     d3=zero
                    tnorm=tnorm+one
                    DO k=1,mf
                        d3=d3+(a(k,i)-a(k,l))**2
                    END DO
                    d21=d21+MIN(zero,d3-dminim(l))
                END IF
            END DO loop1
            fmin1(i)=d21
            fmin=MIN(fmin,d21)
        END DO

        ncand1=0
        gamma0=1.01_prec

        loop2: DO WHILE (ncand1 < nstart0)
            gamma0=gamma0-0.01_prec
            fmin2=(gamma0-0.01_prec)*fmin
            fmin0=gamma0*fmin

            DO i1=1,ncand
                i=lcand(i1)
                IF ((fmin1(i) >= fmin0).AND.(fmin1(i) < fmin2)) THEN
                    ncand1=ncand1+1
                    lcand1(ncand1)=i
                    IF(ncand1>=nstart0) EXIT loop2
                END IF
            END DO
        END DO loop2

        ncand=ncand1
        DO i1=1,ncand
            lcand(i1)=lcand1(i1)
        END DO

        icand=0
        DO i1=1,ncand
            i=lcand(i1)
            nclose=0
            DO j=1,nrecord
                IF(dminim(i)<(4.0_prec*dminim(j))) THEN
                    d1=zero
                    tnorm=tnorm+one
                    DO k=1,mf
                        d1=d1+(a(k,i)-a(k,j))**2
                    END DO
                    IF(d1 < dminim(j)) THEN
                        nclose=nclose+1
                        l4(nclose)=j
                    END IF
                END IF
            END DO

            IF(nclose>0) THEN
                DO k=1,mf
                    d3=zero
                    DO j=1,nclose
                        j1=l4(j)
                        d3=d3+a(k,j1)
                    END DO
                    x2(k+(i1-1)*mf)=d3/REAL(nclose,prec)
                END DO
            END IF
        END DO

        DO i=1,ncand
            dcand(i,i)=zero
        END DO

        DO i=1,ncand
            DO j=i+1,ncand
                dcand(i,j)=zero
                DO k=1,mf
                    dcand(i,j)=dcand(i,j)+(x2(k+(i-1)*mf)-x2(k+(j-1)*mf))**2
                END DO
                dcand(j,i)=dcand(i,j)
            END DO
        END DO

        toler2=0.1_prec*toler
        nstart=0
        loop_ncand: DO j=1,ncand
            DO j1=1,nstart
                j2=l4(j1)
                d4=dcand(j,j2)
                IF(d4 > toler2) CYCLE loop_ncand

            END DO
            nstart=nstart+1
            l4(nstart)=j
            DO k=1,mf
                x2(k+(nstart-1)*mf)=x2(k+(j-1)*mf)
            END DO
        END DO loop_ncand

        RETURN

    END SUBROUTINE step2

END MODULE clusteringmod
