!*************************************************************************
!*                                                                       *
!*     Subroutines for DC optimization based clustering software         *
!*     (last modified 24.02.2016).                                       *
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
                           !    Davies-Bouldin (DB) validity index.
        step1, &           ! S  Computation of the centroid and the value of
                           !    clustering function at centroid.
        step2              ! S  Computation of clusters for each data point.
!    PRIVATE :: &

CONTAINS

    SUBROUTINE check(x,db)

        USE param, ONLY : zero
        USE initclust, ONLY : &
            maxdim, &         ! Maximum number of variables in optimization,
                              !   maxdim >= maxclust * maxnft (parameter).
            maxrec, &         ! maximum number of reports in dataset (parameter).
            maxclust, &       ! maximum number of clusters (parameter).
            ! maxclass, &       ! maximum number of classed (not used).
            a, &              ! Data matrix.
            mf, &             ! Number of used features:
                              !   mf = nft, when no classes,
                              !   mf = nft - 1, when classes.
            nclust, &         ! Maximum number of clusters (given by user).
            nrecord, &        ! Number of reports in data  (given by user).
            nc, &             ! Current number of clusters, loops from 1 to nclust.
            nk, &             ! nk(maxclust,maxrec)
            nel, &            ! nel(maxclust)
            ncand, &          !
            lcand, &          ! lcand(maxrec)
            dminim            ! dminim(maxrec)
            ! nclass, &       ! Number of classes, from user (not used).
            ! plabel, &       ! Index for class output plabel(maxclass), (not used).
            ! nob, &          ! Classes, nob(maxclass) (not used).
            ! npurity, &      ! Switch for classes, from user, (not used):
            !                 !   npurity = 1 if class outputs,
            !                 !   npurity = 2 if no class outputs.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::      &
            db,                 &
            x(maxdim)

        ! Local variables
        REAL(KIND=prec) ::      &
            f12,                &
            f2,                 &
            fdb,                &
            fk2,                &
            fm,                 &
            ! purity,           & ! no classes, not used
            f1(maxclust),       &
            fk(maxclust),       &
            rad(maxclust)

        ! Local integer variables
        INTEGER ::  &
            j,      &
            i,      &
            j1,     &
            k,      &
            k1,     &
            n2,     &
            n3,     &
            ! nob1(maxclass),   &  ! No classes, not used
            list1(maxrec)

        DO j=1,nclust
            nel(j)=0
            rad(j)=zero
        END DO

        outerloop: DO k=1,nrecord
            DO j=1,nc
                f1(j)=zero
                DO k1=1,mf
                    f1(j)=f1(j)+(a(k,k1)-x(k1+(j-1)*mf))**2
                END DO
            END DO
            f2=f1(1)
            DO  j=2,nc
                f2=dmin1(f2,f1(j))
            END DO
            dminim(k)=f2
            DO j=1,nc
                IF (f1(j)==f2) THEN
                    nel(j)=nel(j)+1
                    nk(j,nel(j))=k
                    list1(k)=j
                    rad(j)=rad(j)+f2
                    CYCLE outerloop
                END IF
            END DO
        END DO outerloop

        DO k=1,nc
            IF (nel(k) > 0) THEN
                rad(k)=7.5E-01_prec*rad(k)/REAL(nel(k),prec)
            ELSE
                rad(k)=zero
            END IF
        END DO

        ncand=0
        DO k=1,nc
            DO j=1,nel(k)
                i=nk(k,j)
                IF(dminim(i) > rad(k)) THEN
                    ncand=ncand+1
                    lcand(ncand)=i
                END IF
            END DO
        END DO


        !=====================================================
        ! Calculation of Davies-Bouldin (DB) validity index
        !=====================================================
        DO i=1,nc
            fk(i)=zero
        END DO

        fdb=zero
        DO i=1,nrecord
            k=list1(i)
            fk(k)=fk(k)+dminim(i)
        END DO
        DO i=1,nc
            IF(nel(i) > 0) fk(i)=fk(i)/REAL(nel(i),prec)
        END DO
        DO k=1,nc
            fm=zero
            DO i=1,nc
                IF (i.ne.k) THEN
                    fk2=fk(i)+fk(k)
                    f12=zero
                    DO j=1,mf
                        f12=f12+(x(j+(i-1)*mf)-x(j+(k-1)*mf))**2
                    END DO
                    f2=fk2/f12
                    fm=dmax1(fm,f2)
                END IF
            END DO
            fdb=fdb+fm
        END DO
        db=fdb/REAL(nc,prec)

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
            dminim, &      !
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
            i,          &
            j


        DO i=1,mf
            x(i)=zero
            DO j=1,nrecord
                x(i)=x(i)+a(j,i)
            END DO
            x(i)=x(i)/REAL(nrecord,prec)
        END DO
        f=zero

        DO i=1,nrecord
            f1=zero
            tnorm=tnorm+one
            DO j=1,mf
                f1=f1+(a(i,j)-x(j))*(a(i,j)-x(j))
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

        USE param, ONLY : zero, one
        USE initclust, ONLY : &
            maxsize, &     !
            a, &           ! Data matrix
            mf, &          ! Number of used features:
                           !   mf = nft, when no classes,
                           !   mf = nft - 1, when classes.
            nrecord, &     ! Number of reports in data.
            dminim, &      !
            ncand, &       !
            lcand, &       !
            gamma1, &      ! Parameter affecting the size of the set of starting point,
            gamma2,&       ! Parameter affecting the size of the set of starting point,
            tnorm          ! The total number of norms computed.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            toler,          &
            x2(maxsize)

        INTEGER ::  &
            nstart

        ! Local variables
        REAL(KIND=prec) ::  &
            d1,             &
            d2,             &
            d3,             &
            d4,             &
            d21,            &
            fmin,           &
            fmin2,          &
            x4(mf),         &
            fval(nrecord),  &
            fmin1(nrecord)

        INTEGER ::           &
            ncand1,          &
            nclose,          &
            nstart1,         &
            l4(nrecord),     &
            lcand1(nrecord), &
            i1,              &
            i,               &
            j,               &
            j1,              &
            k,               &
            l

        INTRINSIC :: DMIN1      ! mimimum

        nstart=0
        fmin=zero
        DO i1=1,ncand
            i=lcand(i1)
            d21=zero
            DO l=1,nrecord
                d3=zero
                tnorm=tnorm+one
                DO k=1,mf
                    d3=d3+(a(i,k)-a(l,k))**2
                END DO
                d21=d21+DMIN1(zero,d3-dminim(l))
            END DO
            fmin1(i)=d21
            fmin=dmin1(fmin,d21)
        END DO

        fmin2=gamma1*fmin
        ncand1=0
        DO i1=1,ncand
            i=lcand(i1)
            IF (fmin1(i) <= fmin2) THEN
                ncand1=ncand1+1
                lcand1(ncand1)=i
            END IF
        END DO

        ncand=ncand1
        DO i1=1,ncand
            lcand(i1)=lcand1(i1)
        END DO

        outerloop: DO i1=1,ncand
            i=lcand(i1)
            nclose=0
            DO j=1,nrecord
                d1=zero
                tnorm=tnorm+one
                DO k=1,mf
                    d1=d1+(a(i,k)-a(j,k))**2
                END DO
                IF(d1 < dminim(j)) THEN
                    nclose=nclose+1
                    l4(nclose)=j
                END IF
            END DO

            IF(nclose == 0) CYCLE outerloop
            DO k=1,mf
                d3=zero
                DO j=1,nclose
                    j1=l4(j)
                    d3=d3+a(j1,k)
                END DO
                x4(k)=d3/REAL(nclose,prec)
            END DO
            DO j=1,nstart
                d4=zero
                tnorm=tnorm+one
                DO k=1,mf
                    d4=d4+(x2(k+(j-1)*mf)-x4(k))**2
                END DO
                IF(d4 <= toler) CYCLE outerloop
            END DO

            nstart=nstart+1
            DO k=1,mf
                x2(k+(nstart-1)*mf)=x4(k)
            END DO
        END DO outerloop

        d2=zero
        DO j=1,nstart
            d21=zero
            DO l=1,nrecord
                d3=zero
                tnorm=tnorm+one
                DO k=1,mf
                    d3=d3+(x2(k+(j-1)*mf)-a(l,k))**2
                END DO
                d21=d21+DMIN1(zero,d3-dminim(l))
            END DO
            fval(j)=d21
            d2=dmin1(d2,d21)
        END DO

        d2=gamma2*d2
        nstart1=0

        DO j=1,nstart
            IF (fval(j) <= d2) THEN
                nstart1=nstart1+1
                DO k=1,mf
                    x2(k+(nstart1-1)*mf)=x2(k+(j-1)*mf)
                END DO
            END IF
        END DO

        nstart=nstart1

        RETURN

    END SUBROUTINE step2

END MODULE clusteringmod
