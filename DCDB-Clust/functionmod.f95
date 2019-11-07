!*************************************************************************
!*                                                                       *
!*     Computation of function and (sub)gradients values of              *
!*     minimum sum of squares clustering problem given in DC form        *
!*     (last modified 24.02.2016).                                       *
!*                                                                       *
!*************************************************************************
!*
!*     Modules included:
!*
!*     functionmod         !

MODULE functionmod ! Computation of function and (sub)gradients values for clustering software

    USE r_precision, ONLY : prec      ! Precision for reals.
    IMPLICIT NONE

    ! MODULE functionmod includes the following subroutines (S) and functions (F).

    PUBLIC :: &
        auxfunc, &         ! S Computation of aux clustering problem
        clusterfunc, &     ! S Computation of clustering problem
        fgrad1, &          ! S Computation of the gradient of the first component function
        fgrad2             ! S Computation of the subgradient of the second component function
!    PRIVATE :: &


CONTAINS

    !=============================================
    ! Computation of aux clustering problem
    !=============================================
    SUBROUTINE auxfunc(x,fval)

        USE param, ONLY : zero, one
        USE initclust, ONLY : &
            maxdim, &      ! Maximum number of variables in optimization,
                           !   maxdim >= maxclust * maxnft.
            maxrec, &      ! maximum number of recorts in dataset
            maxclust, &    ! maximum number of clusters
            maxnft, &      ! maximum number of features in dataset
            a, &           ! Data matrix
            mf, &          ! Number of used features:
                           !   mf = nft, when no classes,
                           !   mf = nft - 1, when classes.
            nrecord, &     ! Number of reports in data.
            dminim, &      ! dminim(maxrec)
            tnorm          ! The total number of norms computed.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            x(maxdim), &
            fval

        ! Local variables
        REAL(KIND=prec) ::  &
            f1,f2,f3
        INTEGER :: &
            i,j

        INTRINSIC ::  DMIN1

        fval=zero
        DO i=1,nrecord
            f1=dminim(i)
            tnorm=tnorm+one
            f2=zero
            DO j=1,mf
                f2=f2+(a(i,j)-x(j))**2
            END DO
            f3=DMIN1(f1,f2)
            fval=fval+f3
        END DO
        RETURN

    END SUBROUTINE auxfunc



    !=============================================
    ! Computation of clustering problem
    !=============================================
    SUBROUTINE clusterfunc(x,f)

        USE param, ONLY : zero, one, large
        USE initclust, ONLY : &
            maxdim, &      ! Maximum number of variables in optimization,
                           !   maxdim >= maxclust * maxnft.
            maxrec, &      ! maximum number of recorts in dataset.
            maxclust, &    ! maximum number of clusters.
            maxnft, &      ! maximum number of features in dataset.
            a, &           ! Data matrix.
            mf, &          ! Number of used features:
                           !   mf = nft, when no classes,
                           !   mf = nft - 1, when classes.
            nc, &          ! Current number of clusters, loops from 1 to nclust.
            nrecord, &     ! Number of reports in data.
            dminim, &      ! dminim(maxrec).
            tnorm          ! The total number of norms computed.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            x(maxdim), &
            f

        ! Local variables
        REAL(KIND=prec) ::  &
            f1,f2,f3
        INTEGER :: &
            i,j,k

        INTRINSIC ::  DMIN1

        f=zero
        DO i=1,nrecord
            f2=large
            DO k=1,nc
                tnorm=tnorm+one
                f3=zero
                DO j=1,mf
                    f3=f3+(a(i,j)-x(j+(k-1)*mf))**2
                END DO
                f2=DMIN1(f2,f3)
            END DO
            f=f+f2
        END DO

        RETURN

    END SUBROUTINE clusterfunc


    !==============================================================
    ! Computation of the gradient of the first component function
    !==============================================================
    SUBROUTINE fgrad1(x1,grad1)

        USE param, ONLY : zero, two
        USE initclust, ONLY : &
            maxdim, &      ! Maximum number of variables in optimization,
                           !   maxdim >= maxclust * maxnft.
            maxrec, &      ! maximum number of recorts in dataset.
            maxclust, &    ! maximum number of clusters.
            maxnft, &      ! maximum number of features in dataset.
            a, &           ! Data matrix.
            mf, &          ! Number of used features:
                           !   mf = nft, when no classes,
                           !   mf = nft - 1, when classes.
            nc, &          ! Current number of clusters, loops from 1 to nclust.
            nrecord, &     ! Number of reports in data.
            dminim, &      ! dminim(maxrec).
            ns, &          ! Switch for auxiliary and real clustering problem.
            m              ! Number of variables in optimization:
                           !   m = mf    if ns = 1,
                           !   m = mf*nc if ns = 2.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            x1(maxdim), &
            grad1(maxdim)

        ! Local variables
        INTEGER :: &
            i,j,k

        DO i=1,m
            grad1(i)=zero
        END DO

        IF (ns==1) THEN
            DO i=1,nrecord
                DO j=1,m
                    grad1(j)=grad1(j)+two*(x1(j)-a(i,j))
                END DO
            END DO
        END IF

        IF (ns==2) THEN
            DO i=1,nrecord
                DO j=1,nc
                    DO k=1,mf
                        grad1(k+(j-1)*mf)=grad1(k+(j-1)*mf)+two*(x1(k+(j-1)*mf)-a(i,k))
                    END DO
                END DO
            END DO
        END IF

        RETURN

    END SUBROUTINE fgrad1


    !==============================================================
    ! Computation of the gradient of the second component function
    !==============================================================
    SUBROUTINE fgrad2(x,grad2)

        USE param, ONLY : zero, one, two
        USE initclust, ONLY : &
            maxdim, &      ! Maximum number of variables in optimization,
                           !   maxdim >= maxclust * maxnft.
            maxrec, &      ! maximum number of recorts in dataset.
            maxclust, &    ! maximum number of clusters.
            maxnft, &      ! maximum number of features in dataset.
            a, &           ! Data matrix.
            mf, &          ! Number of used features:
                           !   mf = nft, when no classes,
                           !   mf = nft - 1, when classes.
            nc, &          ! Current number of clusters, loops from 1 to nclust.
            nrecord, &     ! Number of reports in data.
            dminim, &      ! dminim(maxrec).
            tnorm, &       ! The total number of norms computed.
            ns, &          ! Switch for auxiliary and real clustering problem.
            m              ! Number of variables in optimization:
                           !   m = mf    if ns = 1,
                           !   m = mf*nc if ns = 2.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            x(maxdim), &
            grad2(maxdim)

        ! Local variables
        REAL(KIND=prec) ::  &
            d2(maxclust), &
            d3,d4
        INTEGER :: &
            i,j,k,jindex

        DO i=1,m
            grad2(i)=zero
        END DO

        IF (ns==1) THEN
            DO i=1,nrecord
                tnorm=tnorm+one
                d3=zero
                DO j=1,mf
                    d3=d3+(a(i,j)-x(j))**2
                END DO
                IF (d3 > dminim(i)) THEN
                    DO j=1,m
                        grad2(j)=grad2(j)+two*(x(j)-a(i,j))
                    END DO
                END IF
            END DO
        END IF

        IF (ns == 2) THEN
            DO i=1,nrecord
                DO j=1,nc
                    tnorm=tnorm+one
                    d2(j)=zero
                    DO k=1,mf
                        d2(j)=d2(j)+(a(i,k)-x(k+(j-1)*mf))**2
                    END DO
                END DO
                d3=zero
                DO j=1,nc
                    d4=zero
                    DO k=1,nc
                        IF (k.ne.j) THEN
                            d4=d4+d2(k)
                        END IF
                    END DO
                    IF (d3<d4) THEN
                        d3=d4
                        jindex=j
                    END IF
                END DO
                DO j=1,nc
                    IF (j.ne.jindex) THEN
                        DO k=1,mf
                            grad2(k+(j-1)*mf)=grad2(k+(j-1)*mf)+two*(x(k+(j-1)*mf)-a(i,k))
                        END DO
                    END IF
                END DO
            END DO
        END IF

        RETURN

    END SUBROUTINE fgrad2

END MODULE functionmod
