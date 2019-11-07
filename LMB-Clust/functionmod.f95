!*************************************************************************
!*                                                                       *
!*     Computation of function and subgradients values of minimum        *
!*     sum of squares clustering problem (last modified 11.11.2017).     *
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
        auxfunc, &         ! S Computation of aux clustering problem.
        clusterfunc, &     ! S Computation of clustering problem.
        fgrad              ! S Computation of the subgradient of the (aux) clustering problem.
!    PRIVATE :: &


CONTAINS

    !=============================================
    ! Computation of aux clustering problem
    !=============================================
    SUBROUTINE auxfunc(x,fval)

        USE param, ONLY : zero, one, two, large
        USE initclust, ONLY : &
            maxdim, &      ! Maximum number of variables in optimization,
                           !   maxdim >= nclust * nft.
            a, &           ! Data matrix
            mf, &          ! Number of used features:
                           !   mf = nft, when no classes,
                           !   mf = nft - 1, when classes.
            nc, &          ! Current number of clusters.
            nrecord, &     ! Number of reports in data.
            dminim, &      ! dminim(nrecord), Distance of a(i,:) and the nearest centroid.
            xbest, &       ! xbest(nclust*nft), Best solution obtained.
            tnorm, &       ! The total number of norms computed.
            list1          ! list1(i) gives the cluster where point i belongs.


        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            x(maxdim), &
            fval

        ! Local variables
        REAL(KIND=prec) ::  &
            dprev(nc),          &
            f1,f10,f2,f3
        INTEGER :: &
            i,j,k,i1

        INTRINSIC ::  MIN
        
        DO k=1,nc-1
            tnorm=tnorm+one
            dprev(k)=zero
            DO j=1,mf
                dprev(k)=dprev(k)+(x(j)-xbest(j+(k-1)*mf))**2
            END DO
        END DO

        fval=zero
        loop: DO i=1,nrecord
            f1=dminim(i)
            f10=4.0_prec*f1
            i1=list1(i)
            IF(dprev(i1) >= f10) THEN
                fval=fval+f1
                CYCLE loop
            END IF
            tnorm=tnorm+one
            f2=zero
            DO j=1,mf
                f2=f2+(a(j,i)-x(j))**2
            END DO
            f3=MIN(f1,f2)
            fval=fval+f3
        END DO loop
        RETURN

    END SUBROUTINE auxfunc

    !=============================================
    ! Computation of clustering problem
    !=============================================
    SUBROUTINE clusterfunc(x,f)

        USE param, ONLY : zero, one, two, large
        USE initclust, ONLY : &
            maxdim, &      ! Maximum number of variables in optimization,
                           !   maxdim >= nclust * nft.
            a, &           ! Data matrix.
            mf, &          ! Number of used features:
                           !   mf = nft, when no classes,
                           !   mf = nft - 1, when classes.
            nc, &          ! Current number of clusters, loops from 1 to nclust.
            nrecord, &     ! Number of reports in data.
            tnorm          ! The total number of norms computed.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            x(maxdim), &
            dcent1(nc,nc), &
            f

        ! Local variables
        REAL(KIND=prec) ::  &
            f10,f2,f3
        INTEGER :: &
            i,j,k,kmin

        DO i=1,nc
            dcent1(i,i)=zero
        END DO

        DO i=1,nc
            DO j=i+1,nc
                dcent1(i,j)=zero
                DO k=1,mf
                    dcent1(i,j)=dcent1(i,j)+(x(k+(i-1)*mf)-x(k+(j-1)*mf))**2
                END DO
                dcent1(j,i)=dcent1(i,j)
            END DO
        END DO        

        f=zero
        kmin=1
        DO i=1,nrecord
            f2=large
            f10=large
            loop1: DO k=1,nc
                IF(k >= 2) THEN
                    IF(dcent1(k,kmin) >= f10) THEN
                        CYCLE loop1
                    END IF
                END IF                     
                tnorm=tnorm+one
                f3=zero
                DO j=1,mf
                    f3=f3+(a(j,i)-x(j+(k-1)*mf))**2
                END DO
                IF(f3 < f2) THEN
                    f2=f3
                    f10=4.0_prec*f2
                    kmin=k
                END IF
            END DO loop1
            f=f+f2
        END DO

        RETURN

    END SUBROUTINE clusterfunc

    !==============================================================
    ! Computation of the gradient of the first component function
    !==============================================================
    SUBROUTINE fgrad(x,grad)

        USE param, ONLY : zero, one, two, large
        USE initclust, ONLY : &
            maxdim, &      ! Maximum number of variables in optimization,
                           !   maxdim >= nclust * nft.
            a, &           ! Data matrix.
            mf, &          ! Number of used features:
                           !   mf = nft, when no classes,
                           !   mf = nft - 1, when classes.
            nc, &          ! Current number of clusters, loops from 1 to nclust.
            nrecord, &     ! Number of reports in data.
            dminim, &      ! dminim(nrecord), Distance of a(i) and the nearest centroid.
            xbest, &       ! xbest(nc*nft), Best solution obtained.
            list1, &       ! list1(i) gives the cluster where point i belongs.
            tnorm, &       ! The total number of norms computed.
            ns, &          ! Switch for auxiliary and real clustering problem.
            m              ! Number of variables in optimization:
                           !   m = mf    if ns = 1,
                           !   m = mf*nc if ns = 2.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            dcent1(nc,nc), &
            dprev(nc),         &
            x(maxdim), &
            grad(maxdim)

        ! Local variables
        REAL(KIND=prec) ::  &
            f1,f10,f12,f3
        INTEGER :: &
            i,j,k,jmin,i1

        DO i=1,m
            grad(i)=zero
        END DO

        IF (ns == 1) THEN
            DO k=1,nc-1
                tnorm=tnorm+one
                dprev(k)=zero
                DO j=1,mf
                    dprev(k)=dprev(k)+(x(j)-xbest(j+(k-1)*mf))**2
                END DO
            END DO
        
            loop1: DO i=1,nrecord
                f1=dminim(i)
                f10=4.0_prec*f1
                i1=list1(i)
                IF(dprev(i1) >= f10) THEN
                    CYCLE loop1
                END IF

                tnorm=tnorm+one
                f3=zero
                DO j=1,mf
                    f3=f3+(a(j,i)-x(j))**2
                END DO
                IF (f3 < dminim(i)) THEN
                    DO j=1,m
                        grad(j)=grad(j)+two*(x(j)-a(j,i))
                    END DO
                END IF
            END DO loop1

        ELSE
        
            DO i=1,nc
                DO j=i+1,nc
                    dcent1(i,j)=zero
                    DO k=1,mf
                        dcent1(i,j)=dcent1(i,j)+(x(k+(i-1)*mf)-x(k+(j-1)*mf))**2
                    END DO
                    dcent1(j,i)=dcent1(i,j)
                END DO
            END DO        

            jmin=1
            DO i=1,nrecord
                f10=large
                f12=large
                loop2: DO k=1,nc
                    IF(k >= 2) THEN
                        IF(dcent1(k,jmin) >= f10) THEN
                            CYCLE loop2
                        END IF
                    END IF
                    tnorm=tnorm+one
                    f3=zero
                    DO j=1,mf
                        f3=f3+(a(j,i)-x(j+(k-1)*mf))**2
                    END DO
                    IF (f12 > f3) THEN
                        f12=f3
                        f10=4.0_prec*f12
                        jmin=k
                    END IF
                END DO loop2
                DO j=1,mf
                    grad(j+(jmin-1)*mf)=grad(j+(jmin-1)*mf)+two*(x(j+(jmin-1)*mf)-a(j,i))
                END DO
            END DO
        END IF

    END SUBROUTINE fgrad

END MODULE functionmod
