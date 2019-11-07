!*************************************************************************
!*                                                                       *
!*     Computation of the value f(x) = f1(x) - f2(x) and the             *
!*     correbonding (sub)gradients of the component functions f1 and     *
!*     f2. Called from subroutine dcdb (last modified 24.02.2016).       *
!*                                                                       *
!*     This file is specially modified to solve clustering problems.     *
!*                                                                       *
!*************************************************************************
!*
!*     Modules included:
!*
!*     dc_fun         !
!*

MODULE dc_fun  ! Computation of the value f(x) =f1(x) - f2(x) and the
               ! correbonding (sub)gradients of the component functions
               ! f1 and f2.

    USE r_precision, ONLY : prec  ! Precision for reals.
    IMPLICIT NONE

    ! MODULE dc_fun includes the following subroutines (S).

    PUBLIC :: &
        myf,     &   ! S  Computation of the value of the objective.
        myg1,    &   ! S  Computation of the gradient of the first component.
        myg2         ! S  Computation of the subgradient of the second component.

CONTAINS

    !************************************************************************
    !*                                                                      *
    !*     * SUBROUTINE myf *                                               *
    !*                                                                      *
    !*     Computation of the value of the objective f = f1 - f2.           *
    !*                                                                      *
    !************************************************************************
     
    SUBROUTINE myf(n,x,f,iterm)

        USE param, ONLY : large,zero,one          ! Parameters.
        USE initclust, ONLY : ns                  ! Switch for auxiliary and real clustering problem.
        USE functionmod, ONLY : &
            auxfunc, &     ! S Computation of auxiliary clustering problem
            clusterfunc    ! S Computation of clustering problem


        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(n), INTENT(IN) :: &
            x  ! Vector of variables.

        ! Scalar Arguments
        REAL(KIND=prec), INTENT(OUT) :: f  ! Value of the function.
        INTEGER, INTENT(IN) :: n           ! Number of variables.
        INTEGER, INTENT(OUT) :: iterm      ! Cause of termination:
                                           !   0  - Everything is ok.
                                           !  -3  - Failure in function calculations
                                           !        (assigned by the user).

        iterm = 0

        ! Function evaluation
        IF (ns == 1) THEN  ! Auxiliary problem
            CALL auxfunc(x,f)

        ELSE  ! Clustering problem
            CALL clusterfunc(x,f)
        END IF

        ! Error checking.
        !    IF (f > large) iterm = -3  !
        !    IF (f < -large) iterm = -3 !
        RETURN
      
    END SUBROUTINE myf


    !************************************************************************
    !*                                                                      *
    !*     * SUBROUTINE myg1 *                                              *
    !*                                                                      *
    !*     Computation of the gradient of the first DC-component f1.        *
    !*                                                                      *
    !************************************************************************
     
    SUBROUTINE myg1(n,x,g,iterm)

        USE param, ONLY : zero   ! Parameters.
        USE functionmod, ONLY : &
            fgrad1               ! S Computation of the gradient of the
                                 !   first component function

        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(n), INTENT(IN) :: x  ! Vector of variables.
        REAL(KIND=prec), DIMENSION(n), INTENT(OUT) :: g ! Subgradient.

        ! Scalar Arguments
        INTEGER, INTENT(IN) :: n                        ! Number of variables.
        INTEGER, INTENT(OUT) :: iterm                   ! Cause of termination:
                                                        !   0  - Everything is ok.
                                                        !  -3  - Failure in subgradient calculations
                                                        !        (assigned by the user).

        iterm = 0

        ! Gradient evaluation.
        CALL fgrad1(x,g)

        RETURN

    END SUBROUTINE myg1


    !************************************************************************
    !*                                                                      *
    !*     * SUBROUTINE myg2 *                                              *
    !*                                                                      *
    !*     Computation of the subgradient of the second DC-component f2.    *
    !*                                                                      *
    !************************************************************************

    SUBROUTINE myg2(n,x,g,iterm)

        USE param, ONLY : zero, one  ! Parameters.
        USE functionmod, ONLY : &
            fgrad2                   ! S Computation of the subgradient of
                                     !   the second component function

        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(n), INTENT(IN) :: x  ! Vector of variables.
        REAL(KIND=prec), DIMENSION(n), INTENT(OUT) :: g ! Subgradient.

        ! Scalar Arguments
        INTEGER, INTENT(IN) :: n            ! Number of variables.
        INTEGER, INTENT(OUT) :: iterm       ! Cause of termination:
                                            !   0  - Everything is ok.
                                            !  -3  - Failure in subgradient calculations
                                            !        (assigned by the user).


        iterm = 0

        ! Subgradient evaluation.
        CALL fgrad2(x,g)

        RETURN

    END SUBROUTINE myg2

END MODULE dc_fun

