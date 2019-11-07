!*************************************************************************
!*                                                                       *
!*     Computation of the value of the clustering problem and the        *
!*     correbonding subgradients. Called from subroutines lmbm, smdb,    *
!*     and dbundle (last modified 12.05.2016).                           *
!*                                                                       *
!*     This file is specially modified to solve clustering problems.     *
!*                                                                       *
!*************************************************************************
!*
!*     Modules included:
!*
!*     obj_fun         !
!*

MODULE obj_fun  ! Computation of the value and the subgradient of the 
                  ! objective function.

    USE r_precision, ONLY : prec  ! Precision for reals.
    IMPLICIT NONE

    PUBLIC :: &
        myf,    &   ! Computation of the value of the objective.
        myg         ! Computation of the subgradient of the objective.

CONTAINS
    !************************************************************************
    !*                                                                      *
    !*     * SUBROUTINE myf *                                               *
    !*                                                                      *
    !*     Computation of the value of the objective.                       *
    !*                                                                      *
    !************************************************************************
     
    SUBROUTINE myf(n,x,f,iterm)

        USE param, ONLY : large                   ! Parameters.
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
    !*     * SUBROUTINE myg *                                               *
    !*                                                                      *
    !*     Computation of the subgradient of the objective function.        *
    !*                                                                      *
    !************************************************************************
     
    SUBROUTINE myg(n,x,g,iterm)

        USE functionmod, ONLY : &
            fgrad                ! S Computation of the gradient of the
                                 !   clustering function

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
        CALL fgrad(x,g)

        RETURN

    END SUBROUTINE myg

END MODULE obj_fun
