!*************************************************************************
!*                                                                       *
!*     Subprograms for DCDB. Called from subroutine dcdb                 *
!*     (last modified 26.02.2016).                                       *
!*                                                                       *
!*************************************************************************
!*
!*     Modules included:
!*
!*     subpro         !
!*

MODULE subpro              ! Subprograms for dcdb.

    USE r_precision, ONLY : prec      ! Precision for reals.
    IMPLICIT NONE

    ! MODULE subpro includes the following subroutines (S) and functions (F).
    PUBLIC :: &
        copy, &   ! S Copying a vector.
        copy2, &  ! S Copying two vectors.
        xdiffy, & ! S Difference of two vectors z:= x - y.
        vxdiag    ! S Vector is multiplied by a diagonal matrix y:=d*x.

CONTAINS

    SUBROUTINE xdiffy(n,x,y,z)  ! Difference of two vectors z:= x - y.
        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
            x,y         ! Input vectors.
        REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
            z           ! Output vector z:= x - y.

        ! Scalar Arguments
        INTEGER, INTENT(IN) :: &
            n           ! Vectors dimension.
      
        ! Local Scalars
        INTEGER :: i

        DO  i = 1,n
            z(i) = x(i) - y(i)
        END DO
 
    END SUBROUTINE xdiffy

    SUBROUTINE copy(n,x,y)  ! Copying of vector Y:= X.
        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
            x           ! Input vector.
        REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
            y           ! Output vector.

        ! Scalar Arguments
        INTEGER, INTENT(IN) :: &
            n           ! Vectors dimension.

        y = x
      
    END SUBROUTINE copy

    SUBROUTINE copy2(n,x,y,z,v)  ! Copying of two vectors: y:=x, v:=z.
        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
            x,z         ! Input vectors.
        REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
            y,v         ! Output vectors.

        ! Scalar Arguments
        INTEGER, INTENT(IN) :: &
            n           ! Vectors dimension.

        ! Local Scalars
        INTEGER :: i

        DO i = 1,n
            y(i) = x(i)
            v(i) = z(i)
        END DO

    END SUBROUTINE copy2

    SUBROUTINE vxdiag(n,d,x,y)  ! Vector is multiplied by a diagonal matrix y:=d*x.
        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
            x, &        ! Input vector.
            d           ! Diagonal matrix stored as a vector with n elements.
        REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
            y           ! Output vector y:= d*x.

        ! Scalar Arguments
        INTEGER, INTENT(IN) :: &
            n           ! Vectors dimension.
      
        ! Local Scalars
        INTEGER :: i

        DO  i = 1,n
            y(i) = x(i)*d(i)
        END DO
 
    END SUBROUTINE vxdiag



END MODULE subpro
