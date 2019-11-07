!*************************************************************************
!*                                                                       *
!*     Subprograms for LMBM, SMDB and D-Bundle.                          *
!*     (last modified 11.05.2016).                                       *
!*                                                                       *
!*************************************************************************
!*
!*     Modules included:
!*
!*     subpro         !
!*

MODULE subpro  ! Subprograms for lmbm, ldgm and d-bundle

  USE r_precision, ONLY : prec      ! Precision for reals.
  IMPLICIT NONE

! MODULE subpro includes the following subroutines (S) and functions (F).
  PUBLIC :: &
       vdot, &   ! F Dot product of two vectors.
       copy, &   ! S Copying a vector.
       copy2, &  ! S Copying two vectors.
       xdiffy, & ! S Difference of two vectors z:= x - y.
       xsumy, &  ! S Sum of two vectors z:= x + y.
       scdiff, & ! S Difference of the scaled vector and a vector z:= a*x - y.
       scsum, &  ! S Sum of a vector and the scaled vector z:= y + a*x.
       vxdiag, & ! S Vector is multiplied by a diagonal matrix y:=d*x.
       symax, &  ! S Multiplication of a dense symmetric matrix A by a vector x.
       cwmaxv, & ! S Multiplication of a columnwise stored dense rectangular matrix by a vector.
       rwaxv2, & ! S Multiplication of two rowwise stored dense rectangular  
                 !   matrices A and B by vectors X and Y.
       trlieq, & ! S Solving x from linear equation u*x=y or u'*x=y, 
                 !   where u is an upper triangular matrix.
       lineq, &  ! S Solver from linear equation.
       calq      ! S Solving x from linear equation A*x=y. Contains:
                 !     S mxdpgf   Gill-Murray decomposition of a dense symmetric matrix.

CONTAINS

  FUNCTION vdot(n,x,y) RESULT(xty) ! Dot product of two vectors.
    USE param, ONLY : zero
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
         x,y         ! Input vectors.

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n           ! Vectors dimension.

! Local Scalars
    REAL(KIND=prec) xty
    INTEGER :: i

    xty = zero
    DO i = 1,n
       xty = xty + x(i)*y(i)
    END DO

  END FUNCTION vdot


  SUBROUTINE scalex(n,a,x,y)  ! Scaling a vector y:= a*x.
    IMPLICIT NONE
      
! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
         x           ! Input vector.
    REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
         y           ! Output vector y:= a*x.

! Scalar Arguments
    REAL(KIND=prec), INTENT(IN) :: &
         a           ! Scaling parameter.
    INTEGER, INTENT(IN) :: &
         n           ! Vectors dimension.

! Local Scalars
    INTEGER :: i

    DO i = 1,n
       y(i) = a*x(i)
    END DO
      
  END SUBROUTINE scalex

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

  SUBROUTINE xsumy(n,x,y,z)  ! Sum of two vectors z:= x + y.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
         x,y         ! Input vectors.
    REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
         z           ! Output vector z:= x + y.

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n           ! Vectors dimension.
      
! Local Scalars
    INTEGER :: i

    DO  i = 1,n
       z(i) = x(i) + y(i)
    END DO

  END SUBROUTINE xsumy

  SUBROUTINE scdiff(n,a,x,y,z)  ! Difference of the scaled vector and a vector z:= a*x - y.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
         x,y         ! Input vectors.
    REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
         z           ! Output vector z:= a*x - y.

! Scalar Arguments
    REAL(KIND=prec), INTENT(IN) :: &
         a           ! Scaling factor.
    INTEGER, INTENT(IN) :: &
         n           ! Vectors dimension.
      
! Local Scalars
    INTEGER :: i

    DO  i = 1,n
       z(i) = a*x(i) - y(i)
    END DO
 
  END SUBROUTINE scdiff

  SUBROUTINE scsum(n,a,x,y,z)  ! Sum of a vector and the scaled vector z:= y + a*x.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
         x,y         ! Input vectors.
    REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
         z           ! Output vector z:= a*x + y.

! Scalar Arguments
    REAL(KIND=prec), INTENT(IN) :: &
         a           ! Scaling factor.
    INTEGER, INTENT(IN) :: &
         n           ! Vectors dimension.
      
! Local Scalars
    INTEGER :: i

    DO  i = 1,n
       z(i) = a*x(i) + y(i)
    END DO

  END SUBROUTINE scsum

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

! Local Scalars
    INTEGER :: i

    DO i = 1,n
       y(i) = x(i)
    END DO
      
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

  SUBROUTINE symax(n,m,iold,a,x,y)  ! Multiplication of a dense symmetric matrix A by a vector x.
    USE param, ONLY : zero
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
         x, &        ! Input vector stored in a circular order.
         a           ! Dense symmetric matrix stored in the packed form: a(n*(n+1)/2).
    REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
         y           ! Output vector y:= a*x. Vector y has the same circular order than x.

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n, &         ! Order of matrix A.
         m, &         ! Length of vector x, m >= n, note that only n
                      ! components from vector x are used.
         iold         ! Index, which controlls the circular order of
                      ! the vector x.

! Local Scalars
    INTEGER :: i,j,k,l

    DO j=1,n
       l=j+iold-1
       IF (l > m) l=l-m
       y(l) = zero
       k=l
       DO i=j,n
          y(l) = a((i-1)*i/2+j)*x(k)+y(l)
          k=k+1
          IF (k > m) k=k-m
       END DO
    END DO

    DO j=2,n
       l=j+iold-1
       IF (l > m) l=l-m
       k=iold
       DO i=1,j-1
          IF (k > m) k=k-m
          y(l) = a((j-1)*j/2+i)*x(k)+y(l)
          k=k+1
       END DO
    END DO
      
  END SUBROUTINE symax

  SUBROUTINE cwmaxv(n,m,a,x,y)  ! Multiplication of a columnwise stored dense 
                                ! rectangular matrix A by a vector x.
    USE param, ONLY : zero
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
         x, &        ! Input vector (dimension m).
         a           ! Rectangular matrix stored columnwise in the
                     ! one-dimensional array (dimension n*m).
    REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
         y           ! Output vector equal to s*a*x. If m = 0 y is a zero vector. 

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n, &        ! Number of rows of the matrix A.
         m           ! Number of columns of the matrix A.

! Local Scalars
    INTEGER :: i,j,k
      
    DO i = 1,n
       y(i) = zero
    END DO

    k = 1
    DO j = 1,m
       CALL scsum(n,x(j),a(k:),y,y)
       k = k + n
    END DO

  END SUBROUTINE cwmaxv

  SUBROUTINE rwaxv2(n,m,a,b,x,y,v,w)  ! Multiplication of two rowwise stored dense rectangular  
                                      ! matrices A and B by vectors X and Y.
    USE param, ONLY : zero
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
         x,y, &      ! Input vectors (dimension n).
         a,b         ! Rectangular matrices stored rowwise in the
                     ! one-dimensional array (dimension n*m).
    REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
         v,w         ! Output vectors v=a*x and w=b*y. 

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n, &        ! Number of columns of the matrices A and B.
         m           ! Number of rows of the matrices A and B.

! Local Scalars
    REAL(KIND=prec) :: tmp1,tmp2
    INTEGER :: i,j,k
      
    k = 0
    DO i = 1,m
       tmp1 = zero
       tmp2 = zero
       DO j = 1,n
          tmp1 = tmp1 + a(k+j)*x(j)
          tmp2 = tmp2 + b(k+j)*y(j)
       END DO
       v(i) = tmp1
       w(i) = tmp2
       k = k + n
    END DO

  END SUBROUTINE rwaxv2


  SUBROUTINE trlieq(n,m,iold,u,x,y,job,ierr)  ! Solving x from linear equation u*x=y or u'*x=y, 
                                              ! where u is an upper triangular matrix.
    USE param, ONLY : small
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
         y, &        ! Input vector stored in a circular order.
         u           ! Triangular matrix.
    REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
         x           ! Output vector y:= a*x. Vector y has the same circular order than x.
                     ! Note that x may be equal to y in calling sequence.

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n, &         ! Order of matrix U.
         m, &         ! Length of vectors x and y, m >= n, note that only n
                      ! components from vectors are used.
         iold, &      ! Index, which controlls the circular order of
                      ! the vectors x and y.
         job          ! Option:
                      !   0  - x:=(u')**(-1)*y, u upper triangular.
                      !   1  - x:=u**(-1)*y, u upper triangular.
    INTEGER, INTENT(OUT) :: &
         ierr         ! Error indicador: 
                      !   0   - Everything is ok.
                      !  -3   - Error; 0 at diagonal.

! Local Scalars
    INTEGER :: i,ii,ij,j,k,l,ji

! Intrinsic Functions
    INTRINSIC ABS
      
    ierr = -3
      
    DO i=1,m
       x(i)=y(i)
    END DO
      
    IF (job == 0) THEN
     
! x=u'**(-1)*y, u' = [u1         ] is lower triangular.
!                    [u2 u3      ]
!                    [u4 u5 u6   ]
!                    [.  .  .  . ]
         
       ii = 0
       DO  i = 1,n
          ii=ii+i
          l=i+iold-1
          IF (l > m) l=l-m
          IF (ABS(u(ii)) <= small) RETURN
          x(l) = x(l)/u(ii)
          DO j = i+1,n
             ji = (j-1)*j/2+i
             k=j+iold-1
             IF (k > m) k=k-m
             x(k) = x(k) - u(ji)*x(l)
          END DO
       END DO
             
         
    ELSE IF (job == 1) THEN
     
! x=u**(-1)*y, u = [u1 u2 u4 . ] is upper triangular.
!                  [   u3 u5 . ]
!                  [      u6 . ]
!                  [         . ]
         
       ii = n* (n+1)/2
       DO i = n,1,-1
          l=i+iold-1
          IF (l > m) l=l-m
          IF (ABS(u(ii)) <= small) RETURN
          ij = ii
          DO j = i + 1,n
             k=j+iold-1
             IF (k > m) k=k-m
             ij = ij + j - 1
             x(l) = x(l) - u(ij)*x(k)
          END DO
          x(l)=x(l)/u(ii)
          ii = ii - i
       END DO
         
         
    ELSE
         
       RETURN
    END IF
      
    ierr = 0

  END SUBROUTINE trlieq
      
  SUBROUTINE lineq(n,m,iold,a,x,y,ierr)  ! Solving X from linear equation A*X=Y. 
                                         ! Positive definite matrix A+E is given using 
                                         ! the factorization A+E=L*D*L' obtained by the
                                         ! subroutine mxdpgf.
    USE param, ONLY : small
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
         y, &        ! Input vector stored in a circular order (dimension m).
         a           ! Factorization a+e=l*d*l' obtained by the subroutine mxdpgf.
    REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
         x           ! Output vector y:= a*x. Vector x has the same circular order than y.
                     ! Note that x may be equal to y in calling sequence.

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n, &        ! Order of matrix a.
         m, &        ! Length of vectors x and y, m >= n, note that only n
                     ! components from vectors are used.
         iold        ! Index, which controlls the circular order of
                     ! the vectors x and y.
    INTEGER, INTENT(OUT) :: &
         ierr        ! Error indicador: 
                     !   0   - Everything is ok.
                     !  -2   - Error; indefinite matrix.

! Local Scalars
    INTEGER :: i,ii,ij,j,k,l


    ierr = -2
      
! Phase 1: x=l**(-1)*x

    ij = 0
    DO i = 1,n
       l=i+iold-1
       IF (l > m) l=l-m
       x(l) = y(l)
         
       DO j = 1,i - 1
          ij = ij + 1
          k=j+iold-1
          IF (k > m) k=k-m
          x(l) = x(l) - a(ij)*x(k)
       END DO
       ij = ij + 1
    END DO

! Phase 2 : x:=d**(-1)*x

    ii = 0
    DO i = 1,n
       ii = ii + i
       IF (a(ii) <= small) RETURN
       l=i+iold-1
       IF (l > m) l=l-m
       x(l) = x(l)/a(ii)
    END DO

! Phase 3 : x:=trans(l)**(-1)*x

    ii = n* (n-1)/2
    DO i = n - 1,1,-1
       ij = ii
       l=i+iold-1
       IF (l > m) l=l-m
       DO j = i + 1,n
          k=j+iold-1
          IF (k > m) k=k-m
          ij = ij + j - 1
          x(l) = x(l) - a(ij)*x(k)
       END DO
       ii = ii - i
    END DO

    ierr = 0

  END SUBROUTINE lineq

  SUBROUTINE calq(n,m,iold,a,x,y,iprint)  ! Solving x from linear equation A*x=y.
    USE param, ONLY : zero,small,one
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
         y           ! Input vector stored in a circular order (dimension m).
    REAL(KIND=prec), DIMENSION(:), INTENT(INOUT) :: &
         a           ! On input: Dense symmetric matrix stored in the packed form. 
                     ! On output: factorization A+E=L*D*trans(L).
    REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
         x           ! Output vector y:= a*x. Vector x has the same circular order than y.
                     ! Note that x may be equal to y in calling sequence.

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n, &        ! Order of matrix a.
         m, &        ! Length of vectors x and y, m >= n, note that only n
                     ! components from vectors are used.
         iold, &     ! Index, which controlls the circular order of
                     ! the vectors x and y.
         iprint      ! Printout specification.
!    INTEGER, INTENT(OUT) :: &
!         ierr        ! Error indicador: 
!                     !   0   - Everything is ok.
!                     !  -2   - Error; indefinite matrix.

! Local Scalars
    REAL(KIND=prec) :: eta,bet
    INTEGER :: inf,ierr

      
    eta = small+small
      
    CALL mxdpgf(n,a,inf,eta,bet)


    IF (iprint == 2) THEN
       IF (inf < 0) THEN
          WRITE (6,FMT='(1X,''Warning: Insufficiently positive'' &
               '' definite matrix detected. '')')
          WRITE (6,FMT='(1X,''Correction added.'')')
         
       ELSE IF (inf > 0) THEN
          WRITE (6,FMT='(1X,''Warning: Indefinite'' &
            '' matrix detected. Correction added.'')')
       END IF
    END IF
      
    CALL lineq(n,m,iold,a,x,y,ierr)
    IF (ierr /= 0) THEN
       PRINT*,'hihu'
       IF (iprint == 2) THEN
          WRITE (6,FMT='(1X,''Warning: Indefinite matrix detected. '')')
       END IF
    END IF

  CONTAINS

    SUBROUTINE mxdpgf(n,a,inf,alf,tau)  ! Factorization A+E=L*D*trans(L) of a dense symmetric positive
                                        ! definite matrix A+E, where D and E are diagonal positive 
                                        ! definite matrices and L is a lower triangular matrix. 
                                        ! If A is sufficiently positive definite then E=0.
      
! Array Arguments
      REAL(KIND=prec), DIMENSION(:), INTENT(INOUT) :: &
           a         ! On input: Dense symmetric matrix stored in the packed form. 
                     ! On output: factorization A+E=L*D*trans(L).

! Scalar Arguments
      REAL(KIND=prec), INTENT(INOUT) :: &
           alf       ! On input a desired tolerance for positive definiteness. 
                     ! On output the most negative diagonal element used in the factorization
                     ! process (if inf>0).
      REAL(KIND=prec), INTENT(OUT) :: &
           tau       ! Maximum diagonal element of matrix E.

      INTEGER, INTENT(IN) :: &
           n         ! Order of matrix a.
      INTEGER, INTENT(OUT) :: &
           inf       ! An information obtained in the factorization process:
                     !    inf=0  - A is sufficiently positive definite and E=0. 
                     !    inf<0  - A is not sufficiently positive definite and E>0.
                     !    inf>0  - A is indefinite and inf is an index of the most negative 
                     !             diagonal element used in the factorization process.

! Local Scalars
      REAL(KIND=prec) :: bet,del,gam,rho,sig,tol
      INTEGER :: i,ij,ik,j,k,kj,kk,l

! Intrinsic Functions
      INTRINSIC ABS,MAX

      l = 0
      inf = 0
      tol = alf
      

! Estimation of the matrix norm

      alf = zero
      bet = zero
      gam = zero
      tau = zero
      kk = 0

      DO k = 1,n
         kk = kk + k
         bet = MAX(bet,ABS(a(kk)))
         kj = kk
         DO j = k + 1,n
            kj = kj + j - 1
            gam = MAX(gam,ABS(a(kj)))
         END DO
      END DO
      bet = MAX(tol,bet,gam/n)

      del = tol*MAX(bet,one)
      kk = 0
      DO k = 1,n
         kk = kk + k

!     Determination of a diagonal correction

         sig = a(kk)
         IF (alf > sig) THEN
            alf = sig
            l = k
         END IF

         gam = zero
         kj = kk
         DO j = k + 1,n
            kj = kj + j - 1
            gam = MAX(gam,ABS(a(kj)))
         END DO
         gam = gam*gam
         rho = MAX(ABS(sig),gam/bet,del)
         IF (tau < rho-sig) THEN
            tau = rho - sig
            inf = -1
         END IF
         
! Gaussian elimination

         a(kk) = rho
         kj = kk
         DO j = k + 1,n
            kj = kj + j - 1
            gam = a(kj)
            a(kj) = gam/rho
            ik = kk
            ij = kj
            DO i = k + 1,j
               ik = ik + i - 1
               ij = ij + 1
               a(ij) = a(ij) - a(ik)*gam
            END DO
         END DO
      END DO
      IF (l > 0 .AND. ABS(alf) > del) inf = l
      
    END SUBROUTINE mxdpgf
  END SUBROUTINE calq

END MODULE subpro
