!*************************************************************************
!*                                                                       *
!*     DCDB - DCDB method for nonsmooth DC-programming                   *
!*                                                                       *
!*     by Napsu Karmitsa 2016 (last modified 2016)                       *
!*                                                                       *
!*     DCDB is developed to solve problems given in the form             *
!*                                                                       *
!*        f = f1 - f2                                                    *
!*                                                                       *
!*     where f1 and f2 are convex and f1 is a smooth function.           *
!*     This file is specially modified to solve clustering problems.     *
!*                                                                       *
!*************************************************************************
!*
!*     Modules included:
!*
!*     dcdb_mod            ! DCDB Method
!*

MODULE dcdb_mod            ! DCDB Method

    USE r_precision, ONLY : prec  ! Precision for reals.
    IMPLICIT NONE

    ! MODULE dcdb_mod includes the following subroutines (S) and functions (F).

    PUBLIC :: &
        optim2             ! S Initializing DCDB subroutine for clustering.
    PRIVATE :: &
        dcdb, &            ! S DCDB subroutine for nonsmooth large-scale
                           !   optimization. Contains:
                           !     S restar     Initialization and reinitialization.
                           !     S armijo     Armijo line search.
                           !     S nmarmijo   Nonmonotone Armijo line search.
                           !     S aggre1     Simplified subgradient aggregation.
                           !     S aggre2     Subgradient aggregation.
        wprint, &          ! S Printout the error and warning messages.
        rprint             ! S Printout the results.

CONTAINS

    ! Optimization with DCDB
    SUBROUTINE optim2(x0,x1,fvalue)

        USE initclust, ONLY : &
            maxdim, &      ! Maximum number of variables.
            maxrec, &      ! Maximum number of reports in dataset.
            maxclust, &    ! Maximum number of clusters.
            mf, &          ! Number of used features:
                           !   mf = nft, when no classes,
                           !   mf = nft - 1, when classes.
            nc, &          ! Current number of clusters.
            ns, &          ! Switch for auxiliary and real clustering problem.
            m              ! Number of variables for function calculations.
        USE initdcdb, ONLY : &
            n, &           ! Number of variables.
            x, &           ! Vector of variables.
            mcu, &         ! Maximum number of stored corrections.
            epsl, &        ! Line search parameter.
            defaults, &    ! S  Default values for parameters.
            init_par       ! S  Further initialization of parameters.

        REAL(KIND=prec) ::  &
            x0(maxdim),     &
            x1(maxdim),     &
            fvalue

        INTEGER, DIMENSION(4) :: &
            iout           ! Output integer parameters.
                           !   iout(1)   Number of used iterations.
                           !   iout(2)   Number of used function evaluations.
                           !   iout(3)   Number of used subgradient evaluations
                           !   iout(4)   Cause of termination:
                           !               1  - The problem has been solved
                           !                    with desired accuracy.
                           !               2  - Changes in function values < tolf in mtesf
                           !                    subsequent iterations.
                           !               3  - Changes in function value <
                           !                    tolf*small*MAX(|f_k|,|f_(k-1)|,1),
                           !                    where small is the smallest positive number
                           !                    such that 1.0 + small > 1.0.
                           !               4  - Number of function calls > mfe.
                           !               5  - Number of iterations > mit.
                           !               6  - Time limit exceeded.
                           !               7  - f < tolb.
                           !               8  - Failure in attaining the demanded accuracy.
                           !              -1  - Two consecutive restarts.
                           !              -2  - Number of restarts > maximum number
                           !                    of restarts.
                           !              -3  - Failure in function or subgradient
                           !                    calculations (assigned by the user).
                           !              -4  -
                           !              -5  - Invalid input parameters.
                           !              -6  - Unspecified error.
        INTEGER :: i

        ! number of variables
        IF(ns.EQ.1) n = mf
        IF(ns.EQ.2) n = nc*mf
        m=n

        DO i = 1, n
            x(i) = x0(i)
        END DO

        CALL defaults()    !
        CALL init_par()    !

        IF (n <= 0) PRINT*,'n<0'
        IF (n > maxdim) PRINT*,'n is greater than maxdim.'
        IF (epsl >= 0.25_prec) PRINT*,'epsl >= 0.25'
        IF (mcu <= 0) PRINT*,'mcu <= 0'

        CALL dcdb(fvalue,iout(1),iout(2),iout(3),iout(4))

        DO i=1, n
            x1(i) = x(i)
        END DO

        RETURN

    END SUBROUTINE optim2


    !***********************************************************************
    !*                                                                     *
    !*     * SUBROUTINE dcdb *                                             *
    !*                                                                     *
    !*     Diagonal bundle subroutine for nonsmooth DC-optimization.       *
    !*                                                                     *
    !***********************************************************************
      
    SUBROUTINE dcdb(f,nit,nfe,nge,iterm)
        USE param, ONLY : small,large,zero,half,one  ! Parameters.
        USE initdcdb, ONLY : &
            n, &           ! Number of variables.
            x, &           ! Vector of variables
            mcu, &         ! Maximum number of stored corrections, mcu >= 1.
            mnma, &        ! Maximum number of function values used in nonmonotone line search.
            inma, &        ! Selection of line search method:
                           !   nma = 0, Armijo line search,
                           !   nma = 1, nonmonotone Armijo line search.
            iprint, &      ! Printout specification (see initializat for more details).
            tolf, &        ! Tolerance for change of function values.
            tolf2, &       ! Second tolerance for change of function values.
            tolb, &        ! Tolerance for the function value.
            tolg, &        ! Tolerance for the first termination criterion.
            told, &        ! Upper bound for the diagonal elements of diag.
            mintold, &     ! Lower bound for the diagonal elements of diag.
            mintold2, &    ! Lower bound for the "concave" update.
            eta, &         ! Distance measure parameter, eta >= 0.
            epsl, &        ! Line search parameter, 0 < epsl < 0.25.
            mtesf, &       ! Maximum number of iterations with changes of
                           ! function values smaller than tolf.
            mit, &         ! Maximun number of iterations.
            mfe            ! Maximun number of function evaluations.
        USE dc_fun, ONLY : &
            myf, &         ! Computation of the value of the objective f=f1-f2.
            myg1, &        ! Computation of the gradient of the first DC-component.
            myg2           ! Computation of the subgradient of the second DC-component.
        USE subpro, ONLY : &
            vxdiag, &      ! Vector is multiplied by a diagonal matrix y:=diag*x.
            xdiffy, &      ! Difference of two vectors.
            copy, &        ! Copying of a vector.
            copy2          ! Copying of two vectors.


        IMPLICIT NONE

        ! Scalar Arguments
        REAL(KIND=prec), INTENT(OUT) :: &
            f              ! Value of the objective function.
        INTEGER, INTENT(OUT) :: &
            nit,nfe,nge    ! Number of iterations, and function and subgradient evaluations.
        INTEGER, INTENT(OUT) :: &
            iterm          ! Cause of termination:
                           !   1  - The problem has been solved with desired accuracy.
                           !   2  - Changes in function values < tolf in mtesf
                           !        subsequent iterations.
                           !   3  - Changes in function value < tolf*small*MAX(|f_k|,|f_(k-1)|,1),
                           !        where small is the smallest positive number such that
                           !        1.0 + small > 1.0.
                           !   4  - Number of function calls > mfe.
                           !   5  - Number of iterations > mit.
                           !   7  - f < tolb.
                           !   8  - Failure in attaining the demanded accuracy.
                           !  -1  - Two consecutive restarts.
                           !  -2  - Number of restarts > maximum number of restarts.
                           !  -3  - Failure in function or subgradient calculations
                           !        (assigned by the user).
                           !  -4
                           !  -5  - Invalid input parameters.
                           !  -6  - Unspecified error.


        ! Local Arrays
        REAL(KIND=prec), DIMENSION(n) :: &
            g, &           ! Subgradient of the objective function.
            gp, &          ! Previous subgradient of the objective function.
            ga, &          ! Aggregate subgradient.
            g1,g2, &       ! Subgradients of component functions.
            gp1,gp2, &     ! Previous subgradients of component functions.
            xo, &          ! Previous vector of variables.
            s, &           ! Difference of current and previous variables.
            u1, &          ! Difference of current and previous gradients of f1.
            u2, &          ! Difference of current and previous subgradients of f2.
            u, &           ! temporary variable to make things work
            d, &           ! Direction vector.
            diag, &        ! Diagonal matrix.
            cdiag, &       ! "Convex" diagonal matrix.
            ncdiag, &      ! "Concave" diagonal matrix.
            diagtmp        ! Auxiliary diagonal matrix.
        REAL(KIND=prec), DIMENSION(n*mcu) :: &
            ssmatrix, &    ! Matrix whose columns are stored differences of
                           ! variables^2.
            csumatrix, &   ! Matrix whose columns are stored differences of
                           ! variables and subgradients when alfn >= 0.
            ncsumatrix     ! Matrix whose columns are stored differences of
                           ! variables and subgradients when alfn < 0.
        REAL(KIND=prec), DIMENSION(mnma) :: &
            fold           ! Old function values.


        ! Local Scalars
        REAL(KIND=prec) :: &
            f1,f2, &       ! values of component functions. Not needed
            alfn, &        ! Linearization error.
            alfv, &        ! Aggregate linearization error.
            epsr, &        ! Line search parameter.
            dnorm, &       ! Euclidean norm of the direction vector.
            gnorm, &       ! Euclidean norm of the aggregate subgradient.
            xnorm, &       ! Stopping criterion.
            pxnorm, &      ! Previous stopping criterion.
            p, &           ! Directional derivative.
            t, &           ! Stepsize.
            fo, &          ! Previous value of the objective.
            stu, &         ! stu = trans(s)*u.
            utu, &         ! utu = trans(u)*u.
            sts, &         ! sts = trans(s)*s.
            ptmp           ! Scalar for convex combination.
        INTEGER :: i,j, &
            mccc, &        ! Current number of stored "convex" corrections.
            mccnc, &       ! Current number of stored "concave" corrections.
            inewc, &       ! Index for the "convex" circular arrays.
            inewnc, &      ! Index for the "concave" circular arrays.
            inewnma, &     ! Index for the circular arrays containing function values.
            iters, &       ! Null step indicator.
                           !   0  - Null step.
                           !   1  - Serious step.
            nnk, &         ! Consecutive null steps counter.
            nser, &        ! Serious steps counter.
            cnull, &       ! Convex null steps counter.
            ncnull, &      ! Concave null steps counter.
            neps, &        ! Number of consecutive equal stopping criterions.
            ntesf, &       ! Number of tests on function decrease.
            ncres, &       ! Number of restarts.
            nres, &        ! Number of consecutive restarts.
            nout           ! Auxilary printout specification.

        ! Intrinsic Functions
        INTRINSIC ABS,MAX,SQRT,DOT_PRODUCT

       
        ! Parameters
        INTEGER, PARAMETER :: &
            maxeps = 20, & ! Maximum number of consecutive equal stopping criterions.
            maxnrs = 2000  ! Maximum number of restarts.

     
        ! Initialization

        inewnma = 1
        nser    = 0
        cnull   = 0
        ncnull  = 0
        nout    = 0
        nit     = 0
        nfe     = 0
        nge     = 0
        ntesf   = 0
        nres    = 1
        ncres   = -1
        neps    = 0
        iterm   = 0
        iters   = 1
        nnk     = 0
        alfn    = zero
        alfv    = zero
        stu     = one
        sts     = one
        utu     = one
        fold    = -large
      
        xnorm   = large

        epsr   = 0.25_prec+small
        IF (epsl+epsl >= epsr) THEN
            epsr = epsl+epsl + small
            IF (epsr >= half) THEN
                CALL wprint(iterm,-2)
            END IF
        END IF
            
     
        ! Computation of the value and the subgradient of the objective
        ! function and the search direction for the first iteration

        CALL myf(n,x,f,iterm)
        CALL myg1(n,x,g1,iterm)
        CALL myg2(n,x,g2,iterm)
        g=g1-g2

        nfe = nfe + 1
        nge = nge + 1


        IF (iterm /= 0) GO TO 900
        CALL restar(n,mccc,inewc,mccnc,inewnc,iters,gp,gp1,gp2,g,g1,g2,nnk,alfv,alfn,d,diag,cdiag,ncdiag,ncres)
     
        ! Start of the iteration
            
        iteration: DO

            ! Computation of norms

            IF (iters > 0) THEN
                gnorm = DOT_PRODUCT(g,g)
                dnorm = SQRT(DOT_PRODUCT(d,d))

                p = DOT_PRODUCT(g,d)

            ELSE
                gnorm = DOT_PRODUCT(ga,ga)
                dnorm = SQRT(DOT_PRODUCT(d,d))

                p = DOT_PRODUCT(ga,d)
            END IF

   
            ! Test on descent direction

            IF (p+small*SQRT(gnorm)*dnorm <= zero) THEN
                nres = 0
          
            ELSE

                nres = nres + 1
                IF (nres == 1) THEN
                    CALL wprint(iterm,-3)
             
                    CALL restar(n,mccc,inewc,mccnc,inewnc,iters,gp,gp1,gp2,g,g1,g2,nnk,alfv,alfn,d,diag,cdiag,ncdiag,ncres)

                    IF (ncres > maxnrs) THEN
                        nout = maxnrs
                        iterm = -2
                        EXIT iteration
                    END IF

                    CYCLE iteration
                END IF
                nout = -1
                iterm = -1
                EXIT iteration
            END IF
       
       
            ! Stopping criterion

            nit = nit + 1
            pxnorm = xnorm


            IF (iters > 0) THEN
                CALL vxdiag(n,cdiag,g,diag)  !
                xnorm = DOT_PRODUCT(g,diag) + 2.0_prec*alfv

            ELSE
                CALL vxdiag(n,cdiag,ga,diag)  !
                xnorm = DOT_PRODUCT(ga,diag) + 2.0_prec*alfv

            END IF

     
            ! Tests for termination

            IF(xnorm <= tolg) THEN   ! desired accuracy
                iterm = 1
                EXIT iteration
            END IF


            IF (nfe >= mfe) THEN  ! too many function calls
                nout = mfe
                iterm = 4
                EXIT iteration
            END IF

      
            IF (nit >= mit) THEN  ! too many iterations
                nout = mit
                iterm = 5
                EXIT iteration
            END IF

      
            IF (f <= tolb) THEN  ! too small f
                iterm = 7
                EXIT iteration
            END IF
    
      
            IF (iters == 0) THEN
                IF (ABS(xnorm - pxnorm) <= small) THEN
                    neps = neps + 1
          
                    IF (neps > maxeps) THEN
                        iterm = 8
                        EXIT iteration
                    END IF

                ELSE
                    neps = 0
                END IF

            ELSE
                neps = 0
            END IF


            IF (pxnorm < xnorm .AND. nnk > 2) THEN
                CALL wprint(iterm,-4)
            END IF
      

            CALL rprint(nit,nser,cnull,ncnull,nfe,nge,x,f,xnorm,iterm)
      
     
            !     Preparation of line search

            fo = f
       
            IF (iters > 0) THEN
                CALL copy2(n,x,xo,g,gp)
                CALL copy2(n,g1,gp1,g2,gp2)
            END IF


            t = one

            ! Line search

            IF (inma == 0) THEN
                CALL armijo(x,g,g1,g2,d,xo,fo,f,f1,f2,p,alfn,xnorm,epsr,iters,nfe,nge,iterm)
            ELSE
                CALL nmarmijo(x,g,g1,g2,d,xo,fo,f,f1,f2,fold,p,alfn,xnorm,epsr,iters,nfe,nge,inewnma,iterm)
            END IF

            IF (iterm /= 0) EXIT iteration
       
            IF (tolf2 >= 0) THEN
                IF (ABS(fo-f) <= tolf2*small*MAX(ABS(f),ABS(fo),one) &
                    .AND. iters == 1) THEN
             
                    iterm = 3
                    EXIT iteration
                END IF
            END IF

            IF (ABS(fo-f) <= tolf) THEN
                ntesf = ntesf + 1
          
                if (ntesf >= mtesf .AND. iters == 1) THEN
                    iterm = 2
                    EXIT iteration
                END IF
          
            ELSE
                ntesf = 0
            END IF
      


            ! Computation of variables and gradients differences

            CALL xdiffy(n,x,xo,s)
            CALL xdiffy(n,g1,gp1,u1)
            CALL xdiffy(n,g2,gp2,u2)
            CALL xdiffy(n,g,gp,u) ! needed in aggregation
 

            ! Computation of aggregate values

            IF (iters == 0) THEN
                nnk = nnk + 1

                IF (nnk == 1) THEN
                    CALL aggre1(g,gp,ga,u,cdiag,alfn,alfv,dnorm)
             
                ELSE
                    CALL aggre2(g,gp,ga,cdiag,alfn,alfv,dnorm)

                END IF
          
                CALL copy(n,xo,x)
                f = fo
          
            ELSE
                nnk = 0
            END IF

     

      
            ! Direction finding

            IF (iters > 0) THEN  ! Serious step update and "convex" direction determination
                nser = nser + 1
       
                ! Serious step initialization

                alfv = zero

                ! Herskovits diagonal update

                CALL update(mccc,inewc,s,u1,u2,ssmatrix,csumatrix,ncsumatrix)  ! update matrices with sts and stu

                DO i = 1, n
                    stu = zero
                    sts = zero
                    DO j = 1, mccc
                        stu = stu + csumatrix((j-1)*n+i)
                        sts = sts + ssmatrix((j-1)*n+i)
                    END DO
                    IF (stu > small) THEN
                        cdiag(i) = sts/stu
                        IF (cdiag(i) < mintold) cdiag(i) = mintold
                        IF (cdiag(i) > told) cdiag(i) = told
                    ELSE

                        cdiag(i) = told

                    END IF
                END DO

                ! Compute "convex" direction (we have a serious step, no need for "concave" direction).
                CALL vxdiag(n,cdiag,g,d)  !
                d = -one*d



            ELSE  ! Update and direction determination for null steps.
                  ! Herskovits's diagonal update

                CALL update(mccc,inewc,s,u1,u2,ssmatrix,csumatrix,ncsumatrix)  ! update matrices with sts and stu

                IF (alfn >= zero) THEN ! Convex null step
                    IF (nnk == 1) THEN  ! the first null step
                        DO i = 1, n
                            stu = zero
                            sts = zero
                            DO j = 1, mccc
                                stu = stu + csumatrix((j-1)*n+i)
                                sts = sts + ssmatrix((j-1)*n+i)
                            END DO
                            IF (stu > small) THEN
                                cdiag(i) = sts/stu
                                IF (cdiag(i) < mintold) cdiag(i) = mintold
                                IF (cdiag(i) > told) cdiag(i) = told
                            ELSE

                                cdiag(i) = told

                            END IF
                        END DO
                    END IF

                    cnull = cnull+1

                    CALL vxdiag(n,cdiag,ga,d)  !
                    d = -one*d

                ELSE ! Concave null step

                    DO i = 1, n
                        stu = zero
                        sts = zero
                        DO j = 1, mccc
                            stu = stu + ncsumatrix((j-1)*n+i)
                            sts = sts + ssmatrix((j-1)*n+i)
                        END DO
                        IF (stu > small) THEN
                            ncdiag(i) = sts/stu
                            IF (ncdiag(i) < mintold) ncdiag(i) = mintold
                            IF (ncdiag(i) > told) ncdiag(i) = told
                        ELSE
                            ncdiag(i) = told
                        END IF
                    END DO
             
                    ncnull = ncnull+1

!                    ptmp = one
                    ptmp = zero
                    DO I=1,n

!                        IF (ABS(ncdiag(i)) > small) THEN
!                            ptmp = MIN(cdiag(i)/ncdiag(i),ptmp)
!                        END IF
!                        IF (ABS(ncdiag(i)+cdiag(i)) > small) THEN
!                            ptmp = MAX(ncdiag(i)/(ncdiag(i)+cdiag(i)),ptmp)
!                        END IF
                        IF (ABS(ncdiag(i)+cdiag(i)) > small) THEN
                            ptmp = MAX((mintold2+ncdiag(i))/(ncdiag(i)+cdiag(i)),ptmp)
                        END IF

                    END DO

!                    IF (ptmp < zero) ptmp = zero
                    IF (ptmp > one) ptmp = one
                
                    DO I=1,n
!                        diag(i) = cdiag(i)-ptmp*ncdiag(i)
                        diag(i) = ptmp*cdiag(i)-(one-ptmp)*ncdiag(i)
                    END DO

                    CALL vxdiag(n,diag,ga,d)
                    d = -one*d

                END IF
            END IF

        END DO iteration
      

900 CONTINUE
      

    ! Printout the final results

    CALL wprint(iterm,nout)
    CALL rprint(nit,nser,cnull,ncnull,nfe,nge,x,f,xnorm,iterm)

CONTAINS


    ! Initialization
    SUBROUTINE restar(n,mccc,inewc,mccnc,inewnc,iters,gp,gp1,gp2,g,g1,g2,nnk,alfv,alfn,d,diag,cdiag,ncdiag,ncres)
      
        ! USE param, ONLY : zero,one  ! given in host
        ! USE subpro, ONLY : copy     ! given in host
        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
            gp, &     ! Basic subgradient of the objective function.
            gp1, &    ! Basic gradient of f1.
            gp2       ! Basic subgradient of f2.
        REAL(KIND=prec), DIMENSION(:), INTENT(INOUT) :: &
            g, &      ! Current (auxiliary) subgradient of the objective function.
            g1        ! Current (auxiliary) gradient of the f1.
        REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
            d, &      ! Search direction.
            g2, &     ! Current (auxiliary) subgradient of the f2.
            diag, &   ! diag = I.
            cdiag, &  ! cdiag = I.
            ncdiag    ! ncdiag = -I.

        ! Scalar Arguments
        INTEGER, INTENT(IN) :: &
            n         ! Number of variables (given in host)
        INTEGER, INTENT(OUT) :: &
            mccc, &   !
            mccnc, &  !
            inewc, &  !
            inewnc, & !
            nnk       ! Consecutive null steps counter.
        INTEGER, INTENT(INOUT) :: &
            iters, &  ! Null step indicator.
                      !   0  - Null step.
                      !   1  - Serious step.
            ncres     ! Number of restarts.
        REAL(KIND=prec), INTENT(OUT) :: &
            alfn, &   ! Linearization error.
            alfv      ! Aggregate linearization error.


        ! Restart
        mccc   = 0
        mccnc  = 0
        inewc  = 1
        inewnc = 1
        ncres = ncres + 1
        
        IF (iters == 0) THEN
            g = gp
            g1 = gp1
            g2 = gp2
            iters = 1
            nnk = 0
            alfv=zero
            alfn=zero
        END IF

        d = -one*g
        diag = one
        cdiag = one
        ncdiag = -one
      
    END SUBROUTINE restar


    ! Armijo line search
    SUBROUTINE armijo(x,g,g1,g2,d,xo,fo,f,f1,f2,p,alfn,wk,epsr,iters,nfe,nge,iterm)
        ! USE param, ONLY : one,zero  ! given in host
        USE initdcdb, ONLY : &
            ! n, &         ! Number of variables.         ! given in host
            ! epsl, &      ! Descent parameter.           ! given in host
            ! eta, &       ! Distance measure parameter.  ! given in host
            maxnin, &      ! Maximum number of interpolations.
            sigma          ! Tolerance for alfn.
        ! USE dc_fun, ONLY : &
            ! myf          ! Computation of the value f = f1 - f2.
            ! myg1         ! Computation of the gradient of the first DC-component,
            !              ! given in host
            ! myg2         ! Computation of the subgradient of the second DC-component,
                           ! given in host.
        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
            d, &           ! Direction vector.
            xo             ! Previous vector of variables.
        REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
            x, &           ! Vector of variables.
            g, &           ! Subgradient of the objective function.
            g1,g2          ! Subgradients of the component functions.

        ! Scalar Arguments
        REAL(KIND=prec), INTENT(INOUT) :: &
            epsr           ! Linesearch parameter, used only in null step test to print an error message
        REAL(KIND=prec), INTENT(IN) :: &
            fo, &          ! Previous value of the objective function.
            wk             ! Stopping parameter.
        REAL(KIND=prec), INTENT(OUT) :: &
            f, &           ! Value of the objective function.
            f1,f2, &       ! Values of the component functions.
            p, &           ! Directional derivative.
            alfn           ! Linearization error.
        INTEGER, INTENT(INOUT) :: &
            nfe, &         ! Number of function evaluations.
            nge            ! Number of subgradient evaluations.
        INTEGER, INTENT(OUT) :: &
            iters, &       ! Null step indicator.
                           !   0  - Null step.
                           !   1  - Serious step.
            iterm          ! Cause of termination:
                           !   0  - Everything is ok.
                           !  -3  - Failure in function or subgradient calculations.

        ! Local arrays
        REAL(KIND=prec), DIMENSION(n) :: &
            xu             ! trial point

        ! Local Scalars
        REAL(KIND=prec) :: &
            fu,fu1,fu2, &  ! Auxiliary values of the objective functions
            a_orig, &      ! Armijo parameter.
            inv_a_orig, &  ! Inverse of Armijo parameter.
            alpha, &       ! Armijo step size.
            epslwk,epsrwk  ! Auxiliary scalars.
        INTEGER :: i,nin   ! Number of interpolations.

        ! Intrinsic Functions
        INTRINSIC ABS,MAX,DOT_PRODUCT

        ! Initialization
        nin = 0

        epslwk   = epsl*wk
        epsrwk   = epsr*wk

        alpha  = 2.0_prec
        a_orig = alpha
        inv_a_orig = one/a_orig

        ! Function evaluation at a new point

        x = xo + d

        CALL myf(n,x,f,iterm)
        nfe = nfe + 1
        IF (iterm /= 0) THEN
            f=fo
            RETURN
        END IF



        ! Serious step

        iters = 1

        IF (f <= fo - epslwk) THEN ! descent condition

            xu = xo + alpha*d
            CALL myf(n,xu,fu,iterm)
            nfe = nfe + 1
            IF (iterm /= 0) RETURN

            IF(fu > fo - alpha*epslwk .OR. maxnin==0) THEN ! Serious step with t=1

                CALL myg1(n,x,g1,iterm)
                CALL myg2(n,x,g2,iterm)
                g=g1-g2
                nge = nge + 1
                IF (iterm /= 0) RETURN

                p = DOT_PRODUCT(g,d)
                alfn = fo-f+p ! Linearization error
                IF (alfn > -sigma .AND. alfn < zero) alfn=zero

                RETURN
            END IF

            inct: DO ! Increase t

                nin=nin+1
                f=fu
                x=xu
                alpha = alpha * a_orig
                xu = xo + alpha*d

                CALL myf(n,xu,fu,iterm)
                nfe = nfe + 1
                IF (iterm /= 0) RETURN

                IF (fu > fo - alpha*epslwk .OR. nin >= maxnin) THEN ! Serious step

                    CALL myg1(n,x,g1,iterm)
                    CALL myg2(n,x,g2,iterm)
                    g=g1-g2
                    nge = nge + 1
                    IF (iterm /= 0) RETURN

                    p = alpha * inv_a_orig * DOT_PRODUCT(g,d)
                    alfn = fo-f+p ! Linearization error
                    IF (alfn > -sigma .AND. alfn < zero) alfn=zero

                    RETURN
                END IF
            END DO inct

        ELSE

            alpha = one

            dect: DO ! Decrease t
                nin=nin+1
                IF (nin > maxnin) EXIT dect

                alpha = alpha * inv_a_orig !
                x = xo + alpha * d
                CALL myf(n,x,f,iterm)
                nfe = nfe + 1
                IF (iterm /= 0) RETURN

                IF (f <= fo - alpha*epslwk) THEN ! Serious step

                    CALL myg1(n,x,g1,iterm)
                    CALL myg2(n,x,g2,iterm)
                    g=g1-g2
                    nge = nge + 1
                    IF (iterm /= 0) RETURN

                    p = alpha * DOT_PRODUCT(g,d)
                    alfn = fo-f+p ! Linearization error
                    IF (alfn > -sigma .AND. alfn < zero) alfn=zero

                    RETURN
                END IF
            END DO dect

        END IF


        iters = 0 ! Null step

        ! original d
        !       x = xo + d
        !      CALL myf(n,x,f,iterm)
        !      nfe = nfe + 1


        CALL myg1(n,x,g1,iterm) ! x = xo + alpha*d
        CALL myg2(n,x,g2,iterm) ! x = xo + alpha*d
        g=g1-g2
        nge = nge + 1
        IF (iterm /= 0) RETURN

        !       p = DOT_PRODUCT(g,d) ! original d
        p = alpha * DOT_PRODUCT(g,d)
        alfn = fo-f+p ! Linearization error

        IF (alfn > -sigma .AND. alfn < zero) alfn=zero

    ! IF (p-alfn < -epsrwk) PRINT*,'not a good null step',p-alfn, -epsrwk
    ! IF (p-max(abs(alfn),eta*alpha*alpha*DOT_PRODUCT(d,d)) < -epsrwk) PRINT*,'not a good null step',p-max(abs(alfn),eta*alpha*alpha*DOT_PRODUCT(d,d)), -epsrwk

    END SUBROUTINE armijo


    ! Nonmonotone Armijo line search.
    SUBROUTINE nmarmijo(x,g,g1,g2,d,xo,fo,f,f1,f2,fold,p,alfn,wk,epsr,iters,nfe,nge,inewnma,iterm)
        ! USE param, ONLY : one,zero
        USE initdcdb, ONLY : &
            ! n, &         ! Number of variables.         ! given in host
            ! epsl, &      ! Descent parameter.           ! given in host
            ! eta, &       ! Distance measure parameter.  ! given in host
            ! mnma, &      ! Maximum number of function values used in line search.  ! given in host
            maxnin, &      ! Maximum number of interpolations.
            sigma          ! Tolerance for alfn.
        ! USE dc_fun, ONLY : &
            ! myf          ! Computation of the value f = f1 - f2.
            ! myg1         ! Computation of the subgradient of the first DC-component,
            !              ! given in host
            ! myg2         ! Computation of the subgradient of the second DC-component,
            !              ! given in host
        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
            d, &           ! Direction vector.
            xo             ! Previous vector of variables.
        REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
            x, &           ! Vector of variables.
            g, &           ! Subgradient of the objective function.
            g1,g2          ! Subgradients of the component functions.
        REAL(KIND=prec), DIMENSION(:), INTENT(INOUT) :: &
            fold           ! Old function values.

        ! Scalar Arguments
        REAL(KIND=prec), INTENT(INOUT) :: &
            epsr           ! Linesearch parameter, used only in null step test to print error
        REAL(KIND=prec), INTENT(IN) :: &
            fo, &          ! Previous value of the objective function.
            wk             ! Stopping parameter.
        REAL(KIND=prec), INTENT(OUT) :: &
            f, &           ! Value of the objective function.
            f1,f2, &       ! Values of the component functions.
            p, &           ! Directional derivative.
            alfn           ! Linearization error.
        INTEGER, INTENT(INOUT) :: &
            nfe, &         ! Number of function evaluations.
            nge, &         ! Number of subgradient evaluations.
            iters, &       ! Null step indicator.
                           !   0  - Null step.
                           !   1  - Serious step.
            inewnma        ! Index for circular arrays.
        INTEGER, INTENT(OUT) :: &
            iterm          ! Cause of termination:
                           !   0  - Everything is ok.
                           !  -3  - Failure in function or subgradient calculations.

        ! Local arrays
        REAL(KIND=prec), DIMENSION(n) :: &
            xu             ! trial point

        ! Local Scalars
        REAL(KIND=prec) :: &
            fu,fu1,fu2, &  ! Auxiliary values of the objective function
            a_orig, &      ! Armijo parameter.
            inv_a_orig, &  ! Inverse of Armijo parameter.
            alpha, &       ! Armijo step size.
            epslwk,epsrwk  ! Auxiliary scalars.
        INTEGER :: i,nin   ! Number of interpolations.

        ! Intrinsic Functions
        INTRINSIC ABS,MAX,MAXVAL,DOT_PRODUCT

        ! Initialization

        nin = 0

        epslwk   = epsl*wk
        epsrwk   = epsr*wk

        alpha  = 2.0_prec
        a_orig = alpha
        inv_a_orig = one/a_orig

        ! Updating circular array fold

        IF (iters==1) THEN
            fold(inewnma) = fo
            inewnma = inewnma + 1
            IF (inewnma > mnma) inewnma = 1
        END IF


        ! Function evaluation at a new point

        x = xo + d

        CALL myf(n,x,f,iterm)
        nfe = nfe + 1
        IF (iterm /= 0) THEN
            f=fo
            RETURN
        END IF



        ! Serious step

        iters = 1

        IF (f <= MAXVAL(fold) - epslwk) THEN ! Nonmonotone descent condition

            xu = xo + alpha*d
            CALL myf(n,xu,fu,iterm)
            nfe = nfe + 1
            IF (iterm /= 0) RETURN

            IF(fu > fo - alpha*epslwk .OR. maxnin==0) THEN ! Serious step with t=1

                CALL myg1(n,x,g1,iterm)
                CALL myg2(n,x,g2,iterm)
                g=g1-g2
                nge = nge + 1
                IF (iterm /= 0) RETURN

                p = DOT_PRODUCT(g,d)
                alfn = fo-f+p ! Linearization error
                IF (alfn > -sigma .AND. alfn < zero) alfn=zero

                RETURN
            END IF

            inct: DO ! Increase t

                nin=nin+1
                f=fu
                x=xu
                alpha = alpha * a_orig
                xu = xo + alpha*d

                CALL myf(n,xu,fu,iterm)
                nfe = nfe + 1
                IF (iterm /= 0) RETURN

                IF (fu > fo - alpha*epslwk .OR. nin >= maxnin) THEN ! Serious step

                    CALL myg1(n,x,g1,iterm)
                    CALL myg2(n,x,g2,iterm)
                    g=g1-g2
                    nge = nge + 1
                    IF (iterm /= 0) RETURN

                    p = alpha * inv_a_orig * DOT_PRODUCT(g,d)
                    alfn = fo-f+p ! Linearization error
                    IF (alfn > -sigma .AND. alfn < zero) alfn=zero

                    RETURN
                END IF
            END DO inct

        ELSE

            alpha = one

            dect: DO ! Decrease t
                nin=nin+1
                IF (nin > maxnin) EXIT dect

                alpha = alpha * inv_a_orig !
                x = xo + alpha * d
                CALL myf(n,x,f,iterm)
                nfe = nfe + 1
                IF (iterm /= 0) RETURN

                IF (f <= MAXVAL(fold) - alpha*epslwk) THEN ! Serious step

                    CALL myg1(n,x,g1,iterm)
                    CALL myg2(n,x,g2,iterm)
                    g=g1-g2
                    nge = nge + 1
                    IF (iterm /= 0) RETURN

                    p = alpha * DOT_PRODUCT(g,d)
                    alfn = fo-f+p ! Linearization error
                    IF (alfn > -sigma .AND. alfn < zero) alfn=zero

                    RETURN
                END IF
            END DO dect

        END IF


        iters = 0 ! Null step

        ! original d
        !       x = xo + d
        !      CALL myf(n,x,f,iterm)
        !      nfe = nfe + 1


        CALL myg1(n,x,g1,iterm) ! x = xo + alpha*d
        CALL myg2(n,x,g2,iterm) ! x = xo + alpha*d
        g=g1-g2
        nge = nge + 1
        IF (iterm /= 0) RETURN

        !       p = DOT_PRODUCT(g,d) ! original d
        p = alpha * DOT_PRODUCT(g,d)
        alfn = fo-f+p ! Linearization error

        IF (alfn > -sigma .AND. alfn < zero) alfn=zero

    !       IF (p-alfn < -epsrwk) PRINT*,'not a good null step',p-alfn, -epsrwk
    !       IF (p-max(abs(alfn),eta*alpha*alpha*DOT_PRODUCT(d,d)) < -epsrwk) PRINT*,'not a good null step',p-max(abs(alfn),eta*alpha*alpha*DOT_PRODUCT(d,d)), -epsrwk

    END SUBROUTINE nmarmijo


    ! Matrix update for DCDB
    SUBROUTINE update(mcc,inew,s,u1,u2,ssmatrix,csumatrix,ncsumatrix)
      
        ! USE subpro, ONLY : &
            ! copy2        ! Copying of two vectors.  ! Given in host
        ! USE initdcdb, ONLY : &
            ! n, &         ! Number of variables.     ! Given in host
            ! mcu          ! Upper limit for maximum number of stored corrections.  ! Given in host
        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
            s, &           ! Difference of current and previous variables.
            u1, &          ! Difference of current and previous gradients of f1.
            u2             ! Difference of current and previous gradients of f2.
        REAL(KIND=prec), DIMENSION(:), INTENT(INOUT) :: &
            ssmatrix, &    ! Matrix whose columns are stored corrections.
            csumatrix, &   ! Matrix whose columns are stored corrections for f1.
            ncsumatrix     ! Matrix whose columns are stored corrections for f2.

        ! Scalar Arguments
        INTEGER, INTENT(INOUT) :: &
            mcc, &         ! Current number of stored corrections.
            inew           ! Index for circular arrays.
      
        ! Local Arrays
        REAL(KIND=prec), DIMENSION(n) :: &
            ss,su1,su2     ! Componentwise multiplication of s*s and s*u.

        ! Local Scalars
        INTEGER :: i

     
        ! Compute new values

        DO i = 1, n
            ss(i)=s(i)*s(i)
            su1(i)=s(i)*u1(i)
            su2(i)=s(i)*u2(i)
        END DO
         
     
        ! Update ssmatrix and sumatrix

        CALL copy(n,ss,ssmatrix((inew-1)*n+1:))
        CALL copy2(n,su1,csumatrix((inew-1)*n+1:),su2,ncsumatrix((inew-1)*n+1:))
           
        inew = inew + 1
        IF (inew > mcu) inew = 1
        IF (mcc < mcu) mcc = mcc + 1
          
    END SUBROUTINE update


    ! Computation of aggregate values after first null step.
    SUBROUTINE aggre1(g,gp,ga,u,diag,alfn,alfv,dnorm)

        ! USE param, ONLY : zero,half,one     ! Given in host
        ! USE subpro, ONLY : &                ! Given in host
            ! vxdiag       ! Vector multiplied by a diagonal matrix.
        ! USE initdcdb, ONLY : &              ! Given in host
            ! n, &         ! Number of variables.
            ! eta          ! Distance measure parameter, eta >= 0.

        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
            diag, &        ! Diagonal variable metric update.
            g, &           ! Current (auxiliary) subgradient of the objective function.
            gp, &          ! Previous subgradient of the objective function.
            u              ! Difference of trial and aggregate gradients.
        REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
            ga             ! Next aggregate subgradient of the objective function.

        ! Scalar Arguments
        REAL(KIND=prec), INTENT(OUT) :: &
            alfv           ! Aggregate locality measure.
        REAL(KIND=prec), INTENT(IN) :: &
            dnorm          ! sqrt(norm(d)).
        REAL(KIND=prec), INTENT(INOUT) :: &
            alfn           ! Locality measure.

        ! Local arrays
        REAL(KIND=prec), DIMENSION(n) :: tmpn

        ! Local Scalars
        REAL(KIND=prec) :: &
            tmp, &         ! temporary storing of alfn.
            p, &           ! trans(gp)*diag*u - alfn.
            q, &           ! q = trans(u)*diag*u, where diag is the diagonal inverse
                           ! approximation of the Hessian.
            lam            ! Multiplier used to calculate aggregate values.
        INTEGER :: i

        ! Intrinsic Functions
        INTRINSIC ABS,MAX,MIN,SIGN,DOT_PRODUCT

        tmp=alfn
        alfn=MAX(ABS(alfn), eta*dnorm*dnorm)

      
        ! Computation of p=trans(gp)*diag*u - alfn and the product q=trans(u)*diag*u

        CALL vxdiag(n,diag,u,tmpn)
        p = - DOT_PRODUCT(tmpn,gp) - alfn
        q = DOT_PRODUCT(u,tmpn)

    
        lam = half + SIGN(half,p)

        IF (q > zero) lam = MIN(one,MAX(zero,p/q))
      

        ! Computation of the aggregate values

        p = one - lam
        DO i=1,n
            ga(i)=lam*g(i) + p*gp(i)
        END DO
      
        alfv = lam*alfn

        alfn=tmp
      
    END SUBROUTINE aggre1


    ! Computation of aggregate values after consequtive null steps.
    SUBROUTINE aggre2(g,gp,ga,diag,alfn,alfv,dnorm)
      
        ! USE param, ONLY : zero,one,small     ! Given in host
        ! USE subpro, ONLY : &                 ! Given in host
            ! vxdiag, &    ! Vector multiplied by a diagonal matrix.
            ! xdiffy       ! Difference of two vectors.
        ! USE initdcdb, ONLY : &               ! Given in host
            ! n, &         ! Number of variables.
            ! eta          ! Distance measure parameter, eta >= 0.
        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
            diag, &        ! Diagonal variable metric update.
            g, &           ! Current (auxiliary) subgradient of the objective function.
            gp             ! Previous subgradient of the objective function.
        REAL(KIND=prec), DIMENSION(:), INTENT(INOUT) :: &
            ga             ! Aggregate subgradient of the objective function.

        ! Scalar Arguments
        REAL(KIND=prec), INTENT(INOUT) :: &
            alfn, &        ! Locality measure.
            alfv           ! Aggregate locality measure.
        REAL(KIND=prec), INTENT(IN) :: &
            dnorm          ! sqrt(norm(d)).
      
        ! Local arrays
        REAL(KIND=prec), DIMENSION(n) :: tmpn2,tmpn3,tmpn4

        ! Local Scalars
        REAL(KIND=prec) :: &
            tmp, &         ! temporary storage for alfn.
            pr, &          ! pr = trans(gp-ga) diag (gp-ga).
            rrp, &         ! rrp = trans(gp-ga) diag ga - alfv.
            prqr, &        ! prqr = trans(gp-ga) diag (g-ga).
            rrq, &         ! rrq = trans(g-ga) diag ga - alfv + alfn.
            qr, &          ! qr = trans(g-ga) diag (g-ga).
            pq, &          ! pq = trans(g-gp) diag (g-gp).
            qqp, &         ! qqp = trans(g-gp) diag g + alfn.
            lam1, &        ! Multiplier used to calculate aggregate values.
            lam2, &        ! Multiplier used to calculate aggregate values.
            tmp1, &        ! Auxiliary scalar.
            tmp2, &        ! Auxiliary scalar.
            tmp3           ! Auxiliary scalar.
        INTEGER :: i

      
        ! Intrinsic Functions
        INTRINSIC ABS,MIN,MAX,DOT_PRODUCT

        tmp=alfn
        alfn=MAX(ABS(alfn), eta*dnorm*dnorm)

        CALL xdiffy(n,gp,ga,tmpn2)
     
      
        ! Calculation of tmpn3 = trans(gp - ga)dm

        CALL vxdiag(n,diag,tmpn2,tmpn3)

         
        pr = DOT_PRODUCT(tmpn3,tmpn2)
        rrp = DOT_PRODUCT(tmpn3,ga)
        CALL xdiffy(n,g,ga,tmpn4)
        prqr = DOT_PRODUCT(tmpn3,tmpn4)


        ! calculation of qr = trans(g - ga) diag (g - ga) and rrq = trans(g-ga) diag ga + alfn - alfv.
       
        CALL vxdiag(n,diag,tmpn4,tmpn3)
        rrq = DOT_PRODUCT(tmpn3,ga)
        qr = DOT_PRODUCT(tmpn4,tmpn3)
      
        pq = qr - prqr - prqr + pr
        qqp = pq + prqr + rrq - pr - rrp + alfn
        rrp = rrp - alfv
        rrq = rrq + alfn - alfv


        ! computation of multipliers lam1 and lam2

        IF (pr > zero .AND. qr > zero) THEN
            tmp1 = rrq/qr
            tmp2 = prqr/qr
            tmp3 = pr - prqr*tmp2
            IF (tmp3 /= zero) THEN
                lam1 = (tmp1*prqr - rrp)/tmp3
                lam2 = -tmp1 - lam1*tmp2
                IF (lam1*(lam1 - one) < zero .AND. &
                    lam2*(lam1 + lam2 - one) < zero) THEN
                    IF (lam1 == zero .AND. lam2*(lam2 - one) < zero &
                        .AND. -rrp - lam2*prqr > zero .AND. pr > zero) &
                        lam1 = MIN(one - lam2, (-rrp-lam2*prqr)/pr)

                    ! Computation of the aggregate values

                    tmp1 = one - lam1 - lam2
                    DO i=1,n
                        ga(i)=lam1*gp(i)+lam2*g(i)+tmp1*ga(i)
                    END DO

                    alfv = lam2*alfn + tmp1*alfv

                    alfn=tmp

                    RETURN

                END IF

            END IF
        END IF


        ! Minimum on the boundary

        lam1 = zero
        lam2 = zero
        IF (alfn <= alfv) lam2 = one
        IF (qr > zero) lam2 = MIN(one,MAX(zero,-rrq/qr))
        tmp3 = (lam2*qr + rrq+rrq)*lam2
        tmp1 = zero
        IF (alfv >= zero) tmp1 = one
        IF (pr > zero) tmp1 = MIN(one,MAX(zero,-rrp/pr))
        tmp2 = (tmp1*pr + rrp+rrp)*tmp1
        IF (tmp2 < tmp3) THEN
            tmp3 = tmp2
            lam1 = tmp1
            lam2 = zero
        END IF

        IF (qqp*(qqp - pq) < zero) THEN
            IF (qr + rrq + rrq - qqp*qqp/pq < tmp3) THEN
                lam1 = qqp/pq
                lam2 = one - lam1
            END IF
        END IF


        IF (lam1 == zero .AND. lam2*(lam2 - one) < zero &
            .AND. -rrp - lam2*prqr > zero .AND. pr > zero) &
            lam1 = MIN(one - lam2, (-rrp-lam2*prqr)/pr)


        ! Computation of the aggregate values
      
        tmp1 = one - lam1 - lam2
        DO i=1,n
            ga(i)=lam1*gp(i)+lam2*g(i)+tmp1*ga(i)
        END DO
    
        alfv = lam2*alfn + tmp1*alfv

        alfn=tmp

    
    END SUBROUTINE aggre2
      
END SUBROUTINE dcdb

!************************************************************************
!*
!*     * SUBROUTINE wprint *
!*
!*     Printout the warning and error messages.
!*
!************************************************************************
      
SUBROUTINE wprint(iterm,nout)
    USE initdcdb, ONLY : &
        iprint             ! Printout specification:
                           !  -1  - No printout.
                           !   0  - Only the error messages.
                           !   1  - The final values of the objective function.
                           !   2  - The final values of the objective function and the
                           !        most serious warning messages.
                           !   3  - The whole final solution.
                           !   4  - At each iteration values of the objective function.
                           !   5  - At each iteration the whole solution
    IMPLICIT NONE

    ! Scalar Arguments
    INTEGER, INTENT(IN) :: &
        nout, &            ! Auxilary printout specification.
        iterm              ! Cause of termination:
                           !   1  - The problem has been solved with desired accuracy.
                           !   2  - Changes in function values < tolf in mtesf
                           !        subsequent iterations.
                           !   3  - Changes in function value < tolf*small*MAX(|f_k|,|f_(k-1)|,1),
                           !        where small is the smallest positive number such that
                           !        1.0 + small > 1.0.
                           !   4  - Number of function calls > mfe.
                           !   5  - Number of iterations > mit.
                           !   6  - Time limit exceeded.
                           !   7  - f < tolb.
                           !   8  - Failure in attaining the demanded accuracy.
                           !  -1  - Two consecutive restarts.
                           !  -2  - Number of restarts > maximum number of restarts.
                           !  -3  - Failure in function or subgradient calculations
                           !        (assigned by the user).
                           !  -4
                           !  -5  - Invalid input parameters.


    IF (iprint >= 0) THEN

        ! Initial error messages

        IF (iterm == -5) THEN
            IF (nout == 1) WRITE (6,FMT='(1X,''Error: '' &
               ''Number of variables (n) is too small, iterm='',I3)'            ) iterm
            IF (nout == 2) WRITE (6,FMT='(1X,''Error: '' &
               ''The maximum number of stored corrections (mcu) '' &
               ''is too small, iterm='',I3)'            ) iterm
            IF (nout == 4) WRITE (6,FMT='(1X,''Error: '' &
               ''Line search parameter epsl >= 0.25, iterm='',I3)'            ) iterm
            RETURN
        END IF

        
        ! Warning messages

        IF (iprint >= 2) THEN
            IF (iterm == 0) THEN
                IF (nout == -2) WRITE (6,FMT='(1X,''Warning: '' &
                  ''A line search parameter epsr >= 0.5.'')'                )
                IF (nout == -3) WRITE (6,FMT='(1X,''Warning: '' &
                  ''A nondescent search direction occured. Restart.'')'                )
                IF (nout == -4) WRITE (6,FMT='(1X,''Warning: '' &
                  ''Does not converge.'')'                )
                RETURN
            END IF
         

            ! Printout the final results
            
            IF (iterm == 2) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''Too many steps without significant progress.'')'            )
            IF (iterm == 3) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''The value of the function does not change.'')'            )
            IF (iterm == 4) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''Number of function evaluations > '',I5)'            ) nout
            IF (iterm == 5) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''Number of iterations > '',I5)'            ) nout
            IF (iterm == 6) WRITE (6,FMT='(1X,''Abnormal exit: Time is up.'')')
            IF (iterm == 7) WRITE (6,FMT='(1X,''Abnormal exit: f < tolb.'')')
            IF (iterm == 8) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''Failure in attaining the demanded accuracy.'')'            )
            IF (iterm == -1) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''Two consecutive restarts.'')'            )
            IF (iterm == -2) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''Number of restarts > '',I5''.'')'            ) nout
            IF (iterm == -3) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''Failure in function or subgradient calculations.'')'            )
        END IF
    END IF
      
END SUBROUTINE wprint

            
!************************************************************************
!*
!*     * SUBROUTINE rprint *
!*      
!*     Printout the (final) results.
!*
!************************************************************************
      
SUBROUTINE rprint(nit,nser,cnull,ncnull,nfe,nge,x,f,wk,iterm)
    USE initdcdb, ONLY : &
        n, &               ! Number of variables
        iprint             ! Printout specification:
                           !  -1  - No printout.
                           !   0  - Only the error messages.
                           !   1  - The final values of the objective function.
                           !   2  - The final values of the objective function and the
                           !        most serious warning messages.
                           !   3  - The whole final solution.
                           !   4  - At each iteration values of the objective function.
                           !   5  - At each iteration the whole solution
    IMPLICIT NONE

    ! Scalar Arguments
    INTEGER, INTENT(IN) :: & 
        nit, &             ! Number of used iterations.
        nser, &            ! Number of serious steps.
        cnull, &           ! Number of convex null steps.
        ncnull, &          ! Number of nonconvex null steps.
        nfe, &             ! Number of used function evaluations.
        nge, &             ! Number of used subgradient evaluations.
        iterm              ! Cause of termination:
                           !   1  - The problem has been solved with desired accuracy.
                           !   2  - Changes in function values < tolf in mtesf
                           !        subsequent iterations.
                           !   3  - Changes in function value < tolf*small*MAX(|f_k|,|f_(k-1)|,1),
                           !        where small is the smallest positive number such that
                           !        1.0 + small > 1.0.
                           !   4  - Number of function calls > mfe.
                           !   5  - Number of iterations > mit.
                           !   6  - Time limit exceeded.
                           !   7  - f < tolb.
                           !   8  - Failure in attaining the demanded accuracy.
                           !  -1  - Two consecutive restarts.
                           !  -2  - Number of restarts > maximum number of restarts.
                           !  -3  - Failure in function or subgradient calculations
                           !        (assigned by the user).
                           !  -4
                           !  -5  - Invalid input parameters.


    REAL(KIND=prec), INTENT(IN) :: &
        f, &               ! Value of the objective function.
        wk                 ! Value of the stopping criterion.

    ! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
        x                  ! Vector of variables
         
    ! Local Scalars
    INTEGER :: i

    ! Intermediate results
    
    IF (iterm == 0) THEN
        IF (iprint > 3) WRITE (6,FMT='(1X,''nit='',I5,2X, &
            ''nfe='',I5,2X,''nge='',I5,2X,''f='',D15.8,2X,''wk='',D11.4,2X &
            )'        ) nit,nfe,nge,f,wk
        IF (iprint == 5) WRITE (6,FMT='(1X,''x='', &
            5D15.7:/(4X,5D15.7))'        )(x(i),i=1,n)
        RETURN
    END IF
         

    ! Final results

    IF (iprint > 0) WRITE (6,FMT='(1X,''nit='',I6,2X,''nser='',I6,2X,''cnull='',I6,2X, &
         ''ncnull='',I6,2X,''nfe='',I7,2X,''nge='',I7,2X,''f='',D15.8,2X,''wk='',D11.4,2X, &
         ''iterm='',I3)'    ) nit,nser,cnull,ncnull,nfe,nge,f,wk,iterm
    IF (IPRINT .EQ. 3 .OR. IPRINT .EQ. 5) &
        WRITE (6,FMT='(1X,''x='',5D15.7:/(4X,5D15.7))')(x(i),i=1,n)
      
END SUBROUTINE rprint
      
END MODULE dcdb_mod
