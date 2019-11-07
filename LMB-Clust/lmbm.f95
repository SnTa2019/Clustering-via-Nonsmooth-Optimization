!*************************************************************************
!*                                                                       *
!*     LMBM     - Limited Memory Bundle Method for nonsmooth             *
!*                nonconvex optimization                                 *
!*                                                                       *
!*                                                                       *
!*     by Napsu Karmitsa 2004 (last modified 31.05.2016)                 *
!*                                                                       *
!*                                                                       *
!*     This file is specially modified to solve clustering problems.     *
!*                                                                       *
!*************************************************************************
!*
!*     Modules included:
!*
!*     lmbm_mod            ! Limited Memory Bundle Method
!*
MODULE lmbm_mod  ! Limited memory bundle method

    USE r_precision, ONLY : prec  ! Precision for reals.
    IMPLICIT NONE

    ! MODULE lmbm_mod includes the following subroutines (S) and functions (F).
    PUBLIC :: &
        optim2        ! S Initializing LMBM subroutine for clustering.
    PRIVATE :: &
        lmbm, &       ! S Limited memory bundle subroutine for nonsmooth
                      !   large-scale optimization. Contains:
                      !     S restar  Initialization and reinitialization.
                      !     S dobun   Bundle construction.
        indic1, &     ! S Initialization of indices.
        tinit, &      ! S Calculation of initial step size. Contains:
                      !     S destep  Stepsize selection using polyhedral
                      !                 approximation for descent steps.
                      !     S nulstep Stepsize selection using polyhedral
                      !                 approximation for null steps.
        lls, &        ! S Weak Wolfe line search. Contains:
                      !     F qint    Quadratic interpolation.
        nmlls, &      ! S Nonmonotonic Weak Wolfe line search. Contains:
                      !     F qint    Quadratic interpolation.
        armijo, &     ! S Armijo line search.
        nmarmijo, &   ! S Nonmonotone Armijo line search.
        dlbfgs, &     ! S Computing the search direction by limited memory
                      !   BFGS update. Contains:
                      !     F sclpar  Selection of scaling parameter.
        dlsr1, &      ! S Computing the search direction by limited
                      !   memory SR1 update.
        agbfgs, &     ! S Simplified subgradient aggregation.
        aggsr1, &     ! S Subgradient aggregation.
        agskip, &     ! S Subgradient aggregation using BFGS update.
        wprint, &     ! S Printout the error and warning messages.
        rprint        ! S Printout the results.

CONTAINS

    SUBROUTINE optim2(x0,x1,fvalue)

        USE initclust, ONLY : &
            maxdim, &      ! Maximum number of variables.
            mf, &          ! Number of used features:
                           !   mf = nft, when no classes,
            nc, &          ! Current number of clusters.
            ns, &          ! Switch for auxiliary and real clustering problem.
            m              ! Number of variables for function calculations.
        USE initlmbm, ONLY : &
            n, &           ! Number of variables.
            x, &           ! Vector of variables.
            mcu, &         ! Maximum number of stored corrections.
            mcinit, &      ! Initial maximum number of stored corrections.
            epsl, &        ! Line search parameter.
            defaults, &    ! S  Default values for parameters.
            init_lmbmpar   ! S  Further initialization of parameters.

        REAL(KIND=prec) ::  &
            x0(nc*mf),     &
            x1(nc*mf),     &
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
        INTEGER :: i,mc

        ! number of variables
        IF(ns.EQ.1) n = mf
        IF(ns.EQ.2) n = nc*mf
        m=n

        DO i = 1, n
            x(i) = x0(i)
        END DO

        mc = mcinit     ! Initial maximum number of stored corrections

        CALL defaults()    !
        CALL init_lmbmpar()    !

        IF (n <= 0) PRINT*,'n<0'
        IF (n > maxdim) PRINT*,'n is greater than maxdim.'
        IF (epsl >= 0.25_prec) PRINT*,'epsl >= 0.25'
        IF (mcu <= 0) PRINT*,'mcu <= 0'

        CALL lmbm(mc,fvalue,iout(1),iout(2),iout(3),iout(4))

        DO i=1, n
            x1(i) = x(i)
        END DO

        RETURN

    END SUBROUTINE optim2


    !***********************************************************************
    !*                                                                     *
    !*     * SUBROUTINE lmbm *                                             *
    !*                                                                     *
    !*     Limited memory bundle subroutine for nonsmooth optimization.    *
    !*                                                                     *
    !***********************************************************************
      
    SUBROUTINE lmbm(mc,f,nit,nfe,nge,iterm)
        USE param, ONLY : small,large,zero,half,one
        USE initlmbm, ONLY : &
            n, &         ! Number of variables.
            x, &         ! Vector of variables
            na, &        ! Maximum bundle dimension, na >= 2.
            mcu, &       ! Upper limit for maximum number of stored corrections, mcu >= 3.
            inma, &      ! Selection of line search method:
                         !   inma = 0, Armijo line search,
                         !   inma = 1, nonmonotone Armijo line search.
                         !   inma = 2, weak Wolfe line search.
            mnma, &      ! Maximum number of function values used in nonmonotone line search.
            iprint, &    ! Printout specification (see initializat for more details).
            iscale, &    ! Selection of the scaling (see initializat for more details).
            tolf, &      ! Tolerance for change of function values.
            tolf2, &     ! Second tolerance for change of function values.
            tolb, &      ! Tolerance for the function value.
            tolg, &      ! Tolerance for the first termination criterion.
            tolg2, &     ! Tolerance for the second termination criterion.
            xmax, &      ! Maximum stepsize, 1 < XMAX.
            eta, &       ! Distance measure parameter, eta >= 0.
            epsl, &      ! Line search parameter, 0 < epsl < 0.25.
            mtesf, &     ! Maximum number of iterations with changes of
                         ! function values smaller than tolf.
            mit, &       ! Maximun number of iterations.
            mfe          ! Maximun number of function evaluations.
        USE obj_fun, ONLY : &
            myf, &       ! Computation of the value of the objective function.
            myg          ! Computation of the subgradient of the objective function.
        USE subpro, ONLY : &
            xdiffy, &    ! Difference of two vectors.
            copy, &      ! Copying of a vector.
            copy2        ! Copying of two vectors.


        IMPLICIT NONE

        ! Scalar Arguments
        REAL(KIND=prec), INTENT(OUT) :: &
            f            ! Value of the objective function.
        INTEGER, INTENT(OUT) :: &
            nit,nfe,nge  ! Number of iterations, and function and subgradient evaluations.
        INTEGER, INTENT(INOUT) :: &
            mc           ! Maximum number of stored corrections.
        INTEGER, INTENT(OUT) :: &
            iterm        ! Cause of termination:
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
                          !  -1  - Two consecutive restarts.
                          !  -2  - Number of restarts > maximum number of restarts.
                          !  -3  - Failure in function or subgradient calculations
                          !        (assigned by the user).
                          !  -4  - Failure in attaining the demanded accuracy.
                          !  -5  - Invalid input parameters.
                          !  -6  - Unspecified error.


        ! Local Arrays
        REAL(KIND=prec), DIMENSION(n) :: &
            xo, &        ! Previous vector of variables.
            g, &         ! Subgradient of the objective function.
            gp, &        ! Previous subgradient of the ohjective function.
            ga, &        ! Aggregate subgradient.
            s, &         ! Difference of current and previous variables.
            u, &         ! Difference of current and previous subgradients.
            d, &         ! Direction vector.
            tmpn1        ! Auxiliary array.
        REAL(KIND=prec), DIMENSION(n*(mcu+1)) :: &
            sm,um        ! Matrises whose columns are stored differences of
                          ! variables (sm) and subgradients (um).
        REAL(KIND=prec), DIMENSION((mcu+2)*(mcu+1)/2) :: &
            rm, &        ! pper triangular matrix stored columnwise in the
                         ! one-dimensional array.
            umtum        ! Matrix whose columns are stored subgradient differences.
        REAL(KIND=prec), DIMENSION(mcu+1) :: &
            cm, &        ! Diagonal matrix.
            smtgp, &     ! smtgp = trans(sm)*gp.
            umtgp, &     ! umtgp = trans(um)*gp.
            tmpmc1, &    ! Auxiliary array.
            tmpmc2       ! Auxiliary array.
        REAL(KIND=prec), DIMENSION(n*na) :: &
            ax,ag        ! Matrix whose columns are bundle points and subgradients.
        REAL(KIND=prec), DIMENSION(na) :: &
            af           ! Vector of bundle values.
        REAL(KIND=prec), DIMENSION(mnma) :: &
            fold         ! Old function values.

        ! Local Scalars
        REAL(KIND=prec) :: &
            alfn, &      ! Locality measure.
            alfv, &      ! Aggregate locality measure.
            epsr, &      ! Line search parameter.
            dnorm, &     ! Euclidean norm of the direction vector.
            gnorm, &     ! Euclidean norm of the aggregate subgradient.
            xnorm, &     ! Stopping criterion.
            pxnorm, &    ! Previous stopping criterion.
            p, &         ! Directional derivative.
            tmax, &      ! Maximum stepsize.
            t, &         ! Stepsize.
            theta, &     ! Correction parameter for stepsize.
            fo, &        ! Previous value of the objective.
            gamma        ! Scaling parameter.
        INTEGER :: i, &
            mcinit, &    ! Initial maximum number of stored corrections.
            mcc, &       ! Current number of stored corrections.
            inew, &      ! Index for the circular arrays.
            ibfgs, &     ! Index of the type of BFGS update.
            isr1, &      ! Index of the type of SR1 update.
            iters, &     ! Null step indicator.
                         !   0  - Null step.
                         !   1  - Serious step.
            nnk, &       ! Consecutive null steps counter.
            ibun, &      ! Index for the circular arrays in bundle updating.
            mal, &       ! Current size of the bundle.
            inewnma, &   ! Index for the circular arrays containing function values.
            ic, &        ! Correction indicator.
            icn, &       ! Correction indicator for null steps.
            iflag, &     ! Index for adaptive version.
            neps, &      ! Number of consecutive equal stopping criterions.
            ntesf, &     ! Number of tests on function decrease.
            ncres, &     ! Number of restarts.
            nres, &      ! Number of consecutive restarts.
            nress, &     ! Number of consecutive restarts in case of tmax < tmin.
            nout         ! Auxilary printout specification.

        ! Intrinsic Functions
        INTRINSIC ABS,MAX,SQRT,DOT_PRODUCT

       
        ! Parameters
        REAL(KIND=prec), PARAMETER :: &
            fmin    = -large, &        ! Smallest acceptable value of the function.
            tmin    = 1.0E-12_prec, &  ! Minimum stepsize.
            lengthd = 1.0E+20_prec, &  ! Direction vector length.
            rho     = 1.0E-12_prec     ! Correction parameter.
        INTEGER, PARAMETER :: &
            maxeps = 20, &             ! Maximum number of consecutive equal stopping criterions.
            maxnrs = 2000              ! Maximum number of restarts.

        ! Initialization

        inewnma = 1
        nout   = 0
        nit    = 0
        nfe    = 0
        nge    = 0
        ntesf  = 0
        nres   = 1
        ncres  = -1
        nress  = 0
        neps   = 0
        iterm  = 0
        iters  = 1
        nnk    = 0
        isr1   = 0
        alfn   = zero
        alfv   = zero
        mcinit = mc
      
        tmax   = xmax
        xnorm  = large
        fold   = -large

        epsr   = 0.25_prec+small
        IF (epsl+epsl >= epsr) THEN
            epsr = epsl+epsl + small
            IF (epsr >= half) THEN
                CALL wprint(iterm,iprint,-2)
            END IF
        END IF
            
     
        ! Computation of the value and the subgradient of the objective
        ! function and the search direction for the first iteration

        CALL myf(n,x,f,iterm)
        CALL myg(n,x,g,iterm)
        nfe = nfe + 1
        nge = nge + 1

        IF (iterm /= 0) GO TO 900
    
        CALL restar(n,mc,mcc,mcinit,inew,ibun,ibfgs,iters,gp,g,nnk, &
            alfv,alfn,gamma,d,ic,icn,mal,ncres,iflag)
      
        CALL dobun(n,na,mal,x,g,f,ax,ag,af,iters,ibun)

     
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
                    CALL wprint(iterm,iprint,-3)
             
                    CALL restar(n,mc,mcc,mcinit,inew,ibun,ibfgs,iters,gp,g,nnk, &
                        alfv,alfn,gamma,d,ic,icn,mal,ncres,iflag)
                    IF (ncres > maxnrs) THEN
                        nout = maxnrs
                        iterm = -2
                        EXIT iteration
                    END IF
             
                    CALL dobun(n,na,mal,x,g,f,ax,ag,af,iters,ibun)
             
                    CYCLE iteration
                END IF
                nout = -1
                iterm = -1
                EXIT iteration
            END IF
       
       
            ! Stopping criterion

            nit = nit + 1
            pxnorm = xnorm
            xnorm = -p + 2.0_prec*alfv

     
            ! Tests for termination

            IF (xnorm <= 1.0E+03_prec*tolg .AND. &
                (mcc > 0 .OR. ibfgs == 2)) THEN
          
                IF(half*gnorm + alfv <= tolg2 .AND. &
                    xnorm <= tolg) THEN
          
                    iterm = 1
                    EXIT iteration
                END IF
       
                IF (mc < mcu .AND. iflag == 0) THEN
                    mc = mc+1
                    iflag = 1
                END IF
            END IF


            IF (nfe >= mfe) THEN
                nout = mfe
                iterm = 4
                EXIT iteration
            END IF

      
            IF (nit >= mit) THEN
                nout = mit
                iterm = 5
                EXIT iteration
            END IF

      
            IF (f <= tolb) THEN
                iterm = 7
                EXIT iteration
            END IF
    
      
            IF (iters == 0) THEN
                IF (ABS(xnorm - pxnorm) <= small) THEN
                    neps = neps + 1
          
                    IF (neps > maxeps) THEN
                        iterm = -4
                        EXIT iteration
                    END IF

                ELSE
                    neps = 0
                END IF

            ELSE
                neps = 0
            END IF


            ! Correction
    
            IF (-p < rho*gnorm .OR. icn == 1) THEN

                xnorm = xnorm + rho*gnorm
                dnorm = SQRT(dnorm*dnorm-2.0_prec*rho*p+rho*rho*gnorm)
         
                IF (iters > 0) THEN
                    DO i=1,n
                        d(i) = d(i)-rho*g(i)
                    END DO

                ELSE
                    DO i=1,n
                        d(i) = d(i)-rho*ga(i)
                    END DO
                    icn = 1
                END IF

                ic = 1
            
            ELSE
                ic = 0
            END IF


            IF (pxnorm < xnorm .AND. nnk > 2) THEN
                CALL wprint(iterm,iprint,-4)
            END IF
      

            CALL rprint(n,nit,nfe,nge,x,f,xnorm,half*gnorm+alfv,iterm,iprint)
      
     
            !     Preparation of line search

            fo = f
       
            IF (iters > 0) THEN
                CALL copy2(n,x,xo,g,gp)
            END IF

            IF (dnorm > zero) tmax = xmax/dnorm
       
            IF (tmax > tmin) THEN
                nress = 0

            ELSE
                nress = nress + 1
                IF (nress == 1) THEN
                    CALL wprint(iterm,iprint,-5)

                    CALL restar(n,mc,mcc,mcinit,inew,ibun,ibfgs,iters,gp,g,nnk, &
                        alfv,alfn,gamma,d,ic,icn,mal,ncres,iflag)
             
                    IF (ncres > maxnrs) THEN
                        nout = maxnrs
                        iterm = -2
                        EXIT iteration
                    END IF
             
                    CALL dobun(n,na,mal,x,g,f,ax,ag,af,iters,ibun)
                    CYCLE iteration
                END IF
                iterm = -1
                EXIT iteration
            END IF
      

            ! Initial step size

            CALL tinit(n,na,mal,x,af,ag,ax,ibun,d,f,p,t,tmax,tmin, &
                eta,iters,iterm)

            IF (iterm /= 0) EXIT iteration
     
            IF (inma == 0) THEN
                CALL armijo(x,g,d,xo,fo,f,t,p,alfn,xnorm,epsr,iters,nfe,nge,iterm)

            ELSE IF (inma == 1) THEN
                CALL nmarmijo(x,g,d,xo,fo,f,fold,t,p,alfn,xnorm,epsr,iters,nfe,nge,inewnma,iterm)

            ELSE IF (inma == 2) THEN! Line search with directional derivatives which allows null steps
                 ! With this the global convergence can be guaranteed.
                theta = one
                IF (dnorm > lengthd) THEN
                    theta=lengthd/dnorm
                END IF
      
                CALL lls(n,x,g,d,xo,t,fo,f,p,alfn,tmin,dnorm,xnorm, &
                    theta,epsl,epsr,eta,iters,nfe,nge,nnk,iterm)
            ELSE ! Nonmonotone line search with directional derivatives which allows null steps
                 ! With this the global convergence can be guaranteed.
                theta = one
                IF (dnorm > lengthd) THEN
                    theta=lengthd/dnorm
                END IF

                CALL nmlls(n,x,g,d,xo,t,fo,f,fold,p,alfn,tmin,dnorm,xnorm, &
                    theta,epsl,epsr,eta,iters,inewnma,nfe,nge,nnk,iterm)

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
      

            ! Bundle updating

            CALL dobun(n,na,mal,x,g,f,ax,ag,af,iters,ibun)
  

            ! Computation of variables difference

            CALL xdiffy(n,x,xo,s)
  

            ! Computation of aggregate values and gradients difference

            IF (iters == 0) THEN
                nnk = nnk + 1

                IF (nnk == 1) THEN
                    CALL copy(n,gp,tmpn1)
                    CALL xdiffy(n,g,gp,u)
             
                    CALL agbfgs(n,mc,mcc,inew,ibfgs,iflag,g,gp,ga,u,d,sm,um, &
                        rm,cm,umtum,alfn,alfv,gamma,ic,rho)
             
                ELSE
                    CALL copy(n,ga,tmpn1)
                    CALL aggsr1(n,mc,mcc,inew,iflag,g,gp,ga,d,alfn,alfv &
                        ,umtum,rm,gamma,smtgp,umtgp,tmpmc1,tmpmc2,sm &
                        ,um,icn,rho)
                    CALL xdiffy(n,g,gp,u)
                END IF
          
                CALL copy(n,xo,x)
                f = fo
          
            ELSE
                IF (nnk /= 0) THEN
                    CALL copy(n,ga,tmpn1)
                ELSE
                    CALL copy(n,gp,tmpn1)
                END IF
                nnk = 0
                CALL xdiffy(n,g,gp,u)
            END IF

     
            ! Serious step initialization

            IF (iters > 0) THEN
                icn = 0
                alfn = zero
                alfv = zero
            END IF

      
            ! Direction finding
    
            IF (iters > 0) THEN
         
     
                ! BFGS update and direction determination

                CALL dlbfgs(n,mc,mcc,inew,ibfgs,iflag,d,g,gp,s,u,sm,um,rm, &
                    umtum,cm,smtgp,umtgp,gamma,tmpn1,iscale)

            ELSE

                ! SR1 update and direction determination
             
                CALL dlsr1(n,mc,mcc,inew,isr1,iflag,d,gp,ga,s,u,sm,um,rm, &
                    umtum,cm,smtgp,umtgp,gamma,tmpmc1,tmpmc2,tmpn1,nnk,iprint)
                ibfgs=0
            END IF
    
        END DO iteration
      

900 CONTINUE
      

    ! Printout the final results

    CALL wprint(iterm,iprint,nout)
    CALL rprint(n,nit,nfe,nge,x,f,xnorm,half*gnorm+alfv,iterm,iprint)

CONTAINS

    SUBROUTINE restar(n,mc,mcc,mcinit,inew,ibun,ibfgs,iters,gp,g,nnk, &
        alfv,alfn,gamma,d,ic,icn,mal,ncres,iflag)  ! Initialization
      
        !      USE param, ONLY : zero,one  ! given in host
        !      USE lmbm_sub, ONLY : copy  ! given in host
        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
            gp        ! Basic subgradient of the objective function.
        REAL(KIND=prec), DIMENSION(:), INTENT(INOUT) :: &
            g         ! Current (auxiliary) subgradient of the objective function.
        REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
            d         ! Search direction.

        ! Scalar Arguments
        INTEGER, INTENT(IN) :: &
            n, &      ! Number of variables
            mcinit    ! Initial maximum number of stored corrections.
        INTEGER, INTENT(OUT) :: &
            mc, &     ! Current maximum number of stored corrections.
            mcc, &    ! Current number of stored corrections.
            inew, &   ! Index for the circular arrays.
            ibun, &   ! Index for the circular arrays in bundle updating.
            ibfgs, &  ! Index of the type of BFGS update.
            nnk, &    ! Consecutive null steps counter.
            ic, &     ! Correction indicator.
            icn, &    ! Correction indicator for null steps.
            mal, &    ! Current size of the bundle.
            iflag     ! Index for adaptive version.
        INTEGER, INTENT(INOUT) :: &
            iters, &  ! Null step indicator.
                      !   0  - Null step.
                      !   1  - Serious step.
            ncres     ! Number of restarts.
        REAL(KIND=prec), INTENT(OUT) :: &
            alfn, &   ! Locality measure.
            alfv, &   ! Aggregate locality measure.
            gamma     ! Scaling parameter.


        ! Restart
        mc    = mcinit
        mcc   = 0
        inew  = 1
        ibun  = 1
        ibfgs = 0
        ic    = 0
        icn   = 0
        mal   = 0
        ncres = ncres + 1
        iflag = 0
        
        IF (iters == 0) THEN
            CALL copy(n,gp,g)
            iters = 1
            nnk = 0
            alfv=zero
            alfn=zero
        END IF

        gamma = one
        d=-g
      
    END SUBROUTINE restar
    
      
    SUBROUTINE dobun(n,ma,mal,x,g,f,ay,ag,af,iters,ibun)  ! Bundle construction

        !      USE lmbm_sub, ONLY : copy2 ! given in host
        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
            x, &      ! Vector of variables
            g         ! Subgradient of the objective function.
        REAL(KIND=prec), DIMENSION(:), INTENT(INOUT) :: &
            ay, &     ! Matrix whose columns are bundle points
                      ! (stored in one-dimensional n*ma -array).
            ag, &     ! Matrix whose columns are bundle subgradients
                      ! (stored in one-dimensional n*ma -array).
            af        ! Vector of values of bundle functions
                       ! (stored in one-dimensional ma -array).

        ! Scalar Arguments
        INTEGER, INTENT(IN) :: &
            n, &      ! Number of variables
            iters, &  ! Null step indicator.
                      !   0  - Null step.
                      !   1  - Serious step.
            ma        ! Maximum size of the bundle.
        INTEGER, INTENT(INOUT) :: &
            ibun, &   ! Index for the circular arrays in bundle updating.
            mal       ! Current size of the bundle.
        REAL(KIND=prec), INTENT(IN) :: &
            f         ! Value of the objective function.

        ! Local Scalars
        INTEGER :: i,j
      
        IF (iters == 1) THEN
     
            ! Serious step

            af(ibun) = f
            i = (ibun-1)*n + 1
            CALL copy2(n,g,ag(i:),x,ay(i:))
         
        ELSE

            ! Null step

            IF (mal < ma) THEN
            
                af(ibun) = af(mal)
                af(mal) = f
            
                i = mal*n + 1
                CALL copy2(n,ag(i-n:),ag(i:),ay(i-n:),ay(i:))
                CALL copy2(n,g,ag(i-n:),x,ay(i-n:))
            
            ELSE
                i = ibun-1
                IF (i < 1) i = mal
                af(ibun) = af(i)
                af(i) = f
            
                i = (ibun-2)*n + 1
                IF (i < 1) i = (mal-1)*n + 1
                j = (ibun-1)*n + 1
                CALL copy2(n,ag(i:),ag(j:),ay(i:),ay(j:))
                CALL copy2(n,g,ag(i:),x,ay(i:))
            END IF
         
        END IF
      
        mal = mal + 1
        IF (mal > ma) mal = ma
      
        ibun = ibun + 1
        IF (ibun > ma) ibun = 1
      
    END SUBROUTINE dobun
    
END SUBROUTINE lmbm

!************************************************************************
!*                                                                      *
!*     * SUBROUTINE tinit *                                             *
!*                                                                      *
!*     Initial stepsize selection for limited memory bundle method      *
!*                                                                      *
!************************************************************************
SUBROUTINE tinit(n,na,mal,x,af,ag,ay,ibun,d,f,p,t,tmax,tmin,eta,iters,iterm)

    USE param, ONLY : zero,half,one,large
    !    USE param, ONLY : zero,one  ! these are the one needed in tinit itself
                                 ! half and large are used in destep and nulstep
    IMPLICIT NONE

    ! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
        x, &      ! Vector of variables (n array).
        d, &      ! Direction vector (n array).
        ay, &     ! Matrix whose columns are bundle points (stored in one-dimensional n*ma -array).
        ag        ! Matrix whose columns are bundle subgradients (stored in one-dimensional n*ma -array).
    REAL(KIND=prec), DIMENSION(:), INTENT(INOUT) :: &
        af        ! Vector of values of bundle functions (stored in one-dimensional ma -array).

    ! Scalar Arguments
    REAL(KIND=prec), INTENT(OUT) :: & 
        t         ! Initial stepsize
    REAL(KIND=prec), INTENT(IN) :: & 
        p, &      ! Directional derivative.
        eta, &    ! Distance measure parameter.
        f, &      ! Value of the objective function.
        tmax, &   ! Upper limit for stepsize parameter.
        tmin      ! Lower limit for stepsize parameter.
    INTEGER, INTENT(IN) :: & 
        n, &      ! Number of variables
        na, &     ! Maximum size of the bundle.
        mal, &    ! Current size of the bundle.
        ibun, &   ! Index for the circular arrays in bundle updating.
        iters     ! Null step indicator.
                   !   0  - Null step.
                   !   1  - Serious step.
    INTEGER, INTENT(OUT) :: &
        iterm     ! Cause of termination:
                   !   0  - Everything is ok.
                   !  -6  - Error.


    ! Intrinsic Functions
    INTRINSIC MAX,MIN

    t = MIN(one,tmax)

    IF (p == zero) RETURN

    IF (iters == 1) THEN
        CALL destep(n,na,mal,x,af,ag,ay,ibun,d,f,p,t,eta,iterm)
    ELSE
        CALL nulstep(n,na,mal,x,af,ag,ay,ibun,d,f,p,t,eta,iterm)
    END IF

    t = MIN(MAX(t,tmin),tmax)

CONTAINS

    SUBROUTINE destep(n,ma,mal,x,af,ag,ay,ibun,d,f,df,t,eta,iterm)  ! Stepsize selection

        !      USE param, ONLY : zero,half,one,large  ! given in host
        IMPLICIT NONE

        ! Scalar Arguments
        REAL(KIND=prec), INTENT(INOUT) :: &
            t         ! Initial stepsize
        REAL(KIND=prec), INTENT(IN) :: &
            df, &     ! Directional derivative.
            eta, &    ! Distance measure parameter.
            f         ! Value of the objective function.
        INTEGER, INTENT(IN) :: &
            n, &      ! Number of variables
            ma, &     ! Maximum size of the bundle.
            mal, &    ! Current size of the bundle.
            ibun      ! Index for the circular arrays in bundle updating.
        INTEGER, INTENT(OUT) :: &
            iterm     ! Cause of termination:
                       !   0  - Everything is ok.
                       !  -6  - Error.

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
            x, &      ! Vector of variables (n array).
            d, &      ! Direction vector (n array).
            ay, &     ! Matrix whose columns are bundle points (stored in one-dimensional n*ma -array).
            ag, &     ! Matrix whose columns are bundle subgradients (stored in one-dimensional n*ma -array).
            af        ! Vector of values of bundle functions (stored in one-dimensional ma -array).

        ! Local Arrays
        REAL(KIND=prec), DIMENSION(2*ma) :: &
            tmparray  ! Auxiliary array.

        ! Local Scalars
        REAL(KIND=prec) :: alf,alfl,alfr,bet,betl,betr,dx,q,r,w
        INTEGER :: i,j,jn,k,l,lq,ib

        ! Intrinsic Functions
        INTRINSIC ABS,REAL,MAX,MIN,SQRT

        iterm = 0
        alfl = zero
        betl = zero
      
        w = df*t* (one - t*half)
      
     
        ! Initial choice of possibly active lines
      
        k = 0
        l = -1
        jn = (ibun-1)*n
        betr = - large
        DO j=1,mal-1
            ib = ibun - 1 + j
            IF (ib > mal) ib = ib - mal
            IF (jn >= mal*n) jn = jn - mal*n
            r = zero
            bet = zero
            alfl = af(ib) - f
            DO i=1,n
                dx = x(i) - ay(jn+i)
                q = ag(jn+i)
                r = r + dx*dx
                alfl = alfl + dx*q
                bet = bet + d(i)*q
            END DO
            alf = MAX(ABS(alfl),eta*r)
            r = one - bet/df
            IF (r*r + (alf+alf)/df > 1.0E-6_prec) THEN
                k = k + 1
                tmparray(k) = alf
                tmparray(ma+k) = bet
                r = t*bet - alf
                IF (r > w) THEN
                    w = r
                    l = k
                END IF
            END IF
         
            betr = MAX(betr,bet-alf)
            jn = jn + n
        END DO
        lq = -1
        IF (betr <= df*half) RETURN
        lq = 1
        IF (l <= 0) RETURN
        betr = tmparray(ma+l)
        IF (betr <= zero) THEN
            IF (t < one .OR. betr == zero) RETURN
            lq = 2
        END IF
      
        alfr = tmparray(l)


        ! Iteration loop

        ds_iteration: DO
         
            IF (lq >= 1) THEN
                q = one - betr/df
                r = q + SQRT(q*q + (alfr+alfr)/df)
                IF (betr >= zero) r = - (alfr+alfr)/ (df*r)
                r = MIN(1.95_prec,MAX(zero,r))
            ELSE
                IF (ABS(betr-betl)+ABS(alfr-alfl)< -1.0E-4_prec*df) RETURN
                IF (betr-betl  == zero) THEN
                    iterm = -6
                    RETURN
                END IF
                r = (alfr-alfl)/ (betr-betl)
            END IF
       
            IF (ABS(t-r) < 1.0E-4_prec) RETURN
            t = r
            tmparray(l) = - one
            w = t*betr - alfr
            l = -1
            DO j = 1,k
                alf = tmparray(j)
                IF (alf < zero) EXIT
                bet = tmparray(ma+j)
                r = t*bet - alf
                IF (r > w) THEN
                    w = r
                    l = j
                END IF
            END DO
            IF (l < 0) RETURN
       
            bet = tmparray(ma+l)
            IF (bet == zero) RETURN
     

            !     New interval selection

            alf = tmparray(l)
            IF (bet < zero) THEN
                IF (lq == 2) THEN
                    alfr = alf
                    betr = bet
               
                ELSE
                    alfl = alf
                    betl = bet
                    lq = 0
                END IF

            ELSE
                IF (lq == 2) THEN
                    alfl = alfr
                    betl = betr
                    lq = 0
                END IF
       
                alfr = alf
                betr = bet
            END IF
       
        END DO ds_iteration
      
    END SUBROUTINE destep

     
    SUBROUTINE nulstep(n,ma,mal,x,af,ag,ay,ibun,d,f,df,t,eta,iterm)  ! Stepsize selection

        !      USE param, ONLY : zero,half,one,large  ! given in host
        IMPLICIT NONE

        ! Scalar Arguments
        REAL(KIND=prec), INTENT(INOUT) :: &
            t         ! Initial stepsize
        REAL(KIND=prec), INTENT(IN) :: &
            df, &     ! Directional derivative.
            eta, &    ! Distance measure parameter.
            f         ! Value of the objective function.
        INTEGER, INTENT(IN) :: &
            n, &      ! Number of variables
            ma, &     ! Maximum size of the bundle.
            mal, &    ! Current size of the bundle.
            ibun      ! Index for the circular arrays in bundle updating.
        INTEGER, INTENT(OUT) :: &
            iterm     ! Cause of termination:
                       !   0  - Everything is ok.
                       !  -6  - Error.

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
            x, &      ! Vector of variables (n array).
            d, &      ! Direction vector (n array).
            ay, &     ! Matrix whose columns are bundle points (stored in one-dimensional n*ma -array).
            ag, &     ! Matrix whose columns are bundle subgradients (stored in one-dimensional n*ma -array).
            af        ! Vector of values of bundle functions (stored in one-dimensional 4*ma -array).

        ! Local Arrays
        REAL(KIND=prec), DIMENSION(2*ma) :: &
            tmparray  ! Auxiliary array.

        ! Local Scalars
        REAL(KIND=prec) :: alf,alfl,alfr,bet,betl,betr,dx,q,r,w
        INTEGER :: i,j,jn,k,l,ib

        ! Intrinsic Functions
        INTRINSIC ABS,REAL,MAX,MIN,SQRT

        iterm = 0
        w = df*t
     
        ! Initial choice of possibly active parabolas

        k = 0
        l = -1
        jn = (ibun-1)*n
        betr = - large
        DO j = 1,mal - 1
            ib = ibun - 1 + j
            IF (ib > mal) ib = ib - mal
            IF (jn >= mal*n) jn = jn - mal*n
            bet = zero
            r = zero
            alfl = af(ib) - f
            DO i = 1,n
                dx = x(i) - ay(jn+i)
                r = r + dx*dx
                q = ag(jn+i)
                alfl = alfl + dx*q
                bet = bet + d(i)*q
            END DO
            alf = MAX(ABS(alfl),eta*r)
            betr = MAX(betr,bet-alf)
            IF (alf < bet-df) THEN
                k = k + 1
                r = t*bet - alf
                tmparray(k) = alf
                tmparray(ma+k) = bet
                IF (r > w) THEN
                    w = r
                    l = k
                END IF
            END IF
            jn = jn + n
        END DO
        IF (l <= 0) RETURN
        betr = tmparray(ma+l)
        alfr = tmparray(l)
        alf = alfr
        bet = betr
        alfl = zero
        betl = df
    
     
        !     Iteration loop
    
        ns_iteration: DO

            w = bet/df
            IF (ABS(betr-betl)+ABS(alfr-alfl)< -1.0E-4_prec*df) RETURN
            IF (betr-betl  == zero) THEN
                iterm = -6
                RETURN
            END IF
            r = (alfr-alfl)/ (betr-betl)
            IF (ABS(t-w) < ABS(t-r)) r = w
            q = t
            t = r
            IF (ABS(t-q) < 1.0E-3_prec) RETURN
            tmparray(l) = - one
            w = t*bet - alf
            l = -1
            DO j=1,k
                alf = tmparray(j)
                IF (alf < zero) EXIT
                bet = tmparray(ma+j)
                r = t*bet - alf
                IF (r > w) THEN
                    w = r
                    l = j
                END IF
            END DO

            IF (l <= 0) RETURN
            bet = tmparray(ma+l)
            q = bet - t*df
            IF (Q == zero) RETURN

     
            !     New interval selection

            alf = tmparray(l)
            IF (q < zero) THEN
                alfl = alf
                betl = bet
            ELSE
                alfr = alf
                betr = bet
            END IF

        END DO ns_iteration

    END SUBROUTINE nulstep

END SUBROUTINE tinit
      
!************************************************************************
!*                                                                      *
!*     * SUBROUTINE lls *                                               *
!*                                                                      *
!*     Special line search for limited memory bundle method             *
!*                                                                      *
!************************************************************************

SUBROUTINE lls(n,x,g,d,xo,t,fo,f,p,alfn,tmin,dnorm,wk,theta,epsl,epsr,&
    eta,iters,nfe,nge,nnk,iterm)

    USE param, ONLY : zero,half,one
    USE obj_fun, ONLY : &
        myf,myg          ! Computation of the value and the subgradient of
                         ! the objective function.
    USE subpro, ONLY : &
        scsum            ! Sum of a vector and the scaled vector.
    USE initlmbm, ONLY : &
        maxint => maxnin ! Maximum number of interpolations.

    IMPLICIT NONE

    ! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
        d, &      ! Direction vector.
        xo        ! Previous vector of variables.
    REAL(KIND=prec), DIMENSION(:), INTENT(INOUT) :: &
        x         ! Vector of variables.
    REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
        g         ! Subgradient of the objective function.

    ! Scalar Arguments
    REAL(KIND=prec), INTENT(INOUT) :: & 
        epsl,&    ! Linesearch parameter.
        epsr,&    ! Linesearch parameter.
        t, &      ! Stepsize
        p         ! Directional derivative.
    REAL(KIND=prec), INTENT(IN) :: & 
        theta, &  ! Scaling parameter.
        eta, &    ! Distance measure parameter.
        fo, &     ! Previous value of the objective function.
        wk, &     ! Stopping parameter.
        dnorm, &  ! Euclidean norm of the direction vector.
        tmin      ! Lower limit for stepsize parameter.
    REAL(KIND=prec), INTENT(OUT) :: & 
        f, &      ! Value of the objective function.
        alfn      ! Locality measure.
    INTEGER, INTENT(IN) :: & 
        n, &      ! Number of variables
        nnk       ! Number of consequtive null steps.
    INTEGER, INTENT(INOUT) :: & 
        nfe, &    ! Number of function evaluations.
        nge       ! Number of subgradient evaluations.
    INTEGER, INTENT(OUT) :: &
        iters, &  ! Null step indicator.
                  !   0  - Null step.
                  !   1  - Serious step.
        iterm     ! Cause of termination:
                  !   0  - Everything is ok.
                  !  -3  - Failure in function or subgradient calculations.

    ! Local Scalars
    INTEGER :: nin ! Number of interpolations.
    REAL(KIND=prec) :: &
        tl,tu, &  ! Lower and upper limits for t used in interpolation.
        fl,fu, &  ! Values of the objective function for t=tl and t=tu.
        epsa,epst,epslk,epsrk, & ! Line search parameters.
        thdnorm,epsawk,epstwk,epslwk,epsrwk ! Auxiliary scalars.

    ! Intrinsic Functions
    INTRINSIC ABS,MAX

      

    ! Initialization

    nin = 0

    epst = epsl+epsl
    epsa = half*(epsr - epst)
    thdnorm = theta*dnorm

    tl = zero
    tu = t
    fl = fo

    IF (theta < one) THEN
        epst  = theta*epst
        epsa  = theta*epsa
        epslk = epsl
        epsl  = theta*epsl
        epsrk = epsr
        epsr  = theta*epsr
    END IF
 
    epsawk   = epsa*wk
    epslwk   = epsl*wk
    epsrwk   = epsr*wk
    epstwk   = epst*wk


    ! Function evaluation at a new point

    lls_iteration: DO
       
        CALL scsum(n,theta*t,d,xo,x)

        CALL myf(n,x,f,iterm)
        nfe = nfe + 1
        IF (iterm /= 0) RETURN

    
        ! Null/descent step test (ITERS=0/1)

        iters = 1
        IF (f <= fo - t*epstwk) THEN
            tl = t
            fl = f
        ELSE
            tu = t
            fu = f
        END IF
      

        ! Additional interpolation

        IF (f > fo .AND. tu-tl >= tmin*0.1_prec &
            .AND. nnk >= 1 .AND. nin < maxint) THEN

            nin=nin+1
            IF (tl == zero .AND. wk > zero) THEN
                t = qint(tu,fl,fu,wk,one-half/(one-epst))
            ELSE
                t = half*(tu+tl)
            END IF
            CYCLE lls_iteration
        END IF

        CALL myg(n,x,g,iterm)
        nge = nge + 1
        IF (iterm /= 0) RETURN

        p = theta*DOT_PRODUCT(g,d)
        alfn = MAX(ABS(fo-f+p*t),eta*(t*thdnorm)**2)


        ! Serious step

        IF (f <= fo - t*epslwk .AND. (t >= tmin .OR. alfn > epsawk)) EXIT lls_iteration


        ! Null step

        IF (p-alfn >= -epsrwk .OR. tu-tl < tmin*0.1_prec .OR. &
            nin >= maxint) THEN
            ITERS = 0
            EXIT lls_iteration
        END IF


        ! Interpolation

        nin=nin+1
        IF (tl == zero .AND. wk > zero) THEN
            t = qint(tu,fl,fu,wk,one-half/(one-epst))
        ELSE
            t = half*(tu+tl)
        END IF

    END DO lls_iteration

    IF (theta /= one) THEN
        epsl = epslk
        epsr = epsrk
    END IF

CONTAINS

    FUNCTION qint(tu,fl,fu,wk,kappa) RESULT(t) ! Quadratic interpolation

        !      USE param, ONLY : half,one  ! given in host
        IMPLICIT NONE
      
        ! Scalar Arguments
        REAL(KIND=prec), INTENT(IN) :: &
            fl, &  ! Value of the objective function.
            fu, &  ! Value of the objective function for t=tu.
            wk, &  ! Directional derivative.
            tu, &  ! Upper value of the stepsize parameter.
            kappa  ! Interpolation parameter.
        REAL(KIND=prec) :: &
            t      ! Stepsize.

        ! Local Scalars
        REAL(KIND=prec) :: tmp1,tmp2

        ! Intrinsic Functions
        INTRINSIC MAX


        tmp1 = (fu-fl)/ (-wk*tu)

     
        ! Quadratic interpolation with one directional derivative

        tmp2 = 2.0_prec * (one - tmp1)

        IF (tmp2 > one) THEN
   
            ! Interpolation accepted
       
            t = MAX(kappa*tu,tu/tmp2)
            RETURN
        END IF
      
     
        ! Bisection
    
        t = half*tu
      
    END FUNCTION qint

END SUBROUTINE lls


!************************************************************************
!*                                                                      *
!*     * SUBROUTINE nmlls *                                             *
!*                                                                      *
!*     Nonmonotonic weak Wolfe-type line search                         *
!*                                                                      *
!************************************************************************

SUBROUTINE nmlls(n,x,g,d,xo,t,fo,f,fold,p,alfn,tmin,dnorm,wk,theta,epsl,epsr,&
    eta,iters,inewnma,nfe,nge,nnk,iterm)

    USE param, ONLY : zero,half,one
    USE obj_fun, ONLY : &
        myf,myg    ! Computation of the value and the subgradient of
                   ! the objective function.
    USE subpro, ONLY : &
        scsum      ! Sum of a vector and the scaled vector.
    USE initlmbm, ONLY : &
        mnma, &    ! Maximum number of function values used in line search.
        maxint => maxnin   ! Maximum number of interpolations.

    IMPLICIT NONE

    ! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
        d, &      ! Direction vector.
        xo        ! Previous vector of variables.
    REAL(KIND=prec), DIMENSION(:), INTENT(INOUT) :: &
        x         ! Vector of variables.
    REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
        g         ! Subgradient of the objective function.
    REAL(KIND=prec), DIMENSION(:), INTENT(INOUT) :: &
        fold           ! Old function values.

    ! Scalar Arguments
    REAL(KIND=prec), INTENT(INOUT) :: &
        epsl,&    ! Linesearch parameter.
        epsr,&    ! Linesearch parameter.
        t, &      ! Stepsize
        p         ! Directional derivative.
    REAL(KIND=prec), INTENT(IN) :: &
        theta, &  ! Scaling parameter.
        eta, &    ! Distance measure parameter.
        fo, &     ! Previous value of the objective function.
        wk, &     ! Stopping parameter.
        dnorm, &  ! Euclidean norm of the direction vector.
        tmin      ! Lower limit for stepsize parameter.
    REAL(KIND=prec), INTENT(OUT) :: &
        f, &      ! Value of the objective function.
        alfn      ! Locality measure.
    INTEGER, INTENT(IN) :: &
        n, &      ! Number of variables
        nnk       ! Number of consequtive null steps.
    INTEGER, INTENT(INOUT) :: &
        inewnma, &! index for array.
        nfe, &    ! Number of function evaluations.
        nge       ! Number of subgradient evaluations.
    INTEGER, INTENT(OUT) :: &
        iters, &  ! Null step indicator.
                  !   0  - Null step.
                  !   1  - Serious step.
        iterm     ! Cause of termination:
                  !   0  - Everything is ok.
                  !  -3  - Failure in function or subgradient calculations.

    ! Local Scalars
    INTEGER :: nin ! Number of interpolations.
    REAL(KIND=prec) :: &
        tl,tu, &  ! Lower and upper limits for t used in interpolation.
        fl,fu, &  ! Values of the objective function for t=tl and t=tu.
        ftmp, &   ! Maximum function value from mnma last iterations.
        epsa,epst,epslk,epsrk, & ! Line search parameters.
        thdnorm,epsawk,epstwk,epslwk,epsrwk ! Auxiliary scalars.

    ! Intrinsic Functions
    INTRINSIC ABS,MAX,MAXVAL,DOT_PRODUCT



    ! Initialization

    nin = 0

    epst = epsl+epsl
    epsa = half*(epsr - epst)
    thdnorm = theta*dnorm

    tl = zero
    tu = t
    fl = fo

    IF (theta < one) THEN
        epst  = theta*epst
        epsa  = theta*epsa
        epslk = epsl
        epsl  = theta*epsl
        epsrk = epsr
        epsr  = theta*epsr
    END IF

    epsawk   = epsa*wk
    epslwk   = epsl*wk
    epsrwk   = epsr*wk
    epstwk   = epst*wk

    ! Updating circular array fold

    IF (iters==1) THEN
        fold(inewnma) = fo
        inewnma = inewnma + 1
        IF (inewnma > mnma) inewnma = 1
    END IF
    ftmp = MAXVAL(fold)


    ! Function evaluation at a new point

    lls_iteration: DO

        CALL scsum(n,theta*t,d,xo,x)

        CALL myf(n,x,f,iterm)
        nfe = nfe + 1
        IF (iterm /= 0) RETURN


        ! Null/descent step test (ITERS=0/1)

        iters = 1
        IF (f <= ftmp - t*epstwk) THEN
            tl = t
            fl = f
        ELSE
            tu = t
            fu = f
        END IF


        ! Additional interpolation

        IF (f > ftmp .AND. tu-tl >= tmin*0.1_prec &
            .AND. nnk >= 1 .AND. nin < maxint) THEN

            nin=nin+1
            IF (tl == zero .AND. wk > zero) THEN
                t = qint(tu,fl,fu,wk,one-half/(one-epst))
            ELSE
                t = half*(tu+tl)
            END IF
            CYCLE lls_iteration
        END IF

        CALL myg(n,x,g,iterm)
        nge = nge + 1
        IF (iterm /= 0) RETURN

        p = theta*DOT_PRODUCT(g,d)
        alfn = MAX(ABS(fo-f+p*t),eta*(t*thdnorm)**2)


        ! Serious step

        IF (f <= ftmp - t*epslwk .AND. (t >= tmin .OR. alfn > epsawk)) EXIT lls_iteration


        ! Null step

        IF (p-alfn >= -epsrwk .OR. tu-tl < tmin*0.1_prec .OR. &
            nin >= maxint) THEN
            ITERS = 0
            EXIT lls_iteration
        END IF


        ! Interpolation

        nin=nin+1
        IF (tl == zero .AND. wk > zero) THEN
            t = qint(tu,fl,fu,wk,one-half/(one-epst))
        ELSE
            t = half*(tu+tl)
        END IF

    END DO lls_iteration

    IF (theta /= one) THEN
        epsl = epslk
        epsr = epsrk
    END IF

CONTAINS

    FUNCTION qint(tu,fl,fu,wk,kappa) RESULT(t) ! Quadratic interpolation

        !      USE param, ONLY : half,one  ! given in host
        IMPLICIT NONE

        ! Scalar Arguments
        REAL(KIND=prec), INTENT(IN) :: &
            fl, &  ! Value of the objective function.
            fu, &  ! Value of the objective function for t=tu.
            wk, &  ! Directional derivative.
            tu, &  ! Upper value of the stepsize parameter.
            kappa  ! Interpolation parameter.
        REAL(KIND=prec) :: &
            t      ! Stepsize.

        ! Local Scalars
        REAL(KIND=prec) :: tmp1,tmp2

        ! Intrinsic Functions
        INTRINSIC MAX


        tmp1 = (fu-fl)/ (-wk*tu)


        ! Quadratic interpolation with one directional derivative

        tmp2 = 2.0_prec * (one - tmp1)

        IF (tmp2 > one) THEN

            ! Interpolation accepted

            t = MAX(kappa*tu,tu/tmp2)
            RETURN
        END IF


        ! Bisection

        t = half*tu

    END FUNCTION qint

END SUBROUTINE nmlls



!************************************************************************
!*                                                                      *
!*     * SUBROUTINE armijo *                                            *
!*                                                                      *
!*     Armijo line search for limited memory bundle method              *
!*                                                                      *
!************************************************************************


    ! Armijo line search
    SUBROUTINE armijo(x,g,d,xo,fo,f,tin,p,alfn,wk,epsr,iters,nfe,nge,iterm)
        USE param, ONLY : two,one,zero  ! given in host
        USE initlmbm, ONLY : &
            n, &         ! Number of variables.         ! given in host
            epsl, &      ! Descent parameter.           ! given in host
            eta, &       ! Distance measure parameter.  ! given in host
            maxnin         ! Maximum number of interpolations.
        USE obj_fun, ONLY : &
            myf, &       ! Computation of the value f = f1 - f2.
            myg          ! Computation of the gradient of the first DC-component,
                         ! given in host
        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
            d, &           ! Direction vector.
            xo             ! Previous vector of variables.
        REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
            x, &           ! Vector of variables.
            g              ! Subgradient of the objective function.

        ! Scalar Arguments
        REAL(KIND=prec), INTENT(INOUT) :: &
            epsr           ! Linesearch parameter, used only in null step test to print an error message
        REAL(KIND=prec), INTENT(IN) :: &
            fo, &          ! Previous value of the objective function.
            tin, &         ! Initial stepsize.
            wk             ! Stopping parameter.
        REAL(KIND=prec), INTENT(OUT) :: &
            f, &           ! Value of the objective function.
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
            dnorm, &       ! |d|
            epslwk,epsrwk  ! Auxiliary scalars.
        INTEGER :: i,nin   ! Number of interpolations.

        ! Intrinsic Functions
        INTRINSIC ABS,MAX,DOT_PRODUCT

        ! Initialization
        nin = 0

        epslwk   = epsl*wk
        epsrwk   = epsr*wk

        alpha  = two*tin
        a_orig = two
        inv_a_orig = one/a_orig

        dnorm = DOT_PRODUCT(d,d)

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

            IF(fu > fo - alpha*epslwk .OR. maxnin==0) THEN ! Serious step with t=tin

                CALL myg(n,x,g,iterm)
                nge = nge + 1
                IF (iterm /= 0) RETURN

                p = DOT_PRODUCT(g,d)
                alfn = MAX(ABS(fo-f+p*tin),eta*(tin*dnorm)**2)

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

                    CALL myg(n,x,g,iterm)
                    nge = nge + 1
                    IF (iterm /= 0) RETURN

                    p = alpha * inv_a_orig * DOT_PRODUCT(g,d)
                    alfn = MAX(ABS(fo-f+p),eta*(alpha * inv_a_orig*dnorm)**2)

                    RETURN
                END IF
            END DO inct

        ELSE

            alpha = tin

            dect: DO ! Decrease t
                nin=nin+1
                IF (nin > maxnin) EXIT dect

                alpha = alpha * inv_a_orig !
                x = xo + alpha * d
                CALL myf(n,x,f,iterm)
                nfe = nfe + 1
                IF (iterm /= 0) RETURN

                IF (f <= fo - alpha*epslwk) THEN ! Serious step

                    CALL myg(n,x,g,iterm)
                    nge = nge + 1
                    IF (iterm /= 0) RETURN

                    p = alpha * DOT_PRODUCT(g,d)
                    alfn = MAX(ABS(fo-f+p),eta*(alpha*dnorm)**2)

                    RETURN
                END IF
            END DO dect

        END IF


        iters = 0 ! Null step

        ! original d
        !       x = xo + d
        !      CALL myf(n,x,f,iterm)
        !      nfe = nfe + 1


        CALL myg(n,x,g,iterm) ! x = xo + alpha*d
        nge = nge + 1
        IF (iterm /= 0) RETURN

        !       p = DOT_PRODUCT(g,d) ! original d
        p = alpha * DOT_PRODUCT(g,d)
        alfn = MAX(ABS(fo-f+p),eta*(alpha*dnorm)**2)


    END SUBROUTINE armijo


!************************************************************************
!*                                                                      *
!*     * SUBROUTINE nmarmijo *                                          *
!*                                                                      *
!*     Nonmonotone Armijo line search for limited memory bundle method  *
!*                                                                      *
!************************************************************************

    ! Nonmonotone Armijo line search.
    SUBROUTINE nmarmijo(x,g,d,xo,fo,f,fold,tin,p,alfn,wk,epsr,iters,nfe,nge,inewnma,iterm)
        USE param, ONLY : two,one,zero
        USE initlmbm, ONLY : &
            n, &           ! Number of variables.
            epsl, &        ! Descent parameter.
            eta, &         ! Distance measure parameter.
            mnma, &        ! Maximum number of function values used in line search.
            maxnin         ! Maximum number of interpolations.
        USE obj_fun, ONLY : &
            myf, &         ! Computation of the value f = f1 - f2.
            myg            ! Computation of the subgradient,
            !              ! given in host
        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
            d, &           ! Direction vector.
            xo             ! Previous vector of variables.
        REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
            x, &           ! Vector of variables.
            g              ! Subgradient of the objective function.
        REAL(KIND=prec), DIMENSION(:), INTENT(INOUT) :: &
            fold           ! Old function values.

        ! Scalar Arguments
        REAL(KIND=prec), INTENT(INOUT) :: &
            epsr           ! Linesearch parameter, used only in null step test to print error
        REAL(KIND=prec), INTENT(IN) :: &
            fo, &          ! Previous value of the objective function.
            tin, &         ! Initial stepsize.
            wk             ! Stopping parameter.
        REAL(KIND=prec), INTENT(OUT) :: &
            f, &           ! Value of the objective function.
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
            dnorm, &       ! |d|
            epslwk,epsrwk  ! Auxiliary scalars.
        INTEGER :: i,nin   ! Number of interpolations.

        ! Intrinsic Functions
        INTRINSIC ABS,MAX,MAXVAL,DOT_PRODUCT

        ! Initialization

        nin = 0

        epslwk   = epsl*wk
        epsrwk   = epsr*wk
        dnorm    = DOT_PRODUCT(d,d)

        alpha  = two*tin
        a_orig = two
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

            IF(fu > fo - alpha*epslwk .OR. maxnin==0) THEN ! Serious step with t=tin

                CALL myg(n,x,g,iterm)
                nge = nge + 1
                IF (iterm /= 0) RETURN

                p = DOT_PRODUCT(g,d)
                alfn = MAX(ABS(fo-f+p*tin),eta*(tin*dnorm)**2)

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

                    CALL myg(n,x,g,iterm)
                    nge = nge + 1
                    IF (iterm /= 0) RETURN

                    p = alpha * inv_a_orig * DOT_PRODUCT(g,d)
                    alfn = MAX(ABS(fo-f+p),eta*(alpha * inv_a_orig*dnorm)**2)

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

                    CALL myg(n,x,g,iterm)
                    nge = nge + 1
                    IF (iterm /= 0) RETURN

                    p = alpha * DOT_PRODUCT(g,d)
                    alfn = MAX(ABS(fo-f+p),eta*(alpha*dnorm)**2)

                    RETURN
                END IF
            END DO dect

        END IF


        iters = 0 ! Null step

        ! original d
        !       x = xo + d
        !      CALL myf(n,x,f,iterm)
        !      nfe = nfe + 1


        CALL myg(n,x,g,iterm) ! x = xo + alpha*d
        nge = nge + 1
        IF (iterm /= 0) RETURN

        !       p = DOT_PRODUCT(g,d) ! original d
        p = alpha * DOT_PRODUCT(g,d)
        alfn = MAX(ABS(fo-f+p),eta*(alpha*dnorm)**2)


    END SUBROUTINE nmarmijo




!************************************************************************
!*                                                                      *
!*     * SUBROUTINE dlbfgs *                                            *
!*                                                                      *
!*     Matrix update and computation of the search direction d = -dm*g  *
!*     by the limited memory BFGS update.                               *
!*                                                                      *
!************************************************************************
    
SUBROUTINE dlbfgs(n,mc,mcc,inew,ibfgs,iflag,d,g,gp,s,u,sm,um,rm, &
    umtum,cm,smtgp,umtgp,gamma,tmpn1,iscale)
      
    USE param, ONLY : zero,small,one,half ! half is used at sclpar
    USE subpro, ONLY : &
        xdiffy, & ! Difference of two vectors.
        xsumy, &  ! Sum of two vectors.
        scdiff, & ! Difference of the scaled vector and a vector.
        scsum, &  ! Sum of a vector and the scaled vector.
        vxdiag, & ! Multiplication of a vector and a diagonal matrix.
        symax, &  ! Multiplication of a dense symmetric matrix by a vector.
        cwmaxv, & ! Multiplication of a vector by a dense rectangular matrix.
        rwaxv2, & ! Multiplication of two rowwise stored dense rectangular
                  ! matrices A and B by vectors x and y.
        trlieq, & ! Solving x from linear equation L*x=y or trans(L)*x=y.
        copy2     ! Copying of two vectors.
    IMPLICIT NONE

    ! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
        g, &      ! Current subgradient of the objective function.
        gp        ! Previous subgradient of the objective function.
    REAL(KIND=prec), DIMENSION(:), INTENT(INOUT) :: &
        sm, &     ! Matrix whose columns are stored corrections.
        um, &     ! Matrix whose columns are stored subgradient differences.
        rm, &     ! Upper triangular matrix.
        cm, &     ! Diagonal matrix.
        umtum, &  ! Matrix umtum = trans(um) * um.
        smtgp, &  ! Vector smtgp = trans(sm)*gp.
        umtgp, &  ! vector umtgp = trans(um)*gp.
        s, &      ! Difference of current and previous variables.
        u, &      ! Difference of current and previous subgradients.
        tmpn1     ! Auxiliary array. On input: previous aggregate subgradient.
    REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
        d         ! Direction vector.

    ! Scalar Arguments
    REAL(KIND=prec), INTENT(INOUT) :: & 
        gamma     ! Scaling parameter.
    INTEGER, INTENT(IN) :: & 
        n, &      ! Number of variables
        mc, &     ! Declared number of stored corrections.
        iscale    ! Selection of the scaling:
                   !   0  - Scaling at every iteration with STU/UTU.
                   !   1  - Scaling at every iteration with STS/STU.
                   !   2  - Interval scaling with STU/UTU.
                   !   3  - Interval scaling with STS/STU.
                   !   4  - Preliminary scaling with STU/UTU.
                   !   5  - Preliminary scaling with STS/STU.
                   !   6  - No scaling.      
    INTEGER, INTENT(INOUT) :: & 
        mcc, &    ! Current number of stored corrections.
        inew, &   ! Index for circular arrays.
        iflag     ! Index for adaptive version:
                   !   0  - Maximum number of stored corrections
                   !        has not been changed at this iteration.
                   !   1  - Maximum number of stored corrections
                   !        has been changed at this iteration.
    INTEGER, INTENT(OUT) :: & 
        ibfgs     ! Index of the type of BFGS update:
                   !   1  - BFGS update: the corrections are stored.
                   !   2  - BFGS update: the corrections are not stored.
                   !   3  - BFGS update is skipped.
      
    ! Local Arrays
    REAL(KIND=prec), DIMENSION(mcc+1) :: &
        tmpmc1,tmpmc2,tmpmc3,tmpmc4

    ! Local Scalars
    REAL(KIND=prec) :: &
        stu, &    ! stu = trans(s)*u.
        sts       ! sts = trans(s)*s.
    INTEGER :: i,j,k, &
        mcnew, &  ! Current size of vectors.
        iold, &   ! Index of the oldest corrections.
        iflag2, & ! Index for adaptive version.
        ierr      ! Error indicator

    ! Intrinsic Functions
    INTRINSIC SQRT,MIN,MAX


    ierr = 0
    ibfgs = 0
    iflag2 = 0
    stu = DOT_PRODUCT(s,u)
    sts = DOT_PRODUCT(s,s)


    ! Positive definiteness

    IF (stu > zero) THEN
        IF (-DOT_PRODUCT(d,u)-DOT_PRODUCT(tmpn1,s) < -small) THEN
          
     
            ! Update matrices
         
            ibfgs = 1

            ! Initialization of indices.

            CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,ibfgs)
            
     
            ! Update sm and um

            CALL copy2(n,s,sm((inew-1)*n+1:),u,um((inew-1)*n+1:))
           
     
            ! Computation of trans(sm)*g and trans(um)*g

            IF (inew >= mcnew) THEN
                CALL rwaxv2(n,mcnew,sm((inew-mcnew)*n+1:),&
                    um((inew-mcnew)*n+1:),g,g,tmpmc1(iold:),tmpmc2(iold:))
            ELSE
                CALL rwaxv2(n,inew,sm,um,g,g,tmpmc1,tmpmc2)
                CALL rwaxv2(n,mcnew-inew,sm((iold-1)*n+1:),&
                    um((iold-1)*n+1:),g,g,tmpmc1(iold:),tmpmc2(iold:))
            END IF
            
            
            ! Computation of trans(sm)*u and trans(um)*u

            IF (inew >= mcnew) THEN
                DO i=iold,inew-1
                    tmpmc3(i) = tmpmc1(i) - smtgp(i)
                    smtgp(i)  = tmpmc1(i)
                    tmpmc4(i) = tmpmc2(i) - umtgp(i)
                    umtgp(i)  = tmpmc2(i)
                END DO
            ELSE
                DO i=1,inew-1
                    tmpmc3(i) = tmpmc1(i) - smtgp(i)
                    smtgp(i)  = tmpmc1(i)
                    tmpmc4(i) = tmpmc2(i) - umtgp(i)
                    umtgp(i)  = tmpmc2(i)
                END DO
                DO i=iold,mcnew+1
                    tmpmc3(i) = tmpmc1(i) - smtgp(i)
                    smtgp(i)  = tmpmc1(i)
                    tmpmc4(i) = tmpmc2(i) - umtgp(i)
                    umtgp(i)  = tmpmc2(i)
                END DO
            END IF
            tmpmc3(inew) = tmpmc1(inew) - DOT_PRODUCT(s,gp)
            smtgp(inew)  = tmpmc1(inew)
            tmpmc4(inew) = tmpmc2(inew) - DOT_PRODUCT(u,gp)
            umtgp(inew)  = tmpmc2(inew)
            
         
            ! Update rm and umtum

            IF (mcc >= mc .AND. iflag2 /= 1) THEN
                DO i=1,mcnew-1
                    j=(i-1)*i/2+1
                    k=i*(i+1)/2+2
                    CALL copy2(i,rm(k:),rm(j:),umtum(k:),umtum(j:))
                END DO
            END IF
                      
            IF (inew >= mcnew) THEN
                CALL copy2(mcnew,tmpmc3(iold:),rm((mcnew-1)*mcnew/2+1:),&
                    tmpmc4(iold:),umtum((mcnew-1)*mcnew/2+1:))
            ELSE
                CALL copy2(mcnew-inew,tmpmc3(iold:),rm((mcnew-1)*mcnew/2+1:)&
                    ,tmpmc4(iold:),umtum((mcnew-1)*mcnew/2+1:))
                CALL copy2(inew,tmpmc3,rm((mcnew-1)*mcnew/2+mcnew-inew+1:),&
                    tmpmc4,umtum((mcnew-1)*mcnew/2+mcnew-inew+1:))
            END IF
            

            ! Update CM

            cm(inew) = stu
            
            ! Computation of gamma

            gamma = sclpar(mcc,iscale,sts,stu,tmpmc4(inew))
            
            inew = inew + 1
            IF (inew > mc + 1) inew = 1
            IF (iflag == 0 .AND. mcc < mc + 1) mcc = mcc + 1
          
        ELSE

            
            ! BFGS update, corrections are not saved.
     
            ibfgs = 2

            ! Initialization of indices.

            CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,ibfgs)
          
            ! Update sm and um
          
            CALL copy2(n,s,sm((inew-1)*n+1:),u,um((inew-1)*n+1:))
            
            ! Computation of trans(sm)*g and trans(um)*g

            CALL rwaxv2(n,mcnew,sm,um,g,g,tmpmc1,tmpmc2)

            ! Computation of trans(sm)*u and trans(um)*u
          
            IF (iold /= 1) THEN
                DO i=1,inew-1
                    tmpmc3(i) = tmpmc1(i) - smtgp(i)
                    smtgp(i)  = tmpmc1(i)
                    tmpmc4(i) = tmpmc2(i) - umtgp(i)
                    umtgp(i)  = tmpmc2(i)
                END DO
                DO i=iold,mcnew
                    tmpmc3(i) = tmpmc1(i) - smtgp(i)
                    smtgp(i)  = tmpmc1(i)
                    tmpmc4(i) = tmpmc2(i) - umtgp(i)
                    umtgp(i)  = tmpmc2(i)
                END DO
            ELSE
                DO i=1,mcnew-1
                    tmpmc3(i) = tmpmc1(i) - smtgp(i)
                    smtgp(i)  = tmpmc1(i)
                    tmpmc4(i) = tmpmc2(i) - umtgp(i)
                    umtgp(i)  = tmpmc2(i)
                END DO
            END IF
            tmpmc3(inew) = tmpmc1(inew) - DOT_PRODUCT(s,gp)
            smtgp(inew)  = tmpmc1(inew)
            tmpmc4(inew) = tmpmc2(inew) - DOT_PRODUCT(u,gp)
            umtgp(inew)  = tmpmc2(inew)


            ! Update rm and umtum

            IF (iold /= 1) THEN
                CALL copy2(mcnew-inew,tmpmc3(iold:),&
                    rm((mcnew-1)*mcnew/2+1:),tmpmc4(iold:),&
                    umtum((mcnew-1)*mcnew/2+1:))
                CALL copy2(inew,tmpmc3,rm((mcnew-1)*mcnew/2+mcnew-inew+1:),tmpmc4,&
                    umtum((mcnew-1)*mcnew/2+mcnew-inew+1:))
            ELSE
                CALL copy2(mcnew,tmpmc3,rm((mcnew-1)*mcnew/2+1:)&
                    ,tmpmc4,umtum((mcnew-1)*mcnew/2+1:))
            END IF
            

            ! Update cm

            cm(inew) = stu
            

            ! Computation of gamma

            gamma = sclpar(mcc,iscale,sts,stu,tmpmc4(inew))
               
        END IF
    ELSE
         
     
        ! BFGS update is skipped

        ibfgs = 3

        IF (mcc == 0) THEN
            iflag = 0
            d=-g
            RETURN
        END IF
         

        !     Initialization of indices.

        CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,ibfgs)


        !     Computation of gamma

        IF (iscale >= 4) gamma = one
               
         
        !     Computation of trans(sm)*g and trans(um)*g and the two intermediate values

        IF (iold <= 2) THEN
            CALL rwaxv2(n,mcnew,sm((iold-1)*n+1:),um((iold-1)*n+1:),g,g,&
                smtgp(iold:),umtgp(iold:))
        ELSE
            CALL rwaxv2(n,inew-1,sm,um,g,g,smtgp,umtgp)
            CALL rwaxv2(n,mcnew-inew+1,sm((iold-1)*n+1:),&
                um((iold-1)*n+1:),g,g,smtgp(iold:),umtgp(iold:))
        END IF
    END IF


    ! Computation of two intermediate values tmpmc1 and tmpmc2

    IF (iold == 1 .OR. ibfgs == 2) THEN
        CALL trlieq(mcnew,mcnew,iold,rm,tmpmc1,smtgp,1,ierr)
        CALL symax(mcnew,mcnew,iold,umtum,tmpmc1,tmpmc3)
        CALL vxdiag(mcnew,cm,tmpmc1,tmpmc2)
        CALL scsum(mcnew,gamma,tmpmc3,tmpmc2,tmpmc2)
        CALL scsum(mcnew,-gamma,umtgp,tmpmc2,tmpmc3)
        CALL trlieq(mcnew,mcnew,iold,rm,tmpmc2,tmpmc3,0,ierr)

    ELSE IF (iflag == 0) THEN
        CALL trlieq(mcnew,mc+1,iold,rm,tmpmc1,smtgp,1,ierr)
        CALL symax(mcnew,mc+1,iold,umtum,tmpmc1,tmpmc3)
        CALL vxdiag(mc+1,cm,tmpmc1,tmpmc2)
        CALL scsum(mc+1,gamma,tmpmc3,tmpmc2,tmpmc2)
        CALL scsum(mc+1,-gamma,umtgp,tmpmc2,tmpmc3)
        CALL trlieq(mcnew,mc+1,iold,rm,tmpmc2,tmpmc3,0,ierr)

    ELSE
        CALL trlieq(mcnew,mc,iold,rm,tmpmc1,smtgp,1,ierr)
        CALL symax(mcnew,mc,iold,umtum,tmpmc1,tmpmc3)
        CALL vxdiag(mc,cm,tmpmc1,tmpmc2)
        CALL scsum(mc,gamma,tmpmc3,tmpmc2,tmpmc2)
        CALL scsum(mc,-gamma,umtgp,tmpmc2,tmpmc3)
        CALL trlieq(mcnew,mc,iold,rm,tmpmc2,tmpmc3,0,ierr)
    END IF

      
    ! Computation of the search direction d

    IF (iold == 1 .OR. ibfgs == 2) THEN
        CALL cwmaxv(n,mcnew,um,tmpmc1,d)
        CALL xdiffy(n,d,g,d)
        CALL cwmaxv(n,mcnew,sm,tmpmc2,tmpn1)
        CALL scdiff(n,gamma,d,tmpn1,d)
    ELSE 
        CALL cwmaxv(n,inew-1,um,tmpmc1,d)
        CALL cwmaxv(n,mcnew-inew+1,um((iold-1)*n+1:),tmpmc1(iold:),tmpn1)
        CALL xsumy(n,d,tmpn1,d)
        CALL xdiffy(n,d,g,d)
        CALL cwmaxv(n,inew-1,sm,tmpmc2,tmpn1)
        CALL scdiff(n,gamma,d,tmpn1,d)
        CALL cwmaxv(n,mcnew-inew+1,sm((iold-1)*n+1:),tmpmc2(iold:),tmpn1)
        CALL xdiffy(n,d,tmpn1,d)
    END IF

CONTAINS

    FUNCTION sclpar(mcc,iscale,sts,stu,utu) RESULT(spar) ! Calculation of the scaling parameter.
      
        ! USE param, ONLY : small,one,half ! given in host
        IMPLICIT NONE

        ! Scalar Arguments
        REAL(KIND=prec), INTENT(IN) :: &
            sts, &     ! sts = trans(s)*s.
            stu, &     ! stu = trans(s)*u.
            utu        ! utu = trans(u)*u.
        REAL(KIND=prec) :: &
            spar       ! Scaling parameter.
        INTEGER, INTENT(IN) :: &
            mcc, &     ! Current number of stored corrections.
            iscale     ! Selection of the scaling:
                       !   0  - Scaling at every iteration with STU/UTU.
                       !   1  - Scaling at every iteration with STS/STU.
                       !   2  - Interval scaling with STU/UTU.
                       !   3  - Interval scaling with STS/STU.
                       !   4  - Preliminary scaling with STU/UTU.
                       !   5  - Preliminary scaling with STS/STU.
                       !   6  - No scaling.

        ! Intrinsic Functions
        INTRINSIC SQRT


        ! Computation of scaling parameter.

        SELECT CASE(iscale)
     
            ! Scaling parameter = STU/UTU

            CASE(0,2,4)
                IF (utu < SQRT(small)) THEN
                    spar = one
                    RETURN
                ELSE
                    spar = stu/utu
                END IF
    
            ! Scaling parameter = STS/STU

            CASE(1,3,5)
                IF (stu < SQRT(small)) THEN
                    spar = one
                    RETURN
                ELSE
                    spar = sts/stu
                END IF

            ! No scaling

            CASE DEFAULT
                spar = one
                RETURN
        END SELECT

               
        !     Scaling
            
        IF (MCC == 0) THEN
            IF (spar < 0.01_prec) spar=0.01_prec
            IF (spar > 100.0_prec) spar=100.0_prec

        ELSE

            SELECT CASE(iscale)

                ! Interval scaling
                CASE(2)
                    IF (spar < 0.6_prec .OR. spar > 6.0_prec) spar = one

                CASE(3)
                    IF (spar < half .OR. spar > 5.0_prec) spar = one
               
                ! Preliminary scaling
                CASE(4,5)
                    spar = one

                ! Scaling at every iteration
                CASE DEFAULT
                CONTINUE
        END SELECT

    END IF

    IF (spar < 1.0E+03_prec*small) spar = 1.0E+03_prec*small
         
END FUNCTION sclpar

END SUBROUTINE dlbfgs

      
!************************************************************************
!*                                                                      *
!*     * SUBROUTINE dlsr1 *                                             *
!*                                                                      *
!*     Matrix update and computation of the search direction d = -dm*ga *
!*     by the limited memory SR1 update.                                *
!*                                                                      *
!************************************************************************

SUBROUTINE dlsr1(n,mc,mcc,inew,isr1,iflag,d,gp,ga,s,u,sm,um,rm,&
    umtum,cm,smtgp,umtgp,gamma,tmpmc1,tmpmc2,tmpn1,nnk,iprint)
      
    USE param, ONLY : zero,small,one
    USE subpro, ONLY : &
        scalex, & ! Scaling a vector.
        xdiffy, & ! Difference of two vectors.
        scdiff, & ! Difference of the scaled vector and a vector.
        xsumy, &  ! Sum of two vectors.
        cwmaxv, & ! Multiplication of a vector by a dense rectangular matrix.
        rwaxv2, & ! Multiplication of two rowwise stored dense rectangular
                  ! matrices A and B by vectors x and y.
        calq, &   ! Solving x from linear equation A*x=y.
        copy, &   ! Copying of a vector.
        copy2     ! Copying of two vectors.
    IMPLICIT NONE

    ! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
        ga, &      ! Current aggregate subgradient of the objective function.
        gp         ! Basic subgradient of the objective function.
    REAL(KIND=prec), DIMENSION(:), INTENT(INOUT) :: &
        sm, &      ! Matrix whose columns are stored corrections.
        um, &      ! Matrix whose columns are stored subgradient differences.
        rm, &      ! Upper triangular matrix.
        cm, &      ! Diagonal matrix.
        umtum, &   ! Matrix umtum = trans(um) * um.
        smtgp, &   ! Vector smtgp = trans(sm)*gp.
        umtgp, &   ! vector umtgp = trans(um)*gp.
        s, &       ! Difference of current and previous variables.
        u, &       ! Difference of current and previous subgradients.
        tmpn1      ! Auxiliary array. On input: previous aggregate subgradient.
    REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
        tmpmc1, &  ! Auxiliary array. On output: trans(sm)*ga.
        tmpmc2, &  ! Auxiliary array. On output: trans(um)*ga.
        d          ! Direction vector.

    ! Scalar Arguments
    REAL(KIND=prec), INTENT(INOUT) :: & 
        gamma      ! Scaling parameter.
    INTEGER, INTENT(IN) :: & 
        n, &       ! Number of variables
        mc, &      ! Declared number of stored corrections.
        nnk, &     ! Consecutive null steps counter.
        iprint     ! Printout specification.
    INTEGER, INTENT(INOUT) :: & 
        mcc, &     ! Current number of stored corrections.
        inew, &    ! Index for circular arrays.
        iflag      ! Index for adaptive version:
                   !   0  - Maximum number of stored corrections
                   !        has not been changed at this iteration.
                   !   1  - Maximum number of stored corrections
                   !        has been changed at this iteration.
    INTEGER, INTENT(OUT) :: & 
        isr1       ! Index of the type of L-SR1 update:
                   !   1  - SR1 update: the corrections are stored.
                   !   3  - SR1 update is skipped.
      
    ! Local Arrays
    REAL(KIND=prec), DIMENSION(n) :: tmpn2
    REAL(KIND=prec), DIMENSION((mcc+1)*(mcc+2)/2) :: tmpmat
    REAL(KIND=prec), DIMENSION(mcc+1) :: tmpmc3,tmpmc4,tmpmc5,tmpmc6

    ! Local Scalars
    REAL(KIND=prec) :: &
        stu, &    ! stu = trans(s)*u.
        a, &      ! a = trans(ga) dm_(k-1) ga.
        b         ! b = trans(ga) dm_k ga.
    INTEGER :: i,j,k, &
        mcnew, &  ! Current size of vectors.
        iold, &   ! Index of the oldest corrections.
        iflag2    ! Index for adaptive version.


    iflag2 = 0
    isr1 = 0 
      

    ! Initialization of indices
      
    CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,3)
      
         
    ! Computation of gamma

    gamma = one


    ! Computation of trans(sm)*ga and trans(um)*ga

    IF (iold <= 2) THEN
        CALL rwaxv2(n,mcnew,sm((iold-1)*n+1:),um((iold-1)*n+1:),ga,ga,&
            tmpmc1(iold:),tmpmc2(iold:))
    ELSE
        CALL rwaxv2(n,inew-1,sm,um,ga,ga,tmpmc1,tmpmc2)
        CALL rwaxv2(n,mcnew-inew+1,sm((iold-1)*n+1:),um((iold-1)*n+1:),&
            ga,ga,tmpmc1(iold:),tmpmc2(iold:))
    END IF


    ! Positive definiteness

    IF (-DOT_PRODUCT(d,u) - DOT_PRODUCT(tmpn1,s) >= -small) THEN
      
     
        ! SR1 update is skipped

        isr1 = 3
       
        IF (mcc == 0) THEN
            iflag = 0
            d=-ga
            RETURN
        END IF

    ELSE
      
        stu = DOT_PRODUCT(s,u)
        
        tmpmc1(inew) = DOT_PRODUCT(s,ga)
        tmpmc2(inew) = DOT_PRODUCT(u,ga)
      

        !     Convergence conditions

        IF ((nnk == 1 .OR. mcc < mc) .OR. &
            (iflag == 1 .AND. (inew == 1 .OR. inew == mc))) THEN

            ! SR1 Update

            isr1 = 1
 

            ! Initialization of indices

            CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,1)

            if (iflag2 == 1 .and. iold == 2) then
                tmpmc1(inew) = tmpmc1(1)
                tmpmc2(inew) = tmpmc2(1)
            end if


            ! Update sm and um

            CALL copy2(n,s,sm((inew-1)*n+1:),u,um((inew-1)*n+1:))

         
            ! Update trans(sm)*gp and trans(um)*gp

            smtgp(inew) = DOT_PRODUCT(s,gp)
            umtgp(inew) = DOT_PRODUCT(u,gp)

     
            ! Computation of trans(sm)*u and trans(um)*u

            IF (iold <= 2) THEN
                CALL rwaxv2(n,mcnew-1,sm((iold-1)*n+1:),um((iold-1)*n+1:),&
                    u,u,tmpmc3(iold:),tmpmc4(iold:))
            ELSE
                CALL rwaxv2(n,inew-1,sm,um,u,u,tmpmc3,tmpmc4)
                CALL rwaxv2(n,mcnew-inew,sm((iold-1)*n+1:),um((iold-1)*n+1:),&
                    u,u,tmpmc3(iold:),tmpmc4(iold:))
            END IF

            tmpmc3(inew) = stu
            tmpmc4(inew) = DOT_PRODUCT(u,u)

         
            ! Update rm and umtum

            IF (mcc >= mc .AND. iflag2 /= 1) THEN
                DO i=1,mcnew-1
                    j=(i-1)*i/2+1
                    k=i*(i+1)/2+2
                    CALL copy2(i,rm(k:),rm(j:),umtum(k:),umtum(j:))
                END DO
            END IF


            IF (inew >= mcnew) THEN
                CALL copy2(mcnew,tmpmc3(iold:),rm((mcnew-1)*mcnew/2+1:),&
                    tmpmc4(iold:),umtum((mcnew-1)*mcnew/2+1:))
            ELSE
                CALL copy2(mcnew-inew,tmpmc3(iold:),&
                    rm((mcnew-1)*mcnew/2+1:),tmpmc4(iold:),&
                    umtum((mcnew-1)*mcnew/2+1:))
                CALL COPY2(inew,tmpmc3,rm((mcnew-1)*mcnew/2+mcnew-inew+1:),tmpmc4,&
                    umtum((mcnew-1)*mcnew/2+mcnew-inew+1:))
            END IF

      
            ! Update CM

            cm(inew) = stu

            inew = inew + 1
            IF (inew > mc + 1) inew = 1
            IF (iflag == 0 .AND. mcc < mc + 1) mcc = mcc + 1

        ELSE


            ! Calculation of matrix (umtum-rm-trans(rm)+cm) from previous iteration
         
            DO  i=1,mcnew*(mcnew+1)/2
                tmpmat(i)= gamma * umtum(i) - rm(i)
            END DO

     
            ! Computation of tmpmat*tmpmc4 = gamma*trans(um)*ga-trans(sm)*ga

            IF (iold == 1) THEN
                CALL scdiff(mcnew,gamma,tmpmc2,tmpmc1,tmpmc5)
                CALL calq(mcnew,mcnew,iold,tmpmat,tmpmc4,tmpmc5,0)

            ELSE IF (iflag == 0) THEN
                CALL scdiff(mc+1,gamma,tmpmc2,tmpmc1,tmpmc5)
                CALL calq(mcnew,mc+1,iold,tmpmat,tmpmc4,tmpmc5,0)

            ELSE
                CALL scdiff(mc,gamma,tmpmc2,tmpmc1,tmpmc5)
                CALL calq(mcnew,mc,iold,tmpmat,tmpmc4,tmpmc5,0)
            END IF


            ! Computation of a = -trans(ga)*dm_(k-1)*ga
      
            IF (iold <= 2) THEN
                CALL scalex(mcnew,gamma,tmpmc4(iold:),tmpmc3(iold:))
                CALL cwmaxv(n,mcnew,sm((iold-1)*n+1:),tmpmc4(iold:),tmpn1)
                CALL scdiff(n,-gamma,ga,tmpn1,tmpn2)
                CALL cwmaxv(n,mcnew,um((iold-1)*n+1:),tmpmc3(iold:),tmpn1)
                CALL xsumy(n,tmpn2,tmpn1,tmpn2)
             
            ELSE
                CALL scalex(mcc,gamma,tmpmc4,tmpmc3)
                CALL cwmaxv(n,inew-1,sm,tmpmc4,tmpn1)
                CALL scdiff(n,-gamma,ga,tmpn1,tmpn2)
                CALL cwmaxv(n,mcnew-inew+1,sm((iold-1)*n+1:),tmpmc4(iold:),&
                    tmpn1)
                CALL xdiffy(n,tmpn2,tmpn1,tmpn2)
                CALL cwmaxv(n,inew-1,um,tmpmc3,tmpn1)
                CALL xsumy(n,tmpn2,tmpn1,tmpn2)
                CALL cwmaxv(n,mcnew-inew+1,um((iold-1)*n+1:),tmpmc3(iold:),&
                    tmpn1)
                CALL xsumy(n,tmpn2,tmpn1,tmpn2)
            END IF
          
            a = DOT_PRODUCT(ga,tmpn2)
          
            IF (iflag == 0) THEN
                mcnew = mc
                iold = inew + 2
                IF (iold > mc+1) iold = iold - mc - 1
            ELSE
                mcnew = mc - 1
                iold = inew + 2
                IF (iold > mc) iold = iold - mc
            END IF
      

            ! Calculation of the new canditate for search direction
            ! Updates are not necessarily saved

            ! Update sm and um

            CALL copy2(n,s,sm((inew-1)*n+1:),u,um((inew-1)*n+1:))

     
            ! Computation of trans(sm)*u and trans(um)*u

            IF (iold == 1 .OR. iold == 2) THEN
                CALL rwaxv2(n,mcnew-1,sm((iold-1)*n+1:),um((iold-1)*n+1:),u,u,&
                    tmpmc3(iold:),tmpmc4(iold:))
            ELSE
                CALL rwaxv2(n,inew-1,sm,um,u,u,tmpmc3,tmpmc4)
                CALL rwaxv2(n,mcnew-inew,sm((iold-1)*n+1:),um((iold-1)*n+1:),u,u,&
                    tmpmc3(iold:),tmpmc4(iold:))
            END IF
       
            tmpmc3(inew) = stu
            tmpmc4(inew) = DOT_PRODUCT(u,u)


            ! Calculation of matrix (umtum-rm-trans(rm)+cm) without updating
            ! matrices rm, umtum and cm
      
            DO i=1,mcnew*(mcnew+1)/2
                tmpmat(i)= gamma * umtum(i) - rm(i)
            END DO

            DO i=1,mcnew-1
                j=(i-1)*i/2+1
                k=i*(i+1)/2+2
                CALL copy(i,tmpmat(k:),tmpmat(j:))
            END DO
         
            CALL scdiff(mcnew+1,gamma,tmpmc4,tmpmc3,tmpmc5)
         
            IF (inew >= mcnew) THEN
                CALL copy(mcnew,tmpmc5(iold:),tmpmat((mcnew-1)*mcnew/2+1:))
            ELSE
                CALL copy(mcnew-inew,tmpmc5(iold:),tmpmat((mcnew-1)*mcnew/2+1:))
                CALL copy(inew,tmpmc5,tmpmat((mcnew-1)*mcnew/2+mcnew-inew+1:))
            END IF
      
            IF (iflag == 0) THEN
                CALL scdiff(mc+1,gamma,tmpmc2,tmpmc1,tmpmc5)
                CALL calq(mcnew,mc+1,iold,tmpmat,tmpmc5,tmpmc5,iprint)

            ELSE
                CALL scdiff(mc,gamma,tmpmc2,tmpmc1,tmpmc5)
                CALL calq(mcnew,mc,iold,tmpmat,tmpmc5,tmpmc5,iprint)
            END IF


            ! Calculation of the new canditate for search direction d = -dm_k*ga
            ! and computation of b = -trans(ga)*dm_k*ga
      
            IF (iold <= 2) THEN
                CALL scalex(mcnew,gamma,tmpmc5(iold:),tmpmc6(iold:))
                CALL cwmaxv(n,mcnew,sm((iold-1)*n+1:),tmpmc5(iold:),tmpn1)
                CALL scdiff(n,-gamma,ga,tmpn1,d)
                CALL cwmaxv(n,mcnew,um((iold-1)*n+1:),tmpmc6(iold:),tmpn1)
                CALL xsumy(n,d,tmpn1,d)
         
            ELSE
                CALL scalex(mcnew+1,gamma,tmpmc5,tmpmc6)
                CALL cwmaxv(n,inew,sm,tmpmc5,tmpn1)
                CALL scdiff(n,-gamma,ga,tmpn1,d)
                CALL cwmaxv(n,mcnew-inew,sm((iold-1)*n+1:),tmpmc5(iold:),&
                    tmpn1)
                CALL xdiffy(n,d,tmpn1,d)
                CALL cwmaxv(n,inew,um,tmpmc6,tmpn1)
                CALL xsumy(n,d,tmpn1,d)
                CALL cwmaxv(n,mcnew-inew,um((iold-1)*n+1:),tmpmc6(iold:),&
                    tmpn1)
                CALL xsumy(n,d,tmpn1,d)
            END IF

            b = DOT_PRODUCT(ga,d)


            ! Checking the convergence conditions

            IF (b - a < zero) THEN
                isr1 = 3
                CALL copy(n,tmpn2,d)
            
            ELSE

                isr1 = 1
         
     
                ! Update trans(sm)*gp and trans(um)*gp
     
                smtgp(inew) = DOT_PRODUCT(s,gp)
                umtgp(inew) = DOT_PRODUCT(u,gp)

                     
                ! Update rm and umtum

                DO i=1,mcnew-1
                    j=(i-1)*i/2+1
                    k=i*(i+1)/2+2
                    CALL copy2(i,rm(k:),rm(j:),umtum(k:),umtum(j:))
                END DO

                IF (inew >= mcnew) THEN
                    CALL copy2(mcnew,tmpmc3(iold:),rm((mcnew-1)*mcnew/2+1:),tmpmc4(iold:),&
                        umtum((mcnew-1)*mcnew/2+1:))
                ELSE
                    CALL copy2(mcnew-inew,tmpmc3(iold:),rm((mcnew-1)*mcnew/2+1:),tmpmc4(iold:),&
                        umtum((mcnew-1)*mcnew/2+1:))
                    CALL copy2(inew,tmpmc3,rm((mcnew-1)*mcnew/2+mcnew-inew+1:),tmpmc4,&
                        umtum((mcnew-1)*mcnew/2+mcnew-inew+1:))
                END IF
            

                ! Update cm

                cm(inew) = stu
                     
                inew = inew + 1
                IF (inew > mc + 1) inew = 1
                IF (iflag == 0 .AND. mcc < mc + 1) mcc = mcc + 1
            
            END IF
         
            RETURN
         
        END IF
    END IF
      
    DO i=1,mcnew*(mcnew+1)/2
        tmpmat(i)= gamma * umtum(i) - rm(I)
    END DO
      
     
    ! Computation of tmpmat*tmpmc4 = gamma*trans(um)*ga-trans(sm)*ga

    IF (iold == 1) THEN
        CALL scdiff(mcnew,gamma,tmpmc2,tmpmc1,tmpmc4)
        CALL calq(mcnew,mcnew,iold,tmpmat,tmpmc4,tmpmc4,iprint)
    ELSE IF (iflag == 0) THEN
        CALL scdiff(mc+1,gamma,tmpmc2,tmpmc1,tmpmc4)
        CALL calq(mcnew,mc+1,iold,tmpmat,tmpmc4,tmpmc4,iprint)
    ELSE
        CALL scdiff(mc,gamma,tmpmc2,tmpmc1,tmpmc4)
        CALL calq(mcnew,mc,iold,tmpmat,tmpmc4,tmpmc4,iprint)
    END IF
      

    ! Computation of the search direction d
      
    IF (iold <= 2) THEN
        CALL scalex(mcnew,gamma,tmpmc4(iold:),tmpmc3(iold:))
        CALL cwmaxv(n,mcnew,sm((iold-1)*n+1:),tmpmc4(iold:),tmpn1)
        CALL scdiff(n,-gamma,ga,tmpn1,d)
        CALL cwmaxv(n,mcnew,um((iold-1)*n+1:),tmpmc3(iold:),tmpn1)
        CALL xsumy(n,d,tmpn1,d)
    ELSE
        CALL scalex(mcc,gamma,tmpmc4,tmpmc3)
        CALL cwmaxv(n,inew-1,sm,tmpmc4,tmpn1)
        CALL scdiff(n,-gamma,ga,tmpn1,d)
        CALL cwmaxv(n,mcnew-inew+1,sm((iold-1)*n+1:),tmpmc4(iold:),&
            tmpn1)
        CALL xdiffy(n,d,tmpn1,d)
        CALL cwmaxv(n,inew-1,um,tmpmc3,tmpn1)
        CALL xsumy(n,d,tmpn1,d)
        CALL cwmaxv(n,mcnew-inew+1,um((iold-1)*n+1:),tmpmc3(iold:),&
            tmpn1)
        CALL xsumy(n,d,tmpn1,d)
    END IF
      
END SUBROUTINE dlsr1

!************************************************************************
!*                                                                      *
!*     * SUBROUTINE indic1 *                                            *
!*                                                                      *
!*     Initialization of indices.                                       *
!*                                                                      *
!************************************************************************

SUBROUTINE indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,itype)

    IMPLICIT NONE

    ! Scalar Arguments

    INTEGER, INTENT(IN) :: & 
        mc, &     ! Declared number of stored corrections.
        mcc, &    ! Current number of stored corrections.
        itype     ! Type of Initialization:
                   !   1  - corrections are stored,
                   !   2  - corrections are not stored,
                   !   3  - update is skipped.
    INTEGER, INTENT(INOUT) :: & 
        inew, &   ! Index for circular arrays.
        iflag, &  ! Index for adaptive version:
                  !   0  - Maximum number of stored corrections
                  !        has not been changed at this iteration.
                  !   1  - Maximum number of stored corrections
                  !        has been changed at this iteration.
        iflag2    ! Index for adaptive version.
                   !   0  - iflag has not been changed.
                   !   1  - iflag has been changed.
    INTEGER, INTENT(OUT) :: & 
        mcnew, &  ! Current size of vectors.
        iold      ! Index of the oldest corrections.
      
    IF (itype == 1) THEN
        IF (mcc < mc) THEN
            mcnew = mcc + 1
            iold = 1
            iflag = 0
        ELSE
            IF (iflag == 0) THEN
                mcnew = mc
                iold = inew + 2
                IF (iold > mc+1) iold = iold - mc - 1
            ELSE
                IF (inew == 1) THEN
                    inew = mc + 1
                    mcnew = mc
                    iold = 2
                    iflag = 0
                    iflag2 = 1
                ELSE IF (inew == mc) THEN
                    mcnew = mc
                    iold = 1
                    iflag = 0
                    iflag2 = 1
                ELSE
                    mcnew = mc - 1
                    iold = inew + 2
                    IF (iold > mc) iold = iold - mc
                END IF
            END IF
        END IF
      
    ELSE IF (itype == 2) THEN
        IF (mcc < mc) THEN
            mcnew = mcc + 1
            iold = 1
            iflag = 0
        ELSE
            IF (iflag == 0) THEN
                mcnew = mc + 1
                iold = inew + 1
                IF (iold > mc + 1) iold = 1
            ELSE
                mcnew = mc
                iold = inew + 1
                IF (iold > mc) iold = 1
            END IF
        END IF

    ELSE 
        IF (mcc < mc) THEN
            mcnew = mcc
            iold = 1
            iflag = 0
        ELSE
            IF (iflag == 0) THEN
                mcnew = mc
                iold = inew + 1
                IF (iold > mc + 1) iold = 1
            ELSE
                mcnew = mc - 1
                iold = inew + 1
                IF (iold > mc) iold = 1
            END IF
        END IF
    END IF
END SUBROUTINE indic1
  
!************************************************************************
!*
!*     * SUBROUTINE agbfgs *
!*
!*     Computation of aggregate values by the limited memory BFGS update.
!*
!************************************************************************

SUBROUTINE agbfgs(n,mc,mcc,inew,ibfgs,iflag,g,gp,ga,u,d,sm,um, &
    rm,cm,umtum,alfn,alfv,gamma,ic,rho)

    USE param, ONLY : zero,half,one
    USE subpro, ONLY : &
        symax, &  ! Multiplication of a dense symmetric matrix by a vector.
        rwaxv2, & ! Multiplication of two rowwise stored dense rectangular
                  ! matrices A and B by vectors x and y.
        trlieq, & ! Solving x from linear equation L*x=y or trans(L)*x=y.
        vdot      ! Dot product of vectors

    IMPLICIT NONE

    ! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
        d, &      ! Direction vector.
        g, &      ! Current (auxiliary) subgradient of the objective function.
        gp, &     ! Previous subgradient of the objective function.
        u, &      ! Difference of trial and aggregate gradients.
        sm, &     ! Matrix whose columns are stored corrections.
        um, &     ! Matrix whose columns are stored subgradient differences.
        rm, &     ! Upper triangular matrix.
        umtum, &  ! Matrix umtum = trans(um) * um.
        cm        ! Diagonal matrix.
    REAL(KIND=prec), DIMENSION(:), INTENT(OUT) :: &
        ga        ! Next aggregate subgradient of the objective function.

    ! Scalar Arguments
    REAL(KIND=prec), INTENT(OUT) :: & 
        alfv      ! Aggregate locality measure.
    REAL(KIND=prec), INTENT(IN) :: & 
        gamma, &  ! Scaling parameter.
        alfn, &   ! Locality measure.
        rho       ! Correction parameter.
    INTEGER, INTENT(IN) :: & 
        n, &      ! Number of variables
        mc, &     ! Declared number of stored corrections.
        mcc, &    ! Current number of stored corrections.
        inew, &   ! Index for circular arrays.
        ibfgs, &  ! Index of the type of BFGS update.
        ic        ! Correction indicator.
    INTEGER, INTENT(INOUT) :: & 
        iflag     ! Index for adaptive version:
                   !   0  - Maximum number of stored corrections
                   !        has not been changed at this iteration.
                   !   1  - Maximum number of stored corrections
                   !        has been changed at this iteration.

    ! Local arrays
    REAL(KIND=prec), DIMENSION(mcc+1) :: tmpmc1, tmpmc2

    ! Local Scalars
    REAL(KIND=prec) :: &
        p, &      ! p = trans(d)*u - alfn.
        q, &      ! q = trans(u)*dm*u, where dm is the inverse approximation of
                  ! the Hessian calculated by using the L-BFGS formula.
        lam, &    ! Multiplier used to calculate aggregate values.
        w         ! Correction.
    INTEGER :: i, &
        mcnew, &  ! Current size of vectors.
        iold, &   ! Index of the oldest corrections.
        ierr      ! Error indicador.

    ! Intrinsic Functions
    INTRINSIC MAX,MIN,SIGN

    ierr = 0

    IF (mcc < mc) THEN
        IF (ibfgs == 2) THEN
            mcnew = mcc + 1
        ELSE
            mcnew = mcc
        END IF
        iold = 1

    ELSE
        IF (iflag == 0) THEN
            IF (ibfgs == 2) THEN
                mcnew = mc + 1
            ELSE
                mcnew = mc
            END IF
            iold = inew + 1
            IF (iold > mc+1) iold = 1

        ELSE
            IF (ibfgs == 2) THEN
                mcnew = mc
            ELSE
                mcnew = mc - 1
            END IF
            iold = inew + 1
            IF (iold > mc) iold = 1
        END IF
    END IF
      
      
    ! Computation of trans(d) * u - alfn

    p = DOT_PRODUCT(d,u) - alfn
    q = DOT_PRODUCT(u,u)

    IF (ic == 1) THEN
        w = rho * q
    ELSE
        w = zero
    END IF
         
     
    ! Computation of the product trans(u)*dm*u

    IF (mcc > 0 .OR. ibfgs == 2) THEN

        IF (iold == 1 .OR. ibfgs == 2) THEN
            CALL rwaxv2(n,mcnew,sm,um,u,u,tmpmc1,tmpmc2)
            CALL trlieq(mcnew,mcnew,iold,rm,tmpmc1,tmpmc1,1,ierr)

            q = q - 2.0_prec*vdot(mcnew,tmpmc2,tmpmc1)
            q = gamma*q
            
            DO i=1,mcnew
                tmpmc2(i) = cm(i)*tmpmc1(i)
            END DO

            q = q + vdot(mcnew,tmpmc1,tmpmc2)

            CALL symax(mcnew,mcnew,iold,umtum,tmpmc1,tmpmc2)

            q = q + gamma*vdot(mcnew,tmpmc1,tmpmc2)

        ELSE
            CALL rwaxv2(n,inew-1,sm,um,u,u,tmpmc1,tmpmc2)
            CALL rwaxv2(n,mcc-inew,sm((iold-1)*n+1:),um((iold-1)*n+1:), &
                u,u,tmpmc1(iold:),tmpmc2(iold:))
            CALL trlieq(mcnew,mcc,iold,rm,tmpmc1,tmpmc1,1,ierr)

            q = q - 2.0_prec*(vdot(mcc-inew,tmpmc2(iold:),tmpmc1(iold:)) + &
                vdot(inew-1,tmpmc2,tmpmc1))
            q = gamma*q

            DO i=1,mcc
                tmpmc2(i) = cm(i)*tmpmc1(i)
            END DO

            q = q + vdot(mcc-inew,tmpmc1(iold:),tmpmc2(iold:)) + &
                vdot(inew-1,tmpmc1,tmpmc2)

            CALL symax(mcnew,mcc,iold,umtum,tmpmc1,tmpmc2)

            q = q + gamma*(vdot(mcc-inew,tmpmc1(iold:),tmpmc2(iold:)) + &
                vdot(inew-1,tmpmc1,tmpmc2))
        END IF

    END IF
    
    q = q + w
    
    lam = half + SIGN(half,p)

    IF (q > zero) lam = MIN(one,MAX(zero,p/q))
      

    ! Computation of the aggregate values

    p = one - lam
    DO i=1,n
        ga(i)=lam*g(i) + p*gp(i)
    END DO
      
    alfv = lam*alfn
      
END SUBROUTINE agbfgs

!************************************************************************
!*
!*     * SUBROUTINE aggsr1 *
!*
!*     Computation of aggregate values by the limited memory SR1 update.
!*
!************************************************************************
      
SUBROUTINE aggsr1(n,mc,mcc,inew,iflag,g,gp,ga,d,alfn,alfv, &
    umtum,rm,gamma,smtgp,umtgp,smtga,umtga,sm,um,icn,rho)
      
    USE param, ONLY : zero,one,small
    USE subpro, ONLY : &
        vdot, &   ! Dot product.
        scalex, & ! Scaling a vector.
        xsumy, &  ! Sum of two vectors.
        xdiffy, & ! Difference of two vectors.
        scsum, &  ! Sum of a vector and the scaled vector.
        scdiff, & ! Difference of the scaled vector and a vector.
        rwaxv2, & ! Multiplication of two rowwise stored dense rectangular
                  ! matrices A and B by vectors x and y.
        cwmaxv, & ! Multiplication of a vector by a dense rectangular matrix.
        lineq, &  ! Solver for linear equation.
        calq      ! Solving x from linear equation A*x=y.
    IMPLICIT NONE

    ! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
        d, &      ! Direction vector.
        g, &      ! Current (auxiliary) subgradient of the objective function.
        gp, &     ! Previous subgradient of the objective function.
        sm, &     ! Matrix whose columns are stored corrections.
        um, &     ! Matrix whose columns are stored subgradient differences.
        rm, &     ! Upper triangular matrix.
        umtum, &  ! Matrix umtum = trans(um) * um.
        smtgp, &  ! Vector smtgp = trans(sm)*gp.
        umtgp, &  ! vector umtgp = trans(um)*gp.
        smtga, &  ! vector smtga = trans(sm)*ga.
        umtga     ! vector umtga = trans(um)*ga.
    REAL(KIND=prec), DIMENSION(:), INTENT(INOUT) :: &
        ga        ! Aggregate subgradient of the objective function.

    ! Scalar Arguments
    REAL(KIND=prec), INTENT(INOUT) :: & 
        alfv      ! Aggregate locality measure.
    REAL(KIND=prec), INTENT(IN) :: & 
        gamma, &  ! Scaling parameter.
        alfn, &   ! Locality measure.
        rho       ! Correction parameter.
    INTEGER, INTENT(IN) :: & 
        n, &      ! Number of variables
        mc, &     ! Declared number of stored corrections.
        mcc, &    ! Current number of stored corrections.
        inew, &   ! Index for circular arrays.
        icn       ! Correction indicator.
    INTEGER, INTENT(INOUT) :: & 
        iflag     ! Index for adaptive version:
                   !   0  - Maximum number of stored corrections
                   !        has not been changed at this iteration.
                   !   1  - Maximum number of stored corrections
                   !        has been changed at this iteration.
      
    ! Local arrays
    REAL(KIND=prec), DIMENSION(n) :: tmpn2,tmpn3,tmpn4
    REAL(KIND=prec), DIMENSION((mcc+1)*(mcc)/2) :: tmpmat
    REAL(KIND=prec), DIMENSION(mcc+1) :: tmpmc3, tmpmc4

    ! Local Scalars
    REAL(KIND=prec) :: &
        pr, &     ! pr = trans(gp-ga) dm (gp-ga), where dm
                  ! presents the L-SR1- approximation of Hessian.
        rrp, &    ! rrp = trans(gp-ga) dm ga - alfv.
        prqr, &   ! prqr = trans(gp-ga) dm (g-ga).
        rrq, &    ! rrq = trans(g-ga) dm ga - alfv + alfn.
        qr, &     ! qr = trans(g-ga) dm (g-ga).
        pq, &     ! pq = trans(g-gp) dm (g-gp).
        qqp, &    ! qqp = trans(g-gp) dm g + alfn.
        lam1, &   ! Multiplier used to calculate aggregate values.
        lam2, &   ! Multiplier used to calculate aggregate values.
        w, &      ! Correction.
        tmp1, &   ! Auxiliary scalar.
        tmp2      ! Auxiliary scalar.
    INTEGER :: i, &
        mcnew, &  ! Current size of vectors.
        iold, &   ! Index of the oldest corrections.
        ierr      ! Error indicador.

      
    ! Intrinsic Functions
    INTRINSIC MIN,MAX,DOT_PRODUCT


    ierr = 0
      
    IF (mcc < mc) THEN
        iold = 1
        mcnew = mcc
    ELSE IF (iflag == 0) THEN
        mcnew = mc
        iold = inew + 1
        IF (iold > mc+1) iold = 1
    ELSE
        mcnew = mc - 1
        iold = inew + 1
        IF (iold > mc) iold = 1
    END IF
      
    CALL xdiffy(n,gp,ga,tmpn2)
     
      
    ! Calculation of tmpn3 = trans(gp - ga)dm

    IF (mcc > 0) THEN

        DO i=1,mcnew*(mcnew+1)/2
            tmpmat(i)= gamma * umtum(i) - rm(i)
        END DO

        IF (iold == 1) THEN
            CALL xdiffy(mcnew,umtgp,umtga,tmpmc4)
            CALL scdiff(mcnew,gamma,tmpmc4,smtgp,tmpmc4)
            CALL xsumy(mcnew,tmpmc4,smtga,tmpmc4)
          
            CALL calq(mcnew,mcnew,iold,tmpmat,tmpmc3,tmpmc4,0)
            CALL scalex(mcnew,gamma,tmpmc3,tmpmc4)
          
            CALL cwmaxv(n,mcnew,sm,tmpmc3,tmpn4)
            CALL scsum(n,gamma,tmpn2,tmpn4,tmpn3)
            CALL cwmaxv(n,mcnew,um,tmpmc4,tmpn4)
            CALL xdiffy(n,tmpn3,tmpn4,tmpn3)

        ELSE
            CALL xdiffy(mcc,umtgp,umtga,tmpmc4)
            CALL scdiff(mcc,gamma,tmpmc4,smtgp,tmpmc4)
            CALL xsumy(mcc,tmpmc4,smtga,tmpmc4)

            CALL calq(mcnew,mcc,iold,tmpmat,tmpmc3,tmpmc4,0)
            CALL scalex(mcc,gamma,tmpmc3,tmpmc4)
          
            CALL cwmaxv(n,inew-1,sm,tmpmc3,tmpn4)
            CALL scsum(n,gamma,tmpn2,tmpn4,tmpn3)
            CALL cwmaxv(n,mcnew-inew+1,sm((iold-1)*n+1:),tmpmc3(iold:)&
                ,tmpn4)
            CALL xsumy(n,tmpn3,tmpn4,tmpn3)
            CALL cwmaxv(n,inew-1,um,tmpmc4,tmpn4)
            CALL xdiffy(n,tmpn3,tmpn4,tmpn3)
            CALL cwmaxv(n,mcnew-inew+1,um((iold-1)*n+1:),tmpmc4(iold:)&
                ,tmpn4)
            CALL xdiffy(n,tmpn3,tmpn4,tmpn3)
        END IF
       
        IF (icn == 1) THEN
            CALL scsum(n,rho,tmpn2,tmpn3,tmpn3)
        END IF
         
        pr = DOT_PRODUCT(tmpn3,tmpn2)
        rrp = DOT_PRODUCT(tmpn3,ga)
        CALL xdiffy(n,g,ga,tmpn4)
        prqr = DOT_PRODUCT(tmpn3,tmpn4)
        rrq = -DOT_PRODUCT(tmpn4,d)

    ELSE

        pr = DOT_PRODUCT(tmpn2,tmpn2)
        rrp = DOT_PRODUCT(tmpn2,ga)
        CALL xdiffy(n,g,ga,tmpn4)
        prqr = DOT_PRODUCT(tmpn2,tmpn4)
        rrq = -DOT_PRODUCT(tmpn4,d)
    END IF

    ! calculation of qr = trans(g - ga) dm (g - ga)

    qr = DOT_PRODUCT(tmpn4,tmpn4)
    IF (icn == 1) THEN
        w = rho*qr
    ELSE
        w = zero
    END IF
      
    IF (mcc > 0) THEN
        qr = gamma*qr

        IF (iold == 1) THEN
            CALL rwaxv2(n,mcnew,sm,um,tmpn4,tmpn4,tmpmc4,tmpmc3)
            CALL scsum(mcnew,-gamma,tmpmc3,tmpmc4,tmpmc4)
            CALL lineq(mcnew,mcnew,iold,tmpmat,tmpmc3,tmpmc4,ierr)
            
            qr = qr - vdot(mcnew,tmpmc4,tmpmc3) + w

        ELSE
            CALL rwaxv2(n,inew-1,sm,um,tmpn4,tmpn4,tmpmc4,tmpmc3)
            CALL rwaxv2(n,mcnew-inew+1,sm((iold-1)*n+1:),&
                um((iold-1)*n+1:),tmpn4,tmpn4,tmpmc4(iold:),tmpmc3(iold:))
            CALL scsum(mcc,-gamma,tmpmc3,tmpmc4,tmpmc4)
            CALL lineq(mcnew,mcc,iold,tmpmat,tmpmc3,tmpmc4,ierr)
          
            qr = qr - vdot(mcc-inew,tmpmc4(iold:),tmpmc3(iold:)) -&
                vdot(inew-1,tmpmc4,tmpmc3) + w
        END IF

    END IF
      
    pq = qr - prqr - prqr + pr
    qqp = pq + prqr + rrq - pr - rrp + alfn
    rrp = rrp - alfv
    rrq = rrq + alfn - alfv

    ! computation of multipliers lam1 and lam2

    IF (pr > zero .AND. qr > zero) THEN
        tmp1 = rrq/qr
        tmp2 = prqr/qr
        w = pr - prqr*tmp2
        IF (w /= zero) THEN
            lam1 = (tmp1*prqr - rrp)/w
            lam2 = -tmp1 - lam1*tmp2
            IF (lam1*(lam1 - one) < zero .AND. &
                lam2*(lam1 + lam2 - one) < zero) GO TO 200
        END IF
    END IF

! Minimum on the boundary

100 continue
    lam1 = zero
    lam2 = zero
    IF (alfn <= alfv) lam2 = one
    IF (qr > zero) lam2 = MIN(one,MAX(zero,-rrq/qr))
    w = (lam2*qr + rrq+rrq)*lam2
    !    w = (lam2*qr + 2.0_prec*rrq)*lam2
    tmp1 = zero
    IF (alfv >= zero) tmp1 = one
    IF (pr > zero) tmp1 = MIN(one,MAX(zero,-rrp/pr))
    !    tmp2 = (tmp1*pr + 2.0_prec*rrp)*tmp1
    tmp2 = (tmp1*pr + rrp+rrp)*tmp1
    IF (tmp2 < w) THEN
        w = tmp2
        lam1 = tmp1
        lam2 = zero
    END IF
    
    IF (qqp*(qqp - pq) < zero) THEN
        IF (qr + rrq + rrq - qqp*qqp/pq < W) THEN
            lam1 = qqp/pq
            lam2 = one - lam1
        END IF
    END IF
    
200 CONTINUE
    IF (lam1 == zero .AND. lam2*(lam2 - one) < zero &
        .AND. -rrp - lam2*prqr > zero .AND. pr > zero) &
        lam1 = MIN(one - lam2, (-rrp-lam2*prqr)/pr)

    ! Computation of the aggregate values
      
    tmp1 = one - lam1 - lam2
    DO i=1,n
        ga(i)=lam1*gp(i)+lam2*g(i)+tmp1*ga(i)
    END DO
    
    alfv = lam2*alfn + tmp1*alfv
    
END SUBROUTINE aggsr1
      
!************************************************************************
!*
!*     * SUBROUTINE agskip *
!*
!*     Computation of aggregate values after consecutive null steps
!*     by the limited memory BFGS update.
!*
!************************************************************************
      
SUBROUTINE agskip(n,mc,mcc,inew,iflag,g,gp,ga,d,u,alfn,alfv, &
    umtum,rm,cm,gamma,smtgp,umtgp,smtga,umtga,sm,um,icn,rho)
      
    USE param, ONLY : zero,half,one,small
    USE subpro, ONLY : &
        xdiffy, & ! Difference of two vectors.
        symax, &  ! Multiplication of a dense symmetric matrix by a vector.
        rwaxv2, & ! Multiplication of two rowwise stored dense rectangular
                  ! matrices A and B by vectors x and y.
        trlieq, & ! Solving x from linear equation l*x=y or trans(l)*x=y.
        vdot      ! Dot product
    IMPLICIT NONE

    ! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
        d, &      ! Direction vector.
        g, &      ! Current (auxiliary) subgradient of the objective function.
        gp, &     ! Previous subgradient of the objective function.
        u, &      ! Difference of trial and aggregate gradients.
        sm, &     ! Matrix whose columns are stored corrections.
        um, &     ! Matrix whose columns are stored subgradient differences.
        rm, &     ! Upper triangular matrix.
        cm, &     ! Diagonal matrix.
        umtum, &  ! Matrix umtum = trans(um) * um.
        smtgp, &  ! Vector smtgp = trans(sm)*gp.
        umtgp, &  ! vector umtgp = trans(um)*gp.
        smtga, &  ! vector smtga = trans(sm)*ga.
        umtga     ! vector umtga = trans(um)*ga.
    REAL(KIND=prec), DIMENSION(:), INTENT(INOUT) :: &
        ga        ! Aggregate subgradient of the objective function.

    ! Scalar Arguments
    REAL(KIND=prec), INTENT(INOUT) :: & 
        alfv      ! Aggregate locality measure.
    REAL(KIND=prec), INTENT(IN) :: & 
        gamma, &  ! Scaling parameter.
        alfn, &   ! Locality measure.
        rho       ! Correction parameter.
    INTEGER, INTENT(IN) :: & 
        n, &      ! Number of variables
        mc, &     ! Declared number of stored corrections.
        mcc, &    ! Current number of stored corrections.
        inew, &   ! Index for circular arrays.
        icn       ! Correction indicator.
    INTEGER, INTENT(INOUT) :: & 
        iflag     ! Index for adaptive version:
                   !   0  - Maximum number of stored corrections
                   !        has not been changed at this iteration.
                   !   1  - Maximum number of stored corrections
                   !        has been changed at this iteration.
      
    ! Local arrays
    REAL(KIND=prec), DIMENSION(n) :: tmpn2
    REAL(KIND=prec), DIMENSION(mcc+1) :: tmpmc3, tmpmc4

    ! Local Scalars
    REAL(KIND=prec) :: &
        pr, &     ! pr = trans(gp-ga) dm (gp-ga), where dm
                  ! presents the L-SR1- approximation of Hessian.
        rrp, &    ! rrp = trans(gp-ga) dm ga - alfv.
        prqr, &   ! prqr = trans(gp-ga) dm (g-ga).
        rrq, &    ! rrq = trans(g-ga) dm ga - alfv + alfn.
        qr, &     ! qr = trans(g-ga) dm (g-ga).
        pq, &     ! pq = trans(g-gp) dm (g-gp).
        qqp, &    ! qqp = trans(g-gp) dm g + alfn.
        lam1, &   ! Multiplier used to calculate aggregate values.
        lam2, &   ! Multiplier used to calculate aggregate values.
        w, &      ! Correction.
        tmp1, &   ! Auxiliary scalar.
        tmp2      ! Auxiliary scalar.
    INTEGER :: i, &
        mcnew, &  ! Current size of vectors.
        iold, &   ! Index of the oldest corrections.
        ierr      ! Error indicador.
    
    ! Intrinsic Functions
    INTRINSIC MIN,MAX,DOT_PRODUCT

    ierr = 0
           
    IF (mcc < mc) THEN
        iold = 1
        mcnew = mcc
    ELSE
        IF (iflag == 0) THEN
            mcnew = mc
            iold = inew + 1
            IF (iold > mc+1) iold = 1
        ELSE
            mcnew = mc - 1
            iold = inew + 1
            IF (iold > mc) iold = 1
        END IF
    END IF

      
    ! Calculation of pq = trans(g-gp) dm (g-gp) = trans(u) dm u.

    pq = DOT_PRODUCT(u,u)

    IF (icn == 1) THEN
        w = rho * pq
    ELSE
        w = zero
    END IF
    
    IF (mcc > 0) THEN

        IF (iold == 1) THEN
            CALL rwaxv2(n,mcnew,sm,um,u,u,tmpmc3,tmpmc4)
            CALL trlieq(mcnew,mcnew,iold,rm,tmpmc3,tmpmc3,1,ierr)
          
            pq = pq - 2.0_prec*vdot(mcnew,tmpmc4,tmpmc3)
            pq = gamma*pq
          
            DO i=1,mcnew
                tmpmc4(i) = cm(i)*tmpmc3(i)
            END DO

            pq = pq + vdot(mcnew,tmpmc3,tmpmc4)
          
            CALL symax(mcnew,mcnew,iold,umtum,tmpmc3,tmpmc4)
          
            pq = pq + gamma*vdot(mcnew,tmpmc3,tmpmc4)

        ELSE
            CALL rwaxv2(n,inew-1,sm,um,u,u,tmpmc3,tmpmc4)
            CALL rwaxv2(n,mcc-inew,sm((iold-1)*n+1:),um((iold-1)*n+1:),&
                u,u,tmpmc3(iold:),tmpmc4(iold:))
            CALL trlieq(mcnew,mcc,iold,rm,tmpmc3,tmpmc3,1,ierr)

            pq = pq - 2.0_prec*(vdot(mcc-inew,tmpmc4(iold:),tmpmc3(iold:)) + &
                vdot(inew-1,tmpmc4,tmpmc3))
            pq = gamma*pq
          
            DO i=1,mcc
                tmpmc4(i) = cm(i)*tmpmc3(i)
            END DO

            pq = pq + vdot(mcc-inew,tmpmc3(iold:),tmpmc4(iold:)) + &
                vdot(inew-1,tmpmc3,tmpmc4)

            CALL symax(mcnew,mcc,iold,umtum,tmpmc3,tmpmc4)
          
            pq = pq + gamma*(vdot(mcc-inew,tmpmc3(iold:),tmpmc4(iold:)) &
                + vdot(inew-1,tmpmc3,tmpmc4))
        END IF

    END IF

    pq = pq + w
      

    ! Calculation of pr = trans(gp-ga) dm (gp-ga).
      
    CALL xdiffy(n,gp,ga,tmpn2)
    pr = DOT_PRODUCT(tmpn2,tmpn2)

    IF (icn == 1) THEN
        w = rho * pr
    ELSE
        w = zero
    END IF

    IF (mcc > 0) THEN
       
        IF (iold == 1) THEN
            DO i=1, mcnew
                tmpmc3(i)=smtgp(i)-smtga(i)
                tmpmc4(i)=umtgp(i)-umtga(i)
            END DO
            CALL trlieq(mcnew,mcnew,iold,rm,tmpmc3,tmpmc3,1,ierr)
             
            pr = pr - 2.0_prec*vdot(mcnew,tmpmc4,tmpmc3)
            pr = gamma*pr
            
            DO i=1,mcnew
                tmpmc4(i) = cm(i)*tmpmc3(i)
            END DO

            pr = pr + vdot(mcnew,tmpmc3,tmpmc4)

            CALL symax(mcnew,mcnew,iold,umtum,tmpmc3,tmpmc4)

            pr = pr + gamma*vdot(mcnew,tmpmc3,tmpmc4)

        ELSE
            DO i=1, mcc
                tmpmc3(i)=smtgp(i)-smtga(i)
                tmpmc4(i)=umtgp(i)-umtga(i)
            END DO
            CALL trlieq(mcnew,mcc,iold,rm,tmpmc3,tmpmc3,1,ierr)

            pr = pr - 2.0_prec*(vdot(mcc-inew,tmpmc4(iold:),tmpmc3(iold:)) + &
                vdot(inew-1,tmpmc4,tmpmc3))
            pr = gamma*pr

            DO  i=1,mcc
                tmpmc4(i) = cm(i)*tmpmc3(i)
            END DO

            pr = pr + vdot(mcc-inew,tmpmc3(iold:),tmpmc4(iold:)) + &
                vdot(inew-1,tmpmc3,tmpmc4)

            CALL symax(mcnew,mcc,iold,umtum,tmpmc3,tmpmc4)

            pr = pr + gamma*(vdot(mcc-inew,tmpmc3(iold:),tmpmc4(iold:)) &
                + vdot(inew-1,tmpmc3,tmpmc4))
        END IF

    END IF

    pr = pr + w

      
    ! Calculation of rrp = trans(gp-ga) dm ga - alfv.
      
    rrp = - DOT_PRODUCT(tmpn2,d) - alfv
      

    ! Calculation of qr = trans(g-ga) dm (g-ga).

    CALL xdiffy(n,g,ga,tmpn2)
    qr = DOT_PRODUCT(tmpn2,tmpn2)

    IF (icn == 1) THEN
        w = rho * qr
    ELSE
        w = zero
    END IF

    IF (mcc > 0) THEN

        IF (iold == 1) THEN
            CALL rwaxv2(n,mcnew,sm,um,tmpn2,tmpn2,tmpmc3,tmpmc4)
            CALL trlieq(mcnew,mcnew,iold,rm,tmpmc3,tmpmc3,1,ierr)

            qr = qr - 2.0_prec*vdot(mcnew,tmpmc4,tmpmc3)
            qr = gamma*qr
            
            DO i=1,mcnew
                tmpmc4(i) = cm(i)*tmpmc3(i)
            END DO

            qr = qr + vdot(mcnew,tmpmc3,tmpmc4)
          
            CALL symax(mcnew,mcnew,iold,umtum,tmpmc3,tmpmc4)
          
            qr = qr + gamma*vdot(mcnew,tmpmc3,tmpmc4)

        ELSE
            CALL rwaxv2(n,inew-1,sm,um,tmpn2,tmpn2,tmpmc3,tmpmc4)
            CALL rwaxv2(n,mcc-inew,sm((iold-1)*n+1:),um((iold-1)*n+1:),&
                tmpn2,tmpn2,tmpmc3(iold:),tmpmc4(iold:))
            CALL trlieq(mcnew,mcc,iold,rm,tmpmc3,tmpmc3,1,ierr)

            qr = qr - 2.0_prec*(vdot(mcc-inew,tmpmc4(iold:),tmpmc3(iold:)) + &
                vdot(inew-1,tmpmc4,tmpmc3))
            qr = gamma*qr
          
            DO i=1,mcc
                tmpmc4(i) = cm(i)*tmpmc3(i)
            END DO
          
            qr = qr + vdot(mcc-inew,tmpmc3(iold:),tmpmc4(iold:)) + &
                vdot(inew-1,tmpmc3,tmpmc4)
          
            CALL symax(mcnew,mcc,iold,umtum,tmpmc3,tmpmc4)
          
            qr = qr + gamma*(vdot(mcc-inew,tmpmc3(iold:),tmpmc4(iold:)) &
                +vdot(inew-1,tmpmc3,tmpmc4))
        END IF
       
    END IF
    
    qr = qr + w
      

    ! Calculation of rrq = trans(g-ga) dm ga - alfv + alfn.

    rrq = - DOT_PRODUCT(tmpn2,d) - alfv + alfn

     
    ! Calculation of prqr = trans(gp-ga) dm (g-ga).
      
    prqr = half*(qr - pq + pr)

     
    ! Calculation of qqp = trans(g-gp) dm g + alfn.

    qqp = pq + prqr + rrq - pr - rrp

     
    ! Computation of multipliers lam1 and lam2

    IF (pr > zero .AND. qr > zero) THEN
        tmp1 = rrq/qr
        tmp2 = prqr/qr
        w = pr - prqr*tmp2
        IF (w /= zero) THEN

            lam1 = (tmp1*prqr - rrp)/w
            lam2 = -tmp1 - lam1*tmp2

            IF (lam1*(lam1 - one) < zero .AND. &
                lam2*(lam1 + lam2 - one) < zero) GO TO 200
        END IF
    END IF


    ! Minimum on the boundary

    lam1 = zero
    lam2 = zero
    IF (alfn <= alfv) lam2 = one
    IF (qr > zero) lam2 = MIN(one,MAX(zero,-rrq/qr))
    w = (lam2*qr + rrq + rrq)*lam2
    tmp1 = zero
    IF (alfv >= zero) tmp1 = one
    IF (pr > zero) tmp1 = MIN(one,MAX(zero,-rrp/pr))
    tmp2 = (tmp1*pr + rrp + rrp)*tmp1
    IF (tmp2 < w) THEN
        w = tmp2
        lam1 = tmp1
        lam2 = zero
    END IF
      
    IF (qqp*(qqp - pq) < zero) THEN
        IF (qr + rrq + rrq - qqp*qqp/pq < w) THEN
            lam1 = qqp/pq
            lam2 = one - lam1
        END IF
    END IF

200 CONTINUE
    IF (lam1 == zero .AND. lam2*(lam2 - one) < zero &
        .AND. -rrp - lam2*prqr > zero .AND. pr > zero) &
        lam1 = MIN(one - lam2, (-rrp-lam2*prqr)/pr)
      

    ! Computation of the aggregate values
      
    tmp1 = one - lam1 - lam2
    DO i=1,n
        ga(i)=lam1*gp(i)+lam2*g(i)+tmp1*ga(i)
    END DO
    
    alfv = lam2*alfn + tmp1*alfv
      
END SUBROUTINE agskip

!************************************************************************
!*
!*     * SUBROUTINE wprint *
!*
!*     Printout the warning and error messages.
!*
!************************************************************************
      
SUBROUTINE wprint(iterm,iprint,nout)
    IMPLICIT NONE

    ! Scalar Arguments
    INTEGER, INTENT(IN) :: &
        iprint, &    ! Printout specification:
                     !  -1  - No printout.
                     !   0  - Only the error messages.
                     !   1  - The final values of the objective
                     !        function.
                     !   2  - The final values of the objective
                     !        function and the most serious
                     !        warning messages.
                     !   3  - The whole final solution.
                     !   4  - At each iteration values of the
                     !        objective function.
                     !   5  - At each iteration the whole
                     !        solution
        nout, &      ! Auxilary printout specification.
        iterm        ! Cause of termination:
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
                      !  -1  - Two consecutive restarts.
                      !  -2  - Number of restarts > maximum number of restarts.
                      !  -3  - Failure in function or subgradient calculations 
                      !        (assigned by the user).
                      !  -4  - Failure in attaining the demanded accuracy.
                      !  -5  - Invalid input parameters.


    IF (iprint >= 0) THEN

        ! Initial error messages

        IF (iterm == -5) THEN
            IF (nout == 1) WRITE (6,FMT='(1X,''Error: '' &
               ''Number of variables (n) is too small, iterm='',I3)'            ) iterm
            IF (nout == 2) WRITE (6,FMT='(1X,''Error: '' &
               ''The maximum number of stored corrections (mcu) '' &
               ''is too small, iterm='',I3)'            ) iterm
            IF (nout == 3) WRITE (6,FMT='(1X,''Error: '' &
               ''The size of the bundle (na) is too small, iterm='' &
               ,I3)'            ) iterm
            IF (nout == 4) WRITE (6,FMT='(1X,''Error: '' &
               ''Line search parameter epsl >= 0.25, iterm='',I3)'            ) iterm
            RETURN
        END IF

        
        ! Warning messages

        IF (iprint >= 2) THEN
            IF (iterm == 0) THEN
                IF (nout == -1) WRITE (6,FMT='(1X,''Warning: '' &
                  ''mc > mcu. Assigned mc = mcu.'')'                )
                IF (nout == -2) WRITE (6,FMT='(1X,''Warning: '' &
                  ''A line search parameter epsr >= 0.5.'')'                )
                IF (nout == -3) WRITE (6,FMT='(1X,''Warning: '' &
                  ''A nondescent search direction occured. Restart.'')'                )
                IF (nout == -4) WRITE (6,FMT='(1X,''Warning: '' &
                  ''Does not converge.'')'                )
                IF (nout == -5) WRITE (6,FMT='(1X,''Warning: '' &
                  ''tmax < tmin. Restart.'')'                )
                RETURN
            END IF
         

            ! Printout the final results
            
            IF (iterm == 6) WRITE (6,FMT='(1X,''Abnormal exit: Time is up.'')')
            IF (iterm == 7) WRITE (6,FMT='(1X,''Abnormal exit: f < tolb.'')')
            IF (iterm == 2) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''Too many steps without significant progress.'')'            )
            IF (iterm == 3) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''The value of the function does not change.'')'            )
            IF (iterm == 5) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''Number of iterations > '',I5)'            ) nout
            IF (iterm == 4) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''Number of function evaluations > '',I5)'            ) nout
            IF (iterm == -1) THEN
                IF (nout == -1) THEN
                    WRITE (6,FMT='(1X,''Abnormal exit: Two consecutive restarts.'')')
                ELSE
                    WRITE (6,FMT='(1X,''Abnormal exit: tmax < tmin in two subsequent iterations.'')')
                END IF
            END IF
            IF (iterm == -2) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''Number of restarts > '',I5''.'')'            ) nout
            IF (iterm == -3) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''Failure in function or subgradient calculations.'')'            )
            IF (iterm == -4) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               ''Failure in attaining the demanded accuracy.'')'            )
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
      
SUBROUTINE rprint(n,nit,nfe,nge,x,f,wk,qk,iterm,iprint)
    IMPLICIT NONE

    ! Scalar Arguments
    INTEGER, INTENT(IN) :: & 
        n, &         ! Number of variables
        nit, &       ! Number of used iterations.
        nfe, &       ! Number of used function evaluations.
        nge, &       ! Number of used subgradient evaluations.
        iprint, &    ! Printout specification:
                     !  -1  - No printout.
                     !   0  - Only the error messages.
                     !   1  - The final values of the objective
                     !        function.
                     !   2  - The final values of the objective
                     !        function and the most serious
                     !        warning messages.
                     !   3  - The whole final solution.
                     !   4  - At each iteration values of the
                     !        objective function.
                     !   5  - At each iteration the whole
                     !        solution
        iterm        ! Cause of termination:
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
                      !  -1  - Two consecutive restarts.
                      !  -2  - Number of restarts > maximum number of restarts.
                      !  -3  - Failure in function or subgradient calculations 
                      !        (assigned by the user).
                      !  -4  - Failure in attaining the demanded accuracy.
                      !  -5  - Invalid input parameters.


    REAL(KIND=prec), INTENT(IN) :: &
        f, &         ! Value of the objective function.
        wk, &        ! Value of the first stopping criterion.
        qk           ! Value of the second stopping criterion.

    ! Array Arguments
    REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
        x            ! Vector of variables
         
    ! Local Scalars
    INTEGER :: i

    ! Intermediate results
    
    IF (iterm == 0) THEN
        IF (iprint > 3) WRITE (6,FMT='(1X,''nit='',I5,2X, &
            ''nfe='',I5,2X,''nge='',I5,2X,''f='',D15.8,2X,''wk='',D11.4,2X, &
            ''qk='',D11.4,2X)'        ) nit,nfe,nge,f,wk,qk
        IF (iprint == 5) WRITE (6,FMT='(1X,''x='', &
            5D15.7:/(4X,5D15.7))'        )(x(i),i=1,n)
        RETURN
    END IF
         

    ! Final results

    IF (iprint > 0) WRITE (6,FMT='(1X,''nit='',I5,2X, &
         ''nfe='',I5,2X,''nge='',I5,2X,''f='',D15.8,2X,''wk='',D11.4,2X, &
         ''qk='',D11.4,2X,''iterm='',I3)'    ) nit,nfe,nge,f,wk,qk,iterm
    IF (IPRINT .EQ. 3 .OR. IPRINT .EQ. 5) &
        WRITE (6,FMT='(1X,''x='',5D15.7:/(4X,5D15.7))')(x(i),i=1,n)
      
END SUBROUTINE rprint
      
END MODULE lmbm_mod
