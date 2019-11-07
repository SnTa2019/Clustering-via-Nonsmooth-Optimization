!*************************************************************************
!*                                                                       *
!*     Initialization of parameters for DC optimization based            *
!*     clustering software (last modified 24.02.2016)                    *
!*                                                                       *
!*************************************************************************
!*
!*     Modules included:
!*
!*     initclust        ! Initialization of parameters for clustering.
!*     initdcclust      ! Initialization of DC-Clust -solver.
!*     initdcdb         ! Initialization of DCDB -solver.
!*

MODULE initclust  ! Initialization of parameters for clustering codes.

    USE r_precision, ONLY : prec  ! Precision for reals.
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: &
        maxdim   = 50000, &       ! Maximum number of variables in optimization,
                                 !   maxdim >= maxclust * maxnft.
        maxsize  = 1000000, &    ! Maximum number of candidate points of data set,
                                 !   maxsize >> maxdim.
        maxrec   = 500000, &     ! Maximum number of recorts in dataset.
        maxnft   = 250, &        ! Maximum number of features in dataset.
        maxclust = 100, &        ! Maximum number of clusters.
        maxclass = 100           ! Maximum number of classes (not used).


    ! Names of input and output files:
    ! Give names here or uncommend subroutine read_data from dcclustering.f95
    CHARACTER(LEN=80), SAVE :: &
        infile = 'iris.txt', &                        ! Name of dataset file.
        outfile0 = 'cluster_centers_dcdb.txt', &      ! Result file with cluster centers.
        outfile1 = 'cluster_indices_dcdb.txt'         ! Result file with function values,
                                                      ! Davies-Bouldin (DB) validity index,
                                                      ! and cpu-time.
        ! infile, &                ! Name of dataset file.
        ! outfile0, &              ! Result file with cluster centers.
        ! outfile1                 ! Result file with function values, Davies-Bouldin (DB)
                                   ! validity index, and cpu-time.

    ! Input real parameters.
    ! Give values here or uncommend subroutine read_data from dcclustering.f95
    REAL, SAVE :: &
        tlimit    = 72000.0       ! Maximum CPU-time in seconds, from user.

    REAL(KIND=prec), SAVE :: & !
        a(maxrec,maxnft), &      ! Data matrix, from input file, give the name of the file above or
                                 ! uncommend subroutine read_data from dcclustering.f95.
        gamma1 = 0.1_prec, &     ! Parameter affecting the size of the set of starting point,
                                 ! from user, 0 < gamma1 < 1. Defaults: gamma1=0.1 for small problems
                                 ! and gamma1 = 0.9_prec for large problems.
        gamma2 = 0.9_prec        ! Parameter affecting the size of the set of starting point,
                                 ! from user, 0 < gamma2 < 1. Defaults: gamma2=0.9 for small problems
                                 ! and gamma2 = 0.99 for large problems.
        ! gamma1, &                ! Parameter affecting the size of the set of starting point,
        !                          ! from user, 0 < gamma1 < 1.
        ! gamma2                   ! Parameter affecting the size of the set of starting point,
        !                          ! from user, 0 < gamma2 < 1.

    ! Other Real parameters.
    REAL(KIND=prec), SAVE :: & !
        tnorm, &                 ! The total number of norms computed.
        dminim(maxrec), &        !
        plabel(maxclass)         ! Index for class output (not used).

    ! Input integer parameters.
    ! Give values here or uncommend subroutine read_data from dcclustering.f95
    INTEGER, SAVE :: & !
        optmet = 2, &            ! Optimization method:
                                 !   1 = DC-Clust,
                                 !   2 = DCDB.
        nclust = 40, &           ! Maximum number of clusters, from user.
        nft  = 4, &              ! Number of features in data.
        nrecord = 150, &         ! Number of recorts in data.
        nclass = 0, &            ! Number of classes, from user (not used).
        npurity = 2              ! Switch for classes, from user, (not used):
                                 !   npurity = 1 if class outputs,
                                 !   npurity = 2 if no class outputs.
        ! optmet, &               ! Optimization method:
        !                         !   1 = DC-Clust,
        !                         !   2 = DCDB.
        ! nclust, &               ! Maximum number of clusters, from user.
        ! nft, &                  ! Number of features in data, from user.
        ! nrecord, &              ! Number of reports in data, from user.
        ! nclass, &               ! Number of classes, from user (not used).
        ! npurity                 ! Switch for classes, from user, (not used):
        !                         !   npurity = 1 if class outputs,
        !                         !   npurity = 2 if no class outputs.

    ! Other integer parameters.
    INTEGER, SAVE :: & !
        mf, &                    ! Number of used features:
                                 !   mf = nft, when no classes,
                                 !   mf = nft - 1, when classes.
        nc, &                    ! Current number of clusters, loops from 1 to nclust
        ns, &                    ! Switch for auxiliary and real clustering problem.
        m, &                     ! Number of variables in optimization:
                                 !   m = mf    if ns = 1,
                                 !   m = mf*nc if ns = 2.
        nob(maxclass), &         ! Classes, (not used)
        nk(maxclust,maxrec), &   !
        nel(maxclust), &         !
        ncand, &                 ! Number of canditate points.
        lcand(maxrec)            !


CONTAINS

    SUBROUTINE init_clustpar()   ! User supplied subroutine for further initialization of parameters (when needed).
                                 ! May be left empty.
 
        IMPLICIT NONE

        IF (nrecord <  5000) THEN
             gamma1 = 0.1_prec
             gamma2 = 0.9_prec
             ELSE IF (nrecord < 50000) THEN
             gamma1 = 0.9_prec
             gamma2 = 0.95_prec
             ELSE IF (nrecord < 100000) THEN
             gamma1 = 0.95_prec
             gamma2 = 0.99_prec
             ELSE
             gamma1 = 0.975_prec
             gamma2 = 0.995_prec
         END IF

    END SUBROUTINE init_clustpar

END MODULE initclust


MODULE initdcclust  ! Initialization of parameters for DC-Clust method.

    USE r_precision, ONLY : prec  ! Precision for reals.
    USE initclust, ONLY : maxdim  ! Maximum dimension of x.
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: &
        maxdg=1000                ! optimization

    ! Real parameters.
    REAL(KIND=prec), SAVE :: &    ! Previous common variables
        eps, &
        dg3(maxdim), &            ! optimization
        z_opt(maxdg), &
        aa(maxdg,maxdg)

    ! Integer parameters.
    INTEGER, SAVE :: & ! Previous common variables
        nf, &
        multiple, &
        nscale, &
        niter, &
        maxiter, &
        ij(maxdg), &
        jvertex, &
        kmin


CONTAINS

    SUBROUTINE init_dcclustpar()  ! User supplied subroutine for further initialization of parameters (when needed).
                                  ! May be left empty.

        IMPLICIT NONE


    END SUBROUTINE init_dcclustpar

END MODULE initdcclust



MODULE initdcdb  ! Initialization of parameters for DCDB method.

    USE r_precision, ONLY : prec  ! Precision for reals.
    USE param, ONLY : zero, one   ! Parameters.
    USE initclust, ONLY : maxdim  ! Maximum dimension of x.

    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: &
        mcu     = 7, &            ! Maximum number of stored corrections, mcu >=1.
        mnma    = 10, &           ! Maximum number of function values used in nonmonotone line search.
        inma    = 1, &            ! Selection of line search method:
                                  !   inma = 0, Armijo line search,
                                  !   inma = 1, nonmonotone Armijo line search.
        maxnin  =    20           ! Maximum number of interpolations in Armijo line search, maxnin >= 0.
                                  ! The value maxint = 2 is recommented with inma=0 and maxint=20 with inma=1.


    ! Real parameters (if parameter value <= 0.0 the default value of the parameter will be used).
    REAL(KIND=prec), SAVE :: &
        tolb    = zero, &         ! Tolerance for the function value (default = -large).
        tolf    = 1.0E-8_prec, &  ! Tolerance for change of function values (default = 1.0E-8).
        tolf2   = -10.0_prec, &   ! Second tolerance for change of function values.
                                  !   - If tolf2 < 0 the the parameter and the corresponding termination
                                  !   criterion will be ignored (recommended with inma=1).
                                  !   - If tolf2 = 0 the default value 1.0E+4 will be used.
        tolg    = 1.0E-5_prec, &  ! Tolerance for the termination criterion (default = 1.0E-6).
        eta     = 1.0E-4_prec, &  ! Distance measure parameter, eta > 0.
                                  !   - If eta < 0  the default value 0.0001 will be used.
        epsl    = 0.24E+00, &     ! Line search parameter, 0 < epsl < 0.25 (default = 1.0E-4).
        told    = 1.0E+6, &       ! Upper bound for the diagonal elements (default = 1.0E+6).
        mintold = 1.0E-10_prec, & ! Lower bound for the diagonal elements (default = 1.0E-10).
        mintold2 = one, &         ! Lower bound for the "concave" update (default = 1.0).
        sigma   = 1.0E-10_prec    ! -sigma is the smallest value that is considered as positive
                                  ! linearization error (default = 0).

    ! Integer parameters (if value <= 0 the default value of the parameter will be used).
    INTEGER, SAVE :: &
        n       = maxdim,  &      ! Number of variables.
        mit     = 5000, &         ! Maximun number of iterations.
        mfe     = 5000, &         ! Maximun number of function evaluations.
        mtesf   =      0, &       ! Maximum number of iterations with changes of
                                  ! function values smaller than tolf (default = 10).
        iprint  =      1          ! Printout specification:
                                  !    -1  - No printout.
                                  !     0  - Only the error messages.
                                  !     1  - The final values of the objective function
                                  !          (default used if iprint < -1).
                                  !     2  - The final values of the objective function and the
                                  !          most serious warning messages.
                                  !     3  - The whole final solution.
                                  !     4  - At each iteration values of the objective function.
                                  !     5  - At each iteration the whole solution


    REAL(KIND=prec), DIMENSION(maxdim), SAVE :: x  ! Vector of variables


CONTAINS

    SUBROUTINE defaults()  ! Default values for parameters.

        USE param, ONLY: small, large, zero, one
        IMPLICIT NONE


        IF (iprint < -1) iprint  = 1               ! Printout specification.
        IF (mit   <= 0) mit      = 5000            ! Maximum number of iterations.
        IF (mfe   <= 0) mfe      = n*mit           ! Maximum number of function evaluations.
        IF (tolf  <= zero) tolf  = 1.0E-08_prec    ! Tolerance for change of function values.
        IF (tolf2 == zero) tolf2 = 1.0E+04_prec    ! Second tolerance for change of function values.
        IF (tolb  == zero) tolb  = -large + small  ! Tolerance for the function value.
        IF (tolg  <= zero) tolg  = 1.0E-06_prec    ! Tolerance for the termination criterion.
        IF (told  <= zero) told  = 1.0E+06_prec    ! Upper bound for the diagonal elements.
        IF (mintold  <= zero) told  = 1.0E-10_prec ! Lower bound for the diagonal elements.
        IF (mintold2 <= zero) told  = one          ! Lower bound for "concave" update.
        IF (eta   <  zero) eta   = 1.0E-04_prec    ! Distance measure parameter
        IF (epsl  <= zero) epsl  = 1.0E-04_prec    ! Line search parameter,
        IF (mtesf <= 0) mtesf    = 10              ! Maximum number of iterations with changes
                                                   ! of function values smaller than tolf.



    END SUBROUTINE defaults


    SUBROUTINE init_par()  ! User supplied subroutine for further initialization of parameters (when needed).
                           ! May be left empty.

        USE initclust, ONLY : ns,nc
        IMPLICIT NONE

! For big data sets use larger values
!        tolg =1.0E+00
!        IF (ns == 2) THEN
!            tolg = 1.0E+3
!            IF(nc == 5) THEN
!                tolg = 1.0E-3
!            ELSE IF(nc == 10) THEN
!                tolg = 1.0E-3
!            ELSE IF(nc == 15) THEN
!                tolg = 1.0E-3
!            ELSE IF(nc == 20) THEN
!                tolg = 1.0E-3
!            ELSE IF(nc == 25) THEN
!                tolg = 1.0E-3
!            ELSE IF(nc == 30) THEN
!                tolg = 1.0E-3
!            ELSE IF(nc == 35) THEN
!                tolg = 1.0E-3
!            ELSE IF(nc == 40) THEN
!                tolg = 1.0E-3
!            ELSE IF(nc == 45) THEN
!                tolg = 1.0E-3
!            ELSE IF(nc == 50) THEN
!                tolg = 1.0E-3
!            END IF
!        END IF

! For small data sets
        IF(ns == 1) tolg = 1.0E-5
        IF(ns == 2) tolg = 1.0E-4



    END SUBROUTINE init_par

END MODULE initdcdb
