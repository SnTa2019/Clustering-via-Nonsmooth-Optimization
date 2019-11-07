!*************************************************************************
!*                                                                       *
!*     DC Optimization based incremental clustering using L2 norm        *
!*                                                                       *
!*     Original code with DC-Clust algorithm by Adil Bagirov 2015.       *
!*                                                                       *
!*     Fortran 95 version and optimization with DCDB method              *
!*     by Napsu Karmitsa 2016 (last modified 24.02.2016).                *
!*                                                                       *
!*                                                                       *
!*     The software is free for academic teaching and research           *
!*     purposes but I ask you to refer the references given below,       *
!*     if you use it.                                                    *
!*                                                                       *
!*************************************************************************
!*
!*
!*     Codes included:
!*
!*     dcclustering.f95      - Mainprogram for clustering software
!*                             (this file).
!*     parameters.f95        - Parameters. Inludes modules:
!*                               - r_precision - Precision for reals,
!*                               - param - Parameters,
!*                               - exe_time - Execution time.
!*     initclust.f95         - initialization of clustering parameters and
!*                             optimization methods. Includes modules:
!*                               - initclust - Initialization of parameters for clustering.
!*                               - initdcclust - Initialization of DC-Clust -solver.
!*                               - initdcdb - Initialization of DCDB -solver.
!*     clusteringmod.f95     - Subroutines for clustering software.
!*     functionmod.f95       - Computation of function and (sub)gradients values for
!*                             clustering software.
!*     dcclust_method.f95    - DC-Clust method.
!*     dcdb.f95              - DCDB method.
!*     dcfun.f95             - Computation of the function and (sub)gradients
!*                             values for DCDB.
!*     subpro.f95            - subprograms for DCDB.
!*
!*     Makefile              - makefile.
!*
!*
!*     To use the software modify initclust.f95 (and dcclustering.f95) as needed.
!*
!*
!*     References:
!*
!*     "New diagonal bundle method for clustering problems in very large data sets",
!*     Napsu Karmitsa, Adil Bagirov and Sona Taheri, 2016.
!*
!*     "Nonsmooth DC programming approach to the minimum sum-of-squares 
!*     clustering problems", Adil Bagirov, Sona Taheri and Julien Ugon, 
!*     Pattern Recognition, Vol. 53, pp. 12–24, 2016.
!*
!*
!*     Acknowledgements:
!*
!*     The work was financially supported by the Academy of Finland (Project No. 289500)
!*     and Australian Research Counsil's Discovery Projects funding scheme (Project No.
!*     DP140103213).
!*
!*************************************************************************
!*
!*     * PROGRAM dcclustering *
!*
!*     Main program for clustering software.
!*
!*************************************************************************

PROGRAM dcclustering

    USE r_precision, ONLY : prec        ! Precision for reals.
    USE param, ONLY : zero, large       ! Parameters.
    USE initclust, ONLY : &             ! Initialization of clustering parameters.
        infile, &                       ! Dataset file.
        outfile0, &                     ! Result file with cluster centers.
        outfile1, &                     ! Result file with function values,
                                        ! Davies-Bouldin (DB) validity index, and cpu-time.
        maxdim, &                       ! Maximum number of variables for optimization.
        maxsize, &                      ! Maximum number of candidate points of data set.
        maxnft, &                       ! Maximum number of features in dataset.
        maxrec, &                       ! Maximum number of reports in dataset.
        a, &                            ! Data matrix a(maxrec,maxnft), from input file.
        gamma1, &                       ! parameter from user, 0 < gamma1 < 1.
        gamma2, &                       ! parameter from user, 0 < gamma2 < 1.
        tlimit, &                       ! Time limit, from user.
        tnorm, &                        ! The total number of norms computed.
        optmet, &                       ! Optimization method, from user:
                                        !   1 = DC-Clust,
                                        !   2 = DCDB.
        nclust, &                       ! Maximum number of clusters, from user.
        nft, &                          ! Number of features in data, from user.
        nrecord, &                      ! Number of reports in data, from user.
        mf, &                           ! Number of used features:
                                        !   mf = nft, when no classes,
                                        !   mf = nft - 1, when classes.
        nc, &                           ! Current number of clusters, loops from 1 to nclust.
        ns, &                           ! Switch for auxiliary and real clustering problem.
        m, &                            ! Number of variables in optimization:
                                        !   m = mf    if ns = 1,
                                        !   m = mf*nc if ns = 2.
        init_clustpar                   ! Furher initialization of parameters.
    USE clusteringmod                   ! Subprograms for clustering.
    USE dcclust_method                  ! DC-Clust method for optimization.
    USE dcdb_mod                        ! DCDB method for optimization.
    USE exe_time, ONLY : getime         ! Execution time.

    IMPLICIT NONE

    REAL(KIND=prec), DIMENSION(maxdim) :: &
        x, &
        z, &
        x3, &
        xbest, &
        x6
    REAL(KIND=prec), DIMENSION(maxsize) :: &
        x2, &
        x5
    REAL(KIND=prec), DIMENSION(maxnft) :: &
        amed
    REAL(KIND=prec), DIMENSION(maxrec) :: &
        fval, &
        fval1
    REAL(KIND=prec) :: &
        a2, &
        barf, &
        db, &
        f, &
        f31, &
        fbarmin, &
        fbest, &
        fclust1, &
        fcurrent, &
        gamma3, &
        plabel3, &
        toler
    REAL :: &
        time1, &
        time3, &
        time4, &
        time5, &
        timef
    INTEGER :: &
        i,j,k,j1, &
        n2, &
        noutcom, &
        nstart, &
        nstart1, &
        nstart2, &
        limitt

    CALL init_clustpar()

    ! You can either use this or give values directly in initclust.f95
    !    CALL read_data()

    OPEN(40,file=outfile0)
    OPEN(42,file=outfile1)
    OPEN(78,file=infile,status='old',form='formatted')



    IF (nrecord <= 1000) gamma3=2.0E+00_prec
    IF ((nrecord >= 1000) .AND. (nrecord <= 5000)) gamma3=1.25E+00_prec
    IF ((nrecord > 5000)  .AND. (nrecord < 10000)) gamma3=1.10E+00_prec
    IF ((nrecord > 10000)  .AND. (nrecord <= 50000)) gamma3=1.05E+00_prec
    IF (nrecord > 50000) gamma3=1.025E+00_prec
    !===========================================================
    IF (optmet == 1) THEN
        WRITE(42, *) 'Optimization with DC-Clust.'
        WRITE(40, *) 'Optimization with DC-Clust.'
    ELSE
        WRITE(42, *) 'Optimization with DCDB.'
        WRITE(40, *) 'Optimization with DCDB.'
    END IF
    WRITE(42,311) gamma1,gamma2,gamma3
311 FORMAT('Par-ter:',' g1=',f7.4,' g2=',f7.4,' g3=',f7.4)
    WRITE(42,*)
    WRITE(42, *) 'Number of clusters |', ' Value of cluster function |', &
        ' Davies-Bouldin validity index |', '   Number of norms   |', ' CPU-time in seconds '

    !==================================================================
    mf = nft
    tnorm = zero
    CALL getime(time1)

    DO i=1,nrecord
        READ(78,*) (a(i,k),k=1,nft)
    END DO

    !================================================================

    outerloop: DO nc=1,nclust  ! number of clusters

        PRINT 42,nc
42      FORMAT('Cluster No.:',i10)
        IF(nc > 1) THEN
            toler=1.0E-04_prec*fclust1/REAL(nc,prec)
            CALL step2(toler,nstart,x2) !  Step2 computes clusters for each data point
            fbarmin=large
            DO j=1,nstart
                DO k=1,mf
                    z(k)=x2(k+(j-1)*mf)
                END DO
                ns=1
                IF (optmet == 1) THEN  ! DC-Clust
                    CALL optim(z,x6,barf)
                ELSE  ! DCDB
                    CALL optim2(z,x6,barf)
                END IF
                fval(j)=barf
                IF (fbarmin > barf) THEN
                    fbarmin=barf
                END IF
                DO k=1,mf
                    x2(k+(j-1)*mf)=x6(k)
                END DO
            END DO

            fbarmin=gamma3*fbarmin
            nstart1=0
            DO j=1,nstart
                IF (fval(j) <= fbarmin) THEN
                    nstart1=nstart1+1
                    DO k=1,mf
                        x5(k+(nstart1-1)*mf)=x2(k+(j-1)*mf)
                    END DO
                    fval(nstart1)=fval(j)
                END IF
            END DO

            nstart=nstart1
            DO i=1,nstart
                DO k=1,mf
                    x2(k+(i-1)*mf)=x5(k+(i-1)*mf)
                END DO
            END DO

            DO k=1,mf
                x5(k)=x2(k)
            END DO
            fval1(1)=fval(1)
            nstart2=1
            innerloop: DO j=2,nstart
                DO j1=1,nstart2
                    f31=zero
                    DO k=1,mf
                        f31=f31+(x5(k+(j1-1)*mf)-x2(k+(j-1)*mf))**2
                    END DO
                    IF(f31 <= (1.0E-01_prec*toler)) THEN
                        IF(fval1(j1) >= fval(j)) THEN
                            fval1(j1)=fval(j)
                            DO k=1,mf
                                x5(k+(j1-1)*mf)=x2(k+(j-1)*mf)
                            END DO
                        END IF
                        CYCLE innerloop
                    END IF
                END DO
                nstart2=nstart2+1
                DO k=1,mf
                    x5(k+(nstart2-1)*mf)=x2(k+(j-1)*mf)
                END DO
                fval1(nstart2)=fval(j)
            END DO innerloop
            DO i=1,nstart2
                DO k=1,mf
                    x2(k+(i-1)*mf)=x5(k+(i-1)*mf)
                END DO
            END DO
            nstart=nstart2

            m=mf*nc
            fbest=large
            DO j=1,nstart
                DO i=1,mf
                    x(i+(nc-1)*mf)=x2(i+(j-1)*mf)
                END DO
                DO j1=1,m
                    x3(j1)=x(j1)
                END DO
                ns=2
                IF (optmet == 1) THEN
                    CALL optim(x3,x6,fcurrent)
                ELSE
                    CALL optim2(x3,x6,fcurrent)
                END IF
                IF (fcurrent < fbest) THEN
                    fbest=fcurrent
                    DO j1=1,m
                        xbest(j1)=x6(j1)
                    END DO
                END IF
            END DO
            f=fbest
            DO j1=1,m
                x(j1)=xbest(j1)
            END DO

            !================================================================
            WRITE(40,*)
            WRITE(40,*)
            !      WRITE(40,*) 'Final solution:'
            WRITE(40,*)

            WRITE(40,*) '____________________________________________________'
            WRITE(40,43) nc
43          FORMAT('             Total number of clusters:',i8)
            WRITE(40,*)
            WRITE(40,*) '____________________________________________________'
            WRITE(40,*)

!   Print the center of clusters
            DO j=1,nc
                WRITE(40,*)
                WRITE(40,449) j
449             FORMAT('Center of cluster No.',i4)
                WRITE(40,*)
                WRITE(40,49) (x(i+(j-1)*nft),i=1,nft)
            END DO
49          FORMAT(5f16.8)
            WRITE(40,*)
!================================================================

        ELSE  ! nc=1

            CALL step1(f,x) ! Step1 computes the centroid and the value of clustering function at centroid
            fclust1=f

        END IF

        CALL check(x,db)
        CALL getime(time3)
        WRITE(40,*)
        WRITE(40,543) f
543     FORMAT('The value of cluster function:',f28.6)
        WRITE(40,*)
        WRITE(40,142) tnorm
        time4=time3-time1
        IF(time4 > tlimit) EXIT outerloop
        WRITE(40,*)
        WRITE(40,141) time4
        WRITE(42,603) nc,f,db,tnorm,time4
603     FORMAT(i8,f28.8,f15.5,f20.0,f14.4)
    END DO outerloop
    WRITE(40,*) '____________________________________________________'
    CALL getime(time5)
    timef=time5-time1
    WRITE(40,*)
    WRITE(40,141) timef
141 FORMAT('               CPU time:',f12.3)
142 FORMAT('  The total number of norms:',f18.0)
    CLOSE(40)
    CLOSE(42)
    CLOSE(78)
    STOP

CONTAINS

    SUBROUTINE read_data() ! every variable is given by host

        IMPLICIT NONE

        PRINT *,' '
        PRINT *,'Enter the name of dataset file:'
        READ *,infile
        PRINT *,' '
        PRINT *,'Enter the name for first output file:'
        READ *,outfile0
        PRINT *,'Enter the name for second output file:'
        READ *,outfile1


        !========================================================================
        PRINT *,' '
        PRINT *,'Enter number of reports:'
        READ *,nrecord
        PRINT *,' '
        PRINT *,'Enter the number of features:'
        READ *,nft
        PRINT *,' '
        PRINT *,'Enter number of clusters:'
        READ *,nclust
        PRINT *,' '
        PRINT *,'Enter parameter gamma1 between 0 and 1:'
        READ *,gamma1
        PRINT *,' '
        PRINT *,' '
        PRINT *,'Enter parameter gamma2 between 0 and 1:'
        READ *,gamma2
        PRINT *,' '
        PRINT *,'Select optimization method. Enter:'
        PRINT *,' '
        PRINT *,'              1 - DCClust'
        PRINT *,'              2 - DCDB'
        READ *,optmet

        PRINT *,' '
        PRINT *,'Please enter maximum CPU time (in seconds):'
        READ *,limitt
        tlimit=REAL(limitt)
          !===========================================================
        RETURN
    END SUBROUTINE read_data

END PROGRAM dcclustering
