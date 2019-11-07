!*************************************************************************
!*                                                                       *
!*     LMBM-Clust - Nonsmooth Optimization based Incremental             *
!*                  Clustering Software using LMBM (version 2.1)         *
!*                                                                       *
!*     by Napsu Karmitsa 2016 (last modified 10.11.2017).                *
!*     The code is based on clustering algorithms by Adil Bagirov 2015.  *
!*                                                                       *
!*     The software is free for academic teaching and research           *
!*     purposes but I ask you to refer the appropriate references        *
!*     given below, if you use it.                                       *
!*                                                                       *
!*************************************************************************
!*
!*
!*     Codes included:
!*
!*     lmbmclust.f95         - Mainprogram for clustering software
!*                             (this file).
!*     parameters.f95        - Parameters. Inludes modules:
!*                               - r_precision - Precision for reals,
!*                               - param - Parameters,
!*                               - exe_time - Execution time.
!*     initlmbmclust.f95     - initialization of clustering parameters and LMBM.
!*                             Includes modules:
!*                               - initclust - Initialization of parameters for clustering.
!*                               - initlmbm  - Initialization of LMBM.
!*     lmbmclustmod.f95      - Subroutines for clustering software.
!*     functionmod.f95       - Computation of clustering function and (sub)gradients values.
!*     lmbm.f95              - LMBM - limited memory bundle method.
!*     objfun.f95            - Computation of the function and subgradients values.
!*     subpro.f95            - subprograms for LMBM.
!*
!*     Makefile              - makefile.
!*
!*
!*     To use the software modify initlmbmclust.f95 as needed.
!*
!*
!*     References:
!*
!*     for LMBM-Clust:
!*
!*       N. Karmitsa, A. Bagirov and S. Taheri, "Clustering in Large Data Sets with the Limited 
!*       Memory Bundle Method", submitted, 2016.
!*
!*       N. Karmitsa, A. Bagirov and S. Taheri, "MSSC Clustering of Large Data using the Limited 
!*       Memory Bundle Method" TUCS Technical Report, No. 1164, Turku Centre for Computer Science, 
!*       Turku, 2016.
!*
!*     for LMBM:
!*       N. Haarala, K. Miettinen, M.M. Mäkelä, "Globally Convergent Limited Memory Bundle Method  
!*       for Large-Scale Nonsmooth Optimization", Mathematical Programming, Vol. 109, No. 1,
!*       pp. 181-205, 2007. DOI 10.1007/s10107-006-0728-2.
!*
!*       M. Haarala, K. Miettinen, M.M. Mäkelä, "New Limited Memory Bundle Method for Large-Scale 
!*       Nonsmooth Optimization", Optimization Methods and Software, Vol. 19, No. 6, pp. 673-692, 2004. 
!*       DOI 10.1080/10556780410001689225.
!*
!*       N. Karmitsa, "Numerical Methods for Large-Scale Nonsmooth Optimization" in Big Data
!*       Optimization. A. Emrouznejad (eds.), Springer, 2016.
!*
!*     for NSO (and clustering):
!*       A. Bagirov, N. Karmitsa, M.M. Mäkelä, "Introduction to nonsmooth optimization: theory, 
!*       practice and software", Springer, 2014.
!*
!*
!*     Acknowledgements:
!*
!*     The work was financially supported by the Academy of Finland (Project No. 289500).
!*
!*************************************************************************
!*
!*     * PROGRAM lmbmclust *
!*
!*     Main program for nonsmooth clustering software with LMBM.
!*
!*************************************************************************

PROGRAM lmbmclust

    USE r_precision, ONLY : prec        ! Precision for reals.
    USE param, ONLY : zero, one, large  ! Parameters.

    USE initclust, ONLY : &             ! Initialization of clustering parameters.
        infile, &                       ! Dataset file.
        outfile0, &                     ! Result file with cluster centers.
        outfile1, &                     ! Result file with function values,
                                        ! validity indices, and cpu-time.
        maxsize, &                      ! Maximum number of candidate points of data set.
        a, &                            ! Data matrix a(nrecord,nft), from input file.
        xbest, &
        dcent, &                        ! Distance (affinity) matrix for cluster centers
        tlimit, &                       ! Time limit, from user.
        tnorm, &                        ! The total number of norms computed.
        nclust, &                       ! Maximum number of clusters, from user.
        nft, &                          ! Number of features in data, from user.
        nrecord, &                      ! Number of instances in data, from user.
        mf, &                           ! Number of used features:
                                        !   mf = nft, when no classes,
        nc, &                           ! Current number of clusters, loops from 1 to nclust.
        ns, &                           ! Switch for auxiliary and real clustering problem.
        m, &                            ! Number of variables in optimization:
                                        !   m = mf    if ns = 1,
                                        !   m = mf*nc if ns = 2.
        ng2, &
        init_clustpar, &                ! Furher initialization of parameters.
        def_clustpar                    ! Default values of clustering parameters.
    USE clusteringmod                   ! Subprograms for clustering.
    USE lmbm_mod, ONLY : optim2         ! LMBM method for optimization.
    USE exe_time, ONLY : getime         ! Execution time.

    IMPLICIT NONE


    REAL(KIND=prec), DIMENSION(:), allocatable :: &
        x, &
        z, &
        x6
    REAL(KIND=prec), DIMENSION(maxsize) :: &
        x2, &
        x5
    REAL(KIND=prec), DIMENSION(:), allocatable :: &
        amed
    REAL(KIND=prec), DIMENSION(:), allocatable :: &
        fval, &
        fval1, &
        fval2
    REAL(KIND=prec) :: &
        barf, &
        fbarmin1, &
        fbarmin2, &
        db,db2,dn,dn3,sep, & !validity indices
        f, &
        f31, &
        fbarmin, &
        fbest, &
        fclust1, &
        fcurrent, &
        gamma03, &
        gamma04, &
        toler
    REAL :: &
        time1, &
        time3, &
        time4, &
        time5, &
        timef
    INTEGER :: &
        i,j,k,k1,j1, &
        nstart, &
        nstart1, &
        nstart2

    allocate(x(nclust*nft),z(nclust*nft),x6(nclust*nft),amed(nft),fval(nrecord),fval1(nrecord),fval2(nrecord))

    CALL init_clustpar()
    CALL def_clustpar()

    OPEN(40,file=outfile0)
    OPEN(42,file=outfile1)
    OPEN(78,file=infile,status='old',form='formatted')

    !===========================================================
    WRITE(42, *) 'Optimization with LMBM.'
    WRITE(40, *) 'Optimization with LMBM.'
    WRITE(42,*)
    WRITE(42, *) 'nclust |', ' f |',' D-B |',' Dunn |','   norms   |', ' CPU-time '
    !WRITE(42, *) 'nclust |', ' f |',' D-B |',' D-B2 |',' DI_min |', &
    !    ' DI_max |', '   norms   |', ' CPU-time '
    !==================================================================
    mf = nft
    tnorm = zero
    CALL getime(time1)

    DO i=1,nrecord
        READ(78,*) (a(k,i),k=1,nft)
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
                CALL optim2(z,x2((j-1)*mf+1:(j-1)*mf+mf),barf) ! Call for LMBM

                fval(j)=barf
                IF (fbarmin > barf) THEN
                    fbarmin=barf
                END IF
            END DO

            nstart2=ng2*nstart/10
            nstart2=MAX(1,nstart2)
            nstart1=0
            gamma04=one

            loop871: DO
                gamma03=gamma04
                gamma04=gamma04+1.0d-02
                fbarmin1=gamma03*fbarmin
                fbarmin2=gamma04*fbarmin
                DO j=1,nstart
                    IF ((fval(j) >= fbarmin1).AND.(fval(j) < fbarmin2)) THEN
                        nstart1=nstart1+1
                        DO k=1,mf
                            x5(k+(nstart1-1)*mf)=x2(k+(j-1)*mf)
                        END DO
                        fval2(nstart1)=fval(j)
                    END IF
                    IF(nstart1 >= nstart2) EXIT loop871
                END DO
            END DO loop871


872         nstart=nstart1
            DO i=1,nstart
                fval(i)=fval2(i)
            END DO
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
                ns=2
                CALL optim2(x,x6,fcurrent) ! Call for LMBM

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

            DO k=1,nc
                dcent(k,k)=zero
            END DO
            DO k=1,nc
                DO k1=k+1,nc
                    dcent(k,k1)=zero
                    DO j=1,mf
                        dcent(k,k1)=dcent(k,k1)+(xbest(j+(k-1)*mf)-xbest(j+(k1-1)*mf))**2
                    END DO
                    dcent(k1,k)=dcent(k,k1)
                END DO
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

        CALL check(x,db,db2,dn,dn3,sep)
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
        WRITE(42,603) nc,f,db,dn,tnorm,time4
603     FORMAT(i8,f28.4,f15.5,f15.5,f20.0,f14.4)
        !WRITE(42,603) nc,f,db,db2,dn,dn3,tnorm,time4
!603     FORMAT(i8,f28.4,f15.5,f15.5,f15.5,f15.5,f20.0,f14.4)
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
    deallocate(x,z,x6,amed,fval,fval1,fval2)

    STOP

END PROGRAM lmbmclust
