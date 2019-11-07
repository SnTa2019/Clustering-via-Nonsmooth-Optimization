!*************************************************************************
!*                                                                       *
!*     DC-Clust - DC-Clust method for nonsmooth DC Optimization based    *
!*     clustering using L2 norm.                                         *
!*                                                                       *
!*     Original code by Adil Bagirov 2015.                               *
!*                                                                       *
!*     Fortran 95 version (with minor changes)                           *
!*     by Napsu Karmitsa 2016 (last modified 26.02.2016).                *
!*                                                                       *
!*************************************************************************
!*
!*     Modules included:
!*
!*     dcclust_method      !
!*

MODULE dcclust_method      ! DC-Clust Method

    USE r_precision, ONLY : prec      ! Precision for reals.
    IMPLICIT NONE

    ! MODULE dcclust_method includes the following subroutines (S) and functions (F).

    PUBLIC :: &
        optim              ! S  Initialization an calling for DC-Clust subroutine.
    PRIVATE :: &
        func, &            ! S  Computation of the value of an auxiliary of cluster function.
        qsm, &
        optimum, &
        wolfe, &           ! S  Subroutines wolfe and equations solves quadratic
                           !    programming problem, to find descent direction.
        equations, &       ! S
        armijo, &          ! S  Line search
        dgrad              ! S  Computation of (sub)gradients.

CONTAINS

    SUBROUTINE optim(x0,x,fvalue) ! Initialization an calling for DC-Clust subroutine

        USE param, ONLY : zero, one, two
        USE initclust, ONLY : &
            maxdim, &         ! Maximum number of variables in optimization,
                              !   maxdim >= maxclust * maxnft.
            maxrec, &         ! maximum number of reports in dataset
            maxclust, &       ! maximum number of clusters
            mf, &             ! Number of used features:
                              !   mf = nft, when no classes,
                              !   mf = nft - 1, when classes.
            nc, &             ! Current number of clusters, loops from 1 to nclust
            ns                ! Switch for auxiliary and real clustering problem.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            x(maxdim), &
            x0(maxdim), &
            fvalue

        ! Local variables
        INTEGER :: &
            n, &              ! number of variables in optimization
            nbundle, &        ! size of the bundle
            j

        INTRINSIC :: MIN

        !=======================================================
        !  Input data:
        !  n        - number of variables
        !  nbundle  - size of bundle
        !=======================================================
        IF(ns == 1) n = mf
        IF(ns == 2) n = nc*mf
        nbundle=MIN(30,n+3)

        !=======================================================
        !  Calling QSM
        !=======================================================
        DO j=1,n
            x(j)=x0(j)
        END DO

        CALL qsm(n,x0,nbundle,x,fvalue)

        PRINT*,fvalue

        RETURN

    END SUBROUTINE optim

    !======================================================================
    SUBROUTINE func(x,objf) ! Computation of the value of an auxiliary of cluster function.

        USE initclust, ONLY : &
            maxdim, &      !
            ns
        USE initdcclust, ONLY : &
            nf             ! Number of function evaluations.
        USE functionmod, ONLY : &
            auxfunc, &
            clusterfunc

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            x(maxdim), &
            objf

        ! Local variables
        REAL(KIND=prec) ::  &
            f

        IF(ns == 1) CALL auxfunc(x,f)
        IF(ns == 2) CALL clusterfunc(x,f)
        objf=f
        nf=nf+1

        RETURN

    END SUBROUTINE func


    !======================================================================
    ! Without scaling
    !======================================================================
    SUBROUTINE qsm(nvar,x0,nbundle,x,f)

        USE param, ONLY : one
        USE initclust, ONLY : &
            maxdim, &         !
            m, &
            ns
        USE initdcclust, ONLY : &
            niter, &
            nf, &
            maxiter

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            x(maxdim), &
            x0(maxdim), &
            f
        INTEGER :: &
            nvar, &
            nbundle

        ! Local variables
        REAL(KIND=prec):: &
            slinit
        INTEGER :: &
            i,j

        m=nvar
        nf=0
        niter=0
        maxiter=5000
        slinit=one
        DO i=1,m
            x(i)=x0(i)
        END DO
        CALL optimum(x,nbundle,slinit)
        CALL func(x,f)

        RETURN

    END SUBROUTINE qsm

    !============================================================
    SUBROUTINE optimum(x,nbundle,slinit)

        USE param, ONLY : zero, one
        USE initclust, ONLY : &
            maxdim, &         !
            m, &
            ns
        USE initdcclust, ONLY : &
            maxdg, &          ! for optimization
            niter, &
            nf, &
            maxiter, &
            z_opt, &
            kmin, &
            dg3, &
            ij, &
            jvertex


        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            x(maxdim), &
            slinit

        INTEGER :: &
            nbundle

        ! Local variables
        REAL(KIND=prec):: &
            x1(maxdim), &
            g(maxdim), &
            v(maxdim), &
            w(maxdg,maxdim), &
            prod(maxdg,maxdg), &
            fvalues(maxiter), &
            dist1, &
            step0, &
            div, &
            eps0, &
            slmin, &
            sl, &
            sdif, &
            f1,f2,f3,f4,f5, &
            decreas, &
            dotprod, &
            r, &
            ratio1, &
            rmean, &
            rmin, &
            step, &
            toler
        INTEGER :: &
            mturn,mturn2, &
            ndg, &
            nnew, &
            i,j


        !====================================================================
        dist1=1.0E-06_prec
        step0=-2.0E-01_prec
        div=1.0E-01_prec
        eps0=1.0E-06_prec
        slmin=1.0E-10_prec
        sdif=1.0E-05_prec
        mturn=4
        !====================================================================
        sl=slinit/div
        CALL func(x,f2)
        outerloop: DO
            sl=div*sl
            IF(sl < slmin) RETURN
            DO i=1,m
                g(i)=one/dsqrt(REAL(m,prec))
            END DO
            nnew=0
            !================================================================
            innerloop: DO
                niter=niter+1
                IF(niter > maxiter) RETURN
                nnew=nnew+1
                f1=f2
                fvalues(niter)=f1
                !---------------------------------------------------------------
                IF (nnew > mturn) THEN
                    mturn2=niter-mturn+1
                    ratio1=(fvalues(mturn2)-f1)/(dabs(f1)+one)
                    IF(ratio1 < sdif) CYCLE outerloop
                END IF
                IF (nnew >= (2*mturn)) THEN
                    mturn2=niter-2*mturn+1
                    ratio1=(fvalues(mturn2)-f1)/(dabs(f1)+one)
                    IF(ratio1 < (1.0E+01_prec*sdif)) CYCLE outerloop
                END IF
                !--------------------------------------------------------------
                DO ndg=1,nbundle
                    CALL dgrad(ndg,x,sl,g,v)
                    dotprod=zero
                    DO i=1,m
                        dotprod=dotprod+v(i)*v(i)
                    END DO
                    r=dsqrt(dotprod)
                    IF(r < eps0) CYCLE outerloop
                    IF(ndg == 1) THEN
                        rmean=r
                        kmin=1
                        rmin=r
                    END IF
                    IF(ndg > 1) THEN
                        rmin=dmin1(rmin,r)
                        IF(r == rmin) kmin=ndg
                        rmean=((ndg-1)*rmean+r)/ndg
                    END IF
                    toler=dmax1(eps0,dist1*rmean)
                    DO i=1,ndg-1
                        prod(ndg,i)=zero
                        DO j=1,m
                            prod(ndg,i)=prod(ndg,i)+w(i,j)*v(j)
                        END DO
                        prod(i,ndg)=prod(ndg,i)
                    END DO
                    prod(ndg,ndg)=dotprod
                    !====================================================================
                    DO i=1,m
                        w(ndg,i)=v(i)
                    END DO
                    CALL wolfe(ndg,prod)
                    !================================
                    DO i=1,m
                        v(i)=zero
                        DO j=1,jvertex
                            v(i)=v(i)+w(ij(j),i)*z_opt(j)
                        END DO
                    END DO
                    !================================
                    r=zero
                    DO i=1,m
                        r=r+v(i)*v(i)
                    END DO
                    r=dsqrt(r)
                    IF (r < toler) CYCLE outerloop
                    !===========================================================
                    DO i=1,m
                        g(i)=-v(i)/r
                        x1(i)=x(i)+sl*g(i)
                    END DO
                    !===========================================================
                    CALL func(x1,f4)
                    f3=(f4-f1)/sl
                    decreas=step0*r
                    IF (f3 < decreas) THEN
                        CALL armijo(x,g,f1,f5,f4,sl,step,r)
                        f2=f5
                        DO i=1,m
                            x(i)=x(i)+step*g(i)
                        END DO
                        sl=1.2E+00_prec*sl
                        CYCLE innerloop
                    END IF
                END DO
                EXIT innerloop
            END DO innerloop
        !=====================================================
        END DO outerloop

        RETURN

    END SUBROUTINE optimum


    !==============================================================
    !  Subroutines Wolfe and Equations solves quadratic
    !  programming problem, to find
    !  descent direction, Step 3, Algorithm 2.
    !===============================================================

    SUBROUTINE wolfe(ndg,prod)

        USE param, ONLY : zero, one
        USE initclust, ONLY : &
            maxdim, &         !
            m
        USE initdcclust, ONLY : &
            maxdg, &          ! for optimization
            z_opt, &
            kmin, &
            ij, &
            jvertex, &
            aa


        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            prod(maxdg,maxdg)

        INTEGER :: &
            ndg

        ! Local variables
        REAL(KIND=prec):: &
            z1(maxdg), &
            z5, &
            r,r2, &
            rm, &
            teta, &
            t0,t1
        INTEGER :: &
            j2,j9, &
            jmax, &
            kmax, &
            kzero, &
            i,j


        j9=0
        jmax=500*ndg
        jvertex=1
        ij(1)=kmin
        z_opt(1)=one
        !=======================================
        !  To calculate norm of X
        !=======================================
        outerloop: DO
            r=zero
            DO i=1,jvertex
                DO j=1,jvertex
                    r=r+z_opt(i)*z_opt(j)*prod(ij(i),ij(j))
                END DO
            END DO
            IF(ndg == 1) RETURN
            !========================================
            !  To calculate <X,P_J> and J
            !========================================
            t0=1.0E+12_prec
            DO i=1,ndg
                t1=zero
                DO j=1,jvertex
                    t1=t1+z_opt(j)*prod(ij(j),i)
                END DO
                IF (t1 < t0) THEN
                    t0=t1
                    kmax=i
                END IF
            END DO
            !========================================
            !  First stopping criterion
            !========================================
            rm=prod(kmax,kmax)
            DO j=1,jvertex
                rm=dmax1(rm,prod(ij(j),ij(j)))
            END DO
            r2=r-1.0E-12_prec*rm
            IF (t0 > r2) RETURN
            !========================================
            !  Second stopping criterion
            !========================================
            DO i=1,jvertex
                IF (kmax == ij(i)) RETURN
            END DO
            !========================================
            ! Step 1(e) from Wolfe's algorithm
            !========================================
            jvertex=jvertex+1
            ij(jvertex)=kmax
            z_opt(jvertex)=zero
                    !========================================
            innerloop: DO
                DO i=1,jvertex
                    DO j=1,jvertex
                        aa(i,j)=one+prod(ij(i),ij(j))
                    END DO
                END DO
                j9=j9+1
                IF (j9 > jmax) RETURN
                CALL equations(jvertex,z1)
                DO i=1,jvertex
                    IF (z1(i) < 1.0E-10_prec) EXIT outerloop
                END DO
                DO i=1,jvertex
                    z_opt(i)=z1(i)
                END DO
                CYCLE outerloop

                teta=one
                DO i=1,jvertex
                    z5=z_opt(i)-z1(i)
                    IF (z5 > 1.0E-10_prec) teta=dmin1(teta,z_opt(i)/z5)
                END DO
                kzero=0
                DO i=1,jvertex
                    z_opt(i)=(one-teta)*z_opt(i)+teta*z1(i)
                    IF (z_opt(i) <= 1.0E-10_prec) THEN
                        z_opt(i)=zero
                        kzero=i
                    END IF
                END DO
                j2=0
                DO i=1,jvertex
                    IF(i /= kzero) THEN
                        j2=j2+1
                        ij(j2)=ij(i)
                        z_opt(j2)=z_opt(i)
                    END IF
                END DO
                jvertex=j2

            END DO innerloop
        END DO outerloop

        RETURN

    END SUBROUTINE wolfe


    SUBROUTINE equations(n,z1)

        USE param, ONLY : zero, one
        USE initclust, ONLY : &
            maxdim            !
        USE initdcclust, ONLY : &
            maxdg, &          ! for optimization
            aa

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) :: &
            z1(maxdg)

        INTEGER :: &
            n

        ! Local variables
        REAL(KIND=prec):: &
            b(maxdg,maxdg), &
            r, &
            z2
        INTEGER :: &
            i,j,k

        DO i=1,n
            DO j=1,n
                b(i,j)=aa(i,j)
            END DO
            b(i,n+1)=one
        END DO
        DO i=1,n
            r=b(i,i)
            DO j=i,n+1
                b(i,j)=b(i,j)/r
            END DO
            DO j=i+1,n
                DO k=i+1,n+1
                    b(j,k)=b(j,k)-b(i,k)*b(j,i)
                END DO
            END DO
        END DO
        z1(n)=b(n,n+1)
        DO i=1,n-1
            k=n-i
            z1(k)=b(k,n+1)
            DO j=k+1,n
                z1(k)=z1(k)-b(k,j)*z1(j)
            END DO
        END DO
        z2=zero
        DO i=1,n
            z2=z2+z1(i)
        END DO
        DO i=1,n
            z1(i)=z1(i)/z2
        END DO

        RETURN

    END SUBROUTINE equations

    !=====================================================================
    ! Subroutine dgrad calculates subgradients or discrete gradients
    !=====================================================================
    subroutine dgrad(ndg,x,sl,g,dg)

        USE initclust, ONLY : &
            maxdim, &         !
            m
        USE initdcclust, ONLY : &
            maxdg, &          ! for optimization
            dg3
        USE functionmod, ONLY : &
            fgrad1,fgrad2

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            g(maxdim), &
            x(maxdim), &
            dg(maxdim), &
            sl

        INTEGER :: &
            ndg

        ! Local variables
        REAL(KIND=prec):: &
            x1(maxdim), &
            dg1(maxdim), &
            dg2(maxdim), &
            z1(maxdg)

        INTEGER :: &
            i,j,k

        do k=1,m
            x1(k)=x(k)+sl*g(k)
        end do
        call fgrad1(x1,dg1)
        if(ndg.eq.1) then
            call fgrad2(x,dg2)
            do i=1,m
                dg3(i)=dg2(i)
            end do
        end if
        do i=1,m
            dg(i)=dg1(i)-dg3(i)
        end do
    end subroutine dgrad

    !===========================================================
    ! Line search (Armijo-type), Step 5 Algorithm 2.
    !===========================================================
    subroutine armijo(x,g,f1,f5,f4,sl,step,r)

        USE param, ONLY : zero, two
        USE initclust, ONLY : &
            maxdim, &         !
            m
        USE initdcclust, ONLY : &
            maxdg             ! for optimization

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            g(maxdim), &
            x(maxdim), &
            f1,f5,f4, &
            sl,step, &
            r

        ! Local variables
        REAL(KIND=prec):: &
            x1(maxdim), &
            f30,f50, &
            s1

        INTEGER :: &
            i,k

        step=sl
        f5=f4
        s1=sl
        k=0
        loop: DO
            k=k+1
            IF(k.gt.20) RETURN
            s1=two*s1
            do i=1,m
                x1(i)=x(i)+s1*g(i)
            end do
            call func(x1,f50)
            f30=f50-f1+5.0E-02_prec*s1*r
            IF(f30.gt.zero) RETURN
            step=s1
            f5=f50
        END DO loop
    end subroutine armijo


END MODULE dcclust_method
