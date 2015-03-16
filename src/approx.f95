! Subroutine for computation of the approximating gaussian model for non-gaussian models

subroutine approx(yt, ymiss, timevar, zt, tt, rtv, ht, qt, a1, p1,p1inf, p,n,m,r,&
theta, u, ytilde, dist,maxiter,tol,rankp,convtol,diff,lik,info)

    implicit none

    integer, intent(in) ::  p,m, r, n,rankp
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    integer, intent(in), dimension(p) :: dist
    integer, intent(inout) :: maxiter,info
    integer ::  i, k,tvrqr,kk,jt,dt
    double precision, intent(in) :: tol,convtol
    double precision, intent(in), dimension(n,p) :: u
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rtv
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, intent(inout), dimension(n,p) :: theta
    double precision, intent(inout), dimension(n,p) :: ytilde
    double precision, intent(inout), dimension(p,p,n) :: ht
    double precision, intent(inout) :: diff
    double precision, dimension(m,r) :: mr
    double precision, dimension(m,m,(n-1)*max(timevar(4),timevar(5))+1) :: rqr
    double precision, intent(inout) :: lik
    double precision, external :: ddot
    integer, external :: finitex
    double precision dev, devold
    double precision, dimension(n,p) :: thetanew, thetaold
    double precision, dimension(p,n) :: ft,finf
    double precision, dimension(m,p,n) :: kt,kinf

    external dgemm, pytheta, pthetafirst, approxloop, pthetarest

    !compute rqr
    tvrqr = max(timevar(4),timevar(5))
    do i=1, (n-1)*tvrqr+1
        call dgemm('n','n',m,r,r,1.0d0,rtv(:,:,(i-1)*timevar(4)+1),m,&
        qt(:,:,(i-1)*timevar(5)+1),r,0.0d0,mr,m)
        call dgemm('n','t',m,m,r,1.0d0,mr,m,rtv(:,:,(i-1)*timevar(4)+1),m,0.0d0,rqr(:,:,i),m)
    end do


    if(rankp .NE. m) then ! in case of totally diffuse initialization term p(theta) disappears
        call pthetafirst(theta, timevar, zt, tt, rqr, a1, p1, p1inf, p, m, n, devold, tol,rankp,kt,kinf,ft,finf,dt,jt)
    end if
    thetaold = theta
devold = -huge(devold)
    k=0
    do while(k < maxiter)

        k=k+1
        ! compute new guess thetanew
        call approxloop(yt, ymiss, timevar, zt, tt, rtv, ht, qt, rqr, tvrqr, a1, p1,p1inf, p,n,m,r, &
        theta, thetanew, u, ytilde, dist,tol,rankp,lik)
        ! and log(p(theta|y))
        call pytheta(thetanew, dist, u, yt, ymiss, dev, p, n)
        if(rankp .NE. m) then
            call pthetarest(thetanew, timevar, zt, tt, a1, p, m, n, dev, tol,kt,kinf,ft,finf,dt,jt)
        end if
!non-finite value in linear predictor or muhat
        if(finitex(sum(thetanew))==0 .OR. finitex(maxval(exp(thetanew)))==0 ) then 
            if(k>1) then
                kk = 0
                do while(finitex(sum(thetanew))==0 .OR. finitex(maxval(exp(thetanew)))==0)
                kk = kk + 1
                    if(kk>maxiter) then
                        info = 1
                        return
                    end if
                    !backtrack
                    theta = 0.5d0*(thetaold+theta)
                    call approxloop(yt, ymiss, timevar, zt, tt, rtv, ht, qt, rqr, tvrqr, a1, p1,p1inf, p,n,m,r, &
                          theta, thetanew, u, ytilde, dist,tol,rankp,lik)

                    call pytheta(thetanew, dist, u, yt, ymiss, dev, p, n)
                    if(rankp .NE. m) then
                        call pthetarest(thetanew, timevar, zt, tt, a1, p, m, n, dev, tol,kt,kinf,ft,finf,dt,jt)
                    end if
                end do
            else !cannot correct step size as we have just began
                info = 1
                return
            end if
        end if

        if(finitex(dev)==0) then !non-finite value of objective function
            if(k>1) then
                kk = 0
                do while(finitex(dev)==0)
                    kk = kk + 1
                    if(kk>maxiter) then !did not find valid likelihood
                        info = 2
                        return
                    end if

                    theta = 0.5d0*(thetaold+theta)
                    call approxloop(yt, ymiss, timevar, zt, tt, rtv, ht, qt, rqr, tvrqr, a1, p1,p1inf, p,n,m,r, &
                    theta, thetanew, u, ytilde, dist,tol,rankp,lik)

                    call pytheta(thetanew, dist, u, yt, ymiss, dev, p, n)
                    if(rankp .NE. m) then
                        call pthetarest(thetanew, timevar, zt, tt, a1, p, m, n, dev, tol,kt,kinf,ft,finf,dt,jt)
                    end if

                end do
            else !cannot correct step size as we have just began
                info = 2
                return
            end if
        end if


        ! decreasing deviance
        if((dev - devold)/(0.1d0 + abs(dev)) < -convtol .AND. k > 1) then
            kk = 0
            do while((dev - devold)/(0.1d0 + abs(dev)) < convtol .AND. kk < maxiter)
                kk = kk + 1
                ! previous theta produced too 'big' thetanew
                ! new guess by halving the last try
                theta = 0.5d0*(thetaold+theta)
                call approxloop(yt, ymiss, timevar, zt, tt, rtv, ht, qt, rqr, tvrqr, a1, p1,p1inf, p,n,m,r, &
                theta, thetanew, u, ytilde, dist,tol,rankp,lik)

                call pytheta(thetanew, dist, u, yt, ymiss, dev, p, n)
                if(rankp .NE. m) then
                    call pthetarest(thetanew, timevar, zt, tt, a1, p, m, n, dev, tol,kt,kinf,ft,finf,dt,jt)
                end if

            end do
        end if

        diff = abs(dev - devold)/(0.1d0 + abs(dev))
        if(diff < convtol) then !convergence
            theta=thetanew
            info=0
            exit
        else
            thetaold=theta
            theta=thetanew
            devold = dev
        end if
    end do
    if(maxiter==k) then
        info=3
    end if
    maxiter=k
end subroutine approx
