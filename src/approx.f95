! Subroutine for computation of the approximating gaussian model for non-gaussian models

subroutine approx(yt, ymiss, timevar, zt, tt, rtv, ht, qt, a1, p1,p1inf, p,n,m,r,&
theta, u, ytilde, dist,maxiter,tol,rankp,convtol,diff,lik,stepmax,info)

    implicit none

    integer, intent(in) ::  p,m, r, n,rankp
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    integer, intent(in), dimension(p) :: dist
    integer, intent(inout) :: maxiter,info
    integer ::  i, k,tvrqr,kk
    double precision, intent(in) :: tol,convtol,stepmax
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
    double precision, external :: ddot,dnrm2
    integer, external :: finitex
    double precision, dimension(n,p) :: thetaold,thetanew
    double precision np

    np = dble(n*p)

    !compute rqr
    tvrqr = max(timevar(4),timevar(5))
    do i=1, (n-1)*tvrqr+1
        call dgemm('n','n',m,r,r,1.0d0,rtv(:,:,(i-1)*timevar(4)+1),m,&
        qt(:,:,(i-1)*timevar(5)+1),r,0.0d0,mr,m)
        call dgemm('n','t',m,m,r,1.0d0,mr,m,rtv(:,:,(i-1)*timevar(4)+1),m,0.0d0,rqr(:,:,i),m)
    end do

    k=0
    thetaold=theta
    do while(k < maxiter)

        k=k+1
        !compute new guess thetanew
        call approxloop(yt, ymiss, timevar, zt, tt, rtv, ht, qt, rqr, tvrqr, a1, p1,p1inf, p,n,m,r, &
        theta, thetanew, u, ytilde, dist,tol,rankp,lik)

        if(finitex(sum(thetanew))==0 .OR. finitex(maxval(exp(thetanew)))==0) then !non-finite value in linear predictor or in muhat
            if(k>1) then
                kk = 0
                do while(finitex(sum(thetanew))==0 .OR. finitex(maxval(exp(thetanew)))==0)
                    if(kk>maxiter) then !did not find valid likelihood
                        info = 1
                        return
                    end if
                    kk = kk + 1

                    theta = 0.5d0*(thetaold+theta)
                    call approxloop(yt, ymiss, timevar, zt, tt, rtv, ht, qt, rqr, tvrqr, a1, p1,p1inf, p,n,m,r, &
                    theta, thetanew, u, ytilde, dist,tol,rankp,lik)

                end do
            else !cannot correct step size as we have just began
                info = 1
                return
            end if
        end if



        ! relative euclidian distance between the new and old theta
        diff = dnrm2(n*p, thetanew-thetaold, 1)/np


        ! We might be overshooting, let's backtrack
        if(diff > stepmax .AND. k > 1) then
            kk = 0
          do while(diff > stepmax .AND. kk < maxiter)
                kk = kk + 1
                ! previous theta produced too 'big' thetanew
                ! new guess by halving the last try
                theta = 0.5d0*(thetaold+theta)

                call approxloop(yt, ymiss, timevar, zt, tt, rtv, ht, qt, rqr, tvrqr, a1, p1,p1inf, p,n,m,r, &
                theta, thetanew, u, ytilde, dist,tol,rankp,lik)

                diff = dnrm2(n*p, thetanew-thetaold, 1)/np
            end do
        end if

        if(abs(diff) < convtol) then !convergence
            theta=thetanew
            info=0
            exit
        else
            thetaold=theta
            theta=thetanew
        end if
    end do
    if(maxiter==k) then
        info=2
    end if
    maxiter=k
end subroutine approx
