! Non-Gaussian log-likelihood computation
subroutine ngloglik(yt, ymiss, timevar, zt, tt, rtv, qt, a1, p1,p1inf, p,m,&
r, n, lik, theta, u, dist,maxiter,rankp,convtol, &
nnd,nsim,epsplus,etaplus,aplus1,c,tol,info,antit,sim,nsim2,diff,marginal)

    implicit none

    integer, intent(in) ::  p, m, r, n,nnd,antit,nsim,sim,nsim2,rankp
    integer, intent(in), dimension(p,n) :: ymiss
    integer, intent(in), dimension(p) :: dist
    integer, intent(in), dimension(5) :: timevar
    integer, intent(inout) :: maxiter,marginal,info
    integer ::  j,t,info2
    double precision, intent(in) :: convtol,tol
    double precision, intent(in), dimension(p,n) :: u
    double precision, intent(in), dimension(p,n) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rtv
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, intent(inout), dimension(m,nsim) :: aplus1
    double precision, intent(in),dimension(nsim) :: c
    double precision, dimension(p,p,n) :: ht
    double precision, dimension(p,n) :: ytilde
    double precision, intent(inout), dimension(p,n) :: theta
    double precision, intent(inout), dimension(p,n,nsim) :: epsplus
    double precision, intent(inout), dimension(r,n,nsim) :: etaplus
    double precision, intent(inout) :: lik
    double precision, dimension(p,n,nsim2) :: tsim
    double precision, dimension(n) :: tmp
    double precision, dimension(nsim2) :: w
    double precision, intent(inout) :: diff
    double precision, external :: ddot

    external approx, marginalxx, dpoisf, dnormf, dbinomf, dgammaf, dnbinomf, simgaussian
    !approximate
    call approx(yt, ymiss, timevar, zt, tt, rtv, ht, qt, a1, p1,p1inf, p, n, m, r,&
    theta, u, ytilde, dist,maxiter,tol,rankp,convtol,diff,lik, info)

    if(info .ne. 0 .and. info .ne. 3) then
        return
    end if

    if(marginal.EQ.1) then
        j = int(sum(p1inf))
        if(j.GT.0) then
            call marginalxx(p1inf,zt,tt,m,p,n,j,timevar,lik,marginal)
        end if
        if(marginal.EQ.-1) then
            info = 5
            return
        end if
    end if


    do j=1,p
        select case(dist(j))
            case(2)
                do t=1,n
                    if(ymiss(j,t).EQ.0) then
                        call dpoisf(yt(j,t), u(j,t)*exp(theta(j,t)), lik)
                        call dnormf(ytilde(j,t), theta(j,t),sqrt(ht(j,j,t)), lik)
                    end if
                end do
            case(3)
                do t=1,n
                    if(ymiss(j,t).EQ.0) then
                        call dbinomf(yt(j,t), u(j,t), exp(theta(j,t))/(1.0d0+exp(theta(j,t))), lik)
                        call dnormf(ytilde(j,t), theta(j,t),sqrt(ht(j,j,t)), lik)
                    end if
                end do
            case(4)
                do t=1,n
                    if(ymiss(j,t).EQ.0) then
                        call dgammaf(yt(j,t), u(j,t), exp(theta(j,t))/u(j,t), lik)
                        call dnormf(ytilde(j,t), theta(j,t),sqrt(ht(j,j,t)), lik)
                    end if
                end do
            case(5)
                do t=1,n
                    if(ymiss(j,t).EQ.0) then
                        call dnbinomf(yt(j,t), u(j,t), exp(theta(j,t)), lik)
                        call dnormf(ytilde(j,t), theta(j,t),sqrt(ht(j,j,t)), lik)
                    end if
                end do
        end select
    end do

    if(sim .EQ. 1) then

        w=1.0d0
        info2 = 0

        ! simulate signals
        call simgaussian(ymiss,timevar, ytilde, zt, ht, tt, rtv, qt, a1, p1, &
        p1inf, nnd,nsim, epsplus, etaplus, aplus1, p, n, m, r, info2,rankp,&
        tol,tsim,c,5,p,antit)

        if(info2.EQ.0) then
            ! Compute weights
            do j=1,p
                select case(dist(j))
                    case(2)    !poisson
                        tmp = exp(theta(j,:))
                        do t=1,n
                            if(ymiss(j,t) .EQ. 0) then
                                !  do i=1,nsim2
                                w = w*exp(yt(j,t)*(tsim(j,t,:)-theta(j,t))-&
                                u(j,t)*(exp(tsim(j,t,:))-tmp(t)))/&
                                exp(-0.5d0/ht(j,j,t)*( (ytilde(j,t)-tsim(j,t,:))**2 - (ytilde(j,t)-theta(j,t))**2))
                              !  end do
                            end if
                        end do
                    case(3) !binomial
                        tmp = log(1.0d0+exp(theta(j,:)))
                        do t=1,n
                            if(ymiss(j,t) .EQ. 0) then
                                w = w*exp( yt(j,t)*(tsim(j,t,:)-theta(j,t))-&
                                u(j,t)*(log(1.0d0+exp(tsim(j,t,:)))-tmp(t)))/&
                                exp(-0.5d0/ht(j,j,t)*( (ytilde(j,t)-tsim(j,t,:))**2 - (ytilde(j,t)-theta(j,t))**2))
                            end if
                        end do
                    case(4) ! gamma
                        tmp = exp(-theta(j,:))
                        do t=1,n
                            if(ymiss(j,t) .EQ. 0) then
                                w = w*exp( u(j,t)*(yt(j,t)*(tmp(t)-exp(-tsim(j,t,:)))+theta(j,t)-tsim(j,t,:)))/&
                                exp(-0.5d0/ht(j,j,t)*( (ytilde(j,t)-tsim(j,t,:))**2 - (ytilde(j,t)-theta(j,t))**2))
                            end if
                        end do
                    case(5)
                        tmp = exp(theta(j,:))
                        do t=1,n
                            if(ymiss(j,t) .EQ. 0) then
                                w = w*exp(yt(j,t)*(tsim(j,t,:)-theta(j,t)) +&
                                (yt(j,t)+u(j,t))*log((u(j,t)+tmp(t))/(u(j,t)+exp(tsim(j,t,:)))))/&
                                exp(-0.5d0/ht(j,j,t)*( (ytilde(j,t)-tsim(j,t,:))**2 - (ytilde(j,t)-theta(j,t))**2))
                            end if
                        end do
                end select
            end do


            lik= lik+log(sum(w)/dble(nsim2))
        else
            info = info2
            return
        end if

    end if



end subroutine ngloglik
