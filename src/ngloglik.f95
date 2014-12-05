! Non-Gaussian log-likelihood computation
subroutine ngloglik(yt, ymiss, timevar, zt, tt, rtv, qt, a1, p1,p1inf, p,m,&
r, n, lik, theta, u, dist,maxiter,rankp,convtol, &
nnd,nsim,epsplus,etaplus,aplus1,c,tol,info,antit,sim,nsim2,nd,ndl,diff,stepmax,marginal)

    implicit none

    integer, intent(in) ::  p, m, r, n,nnd,antit,nsim,sim,nsim2,ndl,rankp
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(p) :: dist
    integer, intent(in), dimension(ndl) :: nd
    integer, intent(in), dimension(5) :: timevar
    integer, intent(inout) :: maxiter,marginal,info
    integer ::  j,t,info2
    double precision, intent(in) :: convtol,tol,stepmax
    double precision, intent(in), dimension(n,p) :: u
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rtv
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, intent(inout), dimension(m,nsim) :: aplus1
    double precision, intent(in),dimension(nsim) :: c
    double precision, dimension(p,p,n) :: ht
    double precision, dimension(n,p) :: ytilde
    double precision, intent(inout), dimension(n,p) :: theta
    double precision, intent(inout), dimension(p,n,nsim) :: epsplus
    double precision, intent(inout), dimension(r,n,nsim) :: etaplus
    double precision, intent(inout) :: lik
    double precision, dimension(p,n,nsim2) :: tsim
    double precision, dimension(n,p) :: dn
    double precision, dimension(n) :: tmp
    double precision, dimension(nsim2) :: w
    double precision, intent(inout) :: diff
    double precision, external :: ddot


    !approximate
    call approx(yt, ymiss, timevar, zt, tt, rtv, ht, qt, a1, p1,p1inf, p,n,m,r,&
    theta, u, ytilde, dist,maxiter,tol,rankp,convtol,diff,lik,stepmax, info)

    if(info .ne. 0 .and. info .ne. 2) then !check for errors in approximating algorithm
        return
    end if
        
    if(marginal==1) then
        j = int(sum(p1inf))
        if(j>0) then
            call marginalxx(p1inf,zt,tt,m,p,n,j,timevar,lik,marginal)
        end if
        if(marginal==-1) then
            info = 5
            return
        end if
    end if


    do j=1,p
        select case(dist(j))
            case(2)
                do t=1,n
                    if(ymiss(t,j).EQ.0) then
                        call dpoisf(yt(t,j), u(t,j)*exp(theta(t,j)), lik)
                        call dnormf(ytilde(t,j), theta(t,j),sqrt(ht(j,j,t)), lik)
                    end if
                end do
            case(3)
                do t=1,n
                    if(ymiss(t,j).EQ.0) then
                        call dbinomf(yt(t,j), u(t,j), exp(theta(t,j))/(1.0d0+exp(theta(t,j))), lik)
                        call dnormf(ytilde(t,j), theta(t,j),sqrt(ht(j,j,t)), lik)
                    end if
                end do
            case(4)
                do t=1,n
                    if(ymiss(t,j).EQ.0) then
                        call dgammaf(yt(t,j), u(t,j), exp(theta(t,j))/u(t,j), lik)
                        call dnormf(ytilde(t,j), theta(t,j),sqrt(ht(j,j,t)), lik)
                    end if
                end do
            case(5)
                do t=1,n
                    if(ymiss(t,j).EQ.0) then
                        call dnbinomf(yt(t,j), u(t,j), exp(theta(t,j)), lik)
                        call dnormf(ytilde(t,j), theta(t,j),sqrt(ht(j,j,t)), lik)
                    end if
                end do
        end select
    end do

    if(sim .EQ. 1) then

        where(ymiss == 0) dn=(ytilde-theta)**2

        w=1.0d0
        info2 = 0

        ! simulate signals
        call simgaussian(ymiss,timevar, ytilde, zt, ht, tt, rtv, qt, a1, p1, &
        p1inf, nnd,nsim, epsplus, etaplus, aplus1, p, n, m, r, info2,rankp,&
        tol,nd,ndl,tsim,c,5,p,antit)

        if(info2==0) then
            ! Compute weights
            do j=1,p
                select case(dist(j))
                    case(2)    !poisson
                        tmp = exp(theta(:,j))
                        do t=1,n
                            if(ymiss(t,j) .EQ. 0) then
                                !  do i=1,nsim2
                                w = w*exp(yt(t,j)*(tsim(j,t,:)-theta(t,j))-&
                                u(t,j)*(exp(tsim(j,t,:))-tmp(t)))/&
                                exp(-0.5d0/ht(j,j,t)*( (ytilde(t,j)-tsim(j,t,:))**2 - dn(t,j)))
                              !  end do
                            end if
                        end do
                    case(3) !binomial
                        tmp = log(1.0d0+exp(theta(:,j)))
                        do t=1,n
                            if(ymiss(t,j) .EQ. 0) then
                                w = w*exp( yt(t,j)*(tsim(j,t,:)-theta(t,j))-&
                                u(t,j)*(log(1.0d0+exp(tsim(j,t,:)))-tmp(t)))/&
                                exp(-0.5d0/ht(j,j,t)*( (ytilde(t,j)-tsim(j,t,:))**2 -dn(t,j)))
                            end if
                        end do
                    case(4) ! gamma
                        tmp = exp(-theta(:,j))
                        do t=1,n
                            if(ymiss(t,j) .EQ. 0) then
                                w = w*exp( u(t,j)*(yt(t,j)*(tmp(t)-exp(-tsim(j,t,:)))+theta(t,j)-tsim(j,t,:)))/&
                                exp(-0.5d0/ht(j,j,t)*( (ytilde(t,j)-tsim(j,t,:))**2 -dn(t,j)))
                            end if
                        end do
                    case(5)
                        tmp = exp(theta(:,j))
                        do t=1,n
                            if(ymiss(t,j) .EQ. 0) then
                                w = w*exp(yt(t,j)*(tsim(j,t,:)-theta(t,j)) +&
                                (yt(t,j)+u(t,j))*log((u(t,j)+tmp(t))/(u(t,j)+exp(tsim(j,t,:)))))/&
                                exp(-0.5d0/ht(j,j,t)*( (ytilde(t,j)-tsim(j,t,:))**2 -dn(t,j)))
                            end if
                        end do
                end select
            end do


            lik= lik+log(sum(w)/nsim2)
        else
            info = info2
            return
        end if

    end if



end subroutine ngloglik
