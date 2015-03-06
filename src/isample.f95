! Importance sampling of non-gaussian model

subroutine isample(yt, ymiss, timevar, zt, tt, rtv, qt, a1, p1,p1inf, u, dist, &
p, n, m, r, theta, maxiter,rankp,convtol, nnd,nsim,epsplus,etaplus,&
aplus1,c,tol,info,antithetics,w,sim,nd,ndl,simwhat,simdim)

    implicit none

    integer, intent(in) ::  p,m, r, n,nnd,antithetics,nsim, ndl,simwhat,simdim,rankp
    integer, intent(in), dimension(p) :: dist
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(ndl) :: nd
    integer, intent(in), dimension(5) :: timevar
    integer, intent(inout) :: maxiter,info
    integer ::  t, j,i,info2
    double precision, intent(in) :: convtol,tol
    double precision, intent(in), dimension(n,p) :: u
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rtv
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, intent(in),dimension(nsim) :: c
    double precision, intent(inout), dimension(p,n,nsim) :: epsplus
    double precision, intent(inout), dimension(r,n,nsim) :: etaplus
    double precision, intent(inout), dimension(m,nsim) :: aplus1
    double precision, intent(inout), dimension(n,p) :: theta
    double precision, dimension(p,p,n) :: ht
    double precision, intent(inout), dimension(simdim,n,3 * nsim * antithetics + nsim) :: sim
    double precision, dimension(p,(3 * nsim * antithetics + nsim)*(5-simwhat)) :: tsim
    double precision, dimension(n,p) :: ytilde
    double precision, dimension(n,p) :: dn
    double precision, dimension(n) :: tmp
    double precision, dimension(3 * nsim * antithetics + nsim) :: w
    double precision :: diff
    double precision, external :: ddot
    double precision :: lik

    ht=0.0d0

    ! approximate
    call approx(yt, ymiss, timevar, zt, tt, rtv, ht, qt, a1, p1,p1inf, p,n,m,r,&
    theta, u, ytilde, dist,maxiter,tol,rankp,convtol,diff,lik,info)

    if(info .ne. 0 .and. info .ne. 3) then
        return
    end if

    info2 = 0
    ! simulate signals
    call simgaussian(ymiss,timevar, ytilde, zt, ht, tt, rtv, qt, a1, p1, &
    p1inf, nnd,nsim, epsplus, etaplus, aplus1, p, n, m, r, info2,rankp,&
    tol,nd,ndl,sim,c,simwhat,simdim,antithetics)

    if(info2 /= 0) then
        info = info2
        return
    end if

    ! compute importance weights

    where(ymiss .EQ. 0) dn=(ytilde-theta)**2
    w=1.0d0

    if(simwhat==5) then
        do j=1,p
            select case(dist(j))
                case(2)    !poisson
                    tmp = exp(theta(:,j))
                    do t=1,n
                        if(ymiss(t,j) .EQ. 0) then

                            w = w*exp(yt(t,j)*(sim(j,t,:)-theta(t,j))-&
                            u(t,j)*(exp(sim(j,t,:))-tmp(t)))/&
                            exp(-0.5d0/ht(j,j,t)*( (ytilde(t,j)-sim(j,t,:))**2 - dn(t,j)))
                     
                        end if
                    end do
                case(3) !binomial
                    tmp = log(1.0d0+exp(theta(:,j)))
                    do t=1,n
                        if(ymiss(t,j) .EQ. 0) then
                       
                            w = w*exp( yt(t,j)*(sim(j,t,:)-theta(t,j))-&
                            u(t,j)*(log(1.0d0+exp(sim(j,t,:)))-tmp(t)))/&
                            exp(-0.5d0/ht(j,j,t)*( (ytilde(t,j)-sim(j,t,:))**2 -dn(t,j)))
                     
                        end if
                    end do
                case(4) ! gamma
                    tmp = exp(-theta(:,j))
                    do t=1,n
                        if(ymiss(t,j) .EQ. 0) then
                            w = w*exp( u(t,j)*(yt(t,j)*(tmp(t)-exp(-sim(j,t,:)))+theta(t,j)-sim(j,t,:)))/&
                            exp(-0.5d0/ht(j,j,t)*( (ytilde(t,j)-sim(j,t,:))**2 -dn(t,j)))
                        end if
                    end do
                case(5) !negbin
                    tmp = exp(theta(:,j))
                    do t=1,n
                        if(ymiss(t,j) .EQ. 0) then
                            w = w*exp(yt(t,j)*(sim(j,t,:)-theta(t,j)) +&
                            (yt(t,j)+u(t,j))*log((u(t,j)+tmp(t))/(u(t,j)+exp(sim(j,t,:)))))/&
                            exp(-0.5d0/ht(j,j,t)*( (ytilde(t,j)-sim(j,t,:))**2 -dn(t,j)))

                        end if
                    end do
            end select
        end do

    else
        do j=1,p
            select case(dist(j))
                case(2)    !poisson
                    tmp = exp(theta(:,j))
                    do t=1,n
                        if(ymiss(t,j) .EQ. 0) then
                            do i=1,3 * nsim * antithetics + nsim
                                tsim(j,i) = ddot(m,zt(j,:,(t-1)*timevar(1)+1),1,sim(:,t,i),1)
                            end do
                            w = w*exp(yt(t,j)*(tsim(j,:)-theta(t,j))-&
                            u(t,j)*(exp(tsim(j,:))-tmp(t)))/&
                            exp(-0.5d0/ht(j,j,t)*( (ytilde(t,j)-tsim(j,:))**2 - dn(t,j)))

                        end if
                    end do
                case(3) !binomial
                    tmp = log(1.0d0+exp(theta(:,j)))
                    do t=1,n
                        if(ymiss(t,j) .EQ. 0) then
                            do i=1,3 * nsim * antithetics + nsim
                                tsim(j,i) = ddot(m,zt(j,:,(t-1)*timevar(1)+1),1,sim(:,t,i),1)
                            end do
                            w = w*exp( yt(t,j)*(tsim(j,:)-theta(t,j))-&
                            u(t,j)*(log(1.0d0+exp(tsim(j,:)))-tmp(t)))/&
                            exp(-0.5d0/ht(j,j,t)*( (ytilde(t,j)-tsim(j,:))**2 -dn(t,j)))

                        end if
                    end do
                case(4) ! gamma
                    tmp = exp(-theta(:,j))
                    do t=1,n
                        if(ymiss(t,j) .EQ. 0) then
                            do i=1,3 * nsim * antithetics + nsim
                                tsim(j,i) = ddot(m,zt(j,:,(t-1)*timevar(1)+1),1,sim(:,t,i),1)
                            end do
                            w = w*exp( u(t,j)*(yt(t,j)*(tmp(t)-exp(-tsim(j,:)))+theta(t,j)-tsim(j,:)))/&
                            exp(-0.5d0/ht(j,j,t)*( (ytilde(t,j)-tsim(j,:))**2 -dn(t,j)))
                        end if
                    end do
                case(5) !negbin
                    tmp = exp(theta(:,j))
                    do t=1,n
                        if(ymiss(t,j) .EQ. 0) then
                            do i=1,3 * nsim * antithetics + nsim
                                tsim(j,i) = ddot(m,zt(j,:,(t-1)*timevar(1)+1),1,sim(:,t,i),1)
                            end do
                            w = w*exp(yt(t,j)*(tsim(j,:)-theta(t,j)) +&
                            (yt(t,j)+u(t,j))*log((u(t,j)+tmp(t))/(u(t,j)+exp(tsim(j,:)))))/&
                            exp(-0.5d0/ht(j,j,t)*( (ytilde(t,j)-tsim(j,:))**2 -dn(t,j)))
                        end if
                    end do
            end select
        end do


    end if



end subroutine isample
