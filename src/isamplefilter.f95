! Importance sampling filtering of non-gaussian model

subroutine isamplefilter(yt, ymiss, timevar, zt, tt, rtv, qt, a1, p1,p1inf, u, dist, &
p, n, m, r, theta, maxiter,rankp,convtol, nnd,nsim,epsplus,etaplus,&
aplus1,c,tol,info,antithetics,w,sim,simwhat,simdim)

    implicit none

    integer, intent(in) ::  p,m, r, n,nnd,antithetics,nsim&
    ,simwhat,simdim,rankp
    integer, intent(in), dimension(p) :: dist
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    integer, intent(inout) ::info, maxiter
    integer ::  t, j,i,k,maxiter2,maxitermax,info2
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
    double precision, dimension(n) :: tmp
    double precision, dimension(n,3 * nsim * antithetics + nsim) :: w
    double precision, external :: ddot

    integer, dimension(n,p) :: ymiss2
    double precision, dimension(p,n,nsim) :: epsplus2
    double precision, dimension(r,n,nsim) :: etaplus2
    double precision, dimension(m,nsim) :: aplus12
    double precision, dimension(simdim,n,3 * nsim * antithetics + nsim) :: sim2
    double precision :: diff
    double precision :: lik

    external approx, simgaussian

    ht=0.0d0
    ytilde=0.0d0
    w=1.0d0

    epsplus2 = epsplus
    etaplus2 = etaplus
    aplus12  = aplus1
    ymiss2 = ymiss
    ymiss2(1,:) = 1

    call simgaussian(ymiss2(1,:),timevar, ytilde(1,:), zt(:,:,1), &
    ht(:,:,1), tt(:,:,1), rtv(:,:,1), &
    qt(:,:,1), a1, p1, p1inf, nnd,nsim, epsplus2(:,1,:), etaplus2(:,1,:), aplus12(:,:), &
    p, 1, m, r, info,rankp,tol,sim(:,1,:),c,simwhat,simdim,antithetics)


    if(info /= 0) then
        return
    end if
    maxitermax = 0

    do i = 1, (n-1) ! increase time

        ht = 0.0d0
        ytilde = 0.0d0

        maxiter2 = maxiter
        info2 = 0
        ! approximate
        call approx(yt(1:i,:), ymiss(1:i,:), timevar, zt(:,:,1:((i-1)*timevar(1)+1)), &
        tt(:,:,1:((i-1)*timevar(3)+1)), rtv(:,:,1:((i-1)*timevar(4)+1)), ht(:,:,1:i),&
        qt(:,:,1:((i-1)*timevar(5)+1)), a1, p1,p1inf, p,i,m,r,&
        theta(1:i,:), u(1:i,:), ytilde(1:i,:), dist,maxiter2,tol,rankp,convtol,diff,lik,info2)

        if(info2 .ne. 0 .and. info2 .ne. 3) then !check for errors in approximating algorithm
            info = info2
            return
        end if

        if(maxiter2.GT.maxitermax) then
            maxitermax = maxiter2
        end if
        epsplus2 = epsplus
        etaplus2 = etaplus
        aplus12  = aplus1
        ymiss2 = ymiss
        ymiss2(i+1,:) = 1
        sim2=0.0d0
        info2 = 0
        ! simulate signals
        call simgaussian(ymiss2(1:(i+1),:),timevar, ytilde(1:(i+1),:), zt(:,:,1:(i*timevar(1)+1)), &
        ht(:,:,1:(i+1)), tt(:,:,1:(i*timevar(3)+1)), rtv(:,:,1:(i*timevar(4)+1)), &
        qt(:,:,1:(i*timevar(5)+1)), a1, p1, p1inf, nnd,nsim, epsplus2(:,1:(i+1),:), &
        etaplus2(:,1:(i+1),:), aplus12(:,:),p, i+1, m, r, info2,rankp,tol,&
        sim2(:,1:(i+1),:),c,simwhat,simdim,antithetics)

        if(info2 /= 0) then
            info = info2
            return
        end if

        ! compute importance weights

        if(simwhat.EQ.5) then
            do j=1,p
                select case(dist(j))
                    case(2)    !poisson
                        tmp(1:i) = exp(theta(1:i,j))
                        do t=1,i
                            if(ymiss2(t,j) .EQ. 0) then
                                w(i+1,:) = w(i+1,:)*exp(yt(t,j)*(sim2(j,t,:)-theta(t,j))-&
                                u(t,j)*(exp(sim2(j,t,:))-tmp(t)))/&
                                exp(-0.5d0/ht(j,j,t)*( (ytilde(t,j)-sim2(j,t,:))**2 - (ytilde(t,j)-theta(t,j))**2))
                            end if
                        end do
                    case(3) !binomial
                        tmp(1:i) = log(1.0d0+exp(theta(1:i,j)))
                        do t=1,i
                            if(ymiss2(t,j) .EQ. 0) then

                                w(i+1,:) = w(i+1,:)*exp( yt(t,j)*(sim2(j,t,:)-theta(t,j))-&
                                u(t,j)*(log(1.0d0+exp(sim2(j,t,:)))-tmp(t)))/&
                                exp(-0.5d0/ht(j,j,t)*( (ytilde(t,j)-sim2(j,t,:))**2 -(ytilde(t,j)-theta(t,j))**2))

                            end if
                        end do
                    case(4) ! gamma
                        tmp(1:i) = exp(-theta(1:i,j))
                        do t=1,i
                            if(ymiss2(t,j) .EQ. 0) then
                                w(i+1,:) = w(i+1,:)*exp( u(t,j)*(yt(t,j)*(tmp(t)-exp(-sim2(j,t,:)))&
                                +theta(t,j)-sim2(j,t,:)))/&
                                exp(-0.5d0/ht(j,j,t)*( (ytilde(t,j)-sim2(j,t,:))**2 - (ytilde(t,j)-theta(t,j))**2))
                            end if
                        end do
                    case(5) !negbin
                        tmp(1:i) = exp(theta(1:i,j))
                        do t=1,i
                            if(ymiss2(t,j) .EQ. 0) then
                                w(i+1,:) = w(i+1,:)*exp(yt(t,j)*(sim2(j,t,:)-theta(t,j)) +&
                                (yt(t,j)+u(t,j))*log((u(t,j)+tmp(t))/(u(t,j)+exp(sim2(j,t,:)))))/&
                                exp(-0.5d0/ht(j,j,t)*( (ytilde(t,j)-sim2(j,t,:))**2 - (ytilde(t,j)-theta(t,j))**2))
                            end if
                        end do
                end select
            end do

        else
            do j=1,p
                select case(dist(j))
                    case(2)    !poisson
                        tmp(1:i) = exp(theta(1:i,j))
                        do t=1,i
                            if(ymiss2(t,j) .EQ. 0) then
                                do k=1,3 * nsim * antithetics + nsim
                                    tsim(j,k) = ddot(m,zt(j,:,(t-1)*timevar(1)+1),1,sim2(:,t,k),1)
                                end do
                                w(i+1,:) = w(i+1,:)*exp(yt(t,j)*(tsim(j,:)-theta(t,j))-&
                                u(t,j)*(exp(tsim(j,:))-tmp(t)))/&
                                exp(-0.5d0/ht(j,j,t)*( (ytilde(t,j)-tsim(j,:))**2 - (ytilde(t,j)-theta(t,j))**2))

                            end if
                        end do
                    case(3) !binomial
                        tmp(1:i) = log(1.0d0+exp(theta(1:i,j)))
                        do t=1,i
                            if(ymiss2(t,j) .EQ. 0) then
                                do k=1,3 * nsim * antithetics + nsim
                                    tsim(j,k) = ddot(m,zt(j,:,(t-1)*timevar(1)+1),1,sim2(:,t,k),1)
                                end do
                                w(i+1,:) = w(i+1,:)*exp( yt(t,j)*(tsim(j,:)-theta(t,j))-&
                                u(t,j)*(log(1.0d0+exp(tsim(j,:)))-tmp(t)))/&
                                exp(-0.5d0/ht(j,j,t)*( (ytilde(t,j)-tsim(j,:))**2 - (ytilde(t,j)-theta(t,j))**2))

                            end if
                        end do
                    case(4) ! gamma
                        tmp(1:i) = exp(-theta(1:i,j))
                        do t=1,i
                            if(ymiss2(t,j) .EQ. 0) then
                                do k=1,3 * nsim * antithetics + nsim
                                    tsim(j,k) = ddot(m,zt(j,:,(t-1)*timevar(1)+1),1,sim2(:,t,k),1)
                                end do
                                w(i+1,:) = w(i+1,:)*exp( u(t,j)*(yt(t,j)*(tmp(t)-exp(-tsim(j,:)))+theta(t,j)-tsim(j,:)))/&
                                exp(-0.5d0/ht(j,j,t)*( (ytilde(t,j)-tsim(j,:))**2 - (ytilde(t,j)-theta(t,j))**2))
                            end if
                        end do
                    case(5) !negbin
                        tmp(1:i) = exp(theta(1:i,j))
                        do t=1,i
                            if(ymiss2(t,j) .EQ. 0) then
                                do k=1,3 * nsim * antithetics + nsim
                                    tsim(j,k) = ddot(m,zt(j,:,(t-1)*timevar(1)+1),1,sim2(:,t,k),1)
                                end do
                                w(i+1,:) = w(i+1,:)*exp(yt(t,j)*(tsim(j,:)-theta(t,j)) +&
                                (yt(t,j)+u(t,j))*log((u(t,j)+tmp(t))/(u(t,j)+exp(tsim(j,:)))))/&
                                exp(-0.5d0/ht(j,j,t)*( (ytilde(t,j)-tsim(j,:))**2 - (ytilde(t,j)-theta(t,j))**2))
                            end if
                        end do
                end select
            end do
        end if
        sim(:,i+1,:) = sim2(:,i+1,:)
    end do
    maxiter=maxitermax
end subroutine isamplefilter
