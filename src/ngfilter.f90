! Filtering of non-gaussian model

subroutine ngfilter(yt, ymiss, timevar, zt, tt, rtv, qt, a1, p1,p1inf, u, theta,&
dist, p,n, m, r, rankp, nnd,nsim,epsplus,etaplus,aplus1,c,tol,info,maxiter,&
convtol,alphahat,alphavar,thetahat,thetavar,yhat,yvar,smootha,smooths,smoothy,&
expected, htol)

    implicit none

    integer, intent(in) ::  p,m, r, n,nnd,nsim,rankp,smootha,smooths,smoothy&
    , expected
    integer, intent(in), dimension(p) :: dist
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    integer, intent(inout) :: info,maxiter
    integer ::  t, j
    double precision, intent(in) :: tol,convtol
    double precision, intent(inout) :: htol
    double precision, intent(inout), dimension(n,p) :: theta
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
    double precision, intent(inout), dimension((m-1)*smootha+1,(n-1)*smootha+1) :: alphahat
    double precision, intent(inout), dimension((m-1)*smootha+1,(m-1)*smootha+1,(n-1)*smootha+1) :: alphavar
    double precision, intent(inout), dimension((p-1)*smoothy+1,(n-1)*smoothy+1) :: yhat
    double precision, intent(inout), dimension((p-1)*smoothy+1,(p-1)*smoothy+1,(n-1)*smoothy+1) :: yvar
    double precision, intent(inout), dimension((p-1)*smooths+1,(n-1)*smooths+1) :: thetahat
    double precision, intent(inout), dimension((p-1)*smooths+1,(p-1)*smooths+1,(n-1)*smooths+1) :: thetavar
    double precision, dimension(p,m) :: pm
    double precision, dimension(:,:,:), allocatable :: sim
    double precision, dimension(:,:,:), allocatable :: osim
    double precision, dimension(:,:), allocatable :: w

    double precision, external :: ddot
    external isamplefilter, covmeanwprotect, dgemv, dsymm, dgemm, covmeanw
    
    allocate(w(n,4*nsim))
    allocate(sim(smootha*m+(1-smootha)*p,n,4*nsim))

    if(smootha.EQ.1) then

        call isamplefilter(yt, ymiss, timevar, zt, tt, rtv, qt, a1, p1,p1inf, u, dist, &
        p, n, m, r, theta, maxiter,rankp,convtol, nnd,nsim,epsplus,etaplus,&
        aplus1,c,tol,info,1,w,sim,4,m, expected, htol)

        if(info .ne. 0 .and. info .ne. 3) then
            return
        end if

        do t=1,n
            w(t,:) = w(t,:)/sum(w(t,:))
            call covmeanwprotect(sim(:,t,:),w(t,:),m,1,4*nsim,alphahat(:,t),alphavar(:,:,t))
        end do

        if(smooths.EQ.1) then
            do t = 1, n
                call dgemv('n',p,m,1.0d0,zt(:,:,(t-1)*timevar(1)+1),p,alphahat(:,t),1,0.0d0,thetahat(:,t),1)
                call dsymm('r','u',p,m,1.0d0,alphavar(:,:,t),m,zt(:,:,(t-1)*timevar(1)+1),p,0.0d0,pm,p)
                call dgemm('n','t',p,p,m,1.0d0,pm,p,zt(:,:,(t-1)*timevar(1)+1),p,0.0d0,thetavar(:,:,t),p)
            end do
        end if

        if(smoothy.EQ.1) then
            allocate(osim(smootha*p+(1-smootha),n,4*nsim*smootha+(1-smootha)))
            do t = 1, n
                do j = 1,p
                    call dgemv('t',m,4*nsim,1.0d0,sim(:,t,:),m,zt(j,:,(t-1)*timevar(1)+1),1,0.0d0,osim(j,t,:),1)
                end do
            end do

            do j= 1,p
                select case(dist(j))
                    case(1)

                    case(2)
                        do t=1, n
                            osim(j,t,:) = exp(osim(j,t,:))*u(t,j)
                        end do
                    case(3)
                        osim(j,:,:) = exp(osim(j,:,:))/(1.0d0+exp(osim(j,:,:)))
                    case default
                        osim(j,:,:) = exp(osim(j,:,:))
                end select
            end do
            do t=1,n
                call covmeanw(osim(:,t,:),w(t,:),p,1,4*nsim,yhat(:,t),yvar(:,:,t))
            end do
            deallocate(osim)
        end if
    else
        call isamplefilter(yt, ymiss, timevar, zt, tt, rtv, qt, a1, p1,p1inf, u, dist, &
        p, n, m, r, theta, maxiter,rankp,convtol, nnd,nsim,epsplus,etaplus,&
        aplus1,c,tol,info,1,w,sim,5,p, expected, htol)

        if(info .ne. 0 .and. info .ne. 3) then
            return
        end if

        do t=1,n
            w(t,:) = w(t,:)/sum(w(t,:))
        end do

        if(smooths.EQ.1) then
            do t=1,n
                call covmeanwprotect(sim(:,t,:),w(t,:),p,1,4*nsim,thetahat(:,t),thetavar(:,:,t))
            end do
        end if

        if(smoothy.EQ.1) then
            do j= 1,p
                select case(dist(j))
                    case(1)

                    case(2)
                        do t=1, n
                            sim(j,t,:) = exp(sim(j,t,:))*u(t,j)
                        end do
                    case(3)
                        sim(j,:,:) = exp(sim(j,:,:))/(1.0d0+exp(sim(j,:,:)))
                    case default
                        sim(j,:,:) = exp(sim(j,:,:))
                end select
            end do
            do t=1,n
                call covmeanw(sim(:,t,:),w(t,:),p,1,4*nsim,yhat(:,t),yvar(:,:,t))
            end do
        end if
    end if
    deallocate(sim)
    deallocate(w)
end subroutine ngfilter
