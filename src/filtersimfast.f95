!fast filtering algorithm used in simulation filter
subroutine filtersimfast(yt, ymiss, timevar, zt,tt, &
a1, ft,kt,finf, kinf, dt, jt, p, m, n,at)

    implicit none

    integer, intent(in) ::  p, m,n,dt,jt
    integer ::  t
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(p,n) :: ft,finf
    double precision, intent(in), dimension(m,p,n) :: kt,kinf
    double precision, intent(inout), dimension(m,n+1) :: at
    double precision, dimension(p,n) :: vt
    double precision :: lik
    double precision, external :: ddot

    external dgemv

    lik = 0.0d0

    at(:,1) = a1
    if(dt .GT. 0) then
        !diffuse filtering begins
        do t = 1, (dt - 1)
            at(:,t+1) = at(:,t)
            call dfilter1stepnv(ymiss(t,:),yt(t,:),transpose(zt(:,:,(t-1)*timevar(1)+1)),&
            tt(:,:,(t-1)*timevar(3)+1),at(:,t+1),vt(:,t),ft(:,t),kt(:,:,t),&
            finf(:,t),kinf(:,:,t),p,m,p,lik)
        end do

        t = dt
        at(:,t+1) = at(:,t)
        call dfilter1stepnv(ymiss(t,:),yt(t,:),transpose(zt(:,:,(t-1)*timevar(1)+1)),&
        tt(:,:,(t-1)*timevar(3)+1),at(:,t+1),vt(:,t),ft(:,t),kt(:,:,t),&
        finf(:,t),kinf(:,:,t),p,m,jt,lik)
        !non-diffuse filtering begins
        if(jt .LT. p) then
            call filter1stepnv(ymiss(t,:),yt(t,:),transpose(zt(:,:,(t-1)*timevar(1)+1)),&
            tt(:,:,(t-1)*timevar(3)+1),at(:,t+1),vt(:,t),ft(:,t),kt(:,:,t),p,m,jt,lik)
        end if
    end if
    !Non-diffuse filtering continues from t=d+1, i=1
    do t = dt + 1, n
        at(:,t+1) = at(:,t)
        call filter1stepnv(ymiss(t,:),yt(t,:),transpose(zt(:,:,(t-1)*timevar(1)+1)),&
        tt(:,:,(t-1)*timevar(3)+1),at(:,t+1),vt(:,t),ft(:,t),kt(:,:,t),p,m,0,lik)
    end do


end subroutine filtersimfast
