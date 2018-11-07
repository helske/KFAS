!fast disturbance smoother for simulation
subroutine smoothsimfast(yt, ymiss, timevar, zt, ht,tt, rtv,qt,a1, ft,kt,&
finf, kinf, dt, jt, p, m, n,r,epshat,etahat,rt0,rt1,needeps)

    implicit none

    logical, intent(in) :: needeps
    integer, intent(in) ::  p, m, r,n,dt,jt
    integer ::  t, i
    integer, intent(in), dimension(p,n) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    double precision, intent(in), dimension(p,n) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(p,p,(n-1)*timevar(2)+1) :: ht
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rtv
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(p,n) :: ft,finf
    double precision, intent(in), dimension(m,p,n) :: kt,kinf
    double precision, intent(inout), dimension(p,n) :: epshat
    double precision, intent(inout), dimension(r,n) :: etahat
    double precision, dimension(p,n) :: vt
    double precision, dimension(m) :: at
    double precision, dimension(m,m) :: im
    double precision, intent(inout), dimension(m) :: rt0,rt1
    double precision :: lik
    lik = 0.0d0
    at = a1
    if(dt .GT. 0) then
        !diffuse filtering begins
        do t = 1, dt - 1
            call dfilter1stepnv(ymiss(:,t),yt(:,t),transpose(zt(:,:,(t-1)*timevar(1)+1)),&
            tt(:,:,(t-1)*timevar(3)+1),at,vt(:,t),ft(:,t),kt(:,:,t),&
            finf(:,t),kinf(:,:,t),p,m,p,lik)
        end do

        t = dt
        call dfilter1stepnv(ymiss(:,t),yt(:,t),transpose(zt(:,:,(t-1)*timevar(1)+1)),&
        tt(:,:,(t-1)*timevar(3)+1),at,vt(:,t),ft(:,t),kt(:,:,t),&
        finf(:,t),kinf(:,:,t),p,m,jt,lik)
        !non-diffuse filtering begins
        if(jt .LT. p) then
            call filter1stepnv(ymiss(:,t),yt(:,t),transpose(zt(:,:,(t-1)*timevar(1)+1)),&
            tt(:,:,(t-1)*timevar(3)+1),at,vt(:,t),ft(:,t),kt(:,:,t),p,m,jt,lik)
        end if
    end if
    !Non-diffuse filtering continues from t=d+1, i=1
    do t = dt + 1, n
        call filter1stepnv(ymiss(:,t),yt(:,t),transpose(zt(:,:,(t-1)*timevar(1)+1)),&
        tt(:,:,(t-1)*timevar(3)+1),at,vt(:,t),ft(:,t),kt(:,:,t),p,m,0,lik)
    end do


    !smoothing begins

    im = 0.0d0
    do i = 1, m
        im(i,i) = 1.0d0
    end do

    rt0 = 0.0d0

    do t = n, dt+1, -1
        call smooth1step(ymiss(:,t), transpose(zt(:,:,(t-1)*timevar(1)+1)), ht(:,:,(t-1)*timevar(2)+1), &
        tt(:,:,(t-1)*timevar(3)+1), rtv(:,:,(t-1)*timevar(4)+1), qt(:,:,(t-1)*timevar(5)+1), vt(:,t), &
        ft(:,t),kt(:,:,t), im,p,m,r,1,rt0,etahat(:,t),epshat(:,t),needeps)
    end do

    if(dt .GT. 0) then
        t = dt
        if(jt .LT. p) then
            call smooth1step(ymiss(:,t), transpose(zt(:,:,(t-1)*timevar(1)+1)), ht(:,:,(t-1)*timevar(2)+1), &
            tt(:,:,(t-1)*timevar(3)+1), rtv(:,:,(t-1)*timevar(4)+1), qt(:,:,(t-1)*timevar(5)+1), vt(:,t), &
            ft(:,t),kt(:,:,t), im,p,m,r,jt+1,rt0,etahat(:,t),epshat(:,t),needeps)
        end if
        rt1 = 0.0d0
        call dsmooth1step(ymiss(:,t), transpose(zt(:,:,(t-1)*timevar(1)+1)), ht(:,:,(t-1)*timevar(2)+1), &
        tt(:,:,(t-1)*timevar(3)+1), rtv(:,:,(t-1)*timevar(4)+1), qt(:,:,(t-1)*timevar(5)+1), vt(:,t), &
        ft(:,t),kt(:,:,t), im,p,m,r,jt,rt0,rt1,finf(:,t),kinf(:,:,t),etahat(:,t),epshat(:,t),needeps)
        do t = (dt - 1), 1, -1
            call dsmooth1step(ymiss(:,t), transpose(zt(:,:,(t-1)*timevar(1)+1)), ht(:,:,(t-1)*timevar(2)+1), &
            tt(:,:,(t-1)*timevar(3)+1), rtv(:,:,(t-1)*timevar(4)+1), qt(:,:,(t-1)*timevar(5)+1), vt(:,t), &
            ft(:,t),kt(:,:,t), im,p,m,r,p,rt0,rt1,finf(:,t),kinf(:,:,t),etahat(:,t),epshat(:,t),needeps)
        end do
    end if

end subroutine smoothsimfast
