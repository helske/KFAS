! disturbance smoothing algorithm for simulation
subroutine smoothsim(yt, ymiss, timevar, zt, ht,tt, rtv,qt,rqr, a1, p1, p1inf, &
d, j, p, m, n, r,tol,rankp,ft,finf,kt,kinf,epshat,etahat,rt0,rt1,needeps)

    implicit none

    logical, intent(in) :: needeps
    integer, intent(in) ::  p, m, n,r
    integer, intent(inout) :: d, j,rankp
    integer ::  t, i,tv
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(p,p,(n-1)*timevar(2)+1) :: ht
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rtv
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, dimension(m,m,(n-1)*max(timevar(4),timevar(5))+1) :: rqr
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, intent(inout), dimension(p,n) :: ft,finf
    double precision, intent(inout), dimension(m,p,n) :: kt,kinf
    double precision, dimension(p,n) :: vt
    double precision, dimension(m) :: at, arec
    double precision, dimension(m,m) :: pt,pinf,mm
    double precision, intent(in) :: tol
    double precision, dimension(m) :: rrec,rrec1,rhelp,help
    double precision, dimension(m,m) :: im,linf,l0
    double precision, intent(inout), dimension(r,n) :: etahat
    double precision, intent(inout), dimension(p,n) :: epshat
    double precision, intent(inout), dimension(m) :: rt0,rt1
    double precision :: meps, finv,lik
    double precision, external :: ddot

    external dgemm, dsymm, dgemv, dsymv, dsyr, dsyr2, dger
 
    meps = epsilon(meps)
    tv = max(timevar(4),timevar(5))
    lik = 0.0d0
    im = 0.0d0
    do i = 1, m
        im(i,i) = 1.0d0
    end do
    j=0
    d=0
    pinf = p1inf
    pt = p1
    at = a1
    if(rankp .GT. 0) then
        diffuse: do while(d .LT. n .AND. rankp .GT. 0)
            d = d+1
            call diffusefilteronestep(ymiss(d,:),yt(d,:),transpose(zt(:,:,(d-1)*timevar(1)+1)),ht(:,:,(d-1)*timevar(2)+1),&
            tt(:,:,(d-1)*timevar(3)+1),rqr(:,:,(d-1)*tv+1),&
            at,pt,vt(:,d),ft(:,d),kt(:,:,d),pinf,finf(:,d),kinf(:,:,d),rankp,lik,tol,meps,0.0d0,p,m,j)
        end do diffuse
        if(rankp .EQ. 0 .AND. j .LT. p) then
            call filteronestep(ymiss(d,:),yt(d,:),transpose(zt(:,:,(d-1)*timevar(1)+1)),ht(:,:,(d-1)*timevar(2)+1),&
            tt(:,:,(d-1)*timevar(3)+1),rqr(:,:,(d-1)*tv+1),&
            at,pt,vt(:,d),ft(:,d),kt(:,:,d),lik,tol,meps,0.0d0,p,m,j)

        else
            j = p
        end if

    end if

    !Non-diffuse filtering continues from t=d+1, i=1

    do t = d+1, n
        call filteronestep(ymiss(t,:),yt(t,:),transpose(zt(:,:,(t-1)*timevar(1)+1)),ht(:,:,(t-1)*timevar(2)+1),&
        tt(:,:,(t-1)*timevar(3)+1),rqr(:,:,(t-1)*tv+1),&
        at,pt,vt(:,t),ft(:,t),kt(:,:,t),lik,tol,meps,0.0d0,p,m,0)
    end do

    !smoothing begins

    rt0 = 0.0d0

    do t = n, d+1, -1
        call smoothonestep(ymiss(t,:), transpose(zt(:,:,(t-1)*timevar(1)+1)), ht(:,:,(t-1)*timevar(2)+1), &
        tt(:,:,(t-1)*timevar(3)+1), rtv(:,:,(t-1)*timevar(4)+1), qt(:,:,(t-1)*timevar(5)+1), vt(:,t), &
        ft(:,t),kt(:,:,t), im,p,m,r,1,rt0,etahat(:,t),epshat(:,t),needeps)

    end do

    if(d .GT. 0) then
        t = d
        if(j .LT. p) then
            call smoothonestep(ymiss(t,:), transpose(zt(:,:,(t-1)*timevar(1)+1)), ht(:,:,(t-1)*timevar(2)+1), &
            tt(:,:,(t-1)*timevar(3)+1), rtv(:,:,(t-1)*timevar(4)+1), qt(:,:,(t-1)*timevar(5)+1), vt(:,t), &
            ft(:,t),kt(:,:,t), im,p,m,r,j+1,rt0,etahat(:,t),epshat(:,t),needeps)
        end if


        rt1 = 0.0d0
        call diffusesmoothonestep(ymiss(t,:), transpose(zt(:,:,(t-1)*timevar(1)+1)), ht(:,:,(t-1)*timevar(2)+1), &
        tt(:,:,(t-1)*timevar(3)+1), rtv(:,:,(t-1)*timevar(4)+1), qt(:,:,(t-1)*timevar(5)+1), vt(:,t), &
        ft(:,t),kt(:,:,t), im,p,m,r,j,rt0,rt1,finf(:,t),kinf(:,:,t),etahat(:,t),epshat(:,t),needeps)


        do t=(d-1), 1, -1

        call diffusesmoothonestep(ymiss(t,:), transpose(zt(:,:,(t-1)*timevar(1)+1)), ht(:,:,(t-1)*timevar(2)+1), &
        tt(:,:,(t-1)*timevar(3)+1), rtv(:,:,(t-1)*timevar(4)+1), qt(:,:,(t-1)*timevar(5)+1), vt(:,t), &
        ft(:,t),kt(:,:,t), im,p,m,r,p,rt0,rt1,finf(:,t),kinf(:,:,t),etahat(:,t),epshat(:,t),needeps)



        end do
    end if

end subroutine smoothsim
