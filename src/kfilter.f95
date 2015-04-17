  ! Subroutine for Kalman filtering of linear gaussian state space model

subroutine kfilter(yt, ymiss, timevar, zt, ht,tt, rt, qt, a1, p1, p1inf, p,n,m,r,d,j,&
at, pt, vt, ft,kt, pinf, finf, kinf, lik, tol,rankp,theta,thetavar,filtersignal)

    implicit none

    integer, intent(in) :: p, m, r, n,filtersignal
    integer, intent(inout) :: d, j, rankp
    integer :: t,tv
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(p,p,(n-1)*timevar(2)+1) :: ht
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rt
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) :: p1,p1inf
    double precision, intent(in) :: tol
    double precision, intent(inout), dimension(m,n+1) :: at
    double precision, intent(inout), dimension(m,m,n+1) :: pt,pinf
    double precision, intent(inout), dimension(p,n) :: vt,ft,finf
    double precision, intent(inout), dimension(m,p,n) :: kt,kinf
    double precision, intent(inout) :: lik
    double precision, intent(inout), dimension(p,p,n) :: thetavar
    double precision, intent(inout), dimension(n,p) :: theta
    double precision, dimension(m,r) :: mr
    double precision, dimension(p,m) :: pm
    double precision :: c,meps
    double precision, external :: ddot
    double precision, dimension(m,m,(n-1)*max(timevar(4),timevar(5))+1) :: rqr
    external dgemm, dsymm, dgemv, dsymv, dsyr, dsyr2

    meps = epsilon(meps)
    c = 0.5d0*log(8.0d0*atan(1.0d0))

    lik = 0.0d0
    tv= max(timevar(4),timevar(5))
    do t=1, (n-1)*tv+1
        call dsymm('r','l',m,r,1.0d0,qt(:,:,(t-1)*timevar(5)+1),r,rt(:,:,(t-1)*timevar(4)+1),m,0.0d0,mr,m)
        call dgemm('n','t',m,m,r,1.0d0,mr,m,rt(:,:,(t-1)*timevar(4)+1),m,0.0d0,rqr(:,:,t),m)
    end do

    j=0
    d=0
    pinf(:,:,1)=p1inf
    pt(:,:,1) = p1
    at(:,1) = a1
    ! diffuse initialization
    if(rankp .GT. 0) then
        diffuse: do while(d .LT. n .AND. rankp .GT. 0)
            d = d+1
            at(:,d+1) = at(:,d)
            pt(:,:,d+1) = pt(:,:,d)
            pinf(:,:,d+1) = pinf(:,:,d)
            call dfilter1step(ymiss(d,:),yt(d,:),transpose(zt(:,:,(d-1)*timevar(1)+1)),ht(:,:,(d-1)*timevar(2)+1),&
            tt(:,:,(d-1)*timevar(3)+1),rqr(:,:,(d-1)*tv+1),&
            at(:,d+1),pt(:,:,d+1),vt(:,d),ft(:,d),kt(:,:,d),pinf(:,:,d+1),finf(:,d),kinf(:,:,d),rankp,lik,tol,meps,c,p,m,j)
        end do diffuse


        if(rankp .EQ. 0 .AND. j .LT. p) then
                !non-diffuse filtering begins

            call filter1step(ymiss(d,:),yt(d,:),transpose(zt(:,:,(d-1)*timevar(1)+1)),ht(:,:,(d-1)*timevar(2)+1),&
            tt(:,:,(d-1)*timevar(3)+1),rqr(:,:,(d-1)*tv+1),&
            at(:,d+1),pt(:,:,d+1),vt(:,d),ft(:,d),kt(:,:,d),lik,tol,c,p,m,j)

        else
            j = p
        end if
    end if

    !Non-diffuse filtering continues from t=d+1, i=1


    do t = d+1, n
        at(:,t+1) = at(:,t)
        pt(:,:,t+1) = pt(:,:,t)
        call filter1step(ymiss(t,:),yt(t,:),transpose(zt(:,:,(t-1)*timevar(1)+1)),ht(:,:,(t-1)*timevar(2)+1),&
        tt(:,:,(t-1)*timevar(3)+1),rqr(:,:,(t-1)*tv+1),&
        at(:,t+1),pt(:,:,t+1),vt(:,t),ft(:,t),kt(:,:,t),lik,tol,c,p,m,0)
    end do

    if(filtersignal.EQ.1) then
        do t = 1, n
            call dgemv('n',p,m,1.0d0,zt(:,:,(t-1)*timevar(1)+1),p,at(:,t),1,0.0d0,theta(t,:),1)
            call dsymm('r','u',p,m,1.0d0,pt(:,:,t),m,zt(:,:,(t-1)*timevar(1)+1),p,0.0d0,pm,p)
            call dgemm('n','t',p,p,m,1.0d0,pm,p,zt(:,:,(t-1)*timevar(1)+1),p,0.0d0,thetavar(:,:,t),p)
        end do
    end if
end subroutine kfilter
