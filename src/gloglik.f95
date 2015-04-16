! Subroutine for computing the log-Likelihood of general linear gaussian state space model

subroutine gloglik(yt, ymiss, timevar, zt, ht, tt, rt, qt, a1, p1, p1inf,&
p, m, r, n, lik, tol,rankp,marginal)


    implicit none

    integer, intent(in) ::  p, m, r, n
    integer, intent(inout) :: rankp,marginal
    integer ::  t, i,d,j,tv
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(p,p,(n-1)*timevar(2)+1) :: ht
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rt
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, intent(in) :: tol
    double precision, intent(inout) :: lik
    double precision, dimension(m) :: at,arec
    double precision, dimension(p) :: vt,ft,finf
    double precision, dimension(m,p) :: kt,kinf
    double precision, dimension(m,m) :: pt,pinf,mm
    double precision, dimension(m,r) :: mr    
    double precision :: c, meps, finv
    double precision, external :: ddot
    double precision, dimension(m,m,(n-1)*max(timevar(4),timevar(5))+1) :: rqr

    external dgemm, dsymm, dgemv, dsymv, dsyr, dsyr2, marginalxx

    meps = epsilon(meps)

    ! compute RQR'
    tv= max(timevar(4),timevar(5))
    do t=1, (n-1)*tv+1
        call dsymm('r','l',m,r,1.0d0,qt(:,:,(t-1)*timevar(5)+1),r,rt(:,:,(t-1)*timevar(4)+1),m,0.0d0,mr,m)
        call dgemm('n','t',m,m,r,1.0d0,mr,m,rt(:,:,(t-1)*timevar(4)+1),m,0.0d0,rqr(:,:,t),m)
    end do


    ! constant term for log-likelihood
    c = 0.5d0*log(8.0d0*atan(1.0d0))
    lik = 0.0d0

    j=0
    d=0
    pt = p1
    at = a1
    pinf=p1inf
    ! Diffuse initialization
    if(rankp .GT. 0) then
        diffuse: do while(d .LT. n .AND. rankp .GT. 0)
            d = d+1
           
            call diffusefilteronestep(ymiss(d,:),yt(d,:),transpose(zt(:,:,(d-1)*timevar(1)+1)),ht(:,:,(d-1)*timevar(2)+1),&
            tt(:,:,(d-1)*timevar(3)+1),rqr(:,:,(d-1)*tv+1), at,pt,vt,ft,kt,pinf,finf,kinf,rankp,lik,tol,meps,c,p,m,j)

        end do diffuse

        if(rankp .EQ. 0 .AND. j .LT. p) then
            !non-diffuse filtering begins
            call filteronestep(ymiss(d,:),yt(d,:),transpose(zt(:,:,(d-1)*timevar(1)+1)),ht(:,:,(d-1)*timevar(2)+1),&
            tt(:,:,(d-1)*timevar(3)+1),rqr(:,:,(d-1)*tv+1),at,pt,vt,ft,kt,lik,tol,meps,c,p,m,j)

        else
            j = p

        end if
    end if



    !Non-diffuse filtering continues from t=d+1, i=1

    do t = d+1, n
        call filteronestep(ymiss(t,:),yt(t,:),transpose(zt(:,:,(t-1)*timevar(1)+1)),ht(:,:,(t-1)*timevar(2)+1),&
        tt(:,:,(t-1)*timevar(3)+1),rqr(:,:,(t-1)*tv+1),at,pt,vt,ft,kt,lik,tol,meps,c,p,m,0)

    end do

    if(marginal.EQ.1) then
        t = int(sum(p1inf))
        if(t.GT.0) then
            call marginalxx(p1inf,zt,tt,m,p,n,t,timevar,lik,marginal)
        end if
    end if
   

end subroutine gloglik
