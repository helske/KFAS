! subroutine for computing new guess of conditional mode theta for non-Gaussian model given current guess

subroutine approxloop(yt, ymiss, timevar, zt, tt, rtv, ht, qt, rqr, tvrqr, a1, p1,p1inf, p,n,m,r, &
theta, thetanew, u, ytilde, dist,tol,rankp,lik)

    implicit none

    integer, intent(in) ::  p,m, r, n,tvrqr,rankp
    integer, intent(in), dimension(p,n) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    integer, intent(in), dimension(p) :: dist
    integer ::  j,i, rankp2
    double precision, intent(in) :: tol
    double precision, intent(in), dimension(p,n) :: u
    double precision, intent(in), dimension(p,n) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rtv
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, intent(in), dimension(p,n) :: theta
    double precision, intent(inout), dimension(p,n) :: thetanew
    double precision, intent(inout), dimension(p,n) :: ytilde
    double precision, intent(inout), dimension(p,p,n) :: ht
    double precision, intent(in), dimension(m,m,tvrqr) :: rqr
    double precision, intent(inout) :: lik
    double precision, external :: ddot

    external kfstheta

    !construct pseudo-observations and variances for gaussian approximating model
    do j=1,p
        select case(dist(j))
            case(1)
                do i=1,n
                    if(ymiss(j,i).EQ.0) then
                        ht(j,j,i) =  u(j,i)
                        ytilde(j,i) =  yt(j,i)
                    end if
                end do
            case(2)
                do i=1,n
                    if(ymiss(j,i).EQ.0) then
                        ht(j,j,i) =  1.0d0/(exp(theta(j,i))*u(j,i))
                        ytilde(j,i) =  yt(j,i)*ht(j,j,i) + theta(j,i) - 1.0d0
                    end if
                end do
            case(3)
                do i=1,n
                    if(ymiss(j,i).EQ.0) then
                        ht(j,j,i) = (1.0d0+exp(theta(j,i)))**2/(u(j,i)*exp(theta(j,i)))
                        ytilde(j,i) = theta(j,i) + ht(j,j,i)*yt(j,i) - 1.0d0 - exp(theta(j,i))
                    end if
                end do
            case(4)
                do i=1,n
                    if(ymiss(j,i).EQ.0) then
                        ht(j,j,i) =1.0d0/u(j,i)
                        ytilde(j,i) = theta(j,i)+yt(j,i)/exp(theta(j,i))-1.0d0
                    end if
                end do
            case(5)
                do i=1,n
                    if(ymiss(j,i).EQ.0) then
                        ht(j,j,i) = (1.0d0/u(j,i)+1.0d0/exp(theta(j,i)))
                        ytilde(j,i) = theta(j,i)+yt(j,i)/exp(theta(j,i))-1.0d0
                    end if
                end do
        end select
    end do

    ! compute new estimate of thetahat
    rankp2 = rankp
    call kfstheta(ytilde, ymiss, timevar, zt, ht,tt, rtv,qt,rqr, tvrqr, a1, p1, p1inf, &
    p, n, m, r,tol,rankp2,thetanew,lik)



end subroutine approxloop
