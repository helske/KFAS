! Subroutine for computation of the approximating gaussian model for non-gaussian models

subroutine approx(yt, ymiss, timevar, zt, tt, rtv, ht, qt, a1, p1,p1inf, p,n,m,r,&
theta, u, ytilde, dist,maxiter,tol,rankp,convtol,diff)

    implicit none

    integer, intent(in) ::  p,m, r, n
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    integer, intent(in), dimension(p) :: dist
    integer, intent(inout) :: maxiter,rankp
    integer ::  j,i, k,rankp2,tv
    double precision, intent(in) :: tol,convtol
    double precision, intent(in), dimension(n,p) :: u
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rtv
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, intent(inout), dimension(n,p) :: theta
    double precision, intent(inout), dimension(n,p) :: ytilde
    double precision, intent(inout), dimension(p,p,n) :: ht
    double precision, intent(inout) :: diff
    double precision, dimension(m,r) :: mr
    double precision, dimension(m,m,(n-1)*max(timevar(4),timevar(5))+1) :: rqr
    double precision dev, devold
    double precision, dimension(n,p) :: muhat
    double precision, external :: ddot


    !compute rqr
    tv = max(timevar(4),timevar(5))
    do i=1, (n-1)*tv+1
        call dsymm('r','u',m,r,1.0d0,qt(:,:,(i-1)*timevar(5)+1),r,rtv(:,:,(i-1)*timevar(4)+1),m,0.0d0,mr,m)
        call dgemm('n','t',m,m,r,1.0d0,mr,m,rtv(:,:,(i-1)*timevar(4)+1),m,0.0d0,rqr(:,:,i),m)
    end do

    diff = 1000.0d0
    k=0


    do j=1,p
        select case(dist(j))
            case(1)
                do i=1,n
                    if(ymiss(i,j).EQ.0) then
                        ht(j,j,i) =  u(i,j)
                        ytilde(i,j) =  yt(i,j)
                    end if
                end do
            case(2)
                do i=1,n
                    if(ymiss(i,j).EQ.0) then
                        ht(j,j,i) =  1.0d0/(exp(theta(i,j))*u(i,j))
                        ytilde(i,j) =  yt(i,j)*ht(j,j,i) + theta(i,j) - 1.0d0
                    end if
                end do
            case(3)
                do i=1,n
                    if(ymiss(i,j).EQ.0) then
                        ht(j,j,i) = (1.0d0+exp(theta(i,j)))**2/(u(i,j)*exp(theta(i,j)))
                        ytilde(i,j) = theta(i,j) + ht(j,j,i)*yt(i,j) - 1.0d0 - exp(theta(i,j))
                    end if
                end do
            case(4)
                do i=1,n
                    if(ymiss(i,j).EQ.0) then
                        !ht(j,j,i) = exp(theta(i,j))/(u(i,j)*yt(i,j))
                        !ytilde(i,j) = theta(i,j)+1.0d0-exp(theta(i,j))/yt(i,j)
                        ht(j,j,i) =1.0d0/u(i,j) !1.0d0
                        ytilde(i,j) = theta(i,j)+yt(i,j)/exp(theta(i,j))-1.0d0
                    end if
                end do
            case(5)
                do i=1,n
                    if(ymiss(i,j).EQ.0) then
                        !ht(j,j,i) = (exp(theta(i,j))+u(i,j))**2/(u(i,j)*exp(theta(i,j))*(yt(i,j)+u(i,j)))
                        !ytilde(i,j) = theta(i,j) + ht(j,j,i)*u(i,j)*(yt(i,j)-exp(theta(i,j)))/(u(i,j)+exp(theta(i,j)))
                        !ht(j,j,i) = u(i,j)*exp(theta(i,j))/(u(i,j)+exp(theta(i,j)))
                        ht(j,j,i) = (1.0d0/u(i,j)+1.0d0/exp(theta(i,j)))
                        ytilde(i,j) = theta(i,j)+yt(i,j)/exp(theta(i,j))-1.0d0
                    end if
                end do
        end select
    end do


    muhat = theta
        call mu(dist,u,n,p,muhat)
        call deviance(yt,muhat,u,ymiss,n,p,dist,devold)

    do while(diff > convtol .AND. k < maxiter)

        k=k+1


        rankp2 = rankp

        call kfstheta(ytilde, ymiss, timevar, zt, ht,tt, rtv,qt,rqr, a1, p1, p1inf, &
        p, n, m, r,tol,rankp2,theta)

        do j=1,p
            select case(dist(j))
                case(2)
                    do i=1,n
                        if(ymiss(i,j).EQ.0) then
                            ht(j,j,i) =  1.0d0/(exp(theta(i,j))*u(i,j))
                            ytilde(i,j) =  yt(i,j)*ht(j,j,i) + theta(i,j) - 1.0d0
                        end if
                    end do
                case(3)
                    do i=1,n
                        if(ymiss(i,j).EQ.0) then
                            ht(j,j,i) = (1.0d0+exp(theta(i,j)))**2/(u(i,j)*exp(theta(i,j)))
                            ytilde(i,j) = theta(i,j) + ht(j,j,i)*yt(i,j) - 1.0d0 - exp(theta(i,j))
                        end if
                    end do
                case(4)
                    do i=1,n
                        if(ymiss(i,j).EQ.0) then
                               ! ht(j,j,i) = exp(theta(i,j))/(u(i,j)*yt(i,j))
                               ! ytilde(i,j) = theta(i,j)+1.0d0-exp(theta(i,j))/yt(i,j)
                               !ht(j,j,i) = 1.0d0
                            ytilde(i,j) = theta(i,j)+yt(i,j)/exp(theta(i,j))-1.0d0
                        end if
                    end do
                case(5)
                    do i=1,n
                        if(ymiss(i,j).EQ.0) then
                                !ht(j,j,i) = (exp(theta(i,j))+u(i,j))**2/(u(i,j)*exp(theta(i,j))*(yt(i,j)+u(i,j)))
                                !ytilde(i,j) = theta(i,j) + ht(j,j,i)*u(i,j)*(yt(i,j)-exp(theta(i,j)))/(u(i,j)+exp(theta(i,j)))
                            !ht(j,j,i) = u(i,j)*exp(theta(i,j))/(u(i,j)+exp(theta(i,j)))
                            ht(j,j,i) = (1.0d0/u(i,j)+1.0d0/exp(theta(i,j)))
                            ytilde(i,j) = theta(i,j)+yt(i,j)/exp(theta(i,j))-1.0d0
                        end if
                    end do
            end select
        end do
        muhat = theta
        call mu(dist,u,n,p,muhat)
        call deviance(yt,muhat,u,ymiss,n,p,dist,dev)
        diff = abs(dev - devold)/(0.1d0 + abs(dev))
        devold=dev
    end do
    maxiter=k

end subroutine approx
