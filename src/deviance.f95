subroutine mu(dist,u,n,p,theta)
    implicit none

    integer, intent(in) :: n,p
    integer, intent(in), dimension(p) :: dist
    double precision, intent(in), dimension(n,p) :: u
    double precision, intent(inout), dimension(n,p) :: theta
    integer j

    do j=1, p
        select case(dist(j))
            case(2)
                theta(:,j) = exp(theta(:,j))*u(:,j)
            case(3)
                theta(:,j) = exp(theta(:,j))/(1.0d0+exp(theta(:,j)))
            case(4)
                theta(:,j) = exp(theta(:,j))
            case(5)
                theta(:,j) = exp(theta(:,j))
        end select
    end do

end subroutine mu

subroutine deviance(y,mu,u,ymiss,n,p,dist,dev)

    implicit none

    integer, intent(in) :: n,p
    integer, intent(in), dimension(p) :: dist
    integer j
    integer, intent(in), dimension(n,p) :: ymiss
    double precision, intent(in), dimension(n,p) :: y, mu,u
    double precision, intent(inout) :: dev
    double precision, dimension(n,p) :: res
    double precision, dimension(n) :: tmp,tmp2

    res = y
    where(ymiss /= 0) res = 0.0d0

    do j=1, p
        select case(dist(j))
            case(1)
                where(ymiss(:,j) == 0) res(:,j) = (res(:,j) - mu(:,j))**2
            case(2)
                tmp = 1.0d0
                where (res(:,j)/=0.0d0) tmp = res(:,j)/mu(:,j)
                where(ymiss(:,j) == 0) res(:,j) =  2.0d0*(res(:,j)*log(tmp) - res(:,j) + mu(:,j))
            case(3)
                where(ymiss(:,j) == 0) res(:,j) = res(:,j)/u(:,j)

                tmp = 1.0d0
                where (res(:,j)/=0.0d0) tmp = res(:,j)/mu(:,j)
                tmp2 = 1.0d0
                where (res(:,j)/=1.0d0 .and. mu(:,j)/=1.0d0) tmp2 = (1.0d0-res(:,j))/(1.0d0-mu(:,j))

                where(ymiss(:,j) == 0) res(:,j) =  2.0d0*u(:,j)*(res(:,j)*log(tmp)+(1.0d0-res(:,j))*log(tmp2))

            case(4)
                tmp = 1.0d0
                where (res(:,j)/=0.0d0) tmp = res(:,j)/mu(:,j)
                where(ymiss(:,j) == 0) res(:,j) = -2.0d0*(log(tmp)-(res(:,j)-mu(:,j))/mu(:,j))
            case(5)
                tmp = res(:,j)/mu(:,j)
                where (res(:,j)<1.0d0) tmp = 1.0d0/mu(:,j)
                where(ymiss(:,j) == 0) res(:,j) = 2.0d0*(res(:,j)*log(tmp) - &
                (res(:,j)+u(:,j))*log((res(:,j)+u(:,j))/(mu(:,j)+u(:,j))))
        end select
    end do
    dev = sum(res)
end subroutine deviance
