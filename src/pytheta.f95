! Subroutine for computation of p(y|theta)

subroutine pytheta(theta, dist, u, yt, ymiss, dev, p, n)

    implicit none

    integer, intent(in) ::  p,n
    integer, intent(in), dimension(p) :: dist
    integer, intent(in), dimension(n,p) :: ymiss
    integer ::  t,j
    double precision, intent(in), dimension(n,p) :: u
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(n,p) :: theta
    double precision, intent(inout) :: dev

    dev = 0.0d0
    do j=1,p
        select case(dist(j))
            case(2)
                do t=1,n
                    if(ymiss(t,j).EQ.0) then
                        call dpoisf(yt(t,j), u(t,j)*exp(theta(t,j)), dev)
                    end if
                end do
            case(3)
                do t=1,n
                    if(ymiss(t,j).EQ.0) then
                        call dbinomf(yt(t,j), u(t,j), exp(theta(t,j))/(1.0d0+exp(theta(t,j))), dev)
                    end if
                end do
            case(4)
                do t=1,n
                    if(ymiss(t,j).EQ.0) then
                        call dgammaf(yt(t,j), u(t,j), exp(theta(t,j))/u(t,j), dev)
                    end if
                end do
            case(5)
                do t=1,n
                    if(ymiss(t,j).EQ.0) then
                        call dnbinomf(yt(t,j), u(t,j), exp(theta(t,j)), dev)
                    end if
                end do
        end select
    end do
end subroutine pytheta
