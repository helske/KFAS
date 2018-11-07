! Subroutine for computation of p(y|theta)

subroutine pytheta(theta, dist, u, yt, ymiss, dev, p, n)

    implicit none

    integer, intent(in) ::  p,n
    integer, intent(in), dimension(p) :: dist
    integer, intent(in), dimension(p,n) :: ymiss
    integer ::  t,j
    double precision, intent(in), dimension(p,n) :: u
    double precision, intent(in), dimension(p,n) :: yt
    double precision, intent(in), dimension(p,n) :: theta
    double precision, intent(inout) :: dev

    external dpoisf, dbinomf, dgammaf, dnbinomf

    dev = 0.0d0
    do j=1,p
        select case(dist(j))
            case(2)
                do t=1,n
                    if(ymiss(j,t).EQ.0) then
                        call dpoisf(yt(j,t), u(j,t)*exp(theta(j,t)), dev)
                    end if
                end do
            case(3)
                do t=1,n
                    if(ymiss(j,t).EQ.0) then
                        call dbinomf(yt(j,t), u(j,t), exp(theta(j,t))/(1.0d0+exp(theta(j,t))), dev)
                    end if
                end do
            case(4)
                do t=1,n
                    if(ymiss(j,t).EQ.0) then
                        call dgammaf(yt(j,t), u(j,t), exp(theta(j,t))/u(j,t), dev)
                    end if
                end do
            case(5)
                do t=1,n
                    if(ymiss(j,t).EQ.0) then
                        call dnbinomf(yt(j,t), u(j,t), exp(theta(j,t)), dev)
                    end if
                end do
        end select
    end do
end subroutine pytheta
