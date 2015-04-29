!transformation of unconstrained parameters to stationary region
subroutine artransform(phi,p)
    implicit none
 
    integer, intent(in) :: p   
    integer :: i, j
    double precision, intent(inout), dimension(p) :: phi
    double precision, dimension(p,p) :: u

    u = 0.0d0
    do i = 1,p
        u(i,i) = phi(i)
    end do
    do i= 2, p
        do j= 1, i-1
            u(i,j) = u(i-1,j) - phi(i)*u(i-1,i-j)
        end do
    end do
    phi = u(p,:)
end subroutine artransform
