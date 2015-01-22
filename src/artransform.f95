!transformation of unconstrained parameters to stationary region
subroutine artransform(u,phi,p)
    implicit none
 
    integer, intent(in) :: p   
    integer :: i, j
    double precision, intent(inout), dimension(p) :: u
    double precision, intent(inout), dimension(p) :: phi
    double precision :: a

    do i = 2, p
        a = phi(i)
        do j =1, i
            u(j) =  u(j) - a * phi(i - j)
        end do
        do j = 1, i
            phi(j) = u(j)
        end do
    end do


    
end subroutine artransform

