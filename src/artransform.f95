subroutine artransform(u,phi,p)
    implicit none
 
    integer, intent(in) :: p   
    integer :: i, j
    double precision, intent(inout), dimension(p) :: u
    double precision, intent(inout), dimension(p,p) :: phi

    do i= 2, p
        do j= 1, i-1
           phi(i,j) = phi(i-1,j) - u(i)*phi(i-1,i-j)
        end do
    end do
    
end subroutine artransform