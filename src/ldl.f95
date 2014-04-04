! LDL decomposition
subroutine ldl(a,n,tol,info)
  
    implicit none
  
    integer, intent(in) :: n
    integer, intent(inout) :: info
    integer :: i,j,k
    double precision :: di,tmp
    double precision, intent(inout), dimension(n,n) :: a
    double precision, intent(in) :: tol

    do i = 1, n
        di=a(i,i)
        if(abs(di)<=tol) then
            !a(i,:) = 0.0d0
            a(:,i) = 0.0d0
        else
            do j = i+1, n
                tmp = a(j,i)/di
                a(j,i) = tmp
                a(j,j) = a(j,j) - tmp**2*di
                do k = j+1, n
                    a(k,j) = a(k,j) - tmp*a(k,i)
                end do
            end do
        end if
    end do
    do i = 1,n
        a(i,(i+1):n) = 0.0d0
        if(a(i,i)<0.0d0) then
            info = -1
        end if
    end do

end subroutine ldl
