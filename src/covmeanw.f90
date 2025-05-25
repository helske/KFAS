! Subroutines for computation of mean, covariance and variance from weighted sample
! All subroutines assume that the weights are normalized i.e. sum(w)=1
subroutine covmeanw(x,w,m,n,k,meanx,covx)

    implicit none
    integer, intent(in) :: m, n, k
    integer :: t,i
    double precision, intent(inout), dimension(m,n,k) :: x
    double precision, intent(in), dimension(k) :: w
    double precision, intent(inout), dimension(m,n) :: meanx
    double precision, intent(inout), dimension(m,m,n) :: covx

    external dgemm

    do i = 1, k
        meanx = meanx + x(:,:,i)*w(i)
    end do
    do i = 1, k
        x(:,:,i) = sqrt(w(i))*(x(:,:,i) - meanx)
    end do

    do t = 1, n
        call dgemm('n','t',m,m,k,1.0d0,x(:,t,:),m,x(:,t,:),m,0.0d0,covx(:,:,t),m)
    end do

end subroutine covmeanw

subroutine covmeanwprotect(x,w,m,n,k,meanx,covx)

    implicit none
    integer, intent(in) :: m, n, k
    integer :: t,i,j,l
    double precision, intent(in), dimension(m,n,k) :: x
    double precision, intent(in), dimension(k) :: w
    double precision, intent(inout), dimension(m,n) :: meanx
    double precision, intent(inout), dimension(m,m,n) :: covx
    double precision, dimension(:,:), allocatable :: x2

    external dgemm
     do i = 1, k
        do j = 1, m
          do l = 1,n
              meanx(j, l) = meanx(j,l)+x(j,l,i)*w(i)
          end do
        end do
    end do
    allocate(x2(m,k))
    do t = 1, n
        do i = 1, k
          do j = 1, m
              x2(j, i) = sqrt(w(i))*(x(j,t,i) - meanx(j,t))
          end do
        end do
        call dgemm('n','t',m,m,k,1.0d0,x2,m,x2,m,0.0d0,covx(:,:,t),m)
    end do
    deallocate(x2)
end subroutine covmeanwprotect

subroutine varmeanw(x,w,m,n,k,meanx,varx,var)

    implicit none
    integer, intent(in) :: m, n, k,var
    integer :: t,i
    double precision, intent(inout), dimension(n,m,k) :: x
    double precision, intent(in), dimension(k) :: w
    double precision, intent(inout), dimension(n,m) :: meanx
    double precision, intent(inout), dimension(n,m) :: varx


    do i = 1, k
        meanx = meanx + x(:,:,i)*w(i)
    end do
    if(var.EQ.1) then
        do i = 1, m
            do t = 1, n
                varx(t,i) = sum(w*x(t,i,:)**2)-meanx(t,i)**2
            end do
        end do
    end if
end subroutine varmeanw
