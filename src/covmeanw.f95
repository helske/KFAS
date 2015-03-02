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
    integer :: t,i
    double precision, intent(in), dimension(m,n,k) :: x
    double precision, intent(in), dimension(k) :: w
    double precision, intent(inout), dimension(m,n) :: meanx
    double precision, intent(inout), dimension(m,m,n) :: covx

    double precision, dimension(m,n,k) :: x2

    x2 = x
    do i = 1, k
        meanx = meanx + x2(:,:,i)*w(i)
    end do
    do i = 1, k
        x2(:,:,i) = sqrt(w(i))*(x2(:,:,i) - meanx)
    end do

    do t = 1, n
        call dgemm('n','t',m,m,k,1.0d0,x2(:,t,:),m,x2(:,t,:),m,0.0d0,covx(:,:,t),m)
    end do

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
    if(var==1) then
    do i = 1, m
        do t = 1, n
            varx(t,i) = sum(w*x(t,i,:)**2)-meanx(t,i)**2
        end do
    end do
end if
end subroutine varmeanw
