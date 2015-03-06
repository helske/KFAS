! compute the additional term for marginal likelihood
subroutine marginalxx(p1inf,zt,tt,m,p,n,k,timevar,lik,info)

    integer, intent(inout) :: info
    integer, intent(in) :: m,p,n,k
    integer, intent(in), dimension(5) :: timevar
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,m) ::  p1inf
    double precision, intent(inout) :: lik
    integer ::j,i
    double precision, dimension(m,k) :: a,a2
    double precision, dimension(k,k) :: s
    double precision, dimension(p,k) :: v
    external dgemm, dpotrf
    
    a=0.0d0
    j=1
    do i=1, m
        if(sum(p1inf(:,i))>0.0d0) then
            a(i,j) = 1.0d0
            j = j+1
        end if
    end do
    s=0.0d0
    do i = 1, n
        call dgemm('n','n',p,k,m,1.0d0,zt(:,:,(i-1)*timevar(1)+1),p,a,m,0.0d0,v,p)
        call dgemm('n','n',m,k,m,1.0d0,tt(:,:,(i-1)*timevar(3)+1),m,a,m,0.0d0,a2,m)
        a = a2
        call dsyrk('u','t',k,p,1.0d0,v,p,1.0d0,s,k)    
    end do
    call dpotrf('u', k, s, k, info)
    if(info==0) then    
        do i=1, k
            lik = lik + log(s(i,i))
        end do
    else
        info = -1
    end if 
        
end subroutine marginalxx
