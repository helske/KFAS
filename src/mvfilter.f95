! function for computing the multivariate one step ahead prediction errors and their covariances
subroutine mvfilter(tvz, zt, p, m, n, d, at, pt, pinf, vt, ft,finf)


    implicit none

    integer, intent(in) ::  p, m, n,d,tvz
    integer ::  t
    double precision, intent(in), dimension(p,m,(n-1)*tvz+1) :: zt
    double precision, intent(in), dimension(n,m) :: at
    double precision, intent(inout), dimension(n,p) :: vt
    double precision, intent(inout), dimension(p,p,n) :: ft
    double precision, intent(inout), dimension(p,p,d) :: finf
    double precision, intent(in), dimension(m,m,n) :: pt
    double precision, intent(in), dimension(m,m,d) :: pinf
    double precision, dimension(p,m) :: pm
    external dgemv, dgemm, dsymm

    do t = 1,d
        call dgemv('n',p,m,-1.0d0,zt(:,:,(t-1)*tvz+1),p,at(t,:),1,1.0d0,vt(t,:),1)
        call dsymm('r','u',p,m,1.0d0,pt(:,:,t),m,zt(:,:,(t-1)*tvz+1),p,0.0d0,pm,p)
        call dgemm('n','t',p,p,m,1.0d0,pm,p,zt(:,:,(t-1)*tvz+1),p,1.0d0,ft(:,:,t),p)

        call dsymm('r','u',p,m,1.0d0,pinf(:,:,t),m,zt(:,:,(t-1)*tvz+1),p,0.0d0,pm,p)
        call dgemm('n','t',p,p,m,1.0d0,pm,p,zt(:,:,(t-1)*tvz+1),p,0.0d0,finf(:,:,t),p)
    end do

    do t = d+1,n
        call dgemv('n',p,m,-1.0d0,zt(:,:,(t-1)*tvz+1),p,at(t,:),1,1.0d0,vt(t,:),1)
        call dsymm('r','u',p,m,1.0d0,pt(:,:,t),m,zt(:,:,(t-1)*tvz+1),p,0.0d0,pm,p)
        call dgemm('n','t',p,p,m,1.0d0,pm,p,zt(:,:,(t-1)*tvz+1),p,1.0d0,ft(:,:,t),p)
    end do


end subroutine mvfilter
