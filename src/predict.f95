! auxiliary functions for computing Zalpha and ZVZ

subroutine signaltheta(tvz, zt, ahat, vt, p, n, m, theta, thetavar,d,states,m2)

    implicit none

    integer, intent(in) :: p, m, n,tvz,d,m2 !,tvh
    integer, intent(in), dimension(m2) :: states
    integer ::  t
    double precision, intent(in), dimension(p,m,(n-1)*tvz+1) :: zt
    double precision, intent(in), dimension(m,n) :: ahat
    double precision, intent(in), dimension(m,m,n) :: vt
    double precision, intent(inout), dimension(n,p) :: theta
    double precision, intent(inout), dimension(p,p,n) :: thetavar
    double precision, dimension(p,m2) :: pm

    do t = (d+1), n
        call dgemv('n',p,m2,1.0d0,zt(:,states,(t-1)*tvz+1),p,ahat(states,t),1,0.0d0,theta(t,:),1)
        call dsymm('r','u',p,m2,1.0d0,vt(states,states,t),m2,zt(:,states,(t-1)*tvz+1),p,0.0d0,pm,p)
        call dgemm('n','t',p,p,m2,1.0d0,pm,p,zt(:,states,(t-1)*tvz+1),p,0.0d0,thetavar(:,:,t),p)
    end do

end subroutine signaltheta

subroutine zalpha(timevar, zt, alpha,theta,p,m,n,nsim,m2,states)

    implicit none

    integer, intent(in) :: p, m, n, nsim,timevar,m2
    integer, intent(in), dimension(m2) :: states
    integer :: t, i
    double precision, intent(in), dimension(p,m,(n-1)*timevar+1) :: zt
    double precision, intent(in), dimension(n,m,nsim) :: alpha
    double precision, intent(inout), dimension(n,p,nsim) :: theta

   do i=1,nsim
        do t=1,n
            call dgemv('n',p,m2,1.0d0,zt(:,states,(t-1)*timevar+1),p,&
            alpha(t,states,i),1,0.0d0,theta(t,:,i),1)
        end do
    end do

end subroutine zalpha


