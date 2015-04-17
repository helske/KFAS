! simulation filter
subroutine simfilter(ymiss,timevar, yt, zt, ht, tt, rtv, qt, a1, p1, &
p1inf, nnd,nsim, epsplus, etaplus, aplus1, p, n, m, r, info,rankp,&
tol,nd,ndl,sim,c,simwhat,simdim,antithetics)

    implicit none

    integer, intent(in) :: p, m, r, n, nsim,nnd,ndl,simdim,simwhat,antithetics
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    integer, intent(in), dimension(ndl) :: nd
    integer, intent(inout) :: info,rankp
    integer ::  t, i, d, j,k
    double precision, intent(in) :: tol
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(p,p,(n-1)*timevar(2)+1) :: ht
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rtv
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, intent(in), dimension(nsim) :: c
    double precision, intent(inout), dimension(simdim,n,3 * nsim * antithetics + nsim) :: sim
    double precision, intent(inout), dimension(p,n,nsim) :: epsplus
    double precision, intent(inout), dimension(r,n,nsim) :: etaplus
    double precision, intent(inout), dimension(m,nsim) :: aplus1

    double precision, dimension(n,p) :: yplus
    double precision, dimension(m,n+1) :: aplus
    double precision, dimension(p,n) :: ft,finf
    double precision, dimension(m,p,n) :: kt,kinf
    double precision, dimension(r,r,(n-1)*timevar(5)+1) :: cholqt
    double precision, dimension(m,m) :: cholp1
    double precision, dimension(r,r) :: rcholhelp
    double precision, dimension(m,n,4) :: alphatmp

    double precision, dimension(m,n+1) :: at
    double precision, dimension(m,n+1) :: atplus
    double precision, dimension(m,m,n+1) :: pt,pinf
    double precision, dimension(p,n) :: vt
    double precision :: lik
    double precision, dimension(1,p) :: theta
    double precision, dimension(p,p,1) :: thetavar

    double precision, external :: ddot

    external kfilter, filtersimfast, ldl, dtrmv, dgemv


    at=0.0d0
    call kfilter(yt, ymiss, timevar, zt, ht,tt, rtv, qt, a1, p1, p1inf, &
    p,n,m,r,d,j,  at, pt, vt, ft,kt, pinf, finf, kinf, lik, tol,rankp,theta,thetavar,0)

    do t = 1, (n-1)*timevar(5)+1
        if(r.EQ.1) then
            cholqt(1,1,t)=sqrt(qt(1,1,t))
        else
            rcholhelp = qt(:,:,t)
            call ldl(rcholhelp,r,tol,info)
            if(info .NE. 0) then
                info=-2
                return
            end if
            do i=1,r
                cholqt(i,i,t)=sqrt(rcholhelp(i,i))
            end do
            do i=1,r-1
                cholqt((i+1):r,i,t) = rcholhelp((i+1):r,i)*cholqt(i,i,t)
            end do
        end if
    end do
  
  
    if(nnd.GT.0) then
        if(m.EQ.1) then
            cholp1(1,1)=sqrt(p1(1,1))
        else
            cholp1 = p1
            call ldl(cholp1,m,tol,info)
            if(info .NE. 0) then
                info=-3
                return
            end if
            do i=1,m
                cholp1(i,i)=sqrt(cholp1(i,i))
            end do
            do i=1,m-1
                cholp1((i+1):m,i) = cholp1((i+1):m,i)*cholp1(i,i)
            end do
        end if
    end if
  
    do i = 1, nsim
        aplus=0.0d0
        if(ndl.GT.0) then
            aplus(nd,1) = a1(nd)
        end if
        if(nnd.GT.0) then
            call dtrmv('l','n','n',m,cholp1,m,aplus1(:,i),1)
            aplus(:,1) = aplus(:,1)+aplus1(:,i)
        end if

        do t = 1, n
            do k = 1, p
                if(ymiss(t,k).EQ.0) then
                    yplus(t,k) = epsplus(k,t,i)*sqrt(ht(k,k,(t-1)*timevar(2)+1)) + &
                    ddot(m,zt(k,:,(t-1)*timevar(1)+1),1,aplus(:,t),1)
                end if
            end do
            call dtrmv('l','n','n',r,cholqt(:,:,(t-1)*timevar(5)+1),r,etaplus(:,t,i),1)
            call dgemv('n',m,m,1.0d0,tt(:,:,(t-1)*timevar(3)+1),m,aplus(:,t),1,0.0d0,aplus(:,t+1),1)
            call dgemv('n',m,r,1.0d0,rtv(:,:,(t-1)*timevar(4)+1),m,etaplus(:,t,i),1,1.0d0,aplus(:,t+1),1)
        end do


        atplus=0.0d0
        call filtersimfast(yplus, ymiss, timevar, zt,tt, a1, ft,kt,&
        finf, kinf, d, j, p, m, n,atplus)
        if(simwhat.EQ.4) then
            do t = 1, n
                sim(:,t,i) = at(:,t) - atplus(:,t) + aplus(:,t)
                if(antithetics .EQ. 1) then
                    sim(:,t,i+nsim) = at(:,t) + atplus(:,t) - aplus(:,t)
                    sim(:,t,i+2*nsim) = at(:,t)+ c(i)*(sim(:,t,i)-at(:,t))
                    sim(:,t,i+3*nsim) = at(:,t)+ c(i)*(sim(:,t,i+nsim)-at(:,t))
                end if
            end do
        else
            do t = 1, n
                alphatmp(:,t,1) = at(:,t) - atplus(:,t) + aplus(:,t)
                call dgemv('n',p,m,1.0d0,zt(:,:,(t-1)*timevar(1)+1),p,&
                alphatmp(:,t,1),1,0.0d0,sim(:,t,i),1)
                if(antithetics .EQ. 1) then
                    alphatmp(:,t,2) = at(:,t) + atplus(:,t) - aplus(:,t)
                    alphatmp(:,t,3) = at(:,t)+ c(i)*(alphatmp(:,t,1)-at(:,t))
                    alphatmp(:,t,4) = at(:,t)+ c(i)*(alphatmp(:,t,2)-at(:,t))
                    call dgemv('n',p,m,1.0d0,zt(:,:,(t-1)*timevar(1)+1),p,&
                    alphatmp(:,t,2),1,0.0d0,sim(:,t,i+nsim),1)
                    call dgemv('n',p,m,1.0d0,zt(:,:,(t-1)*timevar(1)+1),p,&
                    alphatmp(:,t,3),1,0.0d0,sim(:,t,i+2*nsim),1)
                    call dgemv('n',p,m,1.0d0,zt(:,:,(t-1)*timevar(1)+1),p,&
                    alphatmp(:,t,4),1,0.0d0,sim(:,t,i+3*nsim),1)
                end if
            end do
        end if
    end do
end subroutine simfilter

