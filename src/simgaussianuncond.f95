!simulate gaussian state space model
subroutine simgaussianuncond(timevar, zt, ht, tt, rtv, qt, a1, p1, &
p1inf, nnd,nsim, epsplus, etaplus, aplus1, p, n, m, r, info,rankp,&
tol,nd,ndl,sim,c,simwhat,simdim,antithetics)

    implicit none

    integer, intent(in) :: p, m, r, n, nsim,nnd,ndl,simdim,simwhat,antithetics,rankp
    integer, intent(in), dimension(5) :: timevar
    integer, intent(in), dimension(ndl) :: nd
    integer, intent(inout) :: info
    integer ::  t, i, d, j,k,tv,l,rankp2
    double precision, intent(in) :: tol
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(p,p,(n-1)*timevar(2)+1) :: ht
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rtv
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, intent(in),dimension(nsim) :: c
    double precision, intent(inout), dimension(simdim,n,3 * nsim * antithetics + nsim) :: sim
    double precision, intent(inout), dimension(r,n,nsim) :: etaplus
    double precision, intent(inout), dimension(p,n,nsim) :: epsplus
    double precision, intent(inout), dimension(m,nsim) :: aplus1

    double precision, dimension(n,p) :: yplus
    double precision, dimension(m,n+1) :: aplus
    double precision, dimension(m) :: ahat
    double precision, dimension(m) :: aplushat
    double precision, dimension(p,n) :: ft,finf
    double precision, dimension(m,p,n) :: kt,kinf
    double precision, dimension(r,r,(n-1)*timevar(5)+1) :: cholqt
    double precision, dimension(m,m) :: cholp1
    double precision, dimension(r,r) :: rcholtmp
    double precision, dimension(m) :: rt0,rt1
    double precision, dimension(m,r) :: mr
    double precision, dimension(r,n-1,4) :: etatmp
    double precision, dimension(m,n,4) :: alphatmp
    double precision, external :: ddot

    external dsymm, dgemm, dsymv, ldl, dtrmv, dgemv



    do t = 1, (n-1)*timevar(5)+1
        if(r.EQ.1) then
            cholqt(1,1,t)=sqrt(qt(1,1,t))
        else
            rcholtmp = qt(:,:,t)
            call ldl(rcholtmp,r,tol,info)
            if(info .NE. 0) then
                info = -2
                return
            end if
            do i=1,r
                cholqt(i,i,t)=sqrt(rcholtmp(i,i))
            end do
            do i=1,r-1
                cholqt((i+1):r,i,t) = rcholtmp((i+1):r,i)*cholqt(i,i,t)
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
        aplus = 0.0d0

        aplus(:,1) = a1
        if(nnd.GT.0) then
            call dtrmv('l','n','n',m,cholp1,m,aplus1(:,i),1)
            aplus(:,1) = aplus(:,1)+aplus1(:,i)
        end if

        do t = 1, n
            do k = 1, p
                epsplus(k,t,i) = epsplus(k,t,i)*sqrt(ht(k,k,(t-1)*timevar(2)+1))
            end do
            call dtrmv('l','n','n',r,cholqt(:,:,(t-1)*timevar(5)+1),r,etaplus(:,t,i),1)
            call dgemv('n',m,m,1.0d0,tt(:,:,(t-1)*timevar(3)+1),m,aplus(:,t),1,0.0d0,aplus(:,t+1),1)
            call dgemv('n',m,r,1.0d0,rtv(:,:,(t-1)*timevar(4)+1),m,etaplus(:,t,i),1,1.0d0,aplus(:,t+1),1)
        end do


        !simwhat = 1: epsilon, 2: eta, 3: both, 4: state, 5: signal, 6: observations
        select case(simwhat)
            case(1)
                sim(:,:,i) = epsplus(:,:,i)
                if(antithetics .EQ. 1) then
                    sim(:,:,i+nsim) =  -epsplus(:,:,i)
                    sim(:,:,i+2*nsim) = c(i)*sim(:,:,i)
                    sim(:,:,i+3*nsim) = c(i)*sim(:,:,i+nsim)
                end if
            case(2)
                sim(:,:,i) = etaplus(:,:,i)
                if(antithetics .EQ. 1) then
                    sim(:,:,i+nsim) = -etaplus(:,:,i)
                    sim(:,:,i+2*nsim) = c(i)*sim(:,:,i)
                    sim(:,:,i+3*nsim) = c(i)*sim(:,:,i+nsim)
                end if
            case(3)
                sim(1:p,:,i) = epsplus(:,:,i)
                sim((p+1):,:,i) = etaplus(:,:,i)
                if(antithetics .EQ. 1) then
                    sim(1:p,:,i+nsim) = -epsplus(:,:,i)
                    sim(1:p,:,i+2*nsim) = c(i)*sim(1:p,:,i)
                    sim(1:p,:,i+3*nsim) = c(i)*sim(1:p,:,i+nsim)
                    sim((p+1):,:,i+nsim) = -etaplus(:,:,i)
                    sim((p+1):,:,i+2*nsim) = c(i)*sim((p+1):,:,i)
                    sim((p+1):,:,i+3*nsim) = c(i)*sim((p+1):,:,i+nsim)
                end if
            case(4)
                sim(:,1,i) = aplus(:,1)
                etatmp(:,:,1) = etaplus(:,1:(n-1),i)
                if(antithetics .EQ. 1) then
                    sim(:,1,i+nsim) = -aplus(:,1)
                    sim(:,1,i+2*nsim) = c(i)*sim(:,1,i)
                    sim(:,1,i+3*nsim) = c(i)*sim(:,1,i+nsim)

                    etatmp(:,:,2) = -etaplus(:,1:(n-1),i)
                    etatmp(:,:,3) = c(i)*etatmp(:,:,1)
                    etatmp(:,:,4) = c(i)*etatmp(:,:,2)
                end if
                do k = 1, 3*antithetics+1
                    do t = 2, n
                        call dgemv('n',m,m,1.0d0,tt(:,:,(t-2)*timevar(3)+1),m,sim(:,t-1,i+(k-1)*nsim),&
                        1,0.0d0,sim(:,t,i+(k-1)*nsim),1)
                        call dgemv('n',m,r,1.0d0,rtv(:,:,(t-2)*timevar(4)+1),m,etatmp(:,t-1,k),1,&
                        1.0d0,sim(:,t,i+(k-1)*nsim),1)
                    end do
                end do
            case(5)

                alphatmp(:,1,1) = aplus(:,1)
                etatmp(:,:,1) =  etaplus(:,1:(n-1),i)
                call dgemv('n',p,m,1.0d0,zt(:,:,1),p,alphatmp(:,1,1),1,0.0d0,sim(:,1,i),1)

                if(antithetics .EQ. 1) then
                    alphatmp(:,1,2) = -aplus(:,1)
                    alphatmp(:,1,3) = c(i)*alphatmp(:,1,1)
                    alphatmp(:,1,4) = c(i)*alphatmp(:,1,2)
                    call dgemv('n',p,m,1.0d0,zt(:,:,1),p,alphatmp(:,1,2),1,0.0d0,sim(:,1,i+nsim),1)
                    call dgemv('n',p,m,1.0d0,zt(:,:,1),p,alphatmp(:,1,3),1,0.0d0,sim(:,1,i+2*nsim),1)
                    call dgemv('n',p,m,1.0d0,zt(:,:,1),p,alphatmp(:,1,4),1,0.0d0,sim(:,1,i+3*nsim),1)

                    etatmp(:,:,2) = -etaplus(:,1:(n-1),i)
                    etatmp(:,:,3) = c(i)*etatmp(:,:,1)
                    etatmp(:,:,4) = c(i)*etatmp(:,:,2)

                end if
                do k = 1, 3*antithetics+1
                    do t = 2, n
                        call dgemv('n',m,m,1.0d0,tt(:,:,(t-2)*timevar(3)+1),m,alphatmp(:,t-1,k),&
                        1,0.0d0,alphatmp(:,t,k),1)
                        call dgemv('n',m,r,1.0d0,rtv(:,:,(t-2)*timevar(4)+1),m,etatmp(:,t-1,k),1,&
                        1.0d0,alphatmp(:,t,k),1)
                        call dgemv('n',p,m,1.0d0,zt(:,:,(t-1)*timevar(1)+1),p,&
                        alphatmp(:,t,k),1,0.0d0,sim(:,t,i+(k-1)*nsim),1)
                    end do
                end do
            case(6)


                alphatmp(:,1,1) = aplus(:,1)
                etatmp(:,:,1) = etaplus(:,1:(n-1),i)

                if(antithetics .EQ. 1) then
                    alphatmp(:,1,2) = -aplus(:,1)
                    alphatmp(:,1,3) = c(i)*alphatmp(:,1,1)
                    alphatmp(:,1,4) = c(i)*alphatmp(:,1,2)

                    etatmp(:,:,2) = -etaplus(:,1:(n-1),i)
                    etatmp(:,:,3) = c(i)*etatmp(:,:,1)
                    etatmp(:,:,4) = c(i)*etatmp(:,:,2)

                end if
                do k = 1, 3*antithetics+1
                    do l = 1, p
                          sim(l,1,i+(k-1)*nsim) = ddot(m,zt(l,:,1),1,alphatmp(:,1,k),1)

                    end do
                    do t = 2, n
                        call dgemv('n',m,m,1.0d0,tt(:,:,(t-2)*timevar(3)+1),m,alphatmp(:,t-1,k),&
                        1,0.0d0,alphatmp(:,t,k),1)
                        call dgemv('n',m,r,1.0d0,rtv(:,:,(t-2)*timevar(4)+1),m,etatmp(:,t-1,k),1,&
                        1.0d0,alphatmp(:,t,k),1)
                        do l = 1, p

                                sim(l,t,i+(k-1)*nsim) = ddot(m,zt(l,:,(t-1)*timevar(1)+1),1,alphatmp(:,t,k),1)

                        end do
                    end do
                end do
        end select
    end do

end subroutine simgaussianuncond

