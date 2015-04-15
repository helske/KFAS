  ! Subroutine for Kalman filtering of linear gaussian state space model

subroutine kfilter(yt, ymiss, timevar, zt, ht,tt, rt, qt, a1, p1, p1inf, p,n,m,r,d,j,&
at, pt, vt, ft,kt, pinf, finf, kinf, lik, tol,rankp,theta,thetavar,filtersignal)

    implicit none

    integer, intent(in) :: p, m, r, n,filtersignal
    integer, intent(inout) :: d, j, rankp
    integer :: t, i,tv
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(p,p,(n-1)*timevar(2)+1) :: ht
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rt
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) :: p1,p1inf
    double precision, intent(in) :: tol
    double precision, intent(inout), dimension(m,n+1) :: at
    double precision, intent(inout), dimension(m,m,n+1) :: pt,pinf
    double precision, intent(inout), dimension(p,n) :: vt,ft,finf
    double precision, intent(inout), dimension(m,p,n) :: kt,kinf
    double precision, intent(inout) :: lik
    double precision, intent(inout), dimension(p,p,n) :: thetavar
    double precision, intent(inout), dimension(n,p) :: theta
    double precision, dimension(m) :: arec
    double precision, dimension(m,m) ::pirec,mm,prec
    double precision, dimension(m,r) :: mr
    double precision, dimension(p,m) :: pm
    double precision :: c,meps
    double precision, external :: ddot
    double precision :: finv
    double precision, dimension(m,m,(n-1)*max(timevar(4),timevar(5))+1) :: rqr
    external dgemm, dsymm, dgemv, dsymv, daxpy, dsyr, dsyr2

    meps = epsilon(meps)
    c = 0.5d0*log(8.0d0*atan(1.0d0))

    lik = 0.0d0
    tv= max(timevar(4),timevar(5))
    do t=1, (n-1)*tv+1
        call dsymm('r','l',m,r,1.0d0,qt(:,:,(t-1)*timevar(5)+1),r,rt(:,:,(t-1)*timevar(4)+1),m,0.0d0,mr,m)
        call dgemm('n','t',m,m,r,1.0d0,mr,m,rt(:,:,(t-1)*timevar(4)+1),m,0.0d0,rqr(:,:,t),m)
    end do

    j=0
    d=0
    pinf(:,:,1)=p1inf
    pt(:,:,1) = p1
    prec = p1
    pirec = p1inf
    at(:,1) = a1
    arec = a1
    ! diffuse initialization
    if(rankp > 0) then
        diffuse: do while(d < n .AND. rankp > 0)
            d = d+1
            do j=1, p
                call dsymv('u',m,1.0d0,prec,m,zt(j,:,(d-1)*timevar(1)+1),1,0.0d0,kt(:,j,d),1)
                ft(j,d) = ddot(m,zt(j,:,(d-1)*timevar(1)+1),1,kt(:,j,d),1) + ht(j,j,(d-1)*timevar(2)+1)
                if(ymiss(d,j) .EQ. 0) then
                    call dsymv('u',m,1.0d0,pirec,m,zt(j,:,(d-1)*timevar(1)+1),1,0.0d0,kinf(:,j,d),1) ! kinf_t,i = pinf_t,i*t(z_t,i)
                    finf(j,d) = ddot(m,zt(j,:,(d-1)*timevar(1)+1),1,kinf(:,j,d),1)

                    vt(j,d) = yt(d,j) - ddot(m,zt(j,:,(d-1)*timevar(1)+1),1,arec,1)
                    if (finf(j,d) > tol*maxval(zt(j,:,(d-1)*timevar(1)+1)**2)) then
                        finv = 1.0d0/finf(j,d)
                        arec = arec + kinf(:,j,d)*finv*vt(j,d)
                        call dsyr('u',m,ft(j,d)*finv**2,kinf(:,j,d),1,prec,m) !prec = prec + kinf*kinf'*ft/finf^2
                        call dsyr2('u',m,-finv,kt(:,j,d),1,kinf(:,j,d),1,prec,m) !prec = prec -(kt*kinf'+kinf*kt')/finf
                        call dsyr('u',m,-finv,kinf(:,j,d),1,pirec,m) !pirec = pirec -kinf*kinf'/finf

                        lik = lik - 0.5d0*log(finf(j,d))
                        rankp = rankp -1
                    else
                        finf(j,d) = 0.0d0
                        if(ft(j,d)> tol*maxval(zt(j,:,(d-1)*timevar(1)+1)**2)) then
                            finv = 1.0d0/ft(j,d)
                            call daxpy(m,vt(j,d)*finv,kt(:,j,d),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)/ft(i,t)
                            call dsyr('u',m,-finv,kt(:,j,d),1,prec,m) !prec = prec -kt*kt'/ft
                            lik = lik - c - 0.5d0*(log(ft(j,d)) + vt(j,d)**2*finv)
                        end if
                    end if
                    if (ft(j,d) <= tol*maxval(zt(j,:,(d-1)*timevar(1)+1)**2)) then
                        ft(j,d)=0.0d0
                    end if
                    if(rankp .EQ. 0) then
                        exit diffuse
                    end if
                end if
            end do

            call dgemv('n',m,m,1.0d0,tt(:,:,(d-1)*timevar(3)+1),m,arec,1,0.0d0,at(:,d+1),1)
            arec = at(:,d+1)
            call dsymm('r','u',m,m,1.0d0,prec,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
            call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,pt(:,:,d+1),m)

            pt(:,:,d+1) = pt(:,:,d+1) + rqr(:,:,(d-1)*tv+1)
            prec = pt(:,:,d+1)
            call dsymm('r','u',m,m,1.0d0,pirec,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
            call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,pinf(:,:,d+1),m)

            do i = 1, m
                if(pinf(i,i,d+1) < meps) then
                    pinf(i,:,d+1) = 0.0d0
                    pinf(:,i,d+1) = 0.0d0
                end if
            end do
            pirec = pinf(:,:,d+1)

        end do diffuse


        if(rankp .EQ. 0) then
            !non-diffuse filtering begins
            do i = j+1, p
                call dsymv('u',m,1.0d0,prec,m,zt(i,:,(d-1)*timevar(1)+1),1,0.0d0,kt(:,i,d),1)
                ft(i,d) = ddot(m,zt(i,:,(d-1)*timevar(1)+1),1,kt(:,i,d),1) + ht(i,i,(d-1)*timevar(2)+1)
                if(ymiss(d,i).EQ.0) then
                    vt(i,d) = yt(d,i) - ddot(m,zt(i,:,(d-1)*timevar(1)+1),1,arec,1) !vt
                    if (ft(i,d)> tol*maxval(zt(i,:,(d-1)*timevar(1)+1)**2)) then
                        arec = arec + kt(:,i,d)*vt(i,d)/ft(i,d)
                        call dsyr('u',m,-1.0d0/ft(i,d),kt(:,i,d),1,prec,m) !p_rec = p_rec - kt*kt'*ft(i,t)
                        lik = lik - 0.5d0*(log(ft(i,d)) + vt(i,d)**2/ft(i,d))-c
                    else
                        ft(i,d)=0.0d0
                    end if
                end if
            end do

            call dgemv('n',m,m,1.0d0,tt(:,:,(d-1)*timevar(3)+1),m,arec,1,0.0d0,at(:,d+1),1)
            call dsymm('r','u',m,m,1.0d0,prec,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
            call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,pt(:,:,d+1),m)

            pt(:,:,d+1) = pt(:,:,d+1) + rqr(:,:,(d-1)*tv+1)
            prec = pt(:,:,d+1)
            arec = at(:,d+1)
        end if
    end if

    !Non-diffuse filtering continues from t=d+1, i=1


    if(d .EQ. n .AND. j .EQ. p+1) then
        j = p
    end if

    do t = d+1, n
        do i = 1, p
            call dsymv('u',m,1.0d0,prec,m,zt(i,:,(t-1)*timevar(1)+1),1,0.0d0,kt(:,i,t),1)

            ft(i,t) = ddot(m,zt(i,:,(t-1)*timevar(1)+1),1,kt(:,i,t),1) + ht(i,i,(t-1)*timevar(2)+1)
            if(ymiss(t,i).EQ.0) then
                vt(i,t) = yt(t,i) - ddot(m,zt(i,:,(t-1)*timevar(1)+1),1,arec,1)
                if (ft(i,t)> tol*maxval(zt(i,:,(t-1)*timevar(1)+1)**2)) then
                    arec = arec + kt(:,i,t)/ft(i,t)*vt(i,t)
                    call dsyr('u',m,-1.0d0/ft(i,t),kt(:,i,t),1,prec,m) !p_rec = p_rec - kt*kt'*ft(i,i,t)
                    lik = lik - 0.5d0*(log(ft(i,t)) + vt(i,t)**2/ft(i,t))-c
                else
                    ft(i,t)=0.0d0
                end if
            end if
        end do

        call dgemv('n',m,m,1.0d0,tt(:,:,(t-1)*timevar(3)+1),m,arec,1,0.0d0,at(:,t+1),1)

        call dsymm('r','u',m,m,1.0d0,prec,m,tt(:,:,(t-1)*timevar(3)+1),m,0.0d0,mm,m)
        call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(:,:,(t-1)*timevar(3)+1),m,0.0d0,pt(:,:,t+1),m)

        pt(:,:,t+1) = pt(:,:,t+1) + rqr(:,:,(t-1)*tv+1)
        prec = pt(:,:,t+1)
        arec = at(:,t+1)

    end do

    if(filtersignal==1) then
        do t = 1, n
            call dgemv('n',p,m,1.0d0,zt(:,:,(t-1)*timevar(1)+1),p,at(:,t),1,0.0d0,theta(t,:),1)
            call dsymm('r','u',p,m,1.0d0,pt(:,:,t),m,zt(:,:,(t-1)*timevar(1)+1),p,0.0d0,pm,p)
            call dgemm('n','t',p,p,m,1.0d0,pm,p,zt(:,:,(t-1)*timevar(1)+1),p,0.0d0,thetavar(:,:,t),p)
        end do
    end if
end subroutine kfilter
