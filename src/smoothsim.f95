! disturbance smoothing algorithm for simulation

subroutine smoothsim(yt, ymiss, timevar, zt, ht,tt, rtv,qt,rqr, a1, p1, p1inf, &
d, j, p, m, n, r,tol,rankp,ft,finf,kt,kinf,epshat,etahat,rt0,rt1,needeps)

    implicit none

    logical, intent(in) :: needeps
    integer, intent(in) ::  p, m, n,r
    integer, intent(inout) :: d, j,rankp
    integer ::  t, i,tv
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(p,p,(n-1)*timevar(2)+1) :: ht
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rtv
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, dimension(m,m,(n-1)*max(timevar(4),timevar(5))+1) :: rqr
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, intent(inout), dimension(p,n) :: ft,finf
    double precision, intent(inout), dimension(m,p,n) :: kt,kinf
    double precision, dimension(p,n) :: vt
    double precision, dimension(m) :: at, arec
    double precision, dimension(m,m) :: prec,pirec,mm
    double precision, intent(in) :: tol
    double precision, dimension(m) :: rrec,rrec1,rhelp,help
    double precision, dimension(m,m) :: im,linf,l0,lt
    double precision, intent(inout), dimension(r,n) :: etahat
    double precision, intent(inout), dimension(p,n) :: epshat
    double precision, intent(inout), dimension(m) :: rt0,rt1
    double precision :: meps
    double precision, external :: ddot

    external dgemm, dsymm, dgemv, dsymv, daxpy, dsyr, dsyr2, dger

    meps = epsilon(meps)**0.75d0
    tv = max(timevar(4),timevar(5))

    im = 0.0d0
    do i = 1, m
        im(i,i) = 1.0d0
    end do
    j=0
    d=0
    pirec = p1inf
    prec = p1
    arec = a1
    if(maxval(pirec) .GT.  0.0d0) then


        diffuse: do while(d .LT. n)
            d = d+1
            do j=1, p
                if(ymiss(d,j).EQ.0) then

                    vt(j,d) = yt(d,j) - ddot(m,zt(j,:,(d-1)*timevar(1)+1),1,arec,1) !arec
                    call dsymv('u',m,1.0d0,prec,m,zt(j,:,(d-1)*timevar(1)+1),1,0.0d0,kt(:,j,d),1) ! kt_t,i = pt_t,i*t(z_t,i)
                    ft(j,d) = ddot(m,zt(j,:,(d-1)*timevar(1)+1),1,kt(:,j,d),1)  + ht(j,j,(d-1)*timevar(2)+1)

                    call dsymv('u',m,1.0d0,pirec,m,zt(j,:,(d-1)*timevar(1)+1),1,0.0d0,kinf(:,j,d),1) ! kinf_t,i = pinf_t,i*t(z_t,i)
                    finf(j,d) = ddot(m,zt(j,:,(d-1)*timevar(1)+1),1,kinf(:,j,d),1)! finf


                    if (finf(j,d) > meps*maxval(zt(j,:,(d-1)*timevar(1)+1))**2) then
                        call daxpy(m,vt(j,d)/finf(j,d),kinf(:,j,d),1,arec,1) !a_rec = a_rec + kinf(:,i,t)*vt(:,t)/finf(j,d)
                        call dsyr('u',m,ft(j,d)/(finf(j,d)**2),kinf(:,j,d),1,prec,m) !prec = prec +  kinf*kinf'*ft/finf^2
                        call dsyr2('u',m,-1.0d0/finf(j,d),kt(:,j,d),1,kinf(:,j,d),1,prec,m) !prec = prec -(kt*kinf'+kinf*kt')/finf

                        call dsyr('u',m,(-1.0d0/finf(j,d)),kinf(:,j,d),1,pirec,m) !pirec = pirec -kinf*kinf'/finf

                        rankp = rankp -1

                    else
                        finf(j,d) = 0.0d0
                        if(ft(j,d)> meps*maxval(zt(j,:,(d-1)*timevar(1)+1))**2) then
                            call daxpy(m,vt(j,d)/ft(j,d),kt(:,j,d),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)/ft(i,t)
                            call dsyr('u',m,(-1.0d0)/ft(j,d),kt(:,j,d),1,prec,m) !prec = prec -kt*kt'/ft
                        else
                            ft(j,d)=0.0d0
                        end if
                    end if

                    if(rankp .EQ. 0) then
                        exit diffuse
                    end if
                end if
            end do

            call dgemv('n',m,m,1.0d0,tt(:,:,(d-1)*timevar(3)+1),m,arec,1,0.0d0,at,1)  !at(:,t+1) = matmul(tt,a_rec)
            arec = at
            call dsymm('r','u',m,m,1.0d0,prec,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
            call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,prec,m)

            prec = prec + rqr(:,:,(d-1)*tv+1)

            call dsymm('r','u',m,m,1.0d0,pirec,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
            call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,pirec,m)

            do i = 1, m
                if(pirec(i,i) .LT. meps) then
                    pirec(i,:) = 0.0d0
                    pirec(:,i) = 0.0d0
                end if
            end do

        end do diffuse

        !non-diffuse filtering begins
        if(rankp .EQ. 0) then
            do i = j+1, p
                if(ymiss(d,i).EQ.0) then
                    vt(i,d) = yt(d,i) - ddot(m,zt(i,:,(d-1)*timevar(1)+1),1,arec,1) !vt
                    call dsymv('u',m,1.0d0,prec,m,zt(i,:,(d-1)*timevar(1)+1),1,0.0d0,kt(:,i,d),1) ! p symmetric!
                    ft(i,d) = ddot(m,zt(i,:,(d-1)*timevar(1)+1),1,kt(:,i,d),1)  +  ht(i,i,(d-1)*timevar(2)+1)
                    if (ft(i,d)> meps*maxval(zt(i,:,(d-1)*timevar(1)+1))**2) then
                        call daxpy(m,vt(i,d)/ft(i,d),kt(:,i,d),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)
                        call dsyr('u',m,-1.0d0/ft(i,d),kt(:,i,d),1,prec,m) !p_rec = p_rec - kt*kt'*ft(i,t)
                    else
                        ft(i,d)=0.0d0
                    end if
                end if
            end do

            call dgemv('n',m,m,1.0d0,tt(:,:,(d-1)*timevar(3)+1),m,arec,1,0.0d0,at,1)  !at(:,t+1) = matmul(tt,a_rec)

            call dsymm('r','u',m,m,1.0d0,prec,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
            call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,prec,m)

            prec = prec + rqr(:,:,(d-1)*tv+1)
            arec = at
        end if
    end if

    !Non-diffuse filtering continues from t=d+1, i=1



    if(d .EQ. n .AND. j .EQ. p+1) then
        j = p
    end if

    do t = d+1, n
        do i = 1, p
            if(ymiss(t,i).EQ.0) then
                vt(i,t) = yt(t,i) - ddot(m,zt(i,:,(t-1)*timevar(1)+1),1,arec,1)
                call dsymv('u',m,1.0d0,prec,m,zt(i,:,(t-1)*timevar(1)+1),1,0.0d0,kt(:,i,t),1)
                ft(i,t) = ddot(m,zt(i,:,(t-1)*timevar(1)+1),1,kt(:,i,t),1) +  ht(i,i,(t-1)*timevar(2)+1)
                if (ft(i,t)> meps*maxval(zt(i,:,(t-1)*timevar(1)+1))**2) then
                    call daxpy(m,vt(i,t)/ft(i,t),kt(:,i,t),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)
                    call dsyr('u',m,-1.0d0/ft(i,t),kt(:,i,t),1,prec,m) !p_rec = p_rec - kt*kt'*ft(i,i,t)
                else
                    ft(i,t)=0.0d0
                end if

            end if
        end do

        call dgemv('n',m,m,1.0d0,tt(:,:,(t-1)*timevar(3)+1),m,arec,1,0.0d0,at,1)  !at(:,t+1) = matmul(tt,a_rec)
        call dsymm('r','u',m,m,1.0d0,prec,m,tt(:,:,(t-1)*timevar(3)+1),m,0.0d0,mm,m)
        call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(:,:,(t-1)*timevar(3)+1),m,0.0d0,prec,m)

        prec = prec + rqr(:,:,(t-1)*tv+1)

        arec = at
    end do

    !smoothing begins

    rrec = 0.0d0

    do t = n, d+1, -1
        call dgemv('t',m,r,1.0d0,rtv(:,:,(t-1)*timevar(4)+1),m,rrec,1,0.0d0,help,1)
        call dsymv('l',r,1.0d0,qt(:,:,(t-1)*timevar(5)+1),r,help,1,0.0d0,etahat(:,t),1)
        call dgemv('t',m,m,1.0d0,tt(:,:,(t-1)*timevar(3)+1),m,rrec,1,0.0d0,rhelp,1) !r_t,p=t_t-1'*r_t+1
        rrec = rhelp
        do i = p, 1 , -1
            if(ymiss(t,i)==0) then
                if(ft(i,t) > meps) then
                    if(needeps) then
                        epshat(i,t) = ht(i,i,(t-1)*timevar(2)+1)*(vt(i,t)-ddot(m,kt(:,i,t),1,rrec,1))/ft(i,t)
                    end if
                    lt = im
                    call dger(m,m,-1.0d0/ft(i,t),kt(:,i,t),1,zt(i,:,(t-1)*timevar(1)+1),1,lt,m) !l = I -kz
                    call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1)
                    rrec = rhelp + vt(i,t)/ft(i,t)*zt(i,:,(t-1)*timevar(1)+1)

                end if
            end if
        end do
    end do

    if(d.GT.0) then
        t=d
        call dgemv('t',m,r,1.0d0,rtv(:,:,(t-1)*timevar(4)+1),m,rrec,1,0.0d0,help,1)
        call dsymv('l',r,1.0d0,qt(:,:,(t-1)*timevar(5)+1),r,help,1,0.0d0,etahat(:,t),1)
        call dgemv('t',m,m,1.0d0,tt(:,:,(t-1)*timevar(3)+1),m,rrec,1,0.0d0,rhelp,1)
        rrec = rhelp
        do i = p, (j+1) , -1
            if(ymiss(t,i).EQ.0) then
                if(ft(i,t)  > meps) then
                    if(needeps) then
                        epshat(i,t) = ht(i,i,(t-1)*timevar(2)+1)/ft(i,t)*(vt(i,t)-ddot(m,kt(:,i,t),1,rrec,1))
                    end if
                    lt = im
                    call dger(m,m,-1.0d0/ft(i,t),kt(:,i,t),1,zt(i,:,(t-1)*timevar(1)+1),1,lt,m) !l = i -kz
                    call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1)
                    rrec=rhelp
                    call daxpy(m,vt(i,t)/ft(i,t),zt(i,:,(t-1)*timevar(1)+1),1,rrec,1)
                end if
            end if
        end do
        rrec1 = 0.0d0
        do i = j, 1, -1
            if(ymiss(t,i).EQ.0) then
                if(finf(i,t)> meps) then
                    if(needeps) then
                        epshat(i,t) = -ht(i,i,(t-1)*timevar(2)+1)*ddot(m,kinf(:,i,t),1,rrec,1)/finf(i,t)
                    end if
                    linf = im
                    call dger(m,m,-1.0d0/finf(i,t),kinf(:,i,t),1,zt(i,:,(t-1)*timevar(1)+1),1,linf,m) !linf
                    rhelp = -kt(:,i,t)
                    call daxpy(m,ft(i,t)/finf(i,t),kinf(:,i,t),1,rhelp,1)
                    l0=0.0d0
                    call dger(m,m,(1.0d0/finf(i,t)),rhelp,1,zt(i,:,(t-1)*timevar(1)+1),1,l0,m) !l0
                    call dgemv('t',m,m,1.0d0,linf,m,rrec1,1,0.0d0,rhelp,1) !rt1
                    rrec1 = rhelp
                    call dgemv('t',m,m,1.0d0,l0,m,rrec,1,1.0d0,rrec1,1)
                    call daxpy(m,(vt(i,t)/finf(i,t)),zt(i,:,(t-1)*timevar(1)+1),1,rrec1,1)
                    call dgemv('t',m,m,1.0d0,linf,m,rrec,1,0.0d0,rhelp,1) !rt0
                    rrec = rhelp
                else
                    if(ft(i,t)> meps) then
                        if(needeps) then
                            epshat(i,t) = ht(i,i,(t-1)*timevar(2)+1)*(vt(i,t)-ddot(m,kt(:,i,t),1,rrec,1))/ft(i,t)
                        end if
                        lt= im
                        call dger(m,m,-1.0d0/ft(i,t),kt(:,i,t),1,zt(i,:,(t-1)*timevar(1)+1),1,lt,m) !lt = I -Kt*Z/Ft
                        call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1)
                        rrec = rhelp
                        call daxpy(m,vt(i,t)/ft(i,t),zt(i,:,(t-1)*timevar(1)+1),1,rrec,1) !r0 = Z'vt/Ft - Lt'r0
                        call dgemv('t',m,m,1.0d0,lt,m,rrec1,1,0.0d0,rhelp,1)
                        rrec1=rhelp
                    end if
                end if
            end if
        end do

        do t=(d-1), 1, -1
            call dgemv('t',m,r,1.0d0,rtv(:,:,(t-1)*timevar(4)+1),m,rrec,1,0.0d0,help,1)
            call dsymv('l',r,1.0d0,qt(:,:,(t-1)*timevar(5)+1),r,help,1,0.0d0,etahat(:,t),1)
            call dgemv('t',m,m,1.0d0,tt(:,:,(t-1)*timevar(3)+1),m,rrec,1,0.0d0,rhelp,1)
            rrec = rhelp
            call dgemv('t',m,m,1.0d0,tt(:,:,(t-1)*timevar(3)+1),m,rrec1,1,0.0d0,rhelp,1)
            rrec1 = rhelp

            do i = p, 1, -1
                if(ymiss(t,i).EQ.0) then
                    if(finf(i,t)> meps) then
                        if(needeps) then
                            epshat(i,t) = -ht(i,i,(t-1)*timevar(2)+1)*ddot(m,kinf(:,i,t),1,rrec,1)/finf(i,t)
                        end if
                        linf = im
                        call dger(m,m,-1.0d0/finf(i,t),kinf(:,i,t),1,zt(i,:,(t-1)*timevar(1)+1),1,linf,m) !linf
                        rhelp = -kt(:,i,t)
                        call daxpy(m,ft(i,t)/finf(i,t),kinf(:,i,t),1,rhelp,1)
                        l0=0.0d0
                        call dger(m,m,(1.0d0/finf(i,t)),rhelp,1,zt(i,:,(t-1)*timevar(1)+1),1,l0,m) !l0
                        call dgemv('t',m,m,1.0d0,linf,m,rrec1,1,0.0d0,rhelp,1) !rt1
                        rrec1 = rhelp
                        call dgemv('t',m,m,1.0d0,l0,m,rrec,1,1.0d0,rrec1,1)
                        call daxpy(m,vt(i,t)/finf(i,t),zt(i,:,(t-1)*timevar(1)+1),1,rrec1,1)
                        call dgemv('t',m,m,1.0d0,linf,m,rrec,1,0.0d0,rhelp,1) !rt0
                        rrec = rhelp
                    else
                        if(ft(i,t)> meps) then
                            if(needeps) then
                                epshat(i,t) = ht(i,i,(t-1)*timevar(2)+1)*(vt(i,t)-ddot(m,kt(:,i,t),1,rrec,1))/ft(i,t)
                            end if
                            lt= im
                            call dger(m,m,-1.0d0/ft(i,t),kt(:,i,t),1,zt(i,:,(t-1)*timevar(1)+1),1,lt,m) !lt = I -Kt*Z/Ft
                            call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1)
                            rrec = rhelp
                            call daxpy(m,vt(i,t)/ft(i,t),zt(i,:,(t-1)*timevar(1)+1),1,rrec,1) !r0 = Z'vt/Ft - Lt'r0
                            call dgemv('t',m,m,1.0d0,lt,m,rrec1,1,0.0d0,rhelp,1)
                            rrec1=rhelp

                        end if

                    end if
                end if
            end do

        end do
    end if

    rt0=rrec
    rt1=rrec1

end subroutine smoothsim
