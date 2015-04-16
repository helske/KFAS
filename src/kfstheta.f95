! signal smoothing algorithm for gaussian approximation algorithm
subroutine kfstheta(yt, ymiss, timevar, zt, ht,tt, rtv,qt,rqr, tv, a1, p1, p1inf, &
p, n, m, r,tol,rankp,thetahat,lik)

    implicit none

    integer, intent(in) ::  p, m, n,r,tv
    integer, intent(inout) :: rankp
    integer ::  t, i,d, j
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(p,p,(n-1)*timevar(2)+1) :: ht
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rtv
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, dimension(m,m,(n-1)*tv+1) :: rqr
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, dimension(p,n) :: ft,finf
    double precision, dimension(m,p,n) :: kt,kinf
    double precision, dimension(p,n) :: vt
    double precision, dimension(m) :: at, arec
    double precision, dimension(m,m) :: mm,pinf,pt
    double precision, intent(in) :: tol
    double precision, dimension(m) :: rrec,rrec1,rhelp,help
    double precision, dimension(m,m) :: im,linf,l0
    double precision, dimension(r,n) :: etahat
    double precision :: c, meps, finv
    double precision, external :: ddot
    double precision, intent(inout), dimension(n,p) :: thetahat
    double precision, intent(inout) :: lik
    external dgemm, dsymm, dgemv, dsymv, dsyr, dsyr2, dger

    meps = epsilon(meps)
    lik=0.0d0
    c = 0.5d0*log(8.0d0*atan(1.0d0))

    im = 0.0d0
    do i = 1, m
        im(i,i) = 1.0d0
    end do
    j=0
    d=0
    pt = p1
    pinf = p1inf
    at = a1
    if(rankp .GT. 0) then
        diffuse: do while(d .LT. n .AND. rankp .GT. 0)
            d = d+1
            call diffusefilteronestep(ymiss(d,:),yt(d,:),transpose(zt(:,:,(d-1)*timevar(1)+1)),ht(:,:,(d-1)*timevar(2)+1),&
            tt(:,:,(d-1)*timevar(3)+1),rqr(:,:,(d-1)*tv+1),&
            at,pt,vt(:,d),ft(:,d),kt(:,:,d),pinf,finf(:,d),kinf(:,:,d),rankp,lik,tol,meps,c,p,m,j)
        end do diffuse

        if(rankp .EQ. 0 .AND. j .LT. p) then
            call filteronestep(ymiss(d,:),yt(d,:),transpose(zt(:,:,(d-1)*timevar(1)+1)),ht(:,:,(d-1)*timevar(2)+1),&
            tt(:,:,(d-1)*timevar(3)+1),rqr(:,:,(d-1)*tv+1),&
            at,pt,vt(:,d),ft(:,d),kt(:,:,d),lik,tol,meps,c,p,m,j)

        else
            j = p
        end if

    end if

    !Non-diffuse filtering continues from t=d+1, i=1
    do t = d+1, n
        call filteronestep(ymiss(t,:),yt(t,:),transpose(zt(:,:,(t-1)*timevar(1)+1)),ht(:,:,(t-1)*timevar(2)+1),&
        tt(:,:,(t-1)*timevar(3)+1),rqr(:,:,(t-1)*tv+1),&
        at,pt,vt(:,t),ft(:,t),kt(:,:,t),lik,tol,meps,c,p,m,0)

    end do

    !smoothing begins

    rrec = 0.0d0

    do t = n, d+1, -1
        call dgemv('t',m,r,1.0d0,rtv(:,:,(t-1)*timevar(4)+1),m,rrec,1,0.0d0,help,1)
        call dsymv('l',r,1.0d0,qt(:,:,(t-1)*timevar(5)+1),r,help,1,0.0d0,etahat(:,t),1)
        call dgemv('t',m,m,1.0d0,tt(:,:,(t-1)*timevar(3)+1),m,rrec,1,0.0d0,rhelp,1)
        rrec = rhelp
        do i = p, 1 , -1
            if(ymiss(t,i).EQ.0) then
                if(ft(i,t) .GT. 0.0d0) then
                    finv = 1.0d0/ft(i,t)
                    l0 = im
                    call dger(m,m,-finv,kt(:,i,t),1,zt(i,:,(t-1)*timevar(1)+1),1,l0,m)
                    call dgemv('t',m,m,1.0d0,l0,m,rrec,1,0.0d0,rhelp,1)
                    rrec = rhelp + vt(i,t)*finv*zt(i,:,(t-1)*timevar(1)+1)
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
                if(ft(i,t) .GT.  0.0d0) then
                    finv = 1.0d0/ft(i,t)
                    l0 = im
                    call dger(m,m,-finv,kt(:,i,t),1,zt(i,:,(t-1)*timevar(1)+1),1,l0,m)
                    call dgemv('t',m,m,1.0d0,l0,m,rrec,1,0.0d0,rhelp,1)
                    rrec = rhelp + vt(i,t)*finv*zt(i,:,(t-1)*timevar(1)+1)
                end if
            end if
        end do
        rrec1 = 0.0d0
        do i = j, 1, -1
            if(ymiss(t,i).EQ.0) then
                if(finf(i,t).GT. 0.0d0) then
                    finv = 1.0d0/finf(i,t)

                    linf = im
                    call dger(m,m,-finv,kinf(:,i,t),1,zt(i,:,(t-1)*timevar(1)+1),1,linf,m)

                    rhelp = kinf(:,i,t)*ft(i,t)*finv - kt(:,i,t)
                    l0=0.0d0
                    call dger(m,m,finv,rhelp,1,zt(i,:,(t-1)*timevar(1)+1),1,l0,m)

                    call dgemv('t',m,m,1.0d0,linf,m,rrec1,1,0.0d0,rhelp,1)
                    rrec1 = rhelp
                    call dgemv('t',m,m,1.0d0,l0,m,rrec,1,1.0d0,rrec1,1)
                    rrec1 = rrec1 + vt(i,t)*finv*zt(i,:,(t-1)*timevar(1)+1)

                    call dgemv('t',m,m,1.0d0,linf,m,rrec,1,0.0d0,rhelp,1)
                    rrec = rhelp
                else
                    if(ft(i,t).GT. 0.0d0) then
                        finv = 1.0d0/ft(i,t)
                        l0 = im
                        call dger(m,m,-finv,kt(:,i,t),1,zt(i,:,(t-1)*timevar(1)+1),1,l0,m)

                        call dgemv('t',m,m,1.0d0,l0,m,rrec,1,0.0d0,rhelp,1)
                        rrec = rhelp + vt(i,t)*finv*zt(i,:,(t-1)*timevar(1)+1)

                        call dgemv('t',m,m,1.0d0,l0,m,rrec1,1,0.0d0,rhelp,1)
                        rrec1 = rhelp
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
                    if(finf(i,t).GT. 0.0d0) then
                        finv = 1.0d0/finf(i,t)

                        linf = im
                        call dger(m,m,-finv,kinf(:,i,t),1,zt(i,:,(t-1)*timevar(1)+1),1,linf,m)

                        rhelp = kinf(:,i,t)*ft(i,t)*finv - kt(:,i,t)
                        l0=0.0d0
                        call dger(m,m,finv,rhelp,1,zt(i,:,(t-1)*timevar(1)+1),1,l0,m)

                        call dgemv('t',m,m,1.0d0,linf,m,rrec1,1,0.0d0,rhelp,1)
                        rrec1 = rhelp
                        call dgemv('t',m,m,1.0d0,l0,m,rrec,1,1.0d0,rrec1,1)
                        rrec1 = rrec1 + vt(i,t)*finv*zt(i,:,(t-1)*timevar(1)+1)

                        call dgemv('t',m,m,1.0d0,linf,m,rrec,1,0.0d0,rhelp,1)
                        rrec = rhelp
                    else
                        if(ft(i,t).GT.  0.0d0) then
                            finv = 1.0d0/ft(i,t)
                            l0 = im
                            call dger(m,m,-finv,kt(:,i,t),1,zt(i,:,(t-1)*timevar(1)+1),1,l0,m)

                            call dgemv('t',m,m,1.0d0,l0,m,rrec,1,0.0d0,rhelp,1)
                            rrec = rhelp + vt(i,t)*finv*zt(i,:,(t-1)*timevar(1)+1)

                            call dgemv('t',m,m,1.0d0,l0,m,rrec1,1,0.0d0,rhelp,1)
                            rrec1 = rhelp

                        end if

                    end if
                end if
            end do

        end do
    end if


    at = a1

    call dsymv('l',m,1.0d0,p1,m,rrec,1,1.0d0,at,1)
    if(d .GT. 0) then
        call dsymv('l',m,1.0d0,p1inf,m,rrec1,1,1.0d0,at,1)
    end if
    call dgemv('n',p,m,1.0d0,zt(:,:,1),p,at,1,0.0d0,thetahat(1,:),1)

    do t = 2, n
        call dgemv('n',m,m,1.0d0,tt(:,:,(t-2)*timevar(3)+1),m,at,1,0.0d0,help,1)
        at=help
        call dgemv('n',m,r,1.0d0,rtv(:,:,(t-2)*timevar(4)+1),m,etahat(:,t-1),1,1.0d0,at,1)
        call dgemv('n',p,m,1.0d0,zt(:,:,(t-1)*timevar(1)+1),p,at,1,0.0d0,thetahat(t,:),1)

    end do


end subroutine kfstheta
