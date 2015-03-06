! Subroutine for Kalman smoothing of linear gaussian state space model

subroutine gsmoothall(ymiss, timevar, zt, ht,tt, rtv, qt, p, n, m, r, d,j, at, pt, vt, ft, kt, &
rt, rt0, rt1, nt, nt0, nt1, nt2, pinf, kinf,finf,  tol,ahat, vvt,epshat,epshatvar, &
etahat,etahatvar,thetahat,thetahatvar, ldlsignal,zorig, zorigtv,aug,state,dist,signal)

    implicit none

    integer, intent(in) :: d, j, p, r, m, n,aug,state,dist,signal,ldlsignal,zorigtv
    integer :: t, i
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(p,p,(n-1)*timevar(2)+1) :: ht
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rtv
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m,n+1) :: at
    double precision, intent(in), dimension(m,m,n+1) :: pt
    double precision, intent(in), dimension(p,n) ::  vt,ft
    double precision, intent(in), dimension(m,p,n) :: kt
    double precision, intent(in), dimension(m,m,d+1) ::  pinf
    double precision, intent(in),dimension(m,p,d) ::  kinf
    double precision, intent(in), dimension(p,d) ::  finf
    double precision, intent(in) :: tol
    double precision, intent(inout), dimension(m,m,n+1) :: nt !n_1 = n_0, ..., n_201 = n_200
    double precision, intent(inout), dimension(m,n+1) :: rt !same as n, r_1 = r_0 etc.
    double precision, intent(inout), dimension(m,d+1) :: rt0,rt1
    double precision, intent(inout), dimension(m,m,d+1) :: nt0,nt1,nt2
    double precision, intent(inout), dimension(m*state,n*state) :: ahat
    double precision, intent(inout), dimension(m*state,m*state,n*state) :: vvt
    double precision, intent(inout), dimension(p*dist*aug,n*dist*aug) :: epshat
    double precision, intent(inout), dimension(p*dist*aug,n*dist*aug) :: epshatvar
    double precision, intent(inout), dimension(r*dist,n*dist) :: etahat
    double precision, intent(inout), dimension(r*dist,r*dist,n*dist) :: etahatvar
    double precision, intent(inout), dimension(p*signal,n*signal) :: thetahat
    double precision, intent(inout), dimension(p*signal,p*signal,n*signal) :: thetahatvar
    double precision, intent(in), dimension(ldlsignal*p,ldlsignal*m,ldlsignal*((n-1)*zorigtv+1)) :: zorig
    double precision, dimension(m,m) :: linf,l0
    double precision, dimension(m,m) :: nrec,nrec1,nrec2,im,mm,mm2
    double precision, dimension(m) :: rrec,rrec1,rhelp, help
    double precision, dimension(m,r) :: mr, mr2
    double precision, dimension(p,m) :: pm
    double precision, dimension(p,n) ::  ftinv
    double precision, dimension(p,d) ::  finfinv
    double precision, external :: ddot
    external dgemm, dsymm, dgemv, dsymv, dger, dcopy

    if(aug.EQ.1 .AND. dist.EQ.1) then
        do i = 1, p
            do t = 1, n
                epshatvar(i,t) =  ht(i,i,(t-1)*timevar(2)+1)
            end do
        end do
    end if


    ftinv = 1.0d0/ft
    finfinv = 1.0d0/finf

    im = 0.0d0
    do i = 1, m
        im(i,i) = 1.0d0
    end do

    rrec = 0.0d0
    nrec = 0.0d0
    nt(:,:,n+1) = 0.0d0 !t goes from n+1 to 1, not from n to 0 !
    rt(:,n+1) = 0.0d0


    do t = n, d+1, -1 !do until diffuse starts


        do i = p, 1 , -1
            if(ymiss(t,i).EQ.0) then
                if(ft(i,t) .GT. 0.0d0) then
                    if(aug.EQ.1 .AND. dist.EQ.1) then
                        epshat(i,t) = ht(i,i,(t-1)*timevar(2)+1)/ft(i,t)*(vt(i,t)-ddot(m,kt(:,i,t),1,rrec,1))
                        call dgemv('n',m,m,ftinv(i,t)**2, nrec,m,kt(:,i,t),1,0.0d0,rhelp,1)
                        epshatvar(i,t) = ht(i,i,(t-1)*timevar(2)+1)-(ht(i,i,(t-1)*timevar(2)+1)**2)*&
                        (ftinv(i,t)+ddot(m,kt(:,i,t),1,rhelp,1))
                    end if

                    rhelp = -zt(i,:,(t-1)*timevar(1)+1)*ftinv(i,t)
                    l0 = im
                    call dgemm('n','n',m,m,1,1.0d0,kt(:,i,t),&
                    m,rhelp,1,1.0d0,l0,m)


                    call dgemv('t',m,m,1.0d0,l0,m,rrec,1,0.0d0,rhelp,1)
                    rrec = rhelp + vt(i,t)*ftinv(i,t)*zt(i,:,(t-1)*timevar(1)+1)
                    call dgemm('t','n',m,m,m,1.0d0,l0,m,nrec,m,0.0d0,mm,m) !n = l'nl
                    call dgemm('n','n',m,m,m,1.0d0,mm,m,l0,m,0.0d0,nrec,m)
                    call dgemm('t','n',m,m,1,1.0d0,zt(i,:,(t-1)*timevar(1)+1),&
                    1,zt(i,:,(t-1)*timevar(1)+1),1,0.0d0,mm,m)
                    nrec = nrec+mm*ftinv(i,t)

                end if
            end if
        end do

        rt(:,t) =rrec
        nt(:,:,t) = nrec !n_t-1 = n_t,0

        if(t.GT.1) then
            call dgemv('t',m,m,1.0d0,tt(:,:,(t-2)*timevar(3)+1),m,rrec,1,0.0d0,rhelp,1) !r_t,p=t_t-1'*r_t+1
            rrec = rhelp
            call dgemm('t','n',m,m,m,1.0d0,tt(:,:,(t-2)*timevar(3)+1),m,nrec,m,0.0d0,mm,m) !n*t
            call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(:,:,(t-2)*timevar(3)+1),m,0.0d0,nrec,m) !n_t,p = t'nt
        end if

    end do


    if(d.GT.0) then
        t=d
        rt0(:,d+1)=rt(:,d+1)
        nt0(:,:,d+1) =  nt(:,:,d+1)

        do i = p, (j+1) , -1
            if(ymiss(t,i).EQ.0) then
                if(ft(i,t) .GT. 0.0d0) then
                    if(aug .EQ. 1 .AND. dist.EQ.1) then
                        epshat(i,t) = ht(i,i,(t-1)*timevar(2)+1)/ft(i,t)*(vt(i,t)-ddot(m,kt(:,i,t),1,rrec,1))
                        call dgemv('n',m,m,ftinv(i,t)**2,nrec,m,kt(:,i,t),1,0.0d0,rhelp,1)
                        epshatvar(i,t) = ht(i,i,(t-1)*timevar(2)+1)-(ht(i,i,(t-1)*timevar(2)+1)**2)*&
                        (ftinv(i,t)+ddot(m,kt(:,i,t),1,rhelp,1))
                    end if

                    rhelp = -zt(i,:,(t-1)*timevar(1)+1)*ftinv(i,t)
                    l0 = im
                    call dgemm('n','n',m,m,1,1.0d0,kt(:,i,t),&
                    m,rhelp,1,1.0d0,l0,m)


                    call dgemv('t',m,m,1.0d0,l0,m,rrec,1,0.0d0,rhelp,1)
                    rrec = rhelp + vt(i,t)*ftinv(i,t)*zt(i,:,(t-1)*timevar(1)+1)
                    call dgemm('t','n',m,m,m,1.0d0,l0,m,nrec,m,0.0d0,mm,m) !n = l'nl
                    call dgemm('n','n',m,m,m,1.0d0,mm,m,l0,m,0.0d0,nrec,m)
                    call dgemm('t','n',m,m,1,1.0d0,zt(i,:,(t-1)*timevar(1)+1),&
                    1,zt(i,:,(t-1)*timevar(1)+1),1,0.0d0,mm,m)
                    nrec = mm*ftinv(i,t)+nrec
                end if
            end if
        end do

        rrec1 = 0.0d0
        nrec1 = 0.0d0
        nrec2 = 0.0d0

        do i = j, 1, -1
            if(ymiss(t,i).EQ.0) then
                if(finf(i,t).GT.tol) then
                    if(aug .EQ. 1 .AND. dist.EQ.1) then
                        epshat(i,t) = -ht(i,i,(t-1)*timevar(2)+1)*ddot(m,kinf(:,i,t),1,rrec,1)/finf(i,t)
                        call dgemv('n',m,m,1.0d0,nrec,m,kinf(:,i,t),1,0.0d0,rhelp,1)
                        epshatvar(i,t) = ht(i,i,(t-1)*timevar(2)+1)-(ht(i,i,(t-1)*timevar(2)+1)**2)*&
                        ddot(m,kinf(:,i,t),1,rhelp,1)/finf(i,t)**2
                    end if


                    linf = im
                    rhelp = -zt(i,:,(t-1)*timevar(1)+1)*finfinv(i,t)
                    call dger(m,m,1.0d0,kinf(:,i,t),1,rhelp,1,linf,m)

                    rhelp = (kinf(:,i,t)*ft(i,t)*finfinv(i,t)-kt(:,i,t))*finfinv(i,t)
                    l0=0.0d0
                    call dger(m,m,1.0d0,rhelp,1,zt(i,:,(t-1)*timevar(1)+1),1,l0,m) !l0=  (-kt + ft/finf*kinf)*z/finf

                    call dgemv('t',m,m,1.0d0,linf,m,rrec1,1,0.0d0,rhelp,1) !rt1
                    rrec1 = rhelp + zt(i,:,(t-1)*timevar(1)+1)*vt(i,t)*finfinv(i,t)
                    call dgemv('t',m,m,1.0d0,l0,m,rrec,1,1.0d0,rrec1,1)
                    call dgemv('t',m,m,1.0d0,linf,m,rrec,1,0.0d0,rhelp,1) !rt0
                    rrec = rhelp

                    call dgemm('t','n',m,m,m,1.0d0,linf,m,nrec2,m,0.0d0,mm,m) !mm =linf'*nt2
                    call dgemm('n','n',m,m,m,1.0d0,mm,m,linf,m,0.0d0,nrec2,m) !nt2 = linf'*nt2*linf
                    call dger(m,m,-1.0d0*ft(i,t)*finfinv(i,t)**2.0d0,zt(i,:,(t-1)*timevar(1)+1)&
                    ,1,zt(i,:,(t-1)*timevar(1)+1),1,nrec2,m) !nt2 = linf'nt2'linf + z'z*ft/finf^2
                    call dgemm('t','n',m,m,m,1.0d0,l0,m,nrec,m,0.0d0,mm,m) !mm= nt0*l0
                    call dgemm('n','n',m,m,m,1.0d0,mm,m,l0,m,1.0d0,nrec2,m) !nt2 = linf'nt2'linf + z'z*ft/finf^2 + l0'*nt0*l0

                    call dgemm('t','n',m,m,m,1.0d0,linf,m,nrec1,m,0.0d0,mm,m) !mm = linf'*nt1
                    call dgemm('n','n',m,m,m,1.0d0,mm,m,l0,m,0.0d0,mm2,m) !nt2 = nt2 + linf'*nt1*l0
                    nrec2 = nrec2 +mm2 + transpose(mm2)

                    call dgemm('t','n',m,m,m,1.0d0,linf,m,nrec1,m,0.0d0,mm,m) !mm = linf'*nt1
                    call dgemm('n','n',m,m,m,1.0d0,mm,m,linf,m,0.0d0,nrec1,m) !nt1 = mm*linf

                    call dger(m,m,finfinv(i,t),zt(i,:,(t-1)*timevar(1)+1),1,zt(i,:,(t-1)*timevar(1)+1),1,nrec1,m)
                    !nt1 = linf'nt1'linf + z'z/finf
                    call dgemm('t','n',m,m,m,1.0d0,l0,m,nrec,m,0.0d0,mm,m) !mm= nt0*linf
                    call dgemm('n','n',m,m,m,1.0d0,mm,m,linf,m,1.0d0,nrec1,m) !nt1 = l0'*nt0*linf+ linf'nt1*linf + z'z/finf

                    call dgemm('t','n',m,m,m,1.0d0,linf,m,nrec,m,0.0d0,mm,m) !mm= nt0*linf
                    call dgemm('n','n',m,m,m,1.0d0,mm,m,linf,m,0.0d0,nrec,m)

                else
                    if(ft(i,t).GT.0.0d0) then
                        if(aug .EQ. 1 .AND. dist.EQ.1) then
                            epshat(i,t) = ht(i,i,(t-1)*timevar(2)+1)*(vt(i,t)/ft(i,t)-&
                            ddot(m,kt(:,i,t),1,rrec,1)/ft(i,t))
                            call dgemv('n',m,m,1.0d0,nrec,m,kt(:,i,t),1,0.0d0,rhelp,1)
                            epshatvar(i,t) = ht(i,i,(t-1)*timevar(2)+1)-(ht(i,i,(t-1)*timevar(2)+1)**2)*&
                            (ftinv(i,t)+ddot(m,kt(:,i,t),1,rhelp,1)/ft(i,t)**2)
                        end if
                        l0= im
                        call dger(m,m,-ftinv(i,t),kt(:,i,t),1,zt(i,:,(t-1)*timevar(1)+1),1,l0,m) !l0 = I -Kt*Z/Ft
                        call dgemv('t',m,m,1.0d0,l0,m,rrec,1,0.0d0,rhelp,1)
                        rrec = rhelp+zt(i,:,(t-1)*timevar(1)+1)*vt(i,t)*ftinv(i,t)
                        call dgemv('t',m,m,1.0d0,l0,m,rrec1,1,0.0d0,rhelp,1)
                        rrec1=rhelp

                        call dgemm('t','n',m,m,m,1.0d0,l0,m,nrec,m,0.0d0,mm,m) !mm =l0'*nt0
                        call dgemm('n','n',m,m,m,1.0d0,mm,m,l0,m,0.0d0,nrec,m) !nt0 = l0'*nt0*l0
                        call dger(m,m,ftinv(i,t),zt(i,:,(t-1)*timevar(1)+1),1,&
                        zt(i,:,(t-1)*timevar(1)+1),1,nrec,m)  !nt0 = z'z/ft+l0'*nt0*l0
                        call dgemm('n','n',m,m,m,1.0d0,nrec1,m,l0,m,0.0d0,mm,m) !mm = nt1*l0
                        nrec1 = mm
                        call dgemm('n','n',m,m,m,1.0d0,nrec2,m,l0,m,0.0d0,mm,m) !mm = nt1*l0
                        nrec2 = mm !onko oikein?
                    end if
                end if
            end if
        end do
        rt0(:,t) = rrec
        rt1(:,t) = rrec1
        nt0(:,:,t) = nrec
        nt1(:,:,t) = nrec1
        nt2(:,:,t) = nrec2



        if(t.GT.1) then
            call dgemv('t',m,m,1.0d0,tt(:,:,(t-2)*timevar(3)+1),m,rrec,1,0.0d0,rhelp,1)
            rrec = rhelp
            call dgemv('t',m,m,1.0d0,tt(:,:,(t-2)*timevar(3)+1),m,rrec1,1,0.0d0,rhelp,1)
            rrec1 = rhelp
            call dgemm('t','n',m,m,m,1.0d0,tt(:,:,(t-2)*timevar(3)+1),m,nrec,m,0.0d0,mm,m) !n*t
            call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(:,:,(t-2)*timevar(3)+1),m,0.0d0,nrec,m) !n_t,p = t'nt
            call dgemm('t','n',m,m,m,1.0d0,tt(:,:,(t-2)*timevar(3)+1),m,nrec1,m,0.0d0,mm,m) !n*t
            call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(:,:,(t-2)*timevar(3)+1),m,0.0d0,nrec1,m) !n_t,p = t'nt
            call dgemm('t','n',m,m,m,1.0d0,tt(:,:,(t-2)*timevar(3)+1),m,nrec2,m,0.0d0,mm,m) !n*t
            call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(:,:,(t-2)*timevar(3)+1),m,0.0d0,nrec2,m) !n_t,p = t'nt

        end if

        do t=(d-1), 1, -1

            do i = p, 1, -1
                if(ymiss(t,i).EQ.0) then
                    if(finf(i,t).GT. tol) then
                        if(aug .EQ. 1 .AND. dist.EQ.1) then
                            epshat(i,t) = -ht(i,i,(t-1)*timevar(2)+1)*ddot(m,kinf(:,i,t),1,rrec,1)/finf(i,t)
                            call dgemv('n',m,m,1.0d0,nrec,m,kinf(:,i,t),1,0.0d0,rhelp,1)
                            epshatvar(i,t) = ht(i,i,(t-1)*timevar(2)+1)-(ht(i,i,(t-1)*timevar(2)+1)**2)*&
                            ddot(m,kinf(:,i,t),1,rhelp,1)/finf(i,t)**2
                        end if

                        linf = im
                        rhelp = -zt(i,:,(t-1)*timevar(1)+1)*finfinv(i,t)
                        call dger(m,m,1.0d0,kinf(:,i,t),1,rhelp,1,linf,m)

                        rhelp = (kinf(:,i,t)*ft(i,t)*finfinv(i,t)-kt(:,i,t))*finfinv(i,t)
                        l0=0.0d0
                        call dger(m,m,1.0d0,rhelp,1,zt(i,:,(t-1)*timevar(1)+1),1,l0,m) !l0=  (-kt + ft/finf*kinf)*z/finf

                        call dgemv('t',m,m,1.0d0,linf,m,rrec1,1,0.0d0,rhelp,1) !rt1
                        rrec1 = rhelp + zt(i,:,(t-1)*timevar(1)+1)*vt(i,t)*finfinv(i,t)
                        call dgemv('t',m,m,1.0d0,l0,m,rrec,1,1.0d0,rrec1,1)
                        call dgemv('t',m,m,1.0d0,linf,m,rrec,1,0.0d0,rhelp,1) !rt0
                        rrec = rhelp

                        call dgemm('t','n',m,m,m,1.0d0,linf,m,nrec2,m,0.0d0,mm,m) !mm =linf'*nt2
                        call dgemm('n','n',m,m,m,1.0d0,mm,m,linf,m,0.0d0,nrec2,m) !nt2 = linf'*nt2*linf
                        call dger(m,m,-1.0d0*ft(i,t)*finfinv(i,t)**2.0d0,zt(i,:,(t-1)*timevar(1)+1)&
                        ,1,zt(i,:,(t-1)*timevar(1)+1),1,nrec2,m) !nt2 = linf'nt2'linf + z'z*ft/finf^2
                        call dgemm('t','n',m,m,m,1.0d0,l0,m,nrec,m,0.0d0,mm,m) !mm= nt0*l0
                        call dgemm('n','n',m,m,m,1.0d0,mm,m,l0,m,1.0d0,nrec2,m) !nt2 = linf'nt2'linf + z'z*ft/finf^2 + l0'*nt0*l0

                        call dgemm('t','n',m,m,m,1.0d0,linf,m,nrec1,m,0.0d0,mm,m) !mm = linf'*nt1
                        call dgemm('n','n',m,m,m,1.0d0,mm,m,l0,m,0.0d0,mm2,m) !nt2 = nt2 + linf'*nt1*l0
                        nrec2 = nrec2 +mm2 + transpose(mm2)

                        call dgemm('t','n',m,m,m,1.0d0,linf,m,nrec1,m,0.0d0,mm,m) !mm = linf'*nt1
                        call dgemm('n','n',m,m,m,1.0d0,mm,m,linf,m,0.0d0,nrec1,m) !nt1 = mm*linf

                        call dger(m,m,finfinv(i,t),zt(i,:,(t-1)*timevar(1)+1),1,zt(i,:,(t-1)*timevar(1)+1),1,nrec1,m)
                        !nt1 = linf'nt1'linf + z'z/finf
                        call dgemm('t','n',m,m,m,1.0d0,l0,m,nrec,m,0.0d0,mm,m) !mm= nt0*linf
                        call dgemm('n','n',m,m,m,1.0d0,mm,m,linf,m,1.0d0,nrec1,m) !nt1 = l0'*nt0*linf+ linf'nt1*linf + z'z/finf

                        call dgemm('t','n',m,m,m,1.0d0,linf,m,nrec,m,0.0d0,mm,m) !mm= nt0*linf
                        call dgemm('n','n',m,m,m,1.0d0,mm,m,linf,m,0.0d0,nrec,m)

                    else
                        if(ft(i,t).GT.0.0d0) then !lis?tty 12.1.2012
                            if(aug .EQ. 1 .AND. dist.EQ.1) then
                                epshat(i,t) = ht(i,i,(t-1)*timevar(2)+1)*(vt(i,t)/ft(i,t)-&
                                ddot(m,kt(:,i,t),1,rrec,1)/ft(i,t))
                                call dgemv('n',m,m,1.0d0,nrec,m,kt(:,i,t),1,0.0d0,rhelp,1)
                                epshatvar(i,t) = ht(i,i,(t-1)*timevar(2)+1)-(ht(i,i,(t-1)*timevar(2)+1)**2)*&
                                (ftinv(i,t)+ddot(m,kt(:,i,t),1,rhelp,1)/ft(i,t)**2 )
                            end if
                            l0= im
                            call dger(m,m,-ftinv(i,t),kt(:,i,t),1,zt(i,:,(t-1)*timevar(1)+1),1,l0,m) !l0 = I -Kt*Z/Ft
                            call dgemv('t',m,m,1.0d0,l0,m,rrec,1,0.0d0,rhelp,1)
                            rrec = rhelp+zt(i,:,(t-1)*timevar(1)+1)*vt(i,t)*ftinv(i,t)
                            call dgemv('t',m,m,1.0d0,l0,m,rrec1,1,0.0d0,rhelp,1)
                            rrec1=rhelp

                            call dgemm('t','n',m,m,m,1.0d0,l0,m,nrec,m,0.0d0,mm,m) !mm =l0'*nt0
                            call dgemm('n','n',m,m,m,1.0d0,mm,m,l0,m,0.0d0,nrec,m) !nt0 = l0'*nt0*l0
                            call dger(m,m,ftinv(i,t),zt(i,:,(t-1)*timevar(1)+1),1,&
                            zt(i,:,(t-1)*timevar(1)+1),1,nrec,m)  !nt0 = z'z/ft+l0'*nt0*l0
                            call dgemm('n','n',m,m,m,1.0d0,nrec1,m,l0,m,0.0d0,mm,m) !mm = nt1*l0
                            nrec1 = mm
                            call dgemm('n','n',m,m,m,1.0d0,nrec2,m,l0,m,0.0d0,mm,m) !mm = nt1*l0
                            nrec2 = mm
                        end if
                    end if
                end if
            end do


            rt0(:,t) = rrec
            rt1(:,t) = rrec1
            nt0(:,:,t) = nrec
            nt1(:,:,t) = nrec1
            nt2(:,:,t) = nrec2


            if(t.GT.1) then
                call dgemv('t',m,m,1.0d0,tt(:,:,(t-2)*timevar(3)+1),m,rrec,1,0.0d0,rhelp,1)
                rrec = rhelp
                call dgemv('t',m,m,1.0d0,tt(:,:,(t-2)*timevar(3)+1),m,rrec1,1,0.0d0,rhelp,1)
                rrec1 = rhelp
                call dgemm('t','n',m,m,m,1.0d0,tt(:,:,(t-2)*timevar(3)+1),m,nrec,m,0.0d0,mm,m) !n*t
                call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(:,:,(t-2)*timevar(3)+1),m,0.0d0,nrec,m) !n_t,p = t'nt
                call dgemm('t','n',m,m,m,1.0d0,tt(:,:,(t-2)*timevar(3)+1),m,nrec1,m,0.0d0,mm,m) !n*t
                call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(:,:,(t-2)*timevar(3)+1),m,0.0d0,nrec1,m) !n_t,p = t'nt
                call dgemm('t','n',m,m,m,1.0d0,tt(:,:,(t-2)*timevar(3)+1),m,nrec2,m,0.0d0,mm,m) !n*t
                call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(:,:,(t-2)*timevar(3)+1),m,0.0d0,nrec2,m) !n_t,p = t'nt
            end if


        end do
    end if

    if(state.EQ.1) then
        do t = 1, d
            ahat(:,t) = at(:,t)
            call dsymv('u',m,1.0d0,pt(:,:,t),m,rt0(:,t),1,1.0d0,ahat(:,t),1)
            call dsymv('u',m,1.0d0,pinf(:,:,t),m,rt1(:,t),1,1.0d0,ahat(:,t),1)

            vvt(:,:,t) = pt(:,:,t)
            call dsymm('l','u',m,m,1.0d0,pt(:,:,t),m,nt0(:,:,t),m,0.0d0,mm,m) !mm = pt*nt0
            call dsymm('r','u',m,m,-1.0d0,pt(:,:,t),m,mm,m,1.0d0,vvt(:,:,t),m) !vvt = pt - pt*nt0*pt
            call dsymm('l','u',m,m,1.0d0,pinf(:,:,t),m,nt1(:,:,t),m,0.0d0,mm,m) !mm = pinf*nt1
            call dsymm('r','u',m,m,-1.0d0,pt(:,:,t),m,mm,m,0.0d0,mm2,m) !mm2 = -pinf*nt1*pt
            vvt(:,:,t) = vvt(:,:,t) + mm2 + transpose(mm2) !vvt = pt - pt*nt0*pt  -pinf*nt1*pt - t(pinf*nt1*pt)
            call dsymm('l','u',m,m,1.0d0,pinf(:,:,t),m,nt2(:,:,t),m,0.0d0,mm,m) !mm = pinf*nt2
            call dsymm('r','u',m,m,-1.0d0,pinf(:,:,t),m,mm,m,1.0d0,vvt(:,:,t),m) !vvt = vvt - pinf*nt2*pinf
        end do
        do t = d+1, n
            ahat(:,t) = at(:,t)
            call dsymv('u',m,1.0d0,pt(:,:,t),m,rt(:,t),1,1.0d0,ahat(:,t),1) !ahat = ahat+pt*r_t-1
            call dsymm('l','u',m,m,1.0d0,pt(:,:,t),m,nt(:,:,t),m,0.0d0,mm,m) !pt*n_t-1
            mm = im - mm
            call dsymm('r','u',m,m,1.0d0,pt(:,:,t),m,mm,m,0.0d0,vvt(:,:,t),m) !pt*n_t-1*pt
        end do
        if(m > 1) then
            do t=1, n
                do i=1,m-1
                    vvt((i+1):m,i,t) =vvt(i,(i+1):m,t)
                end do
            end do
        end if
    end if

    if(dist.EQ.1) then
        do t = 1, d
            call dgemv('t',m,r,1.0d0,rtv(:,:,(t-1)*timevar(4)+1),m,rt0(:,t+1),1,0.0d0,help,1)
            call dsymv('u',r,1.0d0,qt(:,:,(t-1)*timevar(5)+1),r,help,1,0.0d0,etahat(:,t),1)
            etahatvar(:,:,t) = qt(:,:,(t-1)*timevar(5)+1)
            call dsymm('r','u',m,r,1.0d0,qt(:,:,(t-1)*timevar(5)+1),r,rtv(:,:,(t-1)*timevar(4)+1),m,0.0d0,mr,m)
            !call dgemm('n','n',m,r,m,1.0d0,nt0(:,:,t+1),m,mr,m,0.0d0,mr2,m)
             call dsymm('l','u',m,r,1.0d0,nt0(:,:,t+1),m,mr,m,0.0d0,mr2,m)
            call dgemm('t','n',r,r,m,-1.0d0,mr,m,mr2,m,1.0d0,etahatvar(:,:,t),r)
        end do
        do t = d+1, n
            call dgemv('t',m,r,1.0d0,rtv(:,:,(t-1)*timevar(4)+1),m,rt(:,t+1),1,0.0d0,help,1)
            call dsymv('u',r,1.0d0,qt(:,:,(t-1)*timevar(5)+1),r,help,1,0.0d0,etahat(:,t),1)
            etahatvar(:,:,t) = qt(:,:,(t-1)*timevar(5)+1)
            call dsymm('r','u',m,r,1.0d0,qt(:,:,(t-1)*timevar(5)+1),r,rtv(:,:,(t-1)*timevar(4)+1),m,0.0d0,mr,m)
            !call dgemm('n','n',m,r,m,1.0d0,nt(:,:,t+1),m,mr,m,0.0d0,mr2,m)
            call dsymm('l','u',m,r,1.0d0,nt(:,:,t+1),m,mr,m,0.0d0,mr2,m)
            call dgemm('t','n',r,r,m,-1.0d0,mr,m,mr2,m,1.0d0,etahatvar(:,:,t),r)
        end do
    end if

    if(signal.EQ.1) then
        if(ldlsignal .EQ. 1) then
            if(state .EQ. 1) then
                do t = 1, n
                    call dgemv('n',p,m,1.0d0,zorig(:,:,(t-1)*zorigtv+1),p,ahat(:,t),1,0.0d0,thetahat(:,t),1)
                    call dsymm('r','u',p,m,1.0d0,vvt(:,:,t),m,zorig(:,:,(t-1)*zorigtv+1),p,0.0d0,pm,p)
                    call dgemm('n','t',p,p,m,1.0d0,pm,p,zorig(:,:,(t-1)*zorigtv+1),p,0.0d0,thetahatvar(:,:,t),p)
                end do
            else
                do t = 1, d
                    call dcopy(m,at(:,t),1,rrec,1) !ahat = at
                    call dsymv('u',m,1.0d0,pt(:,:,t),m,rt0(:,t),1,1.0d0,rrec,1) !ahat = at + pt * rt0_t
                    call dsymv('u',m,1.0d0,pinf(:,:,t),m,rt1(:,t),1,1.0d0,rrec,1) !ahat = at + pt * rt0_t + pinf*rt1_t
                    call dgemv('n',p,m,1.0d0,zorig(:,:,(t-1)*zorigtv+1),p,rrec,1,0.0d0,thetahat(:,t),1)

                    nrec = pt(:,:,t)
                    call dsymm('l','u',m,m,1.0d0,pt(:,:,t),m,nt0(:,:,t),m,0.0d0,mm,m) !mm = pt*nt0
                    call dsymm('r','u',m,m,-1.0d0,pt(:,:,t),m,mm,m,1.0d0,nrec,m) !vvt = pt - pt*nt0*pt
                    call dsymm('l','u',m,m,1.0d0,pinf(:,:,t),m,nt1(:,:,t),m,0.0d0,mm,m) !mm = pinf*nt1
                    call dsymm('r','u',m,m,-1.0d0,pt(:,:,t),m,mm,m,0.0d0,mm2,m) !mm2 = -pinf*nt1*pt
                    nrec = nrec + mm2 + transpose(mm2) !vvt = pt - pt*nt0*pt  -pinf*nt1*pt - t(pinf*nt1*pt)
                    call dsymm('l','u',m,m,1.0d0,pinf(:,:,t),m,nt2(:,:,t),m,0.0d0,mm,m) !mm = pinf*nt2
                    call dsymm('r','u',m,m,-1.0d0,pinf(:,:,t),m,mm,m,1.0d0,nrec,m) !vvt = vvt - pinf*nt2*pinf
                    call dsymm('r','u',p,m,1.0d0,nrec,m,zorig(:,:,(t-1)*zorigtv+1),p,0.0d0,pm,p)
                    call dgemm('n','t',p,p,m,1.0d0,pm,p,zorig(:,:,(t-1)*zorigtv+1),p,0.0d0,thetahatvar(:,:,t),p)
                end do
                do t = d+1, n
                    call dcopy(m,at(:,t),1,rrec,1) !ahat = at
                    call dsymv('u',m,1.0d0,pt(:,:,t),m,rt(:,t),1,1.0d0,rrec,1) !ahat = ahat+pt*r_t-1
                    call dgemv('n',p,m,1.0d0,zorig(:,:,(t-1)*zorigtv+1),p,rrec,1,0.0d0,thetahat(:,t),1)
                    nrec = pt(:,:,t)
                    call dsymm('l','u',m,m,1.0d0,pt(:,:,t),m,nt(:,:,t),m,0.0d0,mm,m) !pt*n_t-1
                    call dsymm('r','u',m,m,-1.0d0,pt(:,:,t),m,mm,m,1.0d0,nrec,m) !pt*n_t-1*pt
                    call dsymm('r','u',p,m,1.0d0,nrec,m,zorig(:,:,(t-1)*zorigtv+1),p,0.0d0,pm,p)
                    call dgemm('n','t',p,p,m,1.0d0,pm,p,zorig(:,:,(t-1)*zorigtv+1),p,0.0d0,thetahatvar(:,:,t),p)
                end do
            end if
        else

            if(state .EQ. 1) then
                do t = 1, n
                    call dgemv('n',p,m,1.0d0,zt(:,:,(t-1)*timevar(1)+1),p,ahat(:,t),1,0.0d0,thetahat(:,t),1)
                    call dsymm('r','u',p,m,1.0d0,vvt(:,:,t),m,zt(:,:,(t-1)*timevar(1)+1),p,0.0d0,pm,p)
                    call dgemm('n','t',p,p,m,1.0d0,pm,p,zt(:,:,(t-1)*timevar(1)+1),p,0.0d0,thetahatvar(:,:,t),p)
                end do
            else
                do t = 1, d
                    call dcopy(m,at(:,t),1,rrec,1) !ahat = at
                    call dsymv('u',m,1.0d0,pt(:,:,t),m,rt0(:,t),1,1.0d0,rrec,1) !ahat = at + pt * rt0_t
                    call dsymv('u',m,1.0d0,pinf(:,:,t),m,rt1(:,t),1,1.0d0,rrec,1) !ahat = at + pt * rt0_t + pinf*rt1_t
                    call dgemv('n',p,m,1.0d0,zt(:,:,(t-1)*timevar(1)+1),p,rrec,1,0.0d0,thetahat(:,t),1)

                    nrec = pt(:,:,t)
                    call dsymm('l','u',m,m,1.0d0,pt(:,:,t),m,nt0(:,:,t),m,0.0d0,mm,m) !mm = pt*nt0
                    call dsymm('r','u',m,m,-1.0d0,pt(:,:,t),m,mm,m,1.0d0,nrec,m) !vvt = pt - pt*nt0*pt
                    call dsymm('l','u',m,m,1.0d0,pinf(:,:,t),m,nt1(:,:,t),m,0.0d0,mm,m) !mm = pinf*nt1
                    call dsymm('r','u',m,m,-1.0d0,pt(:,:,t),m,mm,m,0.0d0,mm2,m) !mm2 = -pinf*nt1*pt
                    nrec = nrec + mm2 + transpose(mm2) !vvt = pt - pt*nt0*pt  -pinf*nt1*pt - t(pinf*nt1*pt)
                    call dsymm('l','u',m,m,1.0d0,pinf(:,:,t),m,nt2(:,:,t),m,0.0d0,mm,m) !mm = pinf*nt2
                    call dsymm('r','u',m,m,-1.0d0,pinf(:,:,t),m,mm,m,1.0d0,nrec,m) !vvt = vvt - pinf*nt2*pinf
                    call dsymm('r','u',p,m,1.0d0,nrec,m,zt(:,:,(t-1)*timevar(1)+1),p,0.0d0,pm,p)
                    call dgemm('n','t',p,p,m,1.0d0,pm,p,zt(:,:,(t-1)*timevar(1)+1),p,0.0d0,thetahatvar(:,:,t),p)
                end do
                do t = d+1, n
                    call dcopy(m,at(:,t),1,rrec,1) !ahat = at
                    call dsymv('u',m,1.0d0,pt(:,:,t),m,rt(:,t),1,1.0d0,rrec,1) !ahat = ahat+pt*r_t-1
                    call dgemv('n',p,m,1.0d0,zt(:,:,(t-1)*timevar(1)+1),p,rrec,1,0.0d0,thetahat(:,t),1)
                    nrec = pt(:,:,t)
                    call dsymm('l','u',m,m,1.0d0,pt(:,:,t),m,nt(:,:,t),m,0.0d0,mm,m) !pt*n_t-1
                    call dsymm('r','u',m,m,-1.0d0,pt(:,:,t),m,mm,m,1.0d0,nrec,m) !pt*n_t-1*pt
                    call dsymm('r','u',p,m,1.0d0,nrec,m,zt(:,:,(t-1)*timevar(1)+1),p,0.0d0,pm,p)
                    call dgemm('n','t',p,p,m,1.0d0,pm,p,zt(:,:,(t-1)*timevar(1)+1),p,0.0d0,thetahatvar(:,:,t),p)
                end do
            end if
        end if
    end if

end subroutine gsmoothall

