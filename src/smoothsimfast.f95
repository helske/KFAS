subroutine smoothsimfast(yt, ymiss, timevar, zt, ht,tt, rtv,qt,a1, ft,kt,&
finf, kinf, dt, jt, p, m, n,r,tol,epshat,etahat,rt0,rt1,needeps)

    implicit none

    logical, intent(in) :: needeps
    integer, intent(in) ::  p, m, r,n,dt,jt
    integer ::  t, i,d,j
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(p,p,(n-1)*timevar(2)+1) :: ht
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rtv
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(p,n) :: ft,finf
    double precision, intent(in), dimension(m,p,n) :: kt,kinf
    double precision, intent(in) :: tol
    double precision, intent(inout), dimension(p,n) :: epshat
    double precision, intent(inout), dimension(r,n) :: etahat
    double precision, dimension(p,n) :: vt
    double precision, dimension(m) :: arec,rrec,rrec1,rhelp,help
    double precision, dimension(m,m) :: im,linf,l0,lt
    double precision, intent(inout), dimension(m) :: rt0,rt1
    double precision :: meps
    double precision, external :: ddot

 meps = tiny(meps)

    j=0
    d=0
    if(dt.GT.0) then
        arec = a1
        diffuse: do while(d .LT. (dt-1))
            d = d+1
            do j=1, p
                if(ymiss(d,j).EQ.0) then
                    vt(j,d) = yt(d,j) - ddot(m,zt(j,:,(d-1)*timevar(1)+1),1,arec,1) !arec
                    if (finf(j,d) .GT. tol) then
                        call daxpy(m,vt(j,d)/finf(j,d),kinf(:,j,d),1,arec,1) !a_rec = a_rec + kinf(:,i,t)*vt(:,t)/finf(j,d)
                    else
                        if(ft(j,d) .GT. meps) then
                            call daxpy(m,vt(j,d)/ft(j,d),kt(:,j,d),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)/ft(i,t)
                        end if
                    end if
                end if
            end do
           
            call dgemv('n',m,m,1.0d0,tt(:,:,(d-1)*timevar(3)+1),m,arec,1,0.0d0,help,1)  !at(:,t+1) = matmul(tt,a_rec)
            arec = help
        end do diffuse

        d = dt
        do j=1, jt
            if(ymiss(d,j).EQ.0) then
                vt(j,d) = yt(d,j) - ddot(m,zt(j,:,(d-1)*timevar(1)+1),1,arec,1) !arec
      
                if (finf(j,d) .GT. tol) then
                    call daxpy(m,vt(j,d)/finf(j,d),kinf(:,j,d),1,arec,1) !a_rec = a_rec + kinf(:,i,t)*vt(:,t)/finf(j,d)
                else
                    if(ft(j,d) .GT. meps ) then
                        call daxpy(m,vt(j,d)/ft(j,d),kt(:,j,d),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)/ft(i,t)
                    end if
                end if
            end if
        end do
   
  
        !non-diffuse filtering begins
 
        do i = jt+1, p
            if(ymiss(d,i).EQ.0) then
                vt(i,d) = yt(d,i) - ddot(m,zt(i,:,(d-1)*timevar(1)+1),1,arec,1) !vt
                if (ft(i,d) .GT.  meps) then !ft.NE.0
                    call daxpy(m,vt(i,d)/ft(i,d),kt(:,i,d),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)
                end if
            end if
        end do
   
        call dgemv('n',m,m,1.0d0,tt(:,:,(d-1)*timevar(3)+1),m,arec,1,0.0d0,help,1)  !at(:,t+1) = matmul(tt,a_rec)
  
   
        arec =help
   
    end if

    if(dt.LT.n) then

        !Non-diffuse filtering continues from t=d+1, i=1
        if(dt.EQ.0) then
            arec = a1
        end if
        do t = dt+1, n
            do i = 1, p
                if(ymiss(t,i).EQ.0) then
                    vt(i,t) = yt(t,i) - ddot(m,zt(i,:,(t-1)*timevar(1)+1),1,arec,1) !variate vt
                    if (ft(i,t) .GT.  meps) then !ft.NE.0
                        call daxpy(m,vt(i,t)/ft(i,t),kt(:,i,t),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)
                    end if
                end if
            end do
   
            call dgemv('n',m,m,1.0d0,tt(:,:,(t-1)*timevar(3)+1),m,arec,1,0.0d0,help,1)  !at(:,t+1) = matmul(tt,a_rec)           
   
            arec =help

        end do

    end if


    !smoothing begins

    im = 0.0d0
    do i = 1, m
        im(i,i) = 1.0d0
    end do

    rrec = 0.0d0

    !rt(:,n+1) = 0.0d0

    do t = n, dt+1, -1 !do until diffuse tts
        call dgemv('t',m,r,1.0d0,rtv(:,:,(t-1)*timevar(4)+1),m,rrec,1,0.0d0,help,1)
        call dsymv('l',r,1.0d0,qt(:,:,(t-1)*timevar(5)+1),r,help,1,0.0d0,etahat(:,t),1)
        call dgemv('t',m,m,1.0d0,tt(:,:,(t-1)*timevar(3)+1),m,rrec,1,0.0d0,rhelp,1) !r_t,p=t_t-1'*r_t+1
        rrec = rhelp
        do i = p, 1 , -1
            if(ymiss(t,i).EQ.0) then
                if(ft(i,t) .GT. meps) then
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

    if(dt.GT.0) then
        t=dt
        call dgemv('t',m,r,1.0d0,rtv(:,:,(t-1)*timevar(4)+1),m,rrec,1,0.0d0,help,1)
        call dsymv('l',r,1.0d0,qt(:,:,(t-1)*timevar(5)+1),r,help,1,0.0d0,etahat(:,t),1)
   
        call dgemv('t',m,m,1.0d0,tt(:,:,(t-1)*timevar(3)+1),m,rrec,1,0.0d0,rhelp,1,1) !tarkiSta tOimivUUS!
        rrec = rhelp       
  
        do i = p, (jt+1) , -1

            if(ymiss(t,i).EQ.0) then
                if(ft(i,t) .GT. meps) then
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
        do i = jt, 1, -1
            if(ymiss(t,i).EQ.0) then
                if(finf(i,t).GT.tol) then
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
                    call dcopy(m,rhelp,1,rrec1,1)
                    call dgemv('t',m,m,1.0d0,l0,m,rrec,1,1.0d0,rrec1,1)
                    call daxpy(m,(vt(i,t)/finf(i,t)),zt(i,:,(t-1)*timevar(1)+1),1,rrec1,1)
                    call dgemv('t',m,m,1.0d0,linf,m,rrec,1,0.0d0,rhelp,1) !rt0
                    rrec = rhelp


                else

                    if(ft(i,t).GT.meps) then
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

        do t=(dt-1), 1, -1
            call dgemv('t',m,r,1.0d0,rtv(:,:,(t-1)*timevar(4)+1),m,rrec,1,0.0d0,help,1)
            call dsymv('l',r,1.0d0,qt(:,:,(t-1)*timevar(5)+1),r,help,1,0.0d0,etahat(:,t),1)
            call dgemv('t',m,m,1.0d0,tt(:,:,(t-1)*timevar(3)+1),m,rrec,1,0.0d0,rhelp,1,1)
            rrec = rhelp
            call dgemv('t',m,m,1.0d0,tt(:,:,(t-1)*timevar(3)+1),m,rrec1,1,0.0d0,rhelp,1,1)
            rrec1 = rhelp

            do i = p, 1, -1
                if(ymiss(t,i).EQ.0) then
                    if(finf(i,t).GT. tol) then
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
                        call dcopy(m,rhelp,1,rrec1,1)
                        call dgemv('t',m,m,1.0d0,l0,m,rrec,1,1.0d0,rrec1,1)
                        call daxpy(m,vt(i,t)/finf(i,t),zt(i,:,(t-1)*timevar(1)+1),1,rrec1,1)
                        call dgemv('t',m,m,1.0d0,linf,m,rrec,1,0.0d0,rhelp,1) !rt0
                        rrec = rhelp

                    else
                        if(ft(i,t).GT. meps) then
                            if(needeps) then
                                epshat(i,t) = ht(i,i,(t-1)*timevar(2)+1)*(vt(i,t)-ddot(m,kt(:,i,t),1,rrec,1))/ft(i,t)
                            end if
                            lt= im
                            call dger(m,m,-1.0d0/ft(i,t),kt(:,i,t),1,zt(i,:,(t-1)*timevar(1)+1),1,lt,m) !lt = I -Kt*Z/Ft
                            call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1) !oli beta 1.0d0!!!!... JA miinusmerkki
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

end subroutine smoothsimfast
