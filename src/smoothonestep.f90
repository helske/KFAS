!non-diffuse smoothing for single time point
subroutine smooth1step(ymiss, zt, ht, tt, rtv, qt, vt, ft,kt,&
im,p,m,r,j,rt,etahat,epshat,needeps)

    implicit none

    logical, intent(in) :: needeps
    integer, intent(in) ::  p, m,j,r
    integer ::  i
    integer, intent(in), dimension(p) :: ymiss
    double precision, intent(in), dimension(m,p) :: zt
    double precision, intent(in), dimension(p,p) :: ht
    double precision, intent(in), dimension(m,m) :: tt
    double precision, intent(in), dimension(m,r) :: rtv
    double precision, intent(in), dimension(r,r) :: qt
    double precision, intent(inout), dimension(m) :: rt
    double precision, intent(inout), dimension(r) :: etahat
    double precision, intent(inout), dimension(p) :: epshat
    double precision, intent(in), dimension(p) :: vt,ft
    double precision, intent(in), dimension(m,p) :: kt
    double precision, intent(in), dimension(m,m) :: im
    double precision, dimension(m,m) :: l0
    double precision, dimension(m) :: mhelp
    double precision, dimension(r) :: rhelp
    double precision :: finv

    double precision, external :: ddot

    external dgemv, dsymv, dger

    call dgemv('t',m,r,1.0d0,rtv,m,rt,1,0.0d0,rhelp,1)
    call dsymv('l',r,1.0d0,qt,r,rhelp,1,0.0d0,etahat,1)
    call dgemv('t',m,m,1.0d0,tt,m,rt,1,0.0d0,mhelp,1)
    rt = mhelp
    do i = p, j , -1
        if(ymiss(i) .EQ. 0) then
            if(ft(i) .GT. 0.0d0) then
                finv = 1.0d0/ft(i)
                if(needeps) then
                    epshat(i) = ht(i,i)*(vt(i)-ddot(m,kt(:,i),1,rt,1))*finv
                end if
                l0 = im
                call dger(m,m,-finv,kt(:,i),1,zt(:,i),1,l0,m)
                call dgemv('t',m,m,1.0d0,l0,m,rt,1,0.0d0,mhelp,1)
                rt = mhelp + vt(i)*finv*zt(:,i)
            end if
        end if
    end do

end subroutine smooth1step

subroutine dsmooth1step(ymiss, zt, ht, tt, rtv, qt, vt, ft,kt,&
im,p,m,r,j,rt,rt1,finf,kinf,etahat,epshat,needeps)

    implicit none

    logical, intent(in) :: needeps
    integer, intent(in) ::  p, m,j,r
    integer ::  i
    integer, intent(in), dimension(p) :: ymiss
    double precision, intent(in), dimension(m,p) :: zt
    double precision, intent(in), dimension(p,p) :: ht
    double precision, intent(in), dimension(m,m) :: tt
    double precision, intent(in), dimension(m,r) :: rtv
    double precision, intent(in), dimension(r,r) :: qt
    double precision, intent(inout), dimension(m) :: rt,rt1
    double precision, intent(inout), dimension(r) :: etahat
    double precision, intent(inout), dimension(p) :: epshat
    double precision, intent(in), dimension(p) :: vt,ft,finf
    double precision, intent(in), dimension(m,p) :: kt,kinf
    double precision, intent(in), dimension(m,m) :: im
    double precision, dimension(m,m) :: l0,linf
    double precision, dimension(m) :: mhelp
    double precision, dimension(r) :: rhelp
    double precision :: finv

    double precision, external :: ddot

    external dgemv, dsymv, dger

    if(j .EQ. p) then
        call dgemv('t',m,r,1.0d0,rtv,m,rt,1,0.0d0,rhelp,1)
        call dsymv('l',r,1.0d0,qt,r,rhelp,1,0.0d0,etahat,1)
        call dgemv('t',m,m,1.0d0,tt,m,rt,1,0.0d0,mhelp,1)
        rt = mhelp
        call dgemv('t',m,m,1.0d0,tt,m,rt1,1,0.0d0,mhelp,1)
        rt1 = mhelp
    end if

    do i = j, 1, -1
        if(ymiss(i) .EQ. 0) then
            if(finf(i) .GT. 0.0d0) then
                finv = 1.0d0/finf(i)
                if(needeps) then
                    epshat(i) = -ht(i,i)*ddot(m,kinf(:,i),1,rt,1)*finv
                end if
                linf = im
                call dger(m,m,-finv,kinf(:,i),1,zt(:,i),1,linf,m)

                mhelp = kinf(:,i)*ft(i)*finv - kt(:,i)
                l0 = 0.0d0
                call dger(m,m,finv,mhelp,1,zt(:,i),1,l0,m)

                call dgemv('t',m,m,1.0d0,linf,m,rt1,1,0.0d0,mhelp,1)
                rt1 = mhelp
                call dgemv('t',m,m,1.0d0,l0,m,rt,1,1.0d0,rt1,1)
                rt1 = rt1 + vt(i)*finv*zt(:,i)

                call dgemv('t',m,m,1.0d0,linf,m,rt,1,0.0d0,mhelp,1)
                rt = mhelp
            else
                if(ft(i) .GT. 0.0d0) then
                    finv = 1.0d0/ft(i)
                    if(needeps) then
                        epshat(i) = ht(i,i)*(vt(i)-ddot(m,kt(:,i),1,rt,1))*finv
                    end if

                    l0 = im
                    call dger(m,m,-finv,kt(:,i),1,zt(:,i),1,l0,m)

                    call dgemv('t',m,m,1.0d0,l0,m,rt,1,0.0d0,mhelp,1)
                    rt = mhelp + vt(i)*finv*zt(:,i)

                    call dgemv('t',m,m,1.0d0,l0,m,rt1,1,0.0d0,mhelp,1)
                    rt1 = mhelp
                end if
            end if
        end if
    end do

end subroutine dsmooth1step
