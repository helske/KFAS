!fast filtering algorithm used in simulation filter
subroutine filtersimfast(yt, ymiss, timevar, zt,tt, &
a1, ft,kt,finf, kinf, dt, jt, p, m, n,tol,at)

    implicit none

    integer, intent(in) ::  p, m,n,dt,jt
    integer ::  t, i,d,j
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(p,n) :: ft,finf
    double precision, intent(in), dimension(m,p,n) :: kt,kinf
    double precision, intent(in) :: tol
    double precision, intent(inout), dimension(m,n+1) :: at
    double precision, dimension(p,n) :: vt
    double precision, dimension(m) :: arec
    double precision, external :: ddot

    external dgemv

    j=0
    d=0
    if(dt.GT.0) then
        !diffuse filtering begins
        arec = a1
        diffuse: do while(d .LT. (dt-1))
            d = d+1
            do j=1, p
                if(ymiss(d,j).EQ.0) then
                    vt(j,d) = yt(d,j) - ddot(m,zt(j,:,(d-1)*timevar(1)+1),1,arec,1)
                    if (finf(j,d) .GT. 0.0d0) then
                        arec = arec + vt(j,d)/finf(j,d)*kinf(:,j,d)
                    else
                        if(ft(j,d) .GT. 0.0d0) then
                            arec = arec + vt(j,d)/ft(j,d)*kt(:,j,d)
                        end if
                    end if
                end if
            end do
            call dgemv('n',m,m,1.0d0,tt(:,:,(d-1)*timevar(3)+1),m,arec,1,0.0d0,at(:,d+1),1)
            arec = at(:,d+1)
        end do diffuse

        d = dt
        do j=1, jt
            if(ymiss(d,j).EQ.0) then
                vt(j,d) = yt(d,j) - ddot(m,zt(j,:,(d-1)*timevar(1)+1),1,arec,1)
      
                if (finf(j,d) .GT. 0.0d0) then
                    arec = arec + vt(j,d)/finf(j,d)*kinf(:,j,d)
                else
                    if(ft(j,d) .GT. 0.0d0) then
                        arec = arec + vt(j,d)/ft(j,d)*kt(:,j,d)
                    end if
                end if
            end if
        end do
        !non-diffuse filtering begins
        do i = jt+1, p
            if(ymiss(d,i).EQ.0) then
                vt(i,d) = yt(d,i) - ddot(m,zt(i,:,(d-1)*timevar(1)+1),1,arec,1)
                if (ft(i,d) .GT. 0.0d0) then
                    arec = arec + vt(i,d)/ft(i,d)*kt(:,i,d)
                end if
            end if
        end do
   
        call dgemv('n',m,m,1.0d0,tt(:,:,(d-1)*timevar(3)+1),m,arec,1,0.0d0,at(:,d+1),1)
        arec = at(:,d+1)
    end if

    if(dt.LT.n) then
        !Non-diffuse filtering continues from t=d+1, i=1
        if(dt.EQ.0) then
            arec = a1
        end if
        do t = dt+1, n
            do i = 1, p
                if(ymiss(t,i).EQ.0) then
                    vt(i,t) = yt(t,i) - ddot(m,zt(i,:,(t-1)*timevar(1)+1),1,arec,1)
                    if (ft(i,t) .GT. 0.0d0) then
                        arec = arec + vt(i,t)/ft(i,t)*kt(:,i,t)
                    end if
                end if
            end do
   
            call dgemv('n',m,m,1.0d0,tt(:,:,(t-1)*timevar(3)+1),m,arec,1,0.0d0,at(:,t+1),1)
            arec = at(:,t+1)
        end do

    end if

end subroutine filtersimfast
