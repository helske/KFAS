! functions for computing p(theta)
subroutine pthetafirst(yt, timevar, zt, tt, rqr, a1, p1, p1inf,&
p, m, r, n, lik, tol,rankp2,kt,kinf,ft,finf,d,j)


    implicit none

    integer, intent(in) ::  p, m, r, n
    integer, intent(inout) :: rankp2,d,j
    integer ::  t, i,tv,rankp
    integer, intent(in), dimension(5) :: timevar
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, intent(in) :: tol
    double precision, intent(inout) :: lik
    double precision, dimension(m) :: at,arec
    double precision, dimension(p) :: vt
    double precision, intent(inout), dimension(p,n) :: ft,finf
    double precision, intent(inout), dimension(m,p,n) :: kt,kinf
    double precision, dimension(m,m) :: pt,pinf,prec,pirec,mm
    double precision, external :: ddot
    double precision :: meps
    double precision, intent(inout), dimension(m,m,(n-1)*max(timevar(4),timevar(5))+1) :: rqr

    tv= max(timevar(4),timevar(5))

    meps = tiny(meps)


    pinf=p1inf
    rankp = rankp2
    j=0
    d=0
    ! Diffuse initialization
    if(maxval(p1inf) .GT.  0.0d0) then

        pt = p1
        prec = pt
        pirec = pinf
        at = a1
        arec = a1
        diffuse: do while(d .LT. n)
            d = d+1
            do j=1, p
                vt(j) = yt(d,j) - ddot(m,zt(j,:,(d-1)*timevar(1)+1),1,arec,1)
                call dsymv('u',m,1.0d0,prec,m,zt(j,:,(d-1)*timevar(1)+1),1,0.0d0,kt(:,j,d),1) ! kt_t,i = pt_t,i*t(z_t,i)
                ft(j,d) = ddot(m,zt(j,:,(d-1)*timevar(1)+1),1,kt(:,j,d),1)
                call dsymv('u',m,1.0d0,pirec,m,zt(j,:,(d-1)*timevar(1)+1),1,0.0d0,kinf(:,j,d),1) ! kinf_t,i = pinf_t,i*t(z_t,i)
                finf(j,d) = ddot(m,zt(j,:,(d-1)*timevar(1)+1),1,kinf(:,j,d),1)
                if (finf(j,d) .GT. tol) then
                    call daxpy(m,vt(j)/finf(j,d),kinf(:,j,d),1,arec,1) !a_rec = a_rec + kinf(:,i,t)*vt(:,t)/finf(j,d)
                    call dsyr('u',m,ft(j,d)/(finf(j,d)**2),kinf(:,j,d),1,prec,m) !prec = prec +  kinf*kinf'*ft/finf^2
                    call dsyr2('u',m,-1.0d0/finf(j,d),kt(:,j,d),1,kinf(:,j,d),1,prec,m) !prec = prec -(kt*kinf'+kinf*kt')/finf
                    call dsyr('u',m,-1.0d0/finf(j,d),kinf(:,j,d),1,pirec,m) !pirec = pirec -kinf*kinf'/finf
                    lik = lik - 0.5d0*log(finf(j,d))
                    rankp = rankp -1
                    do i = 1, m
                        if(pirec(i,i) .LT. meps) then
                            pirec(i,:) = 0.0d0
                            pirec(:,i) = 0.0d0
                        end if
                    end do
                else
                    if (ft(j,d) .GT. meps) then
                        call daxpy(m,vt(j)/ft(j,d),kt(:,j,d),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)/ft(i,t)
                        call dsyr('u',m,(-1.0d0)/ft(j,d),kt(:,j,d),1,prec,m) !prec = prec -kt*kt'/ft
                        lik = lik - 0.5d0*(log(ft(j,d)) + vt(j)**2/ft(j,d))
                    end if
                end if
                if(rankp .EQ. 0) then
                    exit diffuse
                end if

            end do
           
            call dgemv('n',m,m,1.0d0,tt(:,:,(d-1)*timevar(3)+1),m,arec,1,0.0d0,at,1)
            arec = at
            call dsymm('r','u',m,m,1.0d0,prec,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
            call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,pt,m)
            pt = pt + rqr(:,:,(d-1)*tv+1)
            prec = pt
            call dsymm('r','u',m,m,1.0d0,pirec,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
            call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,pinf,m)
            pirec = pinf
     
        end do diffuse
        if(rankp .EQ. 0) then
            !non-diffuse filtering begins
            prec = prec
            do i = j+1, p

                vt(i) = yt(d,i) - ddot(m,zt(i,:,(d-1)*timevar(1)+1),1,arec,1)
                call dsymv('u',m,1.0d0,prec,m,zt(i,:,(d-1)*timevar(1)+1),1,0.0d0,kt(:,i,d),1)
                ft(i,d) = ddot(m,zt(i,:,(d-1)*timevar(1)+1),1,kt(:,i,d),1)
                if (ft(i,d) .GT. meps) then !ft.NE.0
                    call daxpy(m,vt(i)/ft(i,d),kt(:,i,d),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)
                    call dsyr('u',m,-1.0d0/ft(i,d),kt(:,i,d),1,prec,m) !p_rec = p_rec - kt*kt'*ft(i,t)
                    lik = lik - 0.5d0*(log(ft(i,d)) + vt(i)**2/ft(i,d))
                end if

            end do
   
            call dgemv('n',m,m,1.0d0,tt(:,:,(d-1)*timevar(3)+1),m,arec,1,0.0d0,at,1)
            arec = at
            call dsymm('r','u',m,m,1.0d0,prec,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
            call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,pt,m)
 
            pt = pt + rqr(:,:,(d-1)*tv+1)
            prec = pt
        end if
    end if



    !Non-diffuse filtering continues from t=d+1, i=1


    if(d.EQ.0) then
        prec = p1
        arec = a1
    end if
    do t = d+1, n
        do i = 1, p
            vt(i) = yt(t,i) - ddot(m,zt(i,:,(t-1)*timevar(1)+1),1,arec,1)
            call dsymv('u',m,1.0d0,prec,m,zt(i,:,(t-1)*timevar(1)+1),1,0.0d0,kt(:,i,t),1)
            ft(i,t) = ddot(m,zt(i,:,(t-1)*timevar(1)+1),1,kt(:,i,t),1)
            if (ft(i,t) .GT. meps) then
                call daxpy(m,vt(i)/ft(i,t),kt(:,i,t),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)
                call dsyr('u',m,-1.0d0/ft(i,t),kt(:,i,t),1,prec,m) !p_rec = p_rec - kt*kt'*ft(i,i,t)
                lik = lik - 0.5d0*(log(ft(i,t)) + vt(i)**2/ft(i,t))
            end if
        end do
        call dgemv('n',m,m,1.0d0,tt(:,:,(t-1)*timevar(3)+1),m,arec,1,0.0d0,at,1)
        arec = at
        call dsymm('r','u',m,m,1.0d0,prec,m,tt(:,:,(t-1)*timevar(3)+1),m,0.0d0,mm,m)
        call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(:,:,(t-1)*timevar(3)+1),m,0.0d0,pt,m)
        pt = pt + rqr(:,:,(t-1)*tv+1)
        prec = pt
    end do


   

end subroutine pthetafirst

subroutine pthetarest(yt, timevar, zt, tt, a1,&
p, m, n, lik, tol,kt,kinf,ft,finf,dt,jt)


    implicit none

    integer, intent(in) ::  p, m, n,dt,jt
    integer ::  i,j,t,d
    integer, intent(in), dimension(5) :: timevar
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in) :: tol
    double precision, intent(inout) :: lik
    double precision, dimension(m) :: at,arec
    double precision, dimension(p) :: vt
    double precision, intent(inout), dimension(p,n) :: ft,finf
    double precision, intent(inout), dimension(m,p,n) :: kt,kinf
    double precision, external :: ddot
    double precision :: meps

    j=0
    d=0

    if(dt.GT.0) then
        arec = a1
        diffuse: do while(d .LT. (dt-1))
            d = d+1
            do j=1, p

                    vt(j) = yt(d,j) - ddot(m,zt(j,:,(d-1)*timevar(1)+1),1,arec,1) !arec
                    if (finf(j,d) .GT. tol) then
                        call daxpy(m,vt(j)/finf(j,d),kinf(:,j,d),1,arec,1) !a_rec = a_rec + kinf(:,i,t)*vt(:,t)/finf(j,d)
                        lik = lik - 0.5d0*log(finf(j,d))
                    else
                        if(ft(j,d) .GT. meps) then
                            call daxpy(m,vt(j)/ft(j,d),kt(:,j,d),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)/ft(i,t)
                            lik = lik - 0.5d0*(log(ft(j,d)) + vt(j)**2/ft(j,d))
                        end if
                    end if

            end do

        call dgemv('n',m,m,1.0d0,tt(:,:,(d-1)*timevar(3)+1),m,arec,1,0.0d0,at,1)
        arec = at

        end do diffuse

        d = dt
        do j=1, jt

                vt(j) = yt(d,j) - ddot(m,zt(j,:,(d-1)*timevar(1)+1),1,arec,1) !arec

                if (finf(j,d) .GT. tol) then
                    call daxpy(m,vt(j)/finf(j,d),kinf(:,j,d),1,arec,1) !a_rec = a_rec + kinf(:,i,t)*vt(:,t)/finf(j,d)
                    lik = lik - 0.5d0*log(finf(j,d))
                else
                    if(ft(j,d) .GT. meps ) then
                        call daxpy(m,vt(j)/ft(j,d),kt(:,j,d),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)/ft(i,t)
                        lik = lik - 0.5d0*(log(ft(j,d)) + vt(j)**2/ft(j,d))
                    end if
                end if

        end do


        !non-diffuse filtering begins

        do i = jt+1, p

                vt(i) = yt(d,i) - ddot(m,zt(i,:,(d-1)*timevar(1)+1),1,arec,1) !vt
                if (ft(i,d) .GT.  meps) then !ft.NE.0
                    call daxpy(m,vt(i)/ft(i,d),kt(:,i,d),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)
                    lik = lik - 0.5d0*(log(ft(i,d)) + vt(i)**2/ft(i,d))
                end if

        end do

        call dgemv('n',m,m,1.0d0,tt(:,:,(d-1)*timevar(3)+1),m,arec,1,0.0d0,at,1)
        arec = at
    end if

    if(dt.LT.n) then

        !Non-diffuse filtering continues from t=d+1, i=1


        if(dt.EQ.0) then
            arec = a1
        end if
        do t = dt+1, n
            do i = 1, p

                    vt(i) = yt(t,i) - ddot(m,zt(i,:,(t-1)*timevar(1)+1),1,arec,1) !variate vt
                    if (ft(i,t) .GT.  meps) then !ft.NE.0
                        call daxpy(m,vt(i)/ft(i,t),kt(:,i,t),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)
                        lik = lik - 0.5d0*(log(ft(i,t)) + vt(i)**2/ft(i,t))
                    end if

            end do

        call dgemv('n',m,m,1.0d0,tt(:,:,(t-1)*timevar(3)+1),m,arec,1,0.0d0,at,1)
        arec = at
        end do

    end if



end subroutine pthetarest

