! Subroutine for computing the log-Likelihood of univariate linear gaussian state space model

subroutine glogliku(yt, ymiss, timevar, zt, ht,tt, rt, qt, a1, p1, p1inf,&
m, r, n, lik, tol,rankp,marginal)

    implicit none

    integer, intent(in) ::  m, r, n
    integer, intent(in), dimension(n,1) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    integer, intent(inout) :: rankp,marginal
    integer ::  t, i, d
    double precision, intent(in), dimension(n,1) :: yt
    double precision, intent(in), dimension(1,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(1,1,(n-1)*timevar(2)+1) :: ht
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rt
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, intent(in) :: tol
    double precision, intent(inout) :: lik
    double precision, dimension(m) :: at,arec
    double precision  :: vt,ft,finf
    double precision, dimension(m,1) :: kt,kinf
    double precision, dimension(m,m) :: pt, pinf,mm,prec,pirec
    double precision, dimension(m,r) :: mr
    double precision :: c
    double precision, external :: ddot
    double precision :: meps

    meps = tiny(meps)
    c = 0.5d0*log(8.0d0*atan(1.0d0))
    lik = 0.0d0

    pinf=p1inf

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
            if(ymiss(d,1).EQ.0) then
                vt = yt(d,1) - ddot(m,zt(1,:,(d-1)*timevar(1)+1),1,arec,1)
                call dsymv('u',m,1.0d0,prec,m,zt(1,:,(d-1)*timevar(1)+1),1,0.0d0,kt(:,1),1) ! kt_t,i = pt_t,i*t(z_t,i)
                ft = ddot(m,zt(1,:,(d-1)*timevar(1)+1),1,kt(:,1),1)+ ht(1,1,(d-1)*timevar(2)+1)
                call dsymv('u',m,1.0d0,pirec,m,zt(1,:,(d-1)*timevar(1)+1),1,0.0d0,kinf(:,1),1) ! kinf_t,i = pinf_t,i*t(z_t,i)
                finf = ddot(m,zt(1,:,(d-1)*timevar(1)+1),1,kinf(:,1),1)
                if (finf .GT. tol) then
                    call daxpy(m,vt/finf,kinf(:,1),1,arec,1) !a_rec = a_rec + kinf(:,i,t)*vt(:,t)/finf(j,d)
                    call dsyr('u',m,ft/(finf**2),kinf(:,1),1,prec,m) !prec = prec +  kinf*kinf'*ft/finf^2
                    call dsyr2('u',m,-1.0d0/finf,kt(:,1),1,kinf(:,1),1,prec,m) !prec = prec -(kt*kinf'+kinf*kt')/finf
                    call dsyr('u',m,-1.0d0/finf,kinf(:,1),1,pirec,m) !pirec = pirec -kinf*kinf'/finf
                    lik = lik - 0.5d0*log(finf)
                    do i = 1, m
                        if(pirec(i,i) .LT. meps) then
                            pirec(i,:) = 0.0d0
                            pirec(:,i) = 0.0d0
                        end if
                    end do
                    rankp = rankp -1
                else
                    if (ft .GT. meps) then
                        call daxpy(m,vt/ft,kt(:,1),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)/ft(i,t)
                        call dsyr('u',m,(-1.0d0)/ft,kt(:,1),1,prec,m) !prec = prec -kt*kt'/ft
                        lik = lik - c - 0.5d0*(log(ft) + vt**2/ft)
                    end if
                end if
                if(rankp .EQ. 0) then
                    exit diffuse
                end if

            end if

            call dgemv('n',m,m,1.0d0,tt(:,:,(d-1)*timevar(3)+1),m,arec,1,0.0d0,at(:),1)
            call dcopy(m,at(:),1,arec,1) ! a_rec = at(:,t+1)
            call dsymm('r','u',m,m,1.0d0,prec,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
            call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,pt,m)


            call dsymm('r','u',m,r,1.0d0,qt(:,:,(d-1)*timevar(5)+1),r,rt(:,:,(d-1)*timevar(4)+1),m,0.0d0,mr,m)
            call dgemm('n','t',m,m,r,1.0d0,mr,m,rt(:,:,(d-1)*timevar(4)+1),m,1.0d0,pt,m)

            prec = pt
            call dsymm('r','u',m,m,1.0d0,pirec,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
            call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,pinf,m)
            pirec = pinf


        end do diffuse

        !non-diffuse filtering begins
        if(rankp .EQ. 0) then

            call dgemv('n',m,m,1.0d0,tt(:,:,(d-1)*timevar(3)+1),m,arec,1,0.0d0,at(:),1)  !at(:,t+1) = matmul(tt,a_rec)

            call dsymm('r','u',m,m,1.0d0,prec,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
            call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,pt,m)

            call dsymm('r','u',m,r,1.0d0,qt(:,:,(d-1)*timevar(5)+1),r,rt(:,:,(d-1)*timevar(4)+1),m,0.0d0,mr,m)
            call dgemm('n','t',m,m,r,1.0d0,mr,m,rt(:,:,(d-1)*timevar(4)+1),m,1.0d0,pt,m)


            call dcopy(m,at(:),1,arec,1)
            prec = pt
        end if
    end if

   
    !Non-diffuse filtering continues from t=d+1, i=1


    if(d.EQ.0) then
        prec = p1
        arec = a1!   call dcopy(m,a1,1,arec,1)
       !!at(:) = a1 !call dcopy(m,a1,1,at(:),1) !at(:,1) = a1
       !!pt = p1
    end if
    do t = d+1, n
        if(ymiss(t,1).EQ.0) then
            vt = yt(t,1) - ddot(m,zt(1,:,(t-1)*timevar(1)+1),1,arec,1)
            call dsymv('u',m,1.0d0,prec,m,zt(1,:,(t-1)*timevar(1)+1),1,0.0d0,kt(:,1),1)
            ft = ddot(m,zt(1,:,(t-1)*timevar(1)+1),1,kt(:,1),1)+ht(1,1,(t-1)*timevar(2)+1)
            if (ft .GT. meps) then
                call daxpy(m,vt/ft,kt(:,1),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)
                call dsyr('u',m,-1.0d0/ft,kt(:,1),1,prec,m) !p_rec = p_rec - kt*kt'*ft(i,i,t)
                lik = lik - 0.5d0*(log(ft) + vt**2/ft)-c
            end if
        end if
        call dgemv('n',m,m,1.0d0,tt(:,:,(t-1)*timevar(3)+1),m,arec,1,0.0d0,at(:),1)
        call dsymm('r','u',m,m,1.0d0,prec,m,tt(:,:,(t-1)*timevar(3)+1),m,0.0d0,mm,m)
        call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(:,:,(t-1)*timevar(3)+1),m,0.0d0,pt,m)
        call dsymm('r','u',m,r,1.0d0,qt(:,:,(t-1)*timevar(5)+1),r,rt(:,:,(t-1)*timevar(4)+1),m,0.0d0,mr,m)
        call dgemm('n','t',m,m,r,1.0d0,mr,m,rt(:,:,(t-1)*timevar(4)+1),m,1.0d0,pt,m)
        call dcopy(m,at(:),1,arec,1) ! a_rec =at(:,t+1)
        prec = pt
    end do

    if(marginal==1) then
        t = int(sum(p1inf))
        if(t>0) then
            call marginalxx(p1inf,zt,tt,m,1,n,t,timevar,lik,marginal)
        end if
    end if

end subroutine glogliku

