! Subroutine for computing the log-Likelihood of general linear gaussian state space model

subroutine gloglik(yt, ymiss, timevar, zt, ht, tt, rt, qt, a1, p1, p1inf,&
p, m, r, n, lik, tol,rankp,marginal)


    implicit none

    integer, intent(in) ::  p, m, r, n
    integer, intent(inout) :: rankp,marginal
    integer ::  t, i,d,j,tv
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(p,p,(n-1)*timevar(2)+1) :: ht
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rt
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, intent(in) :: tol
    double precision, intent(inout) :: lik
    double precision, dimension(m) :: at,arec
    double precision, dimension(p) :: vt,ft,finf
    double precision, dimension(m,p) :: kt,kinf
    double precision, dimension(m,m) :: pt,pinf,mm
    double precision, dimension(m,r) :: mr    
    double precision :: c
    double precision, external :: ddot
    double precision :: meps
    double precision, dimension(m,m,(n-1)*max(timevar(4),timevar(5))+1) :: rqr

    ! compute RQR'
    tv= max(timevar(4),timevar(5))
    do t=1, (n-1)*tv+1
        call dsymm('r','l',m,r,1.0d0,qt(:,:,(t-1)*timevar(5)+1),r,rt(:,:,(t-1)*timevar(4)+1),m,0.0d0,mr,m)
        call dgemm('n','t',m,m,r,1.0d0,mr,m,rt(:,:,(t-1)*timevar(4)+1),m,0.0d0,rqr(:,:,t),m)
    end do

    ! tolerance for checking whether ft and pinf are > 0
    meps = tiny(meps)

    ! constant term for log-likelihood
    c = 0.5d0*log(8.0d0*atan(1.0d0))
    lik = 0.0d0

    j=0
    d=0
    pt = p1
    arec = a1
    ! Diffuse initialization
    if(maxval(p1inf) .GT.  0.0d0) then
        pinf=p1inf

        diffuse: do while(d .LT. n)
            d = d+1
            do j=1, p !sequential processing
                if(ymiss(d,j).EQ.0) then
                    vt(j) = yt(d,j) - ddot(m,zt(j,:,(d-1)*timevar(1)+1),1,arec,1)
                    call dsymv('u',m,1.0d0,pt,m,zt(j,:,(d-1)*timevar(1)+1),1,0.0d0,kt(:,j),1) ! kt_t,i = pt_t,i*t(z_t,i)
                    ft(j) = ddot(m,zt(j,:,(d-1)*timevar(1)+1),1,kt(:,j),1)+ ht(j,j,(d-1)*timevar(2)+1)
                    call dsymv('u',m,1.0d0,pinf,m,zt(j,:,(d-1)*timevar(1)+1),1,0.0d0,kinf(:,j),1) ! kinf_t,i = pinf_t,i*t(z_t,i)
                    finf(j) = ddot(m,zt(j,:,(d-1)*timevar(1)+1),1,kinf(:,j),1)
                    if (finf(j) .GT. tol) then
                        call daxpy(m,vt(j)/finf(j),kinf(:,j),1,arec,1) !a_rec = a_rec + kinf(:,i,t)*vt(:,t)/finf(j,d)
                        call dsyr('u',m,ft(j)/finf(j)**2,kinf(:,j),1,pt,m) !pt = pt +  kinf*kinf'*ft/finf^2
                        call dsyr2('u',m,-1.0d0/finf(j),kt(:,j),1,kinf(:,j),1,pt,m) !pt = pt -(kt*kinf'+kinf*kt')/finf
                        call dsyr('u',m,-1.0d0/finf(j),kinf(:,j),1,pinf,m) !pirec = pirec -kinf*kinf'/finf
                        lik = lik - 0.5d0*log(finf(j))
                        rankp = rankp -1

                    else
                        if (ft(j) .GT. meps) then
                            call daxpy(m,vt(j)/ft(j),kt(:,j),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)/ft(i,t)
                            call dsyr('u',m,-1.0d0/ft(j),kt(:,j),1,pt,m) !pt = pt -kt*kt'/ft
                            lik = lik - c - 0.5d0*(log(ft(j)) + vt(j)**2/ft(j))
                        end if
                    end if
                    if(rankp .EQ. 0) then
                        exit diffuse
                    end if
                end if
            end do
           
            call dgemv('n',m,m,1.0d0,tt(:,:,(d-1)*timevar(3)+1),m,arec,1,0.0d0,at,1)
            arec = at
            call dsymm('r','u',m,m,1.0d0,pt,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
            call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,pt,m)
            pt = pt + rqr(:,:,(d-1)*tv+1)
            call dsymm('r','u',m,m,1.0d0,pinf,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
            call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,pinf,m)
            do i = 1, m ! try to deal with possible rounding errors, non-diffuse states should have zeros in pinf
                if(pinf(i,i) .LT. meps) then
                    pinf(i,:) = 0.0d0
                    pinf(:,i) = 0.0d0
                end if
            end do
        end do diffuse
        if(rankp .EQ. 0) then
            !non-diffuse filtering begins
            do i = j+1, p
                if(ymiss(d,i).EQ.0) then
                    vt(i) = yt(d,i) - ddot(m,zt(i,:,(d-1)*timevar(1)+1),1,arec,1)
                    call dsymv('u',m,1.0d0,pt,m,zt(i,:,(d-1)*timevar(1)+1),1,0.0d0,kt(:,i),1)
                    ft(i) = ddot(m,zt(i,:,(d-1)*timevar(1)+1),1,kt(:,i),1) + ht(i,i,(d-1)*timevar(2)+1)
                    if (ft(i) .GT. meps) then !ft.NE.0
                        call daxpy(m,vt(i)/ft(i),kt(:,i),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)
                        call dsyr('u',m,-1.0d0/ft(i),kt(:,i),1,pt,m) !p_rec = p_rec - kt*kt'*ft(i,t)
                        lik = lik - 0.5d0*(log(ft(i)) + vt(i)**2/ft(i))-c
                    end if
                end if
            end do
   
            call dgemv('n',m,m,1.0d0,tt(:,:,(d-1)*timevar(3)+1),m,arec,1,0.0d0,at,1)
            arec = at
            call dsymm('r','u',m,m,1.0d0,pt,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
            call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(:,:,(d-1)*timevar(3)+1),m,0.0d0,pt,m)
 
            pt = pt + rqr(:,:,(d-1)*tv+1)
        end if
    end if



    !Non-diffuse filtering continues from t=d+1, i=1
    do t = d+1, n
        do i = 1, p
            if(ymiss(t,i).EQ.0) then
                vt(i) = yt(t,i) - ddot(m,zt(i,:,(t-1)*timevar(1)+1),1,arec,1)
                call dsymv('u',m,1.0d0,pt,m,zt(i,:,(t-1)*timevar(1)+1),1,0.0d0,kt(:,i),1)
                ft(i) = ddot(m,zt(i,:,(t-1)*timevar(1)+1),1,kt(:,i),1) + ht(i,i,(t-1)*timevar(2)+1)
                if (ft(i) .GT. meps) then
                    call daxpy(m,vt(i)/ft(i),kt(:,i),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)
                    call dsyr('u',m,-1.0d0/ft(i),kt(:,i),1,pt,m) !pt = pt - kt*kt'*ft(i,i,t)
                    lik = lik - 0.5d0*(log(ft(i)) + vt(i)**2/ft(i))-c
                end if
            end if
        end do
        call dgemv('n',m,m,1.0d0,tt(:,:,(t-1)*timevar(3)+1),m,arec,1,0.0d0,at,1)
        arec = at
        call dsymm('r','u',m,m,1.0d0,pt,m,tt(:,:,(t-1)*timevar(3)+1),m,0.0d0,mm,m)
        call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(:,:,(t-1)*timevar(3)+1),m,0.0d0,pt,m)
        pt = pt + rqr(:,:,(t-1)*tv+1)


    end do

    if(marginal==1) then
        t = int(sum(p1inf))
        if(t>0) then
            call marginalxx(p1inf,zt,tt,m,p,n,t,timevar,lik,marginal)
        end if
    end if
   

end subroutine gloglik
