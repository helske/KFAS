! functions for computing p(theta)


! Only differences with gloglik is what is stored/returned (ft,finf,kt,kinf)
! Also no missing values as yt is actually the signals, and no ht
subroutine pthetafirst(yt, timevar, zt, tt, rqr, a1, p1, p1inf,&
p, m, n, lik, tol,rankp2,kt,kinf,ft,finf,d,j)


    implicit none

    integer, intent(in) ::  p, m, n
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
    double precision, dimension(m,m) :: pt,pinf,mm
    double precision, external :: ddot
    double precision, intent(inout), dimension(m,m,(n-1)*max(timevar(4),timevar(5))+1) :: rqr
    double precision :: meps, finv
    integer, dimension(p) :: ymiss
    double precision, dimension(p,p) :: ht

    external dgemm, dsymm, dgemv, dsymv, dsyr, dsyr2

    meps = epsilon(meps)
    tv= max(timevar(4),timevar(5))
    ymiss = 0
    ht = 0.0d0

    rankp = rankp2
    j=0
    d=0
    at = a1
    pt = p1
    pinf=p1inf

    ! Diffuse initialization
    if(rankp .GT. 0) then
        diffuse: do while(d .LT. n .AND. rankp .GT. 0)
            d = d+1
            call diffusefilteronestep(ymiss,yt(d,:),transpose(zt(:,:,(d-1)*timevar(1)+1)),ht,&
            tt(:,:,(d-1)*timevar(3)+1),rqr(:,:,(d-1)*tv+1),&
            at,pt,vt,ft(:,d),kt(:,:,d),pinf,finf(:,d),kinf(:,:,d),rankp,lik,tol,meps,0.0d0,p,m,j)


        end do diffuse
        if(rankp .EQ. 0 .AND. j .LT. p) then
            call filteronestep(ymiss,yt(d,:),transpose(zt(:,:,(d-1)*timevar(1)+1)),ht,&
            tt(:,:,(d-1)*timevar(3)+1),rqr(:,:,(d-1)*tv+1),&
            at,pt,vt,ft(:,d),kt(:,:,d),lik,tol,meps,0.0d0,p,m,j)

        else
            j = p
        end if

    end if



    !Non-diffuse filtering continues from t=d+1, i=1

    do t = d+1, n
        call filteronestep(ymiss,yt(t,:),transpose(zt(:,:,(t-1)*timevar(1)+1)),ht,&
        tt(:,:,(t-1)*timevar(3)+1),rqr(:,:,(t-1)*tv+1),&
        at,pt,vt,ft(:,t),kt(:,:,t),lik,tol,meps,0.0d0,p,m,0)
    end do


   

end subroutine pthetafirst


! use output of pthetafirst, kt and ft do not change
subroutine pthetarest(yt, timevar, zt, tt, a1,&
p, m, n, lik, kt,kinf,ft,finf,dt,jt)


    implicit none

    integer, intent(in) ::  p, m, n,dt,jt
    integer ::  i,j,t,d
    integer, intent(in), dimension(5) :: timevar
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(inout) :: lik
    double precision, dimension(m) :: at,arec
    double precision, dimension(p) :: vt
    double precision, intent(in), dimension(p,n) :: ft,finf
    double precision, intent(in), dimension(m,p,n) :: kt,kinf
    double precision, external :: ddot

    external dgemv

    j=0
    d=0
    arec = a1

    if(dt.GT.0) then
        diffuse: do while(d .LT. (dt-1))
            d = d+1
            do j=1, p

                vt(j) = yt(d,j) - ddot(m,zt(j,:,(d-1)*timevar(1)+1),1,arec,1)
                if (finf(j,d) .GT. 0.0d0) then
                    arec = arec + vt(j)/finf(j,d)*kinf(:,j,d)

                    lik = lik - 0.5d0*log(finf(j,d))
                else
                    if(ft(j,d) .GT. 0.0d0) then
                        arec = arec + vt(j)/ft(j,d)*kt(:,j,d)

                        lik = lik - 0.5d0*(log(ft(j,d)) + vt(j)**2/ft(j,d))
                    end if
                end if

            end do

            call dgemv('n',m,m,1.0d0,tt(:,:,(d-1)*timevar(3)+1),m,arec,1,0.0d0,at,1)
            arec = at

        end do diffuse

        d = dt
        do j=1, jt

            vt(j) = yt(d,j) - ddot(m,zt(j,:,(d-1)*timevar(1)+1),1,arec,1)

            if (finf(j,d) .GT. 0.0d0) then
                arec = arec + vt(j)/finf(j,d)*kinf(:,j,d)
                lik = lik - 0.5d0*log(finf(j,d))
            else
                if(ft(j,d) .GT. 0.0d0) then
                    arec = arec + vt(j)/ft(j,d)*kt(:,j,d)
                    lik = lik - 0.5d0*(log(ft(j,d)) + vt(j)**2/ft(j,d))
                end if
            end if

        end do


        !non-diffuse filtering begins

        do i = jt+1, p

            vt(i) = yt(d,i) - ddot(m,zt(i,:,(d-1)*timevar(1)+1),1,arec,1)
            if (ft(i,d) .GT. 0.0d0) then !ft.NE.0
                arec = arec + vt(i)/ft(i,d)*kt(:,i,d)

                lik = lik - 0.5d0*(log(ft(i,d)) + vt(i)**2/ft(i,d))
            end if

        end do

        call dgemv('n',m,m,1.0d0,tt(:,:,(d-1)*timevar(3)+1),m,arec,1,0.0d0,at,1)
        arec = at
    end if

    if(dt.LT.n) then

        !Non-diffuse filtering continues from t=d+1, i=1
        do t = dt+1, n
            do i = 1, p

                vt(i) = yt(t,i) - ddot(m,zt(i,:,(t-1)*timevar(1)+1),1,arec,1)
                if (ft(i,t) .GT. 0.0d0) then
                    arec = arec + vt(i)/ft(i,t)*kt(:,i,t)
                    lik = lik - 0.5d0*(log(ft(i,t)) + vt(i)**2/ft(i,t))
                end if

            end do

            call dgemv('n',m,m,1.0d0,tt(:,:,(t-1)*timevar(3)+1),m,arec,1,0.0d0,at,1)
            arec = at
        end do

    end if



end subroutine pthetarest

