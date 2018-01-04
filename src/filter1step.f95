
!diffuse filtering for single time point
subroutine dfilter1step(ymiss, yt, zt, ht, tt, rqr, at, pt, vt, ft,kt,&
pinf,finf,kinf,rankp,lik,basetol,c,p,m,i)

    implicit none

    integer, intent(in) ::  p, m
    integer, intent(inout) ::  i,rankp
    integer, intent(in), dimension(p) :: ymiss
    double precision, intent(in), dimension(p) :: yt
    double precision, intent(in), dimension(m,p) :: zt
    double precision, intent(in), dimension(p,p) :: ht
    double precision, intent(in), dimension(m,m) :: tt
    double precision, dimension(m,m) :: rqr
    double precision, intent(in) :: basetol,c
    double precision, intent(inout) :: lik
    double precision, intent(inout), dimension(m) :: at
    double precision, intent(inout), dimension(p) :: vt,ft,finf
    double precision, intent(inout), dimension(m,p) :: kt,kinf
    double precision, intent(inout), dimension(m,m) :: pt,pinf
    double precision, dimension(m,m) :: mm
    double precision, dimension(m) :: ahelp
    double precision :: finv, tol

    double precision, external :: ddot

    external dgemm, dsymm, dgemv, dsymv, dsyr, dsyr2
    tol = basetol * minval(abs(zt), mask = abs(zt) .GT. 0.0d0)**2
    
    do i=1, p
        call dsymv('u',m,1.0d0,pt,m,zt(:,i),1,0.0d0,kt(:,i),1)
        ft(i) = ddot(m,zt(:,i),1,kt(:,i),1) + ht(i,i)
        if(ymiss(i) .EQ. 0) then
            call dsymv('u',m,1.0d0,pinf,m,zt(:,i),1,0.0d0,kinf(:,i),1)
            finf(i) = ddot(m,zt(:,i),1,kinf(:,i),1)
            vt(i) = yt(i) - ddot(m,zt(:,i),1,at,1)
            if (finf(i) .GT. tol) then
                finv = 1.0d0/finf(i)
                at = at + vt(i)*finv*kinf(:,i)
                call dsyr('u',m,ft(i)*finv**2,kinf(:,i),1,pt,m)
                call dsyr2('u',m,-finv,kt(:,i),1,kinf(:,i),1,pt,m)
                call dsyr('u',m,-finv,kinf(:,i),1,pinf,m)
                lik = lik - 0.5d0*log(finf(i))
                rankp = rankp -1
            else
                finf(i) = 0.0d0
                if(ft(i) .GT. tol) then
                    finv = 1.0d0/ft(i)
                    at = at + vt(i)*finv*kt(:,i)
                    call dsyr('u',m,-finv,kt(:,i),1,pt,m)
                    lik = lik - c - 0.5d0*(log(ft(i)) + vt(i)**2*finv)
                end if
            end if
            if (ft(i) .LE. tol) then
                ft(i) = 0.0d0
            end if
            if(rankp .EQ. 0 .AND. i .LT. p) then
                return
            end if
        end if
    end do

    call dgemv('n',m,m,1.0d0,tt,m,at,1,0.0d0,ahelp,1)
    at = ahelp
    call dsymm('r','u',m,m,1.0d0,pt,m,tt,m,0.0d0,mm,m)
    call dgemm('n','t',m,m,m,1.0d0,mm,m,tt,m,0.0d0,pt,m)
    pt = pt + rqr

    call dsymm('r','u',m,m,1.0d0,pinf,m,tt,m,0.0d0,mm,m)
    call dgemm('n','t',m,m,m,1.0d0,mm,m,tt,m,0.0d0,pinf,m)

end subroutine dfilter1step

!non-diffuse filtering for single time point
subroutine filter1step(ymiss, yt, zt, ht, tt, rqr, at, pt, vt, & 
    ft, kt, lik, basetol, c, p, m, j)

    implicit none

    integer, intent(in) ::  p, m,j
    integer ::  i
    integer, intent(in), dimension(p) :: ymiss
    double precision, intent(in), dimension(p) :: yt
    double precision, intent(in), dimension(m,p) :: zt
    double precision, intent(in), dimension(p,p) :: ht
    double precision, intent(in), dimension(m,m) :: tt
    double precision, dimension(m,m) :: rqr
    double precision, intent(in) :: basetol,c
    double precision, intent(inout) :: lik
    double precision, intent(inout), dimension(m) :: at
    double precision, intent(inout), dimension(p) :: vt,ft
    double precision, intent(inout), dimension(m,p) :: kt
    double precision, intent(inout), dimension(m,m) :: pt
    double precision, dimension(m,m) :: mm
    double precision, dimension(m) :: ahelp
    double precision :: finv, tol

    double precision, external :: ddot

    external dgemm, dsymm, dgemv, dsymv, dsyr
    tol = basetol * minval(abs(zt), mask = abs(zt) .GT. 0.0d0)**2
    
    do i = j+1, p
        call dsymv('u',m,1.0d0,pt,m,zt(:,i),1,0.0d0,kt(:,i),1)
        ft(i) = ddot(m,zt(:,i),1,kt(:,i),1) + ht(i,i)
        if(ymiss(i).EQ.0) then
            vt(i) = yt(i) - ddot(m,zt(:,i),1,at,1)
            if (ft(i) .GT. tol) then
                finv = 1.0d0/ft(i)
                at = at + vt(i)*finv*kt(:,i)
                call dsyr('u',m,-finv,kt(:,i),1,pt,m)
                lik = lik - c - 0.5d0*(log(ft(i)) + vt(i)**2*finv)
            else
                ft(i) = 0.0d0
            end if
        end if
    end do

    call dgemv('n',m,m,1.0d0,tt,m,at,1,0.0d0,ahelp,1)
    at = ahelp

    call dsymm('r','u',m,m,1.0d0,pt,m,tt,m,0.0d0,mm,m)
    call dgemm('n','t',m,m,m,1.0d0,mm,m,tt,m,0.0d0,pt,m)
    pt = pt + rqr

end subroutine filter1step

!! as above, but returns also att and ptt


!diffuse filtering for single time point
subroutine dfilter1step2(ymiss, yt, zt, ht, tt, rqr, at, pt, vt, ft,kt,&
pinf,finf,kinf,rankp,lik,basetol,c,p,m,i, att, ptt)

    implicit none

    integer, intent(in) ::  p, m
    integer, intent(inout) ::  i,rankp
    integer, intent(in), dimension(p) :: ymiss
    double precision, intent(in), dimension(p) :: yt
    double precision, intent(in), dimension(m,p) :: zt
    double precision, intent(in), dimension(p,p) :: ht
    double precision, intent(in), dimension(m,m) :: tt
    double precision, dimension(m,m) :: rqr
    double precision, intent(in) :: basetol,c
    double precision, intent(inout) :: lik
    double precision, intent(inout), dimension(m) :: at, att
    double precision, intent(inout), dimension(p) :: vt,ft,finf
    double precision, intent(inout), dimension(m,p) :: kt,kinf
    double precision, intent(inout), dimension(m,m) :: pt,pinf, ptt
    double precision, dimension(m,m) :: mm
    double precision, dimension(m) :: ahelp
    double precision :: finv, tol

    double precision, external :: ddot

    external dgemm, dsymm, dgemv, dsymv, dsyr, dsyr2
    tol = basetol * minval(abs(zt), mask = abs(zt) .GT. 0.0d0)**2
    
    do i=1, p
        call dsymv('u',m,1.0d0,pt,m,zt(:,i),1,0.0d0,kt(:,i),1)
        ft(i) = ddot(m,zt(:,i),1,kt(:,i),1) + ht(i,i)
        if(ymiss(i) .EQ. 0) then
            call dsymv('u',m,1.0d0,pinf,m,zt(:,i),1,0.0d0,kinf(:,i),1)
            finf(i) = ddot(m,zt(:,i),1,kinf(:,i),1)
            vt(i) = yt(i) - ddot(m,zt(:,i),1,at,1)
            if (finf(i) .GT. tol) then
                finv = 1.0d0/finf(i)
                at = at + vt(i)*finv*kinf(:,i)
                call dsyr('u',m,ft(i)*finv**2,kinf(:,i),1,pt,m)
                call dsyr2('u',m,-finv,kt(:,i),1,kinf(:,i),1,pt,m)
                call dsyr('u',m,-finv,kinf(:,i),1,pinf,m)
                lik = lik - 0.5d0*log(finf(i))
                rankp = rankp -1
            else
                finf(i) = 0.0d0
                if(ft(i) .GT. tol) then
                    finv = 1.0d0/ft(i)
                    at = at + vt(i)*finv*kt(:,i)
                    call dsyr('u',m,-finv,kt(:,i),1,pt,m)
                    lik = lik - c - 0.5d0*(log(ft(i)) + vt(i)**2*finv)
                end if
            end if
            if (ft(i) .LE. tol) then
                ft(i) = 0.0d0
            end if
            if(rankp .EQ. 0 .AND. i .LT. p) then
                return
            end if
        end if
    end do
    att = at
    ptt = (pt + transpose(ptt)) / 2.0d0
    call dgemv('n',m,m,1.0d0,tt,m,at,1,0.0d0,ahelp,1)
    at = ahelp
    call dsymm('r','u',m,m,1.0d0,pt,m,tt,m,0.0d0,mm,m)
    call dgemm('n','t',m,m,m,1.0d0,mm,m,tt,m,0.0d0,pt,m)
    pt = pt + rqr

    call dsymm('r','u',m,m,1.0d0,pinf,m,tt,m,0.0d0,mm,m)
    call dgemm('n','t',m,m,m,1.0d0,mm,m,tt,m,0.0d0,pinf,m)

end subroutine dfilter1step2

!non-diffuse filtering for single time point
subroutine filter1step2(ymiss, yt, zt, ht, tt, rqr, at, pt, vt, & 
    ft, kt, lik, basetol, c, p, m, j, att, ptt)

    implicit none

    integer, intent(in) ::  p, m,j
    integer ::  i
    integer, intent(in), dimension(p) :: ymiss
    double precision, intent(in), dimension(p) :: yt
    double precision, intent(in), dimension(m,p) :: zt
    double precision, intent(in), dimension(p,p) :: ht
    double precision, intent(in), dimension(m,m) :: tt
    double precision, dimension(m,m) :: rqr
    double precision, intent(in) :: basetol,c
    double precision, intent(inout) :: lik
    double precision, intent(inout), dimension(m) :: at, att
    double precision, intent(inout), dimension(p) :: vt,ft
    double precision, intent(inout), dimension(m,p) :: kt
    double precision, intent(inout), dimension(m,m) :: pt, ptt
    double precision, dimension(m,m) :: mm
    double precision, dimension(m) :: ahelp
    double precision :: finv, tol

    double precision, external :: ddot

    external dgemm, dsymm, dgemv, dsymv, dsyr
    tol = basetol * minval(abs(zt), mask = abs(zt) .GT. 0.0d0)**2
    
    do i = j+1, p
        call dsymv('u',m,1.0d0,pt,m,zt(:,i),1,0.0d0,kt(:,i),1)
        ft(i) = ddot(m,zt(:,i),1,kt(:,i),1) + ht(i,i)
        if(ymiss(i).EQ.0) then
            vt(i) = yt(i) - ddot(m,zt(:,i),1,at,1)
            if (ft(i) .GT. tol) then
                finv = 1.0d0/ft(i)
                at = at + vt(i)*finv*kt(:,i)
                call dsyr('u',m,-finv,kt(:,i),1,pt,m)
                lik = lik - c - 0.5d0*(log(ft(i)) + vt(i)**2*finv)
            else
                ft(i) = 0.0d0
            end if
        end if
    end do

    att = at
    ptt = (pt + transpose(ptt)) / 2.0d0
    call dgemv('n',m,m,1.0d0,tt,m,at,1,0.0d0,ahelp,1)
    at = ahelp

    call dsymm('r','u',m,m,1.0d0,pt,m,tt,m,0.0d0,mm,m)
    call dgemm('n','t',m,m,m,1.0d0,mm,m,tt,m,0.0d0,pt,m)
    pt = pt + rqr

end subroutine filter1step2

