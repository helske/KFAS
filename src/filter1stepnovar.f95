!diffuse filtering for single time point
subroutine dfilter1stepnv(ymiss, yt, zt, tt, at, vt, ft,kt,&
finf,kinf,p,m,j,lik)

    implicit none

    integer, intent(in) ::  p, m, j
    integer :: i
    integer, intent(in), dimension(p) :: ymiss
    double precision, intent(in), dimension(p) :: yt
    double precision, intent(in), dimension(m,p) :: zt
    double precision, intent(in), dimension(m,m) :: tt
    double precision, intent(inout), dimension(m) :: at
    double precision, intent(inout), dimension(p) :: vt
    double precision, intent(in), dimension(p) :: ft,finf
    double precision, intent(in), dimension(m,p) :: kt,kinf
    double precision, intent(inout) :: lik
    double precision, dimension(m) :: ahelp

    double precision, external :: ddot  
    external dgemv

    do i = 1, j
        if(ymiss(i).EQ.0) then
            vt(i) = yt(i) - ddot(m,zt(:,i),1,at,1)
            if (finf(i) .GT. 0.0d0) then
                at = at + vt(i)/finf(i)*kinf(:,i)
                lik = lik - 0.5d0*log(finf(i))
            else
                if (ft(i) .GT. 0.0d0) then
                    at = at + vt(i)/ft(i)*kt(:,i)
                    lik = lik - 0.5d0*(log(ft(i)) + vt(i)**2/ft(i))
                end if
            end if
        end if
    end do
    if(j .EQ. p) then
        call dgemv('n',m,m,1.0d0,tt,m,at,1,0.0d0,ahelp,1)
        at = ahelp
    end if
end subroutine dfilter1stepnv


!non-diffuse filtering for single time point
subroutine filter1stepnv(ymiss, yt, zt, tt, at,vt, ft,kt,p,m,j,lik)

    implicit none

    integer, intent(in) ::  p, m,j
    integer ::  i
    integer, intent(in), dimension(p) :: ymiss
    double precision, intent(in), dimension(p) :: yt
    double precision, intent(in), dimension(m,p) :: zt
    double precision, intent(in), dimension(m,m) :: tt
    double precision, intent(inout), dimension(m) :: at
    double precision, intent(in), dimension(p) :: ft
    double precision, intent(inout), dimension(p) :: vt
    double precision, intent(in), dimension(m,p) :: kt
    double precision, intent(inout) :: lik
    double precision, dimension(m) :: ahelp

    double precision, external :: ddot

    external dgemv

    do i = j+1, p
        if(ymiss(i).EQ.0) then
            vt(i) = yt(i) - ddot(m,zt(:,i),1,at,1)
            if (ft(i) .GT. 0.0d0) then
                at = at + vt(i)/ft(i)*kt(:,i)
                lik = lik - 0.5d0*(log(ft(i)) + vt(i)**2/ft(i))
            end if
        end if
    end do

    call dgemv('n',m,m,1.0d0,tt,m,at,1,0.0d0,ahelp,1)
    at = ahelp

end subroutine filter1stepnv
