    module sreft
      use rocm_real, only: omega => varnf
      use nmprd_int, only: nthes_ => nwtht, netas_ => nweta, nepss_ => nweps
      use sizes, only: dpsize, isize
      implicit none
      integer (kind=isize) :: numbm = 3, numsj = 400
      real (kind=dpsize) :: plotmax = 35.0d0, plotmin = -15.0d0, sum, wsum, ww, ma, mi, err, offsett, m, tc, sigma
      real (kind=dpsize), save :: offt(400)
    contains

      subroutine initialize(id, meanx, meany, coun, a, b, c)
        real (kind=dpsize) :: meanx(*), meany(*), coun(*), a(*), b(*), c(*)
        integer (kind=isize) :: id, i

        sum = 0.0d0
        wsum = 0.0d0
        do i = 1, numbm
          sigma = omega(netas_+i, netas_+i)
          ma = plotmax - meanx(i)
          mi = plotmin + meanx(i)
          if (b(i)==0.0d0) cycle
          if (c(i)==0.0d0) then
            m = (meany(i)-a(i))/b(i)
          else if ((meany(i)-a(i))*c(i)/b(i)+1.0d0<0.0d0) then
            if (c(i)>0.0d0) then
              m = mi
            else
              m = ma
            end if
          else
            m = log((meany(i)-a(i))*c(i)/b(i)+1.0)/c(i)
          end if
          tc = m
          if (tc<mi) then
            tc = mi
          else if (tc>ma) then
            tc = ma
          else
            tc = tc
          end if
          err = sigma / coun(i)
          if (c(i)/=0.0) then
            ww = err / (b(i)*exp(c(i)*tc))**2
          else
            ww = err / b(i)**2
          end if
          ww = 1/ww
          m = tc - meanx(i)
          sum = sum + ww*m
          wsum = wsum + ww
        end do
        if (wsum>0.0d0) then
          offsett = sum/wsum
        else
          offsett = 0.0d0
        end if
        offt(id) = offsett
        return
      end subroutine

      subroutine eval(a, b, c, f, g, datrec, bm, id, meany)
        use prdims, only: gprd
        real (kind=dpsize) :: datrec(*), f, g(gprd, 1), a(*), b(*), c(*), time, meany(*)
        integer (kind=isize) :: bm, id, icall

        time = datrec(2) + offt(id)
        if (c(bm)==0) then
          f = a(bm) + b(bm)*time
          g = 0.0d0
          g(bm, 1) = -1.0d0
          g(bm+numbm, 1) = (a(bm)-meany(bm))/b(bm)
        else
          f = a(bm) + b(bm)/c(bm)*(exp(c(bm)*time)-1.0d0)
          g = 0.0d0
!          if (c(bm)*(meany(bm)-a(bm))+b(bm)==0.0d0) then
            g(bm, 1) = 1.0d0
!          else
!            g(bm, 1) = -b(bm)*exp(c(bm)*time)/(c(bm)*(meany(bm)-a(bm))+b(bm))
!          end if
!
!          if (c(bm)*(meany(bm)-a(bm))+b(bm)==0.0d0) then
            g(bm+numbm, 1) = (exp(c(bm)*time)-1.0d0)/c(bm)
!          else
!            g(bm+numbm, 1) = exp(c(bm)*time)*(a(bm)-meany(bm))/(c(bm)*(meany(bm)-a(bm))+b(bm))
!          end if
!
!          if (c(bm)*(meany(bm)-a(bm))/b(bm)+1.0d0<0.0d0 .or. c(bm)*(meany(bm)-a(bm))+b(bm)==0.0d0) then
            g(bm+2*numbm, 1) = b(bm)*(-exp(c(bm)*time)+1+c(bm)*time*exp(c(bm)*time))/(c(bm)**2)
!          else
!            g(bm+2*numbm, 1) = b(bm)*exp(c(bm)*time)/(c(bm)**2)*(-log(c(bm)*(meany(bm)- &
!              a(bm))/b(bm)+1)+c(bm)*(meany(bm)-a(bm))/(c(bm)*(meany(bm)-a(bm))+b(bm)))
!          end if
        end if
        return
      end subroutine

      subroutine makecsv
        implicit none
        integer (kind=isize) :: i

        open (13, file='offsetT.csv')
        write (13, '(a2, a2, a7)') 'id', ', ', 'offsetT'
        do i = 1, numsj
          write (13, '(i3, a2, e12.4)') i, ', ', offt(i)
        end do
        close (13)
        return
      end subroutine

    end module

    subroutine pred(icall, newind, theta, datrec, indxs, f, g, h)
      use sreft
      use prdims, only: gprd, hprd
      use nmprd_real, only: eta, eps
      use nmprd_int, only: newl2, iquit
      implicit none
      real (kind=dpsize) :: theta(*), datrec(*), f, g(gprd, 1), h(hprd, 1)
      real (kind=dpsize), allocatable :: a(:), b(:), c(:), meanx(:), meany(:), coun(:)
      integer (kind=isize) :: icall, newind, indxs(*), id, bm

      if (icall==4) then
        if (newind/=2) then
          call simeta(eta)
          if (iquit==1) return
        end if
        if (newl2==1) then
          call simeps(eps)
          if (iquit==1) return
        end if
      else
        if (newind/=2) then
          call geteta(eta)
          if (iquit==1) return
          eps = 0.0d0
        end if
      end if
      allocate (a(numbm), b(numbm), c(numbm), meanx(numbm), meany(numbm), coun(numbm))
      id = datrec(5)
      bm = datrec(4)
      meanx = datrec(6:4+numbm)
      meany = datrec(6+numbm:5+2*numbm)
      coun = datrec(6+2*numbm:5+3*numbm)
      a = theta(1:numbm) + eta(1:numbm)
      b = theta(1+numbm:2*numbm) + eta(1+numbm:2*numbm)
      c = theta(1+2*numbm:3*numbm) + eta(1+2*numbm:3*numbm)
      if (newind<=1) call initialize(id, meanx, meany, coun, a, b, c)
      call eval(a, b, c, f, g, datrec, bm, id, meany)
      f = f + eps(bm)
      h = 0.0d0
      h(bm, 1) = 1.0d0
      if (icall==3) call makecsv
      deallocate (a, b, c, meanx, meany, coun)
      return
    end subroutine
