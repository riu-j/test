    module sreft
      use rocm_real, only: omega => varnf
      use nmprd_int, only: nthes_ => nwtht, netas_ => nweta, nepss_ => nweps
      use sizes, only: dpsize, isize
      implicit none
      real (kind=dpsize) :: generalmeanage = 70.0d0, generaltriggerage = 60.0d0, plotmax = 25.0d0, &
        plotmin = -5.0d0, sum, wsum, sigma(3)
      real (kind=dpsize), save :: tc(3), offt(400), ww, ma, mi, err, offsett, m
      integer (kind=isize) :: numbm = 3, numsj = 400
    contains

      subroutine initialize(id, meanx, meany, coun, a, b, c)
        real (kind=dpsize) :: meanx(*), meany(*), coun(*), a(*), b(*), c(*)
        integer (kind=isize) :: id, i, intc

        sum = 0.0d0
        wsum = 0.0d0
        do i = 1, numbm
          sigma(i) = omega(netas_+i, netas_+i)
        end do
        do i = 1, numbm
          ma = plotmax - meanx(i)
          mi = plotmin + meanx(i)
          intc = int(c(i))
          if (b(i)==0.0) cycle
          select case (intc)
          case (0)
            m = (meany(i)-a(i))/b(i)
          case default
            if ((meany(i)-a(i))*c(i)/b(i)+1.0<0) then
              if (c(i)>0) then
                m = mi
              else
                m = ma
              end if
            else
              m = log((meany(i)-a(i))*c(i)/b(i)+1.0)/c(i)
            end if
          end select
          if (m>0.25) then
            if (c(i)>0) then
              tc(i) = mi
            else
              tc(i) = ma
            end if
          else
            tc(i) = (1-sqrt(1-4.0d0*m))/2.0d0
          end if
          if (tc(i)<mi) then
            tc(i) = mi
          else if (tc(i)>ma) then
            tc(i) = ma
          else
            tc(i) = tc(i)
          end if
          m = 1 - (tc(i)-meanx(i))
          b(i) = b(i)*m
          c(i) = c(i)*m
          err = sigma(i)/coun(i)
          if (c(i)/=0.0) then
            if (i==3) then
              b(i) = b(i) - c(i)
            end if
            ww = err/((b(i)*exp(c(i)*tc(i)))*(b(i)*exp(c(i)*tc(i))))
          else
            if (i==3) then
              b(i) = b(i) + m
            end if
            ww = err/(b(i)*b(i))
          end if
          ww = 1/ww
          m = tc(i) - meanx(i)
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

      subroutine eval(a, b, c, f, g, datrec, bm, id)
        use prdims, only: gprd
        real (kind=dpsize) :: datrec(*), f, g(gprd, 1), a(*), b(*), c(*), time
        integer (kind=isize) :: bm, id
        time = datrec(3) + offt(id)
        if (c(bm)==0) then
          f = a(bm) + b(bm)*time
          g = 0.0d0
          g(bm, 1) = 1.0d0
          g(bm+numbm, 1) = time
          g(bm+2*numbm, 1) = 0.0d0
        else
          f = a(bm) + b(bm)/c(bm)*(exp(c(bm)*time)-1.0d0)
          g = 0.0d0
          g(bm, 1) = 1.0d0
          g(bm+numbm, 1) = (exp(c(bm)*time)-1.0d0)/c(bm)
          g(bm+2*numbm, 1) = b(bm)*(exp(c(bm)*time)-1.0d0)*(time-1.0d0/c(bm)/c(bm))
        end if
        return
      end subroutine

      subroutine makecsv
        implicit none
        integer (kind=isize) :: i
        open (10, file='offsetT.csv')
        write (10, '(a2, a2, a7)') 'id', ', ', 'offsetT'
        do i = 1, numsj
          write (10, '(i3, a2, e12.4)') i, ', ', offt(i)
        end do
        close (10)
        return
      end subroutine

    end module

    subroutine pred(icall, newind, theta, datrec, indxs, f, g, h)
      use sreft
      use prdims, only: gprd, hprd
      use nmprd_real, only: eta, eps
      use nmprd_int, only: newl2, iquit
      implicit none
      real (kind=dpsize) :: theta(*), datrec(*), f, g(gprd, 1), h(hprd, 1), a(3), b(3), c(3), meanx(3), &
        meany(3), coun(3)
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
          eps = 0.d0
        end if
      end if
      id = datrec(1)
      bm = datrec(2)
      meanx = datrec(5:4+numbm)
      meany = datrec(5+numbm:4+2*numbm)
      coun = datrec(5+2*numbm:4+2*numbm)
      a = theta(1:numbm) + eta(1:numbm)
      b = theta(1+numbm:2*numbm) + eta(1+numbm:2*numbm)
      c = theta(1+2*numbm:3*numbm) + eta(1+2*numbm:3*numbm)
      if (newind<=1) call initialize(id, meanx, meany, coun, a, b, c)
      call eval(a, b, c, f, g, datrec, bm, id)
      f = f*exp(eps(bm))
      h = 0.0d0
      h(bm, 1) = f
      if (icall==3) call makecsv
      return
    end subroutine
