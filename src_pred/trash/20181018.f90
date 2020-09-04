    module variables
      implicit none
      real (kind=selected_real_kind(8)) :: epss(1:3) = (/ 0.04d0, 0.0025d0, 0.01d0 /), generalmeanage = 70.0d0, &
        generaltriggerage = 60.0d0, plotmax = 25.0d0, plotmin = -5.0d0, sum, wsum
      integer (kind=selected_int_kind(8)) :: numbm = 3, numsj = 400

      real (kind=selected_real_kind(8)), save :: tc(3), offt(400), ww, ma, mi, err, offsett, m
    contains

      subroutine initialize(id, meanx, meany, coun, a, b, c)
        real (kind=selected_real_kind(8)) :: meanx(3), meany(3), coun(3), a(3), b(3), c(3)
        integer (kind=selected_int_kind(8)) :: id, i, intc

        sum = 0.0d0
        wsum = 0.0d0
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
          err = epss(i)/coun(i)
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
        real (kind=selected_real_kind(8)) :: datrec(*), f, g(gprd, 1), a(3), b(3), c(3), time
        integer (kind=selected_int_kind(8)) :: bm, id

        time = datrec(3) + offt(id)

        if (c(bm)==0) then
          f = a(bm) + b(bm)*time
          g = 0.0d0
          g(bm, 1) = 1.0d0
          g(bm+3, 1) = time
          g(bm+6, 1) = 0.0d0
        else
          f = a(bm) + b(bm)/c(bm)*(exp(c(bm)*time)-1.0d0)
          g = 0.0d0
          g(bm, 1) = 1.0d0
          g(bm+3, 1) = (exp(c(bm)*time)-1.0d0)/c(bm)
          g(bm+6, 1) = b(bm)*(exp(c(bm)*time)-1.0d0)*(time-1.0d0/c(bm)/c(bm))
        end if
        return
      end subroutine

      subroutine makecsv
        implicit none
        integer (kind=selected_int_kind(8)) :: i

        open (10, file='offsetT.csv')
        write (10, '(a2, a2, a7)') 'ID', ', ', 'offsetT'
        do i = 1, numsj
          write (10, '(i3, a2, e12.4)') i, ', ', offt(i)
        end do
        close (10)
        return
      end subroutine

    end module

    subroutine pred(icall, newind, theta, datrec, indxs, f, g, h)
      use variables
      use prdims, only: gprd, hprd
      use nmprd_real, only: eta, eps
      implicit none
      real (kind=selected_real_kind(8)) :: theta(*), datrec(*), f, g(gprd, 1), h(hprd, 1), a(3), b(3), c(3), meanx(3), &
        meany(3), coun(3)
      integer (kind=selected_int_kind(8)) :: icall, newind, indxs(*), id, bm

      id = datrec(1)
      bm = datrec(2)
      meanx(1) = datrec(5)
      meanx(2) = datrec(6)
      meanx(3) = datrec(7)
      meany(1) = datrec(8)
      meany(2) = datrec(9)
      meany(3) = datrec(10)
      coun(1) = datrec(11)
      coun(2) = datrec(12)
      coun(3) = datrec(13)
      a(1) = theta(1) + eta(1)
      a(2) = theta(2) + eta(2)
      a(3) = theta(3) + eta(3)
      b(1) = theta(4) + eta(4)
      b(2) = theta(5) + eta(5)
      b(3) = theta(6) + eta(6)
      c(1) = theta(7) + eta(7)
      c(2) = theta(8) + eta(8)
      c(3) = theta(9) + eta(9)

      if (newind<=1) call initialize(id, meanx, meany, coun, a, b, c)

      call eval(a, b, c, f, g, datrec, bm, id)

      f = f*exp(eps(bm))
      h = 0.0d0
      h(bm, 1) = f

      if (icall==3) call makecsv

      return
    end subroutine
