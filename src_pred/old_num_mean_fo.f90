    module sreft
      use :: rocm_real, only: omega => varnf
      use :: nmprd_int, only: nthes_ => nwtht, netas_ => nweta, nepss_ => nweps
      use :: sizes, only: dpsize, isize
      use :: rocm_int, only: nindr => nindobs, indr1 => idxobsf, indr2 => idxobsl
      implicit none
      integer (kind=isize), save :: numbm, numsj
      integer (kind=isize), allocatable, save :: realid(:)
      real (kind=dpsize) :: plotmax = 40.0d0, plotmin = -15.0d0
      real (kind=dpsize), allocatable, save :: offt(:), offt2(:)
    contains

      subroutine initialize(id, meanx, meany, coun, a, b, c, bm, datrec, covt, covy, flag)
        real (kind=dpsize), intent(in) :: meanx(*), meany(*), coun(*), a(*), b(*), c(*), datrec(*), covt, covy(*)
        real (kind=dpsize) :: sum, wsum, ww, ma, mi, err, offsett, m, sigma, omegas(3), p, q, r, mm
        integer (kind=isize), intent(in) :: id, bm, flag
        integer (kind=isize) :: i

        sum = 0.0d0
        wsum = 0.0d0
        do i = 1, numbm
          if (coun(i)==0 .or. b(i)==0.0d0) cycle
          sigma = omega(netas_+i, netas_+i)
          omegas(1) = omega(i, i)
          omegas(2) = omega(i+numbm, i+numbm)
          omegas(3) = omega(i+2*numbm, i+2*numbm)
          ma = plotmax - meanx(i)
          mi = plotmin + meanx(i)
          if (bm==1) then
!main biomarker
            if (c(i)==0.0d0) then
              m = (meany(i)-a(i))/(b(i)+covy(i))/covt
            else if ((meany(i)-a(i))/(b(i)/c(i)+covy(i))+1.0<0.0d0) then
              if (c(i)>0.0d0) then
                m = mi
              else
                m = ma
              end if
            else
              m = log((meany(i)-a(i))/(b(i)/c(i)+covy(i))+1.0)/c(i)/covt
            end if
          else
!other biomarker
            if (c(i)==0.0d0) then
              m = (meany(i)-a(i)-covy(i))/b(i)/covt
            else if ((meany(i)-a(i)-covy(i))*c(i)/b(i)+1.0<0.0d0) then
              if (c(i)>0.0d0) then
                m = mi
              else
                m = ma
              end if
            else
              m = log((meany(i)-a(i)-covy(i))*c(i)/b(i)+1.0)/c(i)/covt
            end if
          end if
          if (m<mi) then
            m = mi
          else if (m>ma) then
            m = ma
          end if

          err = sigma/coun(i)
          if (c(i)==0.0d0) then
            p = 1/b(i)
            q = (a(i)-meany(i))*p
            ww = p**2*err + p**2*omegas(1) + q**2*omegas(2)
          else
            mm = (meany(i)-a(i))*c(i)/b(i) + 1.0
            p = 1/(c(i)*(meany(i)-a(i))+b(i))
            q = (a(i)-b(i))/b(i)*p
            r = (meany(i)-a(i))/b(i)/c(i)/mm - log(mm)/c(i)**2
            ww = p**2*err + p**2*omegas(1) + q**2*omegas(2) + r**2*omegas(3)
          end if
          ww = 1/sqrt(ww)
          m = m - meanx(i)
          sum = sum + ww*m
          wsum = wsum + ww
        end do
        if (wsum>0.0d0) then
          offsett = sum/wsum
        else
          offsett = 0.0d0
        end if
        if(flag==0)then
          offt(id) = offsett
          realid(id) = datrec(1)
        else
          offt2(id) = offsett
        end if
        return
      end subroutine initialize

      subroutine eval(a, b, c, f, time, bm, id, covt, covy, flag)
        real (kind=dpsize), intent(in) :: a(*), b(*), c(*), covt, covy(*), time
        real (kind=dpsize), intent(out) :: f
        real (kind=dpsize) :: diseasetime
        integer (kind=isize), intent(in) :: bm, id, flag

        if(flag==0)then
          diseasetime = time + offt(id)
        else
          diseasetime = time + offt2(id)
        end if
!main biomarker
        if (bm==1) then
          if (c(bm)==0) then
            f = a(bm) + (covy(bm)+b(bm))*covt*diseasetime
          else
            f = a(bm) + (b(bm)/c(bm)-covy(bm))*(exp(c(bm)*covt*diseasetime)-1.0d0)
          end if
        else
!other biomarker
          if (c(bm)==0) then
            f = a(bm) + covy(bm) + b(bm)*covt*diseasetime
          else
            f = a(bm) + covy(bm) + b(bm)/c(bm)*(exp(c(bm)*covt*diseasetime)-1.0d0)
          end if
        end if
        return
      end subroutine eval

      subroutine makecsv
        integer (kind=isize) :: i

        open (13, file='s136_offsetT.csv')
        write (13, '(a2, a1, a7)') 'ID', ',', 'offsetT'
        do i = 1, numsj
          write (13, '(i5, a1, e10.4)') realid(i), ',', offt(i)
        end do
        close (13)
        return
      end subroutine makecsv

    end module sreft

    subroutine pred(icall, newind, theta, datrec, indxs, f, g, h)
      use :: sreft
      use :: prdims, only: gprd, hprd
      use :: nmprd_real, only: eta, eps
      use :: nmprd_int, only: newl2, iquit
      implicit none
      real (kind=dpsize), intent(in) :: theta(*), datrec(*)
      real (kind=dpsize), intent(out) :: f, g(gprd, 1), h(hprd, 1)
      real (kind=dpsize) :: f2, covt, step, time
      real (kind=dpsize), allocatable :: a(:), b(:), c(:), a2(:), b2(:), c2(:), meanx(:), meany(:), coun(:), covy(:)
      integer (kind=isize), intent(in) :: icall, newind, indxs(*)
      integer (kind=isize) :: id, bm, i, flag

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

      if (icall==0 .and. newind==0) then
        numbm = nepss_
        numsj = nindr
        if (.not. allocated(offt)) then
          allocate (offt(numsj), offt2(numsj), realid(numsj))
        end if
      end if

      allocate (a(numbm), b(numbm), c(numbm), a2(numbm), b2(numbm), c2(numbm), meanx(numbm), meany(numbm), coun(numbm), covy(numbm))

      id = datrec(1)
      time = datrec(2)
      bm = datrec(4)
      meanx = datrec(5:4+numbm)
      meany = datrec(5+numbm:4+2*numbm)
      coun = datrec(5+2*numbm:4+3*numbm)

      a = theta(1:numbm) + eta(1:numbm)
      b = theta(1+numbm:2*numbm) + eta(1+numbm:2*numbm)
      c = theta(1+2*numbm:3*numbm) + eta(1+2*numbm:3*numbm)
      covt = 1.0d0
      covy = 0.0d0
      step = 1.0d-6

      flag = 0
      if (newind<=1) call initialize(id, meanx, meany, coun, a, b, c, bm, datrec, covt, covy, flag)
      call eval(a, b, c, f, time, bm, id, covt, covy, flag)
      f = f + eps(bm)

      flag = 1
      do i = 1, numbm
        a2 = a
        b2 = b
        c2 = c

        a2(i) = a2(i) + step
        b2(i) = b2(i) + step
        c2(i) = c2(i) + step

        call initialize(id, meanx, meany, coun, a2, b, c, bm, datrec, covt, covy, flag)
        call eval(a2, b, c, f2, time, bm, id, covt, covy, flag)
        f2 = f2 + eps(bm)
        g(i, 1) = (f2 - f) / step

        call initialize(id, meanx, meany, coun, a, b2, c, bm, datrec, covt, covy, flag)
        call eval(a, b2, c, f2, time, bm, id, covt, covy, flag)
        f2 = f2 + eps(bm)
        g(i+numbm, 1) = (f2 - f) / step

        call initialize(id, meanx, meany, coun, a, b, c2, bm, datrec, covt, covy, flag)
        call eval(a, b, c2, f2, time, bm, id, covt, covy, flag)
        f2 = f2 + eps(bm)
        g(i+2*numbm, 1) = (f2 - f) / step
      end do
      h = 0.0d0
      h(bm, 1) = 1.0d0
      if (icall==3) call makecsv

      deallocate (a, b, c, a2, b2, c2, meanx, meany, coun, covy)

      return
    end subroutine pred
