    module sreft
      use :: rocm_real, only: omega => varnf
      use :: nmprd_int, only: nthes_ => nwtht, netas_ => nweta, nepss_ => nweps
      use :: sizes, only: dpsize, isize
      use :: rocm_int, only: nindr => nindobs, indr1 => idxobsf, indr2 => idxobsl
      implicit none
      integer (kind=isize), save :: numbm, numcof, numsj
      integer (kind=isize), allocatable, save :: realid(:)
      real (kind=dpsize) :: plotmax = 40.0d0, plotmin = -15.0d0
      real (kind=dpsize) :: sum, wsum, ww, ma, mi, err, offsett, m, sigma, covt, covy
      real (kind=dpsize), allocatable, save :: offt(:)
    contains

      subroutine initialize(id, meanx, meany, coun, a, b, c, bm, datrec)
        real (kind=dpsize) :: meanx(*), meany(*), coun(*), a(*), b(*), c(*), datrec(*)
        integer (kind=isize) :: id, i, bm

        sum = 0.0d0
        wsum = 0.0d0
        do i = 1, numbm
          if (coun(i)==0) cycle

          sigma = omega(netas_+i, netas_+i)
          ma = plotmax - meanx(i)
          mi = plotmin + meanx(i)
          if (b(i)==0.0d0) cycle
          if (bm==1) then
!main biomarker
            if (c(i)==0.0d0) then
              m = (meany(i)-a(i))/(b(i)+covy)/covt
            else if ((meany(i)-a(i))/(b(i)/c(i)+covy)+1.0<0.0d0) then
              if (c(i)>0.0d0) then
                m = mi
              else
                m = ma
              end if
            else
              m = log((meany(i)-a(i))/(b(i)/c(i)+covy)+1.0)/c(i)/covt
            end if
          else
!other biomarker
            if (c(i)==0.0d0) then
              m = (meany(i)-a(i)-covy)/b(i)/covt
            else if ((meany(i)-a(i)-covy)*c(i)/b(i)+1.0<0.0d0) then
              if (c(i)>0.0d0) then
                m = mi
              else
                m = ma
              end if
            else
              m = log((meany(i)-a(i)-covy)*c(i)/b(i)+1.0)/c(i)/covt
            end if
          end if
          if (m<mi) then
            m = mi
          else if (m>ma) then
            m = ma
          else
            m = m
          end if

          err = sigma/coun(i)
          if (c(i)/=0.0) then
            ww = err/(b(i)*covt*exp(c(i)*covt*m))**2
          else
            ww = err/(b(i)*covt)**2
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
        offt(id) = offsett
        realid(id) = datrec(1)
        return
      end subroutine

      subroutine eval(a, b, c, f, g, datrec, bm, id, meany)
        use :: prdims, only: gprd
        real (kind=dpsize) :: datrec(*), f, g(gprd, 1), a(*), b(*), c(*), time, meany(*)
        integer (kind=isize) :: bm, id, icall

        time = datrec(2) + offt(id)
!main biomarker
        if (bm==1) then
          if (c(bm)==0) then
            f = a(bm) + (covy+b(bm))*covt*time
            g = 0.0d0
            g(bm, 1) = 1.0d0
            g(bm+numbm, 1) = covt*time
          else
            f = a(bm) + (b(bm)/c(bm)-covy)*(exp(c(bm)*covt*time)-1.0d0)
            g = 0.0d0
            g(bm, 1) = 1.0d0
            g(bm+numbm, 1) = (exp(c(bm)*covt*time)-1.0d0)/c(bm)
            g(bm+2*numbm, 1) = (b(bm)/c(bm)-covy)*covt*time*exp(c(bm)*covt*time) - &
              b(bm)/c(bm)**2*(exp(c(bm)*covt*time)-1.0d0)
          end if
        else
!other biomarker
          if (c(bm)==0) then
            f = a(bm) + covy + b(bm)*covt*time
            g = 0.0d0
            g(bm, 1) = 1.0d0
            g(bm+numbm, 1) = covt*time
          else
            f = a(bm) + covy + b(bm)/c(bm)*(exp(c(bm)*covt*time)-1.0d0)
            g = 0.0d0
            g(bm, 1) = 1.0d0
            g(bm+numbm, 1) = (exp(c(bm)*covt*time)-1.0d0)/c(bm)
            g(bm+2*numbm, 1) = b(bm)/c(bm)*covt*time*exp(c(bm)*covt*time) - &
              b(bm)/c(bm)**2*(exp(c(bm)*covt*time)-1.0d0)
          end if
        end if
        return
      end subroutine

      subroutine makecsv
        implicit none
        integer (kind=isize) :: i

        open (13, file='offsetT.csv')
        write (13, '(a2, a2, a6, a2, a7)') 'id', ', ', 'realid', ',', 'offsetT'
        do i = 1, numsj
          write (13, '(i3, a2, i4, a2, e12.4)') i, ', ', realid(i), ',', offt(i)
        end do
        close (13)
        return
      end subroutine

    end module

    subroutine pred(icall, newind, theta, datrec, indxs, f, g, h)
      use :: sreft
      use :: prdims, only: gprd, hprd
      use :: nmprd_real, only: eta, eps
      use :: nmprd_int, only: newl2, iquit
      implicit none
      real (kind=dpsize) :: theta(*), datrec(*), f, g(gprd, 1), h(hprd, 1)
      real (kind=dpsize), allocatable :: a(:), b(:), c(:), meanx(:), meany(:), coun(:), cof(:)
      integer (kind=isize) :: icall, newind, indxs(*), id, bm, gend

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

      if (newind==0) then
        numbm = nepss_
        numcof = (nthes_-nepss_*3)/numbm
        numsj = nindr
        if (icall==0) then
          allocate (offt(numsj), realid(numsj))
        end if
      end if
      allocate (a(numbm), b(numbm), c(numbm), meanx(numbm), meany(numbm), coun(numbm), &
        cof(numcof*numbm))
      id = datrec(5)
      bm = datrec(4)
      meanx = datrec(6:5+numbm)
      meany = datrec(6+numbm:5+2*numbm)
      coun = datrec(6+2*numbm:5+3*numbm)
      a = theta(1:numbm) + eta(1:numbm)
      b = theta(1+numbm:2*numbm) + eta(1+numbm:2*numbm)
      c = theta(1+2*numbm:3*numbm) + eta(1+2*numbm:3*numbm)
!define covariate
      if (numcof==0) then
        covt = 1.0d0
        covy = 0.0d0
      else
        gend = datrec(6+3*numbm)
        cof = theta(3*numbm+1:3*(numbm+numcof))
        if (gend==1) then
          covt = (1.0d0+cof(bm))
          covy = cof(numbm+bm)
        else
          covt = 1.0d0
          covy = 0.0d0
        end if
      end if
      if (newind<=1) call initialize(id, meanx, meany, coun, a, b, c, bm, datrec)
      call eval(a, b, c, f, g, datrec, bm, id, meany)
      f = f + eps(bm)
      h = 0.0d0
      h(bm, 1) = 1.0d0
      if (icall==3) call makecsv
      deallocate (a, b, c, meanx, meany, coun, cof)
      return
    end subroutine
