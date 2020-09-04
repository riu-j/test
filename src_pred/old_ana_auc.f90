    module sreft
      use :: rocm_real, only: omega => varnf
      use :: nmprd_int, only: nthes_ => nwtht, netas_ => nweta, nepss_ => nweps
      use :: sizes, only: dpsize, isize
      use :: rocm_int, only: nindr => nindobs, indr1 => idxobsf, indr2 => idxobsl
      implicit none
      integer (kind=isize), save :: numbm, numsj
      integer (kind=isize), allocatable, save :: realid(:)
      real (kind=dpsize) :: plotmax = 40.0d0, plotmin = -15.0d0
      real (kind=dpsize), allocatable, save :: offt(:)
    contains

      subroutine initialize(id, med, aucy, inter, coun, a, b, c, bm, datrec, covt, covy)
        real (kind=dpsize), intent(in) :: med(*), aucy(*), inter(*), coun(*), a(*), b(*), c(*), datrec(*), covt, covy(*)
        real (kind=dpsize) :: sum, wsum, ww, ma, mi, err, offsett, m, sigma, check(3)
        integer (kind=isize), intent(in) :: id, bm
        integer (kind=isize) :: i

        sum = 0.0d0
        wsum = 0.0d0
        do i = 1, numbm
          if (coun(i)==0 .or. b(i)==0.0d0) cycle
          sigma = omega(netas_+i, netas_+i)
          ma = plotmax - med(i)
          mi = plotmin + med(i)
          if (coun(i)==1) then
            if (bm==1) then
!main biomarker
              if (c(i)==0.0d0) then
                m = (aucy(i)-a(i))/(b(i)+covy(i))/covt
              else if ((aucy(i)-a(i))/(b(i)/c(i)+covy(i))+1.0<0.0d0) then
                if (c(i)>0.0d0) then
                  m = mi
                else
                  m = ma
                end if
              else
                m = log((aucy(i)-a(i))/(b(i)/c(i)+covy(i))+1.0)/c(i)/covt
              end if
            else
!other biomarker
              if (c(i)==0.0d0) then
                m = (aucy(i)-a(i)-covy(i))/b(i)/covt
              else if ((aucy(i)-a(i)-covy(i))*c(i)/b(i)+1.0<0.0d0) then
                if (c(i)>0.0d0) then
                  m = mi
                else
                  m = ma
                end if
              else
                m = log((aucy(i)-a(i)-covy(i))*c(i)/b(i)+1.0)/c(i)/covt
              end if
            end if
          else
            check(1) = aucy(i)*c(i)**2 - 2*inter(i)*c(i)*(a(i)*c(i)-b(i))
            check(2) = 2*b(i)*sinh(inter(i)*c(i))
            check(3) = check(1)/check(2)
            if (c(i)==0.0d0) then
              m = (aucy(i)/2*inter(i)-a(i))/b(i)
            else if (check(3)<=0.0d0) then
              if (c(i)>0.0d0) then
                m = mi
              else
                m = ma
              end if
            else
              m = log(check(3))/c(i)
            end if
          end if
          if (m<mi) then
            m = mi
          else if (m>ma) then
            m = ma
          end if
          err = sigma/coun(i)
          if (c(i)/=0.0) then
            ww = err/(b(i)*covt*exp(c(i)*covt*m))**2
          else
            ww = err/(b(i)*covt)**2
          end if
          ww = 1/sqrt(ww)
          m = m - med(i)
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
      end subroutine initialize

      subroutine eval(a, b, c, f, g, time, bm, id, covt, covy)
        use :: prdims, only: gprd
        real (kind=dpsize), intent(in) :: a(*), b(*), c(*), covt, covy(*)
        real (kind=dpsize), intent(out) :: f, g(gprd, 1)
        real (kind=dpsize), intent(inout) :: time
        integer (kind=isize) :: bm, id

        time = time + offt(id)
!main biomarker
        if (bm==1) then
          if (c(bm)==0) then
            f = a(bm) + (covy(bm)+b(bm))*covt*time
            g = 0.0d0
            g(bm, 1) = 1.0d0
            g(bm+numbm, 1) = covt*time
          else
            f = a(bm) + (b(bm)/c(bm)-covy(bm))*(exp(c(bm)*covt*time)-1.0d0)
            g = 0.0d0
            g(bm, 1) = 1.0d0
            g(bm+numbm, 1) = (exp(c(bm)*covt*time)-1.0d0)/c(bm)
            g(bm+2*numbm, 1) = (b(bm)/c(bm)-covy(bm))*covt*time*exp(c(bm)*covt*time) - &
              b(bm)/c(bm)**2*(exp(c(bm)*covt*time)-1.0d0)
          end if
        else
!other biomarker
          if (c(bm)==0) then
            f = a(bm) + covy(bm) + b(bm)*covt*time
            g = 0.0d0
            g(bm, 1) = 1.0d0
            g(bm+numbm, 1) = covt*time
          else
            f = a(bm) + covy(bm) + b(bm)/c(bm)*(exp(c(bm)*covt*time)-1.0d0)
            g = 0.0d0
            g(bm, 1) = 1.0d0
            g(bm+numbm, 1) = (exp(c(bm)*covt*time)-1.0d0)/c(bm)
            g(bm+2*numbm, 1) = b(bm)/c(bm)*covt*time*exp(c(bm)*covt*time) - b(bm)/c(bm)**2*(exp(c(bm)*covt*time)-1.0d0)
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
      real (kind=dpsize) :: covt, time
      real (kind=dpsize), allocatable :: a(:), b(:), c(:), med(:), aucy(:), inter(:), coun(:), covy(:)
      integer (kind=isize), intent(in) :: icall, newind, indxs(*)
      integer (kind=isize) :: id, bm

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

      allocate (a(numbm), b(numbm), c(numbm), med(numbm), aucy(numbm), inter(numbm), coun(numbm), covy(numbm))

      id = datrec(1)
      time = datrec(2)
      bm = datrec(4)
      med = datrec(5:4+numbm)
      aucy = datrec(5+numbm:4+2*numbm)
      inter = datrec(5+2*numbm:4+3*numbm)
      coun = datrec(5+3*numbm:4+4*numbm)

      a = theta(1:numbm) + eta(1:numbm)
      b = theta(1+numbm:2*numbm) + eta(1+numbm:2*numbm)
      c = theta(1+2*numbm:3*numbm) + eta(1+2*numbm:3*numbm)
      covt = 1.0d0
      covy = 0.0d0

      if (newind<=1) call initialize(id, med, aucy, inter, coun, a, b, c, bm, datrec, covt, covy)

      call eval(a, b, c, f, g, time, bm, id, covt, covy)

      f = f + eps(bm)
      h = 0.0d0
      h(bm, 1) = 1.0d0

      if (icall==3) call makecsv

      deallocate (a, b, c, med, aucy, inter, coun, covy)

      return
    end subroutine pred
