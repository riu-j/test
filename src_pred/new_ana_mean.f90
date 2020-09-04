    module sreft
      use :: rocm_real, only: omega => varnf
      use :: nmprd_int, only: nthes_ => nwtht, netas_ => nweta, nepss_ => nweps
      use :: sizes, only: dpsize, isize
      use :: rocm_int, only: nindr => nindobs, indr1 => idxobsf, indr2 => idxobsl
      implicit none
      integer (kind=isize), save :: numbm, numsj
      integer (kind=isize), allocatable, save :: realid(:)
      real (kind=dpsize) :: plotmax = 40.0d0, plotmin = -15.0d0, undef = -9999.0d0
      real (kind=dpsize), allocatable, save :: offt(:), df_dv(:, :, :), df_time(:, :, :), meanxs(:, :, :), counts(:, :, :)
    contains

      function mean(data) result(out)
        real(kind=dpsize) :: data(:)
        real(kind=dpsize) :: out
        integer(kind=isize) :: i, N

        out = 0
        N = 0
        do i = 1, size(data)
            if(data(i) /= undef) then
                out = out + data(i)
                N = N + 1
            end if
        end do
        out = out / N
      end function mean

      function length(data) result(out)
      real(kind=dpsize) :: data(:)
      real(kind=dpsize) :: out

      out = count(mask = data /= undef)
      end function length

      subroutine initialize(id, meanx, yval, coun, a, b, c, bm, datrec, covt, covy)
        real (kind=dpsize), intent(in) :: meanx(*), coun(*), a(*), b(*), c(*), datrec(*), covt, covy(*), yval(:, :)
        real (kind=dpsize), allocatable :: obs(:), hoge(:), m(:), now_yval(:)
        real (kind=dpsize) :: sum, wsum, ww, finalm, ma, mi, err, offsett, sigma
        integer (kind=isize), intent(in) :: id, bm
        integer (kind=isize) :: i, now_size

        sum = 0.0d0
        wsum = 0.0d0
        do i = 1, numbm
          if (coun(i)==0) cycle
          now_size = count(mask = yval(:, i)/=undef)
          allocate(obs(now_size), hoge(now_size), m(now_size))
          where (yval(:, i)/=undef)
            obs = yval(:, i)
          end where
          sigma = omega(netas_+i, netas_+i)
          ma = plotmax - meanx(i)
          mi = plotmin + meanx(i)
          if (b(i)==0.0d0) cycle
          if (bm==1) then
!main biomarker
            if (c(i)==0.0d0) then
              m = (obs-a(i))/(b(i)+covy(i))/covt
            else
              hoge = (obs-a(i))/(b(i)/c(i)+covy(i))+1.0
              if (c(i)>0.0d0)then
                where (hoge>0)
                  m = log(hoge)/c(i)/covt
                else where
                  m = mi
                end where
              else
                where (hoge>0)
                  m = log(hoge)/c(i)/covt
                else where
                  m = ma
                end where
              end if
            end if
          else
!other biomarker
            if (c(i)==0.0d0) then
              m = (obs-a(i))/(b(i)+covy(i))/covt
            else
              hoge = (obs-a(i))/(b(i)/c(i)+covy(i))+1.0
              if (c(i)>0.0d0)then
                where (hoge>0)
                  m = log(hoge)/c(i)/covt
                else where
                  m = mi
                end where
              else
                where (hoge>0)
                  m = log(hoge)/c(i)/covt
                else where
                  m = ma
                end where
              end if
            end if
          end if
          where (m<mi)
            m = mi
          else where (m>ma)
            m = ma
          end where
          err = sigma/coun(i)
          if (c(i)/=0.0) then
            ww = err/(b(i)*covt*exp(c(i)*covt*mean(m)))**2
          else
            ww = err/(b(i)*covt)**2
          end if
          ww = 1/sqrt(ww)
          finalm = mean(m) - meanx(i)
          sum = sum + ww*finalm
          wsum = wsum + ww
          deallocate(obs, hoge, m)
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

      subroutine readcsv
        integer (kind=isize) :: i, n, k, read_id, read_bm, read_id2, old_id, old_bm, lengthofdf=15
        real (kind=dpsize) :: read_time, read_dv

        print *, "start read csv file"
        allocate(df_dv(lengthofdf, numbm, numsj), df_time(lengthofdf, numbm, numsj))
        df_dv = undef
        df_time = undef
        open (14, file='s136_data.csv', status='old')
        n = 0
        read(14, '()')
        do
          read (14, *, end=100) read_id, read_time, read_dv, read_bm, read_id2
          n = n + 1
        end do
        100 continue
        rewind (14)
        read(14, '()')
        old_id = undef
        old_bm = undef
        k = 1
        do i = 1, n
          read(14, *) read_id, read_time, read_dv, read_bm, read_id2
          ! print *, read_id, read_time, read_dv, read_bm, read_id2
          if (old_id==read_id .and. old_bm==read_bm)then
            k = k + 1
            if (k>lengthofdf)then
              print *, "PRED_SReFT detect problem!"
              print *, "data length is over limitaion. Please change length of df."
              stop
            end if
          else
            k = 1
          end if
          df_dv(k, read_bm, read_id) = read_dv
          df_time(k, read_bm, read_id) = read_time
          old_id = read_id
          old_bm = read_bm
        end do
        close (14)
        print *, "finish read csv file"
      end subroutine readcsv

      subroutine calcmeancount
        integer (kind=isize) :: i, j

        print *, "start calculate meanx and count"
        allocate(meanxs(1, numbm, numsj), counts(1, numbm, numsj))
        do i = 1, numsj
          do j = 1, numbm
            meanxs(1, j, i) = mean(df_time(:, j, i))
            counts(1, j, i) = length(df_dv(:, j, i))
          end do
        end do
        print *, "finish calculate meanx and count"
      end subroutine calcmeancount

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
      real (kind=dpsize), allocatable :: a(:), b(:), c(:), meanx(:), coun(:), covy(:), yval(:, :)
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
          allocate (offt(numsj), realid(numsj))
        end if
        call readcsv
        call calcmeancount
      end if
      allocate (a(numbm), b(numbm), c(numbm), meanx(numbm), coun(numbm), covy(numbm), yval(15, numbm))
      id = datrec(1)
      time = datrec(2)
      bm = datrec(4)
      meanx = meanxs(1, :, id)
      coun = counts(1, :, id)
      yval = df_dv(:, :, id)
      a = theta(1:numbm) + eta(1:numbm)
      b = theta(1+numbm:2*numbm) + eta(1+numbm:2*numbm)
      c = theta(1+2*numbm:3*numbm) + eta(1+2*numbm:3*numbm)
      covt = 1.0d0
      covy = 0.0d0
      if (newind<=1) call initialize(id, meanx, yval, coun, a, b, c, bm, datrec, covt, covy)
      call eval(a, b, c, f, g, time, bm, id, covt, covy)
      f = f + eps(bm)
      h = 0.0d0
      h(bm, 1) = 1.0d0
      if (icall==3) call makecsv
      deallocate (a, b, c, meanx, coun, covy, yval)
      return
    end subroutine
