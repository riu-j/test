    Module variables
      Implicit None
      Real (Kind=selected_real_kind(8)), Save :: tc(3), offt(400), generalmeanage = 70.0D0, generaltriggerage = 60.0D0, &
        plotmax = 25.0D0, plotmin = 5.0D0, ww, ma, mi, err, offsett, m, sum = 0.0D0, wsum = 0.0D0
      Real (Kind=selected_real_kind(8)), Save :: epss(1:3) = (/ 0.04D0, 0.0025D0, 0.01D0 /)
      Integer (Kind=selected_int_kind(8)), Save :: numbm = 3, numsj = 400
    Contains

      Subroutine initialize(id, bm, meanx, meany, coun, a, b, c)
        Real (Kind=selected_real_kind(8)) :: meanx(3), meany(3), coun(3), a(3), b(3), c(3)
        Integer (Kind=selected_int_kind(8)) :: id, bm, i

        Do i = 1, numbm
          If (b(i)==0.0) Cycle
          Select Case (i)
          Case (3)
            Select Case (c(i))
            Case (0)
              m = (meany(i)-a(i))/b(i)
            Case Default
              If ((meany(i)-a(i))*c(i)/b(i)+1.0<0) Then
                If (c(i)>0) Then
                  m = mi
                Else
                  m = ma
                End If
              Else
                m = log((meany(i)-a(i))*c(i)/b(i)+1.0)/c(i)
              End If
            End Select
          Case (:2)
            Select Case (c(i))
            Case (0)
              m = (meany(i)-a(i))
            Case Default
              if((meany(i)-a(i))*c(i)/b(i)+1.0<0) then
              If(c(i)>0) Then
                m=mi
              Else
                m=ma
              End If
            Else
              m=log((meany(i)-a(i))*c(i)/b(i)+1.0)/c(i)
            End If
          End Select
        End Select
        ma=plotmax-meanx(i)
        mi=plotmin+meanx(i)
        If(m>0.25) Then
          If(c(i)>0) Then
            tc(i)=mi
          Else
            tc(i)=ma
          End If
        Else
          tc(i)=(1-sqrt(1-4.0D0*m))/2.0D0
        End If
        If(tc(i)<mi) Then
          tc(i)=mi
        Else If(tc(i)>ma) Then
          tc(i)=ma
        Else
          tc(i)=tc(i)
        End If
        m=1-(tc(i)-meanx(i))
        b(i)=b(i)*m
        c(i)=c(i)*m
        err=epss(i)/coun(i)
        If(c(i)/=0.0) Then
          If(i==3) Then
            b(i)=b(i)-c(i)
          End If
          ww=err/((b(i)*exp(c(i)*tc(i)))*(b(i)*exp(c(i)*tc(i))))
        Else
          If(i==3) Then
            b(i)=b(i)+m
          End If
          ww=err/(b(i)*b(i))
        End If
        ww=1/ww
        m=tc(i)-meanx(i)
        sum=sum+ww*m
        wsum=wsum+ww
      End Do
      If(wsum>0.0) Then
        offsett=sum/wsum
      Else
        offsett=0.0D0
      End If
      offt(id)=offsett
      Return
    End Subroutine

    Subroutine eval(a,b,c,f,g,datrec,bm)
      Use prdims, Only:gprd
      Real(Kind=selected_real_kind(8)) :: datrec(*), f, g(gprd,1), a(3), b(3), c(3), time
      Integer(Kind=selected_int_kind(8)) :: bm, id

      time=datrec(3)+offt(id)

      If(c(bm)==0) Then
        f=a(bm)+b(bm)*time
        g(bm,1)=1.0D0
        g(bm+3,1)=time
        g(bm+6,1)=0.0D0
      Else
        f=a(bm)+b(bm)/c(bm)*(exp(c(bm)*time)-1.0D0)
        g(bm,1)=1.0D0
        g(bm+3,1)=(exp(c(bm)*time)-1.0D0)/c(bm)
        g(bm+6,1)=b(bm)*(exp(c(bm)*time)-1.0D0)*(time-1.0D0/c(bm)/c(bm))
      End If
      Return
    End Subroutine

    Subroutine makecsv
      Implicit None
      Integer(Kind=selected_int_kind(8)) :: i

      Open(13,File='offsetT.csv')
      Do i=1, numsj
        Write(13,*) i, ', ', offt(i)
      End Do
      Close(13)
      Return
    End Subroutine

  End Module

  Subroutine pred(icall,newind,theta,datrec,indxs,f,g,h)
    Use variables
    Use prdims, Only:gprd, hprd
    Use nmprd_real, Only:eta, eps
    Implicit None
    Real(Kind=selected_real_kind(8)) :: theta(*), datrec(*), f, g(gprd,1), h(hprd,1), a(3), b(3), c(3), meanx(3), &
      meany(3), coun(3)
    Integer(Kind=selected_int_kind(8)) :: icall, newind, indxs(*), id, bm

    id=datrec(1)
    bm=datrec(2)
    meanx(1)=datrec(5)
    meanx(2)=datrec(6)
    meanx(3)=datrec(7)
    meany(1)=datrec(8)
    meany(2)=datrec(9)
    meany(3)=datrec(10)
    coun(1)=datrec(11)
    coun(2)=datrec(12)
    coun(3)=datrec(13)
    a(1)=theta(1)+eta(1)
    a(2)=theta(2)+eta(2)
    a(3)=theta(3)+eta(3)
    b(1)=theta(4)+eta(4)
    b(2)=theta(5)+eta(5)
    b(3)=theta(6)+eta(6)
    c(1)=theta(7)+eta(7)
    c(2)=theta(8)+eta(8)
    c(3)=theta(9)+eta(9)

    If(newind<=1) Call initialize(id,bm,meanx,meany,coun,a,b,c)

    Call eval(a,b,c,f,g,datrec,bm)

    f=f*exp(eps(bm))
    h(bm,1)=f

    If(icall==3) Call makecsv

    Return
  End Subroutine
