MODULE variables
implicit none
REAL(kind = selected_real_kind(8)), SAVE ::  offT(400)
REAL(kind = selected_real_kind(8)) :: generalmeanage = 70.0d0, &
&generaltriggerage = 60.0d0, plotmax = 25.0d0, plotmin = 5.0d0, sum, wsum, ssum, ww, ma, mi, err&
&, tc(3), epss(3), offsetT, m, bb, cc
INTEGER(kind=selected_int_kind(8)):: NumBM = 3
contains

subroutine insertoffT(ID)
implicit none
!REAL(KIND=selected_real_kind(8)) ::
INTEGER(kind = selected_int_kind(8)):: ID
offT(ID) = offsetT
return
end

subroutine makecsv
implicit none
!REAL(KIND=selected_real_kind(8)) ::offtT
INTEGER(kind = selected_int_kind(8)):: i
OPEN(13, file = 'offsetT.csv')
	DO i = 1, 400 !Number of subjects
		write(13, *) i, offT(i)
	END DO
CLOSE (13)
return
end

END MODULE

!program test
!USE variables
!DIMENSION THETA(9), DATREC(13), G(3, 1), H(3, 1)
!REAL(kind= selected_real_kind(8)) THETA, DATREC, F, G, H
!INTEGER(kind= selected_int_kind(8))  ICALL,NEWIND,INDXS(5)
!theta(1) =-1.7516 !a1
!theta(2) =0.0587 !a2
!theta(3) =0.1160 !a3
!theta(4) =1.7498 !b1
!theta(5) =-0.0655 !b2
!theta(6) =0.1009 !b3
!theta(7) =-0.2291 !c1
!theta(8) =-0.0830 !c2
!theta(9) =-0.072 !c3
!DATREC(1) = 1 !ID
!DATREC(2) = 2 !BM
!DATREC(3) = 0 !TIME
!DATREC(5) =2 !MX1
!DATREC(6) =2 !MX2
!DATREC(7) =2 !MX3
!DATREC(8) =-2.038317267 !MY1
!DATREC(9) =0.836190448 !MY2
!DATREC(10) =0.337930273 !MY3
!DATREC(11) =5 !Count1
!DATREC(12) =5 !Count2
!DATREC(13) =5 !Count3
!NEWIND = 0
!CALL PRED(ICALL, NEWIND, THETA, DATREC, INDXS, F, G, H)
!end

SUBROUTINE PRED(ICALL, NEWIND, THETA, DATREC, INDXS, F, G, H)
use variables
use PRDIMS, ONLY: GPRD, HPRD
use NMPRD_REAL, ONLY: ETA, EPS
implicit none
REAL(kind= selected_real_kind(8)) :: THETA(*), DATREC(*), F, G(GPRD, *), H(HPRD, *), a(3) ,b(3), c(3), &
&MeanX(3), MeanY(3), Coun(3), TIME
INTEGER(kind = selected_int_kind(8)) :: ICALL, NEWIND, INDXS(*), ID, BM
!write(6, 80) eta
ID = DATREC(1)
BM = DATREC(2)
MeanX(1) = DATREC(5)
MeanX(2) = DATREC(6)
MeanX(3) = DATREC(7)
MeanY(1) = DATREC(8)
MeanY(2) = DATREC(9)
MeanY(3) = DATREC(10)
Coun(1) = DATREC(11)
Coun(2) = DATREC(12)
Coun(3) = DATREC(13)
a(1) = THETA(1) + ETA(1)
a(2) = THETA(2) + ETA(2)
a(3) = THETA(3) + ETA(3)
b(1) = THETA(4) + ETA(4)
b(2) = THETA(5) + ETA(5)
b(3) = THETA(6) + ETA(6)
c(1) = THETA(7) + ETA(7)
c(2) = THETA(8) + ETA(8)
c(3) = THETA(9) + ETA(9)

IF(NEWIND <= 1)CALL INITIALIZE(ID, BM, MeanX, MeanY, Coun, a, b, c)
TIME = DATREC(3) + offT(ID)

!G = 0.0d0

CALL EVALUATE(TIME, BM, a, b, c, F, G)
!write(6,81) eps
F = F * EXP(EPS(BM))

!H = 0.0d0
H(BM, 1) = F
!write(6,80) H
!print *, F

!IF(ICALL == 3)call makecsv

80	FORMAT ('eta', 10E12.4)
81	FORMAT ('eps', 10E12.4)
RETURN
END

SUBROUTINE EVALUATE(TIME, BM, a, b, c, F, G)!evaluate for compartment
USE variables
use PRDIMS, ONLY: GPRD, HPRD
implicit none
REAL(kind = selected_real_kind(8)) :: F, G(GPRD, *),a(*) ,b(*), c(*), TIME
INTEGER(kind = selected_int_kind(8)) BM
IF(c(BM) == 0)THEN
	F = a(BM) + b(BM) * TIME
	G(BM,1) = 1.0d0
	G(BM+3,1) = TIME
	G(BM+6,1) = 0.0d0
ELSE
	F = a(BM) + b(BM) / c(BM) * (EXP(c(BM) * TIME) - 1.0d0)
	G(BM,1) = 1.0d0
	G(BM+3,1) = (EXP(c(BM) * TIME) - 1.0d0) / c(BM)
	G(BM+6,1) = b(BM) * (EXP(c(BM) * TIME) - 1.0d0) * (TIME - 1.0d0 / c(BM) / c(BM))
END IF
RETURN
END


SUBROUTINE INITIALIZE(ID, BM, MeanX, MeanY, Coun, a, b, c)!preparativeCalculation
USE variables
implicit none
REAL(kind = selected_real_kind(8)) :: MeanX(*), MeanY(*), Coun(*), a(*) ,b(*), c(*)
INTEGER(kind = selected_int_kind(8)) :: ID, BM, i
!define eps
epss(1) = 0.04d0
epss(2) = 0.0025d0
epss(3) = 0.01d0
!end define eps
!write(6, 81) a
!write(6, 82) b
!write(6, 83) c

DO i=1, NumBM!Number of biomarkers
	IF(b(i) == 0.0)CYCLE
	IF(i == 3)THEN
		IF(c(i) == 0)THEN
			m = (MeanY(i) - a(i)) / b(i)
		ELSE
			m = log((MeanY(i) - a(i)) / (b(i) / c(i)) + 1.0) / c(i)
		END IF
	ELSE
		IF(c(i) == 0)THEN
			m = (MeanY(i) - a(i))
		ELSE
			m = log((MeanY(i) - a(i)) * c(i) / b(i) + 1.0) / c(i)
		END IF
	END IF
	ma = plotmax - MeanX(i)
	mi = plotmin + MeanX(i)
	IF(m > 0.25)THEN
		IF(c(i) > 0)THEN
			tc(i) = mi
		ELSE
			tc(i) = ma
		END IF
	ELSE
		tc(i) = (1 - sqrt(1 - 4.0d0 * m)) /  2.0d0
	ENDIF
	IF(tc(i) < mi)THEN
		tc(i) = mi
	ELSEIF(tc(i) > ma)THEN
			tc(i) = ma
		ELSE
			tc(i) = tc(i)
	END IF
	m = 1 - (tc(i) - MeanX(i))
	bb = b(i) * m
	cc = c(i) * m
!	call sim(eps)
	err = epss(i) / Coun(i)
	IF(c(i) /= 0.0)THEN
		IF(i == 3)THEN
			bb = bb - cc
		END IF
		ww = err / ((bb * exp(cc * tc(i))) * (bb * exp(cc * tc(i))))
	ELSE
		IF(i == 3)THEN
			bb = bb + m
		END IF
		ww = err / (bb * bb)
	END IF
	ww = 1 / ww
	m = tc(i) - MeanX(i)
	sum = sum + ww * m
!	ssum = ssum + ww * m * m
	wsum = wsum + ww
END DO
IF(wsum > 0.0)THEN
	offT(ID) = sum / wsum
!	offsetD = sqrt(ssum / wsum - offsetT * offsetT)
ELSE
	offT(ID) = 0.0d0
!	offsetD = 0.0d0
END IF
!call insertoffT(ID)
!write(6, 80) offT
!triggerA = Age - offsetT
80	FORMAT (' offsetT', 10E12.4)
81	FORMAT (' a', 10E12.4)
82	FORMAT (' b', 10E12.4)
83	FORMAT (' c', 10E12.4)
RETURN
END
