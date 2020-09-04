MODULE variables
implicit none
REAL(kind = selected_real_kind(8)), SAVE ::  tc(3), offT(400), generalmeanage = 70.0d0, generaltriggerage = 60.0d0&
&, plotmax = 25.0d0, plotmin = 5.0d0, ww, ma, mi, err, offsetT, m, sum = 0.0d0, wsum = 0.0d0
REAL(kind = selected_real_kind(8)), save :: epss(1:3) = (/ 0.04d0, 0.0025d0, 0.01d0 /)
INTEGER(kind = selected_int_kind(8)), save :: NumBM = 3, NumSJ = 400
contains

SUBROUTINE INITIALIZE(ID, BM, MeanX, MeanY, Coun, a, b, c)
REAL(kind = selected_real_kind(8)) :: MeanX(3), MeanY(3), Coun(3), a(3) ,b(3), c(3)
INTEGER(kind = selected_int_kind(8)) :: ID, BM, i
DO i=1, NumBM!Number of biomarkers
	IF(b(i) == 0.0)CYCLE
	IF(i == 3)THEN
		IF(c(i) == 0)THEN
			m = (MeanY(i) - a(i)) / b(i)
		ELSEif((MeanY(i) - a(i)) / (b(i) / c(i)) + 1.0 < 0)then
				if(c(i) > 0)then
					m = mi
				else
					m = ma
				endif
			else
				m = log((MeanY(i) - a(i)) / (b(i) / c(i)) + 1.0) / c(i)
			endif
		END IF
	endif
	if(i /= 3)then
		IF(c(i) == 0)THEN
			m = (MeanY(i) - a(i))
		ELSEif(((MeanY(i) - a(i)) * c(i) / b(i) + 1.0 < 0)then
				if(c(i) > 0)then
					m = mi
				else
					m = ma
				endif
			else
				m = log((MeanY(i) - a(i)) * c(i) / b(i) + 1.0) / c(i)
			endif
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
	b(i) = b(i) * m
	c(i) = c(i) * m
	err = epss(i) / Coun(i)
	IF(c(i) /= 0.0)THEN
		IF(i == 3)THEN
			b(i) = b(i) - c(i)
		END IF
		ww = err / ((b(i) * exp(c(i) * tc(i))) * (b(i) * exp(c(i) * tc(i))))
	ELSE
		IF(i == 3)THEN
			b(i) = b(i) + m
		END IF
		ww = err / (b(i) * b(i))
	END IF
	ww = 1 / ww
	m = tc(i) - MeanX(i)
	sum = sum + ww * m
	wsum = wsum + ww
END DO
IF(wsum > 0.0)THEN
	offsetT = sum / wsum
ELSE
	offsetT = 0.0d0
END IF
offT(ID) = offsetT
RETURN
END subroutine

subroutine eval(a, b, c, F, G, DATREC, BM)
use PRDIMS, ONLY: GPRD
REAL(kind= selected_real_kind(8)) :: DATREC(*), F, G(GPRD, 1), a(3) ,b(3), c(3), TIME
INTEGER(kind = selected_int_kind(8)) :: BM, ID

TIME = DATREC(3) + offT(ID)

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
return
end subroutine

subroutine makecsv
implicit none
INTEGER(kind = selected_int_kind(8)):: i
OPEN(13, file = 'offsetT.csv')
	DO i = 1, NumSJ
		write(13, *) i, ", ", offT(i)
	END DO
CLOSE (13)
return
end subroutine

END MODULE

SUBROUTINE PRED(ICALL, NEWIND, THETA, DATREC, INDXS, F, G, H)
use variables
use PRDIMS, ONLY: GPRD, HPRD
use NMPRD_REAL, ONLY: ETA, EPS
implicit none
REAL(kind= selected_real_kind(8)) :: THETA(*), DATREC(*), F, G(GPRD, 1), H(HPRD, 1), a(3) ,b(3), c(3), MeanX(3), MeanY(3), Coun(3)
INTEGER(kind = selected_int_kind(8)) :: ICALL, NEWIND, INDXS(*), ID, BM
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

call eval(a, b, c, F, G, DATREC, BM)

F = F * EXP(EPS(BM))
H(BM, 1) = F

IF(ICALL == 3)call makecsv

RETURN
END subroutine
