	PROGRAM CCPHND

        PARAMETER (K=40,N=40,LRWK=1000,LIWK=1000,INCL=19)
	REAL Z(K,N), RWRK(LRWK), RLEVEL(INCL)
	INTEGER M, IWRK(LIWK)

	DATA RLEVEL /-20., -10., -5., -4., -3., -2., -1., 0., 
     +		1., 2., 3., 4., 5., 10., 20., 30., 40., 50., 100. /

	CALL GETDAT (Z, K, M, N) 

C Open GKS
	CALL OPNGKS

C Read in hand picked levels
	CALL CPSETI('CLS - CONTOUR LEVEL SELECTION FLAG',0)
	CALL CPSETI('NCL - NUMBER OF CONTOUR LEVELS',INCL)

	DO 10, I=1,INCL
	  CALL CPSETI('PAI - PARAMETER ARRAY INDEX',I)
	  CALL CPSETR('CLV - CONTOUR LEVEL VALUE',RLEVEL(I))
	  IF (AMOD(RLEVEL(I),5.).EQ.0)
     +		CALL CPSETI('CLU - CONTOUR LEVEL USE FLAG',3)
 10	CONTINUE

C Call conpack normally
	CALL CPRECT(Z,K,M,N,RWRK,LRWK,IWRK,LIWK)
	CALL CPBACK(Z, RWRK, IWRK)
	CALL CPCLDR(Z,RWRK,IWRK)

C Close frame and close GKS
	CALL FRAME
	CALL CLSGKS

	STOP
	END

	SUBROUTINE GETDAT (Z, K, M, N)

	REAL Z(K,N)
	INTEGER I,J,K,M,N

	M=K
	DO 10, I=1,M
	  DO 20, J=1,N
	    Z(I,J)= 10.E-5*(-16.*REAL(I**2*J) +
     +		    34*REAL(I*J**2) - REAL(6*I) + 93.)
  20	  CONTINUE
  10	CONTINUE

	RETURN
	END

