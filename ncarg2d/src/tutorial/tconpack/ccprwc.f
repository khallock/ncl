	PROGRAM CCPRWC

        PARAMETER (K=400,N=400,LRWK=2000,LIWK=2000)
	REAL Z(K,N), RWRK(LRWK)
	INTEGER M, IWRK(LIWK)

	CALL GETDAT (Z, K, M, N) 

C Open GKS
	CALL OPNGKS

C Initialize Conpack
	CALL CPRECT(Z,K,M,N,RWRK,LRWK,IWRK,LIWK)
C Draw perimeter
	CALL CPBACK(Z, RWRK, IWRK)

C Turn on line labels for every line
	CALL CPPKCL(Z, RWRK, IWRK)
	CALL CPGETI('NCL - NUMBER OF CONTOUR LEVELS',NCL)
	DO 10, I=1,NCL
	   CALL CPSETI('PAI - PARAMETER ARRAY INDEX',I)
	   CALL CPSETI('CLU - CONTOUR LEVEL USE FLAG',3)
 10	CONTINUE

C Set RWC so that labels come out on some lines
	CALL CPSETI('RWC - REAL WORKSPACE FOR CONTOURS',125)

C Draw Contours
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
	    Z(I,J)= 10.E-8*(-16.*REAL(I**2*J) +
     +		    34*REAL(I*J**2) - REAL(6*I) + 93.)
  20	  CONTINUE
  10	CONTINUE

	RETURN
	END

