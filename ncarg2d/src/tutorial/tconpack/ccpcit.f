	PROGRAM CCPCIT

        PARAMETER (K=40,N=40,LRWK=1000,LIWK=1000)
	REAL Z(K,N), RWRK(LRWK)
	INTEGER M, IWRK(LIWK)
	DIMENSION CIT(10),LIT(10)

	DATA CIT /1.,2.,3.,4.,5.,6.,7.,8.,9.,0./
	DATA LIT /2, 2, 2, 2, 2, 2, 2, 2, 2, 0 /

	CALL GETDAT (Z, K, M, N) 

C Open GKS
	CALL OPNGKS
C Change nice values to be steps of 1/3. (1/3, 2/3, 3/3...)
C Draw labels at every 5th contour level no matter which contour
C level interval is chosen.
	DO 101 I=1,10
	   CALL CPSETI ('PAI - PARAMETER ARRAY INDEX',I)
	   CALL CPSETR ('CIT - CONTOUR INTERVAL TABLE',CIT(I))
	   CALL CPSETI ('LIT - LABEL INTERVAL TABLE',LIT(I))
  101	CONTINUE

C Initialize Conpack
	CALL CPRECT(Z,K,M,N,RWRK,LRWK,IWRK,LIWK)
C Draw perimeter
	CALL CPBACK(Z, RWRK, IWRK)
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
	    Z(I,J)= 10.E-5*(-16.*REAL(I**2*J) +
     +		    34*REAL(I*J**2) - REAL(6*I) + 93.)
  20	  CONTINUE
  10	CONTINUE

	RETURN
	END

