	PROGRAM CMPXY

	PARAMETER(JX=60,KX=26)
	PARAMETER(LRWK=5000,LIWK=5000)
C
	COMMON /HEIGHT/Z(JX,KX)

	REAL CNTR(JX,KX),RWRK(LRWK)
	INTEGER IWRK(LIWK)

C
C Read arrays containing data 
C
	OPEN(FILE='cpmpxy1.dat',STATUS='OLD',UNIT=13)
	DO 10, I=1,JX
          DO 20 J=1,KX
	    READ (13,*) Z(I,J)
 20       CONTINUE
 10	CONTINUE
	OPEN(FILE='cpmpxy2.dat',STATUS='OLD',UNIT=14)
	DO 30, I=1,JX
          DO 40 J=1,KX
	    READ (14,*) CNTR(I,J)
 40       CONTINUE
 30	CONTINUE
C
C Do a contour plot
C
	CALL OPNGKS
	CALL SET (.1,.95,.25,.85,-110.,-60.,1000.,0.,1)
	CALL CPSETI ('SET - DO-SET-CALL FLAG',0)
	CALL CPSETI ('MAP - MAPPING FLAG',4)
	CALL CPSETR ('XC1 - X COORDINATE AT INDEX 1',-110.)
	CALL CPSETR ('XCM - X COORDINATE AT INDEX M',-60.)
	CALL CPSETR ('SPV - SPECIAL VALUE',-9999.)
	CALL CPRECT (CNTR,JX,JX,KX,RWRK,LRWK,IWRK,LIWK)
	CALL CPCLDR (CNTR,RWRK,IWRK)
	CALL LABMOD ('(F5.0)','(F5.0)',0,0,0,0,0,0,0)
	CALL GRIDAL (5,0,10,0,1,1,5,-110.,1000.)
	CALL SET (0.,1.,0.,1.,0.,1.,0.,1.,1)
	CALL PLCHHQ (.03,.6,'PRESSURE',.012,90.,0.)
	CALL PLCHHQ (.5,.2,'LONGITUDE',.012,0.,0.)
	CALL CLSGKS

	STOP
	END
      SUBROUTINE CPMPXY(IMAP,XINP,YINP,XOTP,YOTP)
C
C Transform contours to overlay various mapping transformations:
C IMAP= 0 - Cartesian data: no transformation necessary
C IMAP= 1 - Lat/Lon transformation
C IMAP=-1 - inverse Lat/Lon transformation
C IMAP= 2 - Rho/Theta transformation
C IMAP=-2 - inverse Rho/Theta transformation
C IMAP= 3 - height in the X direction
C IMAP= 4 - Pressure in the X direction
C IMAP= 5 - height in the X direction
C IMAP= 6 - Pressure in the X direction
C
      PARAMETER(JX=60,KX=26)
C
      COMMON /HEIGHT/Z(JX,KX)
C
C Handle the EZMAP case ...
C
      IF (ABS(IMAP).EQ.1) THEN
        IF (IMAP.GT.0) THEN
          CALL MAPTRA (YINP,XINP,XOTP,YOTP)
        ELSE
          CALL MAPTRI (XINP,YINP,YOTP,XOTP)
        END IF
C
C ... the polar coordinate case ...
C
      ELSE IF (ABS(IMAP).EQ.2) THEN
        IF (IMAP.GT.0) THEN
          XOTP=XINP*COS(.017453292519943*YINP)
          YOTP=XINP*SIN(.017453292519943*YINP)
        ELSE
          XOTP=SQRT(XINP*XINP+YINP*YINP)
          YOTP=57.2957795130823*ATAN2(YINP,XINP)
        END IF
C
C Height & Pressure Data in the X direction
C
C Pressure transformation in the X direction
      ELSEIF(IMAP.EQ.3.OR.IMAP.EQ.4) THEN
C The height transformation in X direction is linear
        XOTP = XINP
C Find next lowest X data point
        X = XINP
C Distance between next lowest data point and contour point
        IIX=INT(X)
        DIFX=X-FLOAT(IIX)
C Find next lowest y data point
        Y = YINP
C Distance between next lowest data point and contour point
        IY=INT(Y)
        DIFY=Y-FLOAT(IY)
C Find next highest X and Y data points, 
C and make sure they are in the domain.
        IXP1 = MIN0(JX,IIX+1)
        IYP1 = MIN0(KX ,IY+1)
C Linear interpolation between points to give height at contour point
        Z1 = Z(IIX ,IY)+DIFY*(Z(IIX ,IYP1)-Z(IIX ,IY))
        Z2 = Z(IXP1,IY)+DIFY*(Z(IXP1,IYP1)-Z(IXP1,IY))
        ZR = Z1 + DIFX*(Z2-Z1)
        YOTP=ZR
C Pressure Data in the X direction
        IF(IMAP.EQ.4) YOTP = 1000. * EXP(-ZR/7.)
C
C If IMAP isn't specified as above, then do an identity transformation.
C
      ELSE
          XOTP = XINP
          YOTP = YINP
      ENDIF
C
      RETURN
      END
