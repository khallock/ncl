	PROGRAM CCPLLP

        PARAMETER (NRAN=30,LRWK=3500,LIWK=4000,LMAP=50000)
        PARAMETER (MREG=50,NREG=50)
	REAL XRAN(NRAN), YRAN(NRAN), ZRAN(NRAN)
	REAL XREG(MREG), YREG(NREG), ZREG(MREG,NREG), RWRK(LRWK)
	INTEGER IWRK(LIWK), MAP(LMAP)

	EXTERNAL CPDRPL

	DATA XRAN /12., 60., 14., 33.,  8., 12., 43., 57., 22., 15.,
     1		   19., 12., 64., 19., 15., 55., 31., 32., 33., 29.,
     2		   18.,  1., 18., 42., 56.,  9.,  6., 12., 44., 19./
	DATA YRAN / 1.,  2.,  3., 53.,  7., 11., 13., 17., 19., 49.,
     1		    1., 31., 37.,  5.,  7., 47., 61., 17.,  5., 23.,
     2		   29.,  3.,  5., 41., 43.,  9., 13., 59.,  1., 67./
	DATA ZRAN /1.0, 1.5, 1.7, 1.4, 1.9, 1.0, 1.5, 1.2, 1.8, 1.4,
     1		   1.8, 1.7, 1.9, 1.5, 1.2, 1.1, 1.3, 1.7, 1.2, 1.6,
     2		   1.9, 1.0, 1.6, 1.3, 1.4, 1.8, 1.7, 1.5, 1.1, 1.0/


C Open GKS
	CALL OPNGKS
	CALL GSCLIP(0)
C
C  Find the min and max data values.
C
      XMIN = 0.0
      XMAX = 65.0
      YMIN =  0.0
      YMAX = 68.0
C
C Choose the X and Y coordinates for interpolation points on the 
C regular grid.
C
      DO 101 I=1,MREG
        XREG(I)=XMIN + (XMAX - XMIN)* REAL(I-1)/MREG
  101 CONTINUE
C
      DO 102 I=1,NREG
        YREG(I)=YMIN + (YMAX - YMIN)* REAL(I-1)/NREG
  102 CONTINUE

C Interpolate data onto a regular grid
	CALL IDSFFT (1,NRAN,XRAN,YRAN,ZRAN,
     +		MREG,NREG,MREG,XREG,YREG,ZREG,IWRK,RWRK)

C ------------Default Labels-----------------
C Set up viewport
	CALL CPSETR('VPT - VIEWPORT TOP',0.95)
	CALL CPSETR('VPB - VIEWPORT BOTTOM',.50)
	CALL CPSETR('VPL - VIEWPORT LEFT',0.0)
	CALL CPSETR('VPR - VIEWPORT RIGHT',.48)
	CALL CPSETC('ILT - INFORMATION LABEL TEXT',' ')
C Initialize Conpack
	CALL CPRECT(ZREG, MREG, MREG, NREG, RWRK, LRWK, IWRK, LIWK)
C Force Conpack to chose contour levels
	CALL CPPKCL(ZREG, RWRK, IWRK)
C Modify Conpack chosen parameters
	CALL CPGETI('NCL - NUMBER OF CONTOUR LEVELS',NCONS)
C choose which labelling scheme will be used.
	CALL CPSETI('LLP - LINE LABEL POSTIONING FLAG',1)
C Turn off high and low labels
	CALL CPSETC('HLT - HIGH/LOW LABEL TEXT',' '' ')
	DO 11, I=1,NCONS
	   CALL CPSETI('PAI - PARAMETER ARRAY INDEX',I)
C Force every line to be labeled.
	   CALL CPSETI('CLU - CONTOUR LEVEL USE FLAG',3)
 11	CONTINUE
C Draw Perimeter
	CALL CPBACK(ZREG, RWRK, IWRK)
C Draw Contours
	CALL CPCLDR(ZREG,RWRK,IWRK)
	CALL CPLBDR(ZREG,RWRK,IWRK)
C Title plot
	CALL GSELNT(0)
	CALL PLCHHQ(.25, .975, 'Default Labels', .017, 0., 0.)

C ------------Medium Labels-----------------
C Set up Areas for drawing boxes
	CALL ARINAM(MAP, LMAP)
C Set up viewport
	CALL CPSETR('VPT - VIEWPORT TOP',0.95)
	CALL CPSETR('VPB - VIEWPORT BOTTOM',.50)
	CALL CPSETR('VPR - VIEWPORT RIGHT',1.0)
	CALL CPSETR('VPL - VIEWPORT LEFT',.52)
	CALL CPSETC('ILT - INFORMATION LABEL TEXT',' ')
C Initialize Conpack
	CALL CPRECT(ZREG, MREG, MREG, NREG, RWRK, LRWK, IWRK, LIWK)
C Force Conpack to chose contour levels
	CALL CPPKCL(ZREG, RWRK, IWRK)
C Modify Conpack chosen parameters
	CALL CPGETI('NCL - NUMBER OF CONTOUR LEVELS',NCONS)
C choose which labelling scheme will be used.
	CALL CPSETI('LLP - LINE LABEL POSTIONING FLAG',2)
C Turn off high and low labels
	CALL CPSETC('HLT - HIGH/LOW LABEL TEXT',' '' ')
	DO 12, I=1,NCONS
	   CALL CPSETI('PAI - PARAMETER ARRAY INDEX',I)
C Force every line to be labeled.
	   CALL CPSETI('CLU - CONTOUR LEVEL USE FLAG',3)
 12	CONTINUE
C Draw Perimeter
	CALL CPBACK(ZREG, RWRK, IWRK)
C Add contours to area map
        CALL CPCLAM(ZREG, RWRK, IWRK, MAP)
C Add label boxes to area map
        CALL CPLBAM(ZREG, RWRK, IWRK, MAP)
C Draw contours
	CALL CPCLDM(ZREG, RWRK, IWRK, MAP, CPDRPL)
C Draw labels
	CALL CPLBDR(ZREG, RWRK, IWRK)
C Title plot
	CALL GSELNT(0)
	CALL PLCHHQ(.75, .975, 'Regular Scheme', .017, 0., 0.)

C ------------High Quality Labels-----------------
C Set up area map
	CALL ARINAM(MAP, LMAP)
C Set up viewport
	CALL CPSETR('VPT - VIEWPORT TOP',.45)
	CALL CPSETR('VPB - VIEWPORT BOTTOM',0.0)
	CALL CPSETR('VPL - VIEWPORT LEFT',0.0)
	CALL CPSETR('VPR - VIEWPORT RIGHT',1.0)
	CALL CPSETC('ILT - INFORMATION LABEL TEXT',' ')
C Initialize Conpack
	CALL CPRECT(ZREG, MREG, MREG, NREG, RWRK, LRWK, IWRK, LIWK)
C Force Conpack to chose contour levels
	CALL CPPKCL(ZREG, RWRK, IWRK)
C Modify Conpack chosen parameters
	CALL CPGETI('NCL - NUMBER OF CONTOUR LEVELS',NCONS)
C choose which labelling scheme will be used.
	CALL CPSETI('LLP - LINE LABEL POSTIONING FLAG',3)
C Turn off high and low labels
	CALL CPSETC('HLT - HIGH/LOW LABEL TEXT',' '' ')
	DO 13, I=1,NCONS
	   CALL CPSETI('PAI - PARAMETER ARRAY INDEX',I)
C Force every line to be labeled.
	   CALL CPSETI('CLU - CONTOUR LEVEL USE FLAG',3)
 13	CONTINUE
C Draw Perimeter
	CALL CPBACK(ZREG, RWRK, IWRK)
C Add contours to area map
        CALL CPCLAM(ZREG, RWRK, IWRK, MAP)
C Add label boxes to area map
        CALL CPLBAM(ZREG, RWRK, IWRK, MAP)
C Draw contours
	CALL CPCLDM(ZREG, RWRK, IWRK, MAP, CPDRPL)
C Draw labels
	CALL CPLBDR(ZREG, RWRK, IWRK)
C Title plot
	CALL GSELNT(0)
	CALL PLCHHQ(.5, .475, 'Penalty Scheme', .017, 0., 0.)

C Close frame and close GKS
	CALL FRAME
	CALL CLSGKS


	WRITE (6,*) 'AREA MAP SIZE =',MAP(1) - MAP(6) + MAP(5)
	STOP
	END

