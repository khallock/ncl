C
C	$Id: cmptra.f,v 1.3 1994-07-11 14:10:53 haley Exp $
C
C
C Define error file, Fortran unit number, and workstation type,
C and workstation ID.
C
      PARAMETER (IERRF=6, LUNIT=2, IWTYPE=SED_WSTYPE, IWKID=1)

      REAL PLIM1(2), PLIM2(2), PLIM3(2), PLIM4(2)

      DATA PLIM1 /0.,0./
      DATA PLIM2 /0.,0./
      DATA PLIM3 /0.,0./
      DATA PLIM4 /0.,0./
C
C Open GKS, Turn Clipping off
C
      CALL GOPKS (IERRF, ISZDM)
      CALL GOPWK (IWKID, LUNIT, IWTYPE)
      CALL GACWK (IWKID)
C
C Invoke demo driver
C
      CALL CMPTRA('OR',35.,-105.,0.,'PO','CO',22.,-120.,47.,-65.)
C
C Deactivate and close workstation, close GKS.
C
      CALL GDAWK (IWKID)
      CALL GCLWK (IWKID)
      CALL GCLKS

      STOP
      END

      SUBROUTINE CMPTRA(PROJ, PLAT, PLON, ROTA, OUTLN,
     +     JLIM, PLIM1, PLIM2, PLIM3, PLIM4)

      CHARACTER*2 PROJ, OUTLN, JLIM
      REAL PLIM1(2), PLIM2(2), PLIM3(2), PLIM4(2)
C
C CMPTRA demonstrates marking points on a map
C
C Draw Continental, political outlines 
C
      CALL MAPSTC ('OU - OUTLINE DATASET SELECTOR',OUTLN)
C
C Set up projection
C
      CALL MAPROJ (PROJ,PLAT,PLON,ROTA)
C
C If it's a satellite projection, choose a satellite distance
C
      IF (PROJ.EQ.'SV') CALL MAPSTR ('SA - SATELLITE DISTANCE',5.)
C
C Set limits of map
C
      CALL MAPSET (JLIM,PLIM1,PLIM2,PLIM3,PLIM4)
C
C Turn off Grid lines
C
      CALL MAPSTR ('GR',0.)
C
C Draw map
C
      CALL MAPDRW
C
C Draw a star over Boulder Colorado
C
      CALL MAPTRA(40.,-105.15,X,Y)
      IF (X.NE.1.E12) CALL POINTS (X, Y, 1, -3, 0)
C
C Draw the state of Colorado in
C
      CALL MAPIT (37.,-109.,0)
      CALL MAPIT (41.,-109.,1)
      CALL MAPIT (41.,-102.,1)
      CALL MAPIT (37.,-102.,1)
      CALL MAPIT (37.,-109.,1)
      CALL MAPIQ
C
C Advance the frame.
C
      CALL FRAME
C
C Done.
C
      RETURN
      END
