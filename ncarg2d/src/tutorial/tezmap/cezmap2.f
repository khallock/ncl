C
C   $Id: cezmap2.f,v 1.3 1994-07-11 14:10:38 haley Exp $
C
C
C Define error file, Fortran unit number, and workstation type,
C and workstation ID.
C
      PARAMETER (IERRF=6, LUNIT=2, IWTYPE=SED_WSTYPE, IWKID=1)
C
C Open GKS, Turn Clipping off
C
      CALL GOPKS (IERRF, ISZDM)
      CALL GOPWK (IWKID, LUNIT, IWTYPE)
      CALL GACWK (IWKID)
C
C Invoke demo driver
C
      CALL CEZMAP('SV',40.,-50.,0.,'PO')
C
C Deactivate and close workstation, close GKS.
C
      CALL GDAWK (IWKID)
      CALL GCLWK (IWKID)
      CALL GCLKS

      STOP
      END

      SUBROUTINE CEZMAP(PROJ,PLAT,PLON,ROTA,OUTLN)

      CHARACTER*2 PROJ, OUTLN
      REAL PLIM1(2), PLIM2(2), PLIM3(2), PLIM4(2)

      DATA PLIM1 /0.,0./
      DATA PLIM2 /0.,0./
      DATA PLIM3 /0.,0./
      DATA PLIM4 /0.,0./
C
C CMPLOT demonstrates MAPLOT drawing continental and political outlines
C
C Set up Maps.
C
C
C Set Continental, political outline options
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
      CALL MAPSET ('MA',PLIM1,PLIM2,PLIM3,PLIM4)
C
C Draw map
C
      CALL MAPDRW
C
C Advance the frame.
C
      CALL FRAME
C
C Done.
C
      RETURN
      END
