C
C	$Id: creb15.f,v 1.2 2000-07-11 23:11:16 haley Exp $
C                                                                      
C                Copyright (C)  2000
C        University Corporation for Atmospheric Research
C                All Rights Reserved
C
C This file is free software; you can redistribute it and/or modify
C it under the terms of the GNU Lesser General Public License as
C published by the Free Software Foundation; either version 2.1 of the
C License, or (at your option) any later version.
C
C This software is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C Lesser General Public License for more details.
C
C You should have received a copy of the GNU Lesser General Public
C License along with this software; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
C USA.
C

      SUBROUTINE CREB15(INUNIT,IOUTUN,INNUM,ITEMP,LENTEM)
C
C ROUTINE TO STORE A SEQUENCE OF INTEGERS INTO AN ARRAY.
C
C ON ENTRY
C     INUNIT IS A FILE WHICH CONTAINS INNUM CARD IMAGES. EACH CARD
C     IMAGE CONTAINS 16 INTEGERS, EACH IN AN I5 FORMAT. EACH INTEGER
C     REPRESENTS A POSITIVE 15 BIT VALUE.
C     INUNIT CAN BE READ BY CREB15.
C     IOUTUN IS A UNIT NUMBER WHERE CREB15 CAN WRITE A RECORD.
C     IOUTUN IS ASSUMED TO BE POSITIONED RIGHT.
C     ITEMP IS AN ARRAY OF LENGTH LENTEM. ITS ELEMENTS MAY HAVE
C     ANY VALUE.
C     LENTEM MUST BE (INNUM*16-1)/NUM15+1 OR BIGGER, WHERE NUM15 IS
C     THE NUMBER OF 15 BIT UNITS WHICH FIT AS A WHOLE INTO 1 WORD.
C     INNUM IS THE NUMBER OF CARD IMAGES ON FILE INUNIT.
C ON EXIT
C     THE VALUES OF INUNIT,IOUTUN, AND INNUM ARE UNCHANGED.
C     FILE INUNIT IS REWOUND BUT OTHERWISE UNCHANGED.
C     ITEMP CONTAINS THE INTEGERS FROM FILE INUNIT STORED LEFT
C     JUSTIFIED AS 15 BIT UNITS WITH AS MANY UNITS PER WORD AS
C     POSSIBLE WITHOUT CROSSING WORD BOUNDARIES.
C     A BINARY RECORD IS WRITTEN ON UNIT IOUTUN WHICH CONTAINS THE FIRST
C     NUMOUT WORDS OF THE ARRAY ITEMP AS ONE LONG BIT STRING,
C     WHERE NUMOUT IS EXACTLY (INNUM*16-1)/NUM15+1 AND NUM15 DEFINED
C     AS ABOVE.
C CALLS
C     IAND,IOR,ISHIFT
C
      DIMENSION ITEMP(LENTEM)
C
C TEMPORARY STORAGE TO CONTAIN 1 CARD IMAGE.
      DIMENSION ICARD(16)
C
C SEE BLOCK DATA DPORT FOR MEANING OF CONSTANTS IN COMMON BLOCK.
      COMMON /IDC1/ NBWD, IZERO, MA15
C
C
C LOCAL VARIABLES
C   IWORD - THE ENTRY IN THE ARRAY ITEMP WHICH IS CURRENTLY BEING
C           FILLED
C   IN15 - CONTAINS 15 BIT UNIT, RIGHT JUSTIFIED.
C   IIN15 - CONTAINS THE 15 BIT UNIT - OR PART OF IT - AT THE POSITION
C           WHICH IT WILL HAVE IN ITEMP(IWORD)
C   IPOS - THE NUMBER OF 15 BIT UNITS ALREADY STORED IN THE CURRENT
C          WORD + 1 .
C   ICARD - TEMPORARY STORAGE FOR 1 CARD IMAGE.
C   NUM15 - THE NUMBER OF 15 BIT UNITS WHICH FIT AS A WHOLE INTO 1 WORD.
C
C
C THE NUMBER OF 15 BIT UNITS WHICH FIT AS A WHOLE INTO ONE WORD.
C
      NUM15 = NBWD/15
C
      REWIND INUNIT
C
C
C
      IPOS = 1
      IWORD = 1
      ITEMP(1) = IZERO
C
C     READ THE INTEGERS FROM FILE INUNIT AND STORE THEM INTO THE
C     ARRAY ITEMP AS WHOLE UNITS, LEFT JUSTIFIED.
C
C LOOP THROUGH ALL CARD IMAGES.
      DO 1 I = 1,INNUM
      READ(INUNIT,100)(ICARD(J),J=1,16)
C LOOP THROUG ALL 15 BIT UNITS ON A CARD IMAGE.
      DO 10 J=1,16
      IN15 = ICARD(J)
C STORE 15 BIT UNIT AT CORRECT POSITION IN WORD.
      IIN15 = ISHIFT(IAND(IN15,MA15),NBWD-IPOS*15)
      ITEMP(IWORD) = IOR(IIN15,ITEMP(IWORD))
      IPOS = IPOS + 1
      IF (IPOS .LE. NUM15) GOTO 11
C START A NEW WORD.
      IPOS = 1
      IWORD = IWORD + 1
      ITEMP(IWORD) = IZERO
   11 CONTINUE
   10 CONTINUE
    1 CONTINUE
C
C
C WRITE 1 RECORD IN A BINARY FILE.
C
      WRITE(IOUTUN) (ITEMP(J),J=1,IWORD)
C
      REWIND INUNIT
C
  100 FORMAT(16I5)
C
      RETURN
      END
