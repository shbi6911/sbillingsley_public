C--------------------------------------------------------------------------
C      TEST DRIVER FOR GTD6 (ATMOSPHERIC MODEL)
      DIMENSION D(8,16),T(2,16),SW(25),APH(7)
      DIMENSION IDAY(15),UT(15),ALT(15),XLAT(15),XLONG(15),XLST(15),
     & F107A(15),F107(15),AP(15)    
      COMMON/GTS3C/DL(16)
      COMMON/DATIME/ISDATE(3),ISTIME(2),NAME(2)
      DATA IDAY/172,81,13*172/
      DATA UT/29000.,29000.,75000.,12*29000./
      DATA ALT/400.,400.,400.,100.,6*400.,0,10.,30.,50.,70./
      DATA XLAT/4*60.,0.,10*60./
      DATA XLONG/5*-70.,0.,9*-70./
      DATA XLST/6*16.,4.,8*16./
      DATA F107A/7*150.,70.,7*150./
      DATA F107/8*150.,180.,6*150./
      DATA AP/9*4.,40.,5*4./
      DATA APH/7*100./,SW/8*1.,-1.,16*1./
      DO I=1,15
         CALL GTD6(IDAY(I),UT(I),ALT(I),XLAT(I),XLONG(I),XLST(I),
     &             F107A(I),F107(I),AP(I),48,D(1,I),T(1,I))
         WRITE(6,100) (D(J,I),J=1,8),T(1,I),T(2,I),DL
      ENDDO
      CALL TSELEC(SW)
      I=16
      CALL GTD6(IDAY(1),UT(1),ALT(1),XLAT(1),XLONG(1),XLST(1),
     &             F107A(1),F107(1),APH,48,D(1,I),T(1,I))
      WRITE(6,100) (D(J,I),J=1,8),T(1,I),T(2,I),DL
      CALL GTD6(IDAY(1),UT(1),ALT(4),XLAT(1),XLONG(1),XLST(1),
     &             F107A(1),F107(1),APH,48,D(1,I),T(1,I))
      WRITE(6,100) (D(J,I),J=1,8),T(1,I),T(2,I),DL
      WRITE(6,300) NAME,ISDATE,ISTIME
      WRITE(6,200) (IDAY(I),I=1,5)
      WRITE(6,201) (UT(I),I=1,5)
      WRITE(6,202) (ALT(I),I=1,5)
      WRITE(6,203) (XLAT(I),I=1,5)
      WRITE(6,204) (XLONG(I),I=1,5)
      WRITE(6,205) (XLST(I),I=1,5)
      WRITE(6,206) (F107A(I),I=1,5)
      WRITE(6,207) (F107(I),I=1,5)
      WRITE(6,210) (T(1,I),I=1,5)
      WRITE(6,211) (T(2,I),I=1,5)
      WRITE(6,212) (D(1,I),I=1,5)
      WRITE(6,213) (D(2,I),I=1,5)
      WRITE(6,214) (D(3,I),I=1,5)
      WRITE(6,215) (D(4,I),I=1,5)
      WRITE(6,216) (D(5,I),I=1,5)
      WRITE(6,217) (D(7,I),I=1,5)
      WRITE(6,219) (D(8,I),I=1,5)
      WRITE(6,218) (D(6,I),I=1,5)
      WRITE(6,200) (IDAY(I),I=6,10)
      WRITE(6,201) (UT(I),I=6,10)
      WRITE(6,202) (ALT(I),I=6,10)
      WRITE(6,203) (XLAT(I),I=6,10)
      WRITE(6,204) (XLONG(I),I=6,10)
      WRITE(6,205) (XLST(I),I=6,10)
      WRITE(6,206) (F107A(I),I=6,10)
      WRITE(6,207) (F107(I),I=6,10)
      WRITE(6,210) (T(1,I),I=6,10)
      WRITE(6,211) (T(2,I),I=6,10)
      WRITE(6,212) (D(1,I),I=6,10)
      WRITE(6,213) (D(2,I),I=6,10)
      WRITE(6,214) (D(3,I),I=6,10)
      WRITE(6,215) (D(4,I),I=6,10)
      WRITE(6,216) (D(5,I),I=6,10)
      WRITE(6,217) (D(7,I),I=6,10)
      WRITE(6,219) (D(8,I),I=6,10)
      WRITE(6,218) (D(6,I),I=6,10)
      WRITE(6,200) (IDAY(I),I=11,15)
      WRITE(6,201) (UT(I),I=11,15)
      WRITE(6,202) (ALT(I),I=11,15)
      WRITE(6,203) (XLAT(I),I=11,15)
      WRITE(6,204) (XLONG(I),I=11,15)
      WRITE(6,205) (XLST(I),I=11,15)
      WRITE(6,206) (F107A(I),I=11,15)
      WRITE(6,207) (F107(I),I=11,15)
      WRITE(6,210) (T(1,I),I=11,15)
      WRITE(6,211) (T(2,I),I=11,15)
      WRITE(6,212) (D(1,I),I=11,15)
      WRITE(6,213) (D(2,I),I=11,15)
      WRITE(6,214) (D(3,I),I=11,15)
      WRITE(6,215) (D(4,I),I=11,15)
      WRITE(6,216) (D(5,I),I=11,15)
      WRITE(6,217) (D(7,I),I=11,15)
      WRITE(6,219) (D(8,I),I=11,15)
      WRITE(6,218) (D(6,I),I=11,15)
  100 FORMAT(1X,1P8E9.2/4X,2E10.3/4X,8E9.2/4X,8E9.2/)
  200 FORMAT(//' DAY  ',5I12)
  201 FORMAT(' UT   ',5F12.0)
  202 FORMAT(' ALT  ',5F12.0)
  203 FORMAT(' LAT  ',5F12.0)
  204 FORMAT(' LONG ',5F12.0)
  205 FORMAT(' LST  ',5F12.0)
  206 FORMAT(' F107A',5F12.0)
  207 FORMAT(' F107 ',5F12.0)
  210 FORMAT(/' TINF ',5F12.2)
  211 FORMAT(' TG   ',5F12.2)
  212 FORMAT(' HE   ',1P5E12.3)
  213 FORMAT(' O    ',1P5E12.3)
  214 FORMAT(' N2   ',1P5E12.3)
  215 FORMAT(' O2   ',1P5E12.3)
  216 FORMAT(' AR   ',1P5E12.3)
  217 FORMAT(' H    ',1P5E12.3)
  219 FORMAT(' N    ',1P5E12.3)
  218 FORMAT(' RHO  ',1P5E12.3)
  300 FORMAT(1X,2A4,2X,3A4,2X,2A4)
      STOP
      END
