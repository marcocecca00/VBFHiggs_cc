C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Identity(2,1)
C     
      SUBROUTINE MP_FFS1_3(F1, F2, COUP, M3, W3,S3)
      IMPLICIT NONE
      COMPLEX*32 CI
      PARAMETER (CI=(0Q0,1Q0))
      COMPLEX*32 COUP
      COMPLEX*32 F1(*)
      COMPLEX*32 F2(*)
      REAL*16 M3
      COMPLEX*32 P3(0:3)
      COMPLEX*32 S3(5)
      COMPLEX*32 TMP0
      REAL*16 W3
      COMPLEX*32 DENOM
      S3(1) = +F1(1)+F2(1)
      S3(2) = +F1(2)+F2(2)
      S3(3) = +F1(3)+F2(3)
      S3(4) = +F1(4)+F2(4)
      P3(0) = -S3(1)
      P3(1) = -S3(2)
      P3(2) = -S3(3)
      P3(3) = -S3(4)
      TMP0 = (F1(5)*F2(5)+F1(6)*F2(6)+F1(7)*F2(7)+F1(8)*F2(8))
      DENOM = COUP/(P3(0)**2-P3(1)**2-P3(2)**2-P3(3)**2 - M3 * (M3 -CI
     $ * W3))
      S3(5)= DENOM*CI * TMP0
      END


