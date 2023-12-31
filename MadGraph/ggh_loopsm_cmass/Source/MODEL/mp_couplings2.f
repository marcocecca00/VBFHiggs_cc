ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE MP_COUP2()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'
      REAL*16 MP__PI, MP__ZERO
      PARAMETER (MP__PI=3.1415926535897932384626433832795E0_16)
      PARAMETER (MP__ZERO=0E0_16)
      INCLUDE 'mp_input.inc'
      INCLUDE 'mp_coupl.inc'

      MP__GC_5 = MP__MDL_COMPLEXI*MP__G
      MP__R2_GGHC = 4.000000E+00_16*(-((MP__MDL_COMPLEXI*MP__MDL_YC)
     $ /MP__MDL_SQRT__2))*(1.000000E+00_16/2.000000E+00_16)
     $ *(MP__MDL_G__EXP__2/(8.000000E+00_16*MP__PI**2))*MP__MDL_MC
      MP__R2_GGHB = 4.000000E+00_16*(-((MP__MDL_COMPLEXI*MP__MDL_YB)
     $ /MP__MDL_SQRT__2))*(1.000000E+00_16/2.000000E+00_16)
     $ *(MP__MDL_G__EXP__2/(8.000000E+00_16*MP__PI**2))*MP__MDL_MB
      MP__R2_GGHT = 4.000000E+00_16*(-((MP__MDL_COMPLEXI*MP__MDL_YT)
     $ /MP__MDL_SQRT__2))*(1.000000E+00_16/2.000000E+00_16)
     $ *(MP__MDL_G__EXP__2/(8.000000E+00_16*MP__PI**2))*MP__MDL_MT
      END
