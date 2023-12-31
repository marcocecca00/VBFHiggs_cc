ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP1()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      REAL*16 MP__PI, MP__ZERO
      PARAMETER (MP__PI=3.1415926535897932384626433832795E0_16)
      PARAMETER (MP__ZERO=0E0_16)
      INCLUDE 'mp_input.inc'
      INCLUDE 'mp_coupl.inc'

      GC_33 = -((MDL_COMPLEXI*MDL_YB)/MDL_SQRT__2)
      MP__GC_33 = -((MP__MDL_COMPLEXI*MP__MDL_YB)/MP__MDL_SQRT__2)
      GC_34 = -((MDL_COMPLEXI*MDL_YC)/MDL_SQRT__2)
      MP__GC_34 = -((MP__MDL_COMPLEXI*MP__MDL_YC)/MP__MDL_SQRT__2)
      GC_37 = -((MDL_COMPLEXI*MDL_YT)/MDL_SQRT__2)
      MP__GC_37 = -((MP__MDL_COMPLEXI*MP__MDL_YT)/MP__MDL_SQRT__2)
      END
