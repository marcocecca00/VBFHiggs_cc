ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP()

      IMPLICIT NONE
      DOUBLE PRECISION PI, ZERO
      LOGICAL READLHA
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'model_functions.inc'
      REAL*16 MP__PI, MP__ZERO
      PARAMETER (MP__PI=3.1415926535897932384626433832795E0_16)
      PARAMETER (MP__ZERO=0E0_16)
      INCLUDE 'mp_input.inc'
      INCLUDE 'mp_coupl.inc'

      LOGICAL UPDATELOOP
      COMMON /TO_UPDATELOOP/UPDATELOOP
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      READLHA = .TRUE.
      INCLUDE 'intparam_definition.inc'
      INCLUDE 'mp_intparam_definition.inc'

      CALL COUP1()
C     
couplings needed to be evaluated points by points
C     
      CALL COUP2()
      CALL MP_COUP2()

      RETURN
      END

      SUBROUTINE UPDATE_AS_PARAM()

      IMPLICIT NONE
      DOUBLE PRECISION PI, ZERO
      LOGICAL READLHA, FIRST
      DATA FIRST /.TRUE./
      SAVE FIRST
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      LOGICAL UPDATELOOP
      COMMON /TO_UPDATELOOP/UPDATELOOP
      INCLUDE 'model_functions.inc'
      DOUBLE PRECISION GOTHER

      DOUBLE PRECISION MODEL_SCALE
      COMMON /MODEL_SCALE/MODEL_SCALE


      INCLUDE '../maxparticles.inc'
      INCLUDE '../cuts.inc'
      INCLUDE '../run.inc'

      DOUBLE PRECISION ALPHAS
      EXTERNAL ALPHAS

      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      READLHA = .FALSE.

      INCLUDE 'intparam_definition.inc'



C     
couplings needed to be evaluated points by points
C     
      CALL COUP2()

      RETURN
      END

      SUBROUTINE UPDATE_AS_PARAM2(MU_R2,AS2)

      IMPLICIT NONE
      DOUBLE PRECISION PI
      PARAMETER  (PI=3.141592653589793D0)
      DOUBLE PRECISION MU_R2, AS2
      INCLUDE 'model_functions.inc'
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      DOUBLE PRECISION MODEL_SCALE
      COMMON /MODEL_SCALE/MODEL_SCALE


      IF (MU_R2.GT.0D0) MU_R = DSQRT(MU_R2)
      MODEL_SCALE = DSQRT(MU_R2)
      G = SQRT(4.0D0*PI*AS2)
      AS = AS2

      CALL UPDATE_AS_PARAM()


      RETURN
      END

      SUBROUTINE MP_UPDATE_AS_PARAM()

      IMPLICIT NONE
      LOGICAL READLHA
      INCLUDE 'model_functions.inc'
      REAL*16 MP__PI, MP__ZERO
      PARAMETER (MP__PI=3.1415926535897932384626433832795E0_16)
      PARAMETER (MP__ZERO=0E0_16)
      INCLUDE 'mp_input.inc'
      INCLUDE 'mp_coupl.inc'

      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      INCLUDE 'actualize_mp_ext_params.inc'
      READLHA = .FALSE.
      INCLUDE 'mp_intparam_definition.inc'


C     
couplings needed to be evaluated points by points
C     
      CALL MP_COUP2()

      RETURN
      END

