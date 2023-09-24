C     Set of subroutines to project the Feynman diagrama amplitudes
C      AMPL onto the color flow space.

      SUBROUTINE ML5_0_0_1_COMPUTE_COLOR_FLOWS(HEL_MULT,DO_CUMULATIVE)
      IMPLICIT NONE
      LOGICAL DO_CUMULATIVE
      INTEGER HEL_MULT
      CALL ML5_0_0_1_REINITIALIZE_JAMPS()
      CALL ML5_0_0_1_DO_COMPUTE_COLOR_FLOWS(.FALSE.)
      END SUBROUTINE

      SUBROUTINE ML5_0_0_1_COMPUTE_COLOR_FLOWS_DERIVED_QUANTITIES(HEL_M
     $ULT)
      IMPLICIT NONE
      INTEGER HEL_MULT
      CALL ML5_0_0_1_COMPUTE_AMPL2(HEL_MULT)
      CALL ML5_0_0_1_COMPUTE_JAMP2(HEL_MULT)
      CONTINUE
      END SUBROUTINE

      SUBROUTINE ML5_0_0_1_DEALLOCATE_COLOR_FLOWS()
      IMPLICIT NONE
      CALL ML5_0_0_1_DO_COMPUTE_COLOR_FLOWS(.TRUE.)
      END SUBROUTINE

      SUBROUTINE ML5_0_0_1_DO_COMPUTE_COLOR_FLOWS(CLEANUP)
      IMPLICIT NONE
C     
C     CONSTANTS 
C     
      CHARACTER*512 PROC_PREFIX
      PARAMETER ( PROC_PREFIX='ML5_0_0_1_')
      CHARACTER*512 LOOPCOLORFLOWCOEFSNAME
      PARAMETER ( LOOPCOLORFLOWCOEFSNAME='LoopColorFlowCoefs.dat')
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(0.0D0,1.0D0))
      INTEGER NCONFIGS
      PARAMETER (NCONFIGS=1)

      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=9)
      INTEGER    NSQUAREDSO, NLOOPAMPSO
      PARAMETER (NSQUAREDSO=1, NLOOPAMPSO=1)
      INTEGER NLOOPFLOWS
      PARAMETER (NLOOPFLOWS=1)
C     
C     LOCAL VARIABLES 
C     
C     When this subroutine is called with CLEANUP=True, it deallocates
C      all its
C     arrays
      LOGICAL CLEANUP
C     Storing concatenated filenames
      CHARACTER*512 TMP
      CHARACTER*512 LOOPCOLORFLOWCOEFSN

      INTEGER I, J, K, SOINDEX, ARRAY_SIZE
      COMPLEX*16 PROJ_COEF

      INTEGER CONFIG_INDEX
      INTEGER CONFIG_MAP(NLOOPAMPS)
      DATA (CONFIG_MAP(I),I=     1,     6) /     1,     1,     1,    
     $  1,     1,     1/
      DATA (CONFIG_MAP(I),I=     7,     9) /     1,     1,     1/

C     
C     FUNCTIONS
C     
      INTEGER ML5_0_0_1_ML5SOINDEX_FOR_BORN_AMP
      INTEGER ML5_0_0_1_ML5SOINDEX_FOR_LOOP_AMP
C     
C     GLOBAL VARIABLES
C     
      CHARACTER(512) MLPATH
      COMMON/MLPATH/MLPATH

      COMPLEX*16 AMPL(3,NLOOPAMPS)
      COMMON/ML5_0_0_1_AMPL/AMPL
      COMPLEX*16 JAMPL(3,NLOOPFLOWS,NLOOPAMPSO)
      COMMON/ML5_0_0_1_JAMPL/JAMPL
      COMPLEX*16 JAMPL_FOR_AMP2(3,NLOOPFLOWS,NLOOPAMPSO,NCONFIGS)
      COMMON/ML5_0_0_1_JAMPL_FOR_AMP2/JAMPL_FOR_AMP2


C     Now a more advanced data structure for storing the projection
C      coefficient
C     This is of course not F77 standard but widely supported by now.
      TYPE PROJCOEFFS
        INTEGER, DIMENSION(:), ALLOCATABLE :: NUM
        INTEGER, DIMENSION(:), ALLOCATABLE :: DENOM
        INTEGER, DIMENSION(:), ALLOCATABLE :: AMPID
      ENDTYPE PROJCOEFFS
      TYPE(PROJCOEFFS), DIMENSION(NLOOPFLOWS), SAVE ::
     $  LOOPCOLORPROJECTOR

C     ----------
C     BEGIN CODE
C     ----------

C     CleanUp duties, i.e. array deallocation
      IF (CLEANUP) THEN
        DO I=1,NLOOPFLOWS
          IF(ALLOCATED(LOOPCOLORPROJECTOR(I)%NUM))
     $      DEALLOCATE(LOOPCOLORPROJECTOR(I)%NUM)
          IF(ALLOCATED(LOOPCOLORPROJECTOR(I)%DENOM))
     $      DEALLOCATE(LOOPCOLORPROJECTOR(I)%DENOM)
          IF(ALLOCATED(LOOPCOLORPROJECTOR(I)%AMPID))
     $      DEALLOCATE(LOOPCOLORPROJECTOR(I)%AMPID)
        ENDDO
        RETURN
      ENDIF

C     Initialization; must allocate array from data files. All arrays
C      are allocated at once, so we only need to check if
C      LoopColorProjector(0)%Num is allocated and all other arrays
C      must share the same status.
      IF(.NOT.ALLOCATED(LOOPCOLORPROJECTOR(1)%NUM)) THEN
        CALL JOINPATH(MLPATH,PROC_PREFIX,TMP)
        CALL JOINPATH(TMP,LOOPCOLORFLOWCOEFSNAME,LOOPCOLORFLOWCOEFSN)

C       Initialize the LoopColorProjector
        OPEN(1, FILE=LOOPCOLORFLOWCOEFSN, ERR=201, STATUS='OLD',      
     $        ACTION='READ')
        DO I=1,NLOOPFLOWS
          READ(1,*,END=998) ARRAY_SIZE
          ALLOCATE(LOOPCOLORPROJECTOR(I)%NUM(ARRAY_SIZE))
          ALLOCATE(LOOPCOLORPROJECTOR(I)%DENOM(ARRAY_SIZE))
          ALLOCATE(LOOPCOLORPROJECTOR(I)%AMPID(ARRAY_SIZE))
          READ(1,*,END=998) (LOOPCOLORPROJECTOR(I)%NUM(J),J=1
     $     ,ARRAY_SIZE)
          READ(1,*,END=998) (LOOPCOLORPROJECTOR(I)%DENOM(J),J=1
     $     ,ARRAY_SIZE)
          READ(1,*,END=998) (LOOPCOLORPROJECTOR(I)%AMPID(J),J=1
     $     ,ARRAY_SIZE)
        ENDDO
        GOTO 203
 201    CONTINUE
        STOP 'Color projection coefficients could not be initialized'
     $   //' from file ML5_0_0_1_LoopColorFlowCoefs.dat.'
 203    CONTINUE
        CLOSE(1)

        GOTO 999
 998    CONTINUE
        STOP 'End of file reached. Should not have happened.'
 999    CONTINUE

      ENDIF

C     Projection of the loop amplitudes
      DO I=1,NLOOPFLOWS
        DO J=1,SIZE(LOOPCOLORPROJECTOR(I)%AMPID)
          SOINDEX = ML5_0_0_1_ML5SOINDEX_FOR_LOOP_AMP(LOOPCOLORPROJECTO
     $R(I)%AMPID(J))
          PROJ_COEF=DCMPLX(LOOPCOLORPROJECTOR(I)%NUM(J)
     $     /DBLE(ABS(LOOPCOLORPROJECTOR(I)%DENOM(J))),0.0D0)
          IF(LOOPCOLORPROJECTOR(I)%DENOM(J).LT.0) PROJ_COEF=PROJ_COEF
     $     *IMAG1
          DO K=1,3
            JAMPL(K,I,SOINDEX) = JAMPL(K,I,SOINDEX) + PROJ_COEF*AMPL(K
     $       ,LOOPCOLORPROJECTOR(I)%AMPID(J))
            CONFIG_INDEX = CONFIG_MAP(LOOPCOLORPROJECTOR(I)%AMPID(J))
            IF (CONFIG_INDEX.NE.0) THEN
              JAMPL_FOR_AMP2(K,I,SOINDEX,CONFIG_INDEX) =
     $          JAMPL_FOR_AMP2(K,I,SOINDEX,CONFIG_INDEX) + PROJ_COEF
     $         *AMPL(K,LOOPCOLORPROJECTOR(I)%AMPID(J))
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE

C     This subroutine initializes the color flow matrix.
      SUBROUTINE ML5_0_0_1_INITIALIZE_FLOW_COLORMATRIX()
      IMPLICIT NONE
C     
C     CONSTANTS 
C     
      CHARACTER*512 PROC_PREFIX
      PARAMETER ( PROC_PREFIX='ML5_0_0_1_')
      CHARACTER*512 LOOPCOLORFLOWMATRIXNAME
      PARAMETER ( LOOPCOLORFLOWMATRIXNAME='LoopColorFlowMatrix.dat')

      INTEGER NLOOPFLOWS
      PARAMETER (NLOOPFLOWS=1)
C     
C     LOCAL VARIABLES 
C     
C     Storing concatenated filenames
      CHARACTER*512 TMP
      CHARACTER*512 LOOPCOLORFLOWMATRIXN
      INTEGER I, J
C     
C     GLOBAL VARIABLES
C     
      CHARACTER(512) MLPATH
      COMMON/MLPATH/MLPATH

C     Now a more advanced data structure for storing the projection
C      coefficient
C     This is of course not F77 standard but widely supported by now.
      TYPE COLORCOEFF
        SEQUENCE
        INTEGER :: NUM
        INTEGER :: DENOM
      ENDTYPE COLORCOEFF
      TYPE(COLORCOEFF), DIMENSION(NLOOPFLOWS,NLOOPFLOWS) ::
     $  LOOPCOLORFLOWMATRIX
      LOGICAL CMINITIALIZED
      DATA CMINITIALIZED/.FALSE./
      COMMON/ML5_0_0_1_FLOW_COLOR_MATRIX/LOOPCOLORFLOWMATRIX,
     $  CMINITIALIZED

C     ----------
C     BEGIN CODE
C     ----------

C     Initialization
      IF(.NOT.CMINITIALIZED) THEN
        CMINITIALIZED = .TRUE.
        CALL JOINPATH(MLPATH,PROC_PREFIX,TMP)
        CALL JOINPATH(TMP,LOOPCOLORFLOWMATRIXNAME,LOOPCOLORFLOWMATRIXN)
C       Initialize the LoopColorFlowMatrix
        OPEN(1, FILE=LOOPCOLORFLOWMATRIXN, ERR=701, STATUS='OLD',     
     $         ACTION='READ')
        DO I=1,NLOOPFLOWS
          READ(1,*,END=898) (LOOPCOLORFLOWMATRIX(I,J)%NUM,J=1
     $     ,NLOOPFLOWS)
          READ(1,*,END=898) (LOOPCOLORFLOWMATRIX(I,J)%DENOM,J=1
     $     ,NLOOPFLOWS)
        ENDDO
        GOTO 703
 701    CONTINUE
        STOP 'Color factors could not be initialized from file'
     $   //' ML5_0_0_1_BornColorFlowMatrix.dat.'
 703    CONTINUE
        CLOSE(1)

        GOTO 899
 898    CONTINUE
        STOP 'End of file reached. Should not have happened.'
 899    CONTINUE

      ENDIF
      END SUBROUTINE

C     For the loop-inducedd MadEvent this subroutine computes the
C      squared of 
C     each individual JAMPL to created the usual JAMP2 array used for
C      the
C     color assignation of the events.
      SUBROUTINE ML5_0_0_1_COMPUTE_JAMP2(HEL_MULT)
      IMPLICIT NONE
C     
C     CONSTANTS 
C     
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(0.0D0,1.0D0))

      INTEGER    NSQUAREDSO, NLOOPAMPSO
      PARAMETER (NSQUAREDSO=1, NLOOPAMPSO=1)
      INTEGER NLOOPFLOWS
      PARAMETER (NLOOPFLOWS=1)
C     
C     ARGUMENT 
C     
C     HEL_MULT is the helicity multiplier which can be more than one
C      if several helicity configuration are mapped onto the one being
C      currently computed.
      INTEGER HEL_MULT
C     
C     LOCAL VARIABLES 
C     
      INTEGER I, J, M, N, K, ISQSO, DOUBLEFACT
      COMPLEX*16 COLOR_COEF
C     
C     FUNCTIONS
C     
      INTEGER ML5_0_0_1_ML5SQSOINDEX
C     
C     GLOBAL VARIABLES
C     
      LOGICAL CHOSEN_SO_CONFIGS(NSQUAREDSO)
      COMMON/ML5_0_0_1_CHOSEN_LOOP_SQSO/CHOSEN_SO_CONFIGS

      COMPLEX*16 JAMPL(3,NLOOPFLOWS,NLOOPAMPSO)
      COMMON/ML5_0_0_1_JAMPL/JAMPL

      REAL*8 JAMP2(NLOOPFLOWS)
      COMMON/ML5_0_0_1_JAMP2/JAMP2

C     Now a more advanced data structure for storing the projection
C      coefficient
C     This is of course not F77 standard but widely supported by now.
      TYPE COLORCOEFF
        SEQUENCE
        INTEGER :: NUM
        INTEGER :: DENOM
      ENDTYPE COLORCOEFF
      TYPE(COLORCOEFF), DIMENSION(NLOOPFLOWS,NLOOPFLOWS) ::
     $  LOOPCOLORFLOWMATRIX
      LOGICAL CMINITIALIZED
      COMMON/ML5_0_0_1_FLOW_COLOR_MATRIX/LOOPCOLORFLOWMATRIX,
     $  CMINITIALIZED

C     ----------
C     BEGIN CODE
C     ----------

C     Initialization
      IF(.NOT.CMINITIALIZED) THEN
        CALL ML5_0_0_1_INITIALIZE_FLOW_COLORMATRIX()
      ENDIF

C     Compute the Loop ME from the loop color flow amplitudes (JAMPL)
      DO I=1,NLOOPFLOWS
        COLOR_COEF=DCMPLX(LOOPCOLORFLOWMATRIX(I,I)%NUM
     $   /DBLE(ABS(LOOPCOLORFLOWMATRIX(I,I)%DENOM)),0.0D0)
        IF (LOOPCOLORFLOWMATRIX(I,I)%DENOM.LT.0) COLOR_COEF=COLOR_COEF
     $   *IMAG1
        DO M=1,NLOOPAMPSO
          DO N=M,NLOOPAMPSO
            ISQSO = ML5_0_0_1_ML5SQSOINDEX(M,N)
            IF(CHOSEN_SO_CONFIGS(ISQSO)) THEN
              IF(M.NE.N) THEN
                DOUBLEFACT=2
              ELSE
                DOUBLEFACT=1
              ENDIF
              JAMP2(I)=JAMP2(I)+DABS(DOUBLEFACT*HEL_MULT
     $         *DBLE(COLOR_COEF*(JAMPL(1,I,M)*DCONJG(JAMPL(1,I,N)))))
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE

C     This routine is used as a crosscheck only to make sure the loop
C      ME computed
C     from the JAMP is equal to the amplitude computed directly. It is
C      a consistency
C     check of the color projection and computation.
      SUBROUTINE ML5_0_0_1_COMPUTE_RES_FROM_JAMP(RES,HEL_MULT)
      IMPLICIT NONE
C     
C     CONSTANTS 
C     
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(0.0D0,1.0D0))

      INTEGER    NSQUAREDSO, NLOOPAMPSO, NSO
      PARAMETER (NSQUAREDSO=1, NLOOPAMPSO=1, NSO=1)
      INTEGER NLOOPFLOWS
      PARAMETER (NLOOPFLOWS=1)
      INTEGER NCONFIGS
      PARAMETER (NCONFIGS=1)
C     
C     ARGUMENT 
C     
      REAL*8 RES(0:3,0:NSQUAREDSO)
C     HEL_MULT is the helicity multiplier which can be more than one
C      if several helicity configuration are mapped onto the one being
C      currently computed.
      INTEGER HEL_MULT
C     
C     LOCAL VARIABLES 
C     
      INTEGER I, J, M, N, K, ISQSO, DOUBLEFACT
      INTEGER LOWERBOUND
      INTEGER CONFIG_I, CONFIG_J
      INTEGER ORDERS_A(NSO), ORDERS_B(NSO)
      COMPLEX*16 COLOR_COEF
      REAL*8 TEMP(3)
      REAL*8 TEMP_AMP2(3,NCONFIGS)
C     
C     FUNCTIONS
C     
C     This function belongs to the loop ME computation (loop_matrix.f)
C      and is prefixed with ML5 
      INTEGER ML5_0_0_1_ML5SQSOINDEX
C     This function belongs to the Born ME computation (born_matrix.f)
C      and is not prefixed at all
      INTEGER ML5_0_0_1_SQSOINDEX, ML5_0_0_1_SOINDEX_FOR_AMPORDERS
C     
C     GLOBAL VARIABLES
C     
      CHARACTER(512) MLPATH
      COMMON/MLPATH/MLPATH
      INTEGER SQSO_TARGET
      COMMON/ML5_0_0_1_SOCHOICE/SQSO_TARGET

      COMPLEX*16 JAMPL(3,NLOOPFLOWS,NLOOPAMPSO)
      COMMON/ML5_0_0_1_JAMPL/JAMPL
      REAL*8 AMP2_ALL(3,NCONFIGS,0:NSQUAREDSO)
      COMMON/ML5_0_0_1_AMP2_ALL/AMP2_ALL
      COMPLEX*16 JAMPL_FOR_AMP2(3,NLOOPFLOWS,NLOOPAMPSO,NCONFIGS)
      COMMON/ML5_0_0_1_JAMPL_FOR_AMP2/JAMPL_FOR_AMP2

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/ML5_0_0_1_SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

C     Now a more advanced data structure for storing the projection
C      coefficient
C     This is of course not F77 standard but widely supported by now.
      TYPE COLORCOEFF
        SEQUENCE
        INTEGER :: NUM
        INTEGER :: DENOM
      ENDTYPE COLORCOEFF
      TYPE(COLORCOEFF), DIMENSION(NLOOPFLOWS,NLOOPFLOWS) ::
     $  LOOPCOLORFLOWMATRIX
      LOGICAL CMINITIALIZED
      COMMON/ML5_0_0_1_FLOW_COLOR_MATRIX/LOOPCOLORFLOWMATRIX,
     $  CMINITIALIZED

C     ----------
C     BEGIN CODE
C     ----------
C     Always reinitialize AMP2_ALL because it is not a cumulative
C      quantity
      CALL ML5_0_0_1_REINITIALIZE_AMP2_ALL()

C     Initialization
      IF(.NOT.CMINITIALIZED) THEN
        CALL ML5_0_0_1_INITIALIZE_FLOW_COLORMATRIX()
      ENDIF

      DO I=0,NSQUAREDSO
        DO K=0,3
          RES(K,I)=0.0D0
        ENDDO
      ENDDO


C     Compute the Loop ME from the loop color flow amplitudes (JAMPL)
      DO I=1,NLOOPFLOWS
        DO J=I,NLOOPFLOWS
          COLOR_COEF=DCMPLX(LOOPCOLORFLOWMATRIX(I,J)%NUM
     $     /DBLE(ABS(LOOPCOLORFLOWMATRIX(I,J)%DENOM)),0.0D0)
          IF (LOOPCOLORFLOWMATRIX(I,J)%DENOM.LT.0)
     $      COLOR_COEF=COLOR_COEF*IMAG1
          DO M=1,NLOOPAMPSO
            IF (I.EQ.J) THEN
              LOWERBOUND = M
            ELSE
              LOWERBOUND = 1
            ENDIF
            DO N=LOWERBOUND,NLOOPAMPSO
              ISQSO = ML5_0_0_1_ML5SQSOINDEX(M,N)
              IF(J.NE.I) THEN
                DOUBLEFACT=2
              ELSEIF (M.NE.N) THEN
                DOUBLEFACT=2
              ELSE
                DOUBLEFACT=1
              ENDIF
              TEMP(1) = DOUBLEFACT*HEL_MULT*DBLE(COLOR_COEF*(JAMPL(1,I
     $         ,M)*DCONJG(JAMPL(1,J,N))))
C             Computing the quantities below is not strictly necessary
C              since the result should be finite
C             It is however a good cross-check.
              TEMP(2) = DOUBLEFACT*HEL_MULT*DBLE(COLOR_COEF*(JAMPL(2,I
     $         ,M)*DCONJG(JAMPL(1,J,N)) + JAMPL(1,I,M)*DCONJG(JAMPL(2
     $         ,J,N))))
              TEMP(3) = DOUBLEFACT*HEL_MULT*DBLE(COLOR_COEF*(JAMPL(3,I
     $         ,M)*DCONJG(JAMPL(1,J,N)) + JAMPL(1,I,M)*DCONJG(JAMPL(3
     $         ,J,N))+JAMPL(2,I,M)*DCONJG(JAMPL(2,J,N))))
              DO CONFIG_I=1,NCONFIGS
                DO K=1,3
                  TEMP_AMP2(K,CONFIG_I) = 0.0D0
                ENDDO
C               For the multi-channeling, to choose to square against
C                only among a given configuration, use the do loop
C                below
                DO CONFIG_J=CONFIG_I,CONFIG_I
C                 For the multi-channeling, to choose to square
C                  against all other configuration, use the do loop
C                  below
C                 DO config_j=1,nconfigs
                  TEMP_AMP2(1,CONFIG_I) = TEMP_AMP2(1,CONFIG_I) +
     $              DOUBLEFACT*HEL_MULT*DBLE(COLOR_COEF
     $             *(JAMPL_FOR_AMP2(1,I,M,CONFIG_I)
     $             *DCONJG(JAMPL_FOR_AMP2(1,J,N,CONFIG_J))))
                  TEMP_AMP2(2,CONFIG_I) = TEMP_AMP2(2,CONFIG_I) +
     $              DOUBLEFACT*HEL_MULT*DBLE(COLOR_COEF
     $             *(JAMPL_FOR_AMP2(2,I,M,CONFIG_I)
     $             *DCONJG(JAMPL_FOR_AMP2(1,J,N,CONFIG_J)) +
     $              JAMPL_FOR_AMP2(1,I,M,CONFIG_I)
     $             *DCONJG(JAMPL_FOR_AMP2(2,J,N,CONFIG_J))))
                  TEMP_AMP2(3,CONFIG_I) = TEMP_AMP2(3,CONFIG_I) +
     $              DOUBLEFACT*HEL_MULT*DBLE(COLOR_COEF
     $             *(JAMPL_FOR_AMP2(3,I,M,CONFIG_I)
     $             *DCONJG(JAMPL_FOR_AMP2(1,J,N,CONFIG_J)) +
     $              JAMPL_FOR_AMP2(1,I,M,CONFIG_I)
     $             *DCONJG(JAMPL_FOR_AMP2(3,J,N,CONFIG_J))
     $             +JAMPL_FOR_AMP2(2,I,M,CONFIG_I)
     $             *DCONJG(JAMPL_FOR_AMP2(2,J,N,CONFIG_J))))

                ENDDO
              ENDDO
              DO K=1,3
                RES(K,ISQSO) = RES(K,ISQSO) + TEMP(K)
                DO CONFIG_I=1,NCONFIGS
                  AMP2_ALL(K,CONFIG_I,ISQSO) = AMP2_ALL(K,CONFIG_I
     $             ,ISQSO) + TEMP_AMP2(K,CONFIG_I)
                ENDDO
              ENDDO
              IF((.NOT.FILTER_SO).OR.SQSO_TARGET.EQ.
     $         -1.OR.SQSO_TARGET.EQ.ISQSO) THEN
                DO K=1,3
                  RES(K,0) = RES(K,0) + TEMP(K)
C                 Notice here that the squared order constraints will
C                  apply to the determination of the AMP2 for
C                  multichanneling. Place it outside the corresponding
C                  statement if this is not desirable.
                  DO CONFIG_I=1,NCONFIGS
                    AMP2_ALL(K,CONFIG_I,0) = AMP2_ALL(K,CONFIG_I,0) +
     $                TEMP_AMP2(K,CONFIG_I)
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE

      SUBROUTINE ML5_0_0_1_COMPUTE_AMPL2(HEL_MULT)
C     
C     CONSTANTS 
C     
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=9)
      INTEGER    NLOOPDIAGRAMS
      PARAMETER (NLOOPDIAGRAMS=6)
      INTEGER    NCONFIGS
      PARAMETER (NCONFIGS=1)
      INTEGER    NSQUAREDSO, NLOOPAMPSO, NSO
      PARAMETER (NSQUAREDSO=1, NLOOPAMPSO=1, NSO=1)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J, SOINDEX
      INTEGER CONFIG_INDEX_MAP(NCONFIGS)
      DATA (CONFIG_INDEX_MAP(I),I=     1,     1) /     1/
C     
C     ARGUMENTS
C     
      INTEGER HEL_MULT
C     
C     GLOBAL VARIABLES
C     
      REAL*8 AMP2(NLOOPDIAGRAMS)
      COMMON/ML5_0_0_1_AMP2/AMP2
      REAL*8  AMP2_ALL(3,NCONFIGS,0:NSQUAREDSO)
      COMMON/ML5_0_0_1_AMP2_ALL/AMP2_ALL
      LOGICAL CHOSEN_SO_CONFIGS(NSQUAREDSO)
      COMMON/ML5_0_0_1_CHOSEN_LOOP_SQSO/CHOSEN_SO_CONFIGS
C     ----------
C     BEGIN CODE
C     ----------

      DO I=1,NCONFIGS
C       The configurations are not necessarily consecutive when
C        grouping subprocesses. So, here, we simply skip those which
C        have no diagrams for this process.
        IF (CONFIG_INDEX_MAP(I).EQ.0) THEN
          CYCLE
        ENDIF
        DO J=1,NSQUAREDSO
C         For the multi-channeling, it is in general better to ignore
C          the squared order constraints (which, if they are any,
C          would restrict the AMP2 quantities to be coming from the
C          interference terms only).  
          IF (.TRUE.) THEN
C           Chosing the if statement below instead would as a result
C            to apply the squared order constraints to the AMP2
C            quantities given to MadEvent for the multichanneling.
C           IF (CHOSEN_SO_CONFIGS(I)) THEN
            AMP2(CONFIG_INDEX_MAP(I))=AMP2(CONFIG_INDEX_MAP(I))
     $       +DABS(AMP2_ALL(1,I,J))
          ENDIF
        ENDDO
      ENDDO

      END SUBROUTINE

C     This subroutine resets to 0 the common arrays JAMPL, JAMPB and
C      possibly JAMPL_FOR_AMP2 if used
      SUBROUTINE ML5_0_0_1_REINITIALIZE_JAMPS()
      IMPLICIT NONE
C     
C     CONSTANTS 
C     
      COMPLEX*16 CMPLXZERO
      PARAMETER (CMPLXZERO=(0.0D0,0.0D0))
      REAL*8 ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER    NLOOPAMPSO
      PARAMETER (NLOOPAMPSO=1)
      INTEGER NLOOPFLOWS
      PARAMETER (NLOOPFLOWS=1)
      INTEGER NCONFIGS
      PARAMETER (NCONFIGS=1)
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K,L
C     
C     GLOBAL VARIABLES
C     
      COMPLEX*16 JAMPL(3,NLOOPFLOWS,NLOOPAMPSO)
      COMMON/ML5_0_0_1_JAMPL/JAMPL
      COMPLEX*16 JAMPL_FOR_AMP2(3,NLOOPFLOWS,NLOOPAMPSO,NCONFIGS)
      COMMON/ML5_0_0_1_JAMPL_FOR_AMP2/JAMPL_FOR_AMP2

C     ----------
C     BEGIN CODE
C     ----------

      DO I=1,NLOOPAMPSO
        DO J=1,NLOOPFLOWS
          DO K=1,3
            JAMPL(K,J,I)=CMPLXZERO
          ENDDO
        ENDDO
      ENDDO
      DO L=1,NCONFIGS
        DO I=1,NLOOPAMPSO
          DO J=1,NLOOPFLOWS
            DO K=1,3
              JAMPL_FOR_AMP2(K,J,I,L)=CMPLXZERO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE

C     This subroutine resets to 0 the common arrays AMP2_ALL
      SUBROUTINE ML5_0_0_1_REINITIALIZE_AMP2_ALL()
      IMPLICIT NONE
C     
C     CONSTANTS 
C     
      INTEGER    NSQUAREDSO, NLOOPAMPSO, NSO
      PARAMETER (NSQUAREDSO=1, NLOOPAMPSO=1, NSO=1)
      REAL*8 ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER NCONFIGS
      PARAMETER (NCONFIGS=1)
C     
C     LOCAL VARIABLES
C     
      INTEGER I,K,L
C     
C     GLOBAL VARIABLES
C     
      REAL*8 AMP2_ALL(3,NCONFIGS,0:NSQUAREDSO)
      COMMON/ML5_0_0_1_AMP2_ALL/AMP2_ALL
C     ----------
C     BEGIN CODE
C     ----------

      DO I=1,NCONFIGS
        DO K=1,3
          DO L=0,NSQUAREDSO
            AMP2_ALL(K,I,L)=ZERO
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE

C     This subroutine resets all cumulative (for each helicity) arrays
C      in common block and deriving from the color flow computation.
C     This subroutine is called by loop_matrix.f at appropriate time
C      when these arrays must be reset.
C     An example of this is the JAMP2 which must be cumulatively
C      computed as loop_matrix.f loops over helicity configs but must
C      be resets when loop_matrix.f starts over with helicity one (for
C      the stability test for example).
      SUBROUTINE ML5_0_0_1_REINITIALIZE_CUMULATIVE_ARRAYS()
      IMPLICIT NONE
C     
C     CONSTANTS 
C     
      COMPLEX*16 CMPLXZERO
      PARAMETER (CMPLXZERO=(0.0D0,0.0D0))
      REAL*8 ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER    NLOOPDIAGRAMS
      PARAMETER (NLOOPDIAGRAMS=6)
      INTEGER NLOOPFLOWS
      PARAMETER (NLOOPFLOWS=1)
      INTEGER NCONFIGS
      PARAMETER (NCONFIGS=1)
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
C     
C     GLOBAL VARIABLES
C     
      REAL*8 AMP2(NLOOPDIAGRAMS)
      COMMON/ML5_0_0_1_AMP2/AMP2
      REAL*8 JAMP2(NLOOPFLOWS)
      COMMON/ML5_0_0_1_JAMP2/JAMP2

C     ----------
C     BEGIN CODE
C     ----------

      DO I=1,NLOOPDIAGRAMS
        AMP2(I)=ZERO
      ENDDO
      DO I=1,NLOOPFLOWS
        JAMP2(I)=ZERO
      ENDDO

      CONTINUE

      END SUBROUTINE

