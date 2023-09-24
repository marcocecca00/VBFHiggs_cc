      SUBROUTINE SMATRIX12(P,ANS)
C     
C     Generated by MadGraph5_aMC@NLO v. 3.5.0, 2023-05-12
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     MadGraph5_aMC@NLO for Madevent Version
C     
C     Returns amplitude squared summed/avg over colors
C     and helicities
C     for the point in phase space P(0:3,NEXTERNAL)
C     
C     Process: d d~ > h d d~ QCD=0 @1
C     Process: s s~ > h s s~ QCD=0 @1
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      Include 'genps.inc'
      Include 'maxconfigs.inc'
      Include 'nexternal.inc'
      Include 'maxamps.inc'
      INTEGER                 NCOMB         
      PARAMETER (             NCOMB=6)
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=2) 
      INTEGER    NDIAGS
      PARAMETER (NDIAGS=2) 
      INTEGER    THEL
      PARAMETER (THEL=2*NCOMB)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS
C     
C     LOCAL VARIABLES 
C     
      INTEGER CONFSUB(MAXSPROC,LMAXCONFIGS)
      INCLUDE 'config_subproc_map.inc'
      INTEGER NHEL(NEXTERNAL,NCOMB)
      REAL*8 T
      REAL*8 R,SUMHEL,TS(NCOMB)
      INTEGER I,IDEN
      INTEGER JC(NEXTERNAL),II
      REAL*8 XTOT
      INTEGER  J, JJ

      double precision get_channel_cut
      external get_channel_cut

C     
C     GLOBAL VARIABLES
C     
      DOUBLE PRECISION AMP2(MAXAMPS), JAMP2(0:MAXFLOW)
      COMMON/TO_AMPS/  AMP2,       JAMP2


C     
C     INFORMATION TO WRITE THE HELICITY IN THE EVENT --not memory
C      efficient--
C     
      CHARACTER*101         HEL_BUFF
      COMMON/TO_HELICITY/  HEL_BUFF

      INTEGER NB_SPIN_STATE_in(2)
      common /nb_hel_state/ nb_spin_state_in

      REAL*8 POL(2)

      COMMON/TO_POLARIZATION/ POL
      double precision tmin_for_channel
      integer sde_strat    !  1 means standard single diagram enhancement strategy,
C     2 means approximation by the	denominator of the propagator
      common/TO_CHANNEL_STRAT/tmin_for_channel,	sde_strat

      INTEGER          ISUM_HEL
      LOGICAL                    MULTI_CHANNEL
      COMMON/TO_MATRIX/ISUM_HEL, MULTI_CHANNEL
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      COMMON/TO_MCONFIGS/MAPCONFIG, ICONFIG
      INTEGER SUBDIAG(MAXSPROC),IB(2)
      COMMON/TO_SUB_DIAG/SUBDIAG,IB
      
      DATA (NHEL(I,1),I=1,5) / 1,-1, 0,-1, 1/
      DATA (NHEL(I,2),I=1,5) / 1,-1, 0, 1,-1/
      DATA (NHEL(I,3),I=1,5) / 1, 1, 0, 1, 1/
      DATA (NHEL(I,4),I=1,5) /-1,-1, 0,-1,-1/
      DATA (NHEL(I,5),I=1,5) /-1, 1, 0,-1, 1/
      DATA (NHEL(I,6),I=1,5) /-1, 1, 0, 1,-1/
      DATA IDEN/36/

C     ----------
C     BEGIN CODE
C     ----------

      DO I=1,NEXTERNAL
        JC(I) = +1
      ENDDO

      IF (multi_channel) THEN
        DO I=1,NDIAGS
          AMP2(I)=0D0
        ENDDO
        JAMP2(0)=2
        DO I=1,INT(JAMP2(0))
          JAMP2(I)=0D0
        ENDDO
      ENDIF
      ANS = 0D0
      WRITE(HEL_BUFF,'(20I5)') (0,I=1,NEXTERNAL)
C     Kiran please check if you need this:    
      DO I=1,NCOMB
        TS(I)=0d0
      ENDDO

      call MATRIX12(P ,JC(1), TS)
      DO I=1,NCOMB     
        T=TS(I)  
        DO JJ=1,nincoming
          IF(POL(JJ).NE.1d0.AND.NHEL(JJ,I).EQ.INT(SIGN(1d0,POL(JJ))))
     $      THEN
            T=T*ABS(POL(JJ))*NB_SPIN_STATE_IN(JJ)/2d0    !  NB_SPIN_STATE(JJ)/2d0 is added for polarised beam
          ELSE IF(POL(JJ).NE.1d0)THEN
            T=T*(2d0-ABS(POL(JJ)))*NB_SPIN_STATE_IN(JJ)/2d0
          ENDIF
        ENDDO
        ANS=ANS+DABS(T)
        TS(I)=T
      ENDDO

      IF (ANS.ne.0d0) THEN
        CALL RANMAR(R)
        SUMHEL=0d0
        DO I=1,NCOMB
          SUMHEL=SUMHEL+DABS(TS(I))/ANS
          IF(R.LT.SUMHEL)THEN
            WRITE(HEL_BUFF,'(20i5)')(NHEL(II,I),II=1,NEXTERNAL)
C           Set right sign for ANS, based on sign of chosen helicity
            ANS=DSIGN(ANS,TS(I))
            GOTO 10
          ENDIF
        ENDDO
 10     CONTINUE   
      ENDIF
      IF (MULTI_CHANNEL) THEN
        XTOT=0D0
        DO I=1,LMAXCONFIGS
          J = CONFSUB(12, I)
          if (J.ne.0)then
            if (sde_strat.eq.1)then
              AMP2(J) = AMP2(J) * GET_CHANNEL_CUT(P, I)
            else
              AMP2(J) = GET_CHANNEL_CUT(P, I)
            endif
            XTOT=XTOT+AMP2(J)

          endif
        ENDDO
        IF (XTOT.NE.0D0) THEN
          ANS=ANS*AMP2(SUBDIAG(12))/XTOT
        ELSE IF(ANS.ne.0d0) THEN
          write(*,*) 'Problem in the multi-channeling. All amp2 are'
     $     //' zero but not the total matrix-element'
          stop 1
        ENDIF
      ENDIF
      ANS=ANS/DBLE(IDEN)
      END


      Subroutine  MATRIX12(P,IC, TS)
C     
C     Generated by MadGraph5_aMC@NLO v. 3.5.0, 2023-05-12
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: d d~ > h d d~ QCD=0 @1
C     Process: s s~ > h s s~ QCD=0 @1
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=2) 
      include 'genps.inc'
      include 'nexternal.inc'
      include 'maxamps.inc'
      INTEGER    NWAVEFUNCS,     NCOLOR
      PARAMETER (NWAVEFUNCS=25, NCOLOR=2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(0D0,1D0))
      INTEGER NAMPSO, NSQAMPSO
      PARAMETER (NAMPSO=1, NSQAMPSO=1)
      LOGICAL CHOSEN_SO_CONFIGS(NSQAMPSO)
      DATA CHOSEN_SO_CONFIGS/.TRUE./
      SAVE CHOSEN_SO_CONFIGS
      INTEGER                 NCOMB         
      PARAMETER (             NCOMB=6)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      REAL*8 TS(NCOMB)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,M,N,K
      COMPLEX*16 ZTEMP,TMP_JAMP(0)
      COMPLEX*16 TMP(6)
      REAL*8 CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NCOMB,NGRAPHS), JAMP(NCOLOR,NAMPSO)
      COMPLEX*16 W(6,NWAVEFUNCS)
C     Needed for v4 models
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0d0, 0d0), (1d0, 0d0)/

      double precision fk_mdl_WZ 
      double precision fk_mdl_WH 
      double precision fk_ZERO 
      save fk_mdl_WZ 
      save fk_mdl_WH 
      save fk_ZERO 

      logical first
      data first /.true./
      save first
C     
C     FUNCTION
C     
      INTEGER SQSOINDEX12

C     
C     GLOBAL VARIABLES
C     
      Double Precision amp2(maxamps), jamp2(0:maxflow)
      common/to_amps/  amp2,       jamp2
      include 'coupl.inc'

      double precision tmin_for_channel
      integer sde_strat    !  1 means standard single diagram enhancement strategy,
C     2 means approximation by the	denominator of the propagator
      common/TO_CHANNEL_STRAT/tmin_for_channel,	sde_strat

      double precision small_width_treatment
      common/narrow_width/small_width_treatment
C     
C     COLOR DATA
C     
      DATA (CF(i,  1),i=  1,  2) /9.000000000000000d+00
     $ ,3.000000000000000d+00/
C     1 T(2,1) T(4,5)
      DATA (CF(i,  2),i=  1,  2) /3.000000000000000d+00
     $ ,9.000000000000000d+00/
C     1 T(2,5) T(4,1)
C     ----------
C     BEGIN CODE
C     ----------
      if (first) then
        first=.false.
        IF(ZERO.ne.0d0) fk_ZERO = SIGN(MAX(ABS(ZERO), ABS(ZERO
     $   *small_width_treatment)), ZERO)
        IF(mdl_WH.ne.0d0) fk_mdl_WH = SIGN(MAX(ABS(mdl_WH), ABS(mdl_MH
     $   *small_width_treatment)), mdl_WH)
        IF(mdl_WZ.ne.0d0) fk_mdl_WZ = SIGN(MAX(ABS(mdl_WZ), ABS(mdl_MZ
     $   *small_width_treatment)), mdl_WZ)
      endif
      AMP(:,:) = (0d0,0d0)
            CALL IXXXXX(P(0,1),ZERO,+1,+1*IC(1),W(1,1))  !  count 4
      CALL IXXXXX(P(0,1),ZERO,-1,+1*IC(1),W(1,2))  !  count 4
      CALL OXXXXX(P(0,2),ZERO,+1,-1*IC(2),W(1,3))  !  count 4
      CALL OXXXXX(P(0,2),ZERO,-1,-1*IC(2),W(1,4))  !  count 4
      CALL SXXXXX(P(0,3),+1*IC(3),W(1,5))  !  count 12
      CALL OXXXXX(P(0,4),ZERO,+1,+1*IC(4),W(1,6))  !  count 4
      CALL OXXXXX(P(0,4),ZERO,-1,+1*IC(4),W(1,7))  !  count 4
      CALL IXXXXX(P(0,5),ZERO,+1,-1*IC(5),W(1,8))  !  count 4
      CALL IXXXXX(P(0,5),ZERO,-1,-1*IC(5),W(1,9))  !  count 4
      CALL FFV2_3_3(W(1,1),W(1,3),GC_50,GC_58,MDL_MZ,FK_MDL_WZ,W(1,10))  !  count 1
      CALL FFV2_3_3(W(1,1),W(1,4),GC_50,GC_58,MDL_MZ,FK_MDL_WZ,W(1,11))  !  count 2
      CALL FFV2_3_3(W(1,2),W(1,3),GC_50,GC_58,MDL_MZ,FK_MDL_WZ,W(1,12))  !  count 2
      CALL FFV2_3_3(W(1,2),W(1,4),GC_50,GC_58,MDL_MZ,FK_MDL_WZ,W(1,13))  !  count 1
      CALL FFV2_3_3(W(1,8),W(1,6),GC_50,GC_58,MDL_MZ,FK_MDL_WZ,W(1,14))  !  count 1
      CALL FFV2_3_3(W(1,8),W(1,7),GC_50,GC_58,MDL_MZ,FK_MDL_WZ,W(1,15))  !  count 2
      CALL FFV2_3_3(W(1,9),W(1,6),GC_50,GC_58,MDL_MZ,FK_MDL_WZ,W(1,16))  !  count 2
      CALL FFV2_3_3(W(1,9),W(1,7),GC_50,GC_58,MDL_MZ,FK_MDL_WZ,W(1,17))  !  count 1
      CALL VVS1_0(W(1,10),W(1,14),W(1,5),GC_81,AMP(3,1))
      CALL VVS1P1N_1(W(1,15), W(1,5), GC_81, TMP(1))
      call CombineAmp(2,
     & (/1,5/), 
     & (/11,12/),
     & TMP, W, AMP(1,1))
      CALL VVS1P1N_1(W(1,16), W(1,5), GC_81, TMP(1))
      call CombineAmp(2,
     & (/2,6/), 
     & (/11,12/),
     & TMP, W, AMP(1,1))
      CALL VVS1_0(W(1,13),W(1,17),W(1,5),GC_81,AMP(4,1))  !  count 1
      CALL FFV2_3_3(W(1,1),W(1,6),GC_50,GC_58,MDL_MZ,ZERO,W(1,18))  !  count 2
      CALL FFV2_3_3(W(1,1),W(1,7),GC_50,GC_58,MDL_MZ,ZERO,W(1,19))  !  count 1
      CALL FFV2_3_3(W(1,2),W(1,6),GC_50,GC_58,MDL_MZ,ZERO,W(1,20))  !  count 1
      CALL FFV2_3_3(W(1,2),W(1,7),GC_50,GC_58,MDL_MZ,ZERO,W(1,21))  !  count 2
      CALL FFV2_3_3(W(1,8),W(1,3),GC_50,GC_58,MDL_MZ,ZERO,W(1,22))  !  count 2
      CALL FFV2_3_3(W(1,8),W(1,4),GC_50,GC_58,MDL_MZ,ZERO,W(1,23))  !  count 1
      CALL FFV2_3_3(W(1,9),W(1,3),GC_50,GC_58,MDL_MZ,ZERO,W(1,24))  !  count 1
      CALL FFV2_3_3(W(1,9),W(1,4),GC_50,GC_58,MDL_MZ,ZERO,W(1,25))  !  count 2
      CALL VVS1P1N_1(W(1,22), W(1,5), GC_81, TMP(1))
      call CombineAmp(2,
     & (/3,5/), 
     & (/18,21/),
     & TMP, W, AMP(1,2))
      CALL VVS1P1N_1(W(1,25), W(1,5), GC_81, TMP(1))
      call CombineAmp(2,
     & (/2,4/), 
     & (/18,21/),
     & TMP, W, AMP(1,2))
      CALL VVS1_0(W(1,19),W(1,23),W(1,5),GC_81,AMP(1,2))
      CALL VVS1_0(W(1,20),W(1,24),W(1,5),GC_81,AMP(6,2))  !  count 1

      JAMP(:,:)  = (0d0,0d0)
      DO K = 1, NCOMB
        
        JAMP(:,:) = (0D0,0D0)
C       JAMPs contributing to orders ALL_ORDERS=1
        JAMP(1,1) = (-1.000000000000000D+00)*AMP( K,1)
        JAMP(2,1) = AMP( K,2)

        TS(K) = 0.D0 
        DO M = 1, NAMPSO
          DO I = 1, NCOLOR
            ZTEMP = (0.D0,0.D0)
            DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J,M)
            ENDDO
            DO N = 1, NAMPSO

              TS(K) = TS(K) + ZTEMP*DCONJG(JAMP(I,N))

            ENDDO
          ENDDO
        ENDDO
        if(sde_strat.eq.1) then
        
        IF(SDE_STRAT.EQ.1)THEN
          AMP2(2)=AMP2(2)+AMP( K,2)*DCONJG(AMP( K,2))
          AMP2(1)=AMP2(1)+AMP( K,1)*DCONJG(AMP( K,1))
        ENDIF

        endif
        Do I = 1, NCOLOR
          DO M = 1, NAMPSO
            DO N = 1, NAMPSO

              Jamp2(i)=Jamp2(i)+DABS(DBLE(Jamp(i,m)*dconjg(Jamp(i,n))))

            enddo
          enddo
        Enddo
      ENDDO

      END


      SUBROUTINE PRINT_ZERO_AMP_12()

      integer i
      i =1
      return
      end
C     Set of functions to handle the array indices of the split orders


      INTEGER FUNCTION SQSOINDEX12(ORDERINDEXA, ORDERINDEXB)
C     
C     This functions plays the role of the interference matrix. It can
C      be hardcoded or 
C     made more elegant using hashtables if its execution speed ever
C      becomes a relevant
C     factor. From two split order indices, it return the
C      corresponding index in the squared 
C     order canonical ordering.
C     
C     CONSTANTS
C     

      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=1, NSQUAREDSO=1, NAMPSO=1)
C     
C     ARGUMENTS
C     
      INTEGER ORDERINDEXA, ORDERINDEXB
C     
C     LOCAL VARIABLES
C     
      INTEGER I, SQORDERS(NSO)
      INTEGER AMPSPLITORDERS(NAMPSO,NSO)
      DATA (AMPSPLITORDERS(  1,i),i=  1,  1) /    1/
      COMMON/AMPSPLITORDERS12/AMPSPLITORDERS
C     
C     FUNCTION
C     
      INTEGER SOINDEX_FOR_SQUARED_ORDERS12
C     
C     BEGIN CODE
C     
      DO I=1,NSO
        SQORDERS(I)=AMPSPLITORDERS(ORDERINDEXA,I)
     $   +AMPSPLITORDERS(ORDERINDEXB,I)
      ENDDO
      SQSOINDEX12=SOINDEX_FOR_SQUARED_ORDERS12(SQORDERS)
      END

      INTEGER FUNCTION SOINDEX_FOR_SQUARED_ORDERS12(ORDERS)
C     
C     This functions returns the integer index identifying the squared
C      split orders list passed in argument which corresponds to the
C      values of the following list of couplings (and in this order).
C     []
C     
C     CONSTANTS
C     
      INTEGER    NSO, NSQSO, NAMPSO
      PARAMETER (NSO=1, NSQSO=1, NAMPSO=1)
C     
C     ARGUMENTS
C     
      INTEGER ORDERS(NSO)
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J
      INTEGER SQSPLITORDERS(NSQSO,NSO)
      DATA (SQSPLITORDERS(  1,i),i=  1,  1) /    2/
      COMMON/SQPLITORDERS12/SQPLITORDERS
C     
C     BEGIN CODE
C     
      DO I=1,NSQSO
        DO J=1,NSO
          IF (ORDERS(J).NE.SQSPLITORDERS(I,J)) GOTO 1009
        ENDDO
        SOINDEX_FOR_SQUARED_ORDERS12 = I
        RETURN
 1009   CONTINUE
      ENDDO

      WRITE(*,*) 'ERROR:: Stopping in function' 
      WRITE(*,*) 'SOINDEX_FOR_SQUARED_ORDERS12'
      WRITE(*,*) 'Could not find squared orders ',(ORDERS(I),I=1,NSO)
      STOP

      END

      SUBROUTINE GET_NSQSO_BORN12(NSQSO)
C     
C     Simple subroutine returning the number of squared split order
C     contributions returned when calling smatrix_split_orders 
C     

      INTEGER    NSQUAREDSO
      PARAMETER  (NSQUAREDSO=1)

      INTEGER NSQSO

      NSQSO=NSQUAREDSO

      END

C     This is the inverse subroutine of SOINDEX_FOR_SQUARED_ORDERS.
C      Not directly useful, but provided nonetheless.
      SUBROUTINE GET_SQUARED_ORDERS_FOR_SOINDEX12(SOINDEX,ORDERS)
C     
C     This functions returns the orders identified by the squared
C      split order index in argument. Order values correspond to
C      following list of couplings (and in this order):
C     []
C     
C     CONSTANTS
C     
      INTEGER    NSO, NSQSO
      PARAMETER (NSO=1, NSQSO=1)
C     
C     ARGUMENTS
C     
      INTEGER SOINDEX, ORDERS(NSO)
C     
C     LOCAL VARIABLES
C     
      INTEGER I
      INTEGER SQPLITORDERS(NSQSO,NSO)
      COMMON/SQPLITORDERS12/SQPLITORDERS	  
C     
C     BEGIN CODE
C     
      IF (SOINDEX.gt.0.and.SOINDEX.le.NSQSO) THEN
        DO I=1,NSO
          ORDERS(I) =  SQPLITORDERS(SOINDEX,I)
        ENDDO
        RETURN
      ENDIF

      WRITE(*,*) 'ERROR:: Stopping function'
     $ //' GET_SQUARED_ORDERS_FOR_SOINDEX12'
      WRITE(*,*) 'Could not find squared orders index ',SOINDEX
      STOP

      END SUBROUTINE

C     This is the inverse subroutine of getting amplitude SO orders.
C      Not directly useful, but provided nonetheless.
      SUBROUTINE GET_ORDERS_FOR_AMPSOINDEX12(SOINDEX,ORDERS)
C     
C     This functions returns the orders identified by the split order
C      index in argument. Order values correspond to following list of
C      couplings (and in this order):
C     []
C     
C     CONSTANTS
C     
      INTEGER    NSO, NAMPSO
      PARAMETER (NSO=1, NAMPSO=1)
C     
C     ARGUMENTS
C     
      INTEGER SOINDEX, ORDERS(NSO)
C     
C     LOCAL VARIABLES
C     
      INTEGER I
      INTEGER AMPSPLITORDERS(NAMPSO,NSO)
      COMMON/AMPSPLITORDERS12/AMPSPLITORDERS
C     
C     BEGIN CODE
C     
      IF (SOINDEX.gt.0.and.SOINDEX.le.NAMPSO) THEN
        DO I=1,NSO
          ORDERS(I) =  AMPSPLITORDERS(SOINDEX,I)
        ENDDO
        RETURN
      ENDIF

      WRITE(*,*) 'ERROR:: Stopping function'
     $ //' GET_ORDERS_FOR_AMPSOINDEX12'
      WRITE(*,*) 'Could not find amplitude split orders index ',SOINDEX
      STOP

      END SUBROUTINE

C     This function is not directly useful, but included for
C      completeness
      INTEGER FUNCTION SOINDEX_FOR_AMPORDERS12(ORDERS)
C     
C     This functions returns the integer index identifying the
C      amplitude split orders passed in argument which correspond to
C      the values of the following list of couplings (and in this
C      order):
C     []
C     
C     CONSTANTS
C     
      INTEGER    NSO, NAMPSO
      PARAMETER (NSO=1, NAMPSO=1)
C     
C     ARGUMENTS
C     
      INTEGER ORDERS(NSO)
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J
      INTEGER AMPSPLITORDERS(NAMPSO,NSO)
      COMMON/AMPSPLITORDERS12/AMPSPLITORDERS
C     
C     BEGIN CODE
C     
      DO I=1,NAMPSO
        DO J=1,NSO
          IF (ORDERS(J).NE.AMPSPLITORDERS(I,J)) GOTO 1009
        ENDDO
        SOINDEX_FOR_AMPORDERS12 = I
        RETURN
 1009   CONTINUE
      ENDDO

      WRITE(*,*) 'ERROR:: Stopping function SOINDEX_FOR_AMPORDERS12'
      WRITE(*,*) 'Could not find squared orders ',(ORDERS(I),I=1,NSO)
      STOP

      END

