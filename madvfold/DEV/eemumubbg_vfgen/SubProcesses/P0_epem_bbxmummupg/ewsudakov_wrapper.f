      SUBROUTINE SUDAKOV_WRAPPER(P)
      IMPLICIT NONE
      INCLUDE 'nexternal.inc'
      INCLUDE 'born_nhel.inc'
      DOUBLE PRECISION P(0:3, NEXTERNAL-1)

      INTEGER I,J
      INTEGER IHEL
      INTEGER     NCOMB
      PARAMETER ( NCOMB=  128 )
      INTEGER NHEL(NEXTERNAL-1,NCOMB)
      DATA (NHEL(I,   1),I=1,7) /-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,   2),I=1,7) /-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,   3),I=1,7) /-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,   4),I=1,7) /-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,   5),I=1,7) /-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,   6),I=1,7) /-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,   7),I=1,7) /-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,   8),I=1,7) /-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,   9),I=1,7) /-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  10),I=1,7) /-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  11),I=1,7) /-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  12),I=1,7) /-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  13),I=1,7) /-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  14),I=1,7) /-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  15),I=1,7) /-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  16),I=1,7) /-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  17),I=1,7) /-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  18),I=1,7) /-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  19),I=1,7) /-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  20),I=1,7) /-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  21),I=1,7) /-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  22),I=1,7) /-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  23),I=1,7) /-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  24),I=1,7) /-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  25),I=1,7) /-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  26),I=1,7) /-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  27),I=1,7) /-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  28),I=1,7) /-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  29),I=1,7) /-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  30),I=1,7) /-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  31),I=1,7) /-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  32),I=1,7) /-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  33),I=1,7) /-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  34),I=1,7) /-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  35),I=1,7) /-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  36),I=1,7) /-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  37),I=1,7) /-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  38),I=1,7) /-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  39),I=1,7) /-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  40),I=1,7) /-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  41),I=1,7) /-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  42),I=1,7) /-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  43),I=1,7) /-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  44),I=1,7) /-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  45),I=1,7) /-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  46),I=1,7) /-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  47),I=1,7) /-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  48),I=1,7) /-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  49),I=1,7) /-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  50),I=1,7) /-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  51),I=1,7) /-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  52),I=1,7) /-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  53),I=1,7) /-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  54),I=1,7) /-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  55),I=1,7) /-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  56),I=1,7) /-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  57),I=1,7) /-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  58),I=1,7) /-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  59),I=1,7) /-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  60),I=1,7) /-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  61),I=1,7) /-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  62),I=1,7) /-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  63),I=1,7) /-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  64),I=1,7) /-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  65),I=1,7) / 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  66),I=1,7) / 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  67),I=1,7) / 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  68),I=1,7) / 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  69),I=1,7) / 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  70),I=1,7) / 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  71),I=1,7) / 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  72),I=1,7) / 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  73),I=1,7) / 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  74),I=1,7) / 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  75),I=1,7) / 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  76),I=1,7) / 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  77),I=1,7) / 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  78),I=1,7) / 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  79),I=1,7) / 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  80),I=1,7) / 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  81),I=1,7) / 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  82),I=1,7) / 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  83),I=1,7) / 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  84),I=1,7) / 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  85),I=1,7) / 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  86),I=1,7) / 1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  87),I=1,7) / 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  88),I=1,7) / 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  89),I=1,7) / 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  90),I=1,7) / 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  91),I=1,7) / 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  92),I=1,7) / 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  93),I=1,7) / 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  94),I=1,7) / 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  95),I=1,7) / 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  96),I=1,7) / 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  97),I=1,7) / 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  98),I=1,7) / 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  99),I=1,7) / 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I, 100),I=1,7) / 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I, 101),I=1,7) / 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I, 102),I=1,7) / 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I, 103),I=1,7) / 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I, 104),I=1,7) / 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(I, 105),I=1,7) / 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I, 106),I=1,7) / 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I, 107),I=1,7) / 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I, 108),I=1,7) / 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I, 109),I=1,7) / 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I, 110),I=1,7) / 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I, 111),I=1,7) / 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I, 112),I=1,7) / 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I, 113),I=1,7) / 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I, 114),I=1,7) / 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I, 115),I=1,7) / 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I, 116),I=1,7) / 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I, 117),I=1,7) / 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I, 118),I=1,7) / 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I, 119),I=1,7) / 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I, 120),I=1,7) / 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I, 121),I=1,7) / 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I, 122),I=1,7) / 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I, 123),I=1,7) / 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I, 124),I=1,7) / 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I, 125),I=1,7) / 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I, 126),I=1,7) / 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I, 127),I=1,7) / 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I, 128),I=1,7) / 1,-1, 1,-1, 1,-1, 1/
      DOUBLE COMPLEX ANS_SUMMED

      INCLUDE 'orders.inc'
      DOUBLE COMPLEX AMP_SPLIT_EWSUD(AMP_SPLIT_SIZE)
      COMMON /TO_AMP_SPLIT_EWSUD/ AMP_SPLIT_EWSUD
      DOUBLE COMPLEX AMP_SPLIT_EWSUD_LO2(AMP_SPLIT_SIZE)
      COMMON /TO_AMP_SPLIT_EWSUD_LO2/ AMP_SPLIT_EWSUD_LO2

      DOUBLE COMPLEX AMP_SPLIT_EWSUD_LSC(AMP_SPLIT_SIZE)
      COMMON /TO_AMP_EWSUD_LSC/AMP_SPLIT_EWSUD_LSC
      DOUBLE COMPLEX AMP_SPLIT_EWSUD_SSC(AMP_SPLIT_SIZE)
      COMMON /TO_AMP_EWSUD_SSC/AMP_SPLIT_EWSUD_SSC
      DOUBLE COMPLEX AMP_SPLIT_EWSUD_XXC(AMP_SPLIT_SIZE)
      COMMON /TO_AMP_EWSUD_XXC/AMP_SPLIT_EWSUD_XXC
      DOUBLE COMPLEX AMP_SPLIT_EWSUD_PAR(AMP_SPLIT_SIZE)
      COMMON /TO_AMP_EWSUD_PAR/AMP_SPLIT_EWSUD_PAR
      DOUBLE COMPLEX AMP_SPLIT_EWSUD_QCD(AMP_SPLIT_SIZE)
      COMMON /TO_AMP_EWSUD_QCD/AMP_SPLIT_EWSUD_QCD
      DOUBLE COMPLEX AMP_SPLIT_BORN_ONEHEL(AMP_SPLIT_SIZE)
      COMMON /TO_AMP_BORN_ONEHEL/AMP_SPLIT_BORN_ONEHEL
      DOUBLE COMPLEX AMP_SPLIT_EWSUD_PARQCD(AMP_SPLIT_SIZE)
      COMMON /TO_AMP_EWSUD_PARQCD/AMP_SPLIT_EWSUD_PARQCD

      INTEGER IFLIST(NEXTERNAL-1), PDGLIST(NEXTERNAL-1),
     $  PDGLIST_OTH(NEXTERNAL-1)

      DOUBLE PRECISION INVARIANTS(NEXTERNAL-1, NEXTERNAL-1)

      DOUBLE COMPLEX GET_LSC_DIAG, GET_LSC_NONDIAG
      EXTERNAL GET_LSC_DIAG, GET_LSC_NONDIAG
      DOUBLE COMPLEX GET_SSC_N_DIAG, GET_SSC_N_NONDIAG_1,
     $  GET_SSC_N_NONDIAG_2, GET_SSC_C
      EXTERNAL GET_SSC_N_DIAG, GET_SSC_N_NONDIAG_1,
     $  GET_SSC_N_NONDIAG_2, GET_SSC_C
      DOUBLE COMPLEX GET_XXC_DIAG, GET_XXC_NONDIAG
      EXTERNAL GET_XXC_DIAG, GET_XXC_NONDIAG
      DOUBLE COMPLEX GET_QCD_LO2
      EXTERNAL GET_QCD_LO2

      LOGICAL CALCULATEDBORN
      COMMON/CCALCULATEDBORN/CALCULATEDBORN

      INTEGER EWSUD_HELSELECT
      COMMON/TO_EWSUD_HELSELECT/EWSUD_HELSELECT

      DOUBLE PRECISION COMPENSATE_IDENTICAL_FACTOR
      EXTERNAL COMPENSATE_IDENTICAL_FACTOR
      DOUBLE PRECISION COMP_IDFAC

      LOGICAL GOOD_HEL(NCOMB), FIRSTTIME_HEL(NCOMB)
      DATA GOOD_HEL  /NCOMB*.TRUE./
      DATA FIRSTTIME_HEL  /NCOMB*.TRUE./
      INTEGER IAMP
      LOGICAL SUD_FILTER_HEL
      COMMON /TO_FILTER_HEL/ SUD_FILTER_HEL
C     stuff for the MC over helicities
      LOGICAL SUD_MC_HEL
      COMMON /TO_MC_HEL/ SUD_MC_HEL
      INTEGER IGOOD_HEL(NCOMB), HEL(0:NCOMB)
      DATA IGOOD_HEL  /NCOMB * 0/
      LOGICAL FIRSTTIME
      DATA FIRSTTIME /.TRUE./
      INTEGER HEL_MIN, HEL_MAX, HEL_PICK
      DOUBLE PRECISION VOLH

      INCLUDE 'ewsudakov_haslo.inc'


      AMP_SPLIT_EWSUD_LSC(:) = DCMPLX(0D0,0D0)
      AMP_SPLIT_EWSUD_SSC(:) = DCMPLX(0D0,0D0)
      AMP_SPLIT_EWSUD_XXC(:) = DCMPLX(0D0,0D0)
      AMP_SPLIT_EWSUD_PAR(:) = DCMPLX(0D0,0D0)

      AMP_SPLIT_EWSUD_QCD(:) = DCMPLX(0D0,0D0)
      AMP_SPLIT_EWSUD_PARQCD(:) = DCMPLX(0D0,0D0)


      IFLIST = (/-1,-1,1,1,1,1,1/)

      CALL SDK_GET_INVARIANTS(P, IFLIST, INVARIANTS)

      IF (FIRSTTIME.OR..NOT.SUD_MC_HEL) THEN
        HEL_MIN = 1
        HEL_MAX = NCOMB
        HEL(:)=0
        VOLH=1D0
        FIRSTTIME=.FALSE.
      ELSE
        CALL PICKHELICITYMC(P,IGOOD_HEL,HEL,HEL_PICK,VOLH)
        HEL_MIN = HEL(HEL_PICK)
        HEL_MAX = HEL_MIN
      ENDIF

      CALCULATEDBORN=.FALSE.

C     Main loop over the helicities
      DO IHEL=HEL_MIN, HEL_MAX
C       if an helicity is selected, then filter it
        IF (EWSUD_HELSELECT.NE.0.AND.EWSUD_HELSELECT.NE.IHEL) CYCLE

        IF (.NOT.GOOD_HEL(IHEL).AND.SUD_FILTER_HEL) CYCLE

C       call the born amplitude and store the result
        CALL SBORN_ONEHEL(P,NHEL(1,IHEL),IHEL,ANS_SUMMED)
        AMP_SPLIT_BORN_ONEHEL(:) = AMP_SPLIT_EWSUD(:)

        IF (FIRSTTIME_HEL(IHEL)) THEN
          GOOD_HEL(IHEL) = .FALSE.
          DO IAMP = 1, AMP_SPLIT_SIZE
            GOOD_HEL(IHEL) = GOOD_HEL(IHEL)
     $       .OR.AMP_SPLIT_BORN_ONEHEL(IAMP)
     $       .NE.0D0.OR.AMP_SPLIT_EWSUD_LO2(IAMP).NE.0D0
          ENDDO

          IF (GOOD_HEL(IHEL)) THEN
            HEL(0) = HEL(0)+1
            HEL(HEL(0)) = IHEL
            IGOOD_HEL(HEL(0)) = 1
          ENDIF
C         write(*,*) 'EW SUDAKOV WRAPPER HEL FILTER', ihel,
C          nhel(:,ihel), good_hel(ihel)
          FIRSTTIME_HEL(IHEL) = .FALSE.
        ENDIF

C       call the amplitudes with the goldstones for the
C        longitudinally-polarised legs and the born for all-transverse

      ENDDO

      AMP_SPLIT_EWSUD_LSC(:) = AMP_SPLIT_EWSUD_LSC(:)/VOLH
      AMP_SPLIT_EWSUD_SSC(:) = AMP_SPLIT_EWSUD_SSC(:)/VOLH
      AMP_SPLIT_EWSUD_XXC(:) = AMP_SPLIT_EWSUD_XXC(:)/VOLH
      AMP_SPLIT_EWSUD_PAR(:) = AMP_SPLIT_EWSUD_PAR(:)/VOLH

      AMP_SPLIT_EWSUD_QCD(:) = AMP_SPLIT_EWSUD_QCD(:)/VOLH
      AMP_SPLIT_EWSUD_PARQCD(:) = AMP_SPLIT_EWSUD_PAR(:)

      DO IAMP = 1, AMP_SPLIT_SIZE

        IF (HAS_LO1) THEN
          IF(IAMP.NE.LO1_POS) AMP_SPLIT_EWSUD_PAR(IAMP)=(0D0,0D0)
        ELSE
          AMP_SPLIT_EWSUD_PAR(IAMP)=(0D0,0D0)
        ENDIF

        IF (HAS_LO2) THEN
          IF(IAMP.NE.LO2_POS) AMP_SPLIT_EWSUD_PARQCD(IAMP)=(0D0,0D0)
        ELSE
          AMP_SPLIT_EWSUD_PARQCD(IAMP)=(0D0,0D0)
        ENDIF

      ENDDO



      CALCULATEDBORN=.FALSE.

      RETURN
      END





      INTEGER FUNCTION SDK_GET_NCOMB()
      IMPLICIT NONE
      INTEGER     NCOMB
      PARAMETER ( NCOMB=  128 )

      SDK_GET_NCOMB = NCOMB

      RETURN
      END


      SUBROUTINE SDK_GET_HELS(IHEL, HELS)
      IMPLICIT NONE
      INCLUDE 'nexternal.inc'
      INTEGER IHEL, I
      INTEGER HELS(NEXTERNAL-1)
      INTEGER     NCOMB
      PARAMETER ( NCOMB=  128 )
      INTEGER NHEL(NEXTERNAL-1,NCOMB)
      DATA (NHEL(I,   1),I=1,7) /-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,   2),I=1,7) /-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,   3),I=1,7) /-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,   4),I=1,7) /-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,   5),I=1,7) /-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,   6),I=1,7) /-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,   7),I=1,7) /-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,   8),I=1,7) /-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,   9),I=1,7) /-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  10),I=1,7) /-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  11),I=1,7) /-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  12),I=1,7) /-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  13),I=1,7) /-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  14),I=1,7) /-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  15),I=1,7) /-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  16),I=1,7) /-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  17),I=1,7) /-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  18),I=1,7) /-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  19),I=1,7) /-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  20),I=1,7) /-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  21),I=1,7) /-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  22),I=1,7) /-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  23),I=1,7) /-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  24),I=1,7) /-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  25),I=1,7) /-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  26),I=1,7) /-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  27),I=1,7) /-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  28),I=1,7) /-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  29),I=1,7) /-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  30),I=1,7) /-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  31),I=1,7) /-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  32),I=1,7) /-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  33),I=1,7) /-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  34),I=1,7) /-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  35),I=1,7) /-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  36),I=1,7) /-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  37),I=1,7) /-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  38),I=1,7) /-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  39),I=1,7) /-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  40),I=1,7) /-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  41),I=1,7) /-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  42),I=1,7) /-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  43),I=1,7) /-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  44),I=1,7) /-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  45),I=1,7) /-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  46),I=1,7) /-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  47),I=1,7) /-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  48),I=1,7) /-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  49),I=1,7) /-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  50),I=1,7) /-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  51),I=1,7) /-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  52),I=1,7) /-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  53),I=1,7) /-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  54),I=1,7) /-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  55),I=1,7) /-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  56),I=1,7) /-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  57),I=1,7) /-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  58),I=1,7) /-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  59),I=1,7) /-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  60),I=1,7) /-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  61),I=1,7) /-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  62),I=1,7) /-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  63),I=1,7) /-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  64),I=1,7) /-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  65),I=1,7) / 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  66),I=1,7) / 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  67),I=1,7) / 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  68),I=1,7) / 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  69),I=1,7) / 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  70),I=1,7) / 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  71),I=1,7) / 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  72),I=1,7) / 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  73),I=1,7) / 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  74),I=1,7) / 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  75),I=1,7) / 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  76),I=1,7) / 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  77),I=1,7) / 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  78),I=1,7) / 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  79),I=1,7) / 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  80),I=1,7) / 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  81),I=1,7) / 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  82),I=1,7) / 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  83),I=1,7) / 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  84),I=1,7) / 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  85),I=1,7) / 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  86),I=1,7) / 1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  87),I=1,7) / 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  88),I=1,7) / 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  89),I=1,7) / 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  90),I=1,7) / 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  91),I=1,7) / 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  92),I=1,7) / 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  93),I=1,7) / 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  94),I=1,7) / 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  95),I=1,7) / 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  96),I=1,7) / 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  97),I=1,7) / 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  98),I=1,7) / 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  99),I=1,7) / 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I, 100),I=1,7) / 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I, 101),I=1,7) / 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I, 102),I=1,7) / 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I, 103),I=1,7) / 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I, 104),I=1,7) / 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(I, 105),I=1,7) / 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I, 106),I=1,7) / 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I, 107),I=1,7) / 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I, 108),I=1,7) / 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I, 109),I=1,7) / 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I, 110),I=1,7) / 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I, 111),I=1,7) / 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I, 112),I=1,7) / 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I, 113),I=1,7) / 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I, 114),I=1,7) / 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I, 115),I=1,7) / 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I, 116),I=1,7) / 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I, 117),I=1,7) / 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I, 118),I=1,7) / 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I, 119),I=1,7) / 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I, 120),I=1,7) / 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I, 121),I=1,7) / 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I, 122),I=1,7) / 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I, 123),I=1,7) / 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I, 124),I=1,7) / 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I, 125),I=1,7) / 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I, 126),I=1,7) / 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I, 127),I=1,7) / 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I, 128),I=1,7) / 1,-1, 1,-1, 1,-1, 1/

      HELS(:) = NHEL(:, IHEL)

      RETURN
      END


      DOUBLE PRECISION FUNCTION COMPENSATE_IDENTICAL_FACTOR(IGOLD_ME)
C     provides the multiplicative factor for the result of the
C     sudakov MEs so that they have the same identical-particle
C     factor as the Born in the specific FKS configuration
      IMPLICIT NONE
      INTEGER IGOLD_ME
C     
      INTEGER SDK_IDENT_GOLDSTONE(1)
      DATA SDK_IDENT_GOLDSTONE / 0 /
C     
      INTEGER NFKSPROCESS
      COMMON/C_NFKSPROCESS/NFKSPROCESS
      INTEGER IDEN_VALUES(5)
      DATA IDEN_VALUES /4, 4, 4, 4, 4/
      INTEGER BORNSPINCOL
      DATA BORNSPINCOL / 4 /

      COMPENSATE_IDENTICAL_FACTOR = DBLE(SDK_IDENT_GOLDSTONE(IGOLD_ME))
     $ /DBLE(IDEN_VALUES(NFKSPROCESS)/BORNSPINCOL)

      RETURN
      END


