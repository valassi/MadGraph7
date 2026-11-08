      SUBROUTINE HELAS_CALLS_AMPB_1(P,NHEL,H,IC)
C     
C     Modules
C     
      USE POLYNOMIAL_CONSTANTS
      USE ALOHA_OBJECT
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NCOMB
      PARAMETER (NCOMB=64)
      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=16)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=264, NLOOPGROUPS=114, NCTAMPS=540)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=804)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=44,NLOOPWAVEFUNCS=413)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      REAL*16     MP__ZERO
      PARAMETER (MP__ZERO=0.0E0_16)
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=2, NSQUAREDSO=1, NAMPSO=2)
C     
C     ARGUMENTS
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      INTEGER H
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
      INTEGER FLAVOR(NEXTERNAL)
      DATA FLAVOR /NEXTERNAL*1/
      COMPLEX*16 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)

      LOGICAL DUMMYFALSE
      DATA DUMMYFALSE/.FALSE./
C     
C     GLOBAL VARIABLES
C     

      INCLUDE 'coupl.inc'
      INCLUDE 'mp_coupl.inc'

      INTEGER HELOFFSET
      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/FILTERS/GOODAMP,GOODHEL,HELOFFSET

      LOGICAL CHECKPHASE
      LOGICAL HELDOUBLECHECKED
      COMMON/INIT/CHECKPHASE, HELDOUBLECHECKED

      INTEGER SQSO_TARGET
      COMMON/SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      INTEGER I_SO
      COMMON/I_SO/I_SO
      INTEGER I_LIB
      COMMON/I_LIB/I_LIB

      COMPLEX*16 AMP(NBORNAMPS)
      COMMON/AMPS/AMP
      TYPE(ALOHA) W(NWAVEFUNCS)
      COMMON/W/W
      COMPLEX*16 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,
     $ -1:NLOOPWAVEFUNCS)
      COMPLEX*16 PL(0:3,-1:NLOOPWAVEFUNCS)
      COMMON/WL/WL,PL

      COMPLEX*16 AMPL(3,NCTAMPS)
      COMMON/AMPL/AMPL

C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.CT_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

      CALL OXXXXX(P(0,1),ZERO,NHEL(1),-1, FLAVOR(1),W(1))
      CALL IXXXXX(P(0,2),ZERO,NHEL(2),+1, FLAVOR(2),W(2))
      CALL OXXXXX(P(0,3),MDL_MB,NHEL(3),+1, FLAVOR(3),W(3))
      CALL OXXXXX(P(0,4),MDL_MB,NHEL(4),+1, FLAVOR(4),W(4))
      CALL IXXXXX(P(0,5),MDL_MB,NHEL(5),-1, FLAVOR(5),W(5))
      CALL IXXXXX(P(0,6),MDL_MB,NHEL(6),-1, FLAVOR(6),W(6))
      CALL FFV1P0_3(W(2),W(1),GC_3,ZERO,ZERO,W(7))
      CALL FFV1P0_3(W(5),W(3),GC_5,ZERO,ZERO,W(8))
      CALL FFV1_1(W(4),W(7),GC_1,MDL_MB,ZERO,W(9))
C     Amplitude(s) for born diagram with ID 1
      CALL FFV1_0(W(6),W(9),W(8),GC_5,AMP(1))
      CALL FFV1_2(W(6),W(7),GC_1,MDL_MB,ZERO,W(10))
C     Amplitude(s) for born diagram with ID 2
      CALL FFV1_0(W(10),W(4),W(8),GC_5,AMP(2))
      CALL FFV2_4_3(W(2),W(1),GC_21,GC_24,MDL_MZ,MDL_WZ,W(11))
      CALL FFV2_3_1(W(4),W(11),GC_21,GC_23,MDL_MB,ZERO,W(12))
C     Amplitude(s) for born diagram with ID 3
      CALL FFV1_0(W(6),W(12),W(8),GC_5,AMP(3))
      CALL FFV2_3_2(W(6),W(11),GC_21,GC_23,MDL_MB,ZERO,W(13))
C     Amplitude(s) for born diagram with ID 4
      CALL FFV1_0(W(13),W(4),W(8),GC_5,AMP(4))
      CALL FFV1P0_3(W(6),W(3),GC_5,ZERO,ZERO,W(14))
C     Amplitude(s) for born diagram with ID 5
      CALL FFV1_0(W(5),W(9),W(14),GC_5,AMP(5))
      CALL FFV1_2(W(5),W(7),GC_1,MDL_MB,ZERO,W(15))
C     Amplitude(s) for born diagram with ID 6
      CALL FFV1_0(W(15),W(4),W(14),GC_5,AMP(6))
C     Amplitude(s) for born diagram with ID 7
      CALL FFV1_0(W(5),W(12),W(14),GC_5,AMP(7))
      CALL FFV2_3_2(W(5),W(11),GC_21,GC_23,MDL_MB,ZERO,W(16))
C     Amplitude(s) for born diagram with ID 8
      CALL FFV1_0(W(16),W(4),W(14),GC_5,AMP(8))
      CALL FFV1P0_3(W(5),W(4),GC_5,ZERO,ZERO,W(17))
      CALL FFV1_1(W(3),W(7),GC_1,MDL_MB,ZERO,W(18))
C     Amplitude(s) for born diagram with ID 9
      CALL FFV1_0(W(6),W(18),W(17),GC_5,AMP(9))
C     Amplitude(s) for born diagram with ID 10
      CALL FFV1_0(W(10),W(3),W(17),GC_5,AMP(10))
      CALL FFV2_3_1(W(3),W(11),GC_21,GC_23,MDL_MB,ZERO,W(19))
C     Amplitude(s) for born diagram with ID 11
      CALL FFV1_0(W(6),W(19),W(17),GC_5,AMP(11))
C     Amplitude(s) for born diagram with ID 12
      CALL FFV1_0(W(13),W(3),W(17),GC_5,AMP(12))
      CALL FFV1P0_3(W(6),W(4),GC_5,ZERO,ZERO,W(20))
C     Amplitude(s) for born diagram with ID 13
      CALL FFV1_0(W(5),W(18),W(20),GC_5,AMP(13))
C     Amplitude(s) for born diagram with ID 14
      CALL FFV1_0(W(15),W(3),W(20),GC_5,AMP(14))
C     Amplitude(s) for born diagram with ID 15
      CALL FFV1_0(W(5),W(19),W(20),GC_5,AMP(15))
C     Amplitude(s) for born diagram with ID 16
      CALL FFV1_0(W(16),W(3),W(20),GC_5,AMP(16))
      CALL FFV1_2(W(6),W(8),GC_5,MDL_MB,ZERO,W(21))
C     Counter-term amplitude(s) for loop diagram number 17
      CALL R2_QQ_1_R2_QQ_2_0(W(21),W(9),R2_QQQ,R2_QQB,AMPL(1,1))
      CALL R2_QQ_2_0(W(21),W(9),UV_BMASS_1EPS,AMPL(2,2))
      CALL R2_QQ_2_0(W(21),W(9),UV_BMASS,AMPL(1,3))
      CALL FFV1P0_3(W(6),W(9),GC_5,ZERO,ZERO,W(22))
C     Counter-term amplitude(s) for loop diagram number 18
      CALL R2_GG_1_R2_GG_2_0(W(8),W(22),R2_GGG_1,R2_GGG_2,AMPL(1,4))
C     Counter-term amplitude(s) for loop diagram number 19
      CALL FFV1_0(W(6),W(9),W(8),R2_GQQ,AMPL(1,5))
      CALL FFV1_0(W(6),W(9),W(8),UV_GQQB_1EPS,AMPL(2,6))
      CALL FFV1_0(W(6),W(9),W(8),UV_GQQB_1EPS,AMPL(2,7))
      CALL FFV1_0(W(6),W(9),W(8),UV_GQQB_1EPS,AMPL(2,8))
      CALL FFV1_0(W(6),W(9),W(8),UV_GQQB_1EPS,AMPL(2,9))
      CALL FFV1_0(W(6),W(9),W(8),UV_GQQB_1EPS,AMPL(2,10))
      CALL FFV1_0(W(6),W(9),W(8),UV_GQQB_1EPS,AMPL(2,11))
      CALL FFV1_0(W(6),W(9),W(8),UV_GQQG_1EPS,AMPL(2,12))
      CALL FFV1_0(W(6),W(9),W(8),UV_GQQB,AMPL(1,13))
      CALL FFV1_0(W(6),W(9),W(8),UV_GQQT,AMPL(1,14))
      CALL FFV1_1(W(4),W(8),GC_5,MDL_MB,ZERO,W(23))
C     Counter-term amplitude(s) for loop diagram number 20
      CALL R2_QQ_1_R2_QQ_2_0(W(10),W(23),R2_QQQ,R2_QQB,AMPL(1,15))
      CALL R2_QQ_2_0(W(10),W(23),UV_BMASS_1EPS,AMPL(2,16))
      CALL R2_QQ_2_0(W(10),W(23),UV_BMASS,AMPL(1,17))
      CALL FFV1P0_3(W(10),W(4),GC_5,ZERO,ZERO,W(24))
C     Counter-term amplitude(s) for loop diagram number 21
      CALL R2_GG_1_R2_GG_2_0(W(8),W(24),R2_GGG_1,R2_GGG_2,AMPL(1,18))
C     Counter-term amplitude(s) for loop diagram number 22
      CALL FFV1_0(W(10),W(4),W(8),R2_GQQ,AMPL(1,19))
      CALL FFV1_0(W(10),W(4),W(8),UV_GQQB_1EPS,AMPL(2,20))
      CALL FFV1_0(W(10),W(4),W(8),UV_GQQB_1EPS,AMPL(2,21))
      CALL FFV1_0(W(10),W(4),W(8),UV_GQQB_1EPS,AMPL(2,22))
      CALL FFV1_0(W(10),W(4),W(8),UV_GQQB_1EPS,AMPL(2,23))
      CALL FFV1_0(W(10),W(4),W(8),UV_GQQB_1EPS,AMPL(2,24))
      CALL FFV1_0(W(10),W(4),W(8),UV_GQQB_1EPS,AMPL(2,25))
      CALL FFV1_0(W(10),W(4),W(8),UV_GQQG_1EPS,AMPL(2,26))
      CALL FFV1_0(W(10),W(4),W(8),UV_GQQB,AMPL(1,27))
      CALL FFV1_0(W(10),W(4),W(8),UV_GQQT,AMPL(1,28))
C     Counter-term amplitude(s) for loop diagram number 23
      CALL R2_QQ_1_R2_QQ_2_0(W(21),W(12),R2_QQQ,R2_QQB,AMPL(1,29))
      CALL R2_QQ_2_0(W(21),W(12),UV_BMASS_1EPS,AMPL(2,30))
      CALL R2_QQ_2_0(W(21),W(12),UV_BMASS,AMPL(1,31))
      CALL FFV1P0_3(W(6),W(12),GC_5,ZERO,ZERO,W(25))
C     Counter-term amplitude(s) for loop diagram number 24
      CALL R2_GG_1_R2_GG_2_0(W(8),W(25),R2_GGG_1,R2_GGG_2,AMPL(1,32))
C     Counter-term amplitude(s) for loop diagram number 25
      CALL FFV1_0(W(6),W(12),W(8),R2_GQQ,AMPL(1,33))
      CALL FFV1_0(W(6),W(12),W(8),UV_GQQB_1EPS,AMPL(2,34))
      CALL FFV1_0(W(6),W(12),W(8),UV_GQQB_1EPS,AMPL(2,35))
      CALL FFV1_0(W(6),W(12),W(8),UV_GQQB_1EPS,AMPL(2,36))
      CALL FFV1_0(W(6),W(12),W(8),UV_GQQB_1EPS,AMPL(2,37))
      CALL FFV1_0(W(6),W(12),W(8),UV_GQQB_1EPS,AMPL(2,38))
      CALL FFV1_0(W(6),W(12),W(8),UV_GQQB_1EPS,AMPL(2,39))
      CALL FFV1_0(W(6),W(12),W(8),UV_GQQG_1EPS,AMPL(2,40))
      CALL FFV1_0(W(6),W(12),W(8),UV_GQQB,AMPL(1,41))
      CALL FFV1_0(W(6),W(12),W(8),UV_GQQT,AMPL(1,42))
C     Counter-term amplitude(s) for loop diagram number 26
      CALL R2_QQ_1_R2_QQ_2_0(W(13),W(23),R2_QQQ,R2_QQB,AMPL(1,43))
      CALL R2_QQ_2_0(W(13),W(23),UV_BMASS_1EPS,AMPL(2,44))
      CALL R2_QQ_2_0(W(13),W(23),UV_BMASS,AMPL(1,45))
      CALL FFV1P0_3(W(13),W(4),GC_5,ZERO,ZERO,W(26))
C     Counter-term amplitude(s) for loop diagram number 27
      CALL R2_GG_1_R2_GG_2_0(W(8),W(26),R2_GGG_1,R2_GGG_2,AMPL(1,46))
C     Counter-term amplitude(s) for loop diagram number 28
      CALL FFV1_0(W(13),W(4),W(8),R2_GQQ,AMPL(1,47))
      CALL FFV1_0(W(13),W(4),W(8),UV_GQQB_1EPS,AMPL(2,48))
      CALL FFV1_0(W(13),W(4),W(8),UV_GQQB_1EPS,AMPL(2,49))
      CALL FFV1_0(W(13),W(4),W(8),UV_GQQB_1EPS,AMPL(2,50))
      CALL FFV1_0(W(13),W(4),W(8),UV_GQQB_1EPS,AMPL(2,51))
      CALL FFV1_0(W(13),W(4),W(8),UV_GQQB_1EPS,AMPL(2,52))
      CALL FFV1_0(W(13),W(4),W(8),UV_GQQB_1EPS,AMPL(2,53))
      CALL FFV1_0(W(13),W(4),W(8),UV_GQQG_1EPS,AMPL(2,54))
      CALL FFV1_0(W(13),W(4),W(8),UV_GQQB,AMPL(1,55))
      CALL FFV1_0(W(13),W(4),W(8),UV_GQQT,AMPL(1,56))
C     Counter-term amplitude(s) for loop diagram number 29
      CALL FFV1_0(W(21),W(4),W(7),R2_DDA,AMPL(1,57))
C     Counter-term amplitude(s) for loop diagram number 32
      CALL FFV2_3_0(W(21),W(4),W(11),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,58))
C     Counter-term amplitude(s) for loop diagram number 40
      CALL FFV1_0(W(6),W(23),W(7),R2_DDA,AMPL(1,59))
C     Counter-term amplitude(s) for loop diagram number 42
      CALL FFV2_3_0(W(6),W(23),W(11),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,60))
      CALL FFV1_2(W(5),W(14),GC_5,MDL_MB,ZERO,W(27))
C     Counter-term amplitude(s) for loop diagram number 43
      CALL R2_QQ_1_R2_QQ_2_0(W(27),W(9),R2_QQQ,R2_QQB,AMPL(1,61))
      CALL R2_QQ_2_0(W(27),W(9),UV_BMASS_1EPS,AMPL(2,62))
      CALL R2_QQ_2_0(W(27),W(9),UV_BMASS,AMPL(1,63))
      CALL FFV1P0_3(W(5),W(9),GC_5,ZERO,ZERO,W(28))
C     Counter-term amplitude(s) for loop diagram number 44
      CALL R2_GG_1_R2_GG_2_0(W(14),W(28),R2_GGG_1,R2_GGG_2,AMPL(1,64))
C     Counter-term amplitude(s) for loop diagram number 45
      CALL FFV1_0(W(5),W(9),W(14),R2_GQQ,AMPL(1,65))
      CALL FFV1_0(W(5),W(9),W(14),UV_GQQB_1EPS,AMPL(2,66))
      CALL FFV1_0(W(5),W(9),W(14),UV_GQQB_1EPS,AMPL(2,67))
      CALL FFV1_0(W(5),W(9),W(14),UV_GQQB_1EPS,AMPL(2,68))
      CALL FFV1_0(W(5),W(9),W(14),UV_GQQB_1EPS,AMPL(2,69))
      CALL FFV1_0(W(5),W(9),W(14),UV_GQQB_1EPS,AMPL(2,70))
      CALL FFV1_0(W(5),W(9),W(14),UV_GQQB_1EPS,AMPL(2,71))
      CALL FFV1_0(W(5),W(9),W(14),UV_GQQG_1EPS,AMPL(2,72))
      CALL FFV1_0(W(5),W(9),W(14),UV_GQQB,AMPL(1,73))
      CALL FFV1_0(W(5),W(9),W(14),UV_GQQT,AMPL(1,74))
      CALL FFV1_1(W(4),W(14),GC_5,MDL_MB,ZERO,W(29))
C     Counter-term amplitude(s) for loop diagram number 46
      CALL R2_QQ_1_R2_QQ_2_0(W(15),W(29),R2_QQQ,R2_QQB,AMPL(1,75))
      CALL R2_QQ_2_0(W(15),W(29),UV_BMASS_1EPS,AMPL(2,76))
      CALL R2_QQ_2_0(W(15),W(29),UV_BMASS,AMPL(1,77))
      CALL FFV1P0_3(W(15),W(4),GC_5,ZERO,ZERO,W(30))
C     Counter-term amplitude(s) for loop diagram number 47
      CALL R2_GG_1_R2_GG_2_0(W(14),W(30),R2_GGG_1,R2_GGG_2,AMPL(1,78))
C     Counter-term amplitude(s) for loop diagram number 48
      CALL FFV1_0(W(15),W(4),W(14),R2_GQQ,AMPL(1,79))
      CALL FFV1_0(W(15),W(4),W(14),UV_GQQB_1EPS,AMPL(2,80))
      CALL FFV1_0(W(15),W(4),W(14),UV_GQQB_1EPS,AMPL(2,81))
      CALL FFV1_0(W(15),W(4),W(14),UV_GQQB_1EPS,AMPL(2,82))
      CALL FFV1_0(W(15),W(4),W(14),UV_GQQB_1EPS,AMPL(2,83))
      CALL FFV1_0(W(15),W(4),W(14),UV_GQQB_1EPS,AMPL(2,84))
      CALL FFV1_0(W(15),W(4),W(14),UV_GQQB_1EPS,AMPL(2,85))
      CALL FFV1_0(W(15),W(4),W(14),UV_GQQG_1EPS,AMPL(2,86))
      CALL FFV1_0(W(15),W(4),W(14),UV_GQQB,AMPL(1,87))
      CALL FFV1_0(W(15),W(4),W(14),UV_GQQT,AMPL(1,88))
C     Counter-term amplitude(s) for loop diagram number 49
      CALL R2_QQ_1_R2_QQ_2_0(W(27),W(12),R2_QQQ,R2_QQB,AMPL(1,89))
      CALL R2_QQ_2_0(W(27),W(12),UV_BMASS_1EPS,AMPL(2,90))
      CALL R2_QQ_2_0(W(27),W(12),UV_BMASS,AMPL(1,91))
      CALL FFV1P0_3(W(5),W(12),GC_5,ZERO,ZERO,W(31))
C     Counter-term amplitude(s) for loop diagram number 50
      CALL R2_GG_1_R2_GG_2_0(W(14),W(31),R2_GGG_1,R2_GGG_2,AMPL(1,92))
C     Counter-term amplitude(s) for loop diagram number 51
      CALL FFV1_0(W(5),W(12),W(14),R2_GQQ,AMPL(1,93))
      CALL FFV1_0(W(5),W(12),W(14),UV_GQQB_1EPS,AMPL(2,94))
      CALL FFV1_0(W(5),W(12),W(14),UV_GQQB_1EPS,AMPL(2,95))
      CALL FFV1_0(W(5),W(12),W(14),UV_GQQB_1EPS,AMPL(2,96))
      CALL FFV1_0(W(5),W(12),W(14),UV_GQQB_1EPS,AMPL(2,97))
      CALL FFV1_0(W(5),W(12),W(14),UV_GQQB_1EPS,AMPL(2,98))
      CALL FFV1_0(W(5),W(12),W(14),UV_GQQB_1EPS,AMPL(2,99))
      CALL FFV1_0(W(5),W(12),W(14),UV_GQQG_1EPS,AMPL(2,100))
      CALL FFV1_0(W(5),W(12),W(14),UV_GQQB,AMPL(1,101))
      CALL FFV1_0(W(5),W(12),W(14),UV_GQQT,AMPL(1,102))
C     Counter-term amplitude(s) for loop diagram number 52
      CALL R2_QQ_1_R2_QQ_2_0(W(16),W(29),R2_QQQ,R2_QQB,AMPL(1,103))
      CALL R2_QQ_2_0(W(16),W(29),UV_BMASS_1EPS,AMPL(2,104))
      CALL R2_QQ_2_0(W(16),W(29),UV_BMASS,AMPL(1,105))
      CALL FFV1P0_3(W(16),W(4),GC_5,ZERO,ZERO,W(32))
C     Counter-term amplitude(s) for loop diagram number 53
      CALL R2_GG_1_R2_GG_2_0(W(14),W(32),R2_GGG_1,R2_GGG_2,AMPL(1,106))
C     Counter-term amplitude(s) for loop diagram number 54
      CALL FFV1_0(W(16),W(4),W(14),R2_GQQ,AMPL(1,107))
      CALL FFV1_0(W(16),W(4),W(14),UV_GQQB_1EPS,AMPL(2,108))
      CALL FFV1_0(W(16),W(4),W(14),UV_GQQB_1EPS,AMPL(2,109))
      CALL FFV1_0(W(16),W(4),W(14),UV_GQQB_1EPS,AMPL(2,110))
      CALL FFV1_0(W(16),W(4),W(14),UV_GQQB_1EPS,AMPL(2,111))
      CALL FFV1_0(W(16),W(4),W(14),UV_GQQB_1EPS,AMPL(2,112))
      CALL FFV1_0(W(16),W(4),W(14),UV_GQQB_1EPS,AMPL(2,113))
      CALL FFV1_0(W(16),W(4),W(14),UV_GQQG_1EPS,AMPL(2,114))
      CALL FFV1_0(W(16),W(4),W(14),UV_GQQB,AMPL(1,115))
      CALL FFV1_0(W(16),W(4),W(14),UV_GQQT,AMPL(1,116))
C     Counter-term amplitude(s) for loop diagram number 55
      CALL FFV1_0(W(27),W(4),W(7),R2_DDA,AMPL(1,117))
C     Counter-term amplitude(s) for loop diagram number 58
      CALL FFV2_3_0(W(27),W(4),W(11),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,118))
C     Counter-term amplitude(s) for loop diagram number 66
      CALL FFV1_0(W(5),W(29),W(7),R2_DDA,AMPL(1,119))
C     Counter-term amplitude(s) for loop diagram number 68
      CALL FFV2_3_0(W(5),W(29),W(11),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,120))
C     Counter-term amplitude(s) for loop diagram number 69
      CALL FFV1_0(W(5),W(3),W(22),R2_GQQ,AMPL(1,121))
      CALL FFV1_0(W(5),W(3),W(22),UV_GQQB_1EPS,AMPL(2,122))
      CALL FFV1_0(W(5),W(3),W(22),UV_GQQB_1EPS,AMPL(2,123))
      CALL FFV1_0(W(5),W(3),W(22),UV_GQQB_1EPS,AMPL(2,124))
      CALL FFV1_0(W(5),W(3),W(22),UV_GQQB_1EPS,AMPL(2,125))
      CALL FFV1_0(W(5),W(3),W(22),UV_GQQB_1EPS,AMPL(2,126))
      CALL FFV1_0(W(5),W(3),W(22),UV_GQQB_1EPS,AMPL(2,127))
      CALL FFV1_0(W(5),W(3),W(22),UV_GQQG_1EPS,AMPL(2,128))
      CALL FFV1_0(W(5),W(3),W(22),UV_GQQB,AMPL(1,129))
      CALL FFV1_0(W(5),W(3),W(22),UV_GQQT,AMPL(1,130))
C     Counter-term amplitude(s) for loop diagram number 71
      CALL FFV1_0(W(6),W(3),W(28),R2_GQQ,AMPL(1,131))
      CALL FFV1_0(W(6),W(3),W(28),UV_GQQB_1EPS,AMPL(2,132))
      CALL FFV1_0(W(6),W(3),W(28),UV_GQQB_1EPS,AMPL(2,133))
      CALL FFV1_0(W(6),W(3),W(28),UV_GQQB_1EPS,AMPL(2,134))
      CALL FFV1_0(W(6),W(3),W(28),UV_GQQB_1EPS,AMPL(2,135))
      CALL FFV1_0(W(6),W(3),W(28),UV_GQQB_1EPS,AMPL(2,136))
      CALL FFV1_0(W(6),W(3),W(28),UV_GQQB_1EPS,AMPL(2,137))
      CALL FFV1_0(W(6),W(3),W(28),UV_GQQG_1EPS,AMPL(2,138))
      CALL FFV1_0(W(6),W(3),W(28),UV_GQQB,AMPL(1,139))
      CALL FFV1_0(W(6),W(3),W(28),UV_GQQT,AMPL(1,140))
C     Counter-term amplitude(s) for loop diagram number 73
      CALL FFV1_0(W(6),W(3),W(30),R2_GQQ,AMPL(1,141))
      CALL FFV1_0(W(6),W(3),W(30),UV_GQQB_1EPS,AMPL(2,142))
      CALL FFV1_0(W(6),W(3),W(30),UV_GQQB_1EPS,AMPL(2,143))
      CALL FFV1_0(W(6),W(3),W(30),UV_GQQB_1EPS,AMPL(2,144))
      CALL FFV1_0(W(6),W(3),W(30),UV_GQQB_1EPS,AMPL(2,145))
      CALL FFV1_0(W(6),W(3),W(30),UV_GQQB_1EPS,AMPL(2,146))
      CALL FFV1_0(W(6),W(3),W(30),UV_GQQB_1EPS,AMPL(2,147))
      CALL FFV1_0(W(6),W(3),W(30),UV_GQQG_1EPS,AMPL(2,148))
      CALL FFV1_0(W(6),W(3),W(30),UV_GQQB,AMPL(1,149))
      CALL FFV1_0(W(6),W(3),W(30),UV_GQQT,AMPL(1,150))
C     Counter-term amplitude(s) for loop diagram number 75
      CALL FFV1_0(W(5),W(3),W(24),R2_GQQ,AMPL(1,151))
      CALL FFV1_0(W(5),W(3),W(24),UV_GQQB_1EPS,AMPL(2,152))
      CALL FFV1_0(W(5),W(3),W(24),UV_GQQB_1EPS,AMPL(2,153))
      CALL FFV1_0(W(5),W(3),W(24),UV_GQQB_1EPS,AMPL(2,154))
      CALL FFV1_0(W(5),W(3),W(24),UV_GQQB_1EPS,AMPL(2,155))
      CALL FFV1_0(W(5),W(3),W(24),UV_GQQB_1EPS,AMPL(2,156))
      CALL FFV1_0(W(5),W(3),W(24),UV_GQQB_1EPS,AMPL(2,157))
      CALL FFV1_0(W(5),W(3),W(24),UV_GQQG_1EPS,AMPL(2,158))
      CALL FFV1_0(W(5),W(3),W(24),UV_GQQB,AMPL(1,159))
      CALL FFV1_0(W(5),W(3),W(24),UV_GQQT,AMPL(1,160))
C     Counter-term amplitude(s) for loop diagram number 77
      CALL FFV1_0(W(5),W(3),W(25),R2_GQQ,AMPL(1,161))
      CALL FFV1_0(W(5),W(3),W(25),UV_GQQB_1EPS,AMPL(2,162))
      CALL FFV1_0(W(5),W(3),W(25),UV_GQQB_1EPS,AMPL(2,163))
      CALL FFV1_0(W(5),W(3),W(25),UV_GQQB_1EPS,AMPL(2,164))
      CALL FFV1_0(W(5),W(3),W(25),UV_GQQB_1EPS,AMPL(2,165))
      CALL FFV1_0(W(5),W(3),W(25),UV_GQQB_1EPS,AMPL(2,166))
      CALL FFV1_0(W(5),W(3),W(25),UV_GQQB_1EPS,AMPL(2,167))
      CALL FFV1_0(W(5),W(3),W(25),UV_GQQG_1EPS,AMPL(2,168))
      CALL FFV1_0(W(5),W(3),W(25),UV_GQQB,AMPL(1,169))
      CALL FFV1_0(W(5),W(3),W(25),UV_GQQT,AMPL(1,170))
C     Counter-term amplitude(s) for loop diagram number 79
      CALL FFV1_0(W(6),W(3),W(31),R2_GQQ,AMPL(1,171))
      CALL FFV1_0(W(6),W(3),W(31),UV_GQQB_1EPS,AMPL(2,172))
      CALL FFV1_0(W(6),W(3),W(31),UV_GQQB_1EPS,AMPL(2,173))
      CALL FFV1_0(W(6),W(3),W(31),UV_GQQB_1EPS,AMPL(2,174))
      CALL FFV1_0(W(6),W(3),W(31),UV_GQQB_1EPS,AMPL(2,175))
      CALL FFV1_0(W(6),W(3),W(31),UV_GQQB_1EPS,AMPL(2,176))
      CALL FFV1_0(W(6),W(3),W(31),UV_GQQB_1EPS,AMPL(2,177))
      CALL FFV1_0(W(6),W(3),W(31),UV_GQQG_1EPS,AMPL(2,178))
      CALL FFV1_0(W(6),W(3),W(31),UV_GQQB,AMPL(1,179))
      CALL FFV1_0(W(6),W(3),W(31),UV_GQQT,AMPL(1,180))
C     Counter-term amplitude(s) for loop diagram number 81
      CALL FFV1_0(W(6),W(3),W(32),R2_GQQ,AMPL(1,181))
      CALL FFV1_0(W(6),W(3),W(32),UV_GQQB_1EPS,AMPL(2,182))
      CALL FFV1_0(W(6),W(3),W(32),UV_GQQB_1EPS,AMPL(2,183))
      CALL FFV1_0(W(6),W(3),W(32),UV_GQQB_1EPS,AMPL(2,184))
      CALL FFV1_0(W(6),W(3),W(32),UV_GQQB_1EPS,AMPL(2,185))
      CALL FFV1_0(W(6),W(3),W(32),UV_GQQB_1EPS,AMPL(2,186))
      CALL FFV1_0(W(6),W(3),W(32),UV_GQQB_1EPS,AMPL(2,187))
      CALL FFV1_0(W(6),W(3),W(32),UV_GQQG_1EPS,AMPL(2,188))
      CALL FFV1_0(W(6),W(3),W(32),UV_GQQB,AMPL(1,189))
      CALL FFV1_0(W(6),W(3),W(32),UV_GQQT,AMPL(1,190))
C     Counter-term amplitude(s) for loop diagram number 83
      CALL FFV1_0(W(5),W(3),W(26),R2_GQQ,AMPL(1,191))
      CALL FFV1_0(W(5),W(3),W(26),UV_GQQB_1EPS,AMPL(2,192))
      CALL FFV1_0(W(5),W(3),W(26),UV_GQQB_1EPS,AMPL(2,193))
      CALL FFV1_0(W(5),W(3),W(26),UV_GQQB_1EPS,AMPL(2,194))
      CALL FFV1_0(W(5),W(3),W(26),UV_GQQB_1EPS,AMPL(2,195))
      CALL FFV1_0(W(5),W(3),W(26),UV_GQQB_1EPS,AMPL(2,196))
      CALL FFV1_0(W(5),W(3),W(26),UV_GQQB_1EPS,AMPL(2,197))
      CALL FFV1_0(W(5),W(3),W(26),UV_GQQG_1EPS,AMPL(2,198))
      CALL FFV1_0(W(5),W(3),W(26),UV_GQQB,AMPL(1,199))
      CALL FFV1_0(W(5),W(3),W(26),UV_GQQT,AMPL(1,200))
      CALL FFV1_2(W(6),W(17),GC_5,MDL_MB,ZERO,W(33))
C     Counter-term amplitude(s) for loop diagram number 85
      CALL FFV1_0(W(33),W(3),W(7),R2_DDA,AMPL(1,201))
C     Counter-term amplitude(s) for loop diagram number 87
      CALL FFV1_0(W(10),W(3),W(17),R2_GQQ,AMPL(1,202))
      CALL FFV1_0(W(10),W(3),W(17),UV_GQQB_1EPS,AMPL(2,203))
      CALL FFV1_0(W(10),W(3),W(17),UV_GQQB_1EPS,AMPL(2,204))
      CALL FFV1_0(W(10),W(3),W(17),UV_GQQB_1EPS,AMPL(2,205))
      CALL FFV1_0(W(10),W(3),W(17),UV_GQQB_1EPS,AMPL(2,206))
      CALL FFV1_0(W(10),W(3),W(17),UV_GQQB_1EPS,AMPL(2,207))
      CALL FFV1_0(W(10),W(3),W(17),UV_GQQB_1EPS,AMPL(2,208))
      CALL FFV1_0(W(10),W(3),W(17),UV_GQQG_1EPS,AMPL(2,209))
      CALL FFV1_0(W(10),W(3),W(17),UV_GQQB,AMPL(1,210))
      CALL FFV1_0(W(10),W(3),W(17),UV_GQQT,AMPL(1,211))
C     Counter-term amplitude(s) for loop diagram number 89
      CALL FFV2_3_0(W(33),W(3),W(11),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,212))
C     Counter-term amplitude(s) for loop diagram number 91
      CALL FFV1_0(W(13),W(3),W(17),R2_GQQ,AMPL(1,213))
      CALL FFV1_0(W(13),W(3),W(17),UV_GQQB_1EPS,AMPL(2,214))
      CALL FFV1_0(W(13),W(3),W(17),UV_GQQB_1EPS,AMPL(2,215))
      CALL FFV1_0(W(13),W(3),W(17),UV_GQQB_1EPS,AMPL(2,216))
      CALL FFV1_0(W(13),W(3),W(17),UV_GQQB_1EPS,AMPL(2,217))
      CALL FFV1_0(W(13),W(3),W(17),UV_GQQB_1EPS,AMPL(2,218))
      CALL FFV1_0(W(13),W(3),W(17),UV_GQQB_1EPS,AMPL(2,219))
      CALL FFV1_0(W(13),W(3),W(17),UV_GQQG_1EPS,AMPL(2,220))
      CALL FFV1_0(W(13),W(3),W(17),UV_GQQB,AMPL(1,221))
      CALL FFV1_0(W(13),W(3),W(17),UV_GQQT,AMPL(1,222))
      CALL FFV1_2(W(5),W(20),GC_5,MDL_MB,ZERO,W(34))
C     Counter-term amplitude(s) for loop diagram number 97
      CALL FFV1_0(W(34),W(3),W(7),R2_DDA,AMPL(1,223))
C     Counter-term amplitude(s) for loop diagram number 99
      CALL FFV1_0(W(15),W(3),W(20),R2_GQQ,AMPL(1,224))
      CALL FFV1_0(W(15),W(3),W(20),UV_GQQB_1EPS,AMPL(2,225))
      CALL FFV1_0(W(15),W(3),W(20),UV_GQQB_1EPS,AMPL(2,226))
      CALL FFV1_0(W(15),W(3),W(20),UV_GQQB_1EPS,AMPL(2,227))
      CALL FFV1_0(W(15),W(3),W(20),UV_GQQB_1EPS,AMPL(2,228))
      CALL FFV1_0(W(15),W(3),W(20),UV_GQQB_1EPS,AMPL(2,229))
      CALL FFV1_0(W(15),W(3),W(20),UV_GQQB_1EPS,AMPL(2,230))
      CALL FFV1_0(W(15),W(3),W(20),UV_GQQG_1EPS,AMPL(2,231))
      CALL FFV1_0(W(15),W(3),W(20),UV_GQQB,AMPL(1,232))
      CALL FFV1_0(W(15),W(3),W(20),UV_GQQT,AMPL(1,233))
C     Counter-term amplitude(s) for loop diagram number 101
      CALL FFV2_3_0(W(34),W(3),W(11),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,234))
C     Counter-term amplitude(s) for loop diagram number 103
      CALL FFV1_0(W(16),W(3),W(20),R2_GQQ,AMPL(1,235))
      CALL FFV1_0(W(16),W(3),W(20),UV_GQQB_1EPS,AMPL(2,236))
      CALL FFV1_0(W(16),W(3),W(20),UV_GQQB_1EPS,AMPL(2,237))
      CALL FFV1_0(W(16),W(3),W(20),UV_GQQB_1EPS,AMPL(2,238))
      CALL FFV1_0(W(16),W(3),W(20),UV_GQQB_1EPS,AMPL(2,239))
      CALL FFV1_0(W(16),W(3),W(20),UV_GQQB_1EPS,AMPL(2,240))
      CALL FFV1_0(W(16),W(3),W(20),UV_GQQB_1EPS,AMPL(2,241))
      CALL FFV1_0(W(16),W(3),W(20),UV_GQQG_1EPS,AMPL(2,242))
      CALL FFV1_0(W(16),W(3),W(20),UV_GQQB,AMPL(1,243))
      CALL FFV1_0(W(16),W(3),W(20),UV_GQQT,AMPL(1,244))
C     Counter-term amplitude(s) for loop diagram number 149
      CALL R2_QQ_1_R2_QQ_2_0(W(33),W(18),R2_QQQ,R2_QQB,AMPL(1,245))
      CALL R2_QQ_2_0(W(33),W(18),UV_BMASS_1EPS,AMPL(2,246))
      CALL R2_QQ_2_0(W(33),W(18),UV_BMASS,AMPL(1,247))
      CALL FFV1P0_3(W(6),W(18),GC_5,ZERO,ZERO,W(35))
C     Counter-term amplitude(s) for loop diagram number 150
      CALL R2_GG_1_R2_GG_2_0(W(17),W(35),R2_GGG_1,R2_GGG_2,AMPL(1,248))
C     Counter-term amplitude(s) for loop diagram number 151
      CALL FFV1_0(W(6),W(18),W(17),R2_GQQ,AMPL(1,249))
      CALL FFV1_0(W(6),W(18),W(17),UV_GQQB_1EPS,AMPL(2,250))
      CALL FFV1_0(W(6),W(18),W(17),UV_GQQB_1EPS,AMPL(2,251))
      CALL FFV1_0(W(6),W(18),W(17),UV_GQQB_1EPS,AMPL(2,252))
      CALL FFV1_0(W(6),W(18),W(17),UV_GQQB_1EPS,AMPL(2,253))
      CALL FFV1_0(W(6),W(18),W(17),UV_GQQB_1EPS,AMPL(2,254))
      CALL FFV1_0(W(6),W(18),W(17),UV_GQQB_1EPS,AMPL(2,255))
      CALL FFV1_0(W(6),W(18),W(17),UV_GQQG_1EPS,AMPL(2,256))
      CALL FFV1_0(W(6),W(18),W(17),UV_GQQB,AMPL(1,257))
      CALL FFV1_0(W(6),W(18),W(17),UV_GQQT,AMPL(1,258))
      CALL FFV1_1(W(3),W(17),GC_5,MDL_MB,ZERO,W(36))
C     Counter-term amplitude(s) for loop diagram number 152
      CALL R2_QQ_1_R2_QQ_2_0(W(10),W(36),R2_QQQ,R2_QQB,AMPL(1,259))
      CALL R2_QQ_2_0(W(10),W(36),UV_BMASS_1EPS,AMPL(2,260))
      CALL R2_QQ_2_0(W(10),W(36),UV_BMASS,AMPL(1,261))
      CALL FFV1P0_3(W(10),W(3),GC_5,ZERO,ZERO,W(37))
C     Counter-term amplitude(s) for loop diagram number 153
      CALL R2_GG_1_R2_GG_2_0(W(17),W(37),R2_GGG_1,R2_GGG_2,AMPL(1,262))
C     Counter-term amplitude(s) for loop diagram number 154
      CALL R2_QQ_1_R2_QQ_2_0(W(33),W(19),R2_QQQ,R2_QQB,AMPL(1,263))
      CALL R2_QQ_2_0(W(33),W(19),UV_BMASS_1EPS,AMPL(2,264))
      CALL R2_QQ_2_0(W(33),W(19),UV_BMASS,AMPL(1,265))
      CALL FFV1P0_3(W(6),W(19),GC_5,ZERO,ZERO,W(38))
C     Counter-term amplitude(s) for loop diagram number 155
      CALL R2_GG_1_R2_GG_2_0(W(17),W(38),R2_GGG_1,R2_GGG_2,AMPL(1,266))
C     Counter-term amplitude(s) for loop diagram number 156
      CALL FFV1_0(W(6),W(19),W(17),R2_GQQ,AMPL(1,267))
      CALL FFV1_0(W(6),W(19),W(17),UV_GQQB_1EPS,AMPL(2,268))
      CALL FFV1_0(W(6),W(19),W(17),UV_GQQB_1EPS,AMPL(2,269))
      CALL FFV1_0(W(6),W(19),W(17),UV_GQQB_1EPS,AMPL(2,270))
      CALL FFV1_0(W(6),W(19),W(17),UV_GQQB_1EPS,AMPL(2,271))
      CALL FFV1_0(W(6),W(19),W(17),UV_GQQB_1EPS,AMPL(2,272))
      CALL FFV1_0(W(6),W(19),W(17),UV_GQQB_1EPS,AMPL(2,273))
      CALL FFV1_0(W(6),W(19),W(17),UV_GQQG_1EPS,AMPL(2,274))
      CALL FFV1_0(W(6),W(19),W(17),UV_GQQB,AMPL(1,275))
      CALL FFV1_0(W(6),W(19),W(17),UV_GQQT,AMPL(1,276))
C     Counter-term amplitude(s) for loop diagram number 157
      CALL R2_QQ_1_R2_QQ_2_0(W(13),W(36),R2_QQQ,R2_QQB,AMPL(1,277))
      CALL R2_QQ_2_0(W(13),W(36),UV_BMASS_1EPS,AMPL(2,278))
      CALL R2_QQ_2_0(W(13),W(36),UV_BMASS,AMPL(1,279))
      CALL FFV1P0_3(W(13),W(3),GC_5,ZERO,ZERO,W(39))
C     Counter-term amplitude(s) for loop diagram number 158
      CALL R2_GG_1_R2_GG_2_0(W(17),W(39),R2_GGG_1,R2_GGG_2,AMPL(1,280))
C     Counter-term amplitude(s) for loop diagram number 160
      CALL FFV1_0(W(6),W(36),W(7),R2_DDA,AMPL(1,281))
C     Counter-term amplitude(s) for loop diagram number 162
      CALL FFV2_3_0(W(6),W(36),W(11),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,282))
C     Counter-term amplitude(s) for loop diagram number 163
      CALL R2_QQ_1_R2_QQ_2_0(W(34),W(18),R2_QQQ,R2_QQB,AMPL(1,283))
      CALL R2_QQ_2_0(W(34),W(18),UV_BMASS_1EPS,AMPL(2,284))
      CALL R2_QQ_2_0(W(34),W(18),UV_BMASS,AMPL(1,285))
      CALL FFV1P0_3(W(5),W(18),GC_5,ZERO,ZERO,W(40))
C     Counter-term amplitude(s) for loop diagram number 164
      CALL R2_GG_1_R2_GG_2_0(W(20),W(40),R2_GGG_1,R2_GGG_2,AMPL(1,286))
C     Counter-term amplitude(s) for loop diagram number 165
      CALL FFV1_0(W(5),W(18),W(20),R2_GQQ,AMPL(1,287))
      CALL FFV1_0(W(5),W(18),W(20),UV_GQQB_1EPS,AMPL(2,288))
      CALL FFV1_0(W(5),W(18),W(20),UV_GQQB_1EPS,AMPL(2,289))
      CALL FFV1_0(W(5),W(18),W(20),UV_GQQB_1EPS,AMPL(2,290))
      CALL FFV1_0(W(5),W(18),W(20),UV_GQQB_1EPS,AMPL(2,291))
      CALL FFV1_0(W(5),W(18),W(20),UV_GQQB_1EPS,AMPL(2,292))
      CALL FFV1_0(W(5),W(18),W(20),UV_GQQB_1EPS,AMPL(2,293))
      CALL FFV1_0(W(5),W(18),W(20),UV_GQQG_1EPS,AMPL(2,294))
      CALL FFV1_0(W(5),W(18),W(20),UV_GQQB,AMPL(1,295))
      CALL FFV1_0(W(5),W(18),W(20),UV_GQQT,AMPL(1,296))
      CALL FFV1_1(W(3),W(20),GC_5,MDL_MB,ZERO,W(41))
C     Counter-term amplitude(s) for loop diagram number 166
      CALL R2_QQ_1_R2_QQ_2_0(W(15),W(41),R2_QQQ,R2_QQB,AMPL(1,297))
      CALL R2_QQ_2_0(W(15),W(41),UV_BMASS_1EPS,AMPL(2,298))
      CALL R2_QQ_2_0(W(15),W(41),UV_BMASS,AMPL(1,299))
      CALL FFV1P0_3(W(15),W(3),GC_5,ZERO,ZERO,W(42))
C     Counter-term amplitude(s) for loop diagram number 167
      CALL R2_GG_1_R2_GG_2_0(W(20),W(42),R2_GGG_1,R2_GGG_2,AMPL(1,300))
C     Counter-term amplitude(s) for loop diagram number 168
      CALL R2_QQ_1_R2_QQ_2_0(W(34),W(19),R2_QQQ,R2_QQB,AMPL(1,301))
      CALL R2_QQ_2_0(W(34),W(19),UV_BMASS_1EPS,AMPL(2,302))
      CALL R2_QQ_2_0(W(34),W(19),UV_BMASS,AMPL(1,303))
      CALL FFV1P0_3(W(5),W(19),GC_5,ZERO,ZERO,W(43))
C     Counter-term amplitude(s) for loop diagram number 169
      CALL R2_GG_1_R2_GG_2_0(W(20),W(43),R2_GGG_1,R2_GGG_2,AMPL(1,304))
C     Counter-term amplitude(s) for loop diagram number 170
      CALL FFV1_0(W(5),W(19),W(20),R2_GQQ,AMPL(1,305))
      CALL FFV1_0(W(5),W(19),W(20),UV_GQQB_1EPS,AMPL(2,306))
      CALL FFV1_0(W(5),W(19),W(20),UV_GQQB_1EPS,AMPL(2,307))
      CALL FFV1_0(W(5),W(19),W(20),UV_GQQB_1EPS,AMPL(2,308))
      CALL FFV1_0(W(5),W(19),W(20),UV_GQQB_1EPS,AMPL(2,309))
      CALL FFV1_0(W(5),W(19),W(20),UV_GQQB_1EPS,AMPL(2,310))
      CALL FFV1_0(W(5),W(19),W(20),UV_GQQB_1EPS,AMPL(2,311))
      CALL FFV1_0(W(5),W(19),W(20),UV_GQQG_1EPS,AMPL(2,312))
      CALL FFV1_0(W(5),W(19),W(20),UV_GQQB,AMPL(1,313))
      CALL FFV1_0(W(5),W(19),W(20),UV_GQQT,AMPL(1,314))
C     Counter-term amplitude(s) for loop diagram number 171
      CALL R2_QQ_1_R2_QQ_2_0(W(16),W(41),R2_QQQ,R2_QQB,AMPL(1,315))
      CALL R2_QQ_2_0(W(16),W(41),UV_BMASS_1EPS,AMPL(2,316))
      CALL R2_QQ_2_0(W(16),W(41),UV_BMASS,AMPL(1,317))
      CALL FFV1P0_3(W(16),W(3),GC_5,ZERO,ZERO,W(44))
C     Counter-term amplitude(s) for loop diagram number 172
      CALL R2_GG_1_R2_GG_2_0(W(20),W(44),R2_GGG_1,R2_GGG_2,AMPL(1,318))
C     Counter-term amplitude(s) for loop diagram number 174
      CALL FFV1_0(W(5),W(41),W(7),R2_DDA,AMPL(1,319))
C     Counter-term amplitude(s) for loop diagram number 176
      CALL FFV2_3_0(W(5),W(41),W(11),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,320))
C     Counter-term amplitude(s) for loop diagram number 177
      CALL FFV1_0(W(5),W(4),W(35),R2_GQQ,AMPL(1,321))
      CALL FFV1_0(W(5),W(4),W(35),UV_GQQB_1EPS,AMPL(2,322))
      CALL FFV1_0(W(5),W(4),W(35),UV_GQQB_1EPS,AMPL(2,323))
      CALL FFV1_0(W(5),W(4),W(35),UV_GQQB_1EPS,AMPL(2,324))
      CALL FFV1_0(W(5),W(4),W(35),UV_GQQB_1EPS,AMPL(2,325))
      CALL FFV1_0(W(5),W(4),W(35),UV_GQQB_1EPS,AMPL(2,326))
      CALL FFV1_0(W(5),W(4),W(35),UV_GQQB_1EPS,AMPL(2,327))
      CALL FFV1_0(W(5),W(4),W(35),UV_GQQG_1EPS,AMPL(2,328))
      CALL FFV1_0(W(5),W(4),W(35),UV_GQQB,AMPL(1,329))
      CALL FFV1_0(W(5),W(4),W(35),UV_GQQT,AMPL(1,330))
C     Counter-term amplitude(s) for loop diagram number 179
      CALL FFV1_0(W(6),W(4),W(40),R2_GQQ,AMPL(1,331))
      CALL FFV1_0(W(6),W(4),W(40),UV_GQQB_1EPS,AMPL(2,332))
      CALL FFV1_0(W(6),W(4),W(40),UV_GQQB_1EPS,AMPL(2,333))
      CALL FFV1_0(W(6),W(4),W(40),UV_GQQB_1EPS,AMPL(2,334))
      CALL FFV1_0(W(6),W(4),W(40),UV_GQQB_1EPS,AMPL(2,335))
      CALL FFV1_0(W(6),W(4),W(40),UV_GQQB_1EPS,AMPL(2,336))
      CALL FFV1_0(W(6),W(4),W(40),UV_GQQB_1EPS,AMPL(2,337))
      CALL FFV1_0(W(6),W(4),W(40),UV_GQQG_1EPS,AMPL(2,338))
      CALL FFV1_0(W(6),W(4),W(40),UV_GQQB,AMPL(1,339))
      CALL FFV1_0(W(6),W(4),W(40),UV_GQQT,AMPL(1,340))
C     Counter-term amplitude(s) for loop diagram number 181
      CALL FFV1_0(W(6),W(4),W(42),R2_GQQ,AMPL(1,341))
      CALL FFV1_0(W(6),W(4),W(42),UV_GQQB_1EPS,AMPL(2,342))
      CALL FFV1_0(W(6),W(4),W(42),UV_GQQB_1EPS,AMPL(2,343))
      CALL FFV1_0(W(6),W(4),W(42),UV_GQQB_1EPS,AMPL(2,344))
      CALL FFV1_0(W(6),W(4),W(42),UV_GQQB_1EPS,AMPL(2,345))
      CALL FFV1_0(W(6),W(4),W(42),UV_GQQB_1EPS,AMPL(2,346))
      CALL FFV1_0(W(6),W(4),W(42),UV_GQQB_1EPS,AMPL(2,347))
      CALL FFV1_0(W(6),W(4),W(42),UV_GQQG_1EPS,AMPL(2,348))
      CALL FFV1_0(W(6),W(4),W(42),UV_GQQB,AMPL(1,349))
      CALL FFV1_0(W(6),W(4),W(42),UV_GQQT,AMPL(1,350))
C     Counter-term amplitude(s) for loop diagram number 182
      CALL FFV1_0(W(5),W(4),W(37),R2_GQQ,AMPL(1,351))
      CALL FFV1_0(W(5),W(4),W(37),UV_GQQB_1EPS,AMPL(2,352))
      CALL FFV1_0(W(5),W(4),W(37),UV_GQQB_1EPS,AMPL(2,353))
      CALL FFV1_0(W(5),W(4),W(37),UV_GQQB_1EPS,AMPL(2,354))
      CALL FFV1_0(W(5),W(4),W(37),UV_GQQB_1EPS,AMPL(2,355))
      CALL FFV1_0(W(5),W(4),W(37),UV_GQQB_1EPS,AMPL(2,356))
      CALL FFV1_0(W(5),W(4),W(37),UV_GQQB_1EPS,AMPL(2,357))
      CALL FFV1_0(W(5),W(4),W(37),UV_GQQG_1EPS,AMPL(2,358))
      CALL FFV1_0(W(5),W(4),W(37),UV_GQQB,AMPL(1,359))
      CALL FFV1_0(W(5),W(4),W(37),UV_GQQT,AMPL(1,360))
C     Counter-term amplitude(s) for loop diagram number 183
      CALL FFV1_0(W(5),W(4),W(38),R2_GQQ,AMPL(1,361))
      CALL FFV1_0(W(5),W(4),W(38),UV_GQQB_1EPS,AMPL(2,362))
      CALL FFV1_0(W(5),W(4),W(38),UV_GQQB_1EPS,AMPL(2,363))
      CALL FFV1_0(W(5),W(4),W(38),UV_GQQB_1EPS,AMPL(2,364))
      CALL FFV1_0(W(5),W(4),W(38),UV_GQQB_1EPS,AMPL(2,365))
      CALL FFV1_0(W(5),W(4),W(38),UV_GQQB_1EPS,AMPL(2,366))
      CALL FFV1_0(W(5),W(4),W(38),UV_GQQB_1EPS,AMPL(2,367))
      CALL FFV1_0(W(5),W(4),W(38),UV_GQQG_1EPS,AMPL(2,368))
      CALL FFV1_0(W(5),W(4),W(38),UV_GQQB,AMPL(1,369))
      CALL FFV1_0(W(5),W(4),W(38),UV_GQQT,AMPL(1,370))
C     Counter-term amplitude(s) for loop diagram number 185
      CALL FFV1_0(W(6),W(4),W(43),R2_GQQ,AMPL(1,371))
      CALL FFV1_0(W(6),W(4),W(43),UV_GQQB_1EPS,AMPL(2,372))
      CALL FFV1_0(W(6),W(4),W(43),UV_GQQB_1EPS,AMPL(2,373))
      CALL FFV1_0(W(6),W(4),W(43),UV_GQQB_1EPS,AMPL(2,374))
      CALL FFV1_0(W(6),W(4),W(43),UV_GQQB_1EPS,AMPL(2,375))
      CALL FFV1_0(W(6),W(4),W(43),UV_GQQB_1EPS,AMPL(2,376))
      CALL FFV1_0(W(6),W(4),W(43),UV_GQQB_1EPS,AMPL(2,377))
      CALL FFV1_0(W(6),W(4),W(43),UV_GQQG_1EPS,AMPL(2,378))
      CALL FFV1_0(W(6),W(4),W(43),UV_GQQB,AMPL(1,379))
      CALL FFV1_0(W(6),W(4),W(43),UV_GQQT,AMPL(1,380))
C     Counter-term amplitude(s) for loop diagram number 187
      CALL FFV1_0(W(6),W(4),W(44),R2_GQQ,AMPL(1,381))
      CALL FFV1_0(W(6),W(4),W(44),UV_GQQB_1EPS,AMPL(2,382))
      CALL FFV1_0(W(6),W(4),W(44),UV_GQQB_1EPS,AMPL(2,383))
      CALL FFV1_0(W(6),W(4),W(44),UV_GQQB_1EPS,AMPL(2,384))
      CALL FFV1_0(W(6),W(4),W(44),UV_GQQB_1EPS,AMPL(2,385))
      CALL FFV1_0(W(6),W(4),W(44),UV_GQQB_1EPS,AMPL(2,386))
      CALL FFV1_0(W(6),W(4),W(44),UV_GQQB_1EPS,AMPL(2,387))
      CALL FFV1_0(W(6),W(4),W(44),UV_GQQG_1EPS,AMPL(2,388))
      CALL FFV1_0(W(6),W(4),W(44),UV_GQQB,AMPL(1,389))
      CALL FFV1_0(W(6),W(4),W(44),UV_GQQT,AMPL(1,390))
C     Counter-term amplitude(s) for loop diagram number 188
      CALL FFV1_0(W(5),W(4),W(39),R2_GQQ,AMPL(1,391))
      CALL FFV1_0(W(5),W(4),W(39),UV_GQQB_1EPS,AMPL(2,392))
      CALL FFV1_0(W(5),W(4),W(39),UV_GQQB_1EPS,AMPL(2,393))
      CALL FFV1_0(W(5),W(4),W(39),UV_GQQB_1EPS,AMPL(2,394))
      CALL FFV1_0(W(5),W(4),W(39),UV_GQQB_1EPS,AMPL(2,395))
      CALL FFV1_0(W(5),W(4),W(39),UV_GQQB_1EPS,AMPL(2,396))
      CALL FFV1_0(W(5),W(4),W(39),UV_GQQB_1EPS,AMPL(2,397))
      CALL FFV1_0(W(5),W(4),W(39),UV_GQQG_1EPS,AMPL(2,398))
      CALL FFV1_0(W(5),W(4),W(39),UV_GQQB,AMPL(1,399))
      CALL FFV1_0(W(5),W(4),W(39),UV_GQQT,AMPL(1,400))
C     Counter-term amplitude(s) for loop diagram number 217
      CALL R2_GG_1_0(W(8),W(22),R2_GGQ,AMPL(1,401))
      CALL R2_GG_1_0(W(8),W(22),R2_GGQ,AMPL(1,402))
      CALL R2_GG_1_0(W(8),W(22),R2_GGQ,AMPL(1,403))
      CALL R2_GG_1_0(W(8),W(22),R2_GGQ,AMPL(1,404))
C     Counter-term amplitude(s) for loop diagram number 218
      CALL R2_GG_1_0(W(8),W(24),R2_GGQ,AMPL(1,405))
      CALL R2_GG_1_0(W(8),W(24),R2_GGQ,AMPL(1,406))
      CALL R2_GG_1_0(W(8),W(24),R2_GGQ,AMPL(1,407))
      CALL R2_GG_1_0(W(8),W(24),R2_GGQ,AMPL(1,408))
C     Counter-term amplitude(s) for loop diagram number 219
      CALL R2_GG_1_0(W(8),W(25),R2_GGQ,AMPL(1,409))
      CALL R2_GG_1_0(W(8),W(25),R2_GGQ,AMPL(1,410))
      CALL R2_GG_1_0(W(8),W(25),R2_GGQ,AMPL(1,411))
      CALL R2_GG_1_0(W(8),W(25),R2_GGQ,AMPL(1,412))
C     Counter-term amplitude(s) for loop diagram number 220
      CALL R2_GG_1_0(W(8),W(26),R2_GGQ,AMPL(1,413))
      CALL R2_GG_1_0(W(8),W(26),R2_GGQ,AMPL(1,414))
      CALL R2_GG_1_0(W(8),W(26),R2_GGQ,AMPL(1,415))
      CALL R2_GG_1_0(W(8),W(26),R2_GGQ,AMPL(1,416))
C     Counter-term amplitude(s) for loop diagram number 221
      CALL R2_GGZ_0(W(8),W(20),W(11),R2_GGZDOWN,AMPL(1,417))
      CALL R2_GGZ_0(W(8),W(20),W(11),R2_GGZDOWN,AMPL(1,418))
C     Counter-term amplitude(s) for loop diagram number 223
      CALL R2_GG_1_0(W(14),W(28),R2_GGQ,AMPL(1,419))
      CALL R2_GG_1_0(W(14),W(28),R2_GGQ,AMPL(1,420))
      CALL R2_GG_1_0(W(14),W(28),R2_GGQ,AMPL(1,421))
      CALL R2_GG_1_0(W(14),W(28),R2_GGQ,AMPL(1,422))
C     Counter-term amplitude(s) for loop diagram number 224
      CALL R2_GG_1_0(W(14),W(30),R2_GGQ,AMPL(1,423))
      CALL R2_GG_1_0(W(14),W(30),R2_GGQ,AMPL(1,424))
      CALL R2_GG_1_0(W(14),W(30),R2_GGQ,AMPL(1,425))
      CALL R2_GG_1_0(W(14),W(30),R2_GGQ,AMPL(1,426))
C     Counter-term amplitude(s) for loop diagram number 225
      CALL R2_GG_1_0(W(14),W(31),R2_GGQ,AMPL(1,427))
      CALL R2_GG_1_0(W(14),W(31),R2_GGQ,AMPL(1,428))
      CALL R2_GG_1_0(W(14),W(31),R2_GGQ,AMPL(1,429))
      CALL R2_GG_1_0(W(14),W(31),R2_GGQ,AMPL(1,430))
C     Counter-term amplitude(s) for loop diagram number 226
      CALL R2_GG_1_0(W(14),W(32),R2_GGQ,AMPL(1,431))
      CALL R2_GG_1_0(W(14),W(32),R2_GGQ,AMPL(1,432))
      CALL R2_GG_1_0(W(14),W(32),R2_GGQ,AMPL(1,433))
      CALL R2_GG_1_0(W(14),W(32),R2_GGQ,AMPL(1,434))
C     Counter-term amplitude(s) for loop diagram number 227
      CALL R2_GGZ_0(W(17),W(14),W(11),R2_GGZDOWN,AMPL(1,435))
      CALL R2_GGZ_0(W(17),W(14),W(11),R2_GGZDOWN,AMPL(1,436))
C     Counter-term amplitude(s) for loop diagram number 229
      CALL R2_GG_1_0(W(17),W(35),R2_GGQ,AMPL(1,437))
      CALL R2_GG_1_0(W(17),W(35),R2_GGQ,AMPL(1,438))
      CALL R2_GG_1_0(W(17),W(35),R2_GGQ,AMPL(1,439))
      CALL R2_GG_1_0(W(17),W(35),R2_GGQ,AMPL(1,440))
C     Counter-term amplitude(s) for loop diagram number 230
      CALL R2_GG_1_0(W(17),W(37),R2_GGQ,AMPL(1,441))
      CALL R2_GG_1_0(W(17),W(37),R2_GGQ,AMPL(1,442))
      CALL R2_GG_1_0(W(17),W(37),R2_GGQ,AMPL(1,443))
      CALL R2_GG_1_0(W(17),W(37),R2_GGQ,AMPL(1,444))
C     Counter-term amplitude(s) for loop diagram number 231
      CALL R2_GG_1_0(W(17),W(38),R2_GGQ,AMPL(1,445))
      CALL R2_GG_1_0(W(17),W(38),R2_GGQ,AMPL(1,446))
      CALL R2_GG_1_0(W(17),W(38),R2_GGQ,AMPL(1,447))
      CALL R2_GG_1_0(W(17),W(38),R2_GGQ,AMPL(1,448))
C     Counter-term amplitude(s) for loop diagram number 232
      CALL R2_GG_1_0(W(17),W(39),R2_GGQ,AMPL(1,449))
      CALL R2_GG_1_0(W(17),W(39),R2_GGQ,AMPL(1,450))
      CALL R2_GG_1_0(W(17),W(39),R2_GGQ,AMPL(1,451))
      CALL R2_GG_1_0(W(17),W(39),R2_GGQ,AMPL(1,452))
C     Counter-term amplitude(s) for loop diagram number 233
      CALL R2_GG_1_0(W(20),W(40),R2_GGQ,AMPL(1,453))
      CALL R2_GG_1_0(W(20),W(40),R2_GGQ,AMPL(1,454))
      CALL R2_GG_1_0(W(20),W(40),R2_GGQ,AMPL(1,455))
      CALL R2_GG_1_0(W(20),W(40),R2_GGQ,AMPL(1,456))
C     Counter-term amplitude(s) for loop diagram number 234
      CALL R2_GG_1_0(W(20),W(42),R2_GGQ,AMPL(1,457))
      CALL R2_GG_1_0(W(20),W(42),R2_GGQ,AMPL(1,458))
      CALL R2_GG_1_0(W(20),W(42),R2_GGQ,AMPL(1,459))
      CALL R2_GG_1_0(W(20),W(42),R2_GGQ,AMPL(1,460))
C     Counter-term amplitude(s) for loop diagram number 235
      CALL R2_GG_1_0(W(20),W(43),R2_GGQ,AMPL(1,461))
      CALL R2_GG_1_0(W(20),W(43),R2_GGQ,AMPL(1,462))
      CALL R2_GG_1_0(W(20),W(43),R2_GGQ,AMPL(1,463))
      CALL R2_GG_1_0(W(20),W(43),R2_GGQ,AMPL(1,464))
C     Counter-term amplitude(s) for loop diagram number 236
      CALL R2_GG_1_0(W(20),W(44),R2_GGQ,AMPL(1,465))
      CALL R2_GG_1_0(W(20),W(44),R2_GGQ,AMPL(1,466))
      CALL R2_GG_1_0(W(20),W(44),R2_GGQ,AMPL(1,467))
      CALL R2_GG_1_0(W(20),W(44),R2_GGQ,AMPL(1,468))
C     Counter-term amplitude(s) for loop diagram number 237
      CALL R2_GGZ_0(W(8),W(20),W(11),R2_GGZUP,AMPL(1,469))
      CALL R2_GGZ_0(W(8),W(20),W(11),R2_GGZUP,AMPL(1,470))
C     Counter-term amplitude(s) for loop diagram number 239
      CALL R2_GGZ_0(W(17),W(14),W(11),R2_GGZUP,AMPL(1,471))
      CALL R2_GGZ_0(W(17),W(14),W(11),R2_GGZUP,AMPL(1,472))
C     Counter-term amplitude(s) for loop diagram number 241
      CALL R2_GG_1_R2_GG_3_0(W(8),W(22),R2_GGQ,R2_GGB,AMPL(1,473))
C     Counter-term amplitude(s) for loop diagram number 242
      CALL R2_GG_1_R2_GG_3_0(W(8),W(24),R2_GGQ,R2_GGB,AMPL(1,474))
C     Counter-term amplitude(s) for loop diagram number 243
      CALL R2_GG_1_R2_GG_3_0(W(8),W(25),R2_GGQ,R2_GGB,AMPL(1,475))
C     Counter-term amplitude(s) for loop diagram number 244
      CALL R2_GG_1_R2_GG_3_0(W(8),W(26),R2_GGQ,R2_GGB,AMPL(1,476))
C     Counter-term amplitude(s) for loop diagram number 245
      CALL R2_GGZ_0(W(8),W(20),W(11),R2_GGZDOWN,AMPL(1,477))
C     Counter-term amplitude(s) for loop diagram number 247
      CALL R2_GG_1_R2_GG_3_0(W(14),W(28),R2_GGQ,R2_GGB,AMPL(1,478))
C     Counter-term amplitude(s) for loop diagram number 248
      CALL R2_GG_1_R2_GG_3_0(W(14),W(30),R2_GGQ,R2_GGB,AMPL(1,479))
C     Counter-term amplitude(s) for loop diagram number 249
      CALL R2_GG_1_R2_GG_3_0(W(14),W(31),R2_GGQ,R2_GGB,AMPL(1,480))
C     Counter-term amplitude(s) for loop diagram number 250
      CALL R2_GG_1_R2_GG_3_0(W(14),W(32),R2_GGQ,R2_GGB,AMPL(1,481))
C     Counter-term amplitude(s) for loop diagram number 251
      CALL R2_GGZ_0(W(17),W(14),W(11),R2_GGZDOWN,AMPL(1,482))
C     Counter-term amplitude(s) for loop diagram number 253
      CALL R2_GG_1_R2_GG_3_0(W(17),W(35),R2_GGQ,R2_GGB,AMPL(1,483))
C     Counter-term amplitude(s) for loop diagram number 254
      CALL R2_GG_1_R2_GG_3_0(W(17),W(37),R2_GGQ,R2_GGB,AMPL(1,484))
C     Counter-term amplitude(s) for loop diagram number 255
      CALL R2_GG_1_R2_GG_3_0(W(17),W(38),R2_GGQ,R2_GGB,AMPL(1,485))
C     Counter-term amplitude(s) for loop diagram number 256
      CALL R2_GG_1_R2_GG_3_0(W(17),W(39),R2_GGQ,R2_GGB,AMPL(1,486))
C     Counter-term amplitude(s) for loop diagram number 257
      CALL R2_GG_1_R2_GG_3_0(W(20),W(40),R2_GGQ,R2_GGB,AMPL(1,487))
C     Counter-term amplitude(s) for loop diagram number 258
      CALL R2_GG_1_R2_GG_3_0(W(20),W(42),R2_GGQ,R2_GGB,AMPL(1,488))
C     Counter-term amplitude(s) for loop diagram number 259
      CALL R2_GG_1_R2_GG_3_0(W(20),W(43),R2_GGQ,R2_GGB,AMPL(1,489))
C     Counter-term amplitude(s) for loop diagram number 260
      CALL R2_GG_1_R2_GG_3_0(W(20),W(44),R2_GGQ,R2_GGB,AMPL(1,490))
C     Counter-term amplitude(s) for loop diagram number 261
      CALL R2_GG_1_R2_GG_3_0(W(8),W(22),R2_GGQ,R2_GGT,AMPL(1,491))
C     Counter-term amplitude(s) for loop diagram number 262
      CALL R2_GG_1_R2_GG_3_0(W(8),W(24),R2_GGQ,R2_GGT,AMPL(1,492))
C     Counter-term amplitude(s) for loop diagram number 263
      CALL R2_GG_1_R2_GG_3_0(W(8),W(25),R2_GGQ,R2_GGT,AMPL(1,493))
C     Counter-term amplitude(s) for loop diagram number 264
      CALL R2_GG_1_R2_GG_3_0(W(8),W(26),R2_GGQ,R2_GGT,AMPL(1,494))
C     Counter-term amplitude(s) for loop diagram number 265
      CALL R2_GGZ_0(W(8),W(20),W(11),R2_GGZUP,AMPL(1,495))
C     Counter-term amplitude(s) for loop diagram number 267
      CALL R2_GG_1_R2_GG_3_0(W(14),W(28),R2_GGQ,R2_GGT,AMPL(1,496))
C     Counter-term amplitude(s) for loop diagram number 268
      CALL R2_GG_1_R2_GG_3_0(W(14),W(30),R2_GGQ,R2_GGT,AMPL(1,497))
C     Counter-term amplitude(s) for loop diagram number 269
      CALL R2_GG_1_R2_GG_3_0(W(14),W(31),R2_GGQ,R2_GGT,AMPL(1,498))
C     Counter-term amplitude(s) for loop diagram number 270
      CALL R2_GG_1_R2_GG_3_0(W(14),W(32),R2_GGQ,R2_GGT,AMPL(1,499))
C     Counter-term amplitude(s) for loop diagram number 271
      CALL R2_GGZ_0(W(17),W(14),W(11),R2_GGZUP,AMPL(1,500))
C     Counter-term amplitude(s) for loop diagram number 273
      CALL R2_GG_1_R2_GG_3_0(W(17),W(35),R2_GGQ,R2_GGT,AMPL(1,501))
C     Counter-term amplitude(s) for loop diagram number 274
      CALL R2_GG_1_R2_GG_3_0(W(17),W(37),R2_GGQ,R2_GGT,AMPL(1,502))
C     Counter-term amplitude(s) for loop diagram number 275
      CALL R2_GG_1_R2_GG_3_0(W(17),W(38),R2_GGQ,R2_GGT,AMPL(1,503))
C     Counter-term amplitude(s) for loop diagram number 276
      CALL R2_GG_1_R2_GG_3_0(W(17),W(39),R2_GGQ,R2_GGT,AMPL(1,504))
C     Counter-term amplitude(s) for loop diagram number 277
      CALL R2_GG_1_R2_GG_3_0(W(20),W(40),R2_GGQ,R2_GGT,AMPL(1,505))
C     Counter-term amplitude(s) for loop diagram number 278
      CALL R2_GG_1_R2_GG_3_0(W(20),W(42),R2_GGQ,R2_GGT,AMPL(1,506))
C     Counter-term amplitude(s) for loop diagram number 279
      CALL R2_GG_1_R2_GG_3_0(W(20),W(43),R2_GGQ,R2_GGT,AMPL(1,507))
C     Counter-term amplitude(s) for loop diagram number 280
      CALL R2_GG_1_R2_GG_3_0(W(20),W(44),R2_GGQ,R2_GGT,AMPL(1,508))
C     At this point, all CT amps needed for (QCD=6 QED=4), i.e. of
C      split order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 2000

      GOTO 1001
 2000 CONTINUE
      CT_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

