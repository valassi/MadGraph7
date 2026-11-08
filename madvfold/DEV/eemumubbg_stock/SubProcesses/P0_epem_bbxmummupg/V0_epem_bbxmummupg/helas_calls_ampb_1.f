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
      PARAMETER (NEXTERNAL=7)
      INTEGER    NCOMB
      PARAMETER (NCOMB=128)
      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=58)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=655, NLOOPGROUPS=125, NCTAMPS=1150)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=1805)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=126,NLOOPWAVEFUNCS=896)
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
      CALL IXXXXX(P(0,4),MDL_MB,NHEL(4),-1, FLAVOR(4),W(4))
      CALL OXXXXX(P(0,5),ZERO,NHEL(5),+1, FLAVOR(5),W(5))
      CALL IXXXXX(P(0,6),ZERO,NHEL(6),-1, FLAVOR(6),W(6))
      CALL VXXXXX(P(0,7),ZERO,NHEL(7),+1,W(7))
      CALL FFV1P0_3(W(2),W(1),GC_3,ZERO,ZERO,W(8))
      CALL FFV1_1(W(3),W(7),GC_5,MDL_MB,ZERO,W(9))
      CALL FFV1_1(W(5),W(8),GC_3,ZERO,ZERO,W(10))
      CALL FFV1P0_3(W(4),W(9),GC_1,ZERO,ZERO,W(11))
C     Amplitude(s) for born diagram with ID 1
      CALL FFV1_0(W(6),W(10),W(11),GC_3,AMP(1))
      CALL FFV2_3_3(W(4),W(9),GC_21,GC_23,MDL_MZ,MDL_WZ,W(12))
C     Amplitude(s) for born diagram with ID 2
      CALL FFV2_4_0(W(6),W(10),W(12),GC_21,GC_24,AMP(2))
      CALL FFV1_2(W(6),W(8),GC_3,ZERO,ZERO,W(13))
C     Amplitude(s) for born diagram with ID 3
      CALL FFV1_0(W(13),W(5),W(11),GC_3,AMP(3))
C     Amplitude(s) for born diagram with ID 4
      CALL FFV2_4_0(W(13),W(5),W(12),GC_21,GC_24,AMP(4))
      CALL FFV2_4_3(W(2),W(1),GC_21,GC_24,MDL_MZ,MDL_WZ,W(14))
      CALL FFV2_4_1(W(5),W(14),GC_21,GC_24,ZERO,ZERO,W(15))
C     Amplitude(s) for born diagram with ID 5
      CALL FFV1_0(W(6),W(15),W(11),GC_3,AMP(5))
C     Amplitude(s) for born diagram with ID 6
      CALL FFV2_4_0(W(6),W(15),W(12),GC_21,GC_24,AMP(6))
      CALL FFV2_4_2(W(6),W(14),GC_21,GC_24,ZERO,ZERO,W(16))
C     Amplitude(s) for born diagram with ID 7
      CALL FFV1_0(W(16),W(5),W(11),GC_3,AMP(7))
C     Amplitude(s) for born diagram with ID 8
      CALL FFV2_4_0(W(16),W(5),W(12),GC_21,GC_24,AMP(8))
      CALL FFV1P0_3(W(6),W(5),GC_3,ZERO,ZERO,W(17))
      CALL FFV1_1(W(9),W(8),GC_1,MDL_MB,ZERO,W(18))
C     Amplitude(s) for born diagram with ID 9
      CALL FFV1_0(W(4),W(18),W(17),GC_1,AMP(9))
      CALL FFV1_2(W(4),W(8),GC_1,MDL_MB,ZERO,W(19))
C     Amplitude(s) for born diagram with ID 10
      CALL FFV1_0(W(19),W(9),W(17),GC_1,AMP(10))
      CALL FFV2_4_3(W(6),W(5),GC_21,GC_24,MDL_MZ,MDL_WZ,W(20))
C     Amplitude(s) for born diagram with ID 11
      CALL FFV2_3_0(W(4),W(18),W(20),GC_21,GC_23,AMP(11))
C     Amplitude(s) for born diagram with ID 12
      CALL FFV2_3_0(W(19),W(9),W(20),GC_21,GC_23,AMP(12))
      CALL FFV2_3_1(W(9),W(14),GC_21,GC_23,MDL_MB,ZERO,W(21))
C     Amplitude(s) for born diagram with ID 13
      CALL FFV1_0(W(4),W(21),W(17),GC_1,AMP(13))
      CALL FFV2_3_2(W(4),W(14),GC_21,GC_23,MDL_MB,ZERO,W(22))
C     Amplitude(s) for born diagram with ID 14
      CALL FFV1_0(W(22),W(9),W(17),GC_1,AMP(14))
C     Amplitude(s) for born diagram with ID 15
      CALL FFV2_3_0(W(4),W(21),W(20),GC_21,GC_23,AMP(15))
C     Amplitude(s) for born diagram with ID 16
      CALL FFV2_3_0(W(22),W(9),W(20),GC_21,GC_23,AMP(16))
      CALL VVS1_3(W(14),W(20),GC_32,MDL_MH,MDL_WH,W(23))
C     Amplitude(s) for born diagram with ID 17
      CALL FFS1_0(W(4),W(9),W(23),GC_33,AMP(17))
      CALL FFV1_2(W(4),W(7),GC_5,MDL_MB,ZERO,W(24))
      CALL FFV1P0_3(W(24),W(3),GC_1,ZERO,ZERO,W(25))
C     Amplitude(s) for born diagram with ID 18
      CALL FFV1_0(W(6),W(10),W(25),GC_3,AMP(18))
      CALL FFV2_3_3(W(24),W(3),GC_21,GC_23,MDL_MZ,MDL_WZ,W(26))
C     Amplitude(s) for born diagram with ID 19
      CALL FFV2_4_0(W(6),W(10),W(26),GC_21,GC_24,AMP(19))
C     Amplitude(s) for born diagram with ID 20
      CALL FFV1_0(W(13),W(5),W(25),GC_3,AMP(20))
C     Amplitude(s) for born diagram with ID 21
      CALL FFV2_4_0(W(13),W(5),W(26),GC_21,GC_24,AMP(21))
C     Amplitude(s) for born diagram with ID 22
      CALL FFV1_0(W(6),W(15),W(25),GC_3,AMP(22))
C     Amplitude(s) for born diagram with ID 23
      CALL FFV2_4_0(W(6),W(15),W(26),GC_21,GC_24,AMP(23))
C     Amplitude(s) for born diagram with ID 24
      CALL FFV1_0(W(16),W(5),W(25),GC_3,AMP(24))
C     Amplitude(s) for born diagram with ID 25
      CALL FFV2_4_0(W(16),W(5),W(26),GC_21,GC_24,AMP(25))
      CALL FFV1_1(W(3),W(8),GC_1,MDL_MB,ZERO,W(27))
C     Amplitude(s) for born diagram with ID 26
      CALL FFV1_0(W(24),W(27),W(17),GC_1,AMP(26))
      CALL FFV1_2(W(24),W(8),GC_1,MDL_MB,ZERO,W(28))
C     Amplitude(s) for born diagram with ID 27
      CALL FFV1_0(W(28),W(3),W(17),GC_1,AMP(27))
C     Amplitude(s) for born diagram with ID 28
      CALL FFV2_3_0(W(24),W(27),W(20),GC_21,GC_23,AMP(28))
C     Amplitude(s) for born diagram with ID 29
      CALL FFV2_3_0(W(28),W(3),W(20),GC_21,GC_23,AMP(29))
      CALL FFV2_3_1(W(3),W(14),GC_21,GC_23,MDL_MB,ZERO,W(29))
C     Amplitude(s) for born diagram with ID 30
      CALL FFV1_0(W(24),W(29),W(17),GC_1,AMP(30))
      CALL FFV2_3_2(W(24),W(14),GC_21,GC_23,MDL_MB,ZERO,W(30))
C     Amplitude(s) for born diagram with ID 31
      CALL FFV1_0(W(30),W(3),W(17),GC_1,AMP(31))
C     Amplitude(s) for born diagram with ID 32
      CALL FFV2_3_0(W(24),W(29),W(20),GC_21,GC_23,AMP(32))
C     Amplitude(s) for born diagram with ID 33
      CALL FFV2_3_0(W(30),W(3),W(20),GC_21,GC_23,AMP(33))
C     Amplitude(s) for born diagram with ID 34
      CALL FFS1_0(W(24),W(3),W(23),GC_33,AMP(34))
      CALL FFV1_2(W(4),W(17),GC_1,MDL_MB,ZERO,W(31))
C     Amplitude(s) for born diagram with ID 35
      CALL FFV1_0(W(31),W(27),W(7),GC_5,AMP(35))
      CALL FFV1_1(W(3),W(17),GC_1,MDL_MB,ZERO,W(32))
C     Amplitude(s) for born diagram with ID 36
      CALL FFV1_0(W(19),W(32),W(7),GC_5,AMP(36))
      CALL FFV2_3_2(W(4),W(20),GC_21,GC_23,MDL_MB,ZERO,W(33))
C     Amplitude(s) for born diagram with ID 37
      CALL FFV1_0(W(33),W(27),W(7),GC_5,AMP(37))
      CALL FFV2_3_1(W(3),W(20),GC_21,GC_23,MDL_MB,ZERO,W(34))
C     Amplitude(s) for born diagram with ID 38
      CALL FFV1_0(W(19),W(34),W(7),GC_5,AMP(38))
C     Amplitude(s) for born diagram with ID 39
      CALL FFV1_0(W(31),W(29),W(7),GC_5,AMP(39))
C     Amplitude(s) for born diagram with ID 40
      CALL FFV1_0(W(22),W(32),W(7),GC_5,AMP(40))
C     Amplitude(s) for born diagram with ID 41
      CALL FFV1_0(W(33),W(29),W(7),GC_5,AMP(41))
C     Amplitude(s) for born diagram with ID 42
      CALL FFV1_0(W(22),W(34),W(7),GC_5,AMP(42))
      CALL FFV1_1(W(1),W(17),GC_3,ZERO,ZERO,W(35))
C     Amplitude(s) for born diagram with ID 43
      CALL FFV1_0(W(2),W(35),W(11),GC_3,AMP(43))
C     Amplitude(s) for born diagram with ID 44
      CALL FFV2_4_0(W(2),W(35),W(12),GC_21,GC_24,AMP(44))
      CALL FFV1_2(W(2),W(17),GC_3,ZERO,ZERO,W(36))
C     Amplitude(s) for born diagram with ID 45
      CALL FFV1_0(W(36),W(1),W(11),GC_3,AMP(45))
C     Amplitude(s) for born diagram with ID 46
      CALL FFV2_4_0(W(36),W(1),W(12),GC_21,GC_24,AMP(46))
      CALL FFV2_4_1(W(1),W(20),GC_21,GC_24,ZERO,ZERO,W(37))
C     Amplitude(s) for born diagram with ID 47
      CALL FFV1_0(W(2),W(37),W(11),GC_3,AMP(47))
C     Amplitude(s) for born diagram with ID 48
      CALL FFV2_4_0(W(2),W(37),W(12),GC_21,GC_24,AMP(48))
      CALL FFV2_4_2(W(2),W(20),GC_21,GC_24,ZERO,ZERO,W(38))
C     Amplitude(s) for born diagram with ID 49
      CALL FFV1_0(W(38),W(1),W(11),GC_3,AMP(49))
C     Amplitude(s) for born diagram with ID 50
      CALL FFV2_4_0(W(38),W(1),W(12),GC_21,GC_24,AMP(50))
C     Amplitude(s) for born diagram with ID 51
      CALL FFV1_0(W(2),W(35),W(25),GC_3,AMP(51))
C     Amplitude(s) for born diagram with ID 52
      CALL FFV2_4_0(W(2),W(35),W(26),GC_21,GC_24,AMP(52))
C     Amplitude(s) for born diagram with ID 53
      CALL FFV1_0(W(36),W(1),W(25),GC_3,AMP(53))
C     Amplitude(s) for born diagram with ID 54
      CALL FFV2_4_0(W(36),W(1),W(26),GC_21,GC_24,AMP(54))
C     Amplitude(s) for born diagram with ID 55
      CALL FFV1_0(W(2),W(37),W(25),GC_3,AMP(55))
C     Amplitude(s) for born diagram with ID 56
      CALL FFV2_4_0(W(2),W(37),W(26),GC_21,GC_24,AMP(56))
C     Amplitude(s) for born diagram with ID 57
      CALL FFV1_0(W(38),W(1),W(25),GC_3,AMP(57))
C     Amplitude(s) for born diagram with ID 58
      CALL FFV2_4_0(W(38),W(1),W(26),GC_21,GC_24,AMP(58))
      CALL FFV1P0_3(W(6),W(10),GC_3,ZERO,ZERO,W(39))
      CALL FFV1_2(W(4),W(39),GC_1,MDL_MB,ZERO,W(40))
C     Counter-term amplitude(s) for loop diagram number 59
      CALL R2_QQ_1_R2_QQ_2_0(W(40),W(9),R2_QQQ,R2_QQB,AMPL(1,1))
      CALL R2_QQ_2_0(W(40),W(9),UV_BMASS_1EPS,AMPL(2,2))
      CALL R2_QQ_2_0(W(40),W(9),UV_BMASS,AMPL(1,3))
      CALL FFV2_4_3(W(6),W(10),GC_21,GC_24,MDL_MZ,MDL_WZ,W(41))
      CALL FFV2_3_2(W(4),W(41),GC_21,GC_23,MDL_MB,ZERO,W(42))
C     Counter-term amplitude(s) for loop diagram number 60
      CALL R2_QQ_1_R2_QQ_2_0(W(42),W(9),R2_QQQ,R2_QQB,AMPL(1,4))
      CALL R2_QQ_2_0(W(42),W(9),UV_BMASS_1EPS,AMPL(2,5))
      CALL R2_QQ_2_0(W(42),W(9),UV_BMASS,AMPL(1,6))
      CALL FFV1P0_3(W(13),W(5),GC_3,ZERO,ZERO,W(43))
      CALL FFV1_2(W(4),W(43),GC_1,MDL_MB,ZERO,W(44))
C     Counter-term amplitude(s) for loop diagram number 61
      CALL R2_QQ_1_R2_QQ_2_0(W(44),W(9),R2_QQQ,R2_QQB,AMPL(1,7))
      CALL R2_QQ_2_0(W(44),W(9),UV_BMASS_1EPS,AMPL(2,8))
      CALL R2_QQ_2_0(W(44),W(9),UV_BMASS,AMPL(1,9))
      CALL FFV2_4_3(W(13),W(5),GC_21,GC_24,MDL_MZ,MDL_WZ,W(45))
      CALL FFV2_3_2(W(4),W(45),GC_21,GC_23,MDL_MB,ZERO,W(46))
C     Counter-term amplitude(s) for loop diagram number 62
      CALL R2_QQ_1_R2_QQ_2_0(W(46),W(9),R2_QQQ,R2_QQB,AMPL(1,10))
      CALL R2_QQ_2_0(W(46),W(9),UV_BMASS_1EPS,AMPL(2,11))
      CALL R2_QQ_2_0(W(46),W(9),UV_BMASS,AMPL(1,12))
      CALL FFV1P0_3(W(6),W(15),GC_3,ZERO,ZERO,W(47))
      CALL FFV1_2(W(4),W(47),GC_1,MDL_MB,ZERO,W(48))
C     Counter-term amplitude(s) for loop diagram number 63
      CALL R2_QQ_1_R2_QQ_2_0(W(48),W(9),R2_QQQ,R2_QQB,AMPL(1,13))
      CALL R2_QQ_2_0(W(48),W(9),UV_BMASS_1EPS,AMPL(2,14))
      CALL R2_QQ_2_0(W(48),W(9),UV_BMASS,AMPL(1,15))
      CALL FFV2_4_3(W(6),W(15),GC_21,GC_24,MDL_MZ,MDL_WZ,W(49))
      CALL FFV2_3_2(W(4),W(49),GC_21,GC_23,MDL_MB,ZERO,W(50))
C     Counter-term amplitude(s) for loop diagram number 64
      CALL R2_QQ_1_R2_QQ_2_0(W(50),W(9),R2_QQQ,R2_QQB,AMPL(1,16))
      CALL R2_QQ_2_0(W(50),W(9),UV_BMASS_1EPS,AMPL(2,17))
      CALL R2_QQ_2_0(W(50),W(9),UV_BMASS,AMPL(1,18))
      CALL FFV1P0_3(W(16),W(5),GC_3,ZERO,ZERO,W(51))
      CALL FFV1_2(W(4),W(51),GC_1,MDL_MB,ZERO,W(52))
C     Counter-term amplitude(s) for loop diagram number 65
      CALL R2_QQ_1_R2_QQ_2_0(W(52),W(9),R2_QQQ,R2_QQB,AMPL(1,19))
      CALL R2_QQ_2_0(W(52),W(9),UV_BMASS_1EPS,AMPL(2,20))
      CALL R2_QQ_2_0(W(52),W(9),UV_BMASS,AMPL(1,21))
      CALL FFV2_4_3(W(16),W(5),GC_21,GC_24,MDL_MZ,MDL_WZ,W(53))
      CALL FFV2_3_2(W(4),W(53),GC_21,GC_23,MDL_MB,ZERO,W(54))
C     Counter-term amplitude(s) for loop diagram number 66
      CALL R2_QQ_1_R2_QQ_2_0(W(54),W(9),R2_QQQ,R2_QQB,AMPL(1,22))
      CALL R2_QQ_2_0(W(54),W(9),UV_BMASS_1EPS,AMPL(2,23))
      CALL R2_QQ_2_0(W(54),W(9),UV_BMASS,AMPL(1,24))
C     Counter-term amplitude(s) for loop diagram number 67
      CALL FFV1_0(W(4),W(9),W(39),R2_DDA,AMPL(1,25))
C     Counter-term amplitude(s) for loop diagram number 68
      CALL FFV2_3_0(W(4),W(9),W(41),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,26))
C     Counter-term amplitude(s) for loop diagram number 69
      CALL FFV1_0(W(4),W(9),W(43),R2_DDA,AMPL(1,27))
C     Counter-term amplitude(s) for loop diagram number 70
      CALL FFV2_3_0(W(4),W(9),W(45),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,28))
C     Counter-term amplitude(s) for loop diagram number 71
      CALL FFV1_0(W(4),W(9),W(47),R2_DDA,AMPL(1,29))
C     Counter-term amplitude(s) for loop diagram number 72
      CALL FFV2_3_0(W(4),W(9),W(49),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,30))
C     Counter-term amplitude(s) for loop diagram number 73
      CALL FFV1_0(W(4),W(9),W(51),R2_DDA,AMPL(1,31))
C     Counter-term amplitude(s) for loop diagram number 74
      CALL FFV2_3_0(W(4),W(9),W(53),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,32))
C     Counter-term amplitude(s) for loop diagram number 75
      CALL FFV1_0(W(4),W(18),W(17),R2_DDA,AMPL(1,33))
      CALL FFV1_1(W(9),W(17),GC_1,MDL_MB,ZERO,W(55))
C     Counter-term amplitude(s) for loop diagram number 76
      CALL FFV1_0(W(4),W(55),W(8),R2_DDA,AMPL(1,34))
C     Counter-term amplitude(s) for loop diagram number 79
      CALL FFV2_3_0(W(4),W(18),W(20),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,35))
      CALL FFV2_3_1(W(9),W(20),GC_21,GC_23,MDL_MB,ZERO,W(56))
C     Counter-term amplitude(s) for loop diagram number 80
      CALL FFV1_0(W(4),W(56),W(8),R2_DDA,AMPL(1,36))
C     Counter-term amplitude(s) for loop diagram number 83
      CALL FFV1_0(W(4),W(21),W(17),R2_DDA,AMPL(1,37))
C     Counter-term amplitude(s) for loop diagram number 84
      CALL FFV2_3_0(W(4),W(55),W(14),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,38))
C     Counter-term amplitude(s) for loop diagram number 87
      CALL FFV2_3_0(W(4),W(21),W(20),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,39))
C     Counter-term amplitude(s) for loop diagram number 88
      CALL FFV2_3_0(W(4),W(56),W(14),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,40))
      CALL VVS1_3(W(20),W(14),GC_32,MDL_MH,MDL_WH,W(57))
C     Counter-term amplitude(s) for loop diagram number 90
      CALL FFS1_0(W(4),W(9),W(57),R2_BBH,AMPL(1,41))
      CALL FFS1_0(W(4),W(9),W(57),UV_HBB_1EPS,AMPL(2,42))
      CALL FFS1_0(W(4),W(9),W(57),UV_HBB,AMPL(1,43))
C     Counter-term amplitude(s) for loop diagram number 92
      CALL R2_QQ_1_R2_QQ_2_0(W(31),W(18),R2_QQQ,R2_QQB,AMPL(1,44))
      CALL R2_QQ_2_0(W(31),W(18),UV_BMASS_1EPS,AMPL(2,45))
      CALL R2_QQ_2_0(W(31),W(18),UV_BMASS,AMPL(1,46))
C     Counter-term amplitude(s) for loop diagram number 93
      CALL R2_QQ_1_R2_QQ_2_0(W(19),W(55),R2_QQQ,R2_QQB,AMPL(1,47))
      CALL R2_QQ_2_0(W(19),W(55),UV_BMASS_1EPS,AMPL(2,48))
      CALL R2_QQ_2_0(W(19),W(55),UV_BMASS,AMPL(1,49))
      CALL FFV1_2(W(19),W(17),GC_1,MDL_MB,ZERO,W(58))
C     Counter-term amplitude(s) for loop diagram number 94
      CALL R2_QQ_1_R2_QQ_2_0(W(58),W(9),R2_QQQ,R2_QQB,AMPL(1,50))
      CALL R2_QQ_2_0(W(58),W(9),UV_BMASS_1EPS,AMPL(2,51))
      CALL R2_QQ_2_0(W(58),W(9),UV_BMASS,AMPL(1,52))
C     Counter-term amplitude(s) for loop diagram number 95
      CALL FFV1_0(W(19),W(9),W(17),R2_DDA,AMPL(1,53))
C     Counter-term amplitude(s) for loop diagram number 96
      CALL FFV1_0(W(31),W(9),W(8),R2_DDA,AMPL(1,54))
      CALL FFV1_2(W(31),W(8),GC_1,MDL_MB,ZERO,W(59))
C     Counter-term amplitude(s) for loop diagram number 97
      CALL R2_QQ_1_R2_QQ_2_0(W(59),W(9),R2_QQQ,R2_QQB,AMPL(1,55))
      CALL R2_QQ_2_0(W(59),W(9),UV_BMASS_1EPS,AMPL(2,56))
      CALL R2_QQ_2_0(W(59),W(9),UV_BMASS,AMPL(1,57))
C     Counter-term amplitude(s) for loop diagram number 98
      CALL R2_QQ_1_R2_QQ_2_0(W(33),W(18),R2_QQQ,R2_QQB,AMPL(1,58))
      CALL R2_QQ_2_0(W(33),W(18),UV_BMASS_1EPS,AMPL(2,59))
      CALL R2_QQ_2_0(W(33),W(18),UV_BMASS,AMPL(1,60))
C     Counter-term amplitude(s) for loop diagram number 99
      CALL R2_QQ_1_R2_QQ_2_0(W(19),W(56),R2_QQQ,R2_QQB,AMPL(1,61))
      CALL R2_QQ_2_0(W(19),W(56),UV_BMASS_1EPS,AMPL(2,62))
      CALL R2_QQ_2_0(W(19),W(56),UV_BMASS,AMPL(1,63))
      CALL FFV2_3_2(W(19),W(20),GC_21,GC_23,MDL_MB,ZERO,W(60))
C     Counter-term amplitude(s) for loop diagram number 100
      CALL R2_QQ_1_R2_QQ_2_0(W(60),W(9),R2_QQQ,R2_QQB,AMPL(1,64))
      CALL R2_QQ_2_0(W(60),W(9),UV_BMASS_1EPS,AMPL(2,65))
      CALL R2_QQ_2_0(W(60),W(9),UV_BMASS,AMPL(1,66))
C     Counter-term amplitude(s) for loop diagram number 101
      CALL FFV2_3_0(W(19),W(9),W(20),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,67))
C     Counter-term amplitude(s) for loop diagram number 102
      CALL FFV1_0(W(33),W(9),W(8),R2_DDA,AMPL(1,68))
      CALL FFV1_2(W(33),W(8),GC_1,MDL_MB,ZERO,W(61))
C     Counter-term amplitude(s) for loop diagram number 103
      CALL R2_QQ_1_R2_QQ_2_0(W(61),W(9),R2_QQQ,R2_QQB,AMPL(1,69))
      CALL R2_QQ_2_0(W(61),W(9),UV_BMASS_1EPS,AMPL(2,70))
      CALL R2_QQ_2_0(W(61),W(9),UV_BMASS,AMPL(1,71))
C     Counter-term amplitude(s) for loop diagram number 104
      CALL R2_QQ_1_R2_QQ_2_0(W(31),W(21),R2_QQQ,R2_QQB,AMPL(1,72))
      CALL R2_QQ_2_0(W(31),W(21),UV_BMASS_1EPS,AMPL(2,73))
      CALL R2_QQ_2_0(W(31),W(21),UV_BMASS,AMPL(1,74))
C     Counter-term amplitude(s) for loop diagram number 105
      CALL R2_QQ_1_R2_QQ_2_0(W(22),W(55),R2_QQQ,R2_QQB,AMPL(1,75))
      CALL R2_QQ_2_0(W(22),W(55),UV_BMASS_1EPS,AMPL(2,76))
      CALL R2_QQ_2_0(W(22),W(55),UV_BMASS,AMPL(1,77))
      CALL FFV1_2(W(22),W(17),GC_1,MDL_MB,ZERO,W(62))
C     Counter-term amplitude(s) for loop diagram number 106
      CALL R2_QQ_1_R2_QQ_2_0(W(62),W(9),R2_QQQ,R2_QQB,AMPL(1,78))
      CALL R2_QQ_2_0(W(62),W(9),UV_BMASS_1EPS,AMPL(2,79))
      CALL R2_QQ_2_0(W(62),W(9),UV_BMASS,AMPL(1,80))
C     Counter-term amplitude(s) for loop diagram number 107
      CALL FFV1_0(W(22),W(9),W(17),R2_DDA,AMPL(1,81))
C     Counter-term amplitude(s) for loop diagram number 108
      CALL FFV2_3_0(W(31),W(9),W(14),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,82))
      CALL FFV2_3_2(W(31),W(14),GC_21,GC_23,MDL_MB,ZERO,W(63))
C     Counter-term amplitude(s) for loop diagram number 109
      CALL R2_QQ_1_R2_QQ_2_0(W(63),W(9),R2_QQQ,R2_QQB,AMPL(1,83))
      CALL R2_QQ_2_0(W(63),W(9),UV_BMASS_1EPS,AMPL(2,84))
      CALL R2_QQ_2_0(W(63),W(9),UV_BMASS,AMPL(1,85))
C     Counter-term amplitude(s) for loop diagram number 110
      CALL R2_QQ_1_R2_QQ_2_0(W(33),W(21),R2_QQQ,R2_QQB,AMPL(1,86))
      CALL R2_QQ_2_0(W(33),W(21),UV_BMASS_1EPS,AMPL(2,87))
      CALL R2_QQ_2_0(W(33),W(21),UV_BMASS,AMPL(1,88))
C     Counter-term amplitude(s) for loop diagram number 111
      CALL R2_QQ_1_R2_QQ_2_0(W(22),W(56),R2_QQQ,R2_QQB,AMPL(1,89))
      CALL R2_QQ_2_0(W(22),W(56),UV_BMASS_1EPS,AMPL(2,90))
      CALL R2_QQ_2_0(W(22),W(56),UV_BMASS,AMPL(1,91))
      CALL FFV2_3_2(W(22),W(20),GC_21,GC_23,MDL_MB,ZERO,W(64))
C     Counter-term amplitude(s) for loop diagram number 112
      CALL R2_QQ_1_R2_QQ_2_0(W(64),W(9),R2_QQQ,R2_QQB,AMPL(1,92))
      CALL R2_QQ_2_0(W(64),W(9),UV_BMASS_1EPS,AMPL(2,93))
      CALL R2_QQ_2_0(W(64),W(9),UV_BMASS,AMPL(1,94))
C     Counter-term amplitude(s) for loop diagram number 113
      CALL FFV2_3_0(W(22),W(9),W(20),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,95))
      CALL FFS1_2(W(4),W(57),GC_33,MDL_MB,ZERO,W(65))
C     Counter-term amplitude(s) for loop diagram number 114
      CALL R2_QQ_1_R2_QQ_2_0(W(65),W(9),R2_QQQ,R2_QQB,AMPL(1,96))
      CALL R2_QQ_2_0(W(65),W(9),UV_BMASS_1EPS,AMPL(2,97))
      CALL R2_QQ_2_0(W(65),W(9),UV_BMASS,AMPL(1,98))
C     Counter-term amplitude(s) for loop diagram number 115
      CALL FFV2_3_0(W(33),W(9),W(14),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,99))
      CALL FFV2_3_2(W(33),W(14),GC_21,GC_23,MDL_MB,ZERO,W(66))
C     Counter-term amplitude(s) for loop diagram number 116
      CALL R2_QQ_1_R2_QQ_2_0(W(66),W(9),R2_QQQ,R2_QQB,AMPL(1,100))
      CALL R2_QQ_2_0(W(66),W(9),UV_BMASS_1EPS,AMPL(2,101))
      CALL R2_QQ_2_0(W(66),W(9),UV_BMASS,AMPL(1,102))
C     Counter-term amplitude(s) for loop diagram number 117
      CALL FFV1_0(W(40),W(3),W(7),R2_GQQ,AMPL(1,103))
      CALL FFV1_0(W(40),W(3),W(7),UV_GQQB_1EPS,AMPL(2,104))
      CALL FFV1_0(W(40),W(3),W(7),UV_GQQB_1EPS,AMPL(2,105))
      CALL FFV1_0(W(40),W(3),W(7),UV_GQQB_1EPS,AMPL(2,106))
      CALL FFV1_0(W(40),W(3),W(7),UV_GQQB_1EPS,AMPL(2,107))
      CALL FFV1_0(W(40),W(3),W(7),UV_GQQB_1EPS,AMPL(2,108))
      CALL FFV1_0(W(40),W(3),W(7),UV_GQQB_1EPS,AMPL(2,109))
      CALL FFV1_0(W(40),W(3),W(7),UV_GQQG_1EPS,AMPL(2,110))
      CALL FFV1_0(W(40),W(3),W(7),UV_GQQB,AMPL(1,111))
      CALL FFV1_0(W(40),W(3),W(7),UV_GQQT,AMPL(1,112))
C     Counter-term amplitude(s) for loop diagram number 118
      CALL FFV1_0(W(42),W(3),W(7),R2_GQQ,AMPL(1,113))
      CALL FFV1_0(W(42),W(3),W(7),UV_GQQB_1EPS,AMPL(2,114))
      CALL FFV1_0(W(42),W(3),W(7),UV_GQQB_1EPS,AMPL(2,115))
      CALL FFV1_0(W(42),W(3),W(7),UV_GQQB_1EPS,AMPL(2,116))
      CALL FFV1_0(W(42),W(3),W(7),UV_GQQB_1EPS,AMPL(2,117))
      CALL FFV1_0(W(42),W(3),W(7),UV_GQQB_1EPS,AMPL(2,118))
      CALL FFV1_0(W(42),W(3),W(7),UV_GQQB_1EPS,AMPL(2,119))
      CALL FFV1_0(W(42),W(3),W(7),UV_GQQG_1EPS,AMPL(2,120))
      CALL FFV1_0(W(42),W(3),W(7),UV_GQQB,AMPL(1,121))
      CALL FFV1_0(W(42),W(3),W(7),UV_GQQT,AMPL(1,122))
C     Counter-term amplitude(s) for loop diagram number 119
      CALL FFV1_0(W(44),W(3),W(7),R2_GQQ,AMPL(1,123))
      CALL FFV1_0(W(44),W(3),W(7),UV_GQQB_1EPS,AMPL(2,124))
      CALL FFV1_0(W(44),W(3),W(7),UV_GQQB_1EPS,AMPL(2,125))
      CALL FFV1_0(W(44),W(3),W(7),UV_GQQB_1EPS,AMPL(2,126))
      CALL FFV1_0(W(44),W(3),W(7),UV_GQQB_1EPS,AMPL(2,127))
      CALL FFV1_0(W(44),W(3),W(7),UV_GQQB_1EPS,AMPL(2,128))
      CALL FFV1_0(W(44),W(3),W(7),UV_GQQB_1EPS,AMPL(2,129))
      CALL FFV1_0(W(44),W(3),W(7),UV_GQQG_1EPS,AMPL(2,130))
      CALL FFV1_0(W(44),W(3),W(7),UV_GQQB,AMPL(1,131))
      CALL FFV1_0(W(44),W(3),W(7),UV_GQQT,AMPL(1,132))
C     Counter-term amplitude(s) for loop diagram number 120
      CALL FFV1_0(W(46),W(3),W(7),R2_GQQ,AMPL(1,133))
      CALL FFV1_0(W(46),W(3),W(7),UV_GQQB_1EPS,AMPL(2,134))
      CALL FFV1_0(W(46),W(3),W(7),UV_GQQB_1EPS,AMPL(2,135))
      CALL FFV1_0(W(46),W(3),W(7),UV_GQQB_1EPS,AMPL(2,136))
      CALL FFV1_0(W(46),W(3),W(7),UV_GQQB_1EPS,AMPL(2,137))
      CALL FFV1_0(W(46),W(3),W(7),UV_GQQB_1EPS,AMPL(2,138))
      CALL FFV1_0(W(46),W(3),W(7),UV_GQQB_1EPS,AMPL(2,139))
      CALL FFV1_0(W(46),W(3),W(7),UV_GQQG_1EPS,AMPL(2,140))
      CALL FFV1_0(W(46),W(3),W(7),UV_GQQB,AMPL(1,141))
      CALL FFV1_0(W(46),W(3),W(7),UV_GQQT,AMPL(1,142))
C     Counter-term amplitude(s) for loop diagram number 121
      CALL FFV1_0(W(48),W(3),W(7),R2_GQQ,AMPL(1,143))
      CALL FFV1_0(W(48),W(3),W(7),UV_GQQB_1EPS,AMPL(2,144))
      CALL FFV1_0(W(48),W(3),W(7),UV_GQQB_1EPS,AMPL(2,145))
      CALL FFV1_0(W(48),W(3),W(7),UV_GQQB_1EPS,AMPL(2,146))
      CALL FFV1_0(W(48),W(3),W(7),UV_GQQB_1EPS,AMPL(2,147))
      CALL FFV1_0(W(48),W(3),W(7),UV_GQQB_1EPS,AMPL(2,148))
      CALL FFV1_0(W(48),W(3),W(7),UV_GQQB_1EPS,AMPL(2,149))
      CALL FFV1_0(W(48),W(3),W(7),UV_GQQG_1EPS,AMPL(2,150))
      CALL FFV1_0(W(48),W(3),W(7),UV_GQQB,AMPL(1,151))
      CALL FFV1_0(W(48),W(3),W(7),UV_GQQT,AMPL(1,152))
C     Counter-term amplitude(s) for loop diagram number 122
      CALL FFV1_0(W(50),W(3),W(7),R2_GQQ,AMPL(1,153))
      CALL FFV1_0(W(50),W(3),W(7),UV_GQQB_1EPS,AMPL(2,154))
      CALL FFV1_0(W(50),W(3),W(7),UV_GQQB_1EPS,AMPL(2,155))
      CALL FFV1_0(W(50),W(3),W(7),UV_GQQB_1EPS,AMPL(2,156))
      CALL FFV1_0(W(50),W(3),W(7),UV_GQQB_1EPS,AMPL(2,157))
      CALL FFV1_0(W(50),W(3),W(7),UV_GQQB_1EPS,AMPL(2,158))
      CALL FFV1_0(W(50),W(3),W(7),UV_GQQB_1EPS,AMPL(2,159))
      CALL FFV1_0(W(50),W(3),W(7),UV_GQQG_1EPS,AMPL(2,160))
      CALL FFV1_0(W(50),W(3),W(7),UV_GQQB,AMPL(1,161))
      CALL FFV1_0(W(50),W(3),W(7),UV_GQQT,AMPL(1,162))
C     Counter-term amplitude(s) for loop diagram number 123
      CALL FFV1_0(W(52),W(3),W(7),R2_GQQ,AMPL(1,163))
      CALL FFV1_0(W(52),W(3),W(7),UV_GQQB_1EPS,AMPL(2,164))
      CALL FFV1_0(W(52),W(3),W(7),UV_GQQB_1EPS,AMPL(2,165))
      CALL FFV1_0(W(52),W(3),W(7),UV_GQQB_1EPS,AMPL(2,166))
      CALL FFV1_0(W(52),W(3),W(7),UV_GQQB_1EPS,AMPL(2,167))
      CALL FFV1_0(W(52),W(3),W(7),UV_GQQB_1EPS,AMPL(2,168))
      CALL FFV1_0(W(52),W(3),W(7),UV_GQQB_1EPS,AMPL(2,169))
      CALL FFV1_0(W(52),W(3),W(7),UV_GQQG_1EPS,AMPL(2,170))
      CALL FFV1_0(W(52),W(3),W(7),UV_GQQB,AMPL(1,171))
      CALL FFV1_0(W(52),W(3),W(7),UV_GQQT,AMPL(1,172))
C     Counter-term amplitude(s) for loop diagram number 124
      CALL FFV1_0(W(54),W(3),W(7),R2_GQQ,AMPL(1,173))
      CALL FFV1_0(W(54),W(3),W(7),UV_GQQB_1EPS,AMPL(2,174))
      CALL FFV1_0(W(54),W(3),W(7),UV_GQQB_1EPS,AMPL(2,175))
      CALL FFV1_0(W(54),W(3),W(7),UV_GQQB_1EPS,AMPL(2,176))
      CALL FFV1_0(W(54),W(3),W(7),UV_GQQB_1EPS,AMPL(2,177))
      CALL FFV1_0(W(54),W(3),W(7),UV_GQQB_1EPS,AMPL(2,178))
      CALL FFV1_0(W(54),W(3),W(7),UV_GQQB_1EPS,AMPL(2,179))
      CALL FFV1_0(W(54),W(3),W(7),UV_GQQG_1EPS,AMPL(2,180))
      CALL FFV1_0(W(54),W(3),W(7),UV_GQQB,AMPL(1,181))
      CALL FFV1_0(W(54),W(3),W(7),UV_GQQT,AMPL(1,182))
C     Counter-term amplitude(s) for loop diagram number 125
      CALL FFV1_0(W(24),W(3),W(39),R2_DDA,AMPL(1,183))
C     Counter-term amplitude(s) for loop diagram number 126
      CALL FFV2_3_0(W(24),W(3),W(41),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,184))
C     Counter-term amplitude(s) for loop diagram number 127
      CALL FFV1_0(W(24),W(3),W(43),R2_DDA,AMPL(1,185))
C     Counter-term amplitude(s) for loop diagram number 128
      CALL FFV2_3_0(W(24),W(3),W(45),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,186))
C     Counter-term amplitude(s) for loop diagram number 129
      CALL FFV1_0(W(24),W(3),W(47),R2_DDA,AMPL(1,187))
C     Counter-term amplitude(s) for loop diagram number 130
      CALL FFV2_3_0(W(24),W(3),W(49),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,188))
C     Counter-term amplitude(s) for loop diagram number 131
      CALL FFV1_0(W(24),W(3),W(51),R2_DDA,AMPL(1,189))
C     Counter-term amplitude(s) for loop diagram number 132
      CALL FFV2_3_0(W(24),W(3),W(53),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,190))
      CALL FFV1_2(W(24),W(17),GC_1,MDL_MB,ZERO,W(67))
C     Counter-term amplitude(s) for loop diagram number 133
      CALL FFV1_0(W(67),W(3),W(8),R2_DDA,AMPL(1,191))
C     Counter-term amplitude(s) for loop diagram number 135
      CALL FFV1_0(W(28),W(3),W(17),R2_DDA,AMPL(1,192))
      CALL FFV2_3_2(W(24),W(20),GC_21,GC_23,MDL_MB,ZERO,W(68))
C     Counter-term amplitude(s) for loop diagram number 137
      CALL FFV1_0(W(68),W(3),W(8),R2_DDA,AMPL(1,193))
C     Counter-term amplitude(s) for loop diagram number 139
      CALL FFV2_3_0(W(28),W(3),W(20),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,194))
C     Counter-term amplitude(s) for loop diagram number 141
      CALL FFV2_3_0(W(67),W(3),W(14),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,195))
C     Counter-term amplitude(s) for loop diagram number 143
      CALL FFV1_0(W(30),W(3),W(17),R2_DDA,AMPL(1,196))
C     Counter-term amplitude(s) for loop diagram number 145
      CALL FFV2_3_0(W(68),W(3),W(14),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,197))
C     Counter-term amplitude(s) for loop diagram number 147
      CALL FFV2_3_0(W(30),W(3),W(20),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,198))
C     Counter-term amplitude(s) for loop diagram number 148
      CALL FFS1_0(W(24),W(3),W(57),R2_BBH,AMPL(1,199))
      CALL FFS1_0(W(24),W(3),W(57),UV_HBB_1EPS,AMPL(2,200))
      CALL FFS1_0(W(24),W(3),W(57),UV_HBB,AMPL(1,201))
      CALL FFV1_2(W(31),W(7),GC_5,MDL_MB,ZERO,W(69))
C     Counter-term amplitude(s) for loop diagram number 193
      CALL FFV1_0(W(69),W(3),W(8),R2_DDA,AMPL(1,202))
      CALL FFV1_2(W(19),W(7),GC_5,MDL_MB,ZERO,W(70))
C     Counter-term amplitude(s) for loop diagram number 194
      CALL FFV1_0(W(70),W(3),W(17),R2_DDA,AMPL(1,203))
C     Counter-term amplitude(s) for loop diagram number 196
      CALL FFV1_0(W(58),W(3),W(7),R2_GQQ,AMPL(1,204))
      CALL FFV1_0(W(58),W(3),W(7),UV_GQQB_1EPS,AMPL(2,205))
      CALL FFV1_0(W(58),W(3),W(7),UV_GQQB_1EPS,AMPL(2,206))
      CALL FFV1_0(W(58),W(3),W(7),UV_GQQB_1EPS,AMPL(2,207))
      CALL FFV1_0(W(58),W(3),W(7),UV_GQQB_1EPS,AMPL(2,208))
      CALL FFV1_0(W(58),W(3),W(7),UV_GQQB_1EPS,AMPL(2,209))
      CALL FFV1_0(W(58),W(3),W(7),UV_GQQB_1EPS,AMPL(2,210))
      CALL FFV1_0(W(58),W(3),W(7),UV_GQQG_1EPS,AMPL(2,211))
      CALL FFV1_0(W(58),W(3),W(7),UV_GQQB,AMPL(1,212))
      CALL FFV1_0(W(58),W(3),W(7),UV_GQQT,AMPL(1,213))
C     Counter-term amplitude(s) for loop diagram number 199
      CALL FFV1_0(W(59),W(3),W(7),R2_GQQ,AMPL(1,214))
      CALL FFV1_0(W(59),W(3),W(7),UV_GQQB_1EPS,AMPL(2,215))
      CALL FFV1_0(W(59),W(3),W(7),UV_GQQB_1EPS,AMPL(2,216))
      CALL FFV1_0(W(59),W(3),W(7),UV_GQQB_1EPS,AMPL(2,217))
      CALL FFV1_0(W(59),W(3),W(7),UV_GQQB_1EPS,AMPL(2,218))
      CALL FFV1_0(W(59),W(3),W(7),UV_GQQB_1EPS,AMPL(2,219))
      CALL FFV1_0(W(59),W(3),W(7),UV_GQQB_1EPS,AMPL(2,220))
      CALL FFV1_0(W(59),W(3),W(7),UV_GQQG_1EPS,AMPL(2,221))
      CALL FFV1_0(W(59),W(3),W(7),UV_GQQB,AMPL(1,222))
      CALL FFV1_0(W(59),W(3),W(7),UV_GQQT,AMPL(1,223))
      CALL FFV1_2(W(33),W(7),GC_5,MDL_MB,ZERO,W(71))
C     Counter-term amplitude(s) for loop diagram number 201
      CALL FFV1_0(W(71),W(3),W(8),R2_DDA,AMPL(1,224))
C     Counter-term amplitude(s) for loop diagram number 202
      CALL FFV2_3_0(W(70),W(3),W(20),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,225))
C     Counter-term amplitude(s) for loop diagram number 204
      CALL FFV1_0(W(60),W(3),W(7),R2_GQQ,AMPL(1,226))
      CALL FFV1_0(W(60),W(3),W(7),UV_GQQB_1EPS,AMPL(2,227))
      CALL FFV1_0(W(60),W(3),W(7),UV_GQQB_1EPS,AMPL(2,228))
      CALL FFV1_0(W(60),W(3),W(7),UV_GQQB_1EPS,AMPL(2,229))
      CALL FFV1_0(W(60),W(3),W(7),UV_GQQB_1EPS,AMPL(2,230))
      CALL FFV1_0(W(60),W(3),W(7),UV_GQQB_1EPS,AMPL(2,231))
      CALL FFV1_0(W(60),W(3),W(7),UV_GQQB_1EPS,AMPL(2,232))
      CALL FFV1_0(W(60),W(3),W(7),UV_GQQG_1EPS,AMPL(2,233))
      CALL FFV1_0(W(60),W(3),W(7),UV_GQQB,AMPL(1,234))
      CALL FFV1_0(W(60),W(3),W(7),UV_GQQT,AMPL(1,235))
C     Counter-term amplitude(s) for loop diagram number 207
      CALL FFV1_0(W(61),W(3),W(7),R2_GQQ,AMPL(1,236))
      CALL FFV1_0(W(61),W(3),W(7),UV_GQQB_1EPS,AMPL(2,237))
      CALL FFV1_0(W(61),W(3),W(7),UV_GQQB_1EPS,AMPL(2,238))
      CALL FFV1_0(W(61),W(3),W(7),UV_GQQB_1EPS,AMPL(2,239))
      CALL FFV1_0(W(61),W(3),W(7),UV_GQQB_1EPS,AMPL(2,240))
      CALL FFV1_0(W(61),W(3),W(7),UV_GQQB_1EPS,AMPL(2,241))
      CALL FFV1_0(W(61),W(3),W(7),UV_GQQB_1EPS,AMPL(2,242))
      CALL FFV1_0(W(61),W(3),W(7),UV_GQQG_1EPS,AMPL(2,243))
      CALL FFV1_0(W(61),W(3),W(7),UV_GQQB,AMPL(1,244))
      CALL FFV1_0(W(61),W(3),W(7),UV_GQQT,AMPL(1,245))
C     Counter-term amplitude(s) for loop diagram number 209
      CALL FFV2_3_0(W(69),W(3),W(14),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,246))
      CALL FFV1_2(W(22),W(7),GC_5,MDL_MB,ZERO,W(72))
C     Counter-term amplitude(s) for loop diagram number 210
      CALL FFV1_0(W(72),W(3),W(17),R2_DDA,AMPL(1,247))
C     Counter-term amplitude(s) for loop diagram number 212
      CALL FFV1_0(W(62),W(3),W(7),R2_GQQ,AMPL(1,248))
      CALL FFV1_0(W(62),W(3),W(7),UV_GQQB_1EPS,AMPL(2,249))
      CALL FFV1_0(W(62),W(3),W(7),UV_GQQB_1EPS,AMPL(2,250))
      CALL FFV1_0(W(62),W(3),W(7),UV_GQQB_1EPS,AMPL(2,251))
      CALL FFV1_0(W(62),W(3),W(7),UV_GQQB_1EPS,AMPL(2,252))
      CALL FFV1_0(W(62),W(3),W(7),UV_GQQB_1EPS,AMPL(2,253))
      CALL FFV1_0(W(62),W(3),W(7),UV_GQQB_1EPS,AMPL(2,254))
      CALL FFV1_0(W(62),W(3),W(7),UV_GQQG_1EPS,AMPL(2,255))
      CALL FFV1_0(W(62),W(3),W(7),UV_GQQB,AMPL(1,256))
      CALL FFV1_0(W(62),W(3),W(7),UV_GQQT,AMPL(1,257))
C     Counter-term amplitude(s) for loop diagram number 215
      CALL FFV1_0(W(63),W(3),W(7),R2_GQQ,AMPL(1,258))
      CALL FFV1_0(W(63),W(3),W(7),UV_GQQB_1EPS,AMPL(2,259))
      CALL FFV1_0(W(63),W(3),W(7),UV_GQQB_1EPS,AMPL(2,260))
      CALL FFV1_0(W(63),W(3),W(7),UV_GQQB_1EPS,AMPL(2,261))
      CALL FFV1_0(W(63),W(3),W(7),UV_GQQB_1EPS,AMPL(2,262))
      CALL FFV1_0(W(63),W(3),W(7),UV_GQQB_1EPS,AMPL(2,263))
      CALL FFV1_0(W(63),W(3),W(7),UV_GQQB_1EPS,AMPL(2,264))
      CALL FFV1_0(W(63),W(3),W(7),UV_GQQG_1EPS,AMPL(2,265))
      CALL FFV1_0(W(63),W(3),W(7),UV_GQQB,AMPL(1,266))
      CALL FFV1_0(W(63),W(3),W(7),UV_GQQT,AMPL(1,267))
C     Counter-term amplitude(s) for loop diagram number 217
      CALL FFV2_3_0(W(71),W(3),W(14),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,268))
C     Counter-term amplitude(s) for loop diagram number 218
      CALL FFV2_3_0(W(72),W(3),W(20),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,269))
C     Counter-term amplitude(s) for loop diagram number 220
      CALL FFV1_0(W(64),W(3),W(7),R2_GQQ,AMPL(1,270))
      CALL FFV1_0(W(64),W(3),W(7),UV_GQQB_1EPS,AMPL(2,271))
      CALL FFV1_0(W(64),W(3),W(7),UV_GQQB_1EPS,AMPL(2,272))
      CALL FFV1_0(W(64),W(3),W(7),UV_GQQB_1EPS,AMPL(2,273))
      CALL FFV1_0(W(64),W(3),W(7),UV_GQQB_1EPS,AMPL(2,274))
      CALL FFV1_0(W(64),W(3),W(7),UV_GQQB_1EPS,AMPL(2,275))
      CALL FFV1_0(W(64),W(3),W(7),UV_GQQB_1EPS,AMPL(2,276))
      CALL FFV1_0(W(64),W(3),W(7),UV_GQQG_1EPS,AMPL(2,277))
      CALL FFV1_0(W(64),W(3),W(7),UV_GQQB,AMPL(1,278))
      CALL FFV1_0(W(64),W(3),W(7),UV_GQQT,AMPL(1,279))
C     Counter-term amplitude(s) for loop diagram number 222
      CALL FFV1_0(W(65),W(3),W(7),R2_GQQ,AMPL(1,280))
      CALL FFV1_0(W(65),W(3),W(7),UV_GQQB_1EPS,AMPL(2,281))
      CALL FFV1_0(W(65),W(3),W(7),UV_GQQB_1EPS,AMPL(2,282))
      CALL FFV1_0(W(65),W(3),W(7),UV_GQQB_1EPS,AMPL(2,283))
      CALL FFV1_0(W(65),W(3),W(7),UV_GQQB_1EPS,AMPL(2,284))
      CALL FFV1_0(W(65),W(3),W(7),UV_GQQB_1EPS,AMPL(2,285))
      CALL FFV1_0(W(65),W(3),W(7),UV_GQQB_1EPS,AMPL(2,286))
      CALL FFV1_0(W(65),W(3),W(7),UV_GQQG_1EPS,AMPL(2,287))
      CALL FFV1_0(W(65),W(3),W(7),UV_GQQB,AMPL(1,288))
      CALL FFV1_0(W(65),W(3),W(7),UV_GQQT,AMPL(1,289))
C     Counter-term amplitude(s) for loop diagram number 224
      CALL FFV1_0(W(66),W(3),W(7),R2_GQQ,AMPL(1,290))
      CALL FFV1_0(W(66),W(3),W(7),UV_GQQB_1EPS,AMPL(2,291))
      CALL FFV1_0(W(66),W(3),W(7),UV_GQQB_1EPS,AMPL(2,292))
      CALL FFV1_0(W(66),W(3),W(7),UV_GQQB_1EPS,AMPL(2,293))
      CALL FFV1_0(W(66),W(3),W(7),UV_GQQB_1EPS,AMPL(2,294))
      CALL FFV1_0(W(66),W(3),W(7),UV_GQQB_1EPS,AMPL(2,295))
      CALL FFV1_0(W(66),W(3),W(7),UV_GQQB_1EPS,AMPL(2,296))
      CALL FFV1_0(W(66),W(3),W(7),UV_GQQG_1EPS,AMPL(2,297))
      CALL FFV1_0(W(66),W(3),W(7),UV_GQQB,AMPL(1,298))
      CALL FFV1_0(W(66),W(3),W(7),UV_GQQT,AMPL(1,299))
      CALL FFV1_1(W(3),W(39),GC_1,MDL_MB,ZERO,W(73))
C     Counter-term amplitude(s) for loop diagram number 267
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(73),R2_QQQ,R2_QQB,AMPL(1,300))
      CALL R2_QQ_2_0(W(24),W(73),UV_BMASS_1EPS,AMPL(2,301))
      CALL R2_QQ_2_0(W(24),W(73),UV_BMASS,AMPL(1,302))
      CALL FFV2_3_1(W(3),W(41),GC_21,GC_23,MDL_MB,ZERO,W(74))
C     Counter-term amplitude(s) for loop diagram number 268
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(74),R2_QQQ,R2_QQB,AMPL(1,303))
      CALL R2_QQ_2_0(W(24),W(74),UV_BMASS_1EPS,AMPL(2,304))
      CALL R2_QQ_2_0(W(24),W(74),UV_BMASS,AMPL(1,305))
      CALL FFV1_1(W(3),W(43),GC_1,MDL_MB,ZERO,W(75))
C     Counter-term amplitude(s) for loop diagram number 269
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(75),R2_QQQ,R2_QQB,AMPL(1,306))
      CALL R2_QQ_2_0(W(24),W(75),UV_BMASS_1EPS,AMPL(2,307))
      CALL R2_QQ_2_0(W(24),W(75),UV_BMASS,AMPL(1,308))
      CALL FFV2_3_1(W(3),W(45),GC_21,GC_23,MDL_MB,ZERO,W(76))
C     Counter-term amplitude(s) for loop diagram number 270
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(76),R2_QQQ,R2_QQB,AMPL(1,309))
      CALL R2_QQ_2_0(W(24),W(76),UV_BMASS_1EPS,AMPL(2,310))
      CALL R2_QQ_2_0(W(24),W(76),UV_BMASS,AMPL(1,311))
      CALL FFV1_1(W(3),W(47),GC_1,MDL_MB,ZERO,W(77))
C     Counter-term amplitude(s) for loop diagram number 271
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(77),R2_QQQ,R2_QQB,AMPL(1,312))
      CALL R2_QQ_2_0(W(24),W(77),UV_BMASS_1EPS,AMPL(2,313))
      CALL R2_QQ_2_0(W(24),W(77),UV_BMASS,AMPL(1,314))
      CALL FFV2_3_1(W(3),W(49),GC_21,GC_23,MDL_MB,ZERO,W(78))
C     Counter-term amplitude(s) for loop diagram number 272
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(78),R2_QQQ,R2_QQB,AMPL(1,315))
      CALL R2_QQ_2_0(W(24),W(78),UV_BMASS_1EPS,AMPL(2,316))
      CALL R2_QQ_2_0(W(24),W(78),UV_BMASS,AMPL(1,317))
      CALL FFV1_1(W(3),W(51),GC_1,MDL_MB,ZERO,W(79))
C     Counter-term amplitude(s) for loop diagram number 273
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(79),R2_QQQ,R2_QQB,AMPL(1,318))
      CALL R2_QQ_2_0(W(24),W(79),UV_BMASS_1EPS,AMPL(2,319))
      CALL R2_QQ_2_0(W(24),W(79),UV_BMASS,AMPL(1,320))
      CALL FFV2_3_1(W(3),W(53),GC_21,GC_23,MDL_MB,ZERO,W(80))
C     Counter-term amplitude(s) for loop diagram number 274
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(80),R2_QQQ,R2_QQB,AMPL(1,321))
      CALL R2_QQ_2_0(W(24),W(80),UV_BMASS_1EPS,AMPL(2,322))
      CALL R2_QQ_2_0(W(24),W(80),UV_BMASS,AMPL(1,323))
C     Counter-term amplitude(s) for loop diagram number 275
      CALL R2_QQ_1_R2_QQ_2_0(W(67),W(27),R2_QQQ,R2_QQB,AMPL(1,324))
      CALL R2_QQ_2_0(W(67),W(27),UV_BMASS_1EPS,AMPL(2,325))
      CALL R2_QQ_2_0(W(67),W(27),UV_BMASS,AMPL(1,326))
      CALL FFV1_1(W(27),W(17),GC_1,MDL_MB,ZERO,W(81))
C     Counter-term amplitude(s) for loop diagram number 276
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(81),R2_QQQ,R2_QQB,AMPL(1,327))
      CALL R2_QQ_2_0(W(24),W(81),UV_BMASS_1EPS,AMPL(2,328))
      CALL R2_QQ_2_0(W(24),W(81),UV_BMASS,AMPL(1,329))
C     Counter-term amplitude(s) for loop diagram number 277
      CALL FFV1_0(W(24),W(27),W(17),R2_DDA,AMPL(1,330))
C     Counter-term amplitude(s) for loop diagram number 278
      CALL R2_QQ_1_R2_QQ_2_0(W(28),W(32),R2_QQQ,R2_QQB,AMPL(1,331))
      CALL R2_QQ_2_0(W(28),W(32),UV_BMASS_1EPS,AMPL(2,332))
      CALL R2_QQ_2_0(W(28),W(32),UV_BMASS,AMPL(1,333))
      CALL FFV1_1(W(32),W(8),GC_1,MDL_MB,ZERO,W(82))
C     Counter-term amplitude(s) for loop diagram number 279
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(82),R2_QQQ,R2_QQB,AMPL(1,334))
      CALL R2_QQ_2_0(W(24),W(82),UV_BMASS_1EPS,AMPL(2,335))
      CALL R2_QQ_2_0(W(24),W(82),UV_BMASS,AMPL(1,336))
C     Counter-term amplitude(s) for loop diagram number 280
      CALL FFV1_0(W(24),W(32),W(8),R2_DDA,AMPL(1,337))
C     Counter-term amplitude(s) for loop diagram number 281
      CALL R2_QQ_1_R2_QQ_2_0(W(68),W(27),R2_QQQ,R2_QQB,AMPL(1,338))
      CALL R2_QQ_2_0(W(68),W(27),UV_BMASS_1EPS,AMPL(2,339))
      CALL R2_QQ_2_0(W(68),W(27),UV_BMASS,AMPL(1,340))
      CALL FFV2_3_1(W(27),W(20),GC_21,GC_23,MDL_MB,ZERO,W(83))
C     Counter-term amplitude(s) for loop diagram number 282
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(83),R2_QQQ,R2_QQB,AMPL(1,341))
      CALL R2_QQ_2_0(W(24),W(83),UV_BMASS_1EPS,AMPL(2,342))
      CALL R2_QQ_2_0(W(24),W(83),UV_BMASS,AMPL(1,343))
C     Counter-term amplitude(s) for loop diagram number 283
      CALL FFV2_3_0(W(24),W(27),W(20),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,344))
C     Counter-term amplitude(s) for loop diagram number 284
      CALL R2_QQ_1_R2_QQ_2_0(W(28),W(34),R2_QQQ,R2_QQB,AMPL(1,345))
      CALL R2_QQ_2_0(W(28),W(34),UV_BMASS_1EPS,AMPL(2,346))
      CALL R2_QQ_2_0(W(28),W(34),UV_BMASS,AMPL(1,347))
      CALL FFV1_1(W(34),W(8),GC_1,MDL_MB,ZERO,W(84))
C     Counter-term amplitude(s) for loop diagram number 285
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(84),R2_QQQ,R2_QQB,AMPL(1,348))
      CALL R2_QQ_2_0(W(24),W(84),UV_BMASS_1EPS,AMPL(2,349))
      CALL R2_QQ_2_0(W(24),W(84),UV_BMASS,AMPL(1,350))
C     Counter-term amplitude(s) for loop diagram number 286
      CALL FFV1_0(W(24),W(34),W(8),R2_DDA,AMPL(1,351))
C     Counter-term amplitude(s) for loop diagram number 287
      CALL R2_QQ_1_R2_QQ_2_0(W(67),W(29),R2_QQQ,R2_QQB,AMPL(1,352))
      CALL R2_QQ_2_0(W(67),W(29),UV_BMASS_1EPS,AMPL(2,353))
      CALL R2_QQ_2_0(W(67),W(29),UV_BMASS,AMPL(1,354))
      CALL FFV1_1(W(29),W(17),GC_1,MDL_MB,ZERO,W(85))
C     Counter-term amplitude(s) for loop diagram number 288
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(85),R2_QQQ,R2_QQB,AMPL(1,355))
      CALL R2_QQ_2_0(W(24),W(85),UV_BMASS_1EPS,AMPL(2,356))
      CALL R2_QQ_2_0(W(24),W(85),UV_BMASS,AMPL(1,357))
C     Counter-term amplitude(s) for loop diagram number 289
      CALL FFV1_0(W(24),W(29),W(17),R2_DDA,AMPL(1,358))
C     Counter-term amplitude(s) for loop diagram number 290
      CALL R2_QQ_1_R2_QQ_2_0(W(30),W(32),R2_QQQ,R2_QQB,AMPL(1,359))
      CALL R2_QQ_2_0(W(30),W(32),UV_BMASS_1EPS,AMPL(2,360))
      CALL R2_QQ_2_0(W(30),W(32),UV_BMASS,AMPL(1,361))
      CALL FFV2_3_1(W(32),W(14),GC_21,GC_23,MDL_MB,ZERO,W(86))
C     Counter-term amplitude(s) for loop diagram number 291
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(86),R2_QQQ,R2_QQB,AMPL(1,362))
      CALL R2_QQ_2_0(W(24),W(86),UV_BMASS_1EPS,AMPL(2,363))
      CALL R2_QQ_2_0(W(24),W(86),UV_BMASS,AMPL(1,364))
C     Counter-term amplitude(s) for loop diagram number 292
      CALL FFV2_3_0(W(24),W(32),W(14),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,365))
C     Counter-term amplitude(s) for loop diagram number 293
      CALL R2_QQ_1_R2_QQ_2_0(W(68),W(29),R2_QQQ,R2_QQB,AMPL(1,366))
      CALL R2_QQ_2_0(W(68),W(29),UV_BMASS_1EPS,AMPL(2,367))
      CALL R2_QQ_2_0(W(68),W(29),UV_BMASS,AMPL(1,368))
      CALL FFV2_3_1(W(29),W(20),GC_21,GC_23,MDL_MB,ZERO,W(87))
C     Counter-term amplitude(s) for loop diagram number 294
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(87),R2_QQQ,R2_QQB,AMPL(1,369))
      CALL R2_QQ_2_0(W(24),W(87),UV_BMASS_1EPS,AMPL(2,370))
      CALL R2_QQ_2_0(W(24),W(87),UV_BMASS,AMPL(1,371))
C     Counter-term amplitude(s) for loop diagram number 295
      CALL FFV2_3_0(W(24),W(29),W(20),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,372))
C     Counter-term amplitude(s) for loop diagram number 296
      CALL R2_QQ_1_R2_QQ_2_0(W(30),W(34),R2_QQQ,R2_QQB,AMPL(1,373))
      CALL R2_QQ_2_0(W(30),W(34),UV_BMASS_1EPS,AMPL(2,374))
      CALL R2_QQ_2_0(W(30),W(34),UV_BMASS,AMPL(1,375))
      CALL FFS1_1(W(3),W(57),GC_33,MDL_MB,ZERO,W(88))
C     Counter-term amplitude(s) for loop diagram number 297
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(88),R2_QQQ,R2_QQB,AMPL(1,376))
      CALL R2_QQ_2_0(W(24),W(88),UV_BMASS_1EPS,AMPL(2,377))
      CALL R2_QQ_2_0(W(24),W(88),UV_BMASS,AMPL(1,378))
      CALL FFV2_3_1(W(34),W(14),GC_21,GC_23,MDL_MB,ZERO,W(89))
C     Counter-term amplitude(s) for loop diagram number 298
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(89),R2_QQQ,R2_QQB,AMPL(1,379))
      CALL R2_QQ_2_0(W(24),W(89),UV_BMASS_1EPS,AMPL(2,380))
      CALL R2_QQ_2_0(W(24),W(89),UV_BMASS,AMPL(1,381))
C     Counter-term amplitude(s) for loop diagram number 299
      CALL FFV2_3_0(W(24),W(34),W(14),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,382))
C     Counter-term amplitude(s) for loop diagram number 300
      CALL FFV1_0(W(4),W(73),W(7),R2_GQQ,AMPL(1,383))
      CALL FFV1_0(W(4),W(73),W(7),UV_GQQB_1EPS,AMPL(2,384))
      CALL FFV1_0(W(4),W(73),W(7),UV_GQQB_1EPS,AMPL(2,385))
      CALL FFV1_0(W(4),W(73),W(7),UV_GQQB_1EPS,AMPL(2,386))
      CALL FFV1_0(W(4),W(73),W(7),UV_GQQB_1EPS,AMPL(2,387))
      CALL FFV1_0(W(4),W(73),W(7),UV_GQQB_1EPS,AMPL(2,388))
      CALL FFV1_0(W(4),W(73),W(7),UV_GQQB_1EPS,AMPL(2,389))
      CALL FFV1_0(W(4),W(73),W(7),UV_GQQG_1EPS,AMPL(2,390))
      CALL FFV1_0(W(4),W(73),W(7),UV_GQQB,AMPL(1,391))
      CALL FFV1_0(W(4),W(73),W(7),UV_GQQT,AMPL(1,392))
C     Counter-term amplitude(s) for loop diagram number 301
      CALL FFV1_0(W(4),W(74),W(7),R2_GQQ,AMPL(1,393))
      CALL FFV1_0(W(4),W(74),W(7),UV_GQQB_1EPS,AMPL(2,394))
      CALL FFV1_0(W(4),W(74),W(7),UV_GQQB_1EPS,AMPL(2,395))
      CALL FFV1_0(W(4),W(74),W(7),UV_GQQB_1EPS,AMPL(2,396))
      CALL FFV1_0(W(4),W(74),W(7),UV_GQQB_1EPS,AMPL(2,397))
      CALL FFV1_0(W(4),W(74),W(7),UV_GQQB_1EPS,AMPL(2,398))
      CALL FFV1_0(W(4),W(74),W(7),UV_GQQB_1EPS,AMPL(2,399))
      CALL FFV1_0(W(4),W(74),W(7),UV_GQQG_1EPS,AMPL(2,400))
      CALL FFV1_0(W(4),W(74),W(7),UV_GQQB,AMPL(1,401))
      CALL FFV1_0(W(4),W(74),W(7),UV_GQQT,AMPL(1,402))
C     Counter-term amplitude(s) for loop diagram number 302
      CALL FFV1_0(W(4),W(75),W(7),R2_GQQ,AMPL(1,403))
      CALL FFV1_0(W(4),W(75),W(7),UV_GQQB_1EPS,AMPL(2,404))
      CALL FFV1_0(W(4),W(75),W(7),UV_GQQB_1EPS,AMPL(2,405))
      CALL FFV1_0(W(4),W(75),W(7),UV_GQQB_1EPS,AMPL(2,406))
      CALL FFV1_0(W(4),W(75),W(7),UV_GQQB_1EPS,AMPL(2,407))
      CALL FFV1_0(W(4),W(75),W(7),UV_GQQB_1EPS,AMPL(2,408))
      CALL FFV1_0(W(4),W(75),W(7),UV_GQQB_1EPS,AMPL(2,409))
      CALL FFV1_0(W(4),W(75),W(7),UV_GQQG_1EPS,AMPL(2,410))
      CALL FFV1_0(W(4),W(75),W(7),UV_GQQB,AMPL(1,411))
      CALL FFV1_0(W(4),W(75),W(7),UV_GQQT,AMPL(1,412))
C     Counter-term amplitude(s) for loop diagram number 303
      CALL FFV1_0(W(4),W(76),W(7),R2_GQQ,AMPL(1,413))
      CALL FFV1_0(W(4),W(76),W(7),UV_GQQB_1EPS,AMPL(2,414))
      CALL FFV1_0(W(4),W(76),W(7),UV_GQQB_1EPS,AMPL(2,415))
      CALL FFV1_0(W(4),W(76),W(7),UV_GQQB_1EPS,AMPL(2,416))
      CALL FFV1_0(W(4),W(76),W(7),UV_GQQB_1EPS,AMPL(2,417))
      CALL FFV1_0(W(4),W(76),W(7),UV_GQQB_1EPS,AMPL(2,418))
      CALL FFV1_0(W(4),W(76),W(7),UV_GQQB_1EPS,AMPL(2,419))
      CALL FFV1_0(W(4),W(76),W(7),UV_GQQG_1EPS,AMPL(2,420))
      CALL FFV1_0(W(4),W(76),W(7),UV_GQQB,AMPL(1,421))
      CALL FFV1_0(W(4),W(76),W(7),UV_GQQT,AMPL(1,422))
C     Counter-term amplitude(s) for loop diagram number 304
      CALL FFV1_0(W(4),W(77),W(7),R2_GQQ,AMPL(1,423))
      CALL FFV1_0(W(4),W(77),W(7),UV_GQQB_1EPS,AMPL(2,424))
      CALL FFV1_0(W(4),W(77),W(7),UV_GQQB_1EPS,AMPL(2,425))
      CALL FFV1_0(W(4),W(77),W(7),UV_GQQB_1EPS,AMPL(2,426))
      CALL FFV1_0(W(4),W(77),W(7),UV_GQQB_1EPS,AMPL(2,427))
      CALL FFV1_0(W(4),W(77),W(7),UV_GQQB_1EPS,AMPL(2,428))
      CALL FFV1_0(W(4),W(77),W(7),UV_GQQB_1EPS,AMPL(2,429))
      CALL FFV1_0(W(4),W(77),W(7),UV_GQQG_1EPS,AMPL(2,430))
      CALL FFV1_0(W(4),W(77),W(7),UV_GQQB,AMPL(1,431))
      CALL FFV1_0(W(4),W(77),W(7),UV_GQQT,AMPL(1,432))
C     Counter-term amplitude(s) for loop diagram number 305
      CALL FFV1_0(W(4),W(78),W(7),R2_GQQ,AMPL(1,433))
      CALL FFV1_0(W(4),W(78),W(7),UV_GQQB_1EPS,AMPL(2,434))
      CALL FFV1_0(W(4),W(78),W(7),UV_GQQB_1EPS,AMPL(2,435))
      CALL FFV1_0(W(4),W(78),W(7),UV_GQQB_1EPS,AMPL(2,436))
      CALL FFV1_0(W(4),W(78),W(7),UV_GQQB_1EPS,AMPL(2,437))
      CALL FFV1_0(W(4),W(78),W(7),UV_GQQB_1EPS,AMPL(2,438))
      CALL FFV1_0(W(4),W(78),W(7),UV_GQQB_1EPS,AMPL(2,439))
      CALL FFV1_0(W(4),W(78),W(7),UV_GQQG_1EPS,AMPL(2,440))
      CALL FFV1_0(W(4),W(78),W(7),UV_GQQB,AMPL(1,441))
      CALL FFV1_0(W(4),W(78),W(7),UV_GQQT,AMPL(1,442))
C     Counter-term amplitude(s) for loop diagram number 306
      CALL FFV1_0(W(4),W(79),W(7),R2_GQQ,AMPL(1,443))
      CALL FFV1_0(W(4),W(79),W(7),UV_GQQB_1EPS,AMPL(2,444))
      CALL FFV1_0(W(4),W(79),W(7),UV_GQQB_1EPS,AMPL(2,445))
      CALL FFV1_0(W(4),W(79),W(7),UV_GQQB_1EPS,AMPL(2,446))
      CALL FFV1_0(W(4),W(79),W(7),UV_GQQB_1EPS,AMPL(2,447))
      CALL FFV1_0(W(4),W(79),W(7),UV_GQQB_1EPS,AMPL(2,448))
      CALL FFV1_0(W(4),W(79),W(7),UV_GQQB_1EPS,AMPL(2,449))
      CALL FFV1_0(W(4),W(79),W(7),UV_GQQG_1EPS,AMPL(2,450))
      CALL FFV1_0(W(4),W(79),W(7),UV_GQQB,AMPL(1,451))
      CALL FFV1_0(W(4),W(79),W(7),UV_GQQT,AMPL(1,452))
C     Counter-term amplitude(s) for loop diagram number 307
      CALL FFV1_0(W(4),W(80),W(7),R2_GQQ,AMPL(1,453))
      CALL FFV1_0(W(4),W(80),W(7),UV_GQQB_1EPS,AMPL(2,454))
      CALL FFV1_0(W(4),W(80),W(7),UV_GQQB_1EPS,AMPL(2,455))
      CALL FFV1_0(W(4),W(80),W(7),UV_GQQB_1EPS,AMPL(2,456))
      CALL FFV1_0(W(4),W(80),W(7),UV_GQQB_1EPS,AMPL(2,457))
      CALL FFV1_0(W(4),W(80),W(7),UV_GQQB_1EPS,AMPL(2,458))
      CALL FFV1_0(W(4),W(80),W(7),UV_GQQB_1EPS,AMPL(2,459))
      CALL FFV1_0(W(4),W(80),W(7),UV_GQQG_1EPS,AMPL(2,460))
      CALL FFV1_0(W(4),W(80),W(7),UV_GQQB,AMPL(1,461))
      CALL FFV1_0(W(4),W(80),W(7),UV_GQQT,AMPL(1,462))
      CALL FFV1_1(W(27),W(7),GC_5,MDL_MB,ZERO,W(90))
C     Counter-term amplitude(s) for loop diagram number 308
      CALL FFV1_0(W(4),W(90),W(17),R2_DDA,AMPL(1,463))
C     Counter-term amplitude(s) for loop diagram number 310
      CALL FFV1_0(W(4),W(81),W(7),R2_GQQ,AMPL(1,464))
      CALL FFV1_0(W(4),W(81),W(7),UV_GQQB_1EPS,AMPL(2,465))
      CALL FFV1_0(W(4),W(81),W(7),UV_GQQB_1EPS,AMPL(2,466))
      CALL FFV1_0(W(4),W(81),W(7),UV_GQQB_1EPS,AMPL(2,467))
      CALL FFV1_0(W(4),W(81),W(7),UV_GQQB_1EPS,AMPL(2,468))
      CALL FFV1_0(W(4),W(81),W(7),UV_GQQB_1EPS,AMPL(2,469))
      CALL FFV1_0(W(4),W(81),W(7),UV_GQQB_1EPS,AMPL(2,470))
      CALL FFV1_0(W(4),W(81),W(7),UV_GQQG_1EPS,AMPL(2,471))
      CALL FFV1_0(W(4),W(81),W(7),UV_GQQB,AMPL(1,472))
      CALL FFV1_0(W(4),W(81),W(7),UV_GQQT,AMPL(1,473))
      CALL FFV1_1(W(32),W(7),GC_5,MDL_MB,ZERO,W(91))
C     Counter-term amplitude(s) for loop diagram number 313
      CALL FFV1_0(W(4),W(91),W(8),R2_DDA,AMPL(1,474))
C     Counter-term amplitude(s) for loop diagram number 314
      CALL FFV1_0(W(4),W(82),W(7),R2_GQQ,AMPL(1,475))
      CALL FFV1_0(W(4),W(82),W(7),UV_GQQB_1EPS,AMPL(2,476))
      CALL FFV1_0(W(4),W(82),W(7),UV_GQQB_1EPS,AMPL(2,477))
      CALL FFV1_0(W(4),W(82),W(7),UV_GQQB_1EPS,AMPL(2,478))
      CALL FFV1_0(W(4),W(82),W(7),UV_GQQB_1EPS,AMPL(2,479))
      CALL FFV1_0(W(4),W(82),W(7),UV_GQQB_1EPS,AMPL(2,480))
      CALL FFV1_0(W(4),W(82),W(7),UV_GQQB_1EPS,AMPL(2,481))
      CALL FFV1_0(W(4),W(82),W(7),UV_GQQG_1EPS,AMPL(2,482))
      CALL FFV1_0(W(4),W(82),W(7),UV_GQQB,AMPL(1,483))
      CALL FFV1_0(W(4),W(82),W(7),UV_GQQT,AMPL(1,484))
C     Counter-term amplitude(s) for loop diagram number 316
      CALL FFV2_3_0(W(4),W(90),W(20),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,485))
C     Counter-term amplitude(s) for loop diagram number 318
      CALL FFV1_0(W(4),W(83),W(7),R2_GQQ,AMPL(1,486))
      CALL FFV1_0(W(4),W(83),W(7),UV_GQQB_1EPS,AMPL(2,487))
      CALL FFV1_0(W(4),W(83),W(7),UV_GQQB_1EPS,AMPL(2,488))
      CALL FFV1_0(W(4),W(83),W(7),UV_GQQB_1EPS,AMPL(2,489))
      CALL FFV1_0(W(4),W(83),W(7),UV_GQQB_1EPS,AMPL(2,490))
      CALL FFV1_0(W(4),W(83),W(7),UV_GQQB_1EPS,AMPL(2,491))
      CALL FFV1_0(W(4),W(83),W(7),UV_GQQB_1EPS,AMPL(2,492))
      CALL FFV1_0(W(4),W(83),W(7),UV_GQQG_1EPS,AMPL(2,493))
      CALL FFV1_0(W(4),W(83),W(7),UV_GQQB,AMPL(1,494))
      CALL FFV1_0(W(4),W(83),W(7),UV_GQQT,AMPL(1,495))
      CALL FFV1_1(W(34),W(7),GC_5,MDL_MB,ZERO,W(92))
C     Counter-term amplitude(s) for loop diagram number 321
      CALL FFV1_0(W(4),W(92),W(8),R2_DDA,AMPL(1,496))
C     Counter-term amplitude(s) for loop diagram number 322
      CALL FFV1_0(W(4),W(84),W(7),R2_GQQ,AMPL(1,497))
      CALL FFV1_0(W(4),W(84),W(7),UV_GQQB_1EPS,AMPL(2,498))
      CALL FFV1_0(W(4),W(84),W(7),UV_GQQB_1EPS,AMPL(2,499))
      CALL FFV1_0(W(4),W(84),W(7),UV_GQQB_1EPS,AMPL(2,500))
      CALL FFV1_0(W(4),W(84),W(7),UV_GQQB_1EPS,AMPL(2,501))
      CALL FFV1_0(W(4),W(84),W(7),UV_GQQB_1EPS,AMPL(2,502))
      CALL FFV1_0(W(4),W(84),W(7),UV_GQQB_1EPS,AMPL(2,503))
      CALL FFV1_0(W(4),W(84),W(7),UV_GQQG_1EPS,AMPL(2,504))
      CALL FFV1_0(W(4),W(84),W(7),UV_GQQB,AMPL(1,505))
      CALL FFV1_0(W(4),W(84),W(7),UV_GQQT,AMPL(1,506))
      CALL FFV1_1(W(29),W(7),GC_5,MDL_MB,ZERO,W(93))
C     Counter-term amplitude(s) for loop diagram number 324
      CALL FFV1_0(W(4),W(93),W(17),R2_DDA,AMPL(1,507))
C     Counter-term amplitude(s) for loop diagram number 326
      CALL FFV1_0(W(4),W(85),W(7),R2_GQQ,AMPL(1,508))
      CALL FFV1_0(W(4),W(85),W(7),UV_GQQB_1EPS,AMPL(2,509))
      CALL FFV1_0(W(4),W(85),W(7),UV_GQQB_1EPS,AMPL(2,510))
      CALL FFV1_0(W(4),W(85),W(7),UV_GQQB_1EPS,AMPL(2,511))
      CALL FFV1_0(W(4),W(85),W(7),UV_GQQB_1EPS,AMPL(2,512))
      CALL FFV1_0(W(4),W(85),W(7),UV_GQQB_1EPS,AMPL(2,513))
      CALL FFV1_0(W(4),W(85),W(7),UV_GQQB_1EPS,AMPL(2,514))
      CALL FFV1_0(W(4),W(85),W(7),UV_GQQG_1EPS,AMPL(2,515))
      CALL FFV1_0(W(4),W(85),W(7),UV_GQQB,AMPL(1,516))
      CALL FFV1_0(W(4),W(85),W(7),UV_GQQT,AMPL(1,517))
C     Counter-term amplitude(s) for loop diagram number 329
      CALL FFV2_3_0(W(4),W(91),W(14),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,518))
C     Counter-term amplitude(s) for loop diagram number 330
      CALL FFV1_0(W(4),W(86),W(7),R2_GQQ,AMPL(1,519))
      CALL FFV1_0(W(4),W(86),W(7),UV_GQQB_1EPS,AMPL(2,520))
      CALL FFV1_0(W(4),W(86),W(7),UV_GQQB_1EPS,AMPL(2,521))
      CALL FFV1_0(W(4),W(86),W(7),UV_GQQB_1EPS,AMPL(2,522))
      CALL FFV1_0(W(4),W(86),W(7),UV_GQQB_1EPS,AMPL(2,523))
      CALL FFV1_0(W(4),W(86),W(7),UV_GQQB_1EPS,AMPL(2,524))
      CALL FFV1_0(W(4),W(86),W(7),UV_GQQB_1EPS,AMPL(2,525))
      CALL FFV1_0(W(4),W(86),W(7),UV_GQQG_1EPS,AMPL(2,526))
      CALL FFV1_0(W(4),W(86),W(7),UV_GQQB,AMPL(1,527))
      CALL FFV1_0(W(4),W(86),W(7),UV_GQQT,AMPL(1,528))
C     Counter-term amplitude(s) for loop diagram number 332
      CALL FFV2_3_0(W(4),W(93),W(20),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,529))
C     Counter-term amplitude(s) for loop diagram number 334
      CALL FFV1_0(W(4),W(87),W(7),R2_GQQ,AMPL(1,530))
      CALL FFV1_0(W(4),W(87),W(7),UV_GQQB_1EPS,AMPL(2,531))
      CALL FFV1_0(W(4),W(87),W(7),UV_GQQB_1EPS,AMPL(2,532))
      CALL FFV1_0(W(4),W(87),W(7),UV_GQQB_1EPS,AMPL(2,533))
      CALL FFV1_0(W(4),W(87),W(7),UV_GQQB_1EPS,AMPL(2,534))
      CALL FFV1_0(W(4),W(87),W(7),UV_GQQB_1EPS,AMPL(2,535))
      CALL FFV1_0(W(4),W(87),W(7),UV_GQQB_1EPS,AMPL(2,536))
      CALL FFV1_0(W(4),W(87),W(7),UV_GQQG_1EPS,AMPL(2,537))
      CALL FFV1_0(W(4),W(87),W(7),UV_GQQB,AMPL(1,538))
      CALL FFV1_0(W(4),W(87),W(7),UV_GQQT,AMPL(1,539))
C     Counter-term amplitude(s) for loop diagram number 337
      CALL FFV2_3_0(W(4),W(92),W(14),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,540))
C     Counter-term amplitude(s) for loop diagram number 338
      CALL FFV1_0(W(4),W(88),W(7),R2_GQQ,AMPL(1,541))
      CALL FFV1_0(W(4),W(88),W(7),UV_GQQB_1EPS,AMPL(2,542))
      CALL FFV1_0(W(4),W(88),W(7),UV_GQQB_1EPS,AMPL(2,543))
      CALL FFV1_0(W(4),W(88),W(7),UV_GQQB_1EPS,AMPL(2,544))
      CALL FFV1_0(W(4),W(88),W(7),UV_GQQB_1EPS,AMPL(2,545))
      CALL FFV1_0(W(4),W(88),W(7),UV_GQQB_1EPS,AMPL(2,546))
      CALL FFV1_0(W(4),W(88),W(7),UV_GQQB_1EPS,AMPL(2,547))
      CALL FFV1_0(W(4),W(88),W(7),UV_GQQG_1EPS,AMPL(2,548))
      CALL FFV1_0(W(4),W(88),W(7),UV_GQQB,AMPL(1,549))
      CALL FFV1_0(W(4),W(88),W(7),UV_GQQT,AMPL(1,550))
C     Counter-term amplitude(s) for loop diagram number 339
      CALL FFV1_0(W(4),W(89),W(7),R2_GQQ,AMPL(1,551))
      CALL FFV1_0(W(4),W(89),W(7),UV_GQQB_1EPS,AMPL(2,552))
      CALL FFV1_0(W(4),W(89),W(7),UV_GQQB_1EPS,AMPL(2,553))
      CALL FFV1_0(W(4),W(89),W(7),UV_GQQB_1EPS,AMPL(2,554))
      CALL FFV1_0(W(4),W(89),W(7),UV_GQQB_1EPS,AMPL(2,555))
      CALL FFV1_0(W(4),W(89),W(7),UV_GQQB_1EPS,AMPL(2,556))
      CALL FFV1_0(W(4),W(89),W(7),UV_GQQB_1EPS,AMPL(2,557))
      CALL FFV1_0(W(4),W(89),W(7),UV_GQQG_1EPS,AMPL(2,558))
      CALL FFV1_0(W(4),W(89),W(7),UV_GQQB,AMPL(1,559))
      CALL FFV1_0(W(4),W(89),W(7),UV_GQQT,AMPL(1,560))
C     Counter-term amplitude(s) for loop diagram number 366
      CALL R2_QQ_1_R2_QQ_2_0(W(31),W(90),R2_QQQ,R2_QQB,AMPL(1,561))
      CALL R2_QQ_2_0(W(31),W(90),UV_BMASS_1EPS,AMPL(2,562))
      CALL R2_QQ_2_0(W(31),W(90),UV_BMASS,AMPL(1,563))
C     Counter-term amplitude(s) for loop diagram number 367
      CALL R2_QQ_1_R2_QQ_2_0(W(69),W(27),R2_QQQ,R2_QQB,AMPL(1,564))
      CALL R2_QQ_2_0(W(69),W(27),UV_BMASS_1EPS,AMPL(2,565))
      CALL R2_QQ_2_0(W(69),W(27),UV_BMASS,AMPL(1,566))
C     Counter-term amplitude(s) for loop diagram number 368
      CALL FFV1_0(W(31),W(27),W(7),R2_GQQ,AMPL(1,567))
      CALL FFV1_0(W(31),W(27),W(7),UV_GQQB_1EPS,AMPL(2,568))
      CALL FFV1_0(W(31),W(27),W(7),UV_GQQB_1EPS,AMPL(2,569))
      CALL FFV1_0(W(31),W(27),W(7),UV_GQQB_1EPS,AMPL(2,570))
      CALL FFV1_0(W(31),W(27),W(7),UV_GQQB_1EPS,AMPL(2,571))
      CALL FFV1_0(W(31),W(27),W(7),UV_GQQB_1EPS,AMPL(2,572))
      CALL FFV1_0(W(31),W(27),W(7),UV_GQQB_1EPS,AMPL(2,573))
      CALL FFV1_0(W(31),W(27),W(7),UV_GQQG_1EPS,AMPL(2,574))
      CALL FFV1_0(W(31),W(27),W(7),UV_GQQB,AMPL(1,575))
      CALL FFV1_0(W(31),W(27),W(7),UV_GQQT,AMPL(1,576))
C     Counter-term amplitude(s) for loop diagram number 369
      CALL R2_QQ_1_R2_QQ_2_0(W(70),W(32),R2_QQQ,R2_QQB,AMPL(1,577))
      CALL R2_QQ_2_0(W(70),W(32),UV_BMASS_1EPS,AMPL(2,578))
      CALL R2_QQ_2_0(W(70),W(32),UV_BMASS,AMPL(1,579))
C     Counter-term amplitude(s) for loop diagram number 370
      CALL R2_QQ_1_R2_QQ_2_0(W(19),W(91),R2_QQQ,R2_QQB,AMPL(1,580))
      CALL R2_QQ_2_0(W(19),W(91),UV_BMASS_1EPS,AMPL(2,581))
      CALL R2_QQ_2_0(W(19),W(91),UV_BMASS,AMPL(1,582))
C     Counter-term amplitude(s) for loop diagram number 371
      CALL FFV1_0(W(19),W(32),W(7),R2_GQQ,AMPL(1,583))
      CALL FFV1_0(W(19),W(32),W(7),UV_GQQB_1EPS,AMPL(2,584))
      CALL FFV1_0(W(19),W(32),W(7),UV_GQQB_1EPS,AMPL(2,585))
      CALL FFV1_0(W(19),W(32),W(7),UV_GQQB_1EPS,AMPL(2,586))
      CALL FFV1_0(W(19),W(32),W(7),UV_GQQB_1EPS,AMPL(2,587))
      CALL FFV1_0(W(19),W(32),W(7),UV_GQQB_1EPS,AMPL(2,588))
      CALL FFV1_0(W(19),W(32),W(7),UV_GQQB_1EPS,AMPL(2,589))
      CALL FFV1_0(W(19),W(32),W(7),UV_GQQG_1EPS,AMPL(2,590))
      CALL FFV1_0(W(19),W(32),W(7),UV_GQQB,AMPL(1,591))
      CALL FFV1_0(W(19),W(32),W(7),UV_GQQT,AMPL(1,592))
C     Counter-term amplitude(s) for loop diagram number 372
      CALL R2_QQ_1_R2_QQ_2_0(W(33),W(90),R2_QQQ,R2_QQB,AMPL(1,593))
      CALL R2_QQ_2_0(W(33),W(90),UV_BMASS_1EPS,AMPL(2,594))
      CALL R2_QQ_2_0(W(33),W(90),UV_BMASS,AMPL(1,595))
C     Counter-term amplitude(s) for loop diagram number 373
      CALL R2_QQ_1_R2_QQ_2_0(W(71),W(27),R2_QQQ,R2_QQB,AMPL(1,596))
      CALL R2_QQ_2_0(W(71),W(27),UV_BMASS_1EPS,AMPL(2,597))
      CALL R2_QQ_2_0(W(71),W(27),UV_BMASS,AMPL(1,598))
C     Counter-term amplitude(s) for loop diagram number 374
      CALL FFV1_0(W(33),W(27),W(7),R2_GQQ,AMPL(1,599))
      CALL FFV1_0(W(33),W(27),W(7),UV_GQQB_1EPS,AMPL(2,600))
      CALL FFV1_0(W(33),W(27),W(7),UV_GQQB_1EPS,AMPL(2,601))
      CALL FFV1_0(W(33),W(27),W(7),UV_GQQB_1EPS,AMPL(2,602))
      CALL FFV1_0(W(33),W(27),W(7),UV_GQQB_1EPS,AMPL(2,603))
      CALL FFV1_0(W(33),W(27),W(7),UV_GQQB_1EPS,AMPL(2,604))
      CALL FFV1_0(W(33),W(27),W(7),UV_GQQB_1EPS,AMPL(2,605))
      CALL FFV1_0(W(33),W(27),W(7),UV_GQQG_1EPS,AMPL(2,606))
      CALL FFV1_0(W(33),W(27),W(7),UV_GQQB,AMPL(1,607))
      CALL FFV1_0(W(33),W(27),W(7),UV_GQQT,AMPL(1,608))
C     Counter-term amplitude(s) for loop diagram number 375
      CALL R2_QQ_1_R2_QQ_2_0(W(70),W(34),R2_QQQ,R2_QQB,AMPL(1,609))
      CALL R2_QQ_2_0(W(70),W(34),UV_BMASS_1EPS,AMPL(2,610))
      CALL R2_QQ_2_0(W(70),W(34),UV_BMASS,AMPL(1,611))
C     Counter-term amplitude(s) for loop diagram number 376
      CALL R2_QQ_1_R2_QQ_2_0(W(19),W(92),R2_QQQ,R2_QQB,AMPL(1,612))
      CALL R2_QQ_2_0(W(19),W(92),UV_BMASS_1EPS,AMPL(2,613))
      CALL R2_QQ_2_0(W(19),W(92),UV_BMASS,AMPL(1,614))
C     Counter-term amplitude(s) for loop diagram number 377
      CALL FFV1_0(W(19),W(34),W(7),R2_GQQ,AMPL(1,615))
      CALL FFV1_0(W(19),W(34),W(7),UV_GQQB_1EPS,AMPL(2,616))
      CALL FFV1_0(W(19),W(34),W(7),UV_GQQB_1EPS,AMPL(2,617))
      CALL FFV1_0(W(19),W(34),W(7),UV_GQQB_1EPS,AMPL(2,618))
      CALL FFV1_0(W(19),W(34),W(7),UV_GQQB_1EPS,AMPL(2,619))
      CALL FFV1_0(W(19),W(34),W(7),UV_GQQB_1EPS,AMPL(2,620))
      CALL FFV1_0(W(19),W(34),W(7),UV_GQQB_1EPS,AMPL(2,621))
      CALL FFV1_0(W(19),W(34),W(7),UV_GQQG_1EPS,AMPL(2,622))
      CALL FFV1_0(W(19),W(34),W(7),UV_GQQB,AMPL(1,623))
      CALL FFV1_0(W(19),W(34),W(7),UV_GQQT,AMPL(1,624))
C     Counter-term amplitude(s) for loop diagram number 378
      CALL R2_QQ_1_R2_QQ_2_0(W(31),W(93),R2_QQQ,R2_QQB,AMPL(1,625))
      CALL R2_QQ_2_0(W(31),W(93),UV_BMASS_1EPS,AMPL(2,626))
      CALL R2_QQ_2_0(W(31),W(93),UV_BMASS,AMPL(1,627))
C     Counter-term amplitude(s) for loop diagram number 379
      CALL R2_QQ_1_R2_QQ_2_0(W(69),W(29),R2_QQQ,R2_QQB,AMPL(1,628))
      CALL R2_QQ_2_0(W(69),W(29),UV_BMASS_1EPS,AMPL(2,629))
      CALL R2_QQ_2_0(W(69),W(29),UV_BMASS,AMPL(1,630))
C     Counter-term amplitude(s) for loop diagram number 380
      CALL FFV1_0(W(31),W(29),W(7),R2_GQQ,AMPL(1,631))
      CALL FFV1_0(W(31),W(29),W(7),UV_GQQB_1EPS,AMPL(2,632))
      CALL FFV1_0(W(31),W(29),W(7),UV_GQQB_1EPS,AMPL(2,633))
      CALL FFV1_0(W(31),W(29),W(7),UV_GQQB_1EPS,AMPL(2,634))
      CALL FFV1_0(W(31),W(29),W(7),UV_GQQB_1EPS,AMPL(2,635))
      CALL FFV1_0(W(31),W(29),W(7),UV_GQQB_1EPS,AMPL(2,636))
      CALL FFV1_0(W(31),W(29),W(7),UV_GQQB_1EPS,AMPL(2,637))
      CALL FFV1_0(W(31),W(29),W(7),UV_GQQG_1EPS,AMPL(2,638))
      CALL FFV1_0(W(31),W(29),W(7),UV_GQQB,AMPL(1,639))
      CALL FFV1_0(W(31),W(29),W(7),UV_GQQT,AMPL(1,640))
C     Counter-term amplitude(s) for loop diagram number 381
      CALL R2_QQ_1_R2_QQ_2_0(W(72),W(32),R2_QQQ,R2_QQB,AMPL(1,641))
      CALL R2_QQ_2_0(W(72),W(32),UV_BMASS_1EPS,AMPL(2,642))
      CALL R2_QQ_2_0(W(72),W(32),UV_BMASS,AMPL(1,643))
C     Counter-term amplitude(s) for loop diagram number 382
      CALL R2_QQ_1_R2_QQ_2_0(W(22),W(91),R2_QQQ,R2_QQB,AMPL(1,644))
      CALL R2_QQ_2_0(W(22),W(91),UV_BMASS_1EPS,AMPL(2,645))
      CALL R2_QQ_2_0(W(22),W(91),UV_BMASS,AMPL(1,646))
C     Counter-term amplitude(s) for loop diagram number 383
      CALL FFV1_0(W(22),W(32),W(7),R2_GQQ,AMPL(1,647))
      CALL FFV1_0(W(22),W(32),W(7),UV_GQQB_1EPS,AMPL(2,648))
      CALL FFV1_0(W(22),W(32),W(7),UV_GQQB_1EPS,AMPL(2,649))
      CALL FFV1_0(W(22),W(32),W(7),UV_GQQB_1EPS,AMPL(2,650))
      CALL FFV1_0(W(22),W(32),W(7),UV_GQQB_1EPS,AMPL(2,651))
      CALL FFV1_0(W(22),W(32),W(7),UV_GQQB_1EPS,AMPL(2,652))
      CALL FFV1_0(W(22),W(32),W(7),UV_GQQB_1EPS,AMPL(2,653))
      CALL FFV1_0(W(22),W(32),W(7),UV_GQQG_1EPS,AMPL(2,654))
      CALL FFV1_0(W(22),W(32),W(7),UV_GQQB,AMPL(1,655))
      CALL FFV1_0(W(22),W(32),W(7),UV_GQQT,AMPL(1,656))
C     Counter-term amplitude(s) for loop diagram number 384
      CALL R2_QQ_1_R2_QQ_2_0(W(33),W(93),R2_QQQ,R2_QQB,AMPL(1,657))
      CALL R2_QQ_2_0(W(33),W(93),UV_BMASS_1EPS,AMPL(2,658))
      CALL R2_QQ_2_0(W(33),W(93),UV_BMASS,AMPL(1,659))
C     Counter-term amplitude(s) for loop diagram number 385
      CALL R2_QQ_1_R2_QQ_2_0(W(71),W(29),R2_QQQ,R2_QQB,AMPL(1,660))
      CALL R2_QQ_2_0(W(71),W(29),UV_BMASS_1EPS,AMPL(2,661))
      CALL R2_QQ_2_0(W(71),W(29),UV_BMASS,AMPL(1,662))
C     Counter-term amplitude(s) for loop diagram number 386
      CALL FFV1_0(W(33),W(29),W(7),R2_GQQ,AMPL(1,663))
      CALL FFV1_0(W(33),W(29),W(7),UV_GQQB_1EPS,AMPL(2,664))
      CALL FFV1_0(W(33),W(29),W(7),UV_GQQB_1EPS,AMPL(2,665))
      CALL FFV1_0(W(33),W(29),W(7),UV_GQQB_1EPS,AMPL(2,666))
      CALL FFV1_0(W(33),W(29),W(7),UV_GQQB_1EPS,AMPL(2,667))
      CALL FFV1_0(W(33),W(29),W(7),UV_GQQB_1EPS,AMPL(2,668))
      CALL FFV1_0(W(33),W(29),W(7),UV_GQQB_1EPS,AMPL(2,669))
      CALL FFV1_0(W(33),W(29),W(7),UV_GQQG_1EPS,AMPL(2,670))
      CALL FFV1_0(W(33),W(29),W(7),UV_GQQB,AMPL(1,671))
      CALL FFV1_0(W(33),W(29),W(7),UV_GQQT,AMPL(1,672))
C     Counter-term amplitude(s) for loop diagram number 387
      CALL R2_QQ_1_R2_QQ_2_0(W(72),W(34),R2_QQQ,R2_QQB,AMPL(1,673))
      CALL R2_QQ_2_0(W(72),W(34),UV_BMASS_1EPS,AMPL(2,674))
      CALL R2_QQ_2_0(W(72),W(34),UV_BMASS,AMPL(1,675))
C     Counter-term amplitude(s) for loop diagram number 388
      CALL R2_QQ_1_R2_QQ_2_0(W(22),W(92),R2_QQQ,R2_QQB,AMPL(1,676))
      CALL R2_QQ_2_0(W(22),W(92),UV_BMASS_1EPS,AMPL(2,677))
      CALL R2_QQ_2_0(W(22),W(92),UV_BMASS,AMPL(1,678))
C     Counter-term amplitude(s) for loop diagram number 389
      CALL FFV1_0(W(22),W(34),W(7),R2_GQQ,AMPL(1,679))
      CALL FFV1_0(W(22),W(34),W(7),UV_GQQB_1EPS,AMPL(2,680))
      CALL FFV1_0(W(22),W(34),W(7),UV_GQQB_1EPS,AMPL(2,681))
      CALL FFV1_0(W(22),W(34),W(7),UV_GQQB_1EPS,AMPL(2,682))
      CALL FFV1_0(W(22),W(34),W(7),UV_GQQB_1EPS,AMPL(2,683))
      CALL FFV1_0(W(22),W(34),W(7),UV_GQQB_1EPS,AMPL(2,684))
      CALL FFV1_0(W(22),W(34),W(7),UV_GQQB_1EPS,AMPL(2,685))
      CALL FFV1_0(W(22),W(34),W(7),UV_GQQG_1EPS,AMPL(2,686))
      CALL FFV1_0(W(22),W(34),W(7),UV_GQQB,AMPL(1,687))
      CALL FFV1_0(W(22),W(34),W(7),UV_GQQT,AMPL(1,688))
      CALL FFV1P0_3(W(2),W(35),GC_3,ZERO,ZERO,W(94))
C     Counter-term amplitude(s) for loop diagram number 398
      CALL FFV1_0(W(4),W(9),W(94),R2_DDA,AMPL(1,689))
      CALL FFV2_4_3(W(2),W(35),GC_21,GC_24,MDL_MZ,MDL_WZ,W(95))
C     Counter-term amplitude(s) for loop diagram number 399
      CALL FFV2_3_0(W(4),W(9),W(95),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,690))
      CALL FFV1P0_3(W(36),W(1),GC_3,ZERO,ZERO,W(96))
C     Counter-term amplitude(s) for loop diagram number 400
      CALL FFV1_0(W(4),W(9),W(96),R2_DDA,AMPL(1,691))
      CALL FFV2_4_3(W(36),W(1),GC_21,GC_24,MDL_MZ,MDL_WZ,W(97))
C     Counter-term amplitude(s) for loop diagram number 401
      CALL FFV2_3_0(W(4),W(9),W(97),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,692))
      CALL FFV1P0_3(W(2),W(37),GC_3,ZERO,ZERO,W(98))
C     Counter-term amplitude(s) for loop diagram number 402
      CALL FFV1_0(W(4),W(9),W(98),R2_DDA,AMPL(1,693))
      CALL FFV2_4_3(W(2),W(37),GC_21,GC_24,MDL_MZ,MDL_WZ,W(99))
C     Counter-term amplitude(s) for loop diagram number 403
      CALL FFV2_3_0(W(4),W(9),W(99),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,694))
      CALL FFV1P0_3(W(38),W(1),GC_3,ZERO,ZERO,W(100))
C     Counter-term amplitude(s) for loop diagram number 404
      CALL FFV1_0(W(4),W(9),W(100),R2_DDA,AMPL(1,695))
      CALL FFV2_4_3(W(38),W(1),GC_21,GC_24,MDL_MZ,MDL_WZ,W(101))
C     Counter-term amplitude(s) for loop diagram number 405
      CALL FFV2_3_0(W(4),W(9),W(101),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,696))
      CALL FFV1_2(W(4),W(94),GC_1,MDL_MB,ZERO,W(102))
C     Counter-term amplitude(s) for loop diagram number 406
      CALL R2_QQ_1_R2_QQ_2_0(W(102),W(9),R2_QQQ,R2_QQB,AMPL(1,697))
      CALL R2_QQ_2_0(W(102),W(9),UV_BMASS_1EPS,AMPL(2,698))
      CALL R2_QQ_2_0(W(102),W(9),UV_BMASS,AMPL(1,699))
      CALL FFV2_3_2(W(4),W(95),GC_21,GC_23,MDL_MB,ZERO,W(103))
C     Counter-term amplitude(s) for loop diagram number 407
      CALL R2_QQ_1_R2_QQ_2_0(W(103),W(9),R2_QQQ,R2_QQB,AMPL(1,700))
      CALL R2_QQ_2_0(W(103),W(9),UV_BMASS_1EPS,AMPL(2,701))
      CALL R2_QQ_2_0(W(103),W(9),UV_BMASS,AMPL(1,702))
      CALL FFV1_2(W(4),W(96),GC_1,MDL_MB,ZERO,W(104))
C     Counter-term amplitude(s) for loop diagram number 408
      CALL R2_QQ_1_R2_QQ_2_0(W(104),W(9),R2_QQQ,R2_QQB,AMPL(1,703))
      CALL R2_QQ_2_0(W(104),W(9),UV_BMASS_1EPS,AMPL(2,704))
      CALL R2_QQ_2_0(W(104),W(9),UV_BMASS,AMPL(1,705))
      CALL FFV2_3_2(W(4),W(97),GC_21,GC_23,MDL_MB,ZERO,W(105))
C     Counter-term amplitude(s) for loop diagram number 409
      CALL R2_QQ_1_R2_QQ_2_0(W(105),W(9),R2_QQQ,R2_QQB,AMPL(1,706))
      CALL R2_QQ_2_0(W(105),W(9),UV_BMASS_1EPS,AMPL(2,707))
      CALL R2_QQ_2_0(W(105),W(9),UV_BMASS,AMPL(1,708))
      CALL FFV1_2(W(4),W(98),GC_1,MDL_MB,ZERO,W(106))
C     Counter-term amplitude(s) for loop diagram number 410
      CALL R2_QQ_1_R2_QQ_2_0(W(106),W(9),R2_QQQ,R2_QQB,AMPL(1,709))
      CALL R2_QQ_2_0(W(106),W(9),UV_BMASS_1EPS,AMPL(2,710))
      CALL R2_QQ_2_0(W(106),W(9),UV_BMASS,AMPL(1,711))
      CALL FFV2_3_2(W(4),W(99),GC_21,GC_23,MDL_MB,ZERO,W(107))
C     Counter-term amplitude(s) for loop diagram number 411
      CALL R2_QQ_1_R2_QQ_2_0(W(107),W(9),R2_QQQ,R2_QQB,AMPL(1,712))
      CALL R2_QQ_2_0(W(107),W(9),UV_BMASS_1EPS,AMPL(2,713))
      CALL R2_QQ_2_0(W(107),W(9),UV_BMASS,AMPL(1,714))
      CALL FFV1_2(W(4),W(100),GC_1,MDL_MB,ZERO,W(108))
C     Counter-term amplitude(s) for loop diagram number 412
      CALL R2_QQ_1_R2_QQ_2_0(W(108),W(9),R2_QQQ,R2_QQB,AMPL(1,715))
      CALL R2_QQ_2_0(W(108),W(9),UV_BMASS_1EPS,AMPL(2,716))
      CALL R2_QQ_2_0(W(108),W(9),UV_BMASS,AMPL(1,717))
      CALL FFV2_3_2(W(4),W(101),GC_21,GC_23,MDL_MB,ZERO,W(109))
C     Counter-term amplitude(s) for loop diagram number 413
      CALL R2_QQ_1_R2_QQ_2_0(W(109),W(9),R2_QQQ,R2_QQB,AMPL(1,718))
      CALL R2_QQ_2_0(W(109),W(9),UV_BMASS_1EPS,AMPL(2,719))
      CALL R2_QQ_2_0(W(109),W(9),UV_BMASS,AMPL(1,720))
C     Counter-term amplitude(s) for loop diagram number 414
      CALL FFV1_0(W(24),W(3),W(94),R2_DDA,AMPL(1,721))
C     Counter-term amplitude(s) for loop diagram number 415
      CALL FFV2_3_0(W(24),W(3),W(95),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,722))
C     Counter-term amplitude(s) for loop diagram number 416
      CALL FFV1_0(W(24),W(3),W(96),R2_DDA,AMPL(1,723))
C     Counter-term amplitude(s) for loop diagram number 417
      CALL FFV2_3_0(W(24),W(3),W(97),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,724))
C     Counter-term amplitude(s) for loop diagram number 418
      CALL FFV1_0(W(24),W(3),W(98),R2_DDA,AMPL(1,725))
C     Counter-term amplitude(s) for loop diagram number 419
      CALL FFV2_3_0(W(24),W(3),W(99),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,726))
C     Counter-term amplitude(s) for loop diagram number 420
      CALL FFV1_0(W(24),W(3),W(100),R2_DDA,AMPL(1,727))
C     Counter-term amplitude(s) for loop diagram number 421
      CALL FFV2_3_0(W(24),W(3),W(101),R2_DDZ_V2,R2_DDZ_V3,AMPL(1,728))
C     Counter-term amplitude(s) for loop diagram number 438
      CALL FFV1_0(W(102),W(3),W(7),R2_GQQ,AMPL(1,729))
      CALL FFV1_0(W(102),W(3),W(7),UV_GQQB_1EPS,AMPL(2,730))
      CALL FFV1_0(W(102),W(3),W(7),UV_GQQB_1EPS,AMPL(2,731))
      CALL FFV1_0(W(102),W(3),W(7),UV_GQQB_1EPS,AMPL(2,732))
      CALL FFV1_0(W(102),W(3),W(7),UV_GQQB_1EPS,AMPL(2,733))
      CALL FFV1_0(W(102),W(3),W(7),UV_GQQB_1EPS,AMPL(2,734))
      CALL FFV1_0(W(102),W(3),W(7),UV_GQQB_1EPS,AMPL(2,735))
      CALL FFV1_0(W(102),W(3),W(7),UV_GQQG_1EPS,AMPL(2,736))
      CALL FFV1_0(W(102),W(3),W(7),UV_GQQB,AMPL(1,737))
      CALL FFV1_0(W(102),W(3),W(7),UV_GQQT,AMPL(1,738))
C     Counter-term amplitude(s) for loop diagram number 439
      CALL FFV1_0(W(103),W(3),W(7),R2_GQQ,AMPL(1,739))
      CALL FFV1_0(W(103),W(3),W(7),UV_GQQB_1EPS,AMPL(2,740))
      CALL FFV1_0(W(103),W(3),W(7),UV_GQQB_1EPS,AMPL(2,741))
      CALL FFV1_0(W(103),W(3),W(7),UV_GQQB_1EPS,AMPL(2,742))
      CALL FFV1_0(W(103),W(3),W(7),UV_GQQB_1EPS,AMPL(2,743))
      CALL FFV1_0(W(103),W(3),W(7),UV_GQQB_1EPS,AMPL(2,744))
      CALL FFV1_0(W(103),W(3),W(7),UV_GQQB_1EPS,AMPL(2,745))
      CALL FFV1_0(W(103),W(3),W(7),UV_GQQG_1EPS,AMPL(2,746))
      CALL FFV1_0(W(103),W(3),W(7),UV_GQQB,AMPL(1,747))
      CALL FFV1_0(W(103),W(3),W(7),UV_GQQT,AMPL(1,748))
C     Counter-term amplitude(s) for loop diagram number 440
      CALL FFV1_0(W(104),W(3),W(7),R2_GQQ,AMPL(1,749))
      CALL FFV1_0(W(104),W(3),W(7),UV_GQQB_1EPS,AMPL(2,750))
      CALL FFV1_0(W(104),W(3),W(7),UV_GQQB_1EPS,AMPL(2,751))
      CALL FFV1_0(W(104),W(3),W(7),UV_GQQB_1EPS,AMPL(2,752))
      CALL FFV1_0(W(104),W(3),W(7),UV_GQQB_1EPS,AMPL(2,753))
      CALL FFV1_0(W(104),W(3),W(7),UV_GQQB_1EPS,AMPL(2,754))
      CALL FFV1_0(W(104),W(3),W(7),UV_GQQB_1EPS,AMPL(2,755))
      CALL FFV1_0(W(104),W(3),W(7),UV_GQQG_1EPS,AMPL(2,756))
      CALL FFV1_0(W(104),W(3),W(7),UV_GQQB,AMPL(1,757))
      CALL FFV1_0(W(104),W(3),W(7),UV_GQQT,AMPL(1,758))
C     Counter-term amplitude(s) for loop diagram number 441
      CALL FFV1_0(W(105),W(3),W(7),R2_GQQ,AMPL(1,759))
      CALL FFV1_0(W(105),W(3),W(7),UV_GQQB_1EPS,AMPL(2,760))
      CALL FFV1_0(W(105),W(3),W(7),UV_GQQB_1EPS,AMPL(2,761))
      CALL FFV1_0(W(105),W(3),W(7),UV_GQQB_1EPS,AMPL(2,762))
      CALL FFV1_0(W(105),W(3),W(7),UV_GQQB_1EPS,AMPL(2,763))
      CALL FFV1_0(W(105),W(3),W(7),UV_GQQB_1EPS,AMPL(2,764))
      CALL FFV1_0(W(105),W(3),W(7),UV_GQQB_1EPS,AMPL(2,765))
      CALL FFV1_0(W(105),W(3),W(7),UV_GQQG_1EPS,AMPL(2,766))
      CALL FFV1_0(W(105),W(3),W(7),UV_GQQB,AMPL(1,767))
      CALL FFV1_0(W(105),W(3),W(7),UV_GQQT,AMPL(1,768))
C     Counter-term amplitude(s) for loop diagram number 442
      CALL FFV1_0(W(106),W(3),W(7),R2_GQQ,AMPL(1,769))
      CALL FFV1_0(W(106),W(3),W(7),UV_GQQB_1EPS,AMPL(2,770))
      CALL FFV1_0(W(106),W(3),W(7),UV_GQQB_1EPS,AMPL(2,771))
      CALL FFV1_0(W(106),W(3),W(7),UV_GQQB_1EPS,AMPL(2,772))
      CALL FFV1_0(W(106),W(3),W(7),UV_GQQB_1EPS,AMPL(2,773))
      CALL FFV1_0(W(106),W(3),W(7),UV_GQQB_1EPS,AMPL(2,774))
      CALL FFV1_0(W(106),W(3),W(7),UV_GQQB_1EPS,AMPL(2,775))
      CALL FFV1_0(W(106),W(3),W(7),UV_GQQG_1EPS,AMPL(2,776))
      CALL FFV1_0(W(106),W(3),W(7),UV_GQQB,AMPL(1,777))
      CALL FFV1_0(W(106),W(3),W(7),UV_GQQT,AMPL(1,778))
C     Counter-term amplitude(s) for loop diagram number 443
      CALL FFV1_0(W(107),W(3),W(7),R2_GQQ,AMPL(1,779))
      CALL FFV1_0(W(107),W(3),W(7),UV_GQQB_1EPS,AMPL(2,780))
      CALL FFV1_0(W(107),W(3),W(7),UV_GQQB_1EPS,AMPL(2,781))
      CALL FFV1_0(W(107),W(3),W(7),UV_GQQB_1EPS,AMPL(2,782))
      CALL FFV1_0(W(107),W(3),W(7),UV_GQQB_1EPS,AMPL(2,783))
      CALL FFV1_0(W(107),W(3),W(7),UV_GQQB_1EPS,AMPL(2,784))
      CALL FFV1_0(W(107),W(3),W(7),UV_GQQB_1EPS,AMPL(2,785))
      CALL FFV1_0(W(107),W(3),W(7),UV_GQQG_1EPS,AMPL(2,786))
      CALL FFV1_0(W(107),W(3),W(7),UV_GQQB,AMPL(1,787))
      CALL FFV1_0(W(107),W(3),W(7),UV_GQQT,AMPL(1,788))
C     Counter-term amplitude(s) for loop diagram number 444
      CALL FFV1_0(W(108),W(3),W(7),R2_GQQ,AMPL(1,789))
      CALL FFV1_0(W(108),W(3),W(7),UV_GQQB_1EPS,AMPL(2,790))
      CALL FFV1_0(W(108),W(3),W(7),UV_GQQB_1EPS,AMPL(2,791))
      CALL FFV1_0(W(108),W(3),W(7),UV_GQQB_1EPS,AMPL(2,792))
      CALL FFV1_0(W(108),W(3),W(7),UV_GQQB_1EPS,AMPL(2,793))
      CALL FFV1_0(W(108),W(3),W(7),UV_GQQB_1EPS,AMPL(2,794))
      CALL FFV1_0(W(108),W(3),W(7),UV_GQQB_1EPS,AMPL(2,795))
      CALL FFV1_0(W(108),W(3),W(7),UV_GQQG_1EPS,AMPL(2,796))
      CALL FFV1_0(W(108),W(3),W(7),UV_GQQB,AMPL(1,797))
      CALL FFV1_0(W(108),W(3),W(7),UV_GQQT,AMPL(1,798))
C     Counter-term amplitude(s) for loop diagram number 445
      CALL FFV1_0(W(109),W(3),W(7),R2_GQQ,AMPL(1,799))
      CALL FFV1_0(W(109),W(3),W(7),UV_GQQB_1EPS,AMPL(2,800))
      CALL FFV1_0(W(109),W(3),W(7),UV_GQQB_1EPS,AMPL(2,801))
      CALL FFV1_0(W(109),W(3),W(7),UV_GQQB_1EPS,AMPL(2,802))
      CALL FFV1_0(W(109),W(3),W(7),UV_GQQB_1EPS,AMPL(2,803))
      CALL FFV1_0(W(109),W(3),W(7),UV_GQQB_1EPS,AMPL(2,804))
      CALL FFV1_0(W(109),W(3),W(7),UV_GQQB_1EPS,AMPL(2,805))
      CALL FFV1_0(W(109),W(3),W(7),UV_GQQG_1EPS,AMPL(2,806))
      CALL FFV1_0(W(109),W(3),W(7),UV_GQQB,AMPL(1,807))
      CALL FFV1_0(W(109),W(3),W(7),UV_GQQT,AMPL(1,808))
      CALL FFV1_1(W(3),W(94),GC_1,MDL_MB,ZERO,W(110))
C     Counter-term amplitude(s) for loop diagram number 462
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(110),R2_QQQ,R2_QQB,AMPL(1,809))
      CALL R2_QQ_2_0(W(24),W(110),UV_BMASS_1EPS,AMPL(2,810))
      CALL R2_QQ_2_0(W(24),W(110),UV_BMASS,AMPL(1,811))
      CALL FFV2_3_1(W(3),W(95),GC_21,GC_23,MDL_MB,ZERO,W(111))
C     Counter-term amplitude(s) for loop diagram number 463
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(111),R2_QQQ,R2_QQB,AMPL(1,812))
      CALL R2_QQ_2_0(W(24),W(111),UV_BMASS_1EPS,AMPL(2,813))
      CALL R2_QQ_2_0(W(24),W(111),UV_BMASS,AMPL(1,814))
      CALL FFV1_1(W(3),W(96),GC_1,MDL_MB,ZERO,W(112))
C     Counter-term amplitude(s) for loop diagram number 464
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(112),R2_QQQ,R2_QQB,AMPL(1,815))
      CALL R2_QQ_2_0(W(24),W(112),UV_BMASS_1EPS,AMPL(2,816))
      CALL R2_QQ_2_0(W(24),W(112),UV_BMASS,AMPL(1,817))
      CALL FFV2_3_1(W(3),W(97),GC_21,GC_23,MDL_MB,ZERO,W(113))
C     Counter-term amplitude(s) for loop diagram number 465
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(113),R2_QQQ,R2_QQB,AMPL(1,818))
      CALL R2_QQ_2_0(W(24),W(113),UV_BMASS_1EPS,AMPL(2,819))
      CALL R2_QQ_2_0(W(24),W(113),UV_BMASS,AMPL(1,820))
      CALL FFV1_1(W(3),W(98),GC_1,MDL_MB,ZERO,W(114))
C     Counter-term amplitude(s) for loop diagram number 466
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(114),R2_QQQ,R2_QQB,AMPL(1,821))
      CALL R2_QQ_2_0(W(24),W(114),UV_BMASS_1EPS,AMPL(2,822))
      CALL R2_QQ_2_0(W(24),W(114),UV_BMASS,AMPL(1,823))
      CALL FFV2_3_1(W(3),W(99),GC_21,GC_23,MDL_MB,ZERO,W(115))
C     Counter-term amplitude(s) for loop diagram number 467
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(115),R2_QQQ,R2_QQB,AMPL(1,824))
      CALL R2_QQ_2_0(W(24),W(115),UV_BMASS_1EPS,AMPL(2,825))
      CALL R2_QQ_2_0(W(24),W(115),UV_BMASS,AMPL(1,826))
      CALL FFV1_1(W(3),W(100),GC_1,MDL_MB,ZERO,W(116))
C     Counter-term amplitude(s) for loop diagram number 468
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(116),R2_QQQ,R2_QQB,AMPL(1,827))
      CALL R2_QQ_2_0(W(24),W(116),UV_BMASS_1EPS,AMPL(2,828))
      CALL R2_QQ_2_0(W(24),W(116),UV_BMASS,AMPL(1,829))
      CALL FFV2_3_1(W(3),W(101),GC_21,GC_23,MDL_MB,ZERO,W(117))
C     Counter-term amplitude(s) for loop diagram number 469
      CALL R2_QQ_1_R2_QQ_2_0(W(24),W(117),R2_QQQ,R2_QQB,AMPL(1,830))
      CALL R2_QQ_2_0(W(24),W(117),UV_BMASS_1EPS,AMPL(2,831))
      CALL R2_QQ_2_0(W(24),W(117),UV_BMASS,AMPL(1,832))
C     Counter-term amplitude(s) for loop diagram number 470
      CALL FFV1_0(W(4),W(110),W(7),R2_GQQ,AMPL(1,833))
      CALL FFV1_0(W(4),W(110),W(7),UV_GQQB_1EPS,AMPL(2,834))
      CALL FFV1_0(W(4),W(110),W(7),UV_GQQB_1EPS,AMPL(2,835))
      CALL FFV1_0(W(4),W(110),W(7),UV_GQQB_1EPS,AMPL(2,836))
      CALL FFV1_0(W(4),W(110),W(7),UV_GQQB_1EPS,AMPL(2,837))
      CALL FFV1_0(W(4),W(110),W(7),UV_GQQB_1EPS,AMPL(2,838))
      CALL FFV1_0(W(4),W(110),W(7),UV_GQQB_1EPS,AMPL(2,839))
      CALL FFV1_0(W(4),W(110),W(7),UV_GQQG_1EPS,AMPL(2,840))
      CALL FFV1_0(W(4),W(110),W(7),UV_GQQB,AMPL(1,841))
      CALL FFV1_0(W(4),W(110),W(7),UV_GQQT,AMPL(1,842))
C     Counter-term amplitude(s) for loop diagram number 471
      CALL FFV1_0(W(4),W(111),W(7),R2_GQQ,AMPL(1,843))
      CALL FFV1_0(W(4),W(111),W(7),UV_GQQB_1EPS,AMPL(2,844))
      CALL FFV1_0(W(4),W(111),W(7),UV_GQQB_1EPS,AMPL(2,845))
      CALL FFV1_0(W(4),W(111),W(7),UV_GQQB_1EPS,AMPL(2,846))
      CALL FFV1_0(W(4),W(111),W(7),UV_GQQB_1EPS,AMPL(2,847))
      CALL FFV1_0(W(4),W(111),W(7),UV_GQQB_1EPS,AMPL(2,848))
      CALL FFV1_0(W(4),W(111),W(7),UV_GQQB_1EPS,AMPL(2,849))
      CALL FFV1_0(W(4),W(111),W(7),UV_GQQG_1EPS,AMPL(2,850))
      CALL FFV1_0(W(4),W(111),W(7),UV_GQQB,AMPL(1,851))
      CALL FFV1_0(W(4),W(111),W(7),UV_GQQT,AMPL(1,852))
C     Counter-term amplitude(s) for loop diagram number 472
      CALL FFV1_0(W(4),W(112),W(7),R2_GQQ,AMPL(1,853))
      CALL FFV1_0(W(4),W(112),W(7),UV_GQQB_1EPS,AMPL(2,854))
      CALL FFV1_0(W(4),W(112),W(7),UV_GQQB_1EPS,AMPL(2,855))
      CALL FFV1_0(W(4),W(112),W(7),UV_GQQB_1EPS,AMPL(2,856))
      CALL FFV1_0(W(4),W(112),W(7),UV_GQQB_1EPS,AMPL(2,857))
      CALL FFV1_0(W(4),W(112),W(7),UV_GQQB_1EPS,AMPL(2,858))
      CALL FFV1_0(W(4),W(112),W(7),UV_GQQB_1EPS,AMPL(2,859))
      CALL FFV1_0(W(4),W(112),W(7),UV_GQQG_1EPS,AMPL(2,860))
      CALL FFV1_0(W(4),W(112),W(7),UV_GQQB,AMPL(1,861))
      CALL FFV1_0(W(4),W(112),W(7),UV_GQQT,AMPL(1,862))
C     Counter-term amplitude(s) for loop diagram number 473
      CALL FFV1_0(W(4),W(113),W(7),R2_GQQ,AMPL(1,863))
      CALL FFV1_0(W(4),W(113),W(7),UV_GQQB_1EPS,AMPL(2,864))
      CALL FFV1_0(W(4),W(113),W(7),UV_GQQB_1EPS,AMPL(2,865))
      CALL FFV1_0(W(4),W(113),W(7),UV_GQQB_1EPS,AMPL(2,866))
      CALL FFV1_0(W(4),W(113),W(7),UV_GQQB_1EPS,AMPL(2,867))
      CALL FFV1_0(W(4),W(113),W(7),UV_GQQB_1EPS,AMPL(2,868))
      CALL FFV1_0(W(4),W(113),W(7),UV_GQQB_1EPS,AMPL(2,869))
      CALL FFV1_0(W(4),W(113),W(7),UV_GQQG_1EPS,AMPL(2,870))
      CALL FFV1_0(W(4),W(113),W(7),UV_GQQB,AMPL(1,871))
      CALL FFV1_0(W(4),W(113),W(7),UV_GQQT,AMPL(1,872))
C     Counter-term amplitude(s) for loop diagram number 474
      CALL FFV1_0(W(4),W(114),W(7),R2_GQQ,AMPL(1,873))
      CALL FFV1_0(W(4),W(114),W(7),UV_GQQB_1EPS,AMPL(2,874))
      CALL FFV1_0(W(4),W(114),W(7),UV_GQQB_1EPS,AMPL(2,875))
      CALL FFV1_0(W(4),W(114),W(7),UV_GQQB_1EPS,AMPL(2,876))
      CALL FFV1_0(W(4),W(114),W(7),UV_GQQB_1EPS,AMPL(2,877))
      CALL FFV1_0(W(4),W(114),W(7),UV_GQQB_1EPS,AMPL(2,878))
      CALL FFV1_0(W(4),W(114),W(7),UV_GQQB_1EPS,AMPL(2,879))
      CALL FFV1_0(W(4),W(114),W(7),UV_GQQG_1EPS,AMPL(2,880))
      CALL FFV1_0(W(4),W(114),W(7),UV_GQQB,AMPL(1,881))
      CALL FFV1_0(W(4),W(114),W(7),UV_GQQT,AMPL(1,882))
C     Counter-term amplitude(s) for loop diagram number 475
      CALL FFV1_0(W(4),W(115),W(7),R2_GQQ,AMPL(1,883))
      CALL FFV1_0(W(4),W(115),W(7),UV_GQQB_1EPS,AMPL(2,884))
      CALL FFV1_0(W(4),W(115),W(7),UV_GQQB_1EPS,AMPL(2,885))
      CALL FFV1_0(W(4),W(115),W(7),UV_GQQB_1EPS,AMPL(2,886))
      CALL FFV1_0(W(4),W(115),W(7),UV_GQQB_1EPS,AMPL(2,887))
      CALL FFV1_0(W(4),W(115),W(7),UV_GQQB_1EPS,AMPL(2,888))
      CALL FFV1_0(W(4),W(115),W(7),UV_GQQB_1EPS,AMPL(2,889))
      CALL FFV1_0(W(4),W(115),W(7),UV_GQQG_1EPS,AMPL(2,890))
      CALL FFV1_0(W(4),W(115),W(7),UV_GQQB,AMPL(1,891))
      CALL FFV1_0(W(4),W(115),W(7),UV_GQQT,AMPL(1,892))
C     Counter-term amplitude(s) for loop diagram number 476
      CALL FFV1_0(W(4),W(116),W(7),R2_GQQ,AMPL(1,893))
      CALL FFV1_0(W(4),W(116),W(7),UV_GQQB_1EPS,AMPL(2,894))
      CALL FFV1_0(W(4),W(116),W(7),UV_GQQB_1EPS,AMPL(2,895))
      CALL FFV1_0(W(4),W(116),W(7),UV_GQQB_1EPS,AMPL(2,896))
      CALL FFV1_0(W(4),W(116),W(7),UV_GQQB_1EPS,AMPL(2,897))
      CALL FFV1_0(W(4),W(116),W(7),UV_GQQB_1EPS,AMPL(2,898))
      CALL FFV1_0(W(4),W(116),W(7),UV_GQQB_1EPS,AMPL(2,899))
      CALL FFV1_0(W(4),W(116),W(7),UV_GQQG_1EPS,AMPL(2,900))
      CALL FFV1_0(W(4),W(116),W(7),UV_GQQB,AMPL(1,901))
      CALL FFV1_0(W(4),W(116),W(7),UV_GQQT,AMPL(1,902))
C     Counter-term amplitude(s) for loop diagram number 477
      CALL FFV1_0(W(4),W(117),W(7),R2_GQQ,AMPL(1,903))
      CALL FFV1_0(W(4),W(117),W(7),UV_GQQB_1EPS,AMPL(2,904))
      CALL FFV1_0(W(4),W(117),W(7),UV_GQQB_1EPS,AMPL(2,905))
      CALL FFV1_0(W(4),W(117),W(7),UV_GQQB_1EPS,AMPL(2,906))
      CALL FFV1_0(W(4),W(117),W(7),UV_GQQB_1EPS,AMPL(2,907))
      CALL FFV1_0(W(4),W(117),W(7),UV_GQQB_1EPS,AMPL(2,908))
      CALL FFV1_0(W(4),W(117),W(7),UV_GQQB_1EPS,AMPL(2,909))
      CALL FFV1_0(W(4),W(117),W(7),UV_GQQG_1EPS,AMPL(2,910))
      CALL FFV1_0(W(4),W(117),W(7),UV_GQQB,AMPL(1,911))
      CALL FFV1_0(W(4),W(117),W(7),UV_GQQT,AMPL(1,912))
      CALL FFV1P0_3(W(4),W(3),GC_5,ZERO,ZERO,W(118))
C     Counter-term amplitude(s) for loop diagram number 486
      CALL R2_GGZ_0(W(7),W(118),W(41),R2_GGZDOWN,AMPL(1,913))
      CALL R2_GGZ_0(W(7),W(118),W(41),R2_GGZDOWN,AMPL(1,914))
C     Counter-term amplitude(s) for loop diagram number 488
      CALL R2_GGZ_0(W(7),W(118),W(45),R2_GGZDOWN,AMPL(1,915))
      CALL R2_GGZ_0(W(7),W(118),W(45),R2_GGZDOWN,AMPL(1,916))
C     Counter-term amplitude(s) for loop diagram number 490
      CALL R2_GGZ_0(W(7),W(118),W(49),R2_GGZDOWN,AMPL(1,917))
      CALL R2_GGZ_0(W(7),W(118),W(49),R2_GGZDOWN,AMPL(1,918))
C     Counter-term amplitude(s) for loop diagram number 492
      CALL R2_GGZ_0(W(7),W(118),W(53),R2_GGZDOWN,AMPL(1,919))
      CALL R2_GGZ_0(W(7),W(118),W(53),R2_GGZDOWN,AMPL(1,920))
C     Counter-term amplitude(s) for loop diagram number 494
      CALL R2_GGVV_0(W(7),W(118),W(8),W(17),R2_GGAADOWN,AMPL(1,921))
      CALL R2_GGVV_0(W(7),W(118),W(8),W(17),R2_GGAADOWN,AMPL(1,922))
C     Counter-term amplitude(s) for loop diagram number 500
      CALL R2_GGVV_0(W(7),W(118),W(20),W(8),R2_GGZADOWN,AMPL(1,923))
      CALL R2_GGVV_0(W(7),W(118),W(20),W(8),R2_GGZADOWN,AMPL(1,924))
C     Counter-term amplitude(s) for loop diagram number 506
      CALL R2_GGVV_0(W(7),W(118),W(14),W(17),R2_GGZADOWN,AMPL(1,925))
      CALL R2_GGVV_0(W(7),W(118),W(14),W(17),R2_GGZADOWN,AMPL(1,926))
C     Counter-term amplitude(s) for loop diagram number 512
      CALL R2_GGVV_0(W(7),W(118),W(14),W(20),R2_GGZZDOWN,AMPL(1,927))
      CALL R2_GGVV_0(W(7),W(118),W(14),W(20),R2_GGZZDOWN,AMPL(1,928))
      CALL FFV1P0_3(W(4),W(27),GC_5,ZERO,ZERO,W(119))
C     Counter-term amplitude(s) for loop diagram number 518
      CALL R2_GGZ_0(W(7),W(119),W(20),R2_GGZDOWN,AMPL(1,929))
      CALL R2_GGZ_0(W(7),W(119),W(20),R2_GGZDOWN,AMPL(1,930))
      CALL FFV1P0_3(W(19),W(3),GC_5,ZERO,ZERO,W(120))
C     Counter-term amplitude(s) for loop diagram number 520
      CALL R2_GGZ_0(W(7),W(120),W(20),R2_GGZDOWN,AMPL(1,931))
      CALL R2_GGZ_0(W(7),W(120),W(20),R2_GGZDOWN,AMPL(1,932))
      CALL FFV1P0_3(W(4),W(32),GC_5,ZERO,ZERO,W(121))
C     Counter-term amplitude(s) for loop diagram number 522
      CALL R2_GGZ_0(W(7),W(121),W(14),R2_GGZDOWN,AMPL(1,933))
      CALL R2_GGZ_0(W(7),W(121),W(14),R2_GGZDOWN,AMPL(1,934))
      CALL FFV1P0_3(W(31),W(3),GC_5,ZERO,ZERO,W(122))
C     Counter-term amplitude(s) for loop diagram number 523
      CALL R2_GGZ_0(W(7),W(122),W(14),R2_GGZDOWN,AMPL(1,935))
      CALL R2_GGZ_0(W(7),W(122),W(14),R2_GGZDOWN,AMPL(1,936))
      CALL FFV1P0_3(W(4),W(29),GC_5,ZERO,ZERO,W(123))
C     Counter-term amplitude(s) for loop diagram number 526
      CALL R2_GGZ_0(W(7),W(123),W(20),R2_GGZDOWN,AMPL(1,937))
      CALL R2_GGZ_0(W(7),W(123),W(20),R2_GGZDOWN,AMPL(1,938))
      CALL FFV1P0_3(W(22),W(3),GC_5,ZERO,ZERO,W(124))
C     Counter-term amplitude(s) for loop diagram number 528
      CALL R2_GGZ_0(W(7),W(124),W(20),R2_GGZDOWN,AMPL(1,939))
      CALL R2_GGZ_0(W(7),W(124),W(20),R2_GGZDOWN,AMPL(1,940))
      CALL FFV1P0_3(W(4),W(34),GC_5,ZERO,ZERO,W(125))
C     Counter-term amplitude(s) for loop diagram number 530
      CALL R2_GGZ_0(W(7),W(125),W(14),R2_GGZDOWN,AMPL(1,941))
      CALL R2_GGZ_0(W(7),W(125),W(14),R2_GGZDOWN,AMPL(1,942))
      CALL FFV1P0_3(W(33),W(3),GC_5,ZERO,ZERO,W(126))
C     Counter-term amplitude(s) for loop diagram number 531
      CALL R2_GGZ_0(W(7),W(126),W(14),R2_GGZDOWN,AMPL(1,943))
      CALL R2_GGZ_0(W(7),W(126),W(14),R2_GGZDOWN,AMPL(1,944))
C     Counter-term amplitude(s) for loop diagram number 534
      CALL R2_GGZ_0(W(7),W(118),W(95),R2_GGZDOWN,AMPL(1,945))
      CALL R2_GGZ_0(W(7),W(118),W(95),R2_GGZDOWN,AMPL(1,946))
C     Counter-term amplitude(s) for loop diagram number 536
      CALL R2_GGZ_0(W(7),W(118),W(97),R2_GGZDOWN,AMPL(1,947))
      CALL R2_GGZ_0(W(7),W(118),W(97),R2_GGZDOWN,AMPL(1,948))
C     Counter-term amplitude(s) for loop diagram number 538
      CALL R2_GGZ_0(W(7),W(118),W(99),R2_GGZDOWN,AMPL(1,949))
      CALL R2_GGZ_0(W(7),W(118),W(99),R2_GGZDOWN,AMPL(1,950))
C     Counter-term amplitude(s) for loop diagram number 540
      CALL R2_GGZ_0(W(7),W(118),W(101),R2_GGZDOWN,AMPL(1,951))
      CALL R2_GGZ_0(W(7),W(118),W(101),R2_GGZDOWN,AMPL(1,952))
C     Counter-term amplitude(s) for loop diagram number 542
      CALL R2_GGZ_0(W(7),W(118),W(41),R2_GGZUP,AMPL(1,953))
      CALL R2_GGZ_0(W(7),W(118),W(41),R2_GGZUP,AMPL(1,954))
C     Counter-term amplitude(s) for loop diagram number 544
      CALL R2_GGZ_0(W(7),W(118),W(45),R2_GGZUP,AMPL(1,955))
      CALL R2_GGZ_0(W(7),W(118),W(45),R2_GGZUP,AMPL(1,956))
C     Counter-term amplitude(s) for loop diagram number 546
      CALL R2_GGZ_0(W(7),W(118),W(49),R2_GGZUP,AMPL(1,957))
      CALL R2_GGZ_0(W(7),W(118),W(49),R2_GGZUP,AMPL(1,958))
C     Counter-term amplitude(s) for loop diagram number 548
      CALL R2_GGZ_0(W(7),W(118),W(53),R2_GGZUP,AMPL(1,959))
      CALL R2_GGZ_0(W(7),W(118),W(53),R2_GGZUP,AMPL(1,960))
C     Counter-term amplitude(s) for loop diagram number 550
      CALL R2_GGVV_0(W(7),W(118),W(8),W(17),R2_GGAAUP,AMPL(1,961))
      CALL R2_GGVV_0(W(7),W(118),W(8),W(17),R2_GGAAUP,AMPL(1,962))
C     Counter-term amplitude(s) for loop diagram number 556
      CALL R2_GGVV_0(W(7),W(118),W(20),W(8),R2_GGZAUP,AMPL(1,963))
      CALL R2_GGVV_0(W(7),W(118),W(20),W(8),R2_GGZAUP,AMPL(1,964))
C     Counter-term amplitude(s) for loop diagram number 562
      CALL R2_GGVV_0(W(7),W(118),W(14),W(17),R2_GGZAUP,AMPL(1,965))
      CALL R2_GGVV_0(W(7),W(118),W(14),W(17),R2_GGZAUP,AMPL(1,966))
C     Counter-term amplitude(s) for loop diagram number 568
      CALL R2_GGVV_0(W(7),W(118),W(14),W(20),R2_GGZZUP,AMPL(1,967))
      CALL R2_GGVV_0(W(7),W(118),W(14),W(20),R2_GGZZUP,AMPL(1,968))
C     Counter-term amplitude(s) for loop diagram number 574
      CALL R2_GGZ_0(W(7),W(119),W(20),R2_GGZUP,AMPL(1,969))
      CALL R2_GGZ_0(W(7),W(119),W(20),R2_GGZUP,AMPL(1,970))
C     Counter-term amplitude(s) for loop diagram number 576
      CALL R2_GGZ_0(W(7),W(120),W(20),R2_GGZUP,AMPL(1,971))
      CALL R2_GGZ_0(W(7),W(120),W(20),R2_GGZUP,AMPL(1,972))
C     Counter-term amplitude(s) for loop diagram number 578
      CALL R2_GGZ_0(W(7),W(121),W(14),R2_GGZUP,AMPL(1,973))
      CALL R2_GGZ_0(W(7),W(121),W(14),R2_GGZUP,AMPL(1,974))
C     Counter-term amplitude(s) for loop diagram number 579
      CALL R2_GGZ_0(W(7),W(122),W(14),R2_GGZUP,AMPL(1,975))
      CALL R2_GGZ_0(W(7),W(122),W(14),R2_GGZUP,AMPL(1,976))
C     Counter-term amplitude(s) for loop diagram number 582
      CALL R2_GGZ_0(W(7),W(123),W(20),R2_GGZUP,AMPL(1,977))
      CALL R2_GGZ_0(W(7),W(123),W(20),R2_GGZUP,AMPL(1,978))
C     Counter-term amplitude(s) for loop diagram number 584
      CALL R2_GGZ_0(W(7),W(124),W(20),R2_GGZUP,AMPL(1,979))
      CALL R2_GGZ_0(W(7),W(124),W(20),R2_GGZUP,AMPL(1,980))
C     Counter-term amplitude(s) for loop diagram number 586
      CALL R2_GGZ_0(W(7),W(125),W(14),R2_GGZUP,AMPL(1,981))
      CALL R2_GGZ_0(W(7),W(125),W(14),R2_GGZUP,AMPL(1,982))
C     Counter-term amplitude(s) for loop diagram number 587
      CALL R2_GGZ_0(W(7),W(126),W(14),R2_GGZUP,AMPL(1,983))
      CALL R2_GGZ_0(W(7),W(126),W(14),R2_GGZUP,AMPL(1,984))
C     Counter-term amplitude(s) for loop diagram number 590
      CALL R2_GGZ_0(W(7),W(118),W(95),R2_GGZUP,AMPL(1,985))
      CALL R2_GGZ_0(W(7),W(118),W(95),R2_GGZUP,AMPL(1,986))
C     Counter-term amplitude(s) for loop diagram number 592
      CALL R2_GGZ_0(W(7),W(118),W(97),R2_GGZUP,AMPL(1,987))
      CALL R2_GGZ_0(W(7),W(118),W(97),R2_GGZUP,AMPL(1,988))
C     Counter-term amplitude(s) for loop diagram number 594
      CALL R2_GGZ_0(W(7),W(118),W(99),R2_GGZUP,AMPL(1,989))
      CALL R2_GGZ_0(W(7),W(118),W(99),R2_GGZUP,AMPL(1,990))
C     Counter-term amplitude(s) for loop diagram number 596
      CALL R2_GGZ_0(W(7),W(118),W(101),R2_GGZUP,AMPL(1,991))
      CALL R2_GGZ_0(W(7),W(118),W(101),R2_GGZUP,AMPL(1,992))
C     Counter-term amplitude(s) for loop diagram number 598
      CALL R2_GGZ_0(W(7),W(118),W(41),R2_GGZDOWN,AMPL(1,993))
C     Counter-term amplitude(s) for loop diagram number 600
      CALL R2_GGZ_0(W(7),W(118),W(45),R2_GGZDOWN,AMPL(1,994))
C     Counter-term amplitude(s) for loop diagram number 602
      CALL R2_GGZ_0(W(7),W(118),W(49),R2_GGZDOWN,AMPL(1,995))
C     Counter-term amplitude(s) for loop diagram number 604
      CALL R2_GGZ_0(W(7),W(118),W(53),R2_GGZDOWN,AMPL(1,996))
C     Counter-term amplitude(s) for loop diagram number 606
      CALL R2_GGVV_0(W(7),W(118),W(8),W(17),R2_GGAADOWN,AMPL(1,997))
C     Counter-term amplitude(s) for loop diagram number 612
      CALL R2_GGVV_0(W(7),W(118),W(20),W(8),R2_GGZADOWN,AMPL(1,998))
C     Counter-term amplitude(s) for loop diagram number 618
      CALL R2_GGVV_0(W(7),W(118),W(14),W(17),R2_GGZADOWN,AMPL(1,999))
C     Counter-term amplitude(s) for loop diagram number 624
      CALL VVS1_0(W(7),W(118),W(57),R2_GGHB,AMPL(1,1000))
C     Counter-term amplitude(s) for loop diagram number 626
      CALL R2_GGVV_0(W(7),W(118),W(14),W(20),R2_GGZZDOWN,AMPL(1,1001))
C     Counter-term amplitude(s) for loop diagram number 632
      CALL R2_GGZ_0(W(7),W(119),W(20),R2_GGZDOWN,AMPL(1,1002))
C     Counter-term amplitude(s) for loop diagram number 634
      CALL R2_GGZ_0(W(7),W(120),W(20),R2_GGZDOWN,AMPL(1,1003))
C     Counter-term amplitude(s) for loop diagram number 636
      CALL R2_GGZ_0(W(7),W(121),W(14),R2_GGZDOWN,AMPL(1,1004))
C     Counter-term amplitude(s) for loop diagram number 637
      CALL R2_GGZ_0(W(7),W(122),W(14),R2_GGZDOWN,AMPL(1,1005))
C     Counter-term amplitude(s) for loop diagram number 640
      CALL R2_GGZ_0(W(7),W(123),W(20),R2_GGZDOWN,AMPL(1,1006))
C     Counter-term amplitude(s) for loop diagram number 642
      CALL R2_GGZ_0(W(7),W(124),W(20),R2_GGZDOWN,AMPL(1,1007))
C     Counter-term amplitude(s) for loop diagram number 644
      CALL R2_GGZ_0(W(7),W(125),W(14),R2_GGZDOWN,AMPL(1,1008))
C     Counter-term amplitude(s) for loop diagram number 645
      CALL R2_GGZ_0(W(7),W(126),W(14),R2_GGZDOWN,AMPL(1,1009))
C     Counter-term amplitude(s) for loop diagram number 648
      CALL R2_GGZ_0(W(7),W(118),W(95),R2_GGZDOWN,AMPL(1,1010))
C     Counter-term amplitude(s) for loop diagram number 650
      CALL R2_GGZ_0(W(7),W(118),W(97),R2_GGZDOWN,AMPL(1,1011))
C     Counter-term amplitude(s) for loop diagram number 652
      CALL R2_GGZ_0(W(7),W(118),W(99),R2_GGZDOWN,AMPL(1,1012))
C     Counter-term amplitude(s) for loop diagram number 654
      CALL R2_GGZ_0(W(7),W(118),W(101),R2_GGZDOWN,AMPL(1,1013))
C     Counter-term amplitude(s) for loop diagram number 656
      CALL R2_GGZ_0(W(7),W(118),W(41),R2_GGZUP,AMPL(1,1014))
C     Counter-term amplitude(s) for loop diagram number 658
      CALL R2_GGZ_0(W(7),W(118),W(45),R2_GGZUP,AMPL(1,1015))
C     Counter-term amplitude(s) for loop diagram number 660
      CALL R2_GGZ_0(W(7),W(118),W(49),R2_GGZUP,AMPL(1,1016))
C     Counter-term amplitude(s) for loop diagram number 662
      CALL R2_GGZ_0(W(7),W(118),W(53),R2_GGZUP,AMPL(1,1017))
C     Counter-term amplitude(s) for loop diagram number 664
      CALL R2_GGVV_0(W(7),W(118),W(8),W(17),R2_GGAAUP,AMPL(1,1018))
C     Counter-term amplitude(s) for loop diagram number 670
      CALL R2_GGVV_0(W(7),W(118),W(20),W(8),R2_GGZAUP,AMPL(1,1019))
C     Counter-term amplitude(s) for loop diagram number 676
      CALL R2_GGVV_0(W(7),W(118),W(14),W(17),R2_GGZAUP,AMPL(1,1020))
C     Counter-term amplitude(s) for loop diagram number 682
      CALL VVS1_0(W(7),W(118),W(57),R2_GGHT,AMPL(1,1021))
C     Counter-term amplitude(s) for loop diagram number 684
      CALL R2_GGVV_0(W(7),W(118),W(14),W(20),R2_GGZZUP,AMPL(1,1022))
C     Counter-term amplitude(s) for loop diagram number 690
      CALL R2_GGZ_0(W(7),W(119),W(20),R2_GGZUP,AMPL(1,1023))
C     Counter-term amplitude(s) for loop diagram number 692
      CALL R2_GGZ_0(W(7),W(120),W(20),R2_GGZUP,AMPL(1,1024))
C     Counter-term amplitude(s) for loop diagram number 694
      CALL R2_GGZ_0(W(7),W(121),W(14),R2_GGZUP,AMPL(1,1025))
C     Counter-term amplitude(s) for loop diagram number 695
      CALL R2_GGZ_0(W(7),W(122),W(14),R2_GGZUP,AMPL(1,1026))
C     Counter-term amplitude(s) for loop diagram number 698
      CALL R2_GGZ_0(W(7),W(123),W(20),R2_GGZUP,AMPL(1,1027))
C     Counter-term amplitude(s) for loop diagram number 700
      CALL R2_GGZ_0(W(7),W(124),W(20),R2_GGZUP,AMPL(1,1028))
C     Counter-term amplitude(s) for loop diagram number 702
      CALL R2_GGZ_0(W(7),W(125),W(14),R2_GGZUP,AMPL(1,1029))
C     Counter-term amplitude(s) for loop diagram number 703
      CALL R2_GGZ_0(W(7),W(126),W(14),R2_GGZUP,AMPL(1,1030))
C     Counter-term amplitude(s) for loop diagram number 706
      CALL R2_GGZ_0(W(7),W(118),W(95),R2_GGZUP,AMPL(1,1031))
C     Counter-term amplitude(s) for loop diagram number 708
      CALL R2_GGZ_0(W(7),W(118),W(97),R2_GGZUP,AMPL(1,1032))
C     Counter-term amplitude(s) for loop diagram number 710
      CALL R2_GGZ_0(W(7),W(118),W(99),R2_GGZUP,AMPL(1,1033))
C     Counter-term amplitude(s) for loop diagram number 712
      CALL R2_GGZ_0(W(7),W(118),W(101),R2_GGZUP,AMPL(1,1034))
C     At this point, all CT amps needed for (QCD=4 QED=8), i.e. of
C      split order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 2000

      GOTO 1001
 2000 CONTINUE
      CT_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

