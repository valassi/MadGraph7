      SUBROUTINE HELAS_CALLS_UVCT_1(P,NHEL,H,IC)
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
      IF (FILTER_SO.AND.UVCT_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

C     Amplitude(s) for UVCT diagram with ID 714
      CALL FFV1_0(W(6),W(10),W(11),GC_3,AMPL(2,1035))
      AMPL(2,1035)=AMPL(2,1035)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 715
      CALL FFV1_0(W(6),W(10),W(11),GC_3,AMPL(1,1036))
      AMPL(1,1036)=AMPL(1,1036)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 716
      CALL FFV2_4_0(W(6),W(10),W(12),GC_21,GC_24,AMPL(2,1037))
      AMPL(2,1037)=AMPL(2,1037)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 717
      CALL FFV2_4_0(W(6),W(10),W(12),GC_21,GC_24,AMPL(1,1038))
      AMPL(1,1038)=AMPL(1,1038)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 718
      CALL FFV1_0(W(13),W(5),W(11),GC_3,AMPL(2,1039))
      AMPL(2,1039)=AMPL(2,1039)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 719
      CALL FFV1_0(W(13),W(5),W(11),GC_3,AMPL(1,1040))
      AMPL(1,1040)=AMPL(1,1040)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 720
      CALL FFV2_4_0(W(13),W(5),W(12),GC_21,GC_24,AMPL(2,1041))
      AMPL(2,1041)=AMPL(2,1041)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 721
      CALL FFV2_4_0(W(13),W(5),W(12),GC_21,GC_24,AMPL(1,1042))
      AMPL(1,1042)=AMPL(1,1042)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 722
      CALL FFV1_0(W(6),W(15),W(11),GC_3,AMPL(2,1043))
      AMPL(2,1043)=AMPL(2,1043)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 723
      CALL FFV1_0(W(6),W(15),W(11),GC_3,AMPL(1,1044))
      AMPL(1,1044)=AMPL(1,1044)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 724
      CALL FFV2_4_0(W(6),W(15),W(12),GC_21,GC_24,AMPL(2,1045))
      AMPL(2,1045)=AMPL(2,1045)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 725
      CALL FFV2_4_0(W(6),W(15),W(12),GC_21,GC_24,AMPL(1,1046))
      AMPL(1,1046)=AMPL(1,1046)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 726
      CALL FFV1_0(W(16),W(5),W(11),GC_3,AMPL(2,1047))
      AMPL(2,1047)=AMPL(2,1047)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 727
      CALL FFV1_0(W(16),W(5),W(11),GC_3,AMPL(1,1048))
      AMPL(1,1048)=AMPL(1,1048)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 728
      CALL FFV2_4_0(W(16),W(5),W(12),GC_21,GC_24,AMPL(2,1049))
      AMPL(2,1049)=AMPL(2,1049)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 729
      CALL FFV2_4_0(W(16),W(5),W(12),GC_21,GC_24,AMPL(1,1050))
      AMPL(1,1050)=AMPL(1,1050)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 730
      CALL FFV1_0(W(4),W(18),W(17),GC_1,AMPL(2,1051))
      AMPL(2,1051)=AMPL(2,1051)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 731
      CALL FFV1_0(W(4),W(18),W(17),GC_1,AMPL(1,1052))
      AMPL(1,1052)=AMPL(1,1052)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 732
      CALL FFV1_0(W(19),W(9),W(17),GC_1,AMPL(2,1053))
      AMPL(2,1053)=AMPL(2,1053)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 733
      CALL FFV1_0(W(19),W(9),W(17),GC_1,AMPL(1,1054))
      AMPL(1,1054)=AMPL(1,1054)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 734
      CALL FFV2_3_0(W(4),W(18),W(20),GC_21,GC_23,AMPL(2,1055))
      AMPL(2,1055)=AMPL(2,1055)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 735
      CALL FFV2_3_0(W(4),W(18),W(20),GC_21,GC_23,AMPL(1,1056))
      AMPL(1,1056)=AMPL(1,1056)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 736
      CALL FFV2_3_0(W(19),W(9),W(20),GC_21,GC_23,AMPL(2,1057))
      AMPL(2,1057)=AMPL(2,1057)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 737
      CALL FFV2_3_0(W(19),W(9),W(20),GC_21,GC_23,AMPL(1,1058))
      AMPL(1,1058)=AMPL(1,1058)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 738
      CALL FFV1_0(W(4),W(21),W(17),GC_1,AMPL(2,1059))
      AMPL(2,1059)=AMPL(2,1059)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 739
      CALL FFV1_0(W(4),W(21),W(17),GC_1,AMPL(1,1060))
      AMPL(1,1060)=AMPL(1,1060)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 740
      CALL FFV1_0(W(22),W(9),W(17),GC_1,AMPL(2,1061))
      AMPL(2,1061)=AMPL(2,1061)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 741
      CALL FFV1_0(W(22),W(9),W(17),GC_1,AMPL(1,1062))
      AMPL(1,1062)=AMPL(1,1062)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 742
      CALL FFV2_3_0(W(4),W(21),W(20),GC_21,GC_23,AMPL(2,1063))
      AMPL(2,1063)=AMPL(2,1063)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 743
      CALL FFV2_3_0(W(4),W(21),W(20),GC_21,GC_23,AMPL(1,1064))
      AMPL(1,1064)=AMPL(1,1064)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 744
      CALL FFV2_3_0(W(22),W(9),W(20),GC_21,GC_23,AMPL(2,1065))
      AMPL(2,1065)=AMPL(2,1065)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 745
      CALL FFV2_3_0(W(22),W(9),W(20),GC_21,GC_23,AMPL(1,1066))
      AMPL(1,1066)=AMPL(1,1066)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 746
      CALL FFS1_0(W(4),W(9),W(23),GC_33,AMPL(2,1067))
      AMPL(2,1067)=AMPL(2,1067)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 747
      CALL FFS1_0(W(4),W(9),W(23),GC_33,AMPL(1,1068))
      AMPL(1,1068)=AMPL(1,1068)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 748
      CALL FFV1_0(W(6),W(10),W(25),GC_3,AMPL(2,1069))
      AMPL(2,1069)=AMPL(2,1069)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 749
      CALL FFV1_0(W(6),W(10),W(25),GC_3,AMPL(1,1070))
      AMPL(1,1070)=AMPL(1,1070)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 750
      CALL FFV2_4_0(W(6),W(10),W(26),GC_21,GC_24,AMPL(2,1071))
      AMPL(2,1071)=AMPL(2,1071)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 751
      CALL FFV2_4_0(W(6),W(10),W(26),GC_21,GC_24,AMPL(1,1072))
      AMPL(1,1072)=AMPL(1,1072)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 752
      CALL FFV1_0(W(13),W(5),W(25),GC_3,AMPL(2,1073))
      AMPL(2,1073)=AMPL(2,1073)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 753
      CALL FFV1_0(W(13),W(5),W(25),GC_3,AMPL(1,1074))
      AMPL(1,1074)=AMPL(1,1074)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 754
      CALL FFV2_4_0(W(13),W(5),W(26),GC_21,GC_24,AMPL(2,1075))
      AMPL(2,1075)=AMPL(2,1075)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 755
      CALL FFV2_4_0(W(13),W(5),W(26),GC_21,GC_24,AMPL(1,1076))
      AMPL(1,1076)=AMPL(1,1076)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 756
      CALL FFV1_0(W(6),W(15),W(25),GC_3,AMPL(2,1077))
      AMPL(2,1077)=AMPL(2,1077)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 757
      CALL FFV1_0(W(6),W(15),W(25),GC_3,AMPL(1,1078))
      AMPL(1,1078)=AMPL(1,1078)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 758
      CALL FFV2_4_0(W(6),W(15),W(26),GC_21,GC_24,AMPL(2,1079))
      AMPL(2,1079)=AMPL(2,1079)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 759
      CALL FFV2_4_0(W(6),W(15),W(26),GC_21,GC_24,AMPL(1,1080))
      AMPL(1,1080)=AMPL(1,1080)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 760
      CALL FFV1_0(W(16),W(5),W(25),GC_3,AMPL(2,1081))
      AMPL(2,1081)=AMPL(2,1081)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 761
      CALL FFV1_0(W(16),W(5),W(25),GC_3,AMPL(1,1082))
      AMPL(1,1082)=AMPL(1,1082)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 762
      CALL FFV2_4_0(W(16),W(5),W(26),GC_21,GC_24,AMPL(2,1083))
      AMPL(2,1083)=AMPL(2,1083)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 763
      CALL FFV2_4_0(W(16),W(5),W(26),GC_21,GC_24,AMPL(1,1084))
      AMPL(1,1084)=AMPL(1,1084)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 764
      CALL FFV1_0(W(24),W(27),W(17),GC_1,AMPL(2,1085))
      AMPL(2,1085)=AMPL(2,1085)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 765
      CALL FFV1_0(W(24),W(27),W(17),GC_1,AMPL(1,1086))
      AMPL(1,1086)=AMPL(1,1086)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 766
      CALL FFV1_0(W(28),W(3),W(17),GC_1,AMPL(2,1087))
      AMPL(2,1087)=AMPL(2,1087)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 767
      CALL FFV1_0(W(28),W(3),W(17),GC_1,AMPL(1,1088))
      AMPL(1,1088)=AMPL(1,1088)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 768
      CALL FFV2_3_0(W(24),W(27),W(20),GC_21,GC_23,AMPL(2,1089))
      AMPL(2,1089)=AMPL(2,1089)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 769
      CALL FFV2_3_0(W(24),W(27),W(20),GC_21,GC_23,AMPL(1,1090))
      AMPL(1,1090)=AMPL(1,1090)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 770
      CALL FFV2_3_0(W(28),W(3),W(20),GC_21,GC_23,AMPL(2,1091))
      AMPL(2,1091)=AMPL(2,1091)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 771
      CALL FFV2_3_0(W(28),W(3),W(20),GC_21,GC_23,AMPL(1,1092))
      AMPL(1,1092)=AMPL(1,1092)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 772
      CALL FFV1_0(W(24),W(29),W(17),GC_1,AMPL(2,1093))
      AMPL(2,1093)=AMPL(2,1093)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 773
      CALL FFV1_0(W(24),W(29),W(17),GC_1,AMPL(1,1094))
      AMPL(1,1094)=AMPL(1,1094)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 774
      CALL FFV1_0(W(30),W(3),W(17),GC_1,AMPL(2,1095))
      AMPL(2,1095)=AMPL(2,1095)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 775
      CALL FFV1_0(W(30),W(3),W(17),GC_1,AMPL(1,1096))
      AMPL(1,1096)=AMPL(1,1096)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 776
      CALL FFV2_3_0(W(24),W(29),W(20),GC_21,GC_23,AMPL(2,1097))
      AMPL(2,1097)=AMPL(2,1097)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 777
      CALL FFV2_3_0(W(24),W(29),W(20),GC_21,GC_23,AMPL(1,1098))
      AMPL(1,1098)=AMPL(1,1098)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 778
      CALL FFV2_3_0(W(30),W(3),W(20),GC_21,GC_23,AMPL(2,1099))
      AMPL(2,1099)=AMPL(2,1099)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 779
      CALL FFV2_3_0(W(30),W(3),W(20),GC_21,GC_23,AMPL(1,1100))
      AMPL(1,1100)=AMPL(1,1100)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 780
      CALL FFS1_0(W(24),W(3),W(23),GC_33,AMPL(2,1101))
      AMPL(2,1101)=AMPL(2,1101)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 781
      CALL FFS1_0(W(24),W(3),W(23),GC_33,AMPL(1,1102))
      AMPL(1,1102)=AMPL(1,1102)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 782
      CALL FFV1_0(W(31),W(27),W(7),GC_5,AMPL(2,1103))
      AMPL(2,1103)=AMPL(2,1103)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 783
      CALL FFV1_0(W(31),W(27),W(7),GC_5,AMPL(1,1104))
      AMPL(1,1104)=AMPL(1,1104)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 784
      CALL FFV1_0(W(19),W(32),W(7),GC_5,AMPL(2,1105))
      AMPL(2,1105)=AMPL(2,1105)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 785
      CALL FFV1_0(W(19),W(32),W(7),GC_5,AMPL(1,1106))
      AMPL(1,1106)=AMPL(1,1106)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 786
      CALL FFV1_0(W(33),W(27),W(7),GC_5,AMPL(2,1107))
      AMPL(2,1107)=AMPL(2,1107)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 787
      CALL FFV1_0(W(33),W(27),W(7),GC_5,AMPL(1,1108))
      AMPL(1,1108)=AMPL(1,1108)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 788
      CALL FFV1_0(W(19),W(34),W(7),GC_5,AMPL(2,1109))
      AMPL(2,1109)=AMPL(2,1109)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 789
      CALL FFV1_0(W(19),W(34),W(7),GC_5,AMPL(1,1110))
      AMPL(1,1110)=AMPL(1,1110)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 790
      CALL FFV1_0(W(31),W(29),W(7),GC_5,AMPL(2,1111))
      AMPL(2,1111)=AMPL(2,1111)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 791
      CALL FFV1_0(W(31),W(29),W(7),GC_5,AMPL(1,1112))
      AMPL(1,1112)=AMPL(1,1112)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 792
      CALL FFV1_0(W(22),W(32),W(7),GC_5,AMPL(2,1113))
      AMPL(2,1113)=AMPL(2,1113)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 793
      CALL FFV1_0(W(22),W(32),W(7),GC_5,AMPL(1,1114))
      AMPL(1,1114)=AMPL(1,1114)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 794
      CALL FFV1_0(W(33),W(29),W(7),GC_5,AMPL(2,1115))
      AMPL(2,1115)=AMPL(2,1115)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 795
      CALL FFV1_0(W(33),W(29),W(7),GC_5,AMPL(1,1116))
      AMPL(1,1116)=AMPL(1,1116)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 796
      CALL FFV1_0(W(22),W(34),W(7),GC_5,AMPL(2,1117))
      AMPL(2,1117)=AMPL(2,1117)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 797
      CALL FFV1_0(W(22),W(34),W(7),GC_5,AMPL(1,1118))
      AMPL(1,1118)=AMPL(1,1118)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 798
      CALL FFV1_0(W(2),W(35),W(11),GC_3,AMPL(2,1119))
      AMPL(2,1119)=AMPL(2,1119)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 799
      CALL FFV1_0(W(2),W(35),W(11),GC_3,AMPL(1,1120))
      AMPL(1,1120)=AMPL(1,1120)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 800
      CALL FFV2_4_0(W(2),W(35),W(12),GC_21,GC_24,AMPL(2,1121))
      AMPL(2,1121)=AMPL(2,1121)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 801
      CALL FFV2_4_0(W(2),W(35),W(12),GC_21,GC_24,AMPL(1,1122))
      AMPL(1,1122)=AMPL(1,1122)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 802
      CALL FFV1_0(W(36),W(1),W(11),GC_3,AMPL(2,1123))
      AMPL(2,1123)=AMPL(2,1123)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 803
      CALL FFV1_0(W(36),W(1),W(11),GC_3,AMPL(1,1124))
      AMPL(1,1124)=AMPL(1,1124)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 804
      CALL FFV2_4_0(W(36),W(1),W(12),GC_21,GC_24,AMPL(2,1125))
      AMPL(2,1125)=AMPL(2,1125)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 805
      CALL FFV2_4_0(W(36),W(1),W(12),GC_21,GC_24,AMPL(1,1126))
      AMPL(1,1126)=AMPL(1,1126)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 806
      CALL FFV1_0(W(2),W(37),W(11),GC_3,AMPL(2,1127))
      AMPL(2,1127)=AMPL(2,1127)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 807
      CALL FFV1_0(W(2),W(37),W(11),GC_3,AMPL(1,1128))
      AMPL(1,1128)=AMPL(1,1128)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 808
      CALL FFV2_4_0(W(2),W(37),W(12),GC_21,GC_24,AMPL(2,1129))
      AMPL(2,1129)=AMPL(2,1129)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 809
      CALL FFV2_4_0(W(2),W(37),W(12),GC_21,GC_24,AMPL(1,1130))
      AMPL(1,1130)=AMPL(1,1130)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 810
      CALL FFV1_0(W(38),W(1),W(11),GC_3,AMPL(2,1131))
      AMPL(2,1131)=AMPL(2,1131)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 811
      CALL FFV1_0(W(38),W(1),W(11),GC_3,AMPL(1,1132))
      AMPL(1,1132)=AMPL(1,1132)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 812
      CALL FFV2_4_0(W(38),W(1),W(12),GC_21,GC_24,AMPL(2,1133))
      AMPL(2,1133)=AMPL(2,1133)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 813
      CALL FFV2_4_0(W(38),W(1),W(12),GC_21,GC_24,AMPL(1,1134))
      AMPL(1,1134)=AMPL(1,1134)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 814
      CALL FFV1_0(W(2),W(35),W(25),GC_3,AMPL(2,1135))
      AMPL(2,1135)=AMPL(2,1135)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 815
      CALL FFV1_0(W(2),W(35),W(25),GC_3,AMPL(1,1136))
      AMPL(1,1136)=AMPL(1,1136)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 816
      CALL FFV2_4_0(W(2),W(35),W(26),GC_21,GC_24,AMPL(2,1137))
      AMPL(2,1137)=AMPL(2,1137)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 817
      CALL FFV2_4_0(W(2),W(35),W(26),GC_21,GC_24,AMPL(1,1138))
      AMPL(1,1138)=AMPL(1,1138)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 818
      CALL FFV1_0(W(36),W(1),W(25),GC_3,AMPL(2,1139))
      AMPL(2,1139)=AMPL(2,1139)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 819
      CALL FFV1_0(W(36),W(1),W(25),GC_3,AMPL(1,1140))
      AMPL(1,1140)=AMPL(1,1140)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 820
      CALL FFV2_4_0(W(36),W(1),W(26),GC_21,GC_24,AMPL(2,1141))
      AMPL(2,1141)=AMPL(2,1141)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 821
      CALL FFV2_4_0(W(36),W(1),W(26),GC_21,GC_24,AMPL(1,1142))
      AMPL(1,1142)=AMPL(1,1142)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 822
      CALL FFV1_0(W(2),W(37),W(25),GC_3,AMPL(2,1143))
      AMPL(2,1143)=AMPL(2,1143)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 823
      CALL FFV1_0(W(2),W(37),W(25),GC_3,AMPL(1,1144))
      AMPL(1,1144)=AMPL(1,1144)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 824
      CALL FFV2_4_0(W(2),W(37),W(26),GC_21,GC_24,AMPL(2,1145))
      AMPL(2,1145)=AMPL(2,1145)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 825
      CALL FFV2_4_0(W(2),W(37),W(26),GC_21,GC_24,AMPL(1,1146))
      AMPL(1,1146)=AMPL(1,1146)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 826
      CALL FFV1_0(W(38),W(1),W(25),GC_3,AMPL(2,1147))
      AMPL(2,1147)=AMPL(2,1147)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 827
      CALL FFV1_0(W(38),W(1),W(25),GC_3,AMPL(1,1148))
      AMPL(1,1148)=AMPL(1,1148)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     Amplitude(s) for UVCT diagram with ID 828
      CALL FFV2_4_0(W(38),W(1),W(26),GC_21,GC_24,AMPL(2,1149))
      AMPL(2,1149)=AMPL(2,1149)*(2.0D0*UVWFCT_B_0_1EPS+2.0D0
     $ *UVWFCT_G_1_1EPS)
C     Amplitude(s) for UVCT diagram with ID 829
      CALL FFV2_4_0(W(38),W(1),W(26),GC_21,GC_24,AMPL(1,1150))
      AMPL(1,1150)=AMPL(1,1150)*(2.0D0*UVWFCT_B_0+1.0D0*UVWFCT_G_1
     $ +1.0D0*UVWFCT_G_2)
C     At this point, all UVCT amps needed for (QCD=4 QED=8), i.e. of
C      split order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 3000

      GOTO 1001
 3000 CONTINUE
      UVCT_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

