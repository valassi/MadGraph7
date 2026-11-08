      SUBROUTINE LOOP_CT_CALLS_1(P,NHEL,H,IC)
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
      PARAMETER (NEXTERNAL=5)
      INTEGER    NCOMB
      PARAMETER (NCOMB=32)
      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=18)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=585, NLOOPGROUPS=272, NCTAMPS=852)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=1437)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=48,NLOOPWAVEFUNCS=866)
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
      IF (FILTER_SO.AND.CTCALL_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

C     CutTools call for loop numbers 1,406,451
      CALL LOOP_2(7,19,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,1)
C     CutTools call for loop numbers 2,407,452
      CALL LOOP_2(6,20,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,2)
C     CutTools call for loop numbers 3,408,409,453,454
      CALL LOOP_3(5,6,7,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO),3,I_SO
     $ ,3)
C     CutTools call for loop numbers 10,11,12,62,63,64
      CALL LOOP_2_3(1,2,5,7,6,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,4)
C     CutTools call for loop numbers 13,88
      CALL LOOP_2(8,21,DCMPLX(ZERO),DCMPLX(MDL_MT),1,I_SO,5)
C     CutTools call for loop numbers 14,410,455,4,5,6
      CALL LOOP_2(6,22,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,6)
C     CutTools call for loop numbers 15
      CALL LOOP_3(4,6,8,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,7)
C     CutTools call for loop numbers 16
      CALL LOOP_3(4,6,8,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,8)
C     CutTools call for loop numbers 17,91
      CALL LOOP_3(3,5,21,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,9)
C     CutTools call for loop numbers 18
      CALL LOOP_3(3,4,19,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO),2
     $ ,I_SO,10)
C     CutTools call for loop numbers 19
      CALL LOOP_4(3,4,5,6,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),3,I_SO,11)
C     CutTools call for loop numbers 20
      CALL LOOP_4(3,5,4,6,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),3,I_SO,12)
C     CutTools call for loop numbers 24
      CALL LOOP_3(3,6,9,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,13)
C     CutTools call for loop numbers 25
      CALL LOOP_3(3,6,9,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,14)
C     CutTools call for loop numbers 26
      CALL LOOP_4(3,4,5,6,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),3,I_SO,15)
C     CutTools call for loop numbers 27
      CALL LOOP_4(3,4,6,5,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),3,I_SO,16)
C     CutTools call for loop numbers 28
      CALL LOOP_3(3,4,19,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT),2
     $ ,I_SO,17)
C     CutTools call for loop numbers 29
      CALL LOOP_4(3,5,4,6,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),3,I_SO,18)
C     CutTools call for loop numbers 30,95
      CALL LOOP_3(3,5,21,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,19)
C     CutTools call for loop numbers 31
      CALL LOOP_4(3,4,6,5,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),3,I_SO,20)
C     CutTools call for loop numbers 32,53
      CALL LOOP_2(9,23,DCMPLX(ZERO),DCMPLX(MDL_MT),1,I_SO,21)
C     CutTools call for loop numbers 33,411,456
      CALL LOOP_2(6,24,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,22)
C     CutTools call for loop numbers 34,56
      CALL LOOP_3(4,5,23,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,23)
C     CutTools call for loop numbers 35,60
      CALL LOOP_3(4,5,23,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,24)
C     CutTools call for loop numbers 36,110
      CALL LOOP_2(11,25,DCMPLX(ZERO),DCMPLX(MDL_MT),1,I_SO,25)
C     CutTools call for loop numbers 37,41
      CALL LOOP_2(10,26,DCMPLX(ZERO),DCMPLX(MDL_MT),1,I_SO,26)
C     CutTools call for loop numbers 38
      CALL LOOP_3(5,10,11,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,27)
C     CutTools call for loop numbers 39
      CALL LOOP_3(5,10,11,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,28)
C     CutTools call for loop numbers 40,412,457,251,252,253
      CALL LOOP_2(12,27,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,29)
C     CutTools call for loop numbers 42
      CALL LOOP_3(4,10,12,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO),2
     $ ,I_SO,30)
C     CutTools call for loop numbers 43
      CALL LOOP_3(4,10,12,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT),2
     $ ,I_SO,31)
C     CutTools call for loop numbers 44,413,414,458,459,274,275,276
      CALL LOOP_3(2,5,27,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO),3,I_SO
     $ ,32)
C     CutTools call for loop numbers 45,46,47,271,272,273
      CALL LOOP_2_3(1,2,2,5,27,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,33)
C     CutTools call for loop numbers 48,114
      CALL LOOP_3(2,4,25,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,34)
C     CutTools call for loop numbers 49
      CALL LOOP_4(2,4,5,10,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,35)
C     CutTools call for loop numbers 50
      CALL LOOP_4(2,5,4,10,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,36)
C     CutTools call for loop numbers 51
      CALL LOOP_3(2,10,9,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,37)
C     CutTools call for loop numbers 52
      CALL LOOP_4(2,4,10,5,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),3,I_SO,38)
C     CutTools call for loop numbers 54
      CALL LOOP_2(10,30,DCMPLX(ZERO),DCMPLX(MDL_MT),1,I_SO,39)
C     CutTools call for loop numbers 55
      CALL LOOP_3(2,10,9,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,40)
C     CutTools call for loop numbers 57,140
      CALL LOOP_3(2,4,25,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,41)
C     CutTools call for loop numbers 58
      CALL LOOP_4(2,4,10,5,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),3,I_SO,42)
C     CutTools call for loop numbers 59
      CALL LOOP_4(2,5,4,10,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,43)
C     CutTools call for loop numbers 61
      CALL LOOP_4(2,4,5,10,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,44)
C     CutTools call for loop numbers
C      65,66,67,391,392,393,394,395,396,397,398,399
      CALL LOOP_2_3(2,1,5,2,27,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,45)
C     CutTools call for loop numbers 68,69,70
      CALL LOOP_3_4(1,1,2,4,10,5,2,DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),1,I_SO,46)
C     CutTools call for loop numbers 71,106
      CALL LOOP_2(14,31,DCMPLX(ZERO),DCMPLX(MDL_MT),1,I_SO,47)
C     CutTools call for loop numbers 72,76
      CALL LOOP_2(13,32,DCMPLX(ZERO),DCMPLX(MDL_MT),1,I_SO,48)
C     CutTools call for loop numbers 73
      CALL LOOP_3(5,14,13,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,49)
C     CutTools call for loop numbers 74
      CALL LOOP_3(5,14,13,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,50)
C     CutTools call for loop numbers 75,415,460
      CALL LOOP_2(12,33,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,51)
C     CutTools call for loop numbers 77
      CALL LOOP_3(3,13,12,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO),2
     $ ,I_SO,52)
C     CutTools call for loop numbers 78
      CALL LOOP_3(3,13,12,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT),2
     $ ,I_SO,53)
C     CutTools call for loop numbers 79,416,417,461,462
      CALL LOOP_3(2,5,33,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO),3,I_SO
     $ ,54)
C     CutTools call for loop numbers 80,81,82
      CALL LOOP_2_3(1,2,2,5,33,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,55)
C     CutTools call for loop numbers 83,115
      CALL LOOP_3(2,3,31,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,56)
C     CutTools call for loop numbers 84
      CALL LOOP_4(2,3,5,13,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,57)
C     CutTools call for loop numbers 85
      CALL LOOP_4(2,5,3,13,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,58)
C     CutTools call for loop numbers 86
      CALL LOOP_3(2,13,8,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,59)
C     CutTools call for loop numbers 87
      CALL LOOP_4(2,3,13,5,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),3,I_SO,60)
C     CutTools call for loop numbers 89
      CALL LOOP_2(13,36,DCMPLX(ZERO),DCMPLX(MDL_MT),1,I_SO,61)
C     CutTools call for loop numbers 90
      CALL LOOP_3(2,13,8,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,62)
C     CutTools call for loop numbers 92,132
      CALL LOOP_3(2,3,31,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,63)
C     CutTools call for loop numbers 93
      CALL LOOP_4(2,3,13,5,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),3,I_SO,64)
C     CutTools call for loop numbers 94
      CALL LOOP_4(2,5,3,13,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,65)
C     CutTools call for loop numbers 96
      CALL LOOP_4(2,3,5,13,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,66)
C     CutTools call for loop numbers 97,98,99
      CALL LOOP_2_3(1,2,5,2,33,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,67)
C     CutTools call for loop numbers 100,101,102
      CALL LOOP_2_3(2,1,5,2,33,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,68)
C     CutTools call for loop numbers 103,104,105
      CALL LOOP_3_4(1,1,2,3,13,5,2,DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),1,I_SO,69)
C     CutTools call for loop numbers 107,418,463,125,126,127
      CALL LOOP_2(15,38,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,70)
C     CutTools call for loop numbers 108
      CALL LOOP_3(4,14,15,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO),2
     $ ,I_SO,71)
C     CutTools call for loop numbers 109
      CALL LOOP_3(4,14,15,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT),2
     $ ,I_SO,72)
C     CutTools call for loop numbers 111,419,464
      CALL LOOP_2(15,40,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,73)
C     CutTools call for loop numbers 112
      CALL LOOP_3(3,11,15,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO),2
     $ ,I_SO,74)
C     CutTools call for loop numbers 113
      CALL LOOP_3(3,11,15,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT),2
     $ ,I_SO,75)
C     CutTools call for loop numbers 116
      CALL LOOP_4(2,3,4,15,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),3,I_SO,76)
C     CutTools call for loop numbers 117
      CALL LOOP_4(2,4,3,15,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),3,I_SO,77)
C     CutTools call for loop numbers 118,119,120
      CALL LOOP_2_3(1,2,2,7,15,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,78)
C     CutTools call for loop numbers 121,420,421,465,466
      CALL LOOP_3(2,7,15,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO),3,I_SO
     $ ,79)
C     CutTools call for loop numbers 122
      CALL LOOP_4(2,3,15,4,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,80)
C     CutTools call for loop numbers 123,422,467,128,129,130
      CALL LOOP_2(7,41,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,81)
C     CutTools call for loop numbers 124,423,468
      CALL LOOP_2(15,42,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,82)
C     CutTools call for loop numbers 131,134,135,136
      CALL LOOP_3(3,4,41,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO),2
     $ ,I_SO,83)
C     CutTools call for loop numbers 133
      CALL LOOP_4(2,3,15,4,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,84)
C     CutTools call for loop numbers 137
      CALL LOOP_3(3,4,41,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT),2
     $ ,I_SO,85)
C     CutTools call for loop numbers 138
      CALL LOOP_4(2,4,3,15,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),3,I_SO,86)
C     CutTools call for loop numbers 139
      CALL LOOP_4(2,3,4,15,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),3,I_SO,87)
C     CutTools call for loop numbers 141,161
      CALL LOOP_3(1,4,32,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,88)
C     CutTools call for loop numbers 142
      CALL LOOP_4(1,4,5,14,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,89)
C     CutTools call for loop numbers 143,424,425,469,470,188,189,190
      CALL LOOP_3(1,5,38,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO),3,I_SO
     $ ,90)
C     CutTools call for loop numbers 144
      CALL LOOP_4(1,5,4,14,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,91)
C     CutTools call for loop numbers 145,146,147,194,195,196
      CALL LOOP_2_3(1,2,1,5,38,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,92)
C     CutTools call for loop numbers 148
      CALL LOOP_3(1,14,9,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,93)
C     CutTools call for loop numbers 149
      CALL LOOP_4(1,4,14,5,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),3,I_SO,94)
C     CutTools call for loop numbers 150,159
      CALL LOOP_3(1,3,26,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,95)
C     CutTools call for loop numbers 151
      CALL LOOP_4(1,3,5,11,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,96)
C     CutTools call for loop numbers 152,426,427,471,472
      CALL LOOP_3(1,5,40,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO),3,I_SO
     $ ,97)
C     CutTools call for loop numbers 153
      CALL LOOP_4(1,5,3,11,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,98)
C     CutTools call for loop numbers 154,155,156
      CALL LOOP_2_3(1,2,1,5,40,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,99)
C     CutTools call for loop numbers 157
      CALL LOOP_3(1,11,8,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,100)
C     CutTools call for loop numbers 158
      CALL LOOP_4(1,3,11,5,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),3,I_SO,101)
C     CutTools call for loop numbers 160
      CALL LOOP_4(1,3,4,12,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),3,I_SO,102)
C     CutTools call for loop numbers 162
      CALL LOOP_4(1,4,3,12,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),3,I_SO,103)
C     CutTools call for loop numbers 163,164,165
      CALL LOOP_2_3(1,2,1,7,12,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,104)
C     CutTools call for loop numbers 166,428,429,473,474
      CALL LOOP_3(1,7,12,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO),3,I_SO
     $ ,105)
C     CutTools call for loop numbers 167
      CALL LOOP_4(1,3,12,4,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,106)
C     CutTools call for loop numbers 168
      CALL LOOP_5(1,2,4,5,3,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),DCMPLX(ZERO),4,I_SO,107)
C     CutTools call for loop numbers 169
      CALL LOOP_5(1,2,5,4,3,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),DCMPLX(ZERO),4,I_SO,108)
C     CutTools call for loop numbers 170
      CALL LOOP_5(1,2,3,5,4,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),DCMPLX(ZERO),4,I_SO,109)
C     CutTools call for loop numbers 171
      CALL LOOP_5(1,2,5,3,4,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),DCMPLX(ZERO),4,I_SO,110)
C     CutTools call for loop numbers 172
      CALL LOOP_5(1,2,3,4,5,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),DCMPLX(ZERO),4,I_SO,111)
C     CutTools call for loop numbers 173
      CALL LOOP_5(1,2,4,3,5,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),DCMPLX(ZERO),4,I_SO,112)
C     CutTools call for loop numbers 177,430,435,475,480
      CALL LOOP_3(1,2,20,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO),3,I_SO
     $ ,113)
C     CutTools call for loop numbers 178,431,436,476,481
      CALL LOOP_4(1,2,5,7,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),4,I_SO,114)
C     CutTools call for loop numbers 179,432,439,477,484
      CALL LOOP_4(1,2,7,5,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),4,I_SO,115)
C     CutTools call for loop numbers 180,433,442,478,487,174,175,176
      CALL LOOP_3(1,2,22,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO),3,I_SO
     $ ,116)
C     CutTools call for loop numbers 181
      CALL LOOP_4(1,2,4,8,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,117)
C     CutTools call for loop numbers 182
      CALL LOOP_4(1,2,8,4,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,118)
C     CutTools call for loop numbers 183,434,443,479,488
      CALL LOOP_3(1,2,24,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO),3,I_SO
     $ ,119)
C     CutTools call for loop numbers 184
      CALL LOOP_4(1,2,9,3,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,120)
C     CutTools call for loop numbers 185
      CALL LOOP_4(1,2,3,9,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,121)
C     CutTools call for loop numbers 186,437,441,482,486
      CALL LOOP_3(1,5,42,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO),3,I_SO
     $ ,122)
C     CutTools call for loop numbers 187,438,440,483,485
      CALL LOOP_4(1,5,2,7,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),4,I_SO,123)
C     CutTools call for loop numbers 191,192,193
      CALL LOOP_3_4(1,1,2,1,7,5,2,DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),2,I_SO,124)
C     CutTools call for loop numbers 197,198,199
      CALL LOOP_2_3(1,2,1,5,42,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,125)
C     CutTools call for loop numbers 200
      CALL LOOP_3(1,4,36,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,126)
C     CutTools call for loop numbers 201
      CALL LOOP_4(1,4,2,8,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,127)
C     CutTools call for loop numbers 202,203,204
      CALL LOOP_2_3(1,2,1,22,2,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,128)
C     CutTools call for loop numbers 205
      CALL LOOP_5(1,3,2,5,4,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),DCMPLX(ZERO),4,I_SO,129)
C     CutTools call for loop numbers 206
      CALL LOOP_5(1,3,5,2,4,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),DCMPLX(ZERO),4,I_SO,130)
C     CutTools call for loop numbers 207
      CALL LOOP_5(1,3,2,4,5,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),DCMPLX(ZERO),4,I_SO,131)
C     CutTools call for loop numbers 208
      CALL LOOP_5(1,3,4,2,5,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),DCMPLX(ZERO),4,I_SO,132)
C     CutTools call for loop numbers 209,210,211
      CALL LOOP_4_5(1,1,1,2,1,3,4,5,2,DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,133)
C     CutTools call for loop numbers 212
      CALL LOOP_3(1,3,30,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,134)
C     CutTools call for loop numbers 213
      CALL LOOP_4(1,3,2,9,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,135)
C     CutTools call for loop numbers 214,215,216
      CALL LOOP_2_3(1,2,1,24,2,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,136)
C     CutTools call for loop numbers 217
      CALL LOOP_5(1,4,2,3,5,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),DCMPLX(ZERO),4,I_SO,137)
C     CutTools call for loop numbers 218
      CALL LOOP_5(1,4,3,2,5,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),DCMPLX(ZERO),4,I_SO,138)
C     CutTools call for loop numbers 219,220,221
      CALL LOOP_4_5(1,1,1,2,1,4,3,5,2,DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,139)
C     CutTools call for loop numbers 222
      CALL LOOP_2(9,43,DCMPLX(ZERO),DCMPLX(MDL_MT),1,I_SO,140)
C     CutTools call for loop numbers 223
      CALL LOOP_2(14,44,DCMPLX(ZERO),DCMPLX(MDL_MT),1,I_SO,141)
C     CutTools call for loop numbers 224
      CALL LOOP_3(1,14,9,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,142)
C     CutTools call for loop numbers 225
      CALL LOOP_3(4,5,43,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,143)
C     CutTools call for loop numbers 226,263
      CALL LOOP_3(1,4,32,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,144)
C     CutTools call for loop numbers 227
      CALL LOOP_4(1,4,14,5,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),3,I_SO,145)
C     CutTools call for loop numbers 228
      CALL LOOP_4(1,5,4,14,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,146)
C     CutTools call for loop numbers 229
      CALL LOOP_3(4,5,43,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,147)
C     CutTools call for loop numbers 230
      CALL LOOP_4(1,4,5,14,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,148)
C     CutTools call for loop numbers 231,232,233
      CALL LOOP_2_3(1,2,5,38,1,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,149)
C     CutTools call for loop numbers 234
      CALL LOOP_2(8,45,DCMPLX(ZERO),DCMPLX(MDL_MT),1,I_SO,150)
C     CutTools call for loop numbers 235
      CALL LOOP_2(11,46,DCMPLX(ZERO),DCMPLX(MDL_MT),1,I_SO,151)
C     CutTools call for loop numbers 236
      CALL LOOP_3(1,11,8,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,152)
C     CutTools call for loop numbers 237
      CALL LOOP_3(3,5,45,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,153)
C     CutTools call for loop numbers 238,255
      CALL LOOP_3(1,3,26,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,154)
C     CutTools call for loop numbers 239
      CALL LOOP_4(1,3,11,5,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),3,I_SO,155)
C     CutTools call for loop numbers 240
      CALL LOOP_4(1,5,3,11,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,156)
C     CutTools call for loop numbers 241
      CALL LOOP_3(3,5,45,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,157)
C     CutTools call for loop numbers 242
      CALL LOOP_4(1,3,5,11,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,158)
C     CutTools call for loop numbers 243,244,245
      CALL LOOP_2_3(1,2,5,40,1,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,159)
C     CutTools call for loop numbers
C      246,328,329,330,444,448,449,450,489,493,494,495,7,8,9
      CALL LOOP_2(7,47,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,160)
C     CutTools call for loop numbers 247,445,490
      CALL LOOP_2(12,48,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,161)
C     CutTools call for loop numbers 248,249,250
      CALL LOOP_2_3(1,2,7,12,1,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,162)
C     CutTools call for loop numbers 254,331,332,333,21,22,23
      CALL LOOP_3(3,4,47,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO),2
     $ ,I_SO,163)
C     CutTools call for loop numbers 256
      CALL LOOP_4(1,3,12,4,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,164)
C     CutTools call for loop numbers 257,258,259
      CALL LOOP_3_4(1,1,2,3,4,12,1,DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),1,I_SO,165)
C     CutTools call for loop numbers 260,334,335,336
      CALL LOOP_3(3,4,47,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT),2
     $ ,I_SO,166)
C     CutTools call for loop numbers 261
      CALL LOOP_4(1,3,4,12,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),3,I_SO,167)
C     CutTools call for loop numbers 262
      CALL LOOP_4(1,4,3,12,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),3,I_SO,168)
C     CutTools call for loop numbers 264,446,447,491,492
      CALL LOOP_3(2,5,48,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO),3,I_SO
     $ ,169)
C     CutTools call for loop numbers 265,266,267
      CALL LOOP_2_3(1,2,2,5,48,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,170)
C     CutTools call for loop numbers 268,269,270
      CALL LOOP_3_4(1,1,2,2,7,5,1,DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),2,I_SO,171)
C     CutTools call for loop numbers 277
      CALL LOOP_3(2,4,46,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,172)
C     CutTools call for loop numbers 278
      CALL LOOP_4(1,4,2,8,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,173)
C     CutTools call for loop numbers 279,280,281
      CALL LOOP_2_3(1,2,2,22,1,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,174)
C     CutTools call for loop numbers 282
      CALL LOOP_5(1,3,2,4,5,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),DCMPLX(MDL_MT),4,I_SO,175)
C     CutTools call for loop numbers 283
      CALL LOOP_5(1,3,2,5,4,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),DCMPLX(MDL_MT),4,I_SO,176)
C     CutTools call for loop numbers 284
      CALL LOOP_5(1,4,2,3,5,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),DCMPLX(MDL_MT),4,I_SO,177)
C     CutTools call for loop numbers 285,286,287
      CALL LOOP_4_5(1,1,1,2,2,3,4,5,1,DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,178)
C     CutTools call for loop numbers 288
      CALL LOOP_4(1,3,2,9,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,179)
C     CutTools call for loop numbers 289
      CALL LOOP_3(2,3,44,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,180)
C     CutTools call for loop numbers 290,291,292
      CALL LOOP_2_3(1,2,2,24,1,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,181)
C     CutTools call for loop numbers 293
      CALL LOOP_5(1,3,5,2,4,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),DCMPLX(MDL_MT),4,I_SO,182)
C     CutTools call for loop numbers 294,295,296
      CALL LOOP_4_5(1,1,1,2,2,4,3,5,1,DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,183)
C     CutTools call for loop numbers 297,298,299
      CALL LOOP_2_3(1,2,5,2,48,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,184)
C     CutTools call for loop numbers 300,301,302
      CALL LOOP_3_4(1,2,1,5,2,1,7,DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),2,I_SO,185)
C     CutTools call for loop numbers 303,304,305
      CALL LOOP_2_3(1,2,5,42,1,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,186)
C     CutTools call for loop numbers 306
      CALL LOOP_3(2,4,46,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,187)
C     CutTools call for loop numbers 307
      CALL LOOP_3(1,4,36,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,188)
C     CutTools call for loop numbers 308
      CALL LOOP_4(1,2,8,4,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,189)
C     CutTools call for loop numbers 309
      CALL LOOP_4(1,2,4,8,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,190)
C     CutTools call for loop numbers 310
      CALL LOOP_3(1,3,30,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,191)
C     CutTools call for loop numbers 311
      CALL LOOP_4(1,2,9,3,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,192)
C     CutTools call for loop numbers 312
      CALL LOOP_3(2,3,44,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,193)
C     CutTools call for loop numbers 313
      CALL LOOP_4(1,2,3,9,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,194)
C     CutTools call for loop numbers 314
      CALL LOOP_5(1,3,4,2,5,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),DCMPLX(MDL_MT),4,I_SO,195)
C     CutTools call for loop numbers 315
      CALL LOOP_5(1,2,5,4,3,DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),4,I_SO,196)
C     CutTools call for loop numbers 316
      CALL LOOP_5(1,4,3,2,5,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),DCMPLX(MDL_MT),4,I_SO,197)
C     CutTools call for loop numbers 317
      CALL LOOP_5(1,2,5,3,4,DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),4,I_SO,198)
C     CutTools call for loop numbers 318
      CALL LOOP_5(1,2,3,4,5,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),DCMPLX(MDL_MT),4,I_SO,199)
C     CutTools call for loop numbers 319
      CALL LOOP_5(1,2,4,3,5,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),DCMPLX(MDL_MT),4,I_SO,200)
C     CutTools call for loop numbers 320
      CALL LOOP_5(1,2,4,5,3,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),DCMPLX(MDL_MT),4,I_SO,201)
C     CutTools call for loop numbers 321
      CALL LOOP_5(1,2,3,5,4,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),DCMPLX(MDL_MT),4,I_SO,202)
C     CutTools call for loop numbers 322,323,324
      CALL LOOP_4_5(1,1,2,1,3,4,2,1,5,DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,203)
C     CutTools call for loop numbers 325,326,327
      CALL LOOP_4_5(1,1,1,2,3,4,5,2,1,DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,204)
C     CutTools call for loop numbers 346,347,348
      CALL LOOP_2_3(2,1,2,1,20,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,205)
C     CutTools call for loop numbers 349,350,351
      CALL LOOP_3_4(1,2,1,4,2,1,8,DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),1,I_SO,206)
C     CutTools call for loop numbers
C      352,353,354,337,338,339,340,341,342,343,344,345
      CALL LOOP_2_3(2,1,2,1,22,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,207)
C     CutTools call for loop numbers 355,356,357
      CALL LOOP_4_5(1,1,1,2,3,5,4,2,1,DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,208)
C     CutTools call for loop numbers 358,359,360
      CALL LOOP_3_4(1,2,1,3,2,1,9,DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),1,I_SO,209)
C     CutTools call for loop numbers 361,362,363
      CALL LOOP_2_3(2,1,2,1,24,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,210)
C     CutTools call for loop numbers 364,365,366
      CALL LOOP_3_4(1,1,2,4,14,5,1,DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),1,I_SO,211)
C     CutTools call for loop numbers
C      367,368,369,376,377,378,379,380,381,382,383,384
      CALL LOOP_2_3(2,1,5,1,38,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,212)
C     CutTools call for loop numbers 370,371,372
      CALL LOOP_3_4(1,1,2,3,11,5,1,DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),1,I_SO,213)
C     CutTools call for loop numbers 373,374,375
      CALL LOOP_2_3(2,1,5,1,40,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,214)
C     CutTools call for loop numbers 385,386,387
      CALL LOOP_2_3(2,1,5,1,42,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,215)
C     CutTools call for loop numbers 388,389,390
      CALL LOOP_4_5(1,1,2,1,2,3,5,1,4,DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),DCMPLX(MDL_MT),2,I_SO,216)
C     CutTools call for loop numbers 400,401,402
      CALL LOOP_2_3(2,1,5,2,48,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,217)
C     CutTools call for loop numbers 403,404,405
      CALL LOOP_4_5(1,1,2,1,1,3,5,2,4,DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),DCMPLX(MDL_MT),2,I_SO,218)
C     CutTools call for loop numbers 496
      CALL LOOP_2(7,19,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,219)
C     CutTools call for loop numbers 497
      CALL LOOP_2(6,20,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,220)
C     CutTools call for loop numbers 498,499
      CALL LOOP_3(5,6,7,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB),3
     $ ,I_SO,221)
C     CutTools call for loop numbers 500
      CALL LOOP_2(6,22,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,222)
C     CutTools call for loop numbers 501
      CALL LOOP_2(6,24,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,223)
C     CutTools call for loop numbers 502
      CALL LOOP_2(12,27,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,224)
C     CutTools call for loop numbers 503,504
      CALL LOOP_3(2,5,27,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,3,I_SO,225)
C     CutTools call for loop numbers 505
      CALL LOOP_2(12,33,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,226)
C     CutTools call for loop numbers 506,507
      CALL LOOP_3(2,5,33,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,3,I_SO,227)
C     CutTools call for loop numbers 508
      CALL LOOP_2(15,38,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,228)
C     CutTools call for loop numbers 509
      CALL LOOP_2(15,40,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,229)
C     CutTools call for loop numbers 510,511
      CALL LOOP_3(2,7,15,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,3,I_SO,230)
C     CutTools call for loop numbers 512
      CALL LOOP_2(7,41,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,231)
C     CutTools call for loop numbers 513
      CALL LOOP_2(15,42,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,232)
C     CutTools call for loop numbers 514,515
      CALL LOOP_3(1,5,38,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,3,I_SO,233)
C     CutTools call for loop numbers 516,517
      CALL LOOP_3(1,5,40,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,3,I_SO,234)
C     CutTools call for loop numbers 518,519
      CALL LOOP_3(1,7,12,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,3,I_SO,235)
C     CutTools call for loop numbers 520,525
      CALL LOOP_3(1,2,20,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,3,I_SO,236)
C     CutTools call for loop numbers 521,526
      CALL LOOP_4(1,2,5,7,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,DCMPLX(MDL_MB),4,I_SO,237)
C     CutTools call for loop numbers 522,529
      CALL LOOP_4(1,2,7,5,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,DCMPLX(MDL_MB),4,I_SO,238)
C     CutTools call for loop numbers 523,532
      CALL LOOP_3(1,2,22,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,3,I_SO,239)
C     CutTools call for loop numbers 524,533
      CALL LOOP_3(1,2,24,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,3,I_SO,240)
C     CutTools call for loop numbers 527,531
      CALL LOOP_3(1,5,42,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,3,I_SO,241)
C     CutTools call for loop numbers 528,530
      CALL LOOP_4(1,5,2,7,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,DCMPLX(MDL_MB),4,I_SO,242)
C     CutTools call for loop numbers 534,538,539,540
      CALL LOOP_2(7,47,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,243)
C     CutTools call for loop numbers 535
      CALL LOOP_2(12,48,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,244)
C     CutTools call for loop numbers 536,537
      CALL LOOP_3(2,5,48,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,3,I_SO,245)
C     CutTools call for loop numbers 541
      CALL LOOP_2(7,19,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,246)
C     CutTools call for loop numbers 542
      CALL LOOP_2(6,20,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,247)
C     CutTools call for loop numbers 543,544
      CALL LOOP_3(5,6,7,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT),3
     $ ,I_SO,248)
C     CutTools call for loop numbers 545
      CALL LOOP_2(6,22,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,249)
C     CutTools call for loop numbers 546
      CALL LOOP_2(6,24,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,250)
C     CutTools call for loop numbers 547
      CALL LOOP_2(12,27,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,251)
C     CutTools call for loop numbers 548,549
      CALL LOOP_3(2,5,27,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,3,I_SO,252)
C     CutTools call for loop numbers 550
      CALL LOOP_2(12,33,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,253)
C     CutTools call for loop numbers 551,552
      CALL LOOP_3(2,5,33,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,3,I_SO,254)
C     CutTools call for loop numbers 553
      CALL LOOP_2(15,38,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,255)
C     CutTools call for loop numbers 554
      CALL LOOP_2(15,40,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,256)
C     CutTools call for loop numbers 555,556
      CALL LOOP_3(2,7,15,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,3,I_SO,257)
C     CutTools call for loop numbers 557
      CALL LOOP_2(7,41,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,258)
C     CutTools call for loop numbers 558
      CALL LOOP_2(15,42,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,259)
C     CutTools call for loop numbers 559,560
      CALL LOOP_3(1,5,38,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,3,I_SO,260)
C     CutTools call for loop numbers 561,562
      CALL LOOP_3(1,5,40,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,3,I_SO,261)
C     CutTools call for loop numbers 563,564
      CALL LOOP_3(1,7,12,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,3,I_SO,262)
C     CutTools call for loop numbers 565,570
      CALL LOOP_3(1,2,20,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,3,I_SO,263)
C     CutTools call for loop numbers 566,571
      CALL LOOP_4(1,2,5,7,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),4,I_SO,264)
C     CutTools call for loop numbers 567,574
      CALL LOOP_4(1,2,7,5,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),4,I_SO,265)
C     CutTools call for loop numbers 568,577
      CALL LOOP_3(1,2,22,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,3,I_SO,266)
C     CutTools call for loop numbers 569,578
      CALL LOOP_3(1,2,24,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,3,I_SO,267)
C     CutTools call for loop numbers 572,576
      CALL LOOP_3(1,5,42,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,3,I_SO,268)
C     CutTools call for loop numbers 573,575
      CALL LOOP_4(1,5,2,7,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),4,I_SO,269)
C     CutTools call for loop numbers 579,583,584,585
      CALL LOOP_2(7,47,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,270)
C     CutTools call for loop numbers 580
      CALL LOOP_2(12,48,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,271)
C     CutTools call for loop numbers 581,582
      CALL LOOP_3(2,5,48,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,3,I_SO,272)
C     At this point, all reductions needed for (QCD=8 QED=0), i.e. of
C      split order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 5000

      GOTO 1001
 5000 CONTINUE
      CTCALL_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

