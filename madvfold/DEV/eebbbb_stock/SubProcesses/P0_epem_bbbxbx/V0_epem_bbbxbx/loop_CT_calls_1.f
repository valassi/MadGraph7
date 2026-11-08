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
      IF (FILTER_SO.AND.CTCALL_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

C     CutTools call for loop numbers 1,7
      CALL LOOP_2(9,21,DCMPLX(ZERO),DCMPLX(MDL_MB),1,I_SO,1)
C     CutTools call for loop numbers 2,8,185,187,201,203
      CALL LOOP_2(8,22,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,2)
C     CutTools call for loop numbers 3,9
      CALL LOOP_3(6,8,9,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MB),2
     $ ,I_SO,3)
C     CutTools call for loop numbers 4,10
      CALL LOOP_2(23,10,DCMPLX(ZERO),DCMPLX(MDL_MB),1,I_SO,4)
C     CutTools call for loop numbers 5,11,186,188,202,204
      CALL LOOP_2(8,24,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,5)
C     CutTools call for loop numbers 6,12
      CALL LOOP_3(4,8,10,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MB),2
     $ ,I_SO,6)
C     CutTools call for loop numbers 13,16
      CALL LOOP_3(4,7,21,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(ZERO),2
     $ ,I_SO,7)
C     CutTools call for loop numbers 14,17
      CALL LOOP_4(4,7,6,8,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),3,I_SO,8)
C     CutTools call for loop numbers 15,18
      CALL LOOP_3(4,8,10,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(ZERO),2
     $ ,I_SO,9)
C     CutTools call for loop numbers 19,21
      CALL LOOP_4(4,6,8,7,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,DCMPLX(MDL_MB),3,I_SO,10)
C     CutTools call for loop numbers 20,22
      CALL LOOP_4(4,6,7,8,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,DCMPLX(MDL_MB),3,I_SO,11)
C     CutTools call for loop numbers 23,25
      CALL LOOP_3(6,8,9,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(ZERO),2
     $ ,I_SO,12)
C     CutTools call for loop numbers 24,26
      CALL LOOP_3(6,7,23,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(ZERO),2
     $ ,I_SO,13)
C     CutTools call for loop numbers 27,33
      CALL LOOP_2(9,27,DCMPLX(ZERO),DCMPLX(MDL_MB),1,I_SO,14)
C     CutTools call for loop numbers 28,34,189,191,207,209
      CALL LOOP_2(14,28,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,15)
C     CutTools call for loop numbers 29,35
      CALL LOOP_3(5,14,9,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MB),2
     $ ,I_SO,16)
C     CutTools call for loop numbers 30,36
      CALL LOOP_2(15,29,DCMPLX(ZERO),DCMPLX(MDL_MB),1,I_SO,17)
C     CutTools call for loop numbers 31,37,190,192,208,210
      CALL LOOP_2(14,30,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,18)
C     CutTools call for loop numbers 32,38
      CALL LOOP_3(4,14,15,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MB),2
     $ ,I_SO,19)
C     CutTools call for loop numbers 39,42
      CALL LOOP_3(4,7,27,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(ZERO),2
     $ ,I_SO,20)
C     CutTools call for loop numbers 40,43
      CALL LOOP_4(4,7,5,14,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),3,I_SO,21)
C     CutTools call for loop numbers 41,44
      CALL LOOP_3(4,14,15,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(ZERO),2
     $ ,I_SO,22)
C     CutTools call for loop numbers 45,47
      CALL LOOP_4(4,5,14,7,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,DCMPLX(MDL_MB),3,I_SO,23)
C     CutTools call for loop numbers 46,48
      CALL LOOP_4(4,5,7,14,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,DCMPLX(MDL_MB),3,I_SO,24)
C     CutTools call for loop numbers 49,51
      CALL LOOP_3(5,14,9,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(ZERO),2
     $ ,I_SO,25)
C     CutTools call for loop numbers 50,52
      CALL LOOP_3(5,7,29,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(ZERO),2
     $ ,I_SO,26)
C     CutTools call for loop numbers 53,61
      CALL LOOP_3(3,5,22,DCMPLX(MDL_MB),DCMPLX(ZERO),DCMPLX(ZERO),2
     $ ,I_SO,27)
C     CutTools call for loop numbers 54,62
      CALL LOOP_4(3,5,6,9,DCMPLX(MDL_MB),DCMPLX(ZERO),DCMPLX(MDL_MB)
     $ ,DCMPLX(ZERO),2,I_SO,28)
C     CutTools call for loop numbers 55,63
      CALL LOOP_3(3,6,28,DCMPLX(MDL_MB),DCMPLX(ZERO),DCMPLX(ZERO),2
     $ ,I_SO,29)
C     CutTools call for loop numbers 56,64
      CALL LOOP_4(3,6,5,9,DCMPLX(MDL_MB),DCMPLX(ZERO),DCMPLX(MDL_MB)
     $ ,DCMPLX(ZERO),2,I_SO,30)
C     CutTools call for loop numbers 57,65
      CALL LOOP_3(3,6,30,DCMPLX(MDL_MB),DCMPLX(ZERO),DCMPLX(ZERO),2
     $ ,I_SO,31)
C     CutTools call for loop numbers 58,66
      CALL LOOP_4(3,6,4,15,DCMPLX(MDL_MB),DCMPLX(ZERO),DCMPLX(MDL_MB)
     $ ,DCMPLX(ZERO),2,I_SO,32)
C     CutTools call for loop numbers 59,67
      CALL LOOP_3(3,5,24,DCMPLX(MDL_MB),DCMPLX(ZERO),DCMPLX(ZERO),2
     $ ,I_SO,33)
C     CutTools call for loop numbers 60,68
      CALL LOOP_4(3,5,4,10,DCMPLX(MDL_MB),DCMPLX(ZERO),DCMPLX(MDL_MB)
     $ ,DCMPLX(ZERO),2,I_SO,34)
C     CutTools call for loop numbers 69,73
      CALL LOOP_3(3,7,33,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(ZERO),2
     $ ,I_SO,35)
C     CutTools call for loop numbers 70,74
      CALL LOOP_4(3,7,6,17,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),3,I_SO,36)
C     CutTools call for loop numbers 71,75
      CALL LOOP_3(3,17,10,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(ZERO),2
     $ ,I_SO,37)
C     CutTools call for loop numbers 72,76
      CALL LOOP_3(3,17,10,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MB),2
     $ ,I_SO,38)
C     CutTools call for loop numbers 77,79
      CALL LOOP_4(3,6,17,7,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,DCMPLX(MDL_MB),3,I_SO,39)
C     CutTools call for loop numbers 78,80
      CALL LOOP_4(3,6,7,17,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,DCMPLX(MDL_MB),3,I_SO,40)
C     CutTools call for loop numbers 81,85
      CALL LOOP_3(3,7,34,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(ZERO),2
     $ ,I_SO,41)
C     CutTools call for loop numbers 82,86
      CALL LOOP_4(3,7,5,20,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),3,I_SO,42)
C     CutTools call for loop numbers 83,87
      CALL LOOP_3(3,20,15,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(ZERO),2
     $ ,I_SO,43)
C     CutTools call for loop numbers 84,88
      CALL LOOP_3(3,20,15,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MB),2
     $ ,I_SO,44)
C     CutTools call for loop numbers 89,91
      CALL LOOP_4(3,5,20,7,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,DCMPLX(MDL_MB),3,I_SO,45)
C     CutTools call for loop numbers 90,92
      CALL LOOP_4(3,5,7,20,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,DCMPLX(MDL_MB),3,I_SO,46)
C     CutTools call for loop numbers 93,101
      CALL LOOP_5(3,4,5,6,7,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MB),DCMPLX(MDL_MB),3,I_SO,47)
C     CutTools call for loop numbers 94,102
      CALL LOOP_5(3,4,6,5,7,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MB),DCMPLX(MDL_MB),3,I_SO,48)
C     CutTools call for loop numbers 95,103
      CALL LOOP_5(3,4,7,6,5,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,DCMPLX(ZERO),DCMPLX(MDL_MB),3,I_SO,49)
C     CutTools call for loop numbers 96,104
      CALL LOOP_5(3,4,7,5,6,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,DCMPLX(ZERO),DCMPLX(MDL_MB),3,I_SO,50)
C     CutTools call for loop numbers 97,105
      CALL LOOP_4(3,4,15,6,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MB),2,I_SO,51)
C     CutTools call for loop numbers 98,106
      CALL LOOP_4(3,4,6,15,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MB),2,I_SO,52)
C     CutTools call for loop numbers 99,107
      CALL LOOP_4(3,4,10,5,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MB),2,I_SO,53)
C     CutTools call for loop numbers 100,108
      CALL LOOP_4(3,4,5,10,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MB),2,I_SO,54)
C     CutTools call for loop numbers 109,115
      CALL LOOP_5(3,5,4,6,7,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MB),DCMPLX(MDL_MB),3,I_SO,55)
C     CutTools call for loop numbers 110,116
      CALL LOOP_3(3,5,22,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(MDL_MB),2
     $ ,I_SO,56)
C     CutTools call for loop numbers 111,117
      CALL LOOP_4(3,5,9,6,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MB),2,I_SO,57)
C     CutTools call for loop numbers 112,118
      CALL LOOP_5(3,5,7,4,6,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,DCMPLX(ZERO),DCMPLX(MDL_MB),3,I_SO,58)
C     CutTools call for loop numbers 113,119
      CALL LOOP_3(3,5,24,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(MDL_MB),2
     $ ,I_SO,59)
C     CutTools call for loop numbers 114,120
      CALL LOOP_4(3,5,4,10,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MB),2,I_SO,60)
C     CutTools call for loop numbers 121,127
      CALL LOOP_5(3,6,4,5,7,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MB),DCMPLX(MDL_MB),3,I_SO,61)
C     CutTools call for loop numbers 122,128
      CALL LOOP_4(3,5,9,6,DCMPLX(MDL_MB),DCMPLX(ZERO),DCMPLX(MDL_MB)
     $ ,DCMPLX(ZERO),2,I_SO,62)
C     CutTools call for loop numbers 123,129
      CALL LOOP_3(3,6,28,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(MDL_MB),2
     $ ,I_SO,63)
C     CutTools call for loop numbers 124,130
      CALL LOOP_3(3,6,30,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(MDL_MB),2
     $ ,I_SO,64)
C     CutTools call for loop numbers 125,131
      CALL LOOP_4(3,6,4,15,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MB),2,I_SO,65)
C     CutTools call for loop numbers 126,132
      CALL LOOP_5(3,5,4,7,6,DCMPLX(MDL_MB),DCMPLX(ZERO),DCMPLX(MDL_MB)
     $ ,DCMPLX(MDL_MB),DCMPLX(ZERO),3,I_SO,66)
C     CutTools call for loop numbers 133,138
      CALL LOOP_2(18,33,DCMPLX(ZERO),DCMPLX(MDL_MB),1,I_SO,67)
C     CutTools call for loop numbers 134,139,193,195,213,215
      CALL LOOP_2(17,35,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,68)
C     CutTools call for loop numbers 135,140
      CALL LOOP_3(6,17,18,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MB),2
     $ ,I_SO,69)
C     CutTools call for loop numbers 136,141
      CALL LOOP_2(36,10,DCMPLX(ZERO),DCMPLX(MDL_MB),1,I_SO,70)
C     CutTools call for loop numbers 137,142,194,196,214,216
      CALL LOOP_2(17,37,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,71)
C     CutTools call for loop numbers 143,145
      CALL LOOP_3(6,17,18,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(ZERO),2
     $ ,I_SO,72)
C     CutTools call for loop numbers 144,146
      CALL LOOP_3(6,7,36,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(ZERO),2
     $ ,I_SO,73)
C     CutTools call for loop numbers 147,152
      CALL LOOP_2(18,34,DCMPLX(ZERO),DCMPLX(MDL_MB),1,I_SO,74)
C     CutTools call for loop numbers 148,153,197,199,217,219
      CALL LOOP_2(20,40,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,75)
C     CutTools call for loop numbers 149,154
      CALL LOOP_3(5,20,18,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MB),2
     $ ,I_SO,76)
C     CutTools call for loop numbers 150,155
      CALL LOOP_2(15,41,DCMPLX(ZERO),DCMPLX(MDL_MB),1,I_SO,77)
C     CutTools call for loop numbers 151,156,198,200,218,220
      CALL LOOP_2(20,42,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,78)
C     CutTools call for loop numbers 157,159
      CALL LOOP_3(5,20,18,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(ZERO),2
     $ ,I_SO,79)
C     CutTools call for loop numbers 158,160
      CALL LOOP_3(5,7,41,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(ZERO),2
     $ ,I_SO,80)
C     CutTools call for loop numbers 161,167
      CALL LOOP_3(4,5,35,DCMPLX(MDL_MB),DCMPLX(ZERO),DCMPLX(ZERO),2
     $ ,I_SO,81)
C     CutTools call for loop numbers 162,168
      CALL LOOP_4(4,5,6,18,DCMPLX(MDL_MB),DCMPLX(ZERO),DCMPLX(MDL_MB)
     $ ,DCMPLX(ZERO),2,I_SO,82)
C     CutTools call for loop numbers 163,169
      CALL LOOP_3(4,6,40,DCMPLX(MDL_MB),DCMPLX(ZERO),DCMPLX(ZERO),2
     $ ,I_SO,83)
C     CutTools call for loop numbers 164,170
      CALL LOOP_4(4,6,5,18,DCMPLX(MDL_MB),DCMPLX(ZERO),DCMPLX(MDL_MB)
     $ ,DCMPLX(ZERO),2,I_SO,84)
C     CutTools call for loop numbers 165,171
      CALL LOOP_3(4,6,42,DCMPLX(MDL_MB),DCMPLX(ZERO),DCMPLX(ZERO),2
     $ ,I_SO,85)
C     CutTools call for loop numbers 166,172
      CALL LOOP_3(4,5,37,DCMPLX(MDL_MB),DCMPLX(ZERO),DCMPLX(ZERO),2
     $ ,I_SO,86)
C     CutTools call for loop numbers 173,176
      CALL LOOP_3(4,5,35,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(MDL_MB),2
     $ ,I_SO,87)
C     CutTools call for loop numbers 174,177
      CALL LOOP_4(4,5,18,6,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MB),2,I_SO,88)
C     CutTools call for loop numbers 175,178
      CALL LOOP_3(4,5,37,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(MDL_MB),2
     $ ,I_SO,89)
C     CutTools call for loop numbers 179,182
      CALL LOOP_4(4,5,18,6,DCMPLX(MDL_MB),DCMPLX(ZERO),DCMPLX(MDL_MB)
     $ ,DCMPLX(ZERO),2,I_SO,90)
C     CutTools call for loop numbers 180,183
      CALL LOOP_3(4,6,40,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(MDL_MB),2
     $ ,I_SO,91)
C     CutTools call for loop numbers 181,184
      CALL LOOP_3(4,6,42,DCMPLX(ZERO),DCMPLX(MDL_MB),DCMPLX(MDL_MB),2
     $ ,I_SO,92)
C     CutTools call for loop numbers 205,206,221,222
      CALL LOOP_3(11,8,20,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO),3
     $ ,I_SO,93)
C     CutTools call for loop numbers 211,212,223,224
      CALL LOOP_3(11,17,14,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO),3
     $ ,I_SO,94)
C     CutTools call for loop numbers 225,227
      CALL LOOP_2(8,22,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,95)
C     CutTools call for loop numbers 226,228
      CALL LOOP_2(8,24,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,96)
C     CutTools call for loop numbers 229,230
      CALL LOOP_3(11,8,20,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,3,I_SO,97)
C     CutTools call for loop numbers 231,233
      CALL LOOP_2(14,28,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,98)
C     CutTools call for loop numbers 232,234
      CALL LOOP_2(14,30,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,99)
C     CutTools call for loop numbers 235,236
      CALL LOOP_3(11,17,14,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,3,I_SO,100)
C     CutTools call for loop numbers 237,239
      CALL LOOP_2(17,35,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,101)
C     CutTools call for loop numbers 238,240
      CALL LOOP_2(17,37,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,102)
C     CutTools call for loop numbers 241,243
      CALL LOOP_2(20,40,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,103)
C     CutTools call for loop numbers 242,244
      CALL LOOP_2(20,42,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,104)
C     CutTools call for loop numbers 245,247
      CALL LOOP_2(8,22,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,105)
C     CutTools call for loop numbers 246,248
      CALL LOOP_2(8,24,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,106)
C     CutTools call for loop numbers 249,250
      CALL LOOP_3(11,8,20,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,3,I_SO,107)
C     CutTools call for loop numbers 251,253
      CALL LOOP_2(14,28,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,108)
C     CutTools call for loop numbers 252,254
      CALL LOOP_2(14,30,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,109)
C     CutTools call for loop numbers 255,256
      CALL LOOP_3(11,17,14,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,3,I_SO,110)
C     CutTools call for loop numbers 257,259
      CALL LOOP_2(17,35,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,111)
C     CutTools call for loop numbers 258,260
      CALL LOOP_2(17,37,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,112)
C     CutTools call for loop numbers 261,263
      CALL LOOP_2(20,40,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,113)
C     CutTools call for loop numbers 262,264
      CALL LOOP_2(20,42,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,114)
C     At this point, all reductions needed for (QCD=6 QED=4), i.e. of
C      split order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 5000

      GOTO 1001
 5000 CONTINUE
      CTCALL_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

