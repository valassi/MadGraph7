      SUBROUTINE MP_HELAS_CALLS_UVCT_1(P,NHEL,H,IC)
C     
      USE POLYNOMIAL_CONSTANTS
      USE ALOHA_OBJECT
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
      REAL*16     ZERO
      PARAMETER (ZERO=0.0E0_16)
      COMPLEX*32     IZERO
      PARAMETER (IZERO=CMPLX(0.0E0_16,0.0E0_16,KIND=16))
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=2, NSQUAREDSO=1, NAMPSO=2)
C     
C     ARGUMENTS
C     
      REAL*16 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      INTEGER H
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
      INTEGER FLAVOR(NEXTERNAL)
      DATA FLAVOR /NEXTERNAL*1/
      COMPLEX*32 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
C     
C     GLOBAL VARIABLES
C     

      INCLUDE 'mp_coupl_same_name.inc'

      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/FILTERS/GOODAMP,GOODHEL

      INTEGER SQSO_TARGET
      COMMON/SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      COMPLEX*32 AMP(NBORNAMPS)
      COMMON/MP_AMPS/AMP
      TYPE(MP_ALOHA) W(NWAVEFUNCS)
      COMMON/MP_W/W

      COMPLEX*32 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,
     $ -1:NLOOPWAVEFUNCS)
      COMPLEX*32 PL(0:3,-1:NLOOPWAVEFUNCS)
      COMMON/MP_WL/WL,PL

      COMPLEX*32 AMPL(3,NCTAMPS)
      COMMON/MP_AMPL/AMPL

C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.MP_UVCT_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

C     Amplitude(s) for UVCT diagram with ID 281
      CALL MP_FFV1_0(W(6),W(9),W(8),GC_5,AMPL(2,509))
      AMPL(2,509)=AMPL(2,509)*(4.0D0*UVWFCT_B_0_1EPS)
C     Amplitude(s) for UVCT diagram with ID 282
      CALL MP_FFV1_0(W(6),W(9),W(8),GC_5,AMPL(1,510))
      AMPL(1,510)=AMPL(1,510)*(4.0D0*UVWFCT_B_0)
C     Amplitude(s) for UVCT diagram with ID 283
      CALL MP_FFV1_0(W(10),W(4),W(8),GC_5,AMPL(2,511))
      AMPL(2,511)=AMPL(2,511)*(4.0D0*UVWFCT_B_0_1EPS)
C     Amplitude(s) for UVCT diagram with ID 284
      CALL MP_FFV1_0(W(10),W(4),W(8),GC_5,AMPL(1,512))
      AMPL(1,512)=AMPL(1,512)*(4.0D0*UVWFCT_B_0)
C     Amplitude(s) for UVCT diagram with ID 285
      CALL MP_FFV1_0(W(6),W(12),W(8),GC_5,AMPL(2,513))
      AMPL(2,513)=AMPL(2,513)*(4.0D0*UVWFCT_B_0_1EPS)
C     Amplitude(s) for UVCT diagram with ID 286
      CALL MP_FFV1_0(W(6),W(12),W(8),GC_5,AMPL(1,514))
      AMPL(1,514)=AMPL(1,514)*(4.0D0*UVWFCT_B_0)
C     Amplitude(s) for UVCT diagram with ID 287
      CALL MP_FFV1_0(W(13),W(4),W(8),GC_5,AMPL(2,515))
      AMPL(2,515)=AMPL(2,515)*(4.0D0*UVWFCT_B_0_1EPS)
C     Amplitude(s) for UVCT diagram with ID 288
      CALL MP_FFV1_0(W(13),W(4),W(8),GC_5,AMPL(1,516))
      AMPL(1,516)=AMPL(1,516)*(4.0D0*UVWFCT_B_0)
C     Amplitude(s) for UVCT diagram with ID 289
      CALL MP_FFV1_0(W(5),W(9),W(14),GC_5,AMPL(2,517))
      AMPL(2,517)=AMPL(2,517)*(4.0D0*UVWFCT_B_0_1EPS)
C     Amplitude(s) for UVCT diagram with ID 290
      CALL MP_FFV1_0(W(5),W(9),W(14),GC_5,AMPL(1,518))
      AMPL(1,518)=AMPL(1,518)*(4.0D0*UVWFCT_B_0)
C     Amplitude(s) for UVCT diagram with ID 291
      CALL MP_FFV1_0(W(15),W(4),W(14),GC_5,AMPL(2,519))
      AMPL(2,519)=AMPL(2,519)*(4.0D0*UVWFCT_B_0_1EPS)
C     Amplitude(s) for UVCT diagram with ID 292
      CALL MP_FFV1_0(W(15),W(4),W(14),GC_5,AMPL(1,520))
      AMPL(1,520)=AMPL(1,520)*(4.0D0*UVWFCT_B_0)
C     Amplitude(s) for UVCT diagram with ID 293
      CALL MP_FFV1_0(W(5),W(12),W(14),GC_5,AMPL(2,521))
      AMPL(2,521)=AMPL(2,521)*(4.0D0*UVWFCT_B_0_1EPS)
C     Amplitude(s) for UVCT diagram with ID 294
      CALL MP_FFV1_0(W(5),W(12),W(14),GC_5,AMPL(1,522))
      AMPL(1,522)=AMPL(1,522)*(4.0D0*UVWFCT_B_0)
C     Amplitude(s) for UVCT diagram with ID 295
      CALL MP_FFV1_0(W(16),W(4),W(14),GC_5,AMPL(2,523))
      AMPL(2,523)=AMPL(2,523)*(4.0D0*UVWFCT_B_0_1EPS)
C     Amplitude(s) for UVCT diagram with ID 296
      CALL MP_FFV1_0(W(16),W(4),W(14),GC_5,AMPL(1,524))
      AMPL(1,524)=AMPL(1,524)*(4.0D0*UVWFCT_B_0)
C     Amplitude(s) for UVCT diagram with ID 297
      CALL MP_FFV1_0(W(6),W(18),W(17),GC_5,AMPL(2,525))
      AMPL(2,525)=AMPL(2,525)*(4.0D0*UVWFCT_B_0_1EPS)
C     Amplitude(s) for UVCT diagram with ID 298
      CALL MP_FFV1_0(W(6),W(18),W(17),GC_5,AMPL(1,526))
      AMPL(1,526)=AMPL(1,526)*(4.0D0*UVWFCT_B_0)
C     Amplitude(s) for UVCT diagram with ID 299
      CALL MP_FFV1_0(W(10),W(3),W(17),GC_5,AMPL(2,527))
      AMPL(2,527)=AMPL(2,527)*(4.0D0*UVWFCT_B_0_1EPS)
C     Amplitude(s) for UVCT diagram with ID 300
      CALL MP_FFV1_0(W(10),W(3),W(17),GC_5,AMPL(1,528))
      AMPL(1,528)=AMPL(1,528)*(4.0D0*UVWFCT_B_0)
C     Amplitude(s) for UVCT diagram with ID 301
      CALL MP_FFV1_0(W(6),W(19),W(17),GC_5,AMPL(2,529))
      AMPL(2,529)=AMPL(2,529)*(4.0D0*UVWFCT_B_0_1EPS)
C     Amplitude(s) for UVCT diagram with ID 302
      CALL MP_FFV1_0(W(6),W(19),W(17),GC_5,AMPL(1,530))
      AMPL(1,530)=AMPL(1,530)*(4.0D0*UVWFCT_B_0)
C     Amplitude(s) for UVCT diagram with ID 303
      CALL MP_FFV1_0(W(13),W(3),W(17),GC_5,AMPL(2,531))
      AMPL(2,531)=AMPL(2,531)*(4.0D0*UVWFCT_B_0_1EPS)
C     Amplitude(s) for UVCT diagram with ID 304
      CALL MP_FFV1_0(W(13),W(3),W(17),GC_5,AMPL(1,532))
      AMPL(1,532)=AMPL(1,532)*(4.0D0*UVWFCT_B_0)
C     Amplitude(s) for UVCT diagram with ID 305
      CALL MP_FFV1_0(W(5),W(18),W(20),GC_5,AMPL(2,533))
      AMPL(2,533)=AMPL(2,533)*(4.0D0*UVWFCT_B_0_1EPS)
C     Amplitude(s) for UVCT diagram with ID 306
      CALL MP_FFV1_0(W(5),W(18),W(20),GC_5,AMPL(1,534))
      AMPL(1,534)=AMPL(1,534)*(4.0D0*UVWFCT_B_0)
C     Amplitude(s) for UVCT diagram with ID 307
      CALL MP_FFV1_0(W(15),W(3),W(20),GC_5,AMPL(2,535))
      AMPL(2,535)=AMPL(2,535)*(4.0D0*UVWFCT_B_0_1EPS)
C     Amplitude(s) for UVCT diagram with ID 308
      CALL MP_FFV1_0(W(15),W(3),W(20),GC_5,AMPL(1,536))
      AMPL(1,536)=AMPL(1,536)*(4.0D0*UVWFCT_B_0)
C     Amplitude(s) for UVCT diagram with ID 309
      CALL MP_FFV1_0(W(5),W(19),W(20),GC_5,AMPL(2,537))
      AMPL(2,537)=AMPL(2,537)*(4.0D0*UVWFCT_B_0_1EPS)
C     Amplitude(s) for UVCT diagram with ID 310
      CALL MP_FFV1_0(W(5),W(19),W(20),GC_5,AMPL(1,538))
      AMPL(1,538)=AMPL(1,538)*(4.0D0*UVWFCT_B_0)
C     Amplitude(s) for UVCT diagram with ID 311
      CALL MP_FFV1_0(W(16),W(3),W(20),GC_5,AMPL(2,539))
      AMPL(2,539)=AMPL(2,539)*(4.0D0*UVWFCT_B_0_1EPS)
C     Amplitude(s) for UVCT diagram with ID 312
      CALL MP_FFV1_0(W(16),W(3),W(20),GC_5,AMPL(1,540))
      AMPL(1,540)=AMPL(1,540)*(4.0D0*UVWFCT_B_0)
C     At this point, all UVCT amps needed for (QCD=6 QED=4), i.e. of
C      split order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 3000

      GOTO 1001
 3000 CONTINUE
      MP_UVCT_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

