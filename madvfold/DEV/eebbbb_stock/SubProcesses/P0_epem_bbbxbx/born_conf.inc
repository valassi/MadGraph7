      INTEGER IFR
      INTEGER LMAXCONFIGSB_USED
      PARAMETER (LMAXCONFIGSB_USED=16)
      INTEGER MAX_BRANCHB_USED
      PARAMETER (MAX_BRANCHB_USED=3)
      INTEGER MAPCONFIG(0 : LMAXCONFIGSB_USED)
      INTEGER IFOREST(2, -MAX_BRANCHB_USED:-1, LMAXCONFIGSB_USED)
      INTEGER SPROP(-MAX_BRANCHB_USED:-1, LMAXCONFIGSB_USED)
      INTEGER TPRID(-MAX_BRANCHB_USED:-1, LMAXCONFIGSB_USED)
      LOGICAL GFORCEBW(-MAX_BRANCHB_USED : -1, LMAXCONFIGSB_USED)

C     Here are the congifurations
C     Process: e+ e- > b b b~ b~ [ all = QCD QED ] QCD^2=6 QED^2=4
C     Diagram 1, Amplitude 1
      DATA MAPCONFIG(   1)/   1/
      DATA (IFOREST(IFR, -1,   1),IFR=1,2)/  5,  3/
      DATA SPROP(  -1,   1)/      21/
      DATA (IFOREST(IFR, -2,   1),IFR=1,2)/  6, -1/
      DATA SPROP(  -2,   1)/      -5/
      DATA (IFOREST(IFR, -3,   1),IFR=1,2)/  4, -2/
      DATA SPROP(  -3,   1)/      22/
C     Diagram 2, Amplitude 2
      DATA MAPCONFIG(   2)/   2/
      DATA (IFOREST(IFR, -1,   2),IFR=1,2)/  5,  3/
      DATA SPROP(  -1,   2)/      21/
      DATA (IFOREST(IFR, -2,   2),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,   2)/       5/
      DATA (IFOREST(IFR, -3,   2),IFR=1,2)/  6, -2/
      DATA SPROP(  -3,   2)/      22/
C     Diagram 3, Amplitude 3
      DATA MAPCONFIG(   3)/   3/
      DATA (IFOREST(IFR, -1,   3),IFR=1,2)/  5,  3/
      DATA SPROP(  -1,   3)/      21/
      DATA (IFOREST(IFR, -2,   3),IFR=1,2)/  6, -1/
      DATA SPROP(  -2,   3)/      -5/
      DATA (IFOREST(IFR, -3,   3),IFR=1,2)/  4, -2/
      DATA SPROP(  -3,   3)/      23/
C     Diagram 4, Amplitude 4
      DATA MAPCONFIG(   4)/   4/
      DATA (IFOREST(IFR, -1,   4),IFR=1,2)/  5,  3/
      DATA SPROP(  -1,   4)/      21/
      DATA (IFOREST(IFR, -2,   4),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,   4)/       5/
      DATA (IFOREST(IFR, -3,   4),IFR=1,2)/  6, -2/
      DATA SPROP(  -3,   4)/      23/
C     Diagram 5, Amplitude 5
      DATA MAPCONFIG(   5)/   5/
      DATA (IFOREST(IFR, -1,   5),IFR=1,2)/  6,  3/
      DATA SPROP(  -1,   5)/      21/
      DATA (IFOREST(IFR, -2,   5),IFR=1,2)/  5, -1/
      DATA SPROP(  -2,   5)/      -5/
      DATA (IFOREST(IFR, -3,   5),IFR=1,2)/  4, -2/
      DATA SPROP(  -3,   5)/      22/
C     Diagram 6, Amplitude 6
      DATA MAPCONFIG(   6)/   6/
      DATA (IFOREST(IFR, -1,   6),IFR=1,2)/  6,  3/
      DATA SPROP(  -1,   6)/      21/
      DATA (IFOREST(IFR, -2,   6),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,   6)/       5/
      DATA (IFOREST(IFR, -3,   6),IFR=1,2)/  5, -2/
      DATA SPROP(  -3,   6)/      22/
C     Diagram 7, Amplitude 7
      DATA MAPCONFIG(   7)/   7/
      DATA (IFOREST(IFR, -1,   7),IFR=1,2)/  6,  3/
      DATA SPROP(  -1,   7)/      21/
      DATA (IFOREST(IFR, -2,   7),IFR=1,2)/  5, -1/
      DATA SPROP(  -2,   7)/      -5/
      DATA (IFOREST(IFR, -3,   7),IFR=1,2)/  4, -2/
      DATA SPROP(  -3,   7)/      23/
C     Diagram 8, Amplitude 8
      DATA MAPCONFIG(   8)/   8/
      DATA (IFOREST(IFR, -1,   8),IFR=1,2)/  6,  3/
      DATA SPROP(  -1,   8)/      21/
      DATA (IFOREST(IFR, -2,   8),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,   8)/       5/
      DATA (IFOREST(IFR, -3,   8),IFR=1,2)/  5, -2/
      DATA SPROP(  -3,   8)/      23/
C     Diagram 9, Amplitude 9
      DATA MAPCONFIG(   9)/   9/
      DATA (IFOREST(IFR, -1,   9),IFR=1,2)/  5,  4/
      DATA SPROP(  -1,   9)/      21/
      DATA (IFOREST(IFR, -2,   9),IFR=1,2)/  6, -1/
      DATA SPROP(  -2,   9)/      -5/
      DATA (IFOREST(IFR, -3,   9),IFR=1,2)/ -2,  3/
      DATA SPROP(  -3,   9)/      22/
C     Diagram 10, Amplitude 10
      DATA MAPCONFIG(  10)/  10/
      DATA (IFOREST(IFR, -1,  10),IFR=1,2)/  5,  4/
      DATA SPROP(  -1,  10)/      21/
      DATA (IFOREST(IFR, -2,  10),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  10)/       5/
      DATA (IFOREST(IFR, -3,  10),IFR=1,2)/  6, -2/
      DATA SPROP(  -3,  10)/      22/
C     Diagram 11, Amplitude 11
      DATA MAPCONFIG(  11)/  11/
      DATA (IFOREST(IFR, -1,  11),IFR=1,2)/  5,  4/
      DATA SPROP(  -1,  11)/      21/
      DATA (IFOREST(IFR, -2,  11),IFR=1,2)/  6, -1/
      DATA SPROP(  -2,  11)/      -5/
      DATA (IFOREST(IFR, -3,  11),IFR=1,2)/ -2,  3/
      DATA SPROP(  -3,  11)/      23/
C     Diagram 12, Amplitude 12
      DATA MAPCONFIG(  12)/  12/
      DATA (IFOREST(IFR, -1,  12),IFR=1,2)/  5,  4/
      DATA SPROP(  -1,  12)/      21/
      DATA (IFOREST(IFR, -2,  12),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  12)/       5/
      DATA (IFOREST(IFR, -3,  12),IFR=1,2)/  6, -2/
      DATA SPROP(  -3,  12)/      23/
C     Diagram 13, Amplitude 13
      DATA MAPCONFIG(  13)/  13/
      DATA (IFOREST(IFR, -1,  13),IFR=1,2)/  6,  4/
      DATA SPROP(  -1,  13)/      21/
      DATA (IFOREST(IFR, -2,  13),IFR=1,2)/  5, -1/
      DATA SPROP(  -2,  13)/      -5/
      DATA (IFOREST(IFR, -3,  13),IFR=1,2)/ -2,  3/
      DATA SPROP(  -3,  13)/      22/
C     Diagram 14, Amplitude 14
      DATA MAPCONFIG(  14)/  14/
      DATA (IFOREST(IFR, -1,  14),IFR=1,2)/  6,  4/
      DATA SPROP(  -1,  14)/      21/
      DATA (IFOREST(IFR, -2,  14),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  14)/       5/
      DATA (IFOREST(IFR, -3,  14),IFR=1,2)/  5, -2/
      DATA SPROP(  -3,  14)/      22/
C     Diagram 15, Amplitude 15
      DATA MAPCONFIG(  15)/  15/
      DATA (IFOREST(IFR, -1,  15),IFR=1,2)/  6,  4/
      DATA SPROP(  -1,  15)/      21/
      DATA (IFOREST(IFR, -2,  15),IFR=1,2)/  5, -1/
      DATA SPROP(  -2,  15)/      -5/
      DATA (IFOREST(IFR, -3,  15),IFR=1,2)/ -2,  3/
      DATA SPROP(  -3,  15)/      23/
C     Diagram 16, Amplitude 16
      DATA MAPCONFIG(  16)/  16/
      DATA (IFOREST(IFR, -1,  16),IFR=1,2)/  6,  4/
      DATA SPROP(  -1,  16)/      21/
      DATA (IFOREST(IFR, -2,  16),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  16)/       5/
      DATA (IFOREST(IFR, -3,  16),IFR=1,2)/  5, -2/
      DATA SPROP(  -3,  16)/      23/
C     Number of configs
      DATA MAPCONFIG(0)/  16/

C     Here are the BWs
      DATA GFORCEBW(-1,1)/.FALSE./
      DATA GFORCEBW(-2,1)/.FALSE./
      DATA GFORCEBW(-3,1)/.FALSE./
      DATA GFORCEBW(-1,2)/.FALSE./
      DATA GFORCEBW(-2,2)/.FALSE./
      DATA GFORCEBW(-3,2)/.FALSE./
      DATA GFORCEBW(-1,3)/.FALSE./
      DATA GFORCEBW(-2,3)/.FALSE./
      DATA GFORCEBW(-3,3)/.FALSE./
      DATA GFORCEBW(-1,4)/.FALSE./
      DATA GFORCEBW(-2,4)/.FALSE./
      DATA GFORCEBW(-3,4)/.FALSE./
      DATA GFORCEBW(-1,5)/.FALSE./
      DATA GFORCEBW(-2,5)/.FALSE./
      DATA GFORCEBW(-3,5)/.FALSE./
      DATA GFORCEBW(-1,6)/.FALSE./
      DATA GFORCEBW(-2,6)/.FALSE./
      DATA GFORCEBW(-3,6)/.FALSE./
      DATA GFORCEBW(-1,7)/.FALSE./
      DATA GFORCEBW(-2,7)/.FALSE./
      DATA GFORCEBW(-3,7)/.FALSE./
      DATA GFORCEBW(-1,8)/.FALSE./
      DATA GFORCEBW(-2,8)/.FALSE./
      DATA GFORCEBW(-3,8)/.FALSE./
      DATA GFORCEBW(-1,9)/.FALSE./
      DATA GFORCEBW(-2,9)/.FALSE./
      DATA GFORCEBW(-3,9)/.FALSE./
      DATA GFORCEBW(-1,10)/.FALSE./
      DATA GFORCEBW(-2,10)/.FALSE./
      DATA GFORCEBW(-3,10)/.FALSE./
      DATA GFORCEBW(-1,11)/.FALSE./
      DATA GFORCEBW(-2,11)/.FALSE./
      DATA GFORCEBW(-3,11)/.FALSE./
      DATA GFORCEBW(-1,12)/.FALSE./
      DATA GFORCEBW(-2,12)/.FALSE./
      DATA GFORCEBW(-3,12)/.FALSE./
      DATA GFORCEBW(-1,13)/.FALSE./
      DATA GFORCEBW(-2,13)/.FALSE./
      DATA GFORCEBW(-3,13)/.FALSE./
      DATA GFORCEBW(-1,14)/.FALSE./
      DATA GFORCEBW(-2,14)/.FALSE./
      DATA GFORCEBW(-3,14)/.FALSE./
      DATA GFORCEBW(-1,15)/.FALSE./
      DATA GFORCEBW(-2,15)/.FALSE./
      DATA GFORCEBW(-3,15)/.FALSE./
      DATA GFORCEBW(-1,16)/.FALSE./
      DATA GFORCEBW(-2,16)/.FALSE./
      DATA GFORCEBW(-3,16)/.FALSE./
