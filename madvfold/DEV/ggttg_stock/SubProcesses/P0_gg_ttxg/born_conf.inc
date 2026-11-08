      INTEGER IFR
      INTEGER LMAXCONFIGSB_USED
      PARAMETER (LMAXCONFIGSB_USED=15)
      INTEGER MAX_BRANCHB_USED
      PARAMETER (MAX_BRANCHB_USED=3)
      INTEGER MAPCONFIG(0 : LMAXCONFIGSB_USED)
      INTEGER IFOREST(2, -MAX_BRANCHB_USED:-1, LMAXCONFIGSB_USED)
      INTEGER SPROP(-MAX_BRANCHB_USED:-1, LMAXCONFIGSB_USED)
      INTEGER TPRID(-MAX_BRANCHB_USED:-1, LMAXCONFIGSB_USED)
      LOGICAL GFORCEBW(-MAX_BRANCHB_USED : -1, LMAXCONFIGSB_USED)

C     Here are the congifurations
C     Process: g g > t t~ g [ all = QCD QED ] QCD^2=8 QED^2=0
C     Diagram 1, Amplitude 1
      DATA MAPCONFIG(   1)/   1/
      DATA (IFOREST(IFR, -1,   1),IFR=1,2)/  4,  3/
      DATA SPROP(  -1,   1)/      21/
      DATA (IFOREST(IFR, -2,   1),IFR=1,2)/  5, -1/
      DATA SPROP(  -2,   1)/      21/
C     Diagram 2, Amplitude 2
      DATA MAPCONFIG(   2)/   2/
      DATA (IFOREST(IFR, -1,   2),IFR=1,2)/  5,  3/
      DATA SPROP(  -1,   2)/       6/
      DATA (IFOREST(IFR, -2,   2),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,   2)/      21/
C     Diagram 3, Amplitude 3
      DATA MAPCONFIG(   3)/   3/
      DATA (IFOREST(IFR, -1,   3),IFR=1,2)/  5,  4/
      DATA SPROP(  -1,   3)/      -6/
      DATA (IFOREST(IFR, -2,   3),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,   3)/      21/
C     Diagram 4, Amplitude 4
      DATA MAPCONFIG(   4)/   4/
      DATA (IFOREST(IFR, -1,   4),IFR=1,2)/  1,  3/
      DATA TPRID(  -1,   4)/       6/
      DATA (IFOREST(IFR, -2,   4),IFR=1,2)/ -1,  5/
      DATA TPRID(  -2,   4)/       6/
      DATA (IFOREST(IFR, -3,   4),IFR=1,2)/ -2,  4/
C     Diagram 5, Amplitude 5
      DATA MAPCONFIG(   5)/   5/
      DATA (IFOREST(IFR, -1,   5),IFR=1,2)/  1,  3/
      DATA TPRID(  -1,   5)/       6/
      DATA (IFOREST(IFR, -2,   5),IFR=1,2)/ -1,  4/
      DATA TPRID(  -2,   5)/      21/
      DATA (IFOREST(IFR, -3,   5),IFR=1,2)/ -2,  5/
C     Diagram 6, Amplitude 6
      DATA MAPCONFIG(   6)/   6/
      DATA (IFOREST(IFR, -1,   6),IFR=1,2)/  5,  4/
      DATA SPROP(  -1,   6)/      -6/
      DATA (IFOREST(IFR, -2,   6),IFR=1,2)/  1,  3/
      DATA TPRID(  -2,   6)/       6/
      DATA (IFOREST(IFR, -3,   6),IFR=1,2)/ -2, -1/
C     Diagram 7, Amplitude 7
      DATA MAPCONFIG(   7)/   7/
      DATA (IFOREST(IFR, -1,   7),IFR=1,2)/  1,  4/
      DATA TPRID(  -1,   7)/       6/
      DATA (IFOREST(IFR, -2,   7),IFR=1,2)/ -1,  5/
      DATA TPRID(  -2,   7)/       6/
      DATA (IFOREST(IFR, -3,   7),IFR=1,2)/ -2,  3/
C     Diagram 8, Amplitude 8
      DATA MAPCONFIG(   8)/   8/
      DATA (IFOREST(IFR, -1,   8),IFR=1,2)/  1,  4/
      DATA TPRID(  -1,   8)/       6/
      DATA (IFOREST(IFR, -2,   8),IFR=1,2)/ -1,  3/
      DATA TPRID(  -2,   8)/      21/
      DATA (IFOREST(IFR, -3,   8),IFR=1,2)/ -2,  5/
C     Diagram 9, Amplitude 9
      DATA MAPCONFIG(   9)/   9/
      DATA (IFOREST(IFR, -1,   9),IFR=1,2)/  5,  3/
      DATA SPROP(  -1,   9)/       6/
      DATA (IFOREST(IFR, -2,   9),IFR=1,2)/  1,  4/
      DATA TPRID(  -2,   9)/       6/
      DATA (IFOREST(IFR, -3,   9),IFR=1,2)/ -2, -1/
C     Diagram 10, Amplitude 10
      DATA MAPCONFIG(  10)/  10/
      DATA (IFOREST(IFR, -1,  10),IFR=1,2)/  1,  5/
      DATA TPRID(  -1,  10)/      21/
      DATA (IFOREST(IFR, -2,  10),IFR=1,2)/ -1,  4/
      DATA TPRID(  -2,  10)/       6/
      DATA (IFOREST(IFR, -3,  10),IFR=1,2)/ -2,  3/
C     Diagram 11, Amplitude 11
      DATA MAPCONFIG(  11)/  11/
      DATA (IFOREST(IFR, -1,  11),IFR=1,2)/  1,  5/
      DATA TPRID(  -1,  11)/      21/
      DATA (IFOREST(IFR, -2,  11),IFR=1,2)/ -1,  3/
      DATA TPRID(  -2,  11)/       6/
      DATA (IFOREST(IFR, -3,  11),IFR=1,2)/ -2,  4/
C     Diagram 12, Amplitude 12
      DATA MAPCONFIG(  12)/  12/
      DATA (IFOREST(IFR, -1,  12),IFR=1,2)/  4,  3/
      DATA SPROP(  -1,  12)/      21/
      DATA (IFOREST(IFR, -2,  12),IFR=1,2)/  1,  5/
      DATA TPRID(  -2,  12)/      21/
      DATA (IFOREST(IFR, -3,  12),IFR=1,2)/ -2, -1/
C     Diagram 13, Amplitude 13
      DATA MAPCONFIG(  13)/  13/
      DATA (IFOREST(IFR, -1,  13),IFR=1,2)/  5,  4/
      DATA SPROP(  -1,  13)/      -6/
      DATA (IFOREST(IFR, -2,  13),IFR=1,2)/  1, -1/
      DATA TPRID(  -2,  13)/       6/
      DATA (IFOREST(IFR, -3,  13),IFR=1,2)/ -2,  3/
C     Diagram 14, Amplitude 14
      DATA MAPCONFIG(  14)/  14/
      DATA (IFOREST(IFR, -1,  14),IFR=1,2)/  5,  3/
      DATA SPROP(  -1,  14)/       6/
      DATA (IFOREST(IFR, -2,  14),IFR=1,2)/  1, -1/
      DATA TPRID(  -2,  14)/       6/
      DATA (IFOREST(IFR, -3,  14),IFR=1,2)/ -2,  4/
C     Diagram 15, Amplitude 15
      DATA MAPCONFIG(  15)/  15/
      DATA (IFOREST(IFR, -1,  15),IFR=1,2)/  4,  3/
      DATA SPROP(  -1,  15)/      21/
      DATA (IFOREST(IFR, -2,  15),IFR=1,2)/  1, -1/
      DATA TPRID(  -2,  15)/      21/
      DATA (IFOREST(IFR, -3,  15),IFR=1,2)/ -2,  5/
C     Number of configs
      DATA MAPCONFIG(0)/  15/

C     Here are the BWs
      DATA GFORCEBW(-1,1)/.FALSE./
      DATA GFORCEBW(-2,1)/.FALSE./
      DATA GFORCEBW(-1,2)/.FALSE./
      DATA GFORCEBW(-2,2)/.FALSE./
      DATA GFORCEBW(-1,3)/.FALSE./
      DATA GFORCEBW(-2,3)/.FALSE./
      DATA GFORCEBW(-1,6)/.FALSE./
      DATA GFORCEBW(-1,9)/.FALSE./
      DATA GFORCEBW(-1,12)/.FALSE./
      DATA GFORCEBW(-1,13)/.FALSE./
      DATA GFORCEBW(-1,14)/.FALSE./
      DATA GFORCEBW(-1,15)/.FALSE./
