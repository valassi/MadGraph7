      INTEGER IFR
      INTEGER LMAXCONFIGSB_USED
      PARAMETER (LMAXCONFIGSB_USED=58)
      INTEGER MAX_BRANCHB_USED
      PARAMETER (MAX_BRANCHB_USED=5)
      INTEGER MAPCONFIG(0 : LMAXCONFIGSB_USED)
      INTEGER IFOREST(2, -MAX_BRANCHB_USED:-1, LMAXCONFIGSB_USED)
      INTEGER SPROP(-MAX_BRANCHB_USED:-1, LMAXCONFIGSB_USED)
      INTEGER TPRID(-MAX_BRANCHB_USED:-1, LMAXCONFIGSB_USED)
      LOGICAL GFORCEBW(-MAX_BRANCHB_USED : -1, LMAXCONFIGSB_USED)

C     Here are the congifurations
C     Process: e+ e- > b b~ mu- mu+ g [ all = QCD QED ] QCD^2=4 QED^2=8
C     Diagram 1, Amplitude 1
      DATA MAPCONFIG(   1)/   1/
      DATA (IFOREST(IFR, -1,   1),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,   1)/       5/
      DATA (IFOREST(IFR, -2,   1),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,   1)/      22/
      DATA (IFOREST(IFR, -3,   1),IFR=1,2)/  6, -2/
      DATA SPROP(  -3,   1)/     -13/
      DATA (IFOREST(IFR, -4,   1),IFR=1,2)/  5, -3/
      DATA SPROP(  -4,   1)/      22/
C     Diagram 2, Amplitude 2
      DATA MAPCONFIG(   2)/   2/
      DATA (IFOREST(IFR, -1,   2),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,   2)/       5/
      DATA (IFOREST(IFR, -2,   2),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,   2)/      23/
      DATA (IFOREST(IFR, -3,   2),IFR=1,2)/  6, -2/
      DATA SPROP(  -3,   2)/     -13/
      DATA (IFOREST(IFR, -4,   2),IFR=1,2)/  5, -3/
      DATA SPROP(  -4,   2)/      22/
C     Diagram 3, Amplitude 3
      DATA MAPCONFIG(   3)/   3/
      DATA (IFOREST(IFR, -1,   3),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,   3)/       5/
      DATA (IFOREST(IFR, -2,   3),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,   3)/      22/
      DATA (IFOREST(IFR, -3,   3),IFR=1,2)/  5, -2/
      DATA SPROP(  -3,   3)/      13/
      DATA (IFOREST(IFR, -4,   3),IFR=1,2)/  6, -3/
      DATA SPROP(  -4,   3)/      22/
C     Diagram 4, Amplitude 4
      DATA MAPCONFIG(   4)/   4/
      DATA (IFOREST(IFR, -1,   4),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,   4)/       5/
      DATA (IFOREST(IFR, -2,   4),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,   4)/      23/
      DATA (IFOREST(IFR, -3,   4),IFR=1,2)/  5, -2/
      DATA SPROP(  -3,   4)/      13/
      DATA (IFOREST(IFR, -4,   4),IFR=1,2)/  6, -3/
      DATA SPROP(  -4,   4)/      22/
C     Diagram 5, Amplitude 5
      DATA MAPCONFIG(   5)/   5/
      DATA (IFOREST(IFR, -1,   5),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,   5)/       5/
      DATA (IFOREST(IFR, -2,   5),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,   5)/      22/
      DATA (IFOREST(IFR, -3,   5),IFR=1,2)/  6, -2/
      DATA SPROP(  -3,   5)/     -13/
      DATA (IFOREST(IFR, -4,   5),IFR=1,2)/  5, -3/
      DATA SPROP(  -4,   5)/      23/
C     Diagram 6, Amplitude 6
      DATA MAPCONFIG(   6)/   6/
      DATA (IFOREST(IFR, -1,   6),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,   6)/       5/
      DATA (IFOREST(IFR, -2,   6),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,   6)/      23/
      DATA (IFOREST(IFR, -3,   6),IFR=1,2)/  6, -2/
      DATA SPROP(  -3,   6)/     -13/
      DATA (IFOREST(IFR, -4,   6),IFR=1,2)/  5, -3/
      DATA SPROP(  -4,   6)/      23/
C     Diagram 7, Amplitude 7
      DATA MAPCONFIG(   7)/   7/
      DATA (IFOREST(IFR, -1,   7),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,   7)/       5/
      DATA (IFOREST(IFR, -2,   7),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,   7)/      22/
      DATA (IFOREST(IFR, -3,   7),IFR=1,2)/  5, -2/
      DATA SPROP(  -3,   7)/      13/
      DATA (IFOREST(IFR, -4,   7),IFR=1,2)/  6, -3/
      DATA SPROP(  -4,   7)/      23/
C     Diagram 8, Amplitude 8
      DATA MAPCONFIG(   8)/   8/
      DATA (IFOREST(IFR, -1,   8),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,   8)/       5/
      DATA (IFOREST(IFR, -2,   8),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,   8)/      23/
      DATA (IFOREST(IFR, -3,   8),IFR=1,2)/  5, -2/
      DATA SPROP(  -3,   8)/      13/
      DATA (IFOREST(IFR, -4,   8),IFR=1,2)/  6, -3/
      DATA SPROP(  -4,   8)/      23/
C     Diagram 9, Amplitude 9
      DATA MAPCONFIG(   9)/   9/
      DATA (IFOREST(IFR, -1,   9),IFR=1,2)/  6,  5/
      DATA SPROP(  -1,   9)/      22/
      DATA (IFOREST(IFR, -2,   9),IFR=1,2)/ -1,  4/
      DATA SPROP(  -2,   9)/      -5/
      DATA (IFOREST(IFR, -3,   9),IFR=1,2)/  7,  3/
      DATA SPROP(  -3,   9)/       5/
      DATA (IFOREST(IFR, -4,   9),IFR=1,2)/ -2, -3/
      DATA SPROP(  -4,   9)/      22/
C     Diagram 10, Amplitude 10
      DATA MAPCONFIG(  10)/  10/
      DATA (IFOREST(IFR, -1,  10),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,  10)/       5/
      DATA (IFOREST(IFR, -2,  10),IFR=1,2)/  6,  5/
      DATA SPROP(  -2,  10)/      22/
      DATA (IFOREST(IFR, -3,  10),IFR=1,2)/ -2, -1/
      DATA SPROP(  -3,  10)/       5/
      DATA (IFOREST(IFR, -4,  10),IFR=1,2)/  4, -3/
      DATA SPROP(  -4,  10)/      22/
C     Diagram 11, Amplitude 11
      DATA MAPCONFIG(  11)/  11/
      DATA (IFOREST(IFR, -1,  11),IFR=1,2)/  6,  5/
      DATA SPROP(  -1,  11)/      23/
      DATA (IFOREST(IFR, -2,  11),IFR=1,2)/ -1,  4/
      DATA SPROP(  -2,  11)/      -5/
      DATA (IFOREST(IFR, -3,  11),IFR=1,2)/  7,  3/
      DATA SPROP(  -3,  11)/       5/
      DATA (IFOREST(IFR, -4,  11),IFR=1,2)/ -2, -3/
      DATA SPROP(  -4,  11)/      22/
C     Diagram 12, Amplitude 12
      DATA MAPCONFIG(  12)/  12/
      DATA (IFOREST(IFR, -1,  12),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,  12)/       5/
      DATA (IFOREST(IFR, -2,  12),IFR=1,2)/  6,  5/
      DATA SPROP(  -2,  12)/      23/
      DATA (IFOREST(IFR, -3,  12),IFR=1,2)/ -2, -1/
      DATA SPROP(  -3,  12)/       5/
      DATA (IFOREST(IFR, -4,  12),IFR=1,2)/  4, -3/
      DATA SPROP(  -4,  12)/      22/
C     Diagram 13, Amplitude 13
      DATA MAPCONFIG(  13)/  13/
      DATA (IFOREST(IFR, -1,  13),IFR=1,2)/  6,  5/
      DATA SPROP(  -1,  13)/      22/
      DATA (IFOREST(IFR, -2,  13),IFR=1,2)/ -1,  4/
      DATA SPROP(  -2,  13)/      -5/
      DATA (IFOREST(IFR, -3,  13),IFR=1,2)/  7,  3/
      DATA SPROP(  -3,  13)/       5/
      DATA (IFOREST(IFR, -4,  13),IFR=1,2)/ -2, -3/
      DATA SPROP(  -4,  13)/      23/
C     Diagram 14, Amplitude 14
      DATA MAPCONFIG(  14)/  14/
      DATA (IFOREST(IFR, -1,  14),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,  14)/       5/
      DATA (IFOREST(IFR, -2,  14),IFR=1,2)/  6,  5/
      DATA SPROP(  -2,  14)/      22/
      DATA (IFOREST(IFR, -3,  14),IFR=1,2)/ -2, -1/
      DATA SPROP(  -3,  14)/       5/
      DATA (IFOREST(IFR, -4,  14),IFR=1,2)/  4, -3/
      DATA SPROP(  -4,  14)/      23/
C     Diagram 15, Amplitude 15
      DATA MAPCONFIG(  15)/  15/
      DATA (IFOREST(IFR, -1,  15),IFR=1,2)/  6,  5/
      DATA SPROP(  -1,  15)/      23/
      DATA (IFOREST(IFR, -2,  15),IFR=1,2)/ -1,  4/
      DATA SPROP(  -2,  15)/      -5/
      DATA (IFOREST(IFR, -3,  15),IFR=1,2)/  7,  3/
      DATA SPROP(  -3,  15)/       5/
      DATA (IFOREST(IFR, -4,  15),IFR=1,2)/ -2, -3/
      DATA SPROP(  -4,  15)/      23/
C     Diagram 16, Amplitude 16
      DATA MAPCONFIG(  16)/  16/
      DATA (IFOREST(IFR, -1,  16),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,  16)/       5/
      DATA (IFOREST(IFR, -2,  16),IFR=1,2)/  6,  5/
      DATA SPROP(  -2,  16)/      23/
      DATA (IFOREST(IFR, -3,  16),IFR=1,2)/ -2, -1/
      DATA SPROP(  -3,  16)/       5/
      DATA (IFOREST(IFR, -4,  16),IFR=1,2)/  4, -3/
      DATA SPROP(  -4,  16)/      23/
C     Diagram 17, Amplitude 17
      DATA MAPCONFIG(  17)/  17/
      DATA (IFOREST(IFR, -1,  17),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,  17)/       5/
      DATA (IFOREST(IFR, -2,  17),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,  17)/      25/
      DATA (IFOREST(IFR, -3,  17),IFR=1,2)/  6,  5/
      DATA SPROP(  -3,  17)/      23/
      DATA (IFOREST(IFR, -4,  17),IFR=1,2)/ -3, -2/
      DATA SPROP(  -4,  17)/      23/
C     Diagram 18, Amplitude 18
      DATA MAPCONFIG(  18)/  18/
      DATA (IFOREST(IFR, -1,  18),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  18)/      -5/
      DATA (IFOREST(IFR, -2,  18),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  18)/      22/
      DATA (IFOREST(IFR, -3,  18),IFR=1,2)/  6, -2/
      DATA SPROP(  -3,  18)/     -13/
      DATA (IFOREST(IFR, -4,  18),IFR=1,2)/  5, -3/
      DATA SPROP(  -4,  18)/      22/
C     Diagram 19, Amplitude 19
      DATA MAPCONFIG(  19)/  19/
      DATA (IFOREST(IFR, -1,  19),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  19)/      -5/
      DATA (IFOREST(IFR, -2,  19),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  19)/      23/
      DATA (IFOREST(IFR, -3,  19),IFR=1,2)/  6, -2/
      DATA SPROP(  -3,  19)/     -13/
      DATA (IFOREST(IFR, -4,  19),IFR=1,2)/  5, -3/
      DATA SPROP(  -4,  19)/      22/
C     Diagram 20, Amplitude 20
      DATA MAPCONFIG(  20)/  20/
      DATA (IFOREST(IFR, -1,  20),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  20)/      -5/
      DATA (IFOREST(IFR, -2,  20),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  20)/      22/
      DATA (IFOREST(IFR, -3,  20),IFR=1,2)/  5, -2/
      DATA SPROP(  -3,  20)/      13/
      DATA (IFOREST(IFR, -4,  20),IFR=1,2)/  6, -3/
      DATA SPROP(  -4,  20)/      22/
C     Diagram 21, Amplitude 21
      DATA MAPCONFIG(  21)/  21/
      DATA (IFOREST(IFR, -1,  21),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  21)/      -5/
      DATA (IFOREST(IFR, -2,  21),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  21)/      23/
      DATA (IFOREST(IFR, -3,  21),IFR=1,2)/  5, -2/
      DATA SPROP(  -3,  21)/      13/
      DATA (IFOREST(IFR, -4,  21),IFR=1,2)/  6, -3/
      DATA SPROP(  -4,  21)/      22/
C     Diagram 22, Amplitude 22
      DATA MAPCONFIG(  22)/  22/
      DATA (IFOREST(IFR, -1,  22),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  22)/      -5/
      DATA (IFOREST(IFR, -2,  22),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  22)/      22/
      DATA (IFOREST(IFR, -3,  22),IFR=1,2)/  6, -2/
      DATA SPROP(  -3,  22)/     -13/
      DATA (IFOREST(IFR, -4,  22),IFR=1,2)/  5, -3/
      DATA SPROP(  -4,  22)/      23/
C     Diagram 23, Amplitude 23
      DATA MAPCONFIG(  23)/  23/
      DATA (IFOREST(IFR, -1,  23),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  23)/      -5/
      DATA (IFOREST(IFR, -2,  23),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  23)/      23/
      DATA (IFOREST(IFR, -3,  23),IFR=1,2)/  6, -2/
      DATA SPROP(  -3,  23)/     -13/
      DATA (IFOREST(IFR, -4,  23),IFR=1,2)/  5, -3/
      DATA SPROP(  -4,  23)/      23/
C     Diagram 24, Amplitude 24
      DATA MAPCONFIG(  24)/  24/
      DATA (IFOREST(IFR, -1,  24),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  24)/      -5/
      DATA (IFOREST(IFR, -2,  24),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  24)/      22/
      DATA (IFOREST(IFR, -3,  24),IFR=1,2)/  5, -2/
      DATA SPROP(  -3,  24)/      13/
      DATA (IFOREST(IFR, -4,  24),IFR=1,2)/  6, -3/
      DATA SPROP(  -4,  24)/      23/
C     Diagram 25, Amplitude 25
      DATA MAPCONFIG(  25)/  25/
      DATA (IFOREST(IFR, -1,  25),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  25)/      -5/
      DATA (IFOREST(IFR, -2,  25),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  25)/      23/
      DATA (IFOREST(IFR, -3,  25),IFR=1,2)/  5, -2/
      DATA SPROP(  -3,  25)/      13/
      DATA (IFOREST(IFR, -4,  25),IFR=1,2)/  6, -3/
      DATA SPROP(  -4,  25)/      23/
C     Diagram 26, Amplitude 26
      DATA MAPCONFIG(  26)/  26/
      DATA (IFOREST(IFR, -1,  26),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  26)/      -5/
      DATA (IFOREST(IFR, -2,  26),IFR=1,2)/  6,  5/
      DATA SPROP(  -2,  26)/      22/
      DATA (IFOREST(IFR, -3,  26),IFR=1,2)/ -2, -1/
      DATA SPROP(  -3,  26)/      -5/
      DATA (IFOREST(IFR, -4,  26),IFR=1,2)/ -3,  3/
      DATA SPROP(  -4,  26)/      22/
C     Diagram 27, Amplitude 27
      DATA MAPCONFIG(  27)/  27/
      DATA (IFOREST(IFR, -1,  27),IFR=1,2)/  6,  5/
      DATA SPROP(  -1,  27)/      22/
      DATA (IFOREST(IFR, -2,  27),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  27)/       5/
      DATA (IFOREST(IFR, -3,  27),IFR=1,2)/  7,  4/
      DATA SPROP(  -3,  27)/      -5/
      DATA (IFOREST(IFR, -4,  27),IFR=1,2)/ -3, -2/
      DATA SPROP(  -4,  27)/      22/
C     Diagram 28, Amplitude 28
      DATA MAPCONFIG(  28)/  28/
      DATA (IFOREST(IFR, -1,  28),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  28)/      -5/
      DATA (IFOREST(IFR, -2,  28),IFR=1,2)/  6,  5/
      DATA SPROP(  -2,  28)/      23/
      DATA (IFOREST(IFR, -3,  28),IFR=1,2)/ -2, -1/
      DATA SPROP(  -3,  28)/      -5/
      DATA (IFOREST(IFR, -4,  28),IFR=1,2)/ -3,  3/
      DATA SPROP(  -4,  28)/      22/
C     Diagram 29, Amplitude 29
      DATA MAPCONFIG(  29)/  29/
      DATA (IFOREST(IFR, -1,  29),IFR=1,2)/  6,  5/
      DATA SPROP(  -1,  29)/      23/
      DATA (IFOREST(IFR, -2,  29),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  29)/       5/
      DATA (IFOREST(IFR, -3,  29),IFR=1,2)/  7,  4/
      DATA SPROP(  -3,  29)/      -5/
      DATA (IFOREST(IFR, -4,  29),IFR=1,2)/ -3, -2/
      DATA SPROP(  -4,  29)/      22/
C     Diagram 30, Amplitude 30
      DATA MAPCONFIG(  30)/  30/
      DATA (IFOREST(IFR, -1,  30),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  30)/      -5/
      DATA (IFOREST(IFR, -2,  30),IFR=1,2)/  6,  5/
      DATA SPROP(  -2,  30)/      22/
      DATA (IFOREST(IFR, -3,  30),IFR=1,2)/ -2, -1/
      DATA SPROP(  -3,  30)/      -5/
      DATA (IFOREST(IFR, -4,  30),IFR=1,2)/ -3,  3/
      DATA SPROP(  -4,  30)/      23/
C     Diagram 31, Amplitude 31
      DATA MAPCONFIG(  31)/  31/
      DATA (IFOREST(IFR, -1,  31),IFR=1,2)/  6,  5/
      DATA SPROP(  -1,  31)/      22/
      DATA (IFOREST(IFR, -2,  31),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  31)/       5/
      DATA (IFOREST(IFR, -3,  31),IFR=1,2)/  7,  4/
      DATA SPROP(  -3,  31)/      -5/
      DATA (IFOREST(IFR, -4,  31),IFR=1,2)/ -3, -2/
      DATA SPROP(  -4,  31)/      23/
C     Diagram 32, Amplitude 32
      DATA MAPCONFIG(  32)/  32/
      DATA (IFOREST(IFR, -1,  32),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  32)/      -5/
      DATA (IFOREST(IFR, -2,  32),IFR=1,2)/  6,  5/
      DATA SPROP(  -2,  32)/      23/
      DATA (IFOREST(IFR, -3,  32),IFR=1,2)/ -2, -1/
      DATA SPROP(  -3,  32)/      -5/
      DATA (IFOREST(IFR, -4,  32),IFR=1,2)/ -3,  3/
      DATA SPROP(  -4,  32)/      23/
C     Diagram 33, Amplitude 33
      DATA MAPCONFIG(  33)/  33/
      DATA (IFOREST(IFR, -1,  33),IFR=1,2)/  6,  5/
      DATA SPROP(  -1,  33)/      23/
      DATA (IFOREST(IFR, -2,  33),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  33)/       5/
      DATA (IFOREST(IFR, -3,  33),IFR=1,2)/  7,  4/
      DATA SPROP(  -3,  33)/      -5/
      DATA (IFOREST(IFR, -4,  33),IFR=1,2)/ -3, -2/
      DATA SPROP(  -4,  33)/      23/
C     Diagram 34, Amplitude 34
      DATA MAPCONFIG(  34)/  34/
      DATA (IFOREST(IFR, -1,  34),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  34)/      -5/
      DATA (IFOREST(IFR, -2,  34),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  34)/      25/
      DATA (IFOREST(IFR, -3,  34),IFR=1,2)/  6,  5/
      DATA SPROP(  -3,  34)/      23/
      DATA (IFOREST(IFR, -4,  34),IFR=1,2)/ -3, -2/
      DATA SPROP(  -4,  34)/      23/
C     Diagram 35, Amplitude 35
      DATA MAPCONFIG(  35)/  35/
      DATA (IFOREST(IFR, -1,  35),IFR=1,2)/  6,  5/
      DATA SPROP(  -1,  35)/      22/
      DATA (IFOREST(IFR, -2,  35),IFR=1,2)/ -1,  4/
      DATA SPROP(  -2,  35)/      -5/
      DATA (IFOREST(IFR, -3,  35),IFR=1,2)/  7, -2/
      DATA SPROP(  -3,  35)/      -5/
      DATA (IFOREST(IFR, -4,  35),IFR=1,2)/ -3,  3/
      DATA SPROP(  -4,  35)/      22/
C     Diagram 36, Amplitude 36
      DATA MAPCONFIG(  36)/  36/
      DATA (IFOREST(IFR, -1,  36),IFR=1,2)/  6,  5/
      DATA SPROP(  -1,  36)/      22/
      DATA (IFOREST(IFR, -2,  36),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  36)/       5/
      DATA (IFOREST(IFR, -3,  36),IFR=1,2)/  7, -2/
      DATA SPROP(  -3,  36)/       5/
      DATA (IFOREST(IFR, -4,  36),IFR=1,2)/  4, -3/
      DATA SPROP(  -4,  36)/      22/
C     Diagram 37, Amplitude 37
      DATA MAPCONFIG(  37)/  37/
      DATA (IFOREST(IFR, -1,  37),IFR=1,2)/  6,  5/
      DATA SPROP(  -1,  37)/      23/
      DATA (IFOREST(IFR, -2,  37),IFR=1,2)/ -1,  4/
      DATA SPROP(  -2,  37)/      -5/
      DATA (IFOREST(IFR, -3,  37),IFR=1,2)/  7, -2/
      DATA SPROP(  -3,  37)/      -5/
      DATA (IFOREST(IFR, -4,  37),IFR=1,2)/ -3,  3/
      DATA SPROP(  -4,  37)/      22/
C     Diagram 38, Amplitude 38
      DATA MAPCONFIG(  38)/  38/
      DATA (IFOREST(IFR, -1,  38),IFR=1,2)/  6,  5/
      DATA SPROP(  -1,  38)/      23/
      DATA (IFOREST(IFR, -2,  38),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  38)/       5/
      DATA (IFOREST(IFR, -3,  38),IFR=1,2)/  7, -2/
      DATA SPROP(  -3,  38)/       5/
      DATA (IFOREST(IFR, -4,  38),IFR=1,2)/  4, -3/
      DATA SPROP(  -4,  38)/      22/
C     Diagram 39, Amplitude 39
      DATA MAPCONFIG(  39)/  39/
      DATA (IFOREST(IFR, -1,  39),IFR=1,2)/  6,  5/
      DATA SPROP(  -1,  39)/      22/
      DATA (IFOREST(IFR, -2,  39),IFR=1,2)/ -1,  4/
      DATA SPROP(  -2,  39)/      -5/
      DATA (IFOREST(IFR, -3,  39),IFR=1,2)/  7, -2/
      DATA SPROP(  -3,  39)/      -5/
      DATA (IFOREST(IFR, -4,  39),IFR=1,2)/ -3,  3/
      DATA SPROP(  -4,  39)/      23/
C     Diagram 40, Amplitude 40
      DATA MAPCONFIG(  40)/  40/
      DATA (IFOREST(IFR, -1,  40),IFR=1,2)/  6,  5/
      DATA SPROP(  -1,  40)/      22/
      DATA (IFOREST(IFR, -2,  40),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  40)/       5/
      DATA (IFOREST(IFR, -3,  40),IFR=1,2)/  7, -2/
      DATA SPROP(  -3,  40)/       5/
      DATA (IFOREST(IFR, -4,  40),IFR=1,2)/  4, -3/
      DATA SPROP(  -4,  40)/      23/
C     Diagram 41, Amplitude 41
      DATA MAPCONFIG(  41)/  41/
      DATA (IFOREST(IFR, -1,  41),IFR=1,2)/  6,  5/
      DATA SPROP(  -1,  41)/      23/
      DATA (IFOREST(IFR, -2,  41),IFR=1,2)/ -1,  4/
      DATA SPROP(  -2,  41)/      -5/
      DATA (IFOREST(IFR, -3,  41),IFR=1,2)/  7, -2/
      DATA SPROP(  -3,  41)/      -5/
      DATA (IFOREST(IFR, -4,  41),IFR=1,2)/ -3,  3/
      DATA SPROP(  -4,  41)/      23/
C     Diagram 42, Amplitude 42
      DATA MAPCONFIG(  42)/  42/
      DATA (IFOREST(IFR, -1,  42),IFR=1,2)/  6,  5/
      DATA SPROP(  -1,  42)/      23/
      DATA (IFOREST(IFR, -2,  42),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  42)/       5/
      DATA (IFOREST(IFR, -3,  42),IFR=1,2)/  7, -2/
      DATA SPROP(  -3,  42)/       5/
      DATA (IFOREST(IFR, -4,  42),IFR=1,2)/  4, -3/
      DATA SPROP(  -4,  42)/      23/
C     Diagram 43, Amplitude 43
      DATA MAPCONFIG(  43)/  43/
      DATA (IFOREST(IFR, -1,  43),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,  43)/       5/
      DATA (IFOREST(IFR, -2,  43),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,  43)/      22/
      DATA (IFOREST(IFR, -3,  43),IFR=1,2)/  6,  5/
      DATA SPROP(  -3,  43)/      22/
      DATA (IFOREST(IFR, -4,  43),IFR=1,2)/  1, -3/
      DATA TPRID(  -4,  43)/      11/
      DATA (IFOREST(IFR, -5,  43),IFR=1,2)/ -4, -2/
C     Diagram 44, Amplitude 44
      DATA MAPCONFIG(  44)/  44/
      DATA (IFOREST(IFR, -1,  44),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,  44)/       5/
      DATA (IFOREST(IFR, -2,  44),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,  44)/      23/
      DATA (IFOREST(IFR, -3,  44),IFR=1,2)/  6,  5/
      DATA SPROP(  -3,  44)/      22/
      DATA (IFOREST(IFR, -4,  44),IFR=1,2)/  1, -3/
      DATA TPRID(  -4,  44)/      11/
      DATA (IFOREST(IFR, -5,  44),IFR=1,2)/ -4, -2/
C     Diagram 45, Amplitude 45
      DATA MAPCONFIG(  45)/  45/
      DATA (IFOREST(IFR, -1,  45),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,  45)/       5/
      DATA (IFOREST(IFR, -2,  45),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,  45)/      22/
      DATA (IFOREST(IFR, -3,  45),IFR=1,2)/  6,  5/
      DATA SPROP(  -3,  45)/      22/
      DATA (IFOREST(IFR, -4,  45),IFR=1,2)/  1, -2/
      DATA TPRID(  -4,  45)/      11/
      DATA (IFOREST(IFR, -5,  45),IFR=1,2)/ -4, -3/
C     Diagram 46, Amplitude 46
      DATA MAPCONFIG(  46)/  46/
      DATA (IFOREST(IFR, -1,  46),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,  46)/       5/
      DATA (IFOREST(IFR, -2,  46),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,  46)/      23/
      DATA (IFOREST(IFR, -3,  46),IFR=1,2)/  6,  5/
      DATA SPROP(  -3,  46)/      22/
      DATA (IFOREST(IFR, -4,  46),IFR=1,2)/  1, -2/
      DATA TPRID(  -4,  46)/      11/
      DATA (IFOREST(IFR, -5,  46),IFR=1,2)/ -4, -3/
C     Diagram 47, Amplitude 47
      DATA MAPCONFIG(  47)/  47/
      DATA (IFOREST(IFR, -1,  47),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,  47)/       5/
      DATA (IFOREST(IFR, -2,  47),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,  47)/      22/
      DATA (IFOREST(IFR, -3,  47),IFR=1,2)/  6,  5/
      DATA SPROP(  -3,  47)/      23/
      DATA (IFOREST(IFR, -4,  47),IFR=1,2)/  1, -3/
      DATA TPRID(  -4,  47)/      11/
      DATA (IFOREST(IFR, -5,  47),IFR=1,2)/ -4, -2/
C     Diagram 48, Amplitude 48
      DATA MAPCONFIG(  48)/  48/
      DATA (IFOREST(IFR, -1,  48),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,  48)/       5/
      DATA (IFOREST(IFR, -2,  48),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,  48)/      23/
      DATA (IFOREST(IFR, -3,  48),IFR=1,2)/  6,  5/
      DATA SPROP(  -3,  48)/      23/
      DATA (IFOREST(IFR, -4,  48),IFR=1,2)/  1, -3/
      DATA TPRID(  -4,  48)/      11/
      DATA (IFOREST(IFR, -5,  48),IFR=1,2)/ -4, -2/
C     Diagram 49, Amplitude 49
      DATA MAPCONFIG(  49)/  49/
      DATA (IFOREST(IFR, -1,  49),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,  49)/       5/
      DATA (IFOREST(IFR, -2,  49),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,  49)/      22/
      DATA (IFOREST(IFR, -3,  49),IFR=1,2)/  6,  5/
      DATA SPROP(  -3,  49)/      23/
      DATA (IFOREST(IFR, -4,  49),IFR=1,2)/  1, -2/
      DATA TPRID(  -4,  49)/      11/
      DATA (IFOREST(IFR, -5,  49),IFR=1,2)/ -4, -3/
C     Diagram 50, Amplitude 50
      DATA MAPCONFIG(  50)/  50/
      DATA (IFOREST(IFR, -1,  50),IFR=1,2)/  7,  3/
      DATA SPROP(  -1,  50)/       5/
      DATA (IFOREST(IFR, -2,  50),IFR=1,2)/  4, -1/
      DATA SPROP(  -2,  50)/      23/
      DATA (IFOREST(IFR, -3,  50),IFR=1,2)/  6,  5/
      DATA SPROP(  -3,  50)/      23/
      DATA (IFOREST(IFR, -4,  50),IFR=1,2)/  1, -2/
      DATA TPRID(  -4,  50)/      11/
      DATA (IFOREST(IFR, -5,  50),IFR=1,2)/ -4, -3/
C     Diagram 51, Amplitude 51
      DATA MAPCONFIG(  51)/  51/
      DATA (IFOREST(IFR, -1,  51),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  51)/      -5/
      DATA (IFOREST(IFR, -2,  51),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  51)/      22/
      DATA (IFOREST(IFR, -3,  51),IFR=1,2)/  6,  5/
      DATA SPROP(  -3,  51)/      22/
      DATA (IFOREST(IFR, -4,  51),IFR=1,2)/  1, -3/
      DATA TPRID(  -4,  51)/      11/
      DATA (IFOREST(IFR, -5,  51),IFR=1,2)/ -4, -2/
C     Diagram 52, Amplitude 52
      DATA MAPCONFIG(  52)/  52/
      DATA (IFOREST(IFR, -1,  52),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  52)/      -5/
      DATA (IFOREST(IFR, -2,  52),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  52)/      23/
      DATA (IFOREST(IFR, -3,  52),IFR=1,2)/  6,  5/
      DATA SPROP(  -3,  52)/      22/
      DATA (IFOREST(IFR, -4,  52),IFR=1,2)/  1, -3/
      DATA TPRID(  -4,  52)/      11/
      DATA (IFOREST(IFR, -5,  52),IFR=1,2)/ -4, -2/
C     Diagram 53, Amplitude 53
      DATA MAPCONFIG(  53)/  53/
      DATA (IFOREST(IFR, -1,  53),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  53)/      -5/
      DATA (IFOREST(IFR, -2,  53),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  53)/      22/
      DATA (IFOREST(IFR, -3,  53),IFR=1,2)/  6,  5/
      DATA SPROP(  -3,  53)/      22/
      DATA (IFOREST(IFR, -4,  53),IFR=1,2)/  1, -2/
      DATA TPRID(  -4,  53)/      11/
      DATA (IFOREST(IFR, -5,  53),IFR=1,2)/ -4, -3/
C     Diagram 54, Amplitude 54
      DATA MAPCONFIG(  54)/  54/
      DATA (IFOREST(IFR, -1,  54),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  54)/      -5/
      DATA (IFOREST(IFR, -2,  54),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  54)/      23/
      DATA (IFOREST(IFR, -3,  54),IFR=1,2)/  6,  5/
      DATA SPROP(  -3,  54)/      22/
      DATA (IFOREST(IFR, -4,  54),IFR=1,2)/  1, -2/
      DATA TPRID(  -4,  54)/      11/
      DATA (IFOREST(IFR, -5,  54),IFR=1,2)/ -4, -3/
C     Diagram 55, Amplitude 55
      DATA MAPCONFIG(  55)/  55/
      DATA (IFOREST(IFR, -1,  55),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  55)/      -5/
      DATA (IFOREST(IFR, -2,  55),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  55)/      22/
      DATA (IFOREST(IFR, -3,  55),IFR=1,2)/  6,  5/
      DATA SPROP(  -3,  55)/      23/
      DATA (IFOREST(IFR, -4,  55),IFR=1,2)/  1, -3/
      DATA TPRID(  -4,  55)/      11/
      DATA (IFOREST(IFR, -5,  55),IFR=1,2)/ -4, -2/
C     Diagram 56, Amplitude 56
      DATA MAPCONFIG(  56)/  56/
      DATA (IFOREST(IFR, -1,  56),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  56)/      -5/
      DATA (IFOREST(IFR, -2,  56),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  56)/      23/
      DATA (IFOREST(IFR, -3,  56),IFR=1,2)/  6,  5/
      DATA SPROP(  -3,  56)/      23/
      DATA (IFOREST(IFR, -4,  56),IFR=1,2)/  1, -3/
      DATA TPRID(  -4,  56)/      11/
      DATA (IFOREST(IFR, -5,  56),IFR=1,2)/ -4, -2/
C     Diagram 57, Amplitude 57
      DATA MAPCONFIG(  57)/  57/
      DATA (IFOREST(IFR, -1,  57),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  57)/      -5/
      DATA (IFOREST(IFR, -2,  57),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  57)/      22/
      DATA (IFOREST(IFR, -3,  57),IFR=1,2)/  6,  5/
      DATA SPROP(  -3,  57)/      23/
      DATA (IFOREST(IFR, -4,  57),IFR=1,2)/  1, -2/
      DATA TPRID(  -4,  57)/      11/
      DATA (IFOREST(IFR, -5,  57),IFR=1,2)/ -4, -3/
C     Diagram 58, Amplitude 58
      DATA MAPCONFIG(  58)/  58/
      DATA (IFOREST(IFR, -1,  58),IFR=1,2)/  7,  4/
      DATA SPROP(  -1,  58)/      -5/
      DATA (IFOREST(IFR, -2,  58),IFR=1,2)/ -1,  3/
      DATA SPROP(  -2,  58)/      23/
      DATA (IFOREST(IFR, -3,  58),IFR=1,2)/  6,  5/
      DATA SPROP(  -3,  58)/      23/
      DATA (IFOREST(IFR, -4,  58),IFR=1,2)/  1, -2/
      DATA TPRID(  -4,  58)/      11/
      DATA (IFOREST(IFR, -5,  58),IFR=1,2)/ -4, -3/
C     Number of configs
      DATA MAPCONFIG(0)/  58/

C     Here are the BWs
      DATA GFORCEBW(-1,1)/.FALSE./
      DATA GFORCEBW(-2,1)/.FALSE./
      DATA GFORCEBW(-3,1)/.FALSE./
      DATA GFORCEBW(-4,1)/.FALSE./
      DATA GFORCEBW(-1,2)/.FALSE./
      DATA GFORCEBW(-2,2)/.FALSE./
      DATA GFORCEBW(-3,2)/.FALSE./
      DATA GFORCEBW(-4,2)/.FALSE./
      DATA GFORCEBW(-1,3)/.FALSE./
      DATA GFORCEBW(-2,3)/.FALSE./
      DATA GFORCEBW(-3,3)/.FALSE./
      DATA GFORCEBW(-4,3)/.FALSE./
      DATA GFORCEBW(-1,4)/.FALSE./
      DATA GFORCEBW(-2,4)/.FALSE./
      DATA GFORCEBW(-3,4)/.FALSE./
      DATA GFORCEBW(-4,4)/.FALSE./
      DATA GFORCEBW(-1,5)/.FALSE./
      DATA GFORCEBW(-2,5)/.FALSE./
      DATA GFORCEBW(-3,5)/.FALSE./
      DATA GFORCEBW(-4,5)/.FALSE./
      DATA GFORCEBW(-1,6)/.FALSE./
      DATA GFORCEBW(-2,6)/.FALSE./
      DATA GFORCEBW(-3,6)/.FALSE./
      DATA GFORCEBW(-4,6)/.FALSE./
      DATA GFORCEBW(-1,7)/.FALSE./
      DATA GFORCEBW(-2,7)/.FALSE./
      DATA GFORCEBW(-3,7)/.FALSE./
      DATA GFORCEBW(-4,7)/.FALSE./
      DATA GFORCEBW(-1,8)/.FALSE./
      DATA GFORCEBW(-2,8)/.FALSE./
      DATA GFORCEBW(-3,8)/.FALSE./
      DATA GFORCEBW(-4,8)/.FALSE./
      DATA GFORCEBW(-1,9)/.FALSE./
      DATA GFORCEBW(-2,9)/.FALSE./
      DATA GFORCEBW(-3,9)/.FALSE./
      DATA GFORCEBW(-4,9)/.FALSE./
      DATA GFORCEBW(-1,10)/.FALSE./
      DATA GFORCEBW(-2,10)/.FALSE./
      DATA GFORCEBW(-3,10)/.FALSE./
      DATA GFORCEBW(-4,10)/.FALSE./
      DATA GFORCEBW(-1,11)/.FALSE./
      DATA GFORCEBW(-2,11)/.FALSE./
      DATA GFORCEBW(-3,11)/.FALSE./
      DATA GFORCEBW(-4,11)/.FALSE./
      DATA GFORCEBW(-1,12)/.FALSE./
      DATA GFORCEBW(-2,12)/.FALSE./
      DATA GFORCEBW(-3,12)/.FALSE./
      DATA GFORCEBW(-4,12)/.FALSE./
      DATA GFORCEBW(-1,13)/.FALSE./
      DATA GFORCEBW(-2,13)/.FALSE./
      DATA GFORCEBW(-3,13)/.FALSE./
      DATA GFORCEBW(-4,13)/.FALSE./
      DATA GFORCEBW(-1,14)/.FALSE./
      DATA GFORCEBW(-2,14)/.FALSE./
      DATA GFORCEBW(-3,14)/.FALSE./
      DATA GFORCEBW(-4,14)/.FALSE./
      DATA GFORCEBW(-1,15)/.FALSE./
      DATA GFORCEBW(-2,15)/.FALSE./
      DATA GFORCEBW(-3,15)/.FALSE./
      DATA GFORCEBW(-4,15)/.FALSE./
      DATA GFORCEBW(-1,16)/.FALSE./
      DATA GFORCEBW(-2,16)/.FALSE./
      DATA GFORCEBW(-3,16)/.FALSE./
      DATA GFORCEBW(-4,16)/.FALSE./
      DATA GFORCEBW(-1,17)/.FALSE./
      DATA GFORCEBW(-2,17)/.FALSE./
      DATA GFORCEBW(-3,17)/.FALSE./
      DATA GFORCEBW(-4,17)/.FALSE./
      DATA GFORCEBW(-1,18)/.FALSE./
      DATA GFORCEBW(-2,18)/.FALSE./
      DATA GFORCEBW(-3,18)/.FALSE./
      DATA GFORCEBW(-4,18)/.FALSE./
      DATA GFORCEBW(-1,19)/.FALSE./
      DATA GFORCEBW(-2,19)/.FALSE./
      DATA GFORCEBW(-3,19)/.FALSE./
      DATA GFORCEBW(-4,19)/.FALSE./
      DATA GFORCEBW(-1,20)/.FALSE./
      DATA GFORCEBW(-2,20)/.FALSE./
      DATA GFORCEBW(-3,20)/.FALSE./
      DATA GFORCEBW(-4,20)/.FALSE./
      DATA GFORCEBW(-1,21)/.FALSE./
      DATA GFORCEBW(-2,21)/.FALSE./
      DATA GFORCEBW(-3,21)/.FALSE./
      DATA GFORCEBW(-4,21)/.FALSE./
      DATA GFORCEBW(-1,22)/.FALSE./
      DATA GFORCEBW(-2,22)/.FALSE./
      DATA GFORCEBW(-3,22)/.FALSE./
      DATA GFORCEBW(-4,22)/.FALSE./
      DATA GFORCEBW(-1,23)/.FALSE./
      DATA GFORCEBW(-2,23)/.FALSE./
      DATA GFORCEBW(-3,23)/.FALSE./
      DATA GFORCEBW(-4,23)/.FALSE./
      DATA GFORCEBW(-1,24)/.FALSE./
      DATA GFORCEBW(-2,24)/.FALSE./
      DATA GFORCEBW(-3,24)/.FALSE./
      DATA GFORCEBW(-4,24)/.FALSE./
      DATA GFORCEBW(-1,25)/.FALSE./
      DATA GFORCEBW(-2,25)/.FALSE./
      DATA GFORCEBW(-3,25)/.FALSE./
      DATA GFORCEBW(-4,25)/.FALSE./
      DATA GFORCEBW(-1,26)/.FALSE./
      DATA GFORCEBW(-2,26)/.FALSE./
      DATA GFORCEBW(-3,26)/.FALSE./
      DATA GFORCEBW(-4,26)/.FALSE./
      DATA GFORCEBW(-1,27)/.FALSE./
      DATA GFORCEBW(-2,27)/.FALSE./
      DATA GFORCEBW(-3,27)/.FALSE./
      DATA GFORCEBW(-4,27)/.FALSE./
      DATA GFORCEBW(-1,28)/.FALSE./
      DATA GFORCEBW(-2,28)/.FALSE./
      DATA GFORCEBW(-3,28)/.FALSE./
      DATA GFORCEBW(-4,28)/.FALSE./
      DATA GFORCEBW(-1,29)/.FALSE./
      DATA GFORCEBW(-2,29)/.FALSE./
      DATA GFORCEBW(-3,29)/.FALSE./
      DATA GFORCEBW(-4,29)/.FALSE./
      DATA GFORCEBW(-1,30)/.FALSE./
      DATA GFORCEBW(-2,30)/.FALSE./
      DATA GFORCEBW(-3,30)/.FALSE./
      DATA GFORCEBW(-4,30)/.FALSE./
      DATA GFORCEBW(-1,31)/.FALSE./
      DATA GFORCEBW(-2,31)/.FALSE./
      DATA GFORCEBW(-3,31)/.FALSE./
      DATA GFORCEBW(-4,31)/.FALSE./
      DATA GFORCEBW(-1,32)/.FALSE./
      DATA GFORCEBW(-2,32)/.FALSE./
      DATA GFORCEBW(-3,32)/.FALSE./
      DATA GFORCEBW(-4,32)/.FALSE./
      DATA GFORCEBW(-1,33)/.FALSE./
      DATA GFORCEBW(-2,33)/.FALSE./
      DATA GFORCEBW(-3,33)/.FALSE./
      DATA GFORCEBW(-4,33)/.FALSE./
      DATA GFORCEBW(-1,34)/.FALSE./
      DATA GFORCEBW(-2,34)/.FALSE./
      DATA GFORCEBW(-3,34)/.FALSE./
      DATA GFORCEBW(-4,34)/.FALSE./
      DATA GFORCEBW(-1,35)/.FALSE./
      DATA GFORCEBW(-2,35)/.FALSE./
      DATA GFORCEBW(-3,35)/.FALSE./
      DATA GFORCEBW(-4,35)/.FALSE./
      DATA GFORCEBW(-1,36)/.FALSE./
      DATA GFORCEBW(-2,36)/.FALSE./
      DATA GFORCEBW(-3,36)/.FALSE./
      DATA GFORCEBW(-4,36)/.FALSE./
      DATA GFORCEBW(-1,37)/.FALSE./
      DATA GFORCEBW(-2,37)/.FALSE./
      DATA GFORCEBW(-3,37)/.FALSE./
      DATA GFORCEBW(-4,37)/.FALSE./
      DATA GFORCEBW(-1,38)/.FALSE./
      DATA GFORCEBW(-2,38)/.FALSE./
      DATA GFORCEBW(-3,38)/.FALSE./
      DATA GFORCEBW(-4,38)/.FALSE./
      DATA GFORCEBW(-1,39)/.FALSE./
      DATA GFORCEBW(-2,39)/.FALSE./
      DATA GFORCEBW(-3,39)/.FALSE./
      DATA GFORCEBW(-4,39)/.FALSE./
      DATA GFORCEBW(-1,40)/.FALSE./
      DATA GFORCEBW(-2,40)/.FALSE./
      DATA GFORCEBW(-3,40)/.FALSE./
      DATA GFORCEBW(-4,40)/.FALSE./
      DATA GFORCEBW(-1,41)/.FALSE./
      DATA GFORCEBW(-2,41)/.FALSE./
      DATA GFORCEBW(-3,41)/.FALSE./
      DATA GFORCEBW(-4,41)/.FALSE./
      DATA GFORCEBW(-1,42)/.FALSE./
      DATA GFORCEBW(-2,42)/.FALSE./
      DATA GFORCEBW(-3,42)/.FALSE./
      DATA GFORCEBW(-4,42)/.FALSE./
      DATA GFORCEBW(-1,43)/.FALSE./
      DATA GFORCEBW(-2,43)/.FALSE./
      DATA GFORCEBW(-3,43)/.FALSE./
      DATA GFORCEBW(-1,44)/.FALSE./
      DATA GFORCEBW(-2,44)/.FALSE./
      DATA GFORCEBW(-3,44)/.FALSE./
      DATA GFORCEBW(-1,45)/.FALSE./
      DATA GFORCEBW(-2,45)/.FALSE./
      DATA GFORCEBW(-3,45)/.FALSE./
      DATA GFORCEBW(-1,46)/.FALSE./
      DATA GFORCEBW(-2,46)/.FALSE./
      DATA GFORCEBW(-3,46)/.FALSE./
      DATA GFORCEBW(-1,47)/.FALSE./
      DATA GFORCEBW(-2,47)/.FALSE./
      DATA GFORCEBW(-3,47)/.FALSE./
      DATA GFORCEBW(-1,48)/.FALSE./
      DATA GFORCEBW(-2,48)/.FALSE./
      DATA GFORCEBW(-3,48)/.FALSE./
      DATA GFORCEBW(-1,49)/.FALSE./
      DATA GFORCEBW(-2,49)/.FALSE./
      DATA GFORCEBW(-3,49)/.FALSE./
      DATA GFORCEBW(-1,50)/.FALSE./
      DATA GFORCEBW(-2,50)/.FALSE./
      DATA GFORCEBW(-3,50)/.FALSE./
      DATA GFORCEBW(-1,51)/.FALSE./
      DATA GFORCEBW(-2,51)/.FALSE./
      DATA GFORCEBW(-3,51)/.FALSE./
      DATA GFORCEBW(-1,52)/.FALSE./
      DATA GFORCEBW(-2,52)/.FALSE./
      DATA GFORCEBW(-3,52)/.FALSE./
      DATA GFORCEBW(-1,53)/.FALSE./
      DATA GFORCEBW(-2,53)/.FALSE./
      DATA GFORCEBW(-3,53)/.FALSE./
      DATA GFORCEBW(-1,54)/.FALSE./
      DATA GFORCEBW(-2,54)/.FALSE./
      DATA GFORCEBW(-3,54)/.FALSE./
      DATA GFORCEBW(-1,55)/.FALSE./
      DATA GFORCEBW(-2,55)/.FALSE./
      DATA GFORCEBW(-3,55)/.FALSE./
      DATA GFORCEBW(-1,56)/.FALSE./
      DATA GFORCEBW(-2,56)/.FALSE./
      DATA GFORCEBW(-3,56)/.FALSE./
      DATA GFORCEBW(-1,57)/.FALSE./
      DATA GFORCEBW(-2,57)/.FALSE./
      DATA GFORCEBW(-3,57)/.FALSE./
      DATA GFORCEBW(-1,58)/.FALSE./
      DATA GFORCEBW(-2,58)/.FALSE./
      DATA GFORCEBW(-3,58)/.FALSE./
