      INTEGER IPOS, JPOS
      INTEGER FKS_I_D(5), FKS_J_D(5)
      INTEGER EXTRA_CNT_D(5), ISPLITORDER_BORN_D(5),
     $  ISPLITORDER_CNT_D(5)
      INTEGER FKS_J_FROM_I_D(5, NEXTERNAL, 0:NEXTERNAL)
      INTEGER PARTICLE_TYPE_D(5, NEXTERNAL), PDG_TYPE_D(5, NEXTERNAL)
      LOGICAL PARTICLE_TAG_D(5, NEXTERNAL)
      REAL*8 PARTICLE_CHARGE_D(5, NEXTERNAL)
      LOGICAL SPLIT_TYPE_D(5, 2)
      LOGICAL NEED_COLOR_LINKS_D(5), NEED_CHARGE_LINKS_D(5)

      DATA FKS_I_D / 8, 8, 8, 8, 8 /
      DATA FKS_J_D / 3, 4, 7, 7, 7 /

      DATA NEED_COLOR_LINKS_D / .TRUE., .TRUE., .TRUE., .FALSE.,
     $  .FALSE. /
      DATA NEED_CHARGE_LINKS_D / .FALSE., .FALSE., .FALSE., .FALSE.,
     $  .FALSE. /

      DATA (FKS_J_FROM_I_D(1, 7, JPOS), JPOS = 0, 3)  / 3, 3, 4, 8 /
      DATA (FKS_J_FROM_I_D(1, 8, JPOS), JPOS = 0, 3)  / 3, 3, 4, 7 /

      DATA (FKS_J_FROM_I_D(2, 7, JPOS), JPOS = 0, 3)  / 3, 3, 4, 8 /
      DATA (FKS_J_FROM_I_D(2, 8, JPOS), JPOS = 0, 3)  / 3, 3, 4, 7 /

      DATA (FKS_J_FROM_I_D(3, 7, JPOS), JPOS = 0, 3)  / 3, 3, 4, 8 /
      DATA (FKS_J_FROM_I_D(3, 8, JPOS), JPOS = 0, 3)  / 3, 3, 4, 7 /

      DATA (FKS_J_FROM_I_D(4, 8, JPOS), JPOS = 0, 1)  / 1, 7 /

      DATA (FKS_J_FROM_I_D(5, 8, JPOS), JPOS = 0, 1)  / 1, 7 /


C     Split type: the order of elements is QCD, QED
      DATA (SPLIT_TYPE_D (1, IPOS), IPOS=1, 2) / .TRUE., .FALSE. /
      DATA (SPLIT_TYPE_D (2, IPOS), IPOS=1, 2) / .TRUE., .FALSE. /
      DATA (SPLIT_TYPE_D (3, IPOS), IPOS=1, 2) / .TRUE., .FALSE. /
      DATA (SPLIT_TYPE_D (4, IPOS), IPOS=1, 2) / .TRUE., .FALSE. /
      DATA (SPLIT_TYPE_D (5, IPOS), IPOS=1, 2) / .TRUE., .FALSE. /

C     Extra cnt: if the real emission needs an extra cnt (g/a >
C      qqbar), then
C     this is the id of the extra cnt. Otherwise it should be set to 0.
C     If not 0, store also the position in the splitorders array of
C      the 
C     coupling corresponding to the splitting type
      DATA EXTRA_CNT_D / 0, 0, 0, 0, 0 /
      DATA ISPLITORDER_BORN_D / 0, 0, 0, 0, 0 /
      DATA ISPLITORDER_CNT_D / 0, 0, 0, 0, 0 /

C     
C     Particle type:
C     octet = 8, triplet = 3, singlet = 1
      DATA (PARTICLE_TYPE_D(1, IPOS), IPOS=1, NEXTERNAL) / 1, 1, 3, -3
     $ , 1, 1, 8, 8 /
      DATA (PARTICLE_TYPE_D(2, IPOS), IPOS=1, NEXTERNAL) / 1, 1, 3, -3
     $ , 1, 1, 8, 8 /
      DATA (PARTICLE_TYPE_D(3, IPOS), IPOS=1, NEXTERNAL) / 1, 1, 3, -3
     $ , 1, 1, 8, 8 /
      DATA (PARTICLE_TYPE_D(4, IPOS), IPOS=1, NEXTERNAL) / 1, 1, 3, -3
     $ , 1, 1, 3, -3 /
      DATA (PARTICLE_TYPE_D(5, IPOS), IPOS=1, NEXTERNAL) / 1, 1, 3, -3
     $ , 1, 1, 3, -3 /

C     
C     Particle type according to PDG:
C     
      DATA (PDG_TYPE_D(1, IPOS), IPOS=1, NEXTERNAL) / -11, 11, 5, -5,
     $  13, -13, 21, 21 /
      DATA (PDG_TYPE_D(2, IPOS), IPOS=1, NEXTERNAL) / -11, 11, 5, -5,
     $  13, -13, 21, 21 /
      DATA (PDG_TYPE_D(3, IPOS), IPOS=1, NEXTERNAL) / -11, 11, 5, -5,
     $  13, -13, 21, 21 /
      DATA (PDG_TYPE_D(4, IPOS), IPOS=1, NEXTERNAL) / -11, 11, 5, -5,
     $  13, -13, 1, -1 /
      DATA (PDG_TYPE_D(5, IPOS), IPOS=1, NEXTERNAL) / -11, 11, 5, -5,
     $  13, -13, 2, -2 /

C     
C     Particle charge:
C     charge is set 0. with QCD corrections, which is irrelevant
      DATA (PARTICLE_CHARGE_D(1, IPOS), IPOS=1, NEXTERNAL) /  
     $  1.000000000000000D0,  -1.000000000000000D0,  
     $ -0.333333333333333D0,   0.333333333333333D0,  
     $ -1.000000000000000D0,   1.000000000000000D0,  
     $  0.000000000000000D0,   0.000000000000000D0 /
      DATA (PARTICLE_CHARGE_D(2, IPOS), IPOS=1, NEXTERNAL) /  
     $  1.000000000000000D0,  -1.000000000000000D0,  
     $ -0.333333333333333D0,   0.333333333333333D0,  
     $ -1.000000000000000D0,   1.000000000000000D0,  
     $  0.000000000000000D0,   0.000000000000000D0 /
      DATA (PARTICLE_CHARGE_D(3, IPOS), IPOS=1, NEXTERNAL) /  
     $  1.000000000000000D0,  -1.000000000000000D0,  
     $ -0.333333333333333D0,   0.333333333333333D0,  
     $ -1.000000000000000D0,   1.000000000000000D0,  
     $  0.000000000000000D0,   0.000000000000000D0 /
      DATA (PARTICLE_CHARGE_D(4, IPOS), IPOS=1, NEXTERNAL) /  
     $  1.000000000000000D0,  -1.000000000000000D0,  
     $ -0.333333333333333D0,   0.333333333333333D0,  
     $ -1.000000000000000D0,   1.000000000000000D0,  
     $ -0.333333333333333D0,   0.333333333333333D0 /
      DATA (PARTICLE_CHARGE_D(5, IPOS), IPOS=1, NEXTERNAL) /  
     $  1.000000000000000D0,  -1.000000000000000D0,  
     $ -0.333333333333333D0,   0.333333333333333D0,  
     $ -1.000000000000000D0,   1.000000000000000D0,  
     $  0.666666666666667D0,  -0.666666666666667D0 /

C     
C     Tagged particles:
C     
      DATA (PARTICLE_TAG_D(1, IPOS), IPOS=1, NEXTERNAL) / .FALSE.,
     $  .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /
      DATA (PARTICLE_TAG_D(2, IPOS), IPOS=1, NEXTERNAL) / .FALSE.,
     $  .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /
      DATA (PARTICLE_TAG_D(3, IPOS), IPOS=1, NEXTERNAL) / .FALSE.,
     $  .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /
      DATA (PARTICLE_TAG_D(4, IPOS), IPOS=1, NEXTERNAL) / .FALSE.,
     $  .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /
      DATA (PARTICLE_TAG_D(5, IPOS), IPOS=1, NEXTERNAL) / .FALSE.,
     $  .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /


