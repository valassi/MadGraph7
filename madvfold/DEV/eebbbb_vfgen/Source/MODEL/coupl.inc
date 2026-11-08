ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     
C     NB: VECSIZE_MEMMAX is defined in vector.inc
C     NB: vector.inc must be included before coupl.inc
C     

      DOUBLE PRECISION G, ALL_G
      COMMON/STRONG/ G, ALL_G

      DOUBLE COMPLEX GAL(2)
      COMMON/WEAK/ GAL

      DOUBLE PRECISION MU_R, ALL_MU_R
      COMMON/RSCALE/ MU_R, ALL_MU_R

      DOUBLE PRECISION NF
      PARAMETER(NF=4D0)
      DOUBLE PRECISION NL
      PARAMETER(NL=2D0)

      DOUBLE PRECISION MDL_MB,MDL_MH,MDL_MT,MDL_MTA,MDL_MW,MDL_MZ

      COMMON/MASSES/ MDL_MB,MDL_MH,MDL_MT,MDL_MTA,MDL_MW,MDL_MZ


      DOUBLE PRECISION MDL_WH,MDL_WT,MDL_WW,MDL_WZ

      COMMON/WIDTHS/ MDL_WH,MDL_WT,MDL_WW,MDL_WZ


      DOUBLE COMPLEX, TARGET :: GC_1, GC_3, GC_21, GC_22, GC_23, GC_24

      DOUBLE COMPLEX, TARGET :: GC_4, GC_5, R2_GQQ, R2_GGQ, R2_GGB,
     $  R2_GGT, R2_GGG_1, R2_GGG_2, R2_QQQ, R2_QQB, UV_GQQG_1EPS,
     $  UV_GQQB_1EPS, UV_BMASS_1EPS, UVWFCT_B_0_1EPS, R2_DDA,
     $  R2_DDZ_V2, R2_DDZ_V3, R2_GGZUP, R2_GGZDOWN, UV_GQQB, UV_GQQT,
     $  UV_BMASS, UVWFCT_B_0

      COMMON/COUPLINGS/ GC_4, GC_5, R2_GQQ, R2_GGQ, R2_GGB, R2_GGT,
     $  R2_GGG_1, R2_GGG_2, R2_QQQ, R2_QQB, UV_GQQG_1EPS, UV_GQQB_1EPS
     $ , UV_BMASS_1EPS, UVWFCT_B_0_1EPS, R2_DDA, R2_DDZ_V2, R2_DDZ_V3,
     $  R2_GGZUP, R2_GGZDOWN, UV_GQQB, UV_GQQT, UV_BMASS, UVWFCT_B_0,
     $  GC_1, GC_3, GC_21, GC_22, GC_23, GC_24

