ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      REAL*16 MP__G
      COMMON/MP_STRONG/ MP__G

      COMPLEX*32 MP__GAL(2)
      COMMON/MP_WEAK/ MP__GAL

      REAL*16 MP__MU_R
      COMMON/MP_RSCALE/ MP__MU_R


      REAL*16 MP__MDL_MB,MP__MDL_MH,MP__MDL_MT,MP__MDL_MTA,MP__MDL_MW
     $ ,MP__MDL_MZ

      COMMON/MP_MASSES/ MP__MDL_MB,MP__MDL_MH,MP__MDL_MT,MP__MDL_MTA
     $ ,MP__MDL_MW,MP__MDL_MZ


      REAL*16 MP__MDL_WH,MP__MDL_WT,MP__MDL_WW,MP__MDL_WZ

      COMMON/MP_WIDTHS/ MP__MDL_WH,MP__MDL_WT,MP__MDL_WW,MP__MDL_WZ


      COMPLEX*32 MP__GC_1,MP__GC_3,MP__GC_21,MP__GC_22,MP__GC_23
     $ ,MP__GC_24

      COMPLEX*32 MP__GC_4,MP__GC_5,MP__R2_GQQ,MP__R2_GGQ,MP__R2_GGB
     $ ,MP__R2_GGT,MP__R2_GGG_1,MP__R2_GGG_2,MP__R2_QQQ,MP__R2_QQB
     $ ,MP__UV_GQQG_1EPS,MP__UV_GQQB_1EPS,MP__UV_BMASS_1EPS
     $ ,MP__UVWFCT_B_0_1EPS,MP__R2_DDA,MP__R2_DDZ_V2,MP__R2_DDZ_V3
     $ ,MP__R2_GGZUP,MP__R2_GGZDOWN,MP__UV_GQQB,MP__UV_GQQT
     $ ,MP__UV_BMASS,MP__UVWFCT_B_0

      COMMON/MP_COUPLINGS/ MP__GC_4,MP__GC_5,MP__R2_GQQ,MP__R2_GGQ
     $ ,MP__R2_GGB,MP__R2_GGT,MP__R2_GGG_1,MP__R2_GGG_2,MP__R2_QQQ
     $ ,MP__R2_QQB,MP__UV_GQQG_1EPS,MP__UV_GQQB_1EPS,MP__UV_BMASS_1EPS
     $ ,MP__UVWFCT_B_0_1EPS,MP__R2_DDA,MP__R2_DDZ_V2,MP__R2_DDZ_V3
     $ ,MP__R2_GGZUP,MP__R2_GGZDOWN,MP__UV_GQQB,MP__UV_GQQT
     $ ,MP__UV_BMASS,MP__UVWFCT_B_0,MP__GC_1,MP__GC_3,MP__GC_21
     $ ,MP__GC_22,MP__GC_23,MP__GC_24


