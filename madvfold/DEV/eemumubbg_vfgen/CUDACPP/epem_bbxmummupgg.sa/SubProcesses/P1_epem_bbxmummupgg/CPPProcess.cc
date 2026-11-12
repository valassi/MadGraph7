// Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
// Created by: J. Alwall (Oct 2010) for the MG5aMC CPP backend.
//==========================================================================
// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: S. Roiser (Feb 2020) for the MG5aMC CUDACPP plugin.
// Further modified by: S. Hageboeck, O. Mattelaer, S. Roiser, J. Teig, A. Valassi, Z. Wettersten (2020-2024) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2026) for the MadVfold project (arXiv:2609.23500).
//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
//  MadGraph7 v. 0.2.1, 2026-09-24
//  By the MadGraph7 Development Team
//  Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "CPPProcess.h"

#include "mgOnGpuConfig.h"

#include "GpuRuntime.h"
#include "HelAmps_sm.h"
#include "MemoryAccessAmplitudes.h"
#include "MemoryAccessChannelIds.h"
#include "MemoryAccessCouplings.h"
#include "MemoryAccessCouplingsFixed.h"
#include "MemoryAccessGs.h"
#include "MemoryAccessMatrixElements.h"
#include "MemoryAccessMomenta.h"
#include "MemoryAccessWavefunctions.h"
#include "color_sum.h"

#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
#include "MemoryAccessDenominators.h"
#include "MemoryAccessNumerators.h"
#include "coloramps.h"
#endif

#include <algorithm>
#include <array>
#include <cfenv>  // for feenableexcept, fegetexcept and FE_XXX
#include <cfloat> // for FLT_MIN
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>

// Test ncu metrics for CUDA thread divergence
#undef MGONGPU_TEST_DIVERGENCE
//#define MGONGPU_TEST_DIVERGENCE 1

//--------------------------------------------------------------------------

// Enable FPE traps (see #701, #733, #831 - except on MacOS where feenableexcept is not defined #730)
// [NB1: Fortran default is -ffpe-trap=none, i.e. FPE traps are not enabled, https://gcc.gnu.org/onlinedocs/gfortran/Debugging-Options.html]
// [NB2: Fortran default is -ffpe-summary=invalid,zero,overflow,underflow,denormal, i.e. warn at the end on STOP]
inline void
fpeEnable()
{
  static bool first = true; // FIXME: quick and dirty hack to do this only once (can be removed when separate C++/CUDA builds are implemented)
  if( !first ) return;
  first = false;
#ifndef __APPLE__ // on MacOS feenableexcept is not defined #730
  //int fpes = fegetexcept();
  //std::cout << "fpeEnable: analyse fegetexcept()=" << fpes << std::endl;
  //std::cout << "fpeEnable:     FE_DIVBYZERO is" << ( ( fpes & FE_DIVBYZERO ) ? " " : " NOT " ) << "enabled" << std::endl;
  //std::cout << "fpeEnable:     FE_INEXACT is" << ( ( fpes & FE_INEXACT ) ? " " : " NOT " ) << "enabled" << std::endl;
  //std::cout << "fpeEnable:     FE_INVALID is" << ( ( fpes & FE_INVALID ) ? " " : " NOT " ) << "enabled" << std::endl;
  //std::cout << "fpeEnable:     FE_OVERFLOW is" << ( ( fpes & FE_OVERFLOW ) ? " " : " NOT " ) << "enabled" << std::endl;
  //std::cout << "fpeEnable:     FE_UNDERFLOW is" << ( ( fpes & FE_UNDERFLOW ) ? " " : " NOT " ) << "enabled" << std::endl;
  constexpr bool enableFPE = true; // this is hardcoded and no longer controlled by getenv( "CUDACPP_RUNTIME_ENABLEFPE" )
  if( enableFPE )
  {
    std::cout << "INFO: The following Floating Point Exceptions will cause SIGFPE program aborts: FE_DIVBYZERO, FE_INVALID, FE_OVERFLOW" << std::endl;
    feenableexcept( FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW ); // new strategy #831 (do not enable FE_UNDERFLOW)
    //fpes = fegetexcept();
    //std::cout << "fpeEnable: analyse fegetexcept()=" << fpes << std::endl;
    //std::cout << "fpeEnable:     FE_DIVBYZERO is" << ( ( fpes & FE_DIVBYZERO ) ? " " : " NOT " ) << "enabled" << std::endl;
    //std::cout << "fpeEnable:     FE_INEXACT is" << ( ( fpes & FE_INEXACT ) ? " " : " NOT " ) << "enabled" << std::endl;
    //std::cout << "fpeEnable:     FE_INVALID is" << ( ( fpes & FE_INVALID ) ? " " : " NOT " ) << "enabled" << std::endl;
    //std::cout << "fpeEnable:     FE_OVERFLOW is" << ( ( fpes & FE_OVERFLOW ) ? " " : " NOT " ) << "enabled" << std::endl;
    //std::cout << "fpeEnable:     FE_UNDERFLOW is" << ( ( fpes & FE_UNDERFLOW ) ? " " : " NOT " ) << "enabled" << std::endl;
  }
  else
  {
    //std::cout << "INFO: Do not enable SIGFPE traps for Floating Point Exceptions" << std::endl;
  }
#else
  //std::cout << "INFO: Keep default SIGFPE settings because feenableexcept is not available on MacOS" << std::endl;
#endif
}

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: e+ e- > b b~ mu- mu+ g g WEIGHTED<=10 @1

#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  constexpr int nw6 = CPPProcess::nw6;       // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)
  constexpr int npar = CPPProcess::npar;     // #particles in total (external = initial + final): e.g. 4 for e+ e- -> mu+ mu-
  constexpr int ncomb = CPPProcess::ncomb;   // #helicity combinations: e.g. 16 for e+ e- -> mu+ mu- (2**4 = fermion spin up/down ** npar)
  constexpr int ncolor = CPPProcess::ncolor; // the number of leading colors

  // [NB: I am currently unable to get the right value of nwf in CPPProcess.h - will hardcode it in CPPProcess.cc instead (#644)]
  //using CPPProcess::nwf; // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)

  using Parameters_sm_dependentCouplings::ndcoup;   // #couplings that vary event by event (depend on running alphas QCD)
  using Parameters_sm_independentCouplings::nicoup; // #couplings that are fixed for all events (do not depend on running alphas QCD)

  // The number of SIMD vectors of events processed by calculate_jamps
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
  constexpr int nParity = 2;
#else
  constexpr int nParity = 1;
#endif

  // Physics parameters (masses, coupling, etc...)
  // For CUDA performance, hardcoded constexpr's would be better: fewer registers and a tiny throughput increase
  // However, physics parameters are user-defined through card files: use CUDA constant memory instead (issue #39)
  // [NB if hardcoded parameters are used, it's better to define them here to avoid silent shadowing (issue #263)]
  constexpr int nIPD = 5; // SM independent parameters used in this CPPProcess.cc (FIXME? rename as sm_IndepParam?)
  // Note: in the Python code generator, nIPD == nparam, while nIPC <= nicoup, because (see #823)
  // nIPC may vary from one P*/CPPProcess.cc to another, while nicoup is defined in src/Param.h and is common to all P*
  constexpr int nIPC = 7; // SM independent couplings used in this CPPProcess.cc (FIXME? rename as sm_IndepCoupl?)
  static_assert( nIPC <= nicoup );
  static_assert( nIPD >= 0 ); // Hack to avoid build warnings when nIPD==0 is unused
  static_assert( nIPC >= 0 ); // Hack to avoid build warnings when nIPC==0 is unused
#ifdef MGONGPU_HARDCODE_PARAM
  __device__ const fptype cIPD[nIPD] = { (fptype)Parameters_sm::mdl_MB, (fptype)Parameters_sm::mdl_MZ, (fptype)Parameters_sm::mdl_WZ, (fptype)Parameters_sm::mdl_MH, (fptype)Parameters_sm::mdl_WH };
  __device__ const fptype cIPC[nIPC * 2] = { (fptype)Parameters_sm::GC_3.real(), (fptype)Parameters_sm::GC_3.imag(), (fptype)Parameters_sm::GC_1.real(), (fptype)Parameters_sm::GC_1.imag(), (fptype)Parameters_sm::GC_50.real(), (fptype)Parameters_sm::GC_50.imag(), (fptype)Parameters_sm::GC_59.real(), (fptype)Parameters_sm::GC_59.imag(), (fptype)Parameters_sm::GC_58.real(), (fptype)Parameters_sm::GC_58.imag(), (fptype)Parameters_sm::GC_81.real(), (fptype)Parameters_sm::GC_81.imag(), (fptype)Parameters_sm::GC_83.real(), (fptype)Parameters_sm::GC_83.imag() };
#else
#ifdef MGONGPUCPP_GPUIMPL
  __device__ __constant__ fptype cIPD[nIPD];
  __device__ __constant__ fptype cIPC[nIPC * 2];
#else
  static fptype cIPD[nIPD];
  static fptype cIPC[nIPC * 2];
#endif
#endif

  // AV Jan 2024 (PR #625): this ugly #define was the only way I found to avoid creating arrays[nBsm] in CPPProcess.cc if nBsm is 0
  // The problem is that nBsm is determined when generating Parameters.h, which happens after CPPProcess.cc has already been generated
  // For simplicity, keep this code hardcoded also for SM processes (a nullptr is needed as in the case nBsm == 0)
#ifdef MGONGPUCPP_NBSMINDEPPARAM_GT_0
#ifdef MGONGPU_HARDCODE_PARAM
  __device__ const double* bsmIndepParam = Parameters_sm::mdl_bsmIndepParam;
#else
#ifdef MGONGPUCPP_GPUIMPL
  __device__ __constant__ double bsmIndepParam[Parameters_sm::nBsmIndepParam];
#else
  static double bsmIndepParam[Parameters_sm::nBsmIndepParam];
#endif
#endif
#else
#ifdef MGONGPU_HARDCODE_PARAM
  __device__ const double* bsmIndepParam = nullptr;
#else
#ifdef MGONGPUCPP_GPUIMPL
  __device__ __constant__ double* bsmIndepParam = nullptr;
#else
  static double* bsmIndepParam = nullptr;
#endif
#endif
#endif

  // Helicity combinations (and filtering of "good" helicity combinations)
#ifdef MGONGPUCPP_GPUIMPL
  __device__ __constant__ short cHel[ncomb][npar];
  __device__ __constant__ int dcNGoodHel;
  __device__ __constant__ int dcGoodHel[ncomb];
#else
  static short cHel[ncomb][npar];
#endif
  static int cNGoodHel;
  static int cGoodHel[ncomb];

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  class DeviceAccessJamp2
  {
  public:
    static __device__ inline fptype&
    kernelAccessIcol( fptype* buffer, const int icol )
    {
      const int nevt = gridDim.x * blockDim.x;
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;
      return buffer[icol * nevt + ievt];
    }
    static __device__ inline const fptype&
    kernelAccessIcolConst( const fptype* buffer, const int icol )
    {
      const int nevt = gridDim.x * blockDim.x;
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;
      return buffer[icol * nevt + ievt];
    }
  };
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  __device__ INLINE unsigned int
  gpu_channelId( const unsigned int* allChannelIds )
  {
    unsigned int channelId = 0; // disable multichannel single-diagram enhancement unless allChannelIds != nullptr
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    using CID_ACCESS = DeviceAccessChannelIds; // non-trivial access: buffer includes all events
    // SCALAR channelId for the current event (CUDA)
    if( allChannelIds != nullptr )
    {
      const unsigned int* channelIds = allChannelIds;                            // fix #899 (distinguish channelIds and allChannelIds)
      const uint_sv channelIds_sv = CID_ACCESS::kernelAccessConst( channelIds ); // fix #895 (compute this only once for all diagrams)
      // NB: channelIds_sv is a scalar in CUDA
      channelId = channelIds_sv;
      assert( channelId > 0 ); // SANITY CHECK: scalar channelId must be > 0 if multichannel is enabled (allChannelIds != nullptr)
    }
#endif
    return channelId;
  }
#endif

  //--------------------------------------------------------------------------

  // Evaluate QCD partial amplitudes jamps for this given helicity from Feynman diagrams
  // Also compute running sums over helicities adding jamp2, numerator, denominator
  // (NB: this function no longer handles matrix elements as the color sum has now been moved to a separate function/kernel)
  // In CUDA, this function processes a single event
  // ** NB1: NEW Nov2024! In CUDA this is now a kernel function (it used to be a device function)
  // ** NB2: NEW Nov2024! in CUDA this now takes a channelId array as input (it used to take a scalar channelId as input)
  // In C++, this function processes a single event "page" or SIMD vector (or for two in "mixed" precision mode, nParity=2)
  // *** NB: in C++, calculate_jamps accepts a SCALAR channelId because it is GUARANTEED that all events in a SIMD vector have the same channelId #898
  __global__ void /* clang-format off */
  calculate_jamps( int ihel,
                   const fptype* allmomenta,          // input: momenta[nevt*npar*4]
                   const fptype* allcouplings,        // input: couplings[nevt*ndcoup*2]
#ifdef MGONGPUCPP_GPUIMPL
                   fptype* allJamps,                  // output: jamp[2*ncolor*nevt] buffer for one helicity _within a super-buffer for dcNGoodHel helicities_
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                   const unsigned int* allChannelIds, // input: multichannel channelIds[nevt] (1 to #diagrams); nullptr to disable SDE (#899/#911)
                   fptype* allNumerators,             // input/output: multichannel numerators[nevt], add helicity ihel
                   fptype* allDenominators,           // input/output: multichannel denominators[nevt], add helicity ihel
                   fptype* colAllJamp2s,              // output: allJamp2s[ncolor][nevt] super-buffer, sum over col/hel (nullptr to disable)
#endif
                   const int nevt                     // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#else
                   cxtype_sv* allJamp_sv,             // output: jamp_sv[ncolor] (float/double) or jamp_sv[2*ncolor] (mixed) for this helicity
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                   const unsigned int channelId,      // input: SCALAR channelId (1 to #diagrams, 0 to disable SDE) for this event or SIMD vector
                   fptype* allNumerators,             // input/output: multichannel numerators[nevt], add helicity ihel
                   fptype* allDenominators,           // input/output: multichannel denominators[nevt], add helicity ihel
                   fptype_sv* jamp2_sv,               // output: jamp2[nParity][ncolor][neppV] for color choice (nullptr if disabled)
#endif
                   const int ievt00                   // input: first event number in current C++ event page (for CUDA, ievt depends on threadid)
#endif
                   )
  //ALWAYS_INLINE // attributes are not permitted in a function definition
  {
#ifdef MGONGPUCPP_GPUIMPL
    using namespace mg5amcGpu;
    using M_ACCESS = DeviceAccessMomenta;         // non-trivial access: buffer includes all events
    using W_ACCESS = DeviceAccessWavefunctions;   // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
    using A_ACCESS = DeviceAccessAmplitudes;      // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
    using CD_ACCESS = DeviceAccessCouplings;      // non-trivial access (dependent couplings): buffer includes all events
    using CI_ACCESS = DeviceAccessCouplingsFixed; // TRIVIAL access (independent couplings): buffer for one event
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    using NUM_ACCESS = DeviceAccessNumerators;    // non-trivial access: buffer includes all events
    using DEN_ACCESS = DeviceAccessDenominators;  // non-trivial access: buffer includes all events
#endif
#else
    using namespace mg5amcCpu;
    using M_ACCESS = HostAccessMomenta;         // non-trivial access: buffer includes all events
    using W_ACCESS = HostAccessWavefunctions;   // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
    using A_ACCESS = HostAccessAmplitudes;      // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
    using CD_ACCESS = HostAccessCouplings;      // non-trivial access (dependent couplings): buffer includes all events
    using CI_ACCESS = HostAccessCouplingsFixed; // TRIVIAL access (independent couplings): buffer for one event
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    using NUM_ACCESS = HostAccessNumerators;    // non-trivial access: buffer includes all events
    using DEN_ACCESS = HostAccessDenominators;  // non-trivial access: buffer includes all events
#endif
#endif
    mgDebug( 0, __FUNCTION__ );
    //bool debug = true;
#ifndef MGONGPUCPP_GPUIMPL
    //debug = ( ievt00 >= 64 && ievt00 < 80 && ihel == 3 ); // example: debug #831
    //if( debug ) printf( "calculate_jamps: ievt00=%d ihel=%2d\n", ievt00, ihel );
#else
    //const int ievt = blockDim.x * blockIdx.x + threadIdx.x;
    //debug = ( ievt == 0 );
    //if( debug ) printf( "calculate_jamps: ievt=%6d ihel=%2d\n", ievt, ihel );
#endif /* clang-format on */

    // The variable nwf (which is specific to each P1 subdirectory, #644) is only used here
    // It is hardcoded here because various attempts to hardcode it in CPPProcess.h at generation time gave the wrong result...
    static const int nwf = 40; // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)

    // Local TEMPORARY variables for a subset of Feynman diagrams in the given CUDA event (ievt) or C++ event page (ipagV)
    // [NB these variables are reused several times (and re-initialised each time) within the same event or event page]
    // ** NB: in other words, amplitudes and wavefunctions still have TRIVIAL ACCESS: there is currently no need
    // ** NB: to have large memory structurs for wavefunctions/amplitudes in all events (no kernel splitting yet)!
    //MemoryBufferWavefunctions w_buffer[nwf]{ neppV };
    cxtype_sv w_sv[nwf][nw6]; // particle wavefunctions within Feynman diagrams (nw6 is often 6, the dimension of spin 1/2 or spin 1 particles)
    cxtype_sv amp_sv[1];      // invariant amplitude for one given Feynman diagram

    // Proof of concept for using fptype* in the interface
    fptype* w_fp[nwf];
    for( int iwf = 0; iwf < nwf; iwf++ ) w_fp[iwf] = reinterpret_cast<fptype*>( w_sv[iwf] );
    fptype* amp_fp;
    amp_fp = reinterpret_cast<fptype*>( amp_sv );

    // Local variables for the given CUDA event (ievt) or C++ event page (ipagV)
    // [jamp: sum (for one event or event page) of the invariant amplitudes for all Feynman diagrams in a given color combination]
    cxtype_sv jamp_sv[ncolor] = {}; // all zeros (NB: vector cxtype_v IS initialized to 0, but scalar cxtype is NOT, if "= {}" is missing!)

    // === Calculate wavefunctions and amplitudes for all diagrams in all processes         ===
    // === (for one event in CUDA, for one - or two in mixed mode - SIMD event pages in C++ ===

    // START LOOP ON IPARITY
    for( int iParity = 0; iParity < nParity; ++iParity )
    {
#ifndef MGONGPUCPP_GPUIMPL
      const int ievt0 = ievt00 + iParity * neppV;
#endif
      //constexpr size_t nxcoup = ndcoup + nicoup; // both dependent and independent couplings (BUG #823)
      constexpr size_t nxcoup = ndcoup + nIPC; // both dependent and independent couplings (FIX #823)
      const fptype* allCOUPs[nxcoup];
#ifdef __CUDACC__ // this must be __CUDACC__ (not MGONGPUCPP_GPUIMPL)
#pragma nv_diagnostic push
#pragma nv_diag_suppress 186 // e.g. <<warning #186-D: pointless comparison of unsigned integer with zero>>
#endif
      for( size_t idcoup = 0; idcoup < ndcoup; idcoup++ )
        allCOUPs[idcoup] = CD_ACCESS::idcoupAccessBufferConst( allcouplings, idcoup ); // dependent couplings, vary event-by-event
      //for( size_t iicoup = 0; iicoup < nicoup; iicoup++ )                             // BUG #823
      for( size_t iicoup = 0; iicoup < nIPC; iicoup++ )                                 // FIX #823
        allCOUPs[ndcoup + iicoup] = CI_ACCESS::iicoupAccessBufferConst( cIPC, iicoup ); // independent couplings, fixed for all events
#ifdef MGONGPUCPP_GPUIMPL
#ifdef __CUDACC__ // this must be __CUDACC__ (not MGONGPUCPP_GPUIMPL)
#pragma nv_diagnostic pop
#endif
      // CUDA kernels take input/output buffers with momenta/MEs for all events
      const fptype* momenta = allmomenta;
      const fptype* COUPs[nxcoup];
      for( size_t ixcoup = 0; ixcoup < nxcoup; ixcoup++ ) COUPs[ixcoup] = allCOUPs[ixcoup];
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      fptype* numerators = allNumerators;
      fptype* denominators = allDenominators;
#endif
#else
      // C++ kernels take input/output buffers with momenta/MEs for one specific event (the first in the current event page)
      const fptype* momenta = M_ACCESS::ieventAccessRecordConst( allmomenta, ievt0 );
      const fptype* COUPs[nxcoup];
      for( size_t idcoup = 0; idcoup < ndcoup; idcoup++ )
        COUPs[idcoup] = CD_ACCESS::ieventAccessRecordConst( allCOUPs[idcoup], ievt0 ); // dependent couplings, vary event-by-event
      //for( size_t iicoup = 0; iicoup < nicoup; iicoup++ ) // BUG #823
      for( size_t iicoup = 0; iicoup < nIPC; iicoup++ )     // FIX #823
        COUPs[ndcoup + iicoup] = allCOUPs[ndcoup + iicoup]; // independent couplings, fixed for all events
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      fptype* numerators = NUM_ACCESS::ieventAccessRecord( allNumerators, ievt0 );
      fptype* denominators = DEN_ACCESS::ieventAccessRecord( allDenominators, ievt0 );
#endif
#endif

      // Reset color flows (reset jamp_sv) at the beginning of a new event or event page
      for( int i = 0; i < ncolor; i++ ) { jamp_sv[i] = cxzero_sv(); }

#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
#ifdef MGONGPUCPP_GPUIMPL
      // SCALAR channelId for the current event (CUDA)
      unsigned int channelId = gpu_channelId( allChannelIds );
#endif
      // Numerators and denominators for the current event (CUDA) or SIMD event page (C++)
      fptype_sv& numerators_sv = NUM_ACCESS::kernelAccess( numerators );
      fptype_sv& denominators_sv = DEN_ACCESS::kernelAccess( denominators );
#endif

      // *** DIAGRAM 1 OF 256 ***

      // Wavefunction(s) for diagram number 1
      oxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][0], -1, w_fp[0], 0 );

      ixxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][1], +1, w_fp[1], 1 );

      oxxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][2], +1, w_fp[2], 2 );

      ixxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][3], -1, w_fp[3], 3 );

      oxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][4], +1, w_fp[4], 4 );

      ixxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][5], -1, w_fp[5], 5 );

      vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][6], +1, w_fp[6], 6 );

      vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][7], +1, w_fp[7], 7 );

      FFV1P0_3<W_ACCESS, CI_ACCESS>( w_fp[1], w_fp[0], COUPs[ndcoup + 0], 1.0, 0., 0., w_fp[8] );
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[6], COUPs[0], 1.0, cIPD[0], 0., w_fp[9] );
      FFV1_1<W_ACCESS, CI_ACCESS>( w_fp[4], w_fp[8], COUPs[ndcoup + 0], 1.0, 0., 0., w_fp[10] );
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[9], w_fp[7], COUPs[0], 1.0, cIPD[0], 0., w_fp[11] );
      FFV1P0_3<W_ACCESS, CI_ACCESS>( w_fp[5], w_fp[10], COUPs[ndcoup + 0], 1.0, 0., 0., w_fp[12] );

      // Amplitude(s) for diagram number 1
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[11], w_fp[12], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 2 OF 256 ***

      // Wavefunction(s) for diagram number 2
      FFV2_4_3<W_ACCESS, CI_ACCESS>( w_fp[5], w_fp[10], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, cIPD[1], cIPD[2], w_fp[13] );

      // Amplitude(s) for diagram number 2
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[11], w_fp[13], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 3 OF 256 ***

      // Wavefunction(s) for diagram number 3
      FFV1_2<W_ACCESS, CI_ACCESS>( w_fp[5], w_fp[8], COUPs[ndcoup + 0], 1.0, 0., 0., w_fp[14] );
      FFV1P0_3<W_ACCESS, CI_ACCESS>( w_fp[14], w_fp[4], COUPs[ndcoup + 0], 1.0, 0., 0., w_fp[15] );

      // Amplitude(s) for diagram number 3
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[11], w_fp[15], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 4 OF 256 ***

      // Wavefunction(s) for diagram number 4
      FFV2_4_3<W_ACCESS, CI_ACCESS>( w_fp[14], w_fp[4], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, cIPD[1], cIPD[2], w_fp[16] );

      // Amplitude(s) for diagram number 4
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[11], w_fp[16], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 5 OF 256 ***

      // Wavefunction(s) for diagram number 5
      FFV2_4_3<W_ACCESS, CI_ACCESS>( w_fp[1], w_fp[0], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, cIPD[1], cIPD[2], w_fp[17] );
      FFV2_4_1<W_ACCESS, CI_ACCESS>( w_fp[4], w_fp[17], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, 0., 0., w_fp[18] );
      FFV1P0_3<W_ACCESS, CI_ACCESS>( w_fp[5], w_fp[18], COUPs[ndcoup + 0], 1.0, 0., 0., w_fp[19] );

      // Amplitude(s) for diagram number 5
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[11], w_fp[19], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 6 OF 256 ***

      // Wavefunction(s) for diagram number 6
      FFV2_4_3<W_ACCESS, CI_ACCESS>( w_fp[5], w_fp[18], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, cIPD[1], cIPD[2], w_fp[20] );

      // Amplitude(s) for diagram number 6
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[11], w_fp[20], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 7 OF 256 ***

      // Wavefunction(s) for diagram number 7
      FFV2_4_2<W_ACCESS, CI_ACCESS>( w_fp[5], w_fp[17], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, 0., 0., w_fp[21] );
      FFV1P0_3<W_ACCESS, CI_ACCESS>( w_fp[21], w_fp[4], COUPs[ndcoup + 0], 1.0, 0., 0., w_fp[22] );

      // Amplitude(s) for diagram number 7
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[11], w_fp[22], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 8 OF 256 ***

      // Wavefunction(s) for diagram number 8
      FFV2_4_3<W_ACCESS, CI_ACCESS>( w_fp[21], w_fp[4], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, cIPD[1], cIPD[2], w_fp[23] );

      // Amplitude(s) for diagram number 8
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[11], w_fp[23], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 9 OF 256 ***

      // Wavefunction(s) for diagram number 9
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[7], COUPs[0], 1.0, cIPD[0], 0., w_fp[24] );
      FFV1P0_3<W_ACCESS, CI_ACCESS>( w_fp[24], w_fp[9], COUPs[ndcoup + 1], 1.0, 0., 0., w_fp[25] );

      // Amplitude(s) for diagram number 9
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[10], w_fp[25], COUPs[ndcoup + 0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 10 OF 256 ***

      // Wavefunction(s) for diagram number 10
      FFV2_3_3<W_ACCESS, CI_ACCESS>( w_fp[24], w_fp[9], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, cIPD[1], cIPD[2], w_fp[26] );

      // Amplitude(s) for diagram number 10
      FFV2_4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[10], w_fp[26], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 11 OF 256 ***

      // Wavefunction(s) for diagram number 11
      // (none)

      // Amplitude(s) for diagram number 11
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[14], w_fp[4], w_fp[25], COUPs[ndcoup + 0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 12 OF 256 ***

      // Wavefunction(s) for diagram number 12
      // (none)

      // Amplitude(s) for diagram number 12
      FFV2_4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[14], w_fp[4], w_fp[26], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 13 OF 256 ***

      // Wavefunction(s) for diagram number 13
      // (none)

      // Amplitude(s) for diagram number 13
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[18], w_fp[25], COUPs[ndcoup + 0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 14 OF 256 ***

      // Wavefunction(s) for diagram number 14
      // (none)

      // Amplitude(s) for diagram number 14
      FFV2_4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[18], w_fp[26], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 15 OF 256 ***

      // Wavefunction(s) for diagram number 15
      // (none)

      // Amplitude(s) for diagram number 15
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[4], w_fp[25], COUPs[ndcoup + 0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 16 OF 256 ***

      // Wavefunction(s) for diagram number 16
      // (none)

      // Amplitude(s) for diagram number 16
      FFV2_4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[4], w_fp[26], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 17 OF 256 ***

      // Wavefunction(s) for diagram number 17
      FFV1P0_3<W_ACCESS, CI_ACCESS>( w_fp[5], w_fp[4], COUPs[ndcoup + 0], 1.0, 0., 0., w_fp[27] );
      FFV1_1<W_ACCESS, CI_ACCESS>( w_fp[9], w_fp[8], COUPs[ndcoup + 1], 1.0, cIPD[0], 0., w_fp[28] );

      // Amplitude(s) for diagram number 17
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[24], w_fp[28], w_fp[27], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 18 OF 256 ***

      // Wavefunction(s) for diagram number 18
      FFV1_2<W_ACCESS, CI_ACCESS>( w_fp[24], w_fp[8], COUPs[ndcoup + 1], 1.0, cIPD[0], 0., w_fp[29] );

      // Amplitude(s) for diagram number 18
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[29], w_fp[9], w_fp[27], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 19 OF 256 ***

      // Wavefunction(s) for diagram number 19
      FFV2_4_3<W_ACCESS, CI_ACCESS>( w_fp[5], w_fp[4], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, cIPD[1], cIPD[2], w_fp[30] );

      // Amplitude(s) for diagram number 19
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[24], w_fp[28], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 20 OF 256 ***

      // Wavefunction(s) for diagram number 20
      // (none)

      // Amplitude(s) for diagram number 20
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[29], w_fp[9], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 21 OF 256 ***

      // Wavefunction(s) for diagram number 21
      FFV2_3_1<W_ACCESS, CI_ACCESS>( w_fp[9], w_fp[17], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, cIPD[0], 0., w_fp[31] );

      // Amplitude(s) for diagram number 21
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[24], w_fp[31], w_fp[27], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 22 OF 256 ***

      // Wavefunction(s) for diagram number 22
      FFV2_3_2<W_ACCESS, CI_ACCESS>( w_fp[24], w_fp[17], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, cIPD[0], 0., w_fp[32] );

      // Amplitude(s) for diagram number 22
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[32], w_fp[9], w_fp[27], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 23 OF 256 ***

      // Wavefunction(s) for diagram number 23
      // (none)

      // Amplitude(s) for diagram number 23
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[24], w_fp[31], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 24 OF 256 ***

      // Wavefunction(s) for diagram number 24
      // (none)

      // Amplitude(s) for diagram number 24
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[32], w_fp[9], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 25 OF 256 ***

      // Wavefunction(s) for diagram number 25
      VVS1_3<W_ACCESS, CI_ACCESS>( w_fp[17], w_fp[30], COUPs[ndcoup + 5], 1.0, cIPD[3], cIPD[4], w_fp[33] );

      // Amplitude(s) for diagram number 25
      FFS4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[24], w_fp[9], w_fp[33], COUPs[ndcoup + 6], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 26 OF 256 ***

      // Wavefunction(s) for diagram number 26
      FFV1_2<W_ACCESS, CI_ACCESS>( w_fp[3], w_fp[27], COUPs[ndcoup + 1], 1.0, cIPD[0], 0., w_fp[34] );

      // Amplitude(s) for diagram number 26
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[28], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 27 OF 256 ***

      // Wavefunction(s) for diagram number 27
      FFV1_2<W_ACCESS, CI_ACCESS>( w_fp[3], w_fp[8], COUPs[ndcoup + 1], 1.0, cIPD[0], 0., w_fp[35] );
      FFV1_1<W_ACCESS, CI_ACCESS>( w_fp[9], w_fp[27], COUPs[ndcoup + 1], 1.0, cIPD[0], 0., w_fp[36] );

      // Amplitude(s) for diagram number 27
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[35], w_fp[36], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 28 OF 256 ***

      // Wavefunction(s) for diagram number 28
      // (none)

      // Amplitude(s) for diagram number 28
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[35], w_fp[11], w_fp[27], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 29 OF 256 ***

      // Wavefunction(s) for diagram number 29
      // (none)

      // Amplitude(s) for diagram number 29
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[34], w_fp[11], w_fp[8], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 30 OF 256 ***

      // Wavefunction(s) for diagram number 30
      FFV2_3_2<W_ACCESS, CI_ACCESS>( w_fp[3], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, cIPD[0], 0., w_fp[37] );

      // Amplitude(s) for diagram number 30
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[37], w_fp[28], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 31 OF 256 ***

      // Wavefunction(s) for diagram number 31
      FFV2_3_1<W_ACCESS, CI_ACCESS>( w_fp[9], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, cIPD[0], 0., w_fp[28] );

      // Amplitude(s) for diagram number 31
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[35], w_fp[28], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 32 OF 256 ***

      // Wavefunction(s) for diagram number 32
      // (none)

      // Amplitude(s) for diagram number 32
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[35], w_fp[11], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 33 OF 256 ***

      // Wavefunction(s) for diagram number 33
      // (none)

      // Amplitude(s) for diagram number 33
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[37], w_fp[11], w_fp[8], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 34 OF 256 ***

      // Wavefunction(s) for diagram number 34
      // (none)

      // Amplitude(s) for diagram number 34
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[31], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 35 OF 256 ***

      // Wavefunction(s) for diagram number 35
      FFV2_3_2<W_ACCESS, CI_ACCESS>( w_fp[3], w_fp[17], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, cIPD[0], 0., w_fp[9] );

      // Amplitude(s) for diagram number 35
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[36], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 36 OF 256 ***

      // Wavefunction(s) for diagram number 36
      // (none)

      // Amplitude(s) for diagram number 36
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[9], w_fp[11], w_fp[27], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 37 OF 256 ***

      // Wavefunction(s) for diagram number 37
      // (none)

      // Amplitude(s) for diagram number 37
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[34], w_fp[11], w_fp[17], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 38 OF 256 ***

      // Wavefunction(s) for diagram number 38
      // (none)

      // Amplitude(s) for diagram number 38
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[37], w_fp[31], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 39 OF 256 ***

      // Wavefunction(s) for diagram number 39
      // (none)

      // Amplitude(s) for diagram number 39
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[28], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 40 OF 256 ***

      // Wavefunction(s) for diagram number 40
      // (none)

      // Amplitude(s) for diagram number 40
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[9], w_fp[11], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 41 OF 256 ***

      // Wavefunction(s) for diagram number 41
      // (none)

      // Amplitude(s) for diagram number 41
      FFS4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[11], w_fp[33], COUPs[ndcoup + 6], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 42 OF 256 ***

      // Wavefunction(s) for diagram number 42
      // (none)

      // Amplitude(s) for diagram number 42
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[37], w_fp[11], w_fp[17], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 43 OF 256 ***

      // Wavefunction(s) for diagram number 43
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[7], COUPs[0], 1.0, cIPD[0], 0., w_fp[28] );
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[28], w_fp[6], COUPs[0], 1.0, cIPD[0], 0., w_fp[31] );

      // Amplitude(s) for diagram number 43
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[31], w_fp[12], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 44 OF 256 ***

      // Wavefunction(s) for diagram number 44
      // (none)

      // Amplitude(s) for diagram number 44
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[31], w_fp[13], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 45 OF 256 ***

      // Wavefunction(s) for diagram number 45
      // (none)

      // Amplitude(s) for diagram number 45
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[31], w_fp[15], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 46 OF 256 ***

      // Wavefunction(s) for diagram number 46
      // (none)

      // Amplitude(s) for diagram number 46
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[31], w_fp[16], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 47 OF 256 ***

      // Wavefunction(s) for diagram number 47
      // (none)

      // Amplitude(s) for diagram number 47
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[31], w_fp[19], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 48 OF 256 ***

      // Wavefunction(s) for diagram number 48
      // (none)

      // Amplitude(s) for diagram number 48
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[31], w_fp[20], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 49 OF 256 ***

      // Wavefunction(s) for diagram number 49
      // (none)

      // Amplitude(s) for diagram number 49
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[31], w_fp[22], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 50 OF 256 ***

      // Wavefunction(s) for diagram number 50
      // (none)

      // Amplitude(s) for diagram number 50
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[31], w_fp[23], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 51 OF 256 ***

      // Wavefunction(s) for diagram number 51
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[6], COUPs[0], 1.0, cIPD[0], 0., w_fp[36] );
      FFV1P0_3<W_ACCESS, CI_ACCESS>( w_fp[36], w_fp[28], COUPs[ndcoup + 1], 1.0, 0., 0., w_fp[38] );

      // Amplitude(s) for diagram number 51
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[10], w_fp[38], COUPs[ndcoup + 0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 52 OF 256 ***

      // Wavefunction(s) for diagram number 52
      FFV2_3_3<W_ACCESS, CI_ACCESS>( w_fp[36], w_fp[28], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, cIPD[1], cIPD[2], w_fp[39] );

      // Amplitude(s) for diagram number 52
      FFV2_4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[10], w_fp[39], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 53 OF 256 ***

      // Wavefunction(s) for diagram number 53
      // (none)

      // Amplitude(s) for diagram number 53
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[14], w_fp[4], w_fp[38], COUPs[ndcoup + 0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 54 OF 256 ***

      // Wavefunction(s) for diagram number 54
      // (none)

      // Amplitude(s) for diagram number 54
      FFV2_4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[14], w_fp[4], w_fp[39], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 55 OF 256 ***

      // Wavefunction(s) for diagram number 55
      // (none)

      // Amplitude(s) for diagram number 55
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[18], w_fp[38], COUPs[ndcoup + 0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 56 OF 256 ***

      // Wavefunction(s) for diagram number 56
      // (none)

      // Amplitude(s) for diagram number 56
      FFV2_4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[18], w_fp[39], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 57 OF 256 ***

      // Wavefunction(s) for diagram number 57
      // (none)

      // Amplitude(s) for diagram number 57
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[4], w_fp[38], COUPs[ndcoup + 0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 58 OF 256 ***

      // Wavefunction(s) for diagram number 58
      // (none)

      // Amplitude(s) for diagram number 58
      FFV2_4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[4], w_fp[39], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 59 OF 256 ***

      // Wavefunction(s) for diagram number 59
      FFV1_1<W_ACCESS, CI_ACCESS>( w_fp[28], w_fp[8], COUPs[ndcoup + 1], 1.0, cIPD[0], 0., w_fp[21] );

      // Amplitude(s) for diagram number 59
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[36], w_fp[21], w_fp[27], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 60 OF 256 ***

      // Wavefunction(s) for diagram number 60
      FFV1_2<W_ACCESS, CI_ACCESS>( w_fp[36], w_fp[8], COUPs[ndcoup + 1], 1.0, cIPD[0], 0., w_fp[4] );

      // Amplitude(s) for diagram number 60
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[4], w_fp[28], w_fp[27], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 61 OF 256 ***

      // Wavefunction(s) for diagram number 61
      // (none)

      // Amplitude(s) for diagram number 61
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[36], w_fp[21], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 62 OF 256 ***

      // Wavefunction(s) for diagram number 62
      // (none)

      // Amplitude(s) for diagram number 62
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[4], w_fp[28], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 63 OF 256 ***

      // Wavefunction(s) for diagram number 63
      FFV2_3_1<W_ACCESS, CI_ACCESS>( w_fp[28], w_fp[17], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, cIPD[0], 0., w_fp[18] );

      // Amplitude(s) for diagram number 63
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[36], w_fp[18], w_fp[27], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 64 OF 256 ***

      // Wavefunction(s) for diagram number 64
      FFV2_3_2<W_ACCESS, CI_ACCESS>( w_fp[36], w_fp[17], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, cIPD[0], 0., w_fp[5] );

      // Amplitude(s) for diagram number 64
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[28], w_fp[27], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 65 OF 256 ***

      // Wavefunction(s) for diagram number 65
      // (none)

      // Amplitude(s) for diagram number 65
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[36], w_fp[18], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 66 OF 256 ***

      // Wavefunction(s) for diagram number 66
      // (none)

      // Amplitude(s) for diagram number 66
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[28], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 67 OF 256 ***

      // Wavefunction(s) for diagram number 67
      // (none)

      // Amplitude(s) for diagram number 67
      FFS4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[36], w_fp[28], w_fp[33], COUPs[ndcoup + 6], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 68 OF 256 ***

      // Wavefunction(s) for diagram number 68
      // (none)

      // Amplitude(s) for diagram number 68
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[21], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 69 OF 256 ***

      // Wavefunction(s) for diagram number 69
      FFV1_1<W_ACCESS, CI_ACCESS>( w_fp[28], w_fp[27], COUPs[ndcoup + 1], 1.0, cIPD[0], 0., w_fp[14] );

      // Amplitude(s) for diagram number 69
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[35], w_fp[14], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 70 OF 256 ***

      // Wavefunction(s) for diagram number 70
      // (none)

      // Amplitude(s) for diagram number 70
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[35], w_fp[31], w_fp[27], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 71 OF 256 ***

      // Wavefunction(s) for diagram number 71
      // (none)

      // Amplitude(s) for diagram number 71
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[34], w_fp[31], w_fp[8], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 72 OF 256 ***

      // Wavefunction(s) for diagram number 72
      // (none)

      // Amplitude(s) for diagram number 72
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[37], w_fp[21], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 73 OF 256 ***

      // Wavefunction(s) for diagram number 73
      FFV2_3_1<W_ACCESS, CI_ACCESS>( w_fp[28], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, cIPD[0], 0., w_fp[21] );

      // Amplitude(s) for diagram number 73
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[35], w_fp[21], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 74 OF 256 ***

      // Wavefunction(s) for diagram number 74
      // (none)

      // Amplitude(s) for diagram number 74
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[35], w_fp[31], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 75 OF 256 ***

      // Wavefunction(s) for diagram number 75
      // (none)

      // Amplitude(s) for diagram number 75
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[37], w_fp[31], w_fp[8], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 76 OF 256 ***

      // Wavefunction(s) for diagram number 76
      // (none)

      // Amplitude(s) for diagram number 76
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[18], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 77 OF 256 ***

      // Wavefunction(s) for diagram number 77
      // (none)

      // Amplitude(s) for diagram number 77
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[14], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 78 OF 256 ***

      // Wavefunction(s) for diagram number 78
      // (none)

      // Amplitude(s) for diagram number 78
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[9], w_fp[31], w_fp[27], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 79 OF 256 ***

      // Wavefunction(s) for diagram number 79
      // (none)

      // Amplitude(s) for diagram number 79
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[34], w_fp[31], w_fp[17], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 80 OF 256 ***

      // Wavefunction(s) for diagram number 80
      // (none)

      // Amplitude(s) for diagram number 80
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[37], w_fp[18], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 81 OF 256 ***

      // Wavefunction(s) for diagram number 81
      // (none)

      // Amplitude(s) for diagram number 81
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[21], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 82 OF 256 ***

      // Wavefunction(s) for diagram number 82
      // (none)

      // Amplitude(s) for diagram number 82
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[9], w_fp[31], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 83 OF 256 ***

      // Wavefunction(s) for diagram number 83
      // (none)

      // Amplitude(s) for diagram number 83
      FFS4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[31], w_fp[33], COUPs[ndcoup + 6], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 84 OF 256 ***

      // Wavefunction(s) for diagram number 84
      // (none)

      // Amplitude(s) for diagram number 84
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[37], w_fp[31], w_fp[17], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 85 OF 256 ***

      // Wavefunction(s) for diagram number 85
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[36], w_fp[7], COUPs[0], 1.0, cIPD[0], 0., w_fp[21] );

      // Amplitude(s) for diagram number 85
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[2], w_fp[12], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 86 OF 256 ***

      // Wavefunction(s) for diagram number 86
      // (none)

      // Amplitude(s) for diagram number 86
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[2], w_fp[13], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 87 OF 256 ***

      // Wavefunction(s) for diagram number 87
      // (none)

      // Amplitude(s) for diagram number 87
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[2], w_fp[15], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 88 OF 256 ***

      // Wavefunction(s) for diagram number 88
      // (none)

      // Amplitude(s) for diagram number 88
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[2], w_fp[16], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 89 OF 256 ***

      // Wavefunction(s) for diagram number 89
      // (none)

      // Amplitude(s) for diagram number 89
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[2], w_fp[19], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 90 OF 256 ***

      // Wavefunction(s) for diagram number 90
      // (none)

      // Amplitude(s) for diagram number 90
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[2], w_fp[20], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 91 OF 256 ***

      // Wavefunction(s) for diagram number 91
      // (none)

      // Amplitude(s) for diagram number 91
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[2], w_fp[22], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 92 OF 256 ***

      // Wavefunction(s) for diagram number 92
      // (none)

      // Amplitude(s) for diagram number 92
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[2], w_fp[23], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 93 OF 256 ***

      // Wavefunction(s) for diagram number 93
      FFV1_1<W_ACCESS, CI_ACCESS>( w_fp[2], w_fp[8], COUPs[ndcoup + 1], 1.0, cIPD[0], 0., w_fp[18] );
      FFV1_2<W_ACCESS, CI_ACCESS>( w_fp[36], w_fp[27], COUPs[ndcoup + 1], 1.0, cIPD[0], 0., w_fp[14] );

      // Amplitude(s) for diagram number 93
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[18], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 94 OF 256 ***

      // Wavefunction(s) for diagram number 94
      // (none)

      // Amplitude(s) for diagram number 94
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[18], w_fp[27], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 95 OF 256 ***

      // Wavefunction(s) for diagram number 95
      FFV1_1<W_ACCESS, CI_ACCESS>( w_fp[2], w_fp[27], COUPs[ndcoup + 1], 1.0, cIPD[0], 0., w_fp[28] );

      // Amplitude(s) for diagram number 95
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[28], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 96 OF 256 ***

      // Wavefunction(s) for diagram number 96
      // (none)

      // Amplitude(s) for diagram number 96
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[28], w_fp[8], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 97 OF 256 ***

      // Wavefunction(s) for diagram number 97
      FFV2_3_2<W_ACCESS, CI_ACCESS>( w_fp[36], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, cIPD[0], 0., w_fp[10] );

      // Amplitude(s) for diagram number 97
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[18], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 98 OF 256 ***

      // Wavefunction(s) for diagram number 98
      // (none)

      // Amplitude(s) for diagram number 98
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[18], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 99 OF 256 ***

      // Wavefunction(s) for diagram number 99
      FFV2_3_1<W_ACCESS, CI_ACCESS>( w_fp[2], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, cIPD[0], 0., w_fp[36] );

      // Amplitude(s) for diagram number 99
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[36], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 100 OF 256 ***

      // Wavefunction(s) for diagram number 100
      // (none)

      // Amplitude(s) for diagram number 100
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[36], w_fp[8], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 101 OF 256 ***

      // Wavefunction(s) for diagram number 101
      FFV2_3_1<W_ACCESS, CI_ACCESS>( w_fp[2], w_fp[17], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, cIPD[0], 0., w_fp[4] );

      // Amplitude(s) for diagram number 101
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[4], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 102 OF 256 ***

      // Wavefunction(s) for diagram number 102
      // (none)

      // Amplitude(s) for diagram number 102
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[4], w_fp[27], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 103 OF 256 ***

      // Wavefunction(s) for diagram number 103
      // (none)

      // Amplitude(s) for diagram number 103
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[28], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 104 OF 256 ***

      // Wavefunction(s) for diagram number 104
      // (none)

      // Amplitude(s) for diagram number 104
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[28], w_fp[17], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 105 OF 256 ***

      // Wavefunction(s) for diagram number 105
      // (none)

      // Amplitude(s) for diagram number 105
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[4], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 106 OF 256 ***

      // Wavefunction(s) for diagram number 106
      // (none)

      // Amplitude(s) for diagram number 106
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[4], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 107 OF 256 ***

      // Wavefunction(s) for diagram number 107
      // (none)

      // Amplitude(s) for diagram number 107
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[36], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 108 OF 256 ***

      // Wavefunction(s) for diagram number 108
      // (none)

      // Amplitude(s) for diagram number 108
      FFS4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[2], w_fp[33], COUPs[ndcoup + 6], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 109 OF 256 ***

      // Wavefunction(s) for diagram number 109
      // (none)

      // Amplitude(s) for diagram number 109
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[36], w_fp[17], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 110 OF 256 ***

      // Wavefunction(s) for diagram number 110
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[24], w_fp[6], COUPs[0], 1.0, cIPD[0], 0., w_fp[5] );

      // Amplitude(s) for diagram number 110
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[2], w_fp[12], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 111 OF 256 ***

      // Wavefunction(s) for diagram number 111
      // (none)

      // Amplitude(s) for diagram number 111
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[2], w_fp[13], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 112 OF 256 ***

      // Wavefunction(s) for diagram number 112
      // (none)

      // Amplitude(s) for diagram number 112
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[2], w_fp[15], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 113 OF 256 ***

      // Wavefunction(s) for diagram number 113
      // (none)

      // Amplitude(s) for diagram number 113
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[2], w_fp[16], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 114 OF 256 ***

      // Wavefunction(s) for diagram number 114
      // (none)

      // Amplitude(s) for diagram number 114
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[2], w_fp[19], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 115 OF 256 ***

      // Wavefunction(s) for diagram number 115
      // (none)

      // Amplitude(s) for diagram number 115
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[2], w_fp[20], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 116 OF 256 ***

      // Wavefunction(s) for diagram number 116
      // (none)

      // Amplitude(s) for diagram number 116
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[2], w_fp[22], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 117 OF 256 ***

      // Wavefunction(s) for diagram number 117
      // (none)

      // Amplitude(s) for diagram number 117
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[2], w_fp[23], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 118 OF 256 ***

      // Wavefunction(s) for diagram number 118
      FFV1_2<W_ACCESS, CI_ACCESS>( w_fp[24], w_fp[27], COUPs[ndcoup + 1], 1.0, cIPD[0], 0., w_fp[10] );

      // Amplitude(s) for diagram number 118
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[18], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 119 OF 256 ***

      // Wavefunction(s) for diagram number 119
      // (none)

      // Amplitude(s) for diagram number 119
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[18], w_fp[27], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 120 OF 256 ***

      // Wavefunction(s) for diagram number 120
      // (none)

      // Amplitude(s) for diagram number 120
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[28], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 121 OF 256 ***

      // Wavefunction(s) for diagram number 121
      // (none)

      // Amplitude(s) for diagram number 121
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[28], w_fp[8], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 122 OF 256 ***

      // Wavefunction(s) for diagram number 122
      FFV2_3_2<W_ACCESS, CI_ACCESS>( w_fp[24], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, cIPD[0], 0., w_fp[14] );

      // Amplitude(s) for diagram number 122
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[18], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 123 OF 256 ***

      // Wavefunction(s) for diagram number 123
      // (none)

      // Amplitude(s) for diagram number 123
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[18], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 124 OF 256 ***

      // Wavefunction(s) for diagram number 124
      // (none)

      // Amplitude(s) for diagram number 124
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[36], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 125 OF 256 ***

      // Wavefunction(s) for diagram number 125
      // (none)

      // Amplitude(s) for diagram number 125
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[36], w_fp[8], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 126 OF 256 ***

      // Wavefunction(s) for diagram number 126
      // (none)

      // Amplitude(s) for diagram number 126
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[4], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 127 OF 256 ***

      // Wavefunction(s) for diagram number 127
      // (none)

      // Amplitude(s) for diagram number 127
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[4], w_fp[27], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 128 OF 256 ***

      // Wavefunction(s) for diagram number 128
      // (none)

      // Amplitude(s) for diagram number 128
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[32], w_fp[28], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 129 OF 256 ***

      // Wavefunction(s) for diagram number 129
      // (none)

      // Amplitude(s) for diagram number 129
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[28], w_fp[17], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 130 OF 256 ***

      // Wavefunction(s) for diagram number 130
      // (none)

      // Amplitude(s) for diagram number 130
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[4], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 131 OF 256 ***

      // Wavefunction(s) for diagram number 131
      // (none)

      // Amplitude(s) for diagram number 131
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[4], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 132 OF 256 ***

      // Wavefunction(s) for diagram number 132
      // (none)

      // Amplitude(s) for diagram number 132
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[32], w_fp[36], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 133 OF 256 ***

      // Wavefunction(s) for diagram number 133
      // (none)

      // Amplitude(s) for diagram number 133
      FFS4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[2], w_fp[33], COUPs[ndcoup + 6], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 134 OF 256 ***

      // Wavefunction(s) for diagram number 134
      // (none)

      // Amplitude(s) for diagram number 134
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[36], w_fp[17], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 135 OF 256 ***

      // Wavefunction(s) for diagram number 135
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[18], w_fp[6], COUPs[0], 1.0, cIPD[0], 0., w_fp[32] );

      // Amplitude(s) for diagram number 135
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[32], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 136 OF 256 ***

      // Wavefunction(s) for diagram number 136
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[18], w_fp[7], COUPs[0], 1.0, cIPD[0], 0., w_fp[14] );

      // Amplitude(s) for diagram number 136
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[14], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 137 OF 256 ***

      // Wavefunction(s) for diagram number 137
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[35], w_fp[6], COUPs[0], 1.0, cIPD[0], 0., w_fp[10] );

      // Amplitude(s) for diagram number 137
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[28], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 138 OF 256 ***

      // Wavefunction(s) for diagram number 138
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[35], w_fp[7], COUPs[0], 1.0, cIPD[0], 0., w_fp[29] );

      // Amplitude(s) for diagram number 138
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[28], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 139 OF 256 ***

      // Wavefunction(s) for diagram number 139
      // (none)

      // Amplitude(s) for diagram number 139
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[37], w_fp[32], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 140 OF 256 ***

      // Wavefunction(s) for diagram number 140
      // (none)

      // Amplitude(s) for diagram number 140
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[37], w_fp[14], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 141 OF 256 ***

      // Wavefunction(s) for diagram number 141
      // (none)

      // Amplitude(s) for diagram number 141
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[36], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 142 OF 256 ***

      // Wavefunction(s) for diagram number 142
      // (none)

      // Amplitude(s) for diagram number 142
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[36], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 143 OF 256 ***

      // Wavefunction(s) for diagram number 143
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[4], w_fp[6], COUPs[0], 1.0, cIPD[0], 0., w_fp[29] );

      // Amplitude(s) for diagram number 143
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[29], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 144 OF 256 ***

      // Wavefunction(s) for diagram number 144
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[4], w_fp[7], COUPs[0], 1.0, cIPD[0], 0., w_fp[10] );

      // Amplitude(s) for diagram number 144
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[10], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 145 OF 256 ***

      // Wavefunction(s) for diagram number 145
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[9], w_fp[6], COUPs[0], 1.0, cIPD[0], 0., w_fp[14] );

      // Amplitude(s) for diagram number 145
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[28], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 146 OF 256 ***

      // Wavefunction(s) for diagram number 146
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[9], w_fp[7], COUPs[0], 1.0, cIPD[0], 0., w_fp[32] );

      // Amplitude(s) for diagram number 146
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[32], w_fp[28], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 147 OF 256 ***

      // Wavefunction(s) for diagram number 147
      // (none)

      // Amplitude(s) for diagram number 147
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[37], w_fp[29], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 148 OF 256 ***

      // Wavefunction(s) for diagram number 148
      // (none)

      // Amplitude(s) for diagram number 148
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[37], w_fp[10], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 149 OF 256 ***

      // Wavefunction(s) for diagram number 149
      // (none)

      // Amplitude(s) for diagram number 149
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[36], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 150 OF 256 ***

      // Wavefunction(s) for diagram number 150
      // (none)

      // Amplitude(s) for diagram number 150
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[32], w_fp[36], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 151 OF 256 ***

      // Wavefunction(s) for diagram number 151
      VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[6], w_fp[7], COUPs[1], 1.0, 0., 0., w_fp[32] );

      // Amplitude(s) for diagram number 151
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[18], w_fp[32], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 152 OF 256 ***

      // Wavefunction(s) for diagram number 152
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[32], COUPs[0], 1.0, cIPD[0], 0., w_fp[7] );

      // Amplitude(s) for diagram number 152
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[18], w_fp[27], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 153 OF 256 ***

      // Wavefunction(s) for diagram number 153
      // (none)

      // Amplitude(s) for diagram number 153
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[35], w_fp[28], w_fp[32], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 154 OF 256 ***

      // Wavefunction(s) for diagram number 154
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[32], COUPs[0], 1.0, cIPD[0], 0., w_fp[6] );

      // Amplitude(s) for diagram number 154
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[35], w_fp[6], w_fp[27], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 155 OF 256 ***

      // Wavefunction(s) for diagram number 155
      // (none)

      // Amplitude(s) for diagram number 155
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[28], w_fp[8], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 156 OF 256 ***

      // Wavefunction(s) for diagram number 156
      // (none)

      // Amplitude(s) for diagram number 156
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[34], w_fp[6], w_fp[8], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 157 OF 256 ***

      // Wavefunction(s) for diagram number 157
      // (none)

      // Amplitude(s) for diagram number 157
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[37], w_fp[18], w_fp[32], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 158 OF 256 ***

      // Wavefunction(s) for diagram number 158
      // (none)

      // Amplitude(s) for diagram number 158
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[18], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 159 OF 256 ***

      // Wavefunction(s) for diagram number 159
      // (none)

      // Amplitude(s) for diagram number 159
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[35], w_fp[36], w_fp[32], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 160 OF 256 ***

      // Wavefunction(s) for diagram number 160
      // (none)

      // Amplitude(s) for diagram number 160
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[35], w_fp[6], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 161 OF 256 ***

      // Wavefunction(s) for diagram number 161
      // (none)

      // Amplitude(s) for diagram number 161
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[36], w_fp[8], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 162 OF 256 ***

      // Wavefunction(s) for diagram number 162
      // (none)

      // Amplitude(s) for diagram number 162
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[37], w_fp[6], w_fp[8], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 163 OF 256 ***

      // Wavefunction(s) for diagram number 163
      // (none)

      // Amplitude(s) for diagram number 163
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[4], w_fp[32], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 164 OF 256 ***

      // Wavefunction(s) for diagram number 164
      // (none)

      // Amplitude(s) for diagram number 164
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[4], w_fp[27], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 165 OF 256 ***

      // Wavefunction(s) for diagram number 165
      // (none)

      // Amplitude(s) for diagram number 165
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[28], w_fp[32], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 166 OF 256 ***

      // Wavefunction(s) for diagram number 166
      // (none)

      // Amplitude(s) for diagram number 166
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[9], w_fp[6], w_fp[27], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 167 OF 256 ***

      // Wavefunction(s) for diagram number 167
      // (none)

      // Amplitude(s) for diagram number 167
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[28], w_fp[17], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 168 OF 256 ***

      // Wavefunction(s) for diagram number 168
      // (none)

      // Amplitude(s) for diagram number 168
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[34], w_fp[6], w_fp[17], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 169 OF 256 ***

      // Wavefunction(s) for diagram number 169
      // (none)

      // Amplitude(s) for diagram number 169
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[37], w_fp[4], w_fp[32], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 170 OF 256 ***

      // Wavefunction(s) for diagram number 170
      // (none)

      // Amplitude(s) for diagram number 170
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[4], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 171 OF 256 ***

      // Wavefunction(s) for diagram number 171
      // (none)

      // Amplitude(s) for diagram number 171
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[36], w_fp[32], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 172 OF 256 ***

      // Wavefunction(s) for diagram number 172
      // (none)

      // Amplitude(s) for diagram number 172
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[9], w_fp[6], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 173 OF 256 ***

      // Wavefunction(s) for diagram number 173
      // (none)

      // Amplitude(s) for diagram number 173
      FFS4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[6], w_fp[33], COUPs[ndcoup + 6], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 174 OF 256 ***

      // Wavefunction(s) for diagram number 174
      // (none)

      // Amplitude(s) for diagram number 174
      FFS4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[2], w_fp[33], COUPs[ndcoup + 6], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 175 OF 256 ***

      // Wavefunction(s) for diagram number 175
      // (none)

      // Amplitude(s) for diagram number 175
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[36], w_fp[17], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 176 OF 256 ***

      // Wavefunction(s) for diagram number 176
      // (none)

      // Amplitude(s) for diagram number 176
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[37], w_fp[6], w_fp[17], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 177 OF 256 ***

      // Wavefunction(s) for diagram number 177
      // (none)

      // Amplitude(s) for diagram number 177
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[6], w_fp[12], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 178 OF 256 ***

      // Wavefunction(s) for diagram number 178
      // (none)

      // Amplitude(s) for diagram number 178
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[6], w_fp[13], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 179 OF 256 ***

      // Wavefunction(s) for diagram number 179
      // (none)

      // Amplitude(s) for diagram number 179
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[2], w_fp[12], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 180 OF 256 ***

      // Wavefunction(s) for diagram number 180
      // (none)

      // Amplitude(s) for diagram number 180
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[2], w_fp[13], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 181 OF 256 ***

      // Wavefunction(s) for diagram number 181
      // (none)

      // Amplitude(s) for diagram number 181
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[6], w_fp[15], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 182 OF 256 ***

      // Wavefunction(s) for diagram number 182
      // (none)

      // Amplitude(s) for diagram number 182
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[6], w_fp[16], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 183 OF 256 ***

      // Wavefunction(s) for diagram number 183
      // (none)

      // Amplitude(s) for diagram number 183
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[2], w_fp[15], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 184 OF 256 ***

      // Wavefunction(s) for diagram number 184
      // (none)

      // Amplitude(s) for diagram number 184
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[2], w_fp[16], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 185 OF 256 ***

      // Wavefunction(s) for diagram number 185
      // (none)

      // Amplitude(s) for diagram number 185
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[6], w_fp[19], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 186 OF 256 ***

      // Wavefunction(s) for diagram number 186
      // (none)

      // Amplitude(s) for diagram number 186
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[6], w_fp[20], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 187 OF 256 ***

      // Wavefunction(s) for diagram number 187
      // (none)

      // Amplitude(s) for diagram number 187
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[2], w_fp[19], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 188 OF 256 ***

      // Wavefunction(s) for diagram number 188
      // (none)

      // Amplitude(s) for diagram number 188
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[2], w_fp[20], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 189 OF 256 ***

      // Wavefunction(s) for diagram number 189
      // (none)

      // Amplitude(s) for diagram number 189
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[6], w_fp[22], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 190 OF 256 ***

      // Wavefunction(s) for diagram number 190
      // (none)

      // Amplitude(s) for diagram number 190
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[6], w_fp[23], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 191 OF 256 ***

      // Wavefunction(s) for diagram number 191
      // (none)

      // Amplitude(s) for diagram number 191
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[2], w_fp[22], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 192 OF 256 ***

      // Wavefunction(s) for diagram number 192
      // (none)

      // Amplitude(s) for diagram number 192
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[2], w_fp[23], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 193 OF 256 ***

      // Wavefunction(s) for diagram number 193
      FFV1_1<W_ACCESS, CI_ACCESS>( w_fp[0], w_fp[27], COUPs[ndcoup + 0], 1.0, 0., 0., w_fp[23] );

      // Amplitude(s) for diagram number 193
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[1], w_fp[23], w_fp[25], COUPs[ndcoup + 0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 194 OF 256 ***

      // Wavefunction(s) for diagram number 194
      // (none)

      // Amplitude(s) for diagram number 194
      FFV2_4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[1], w_fp[23], w_fp[26], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 195 OF 256 ***

      // Wavefunction(s) for diagram number 195
      FFV1_2<W_ACCESS, CI_ACCESS>( w_fp[1], w_fp[27], COUPs[ndcoup + 0], 1.0, 0., 0., w_fp[22] );

      // Amplitude(s) for diagram number 195
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[22], w_fp[0], w_fp[25], COUPs[ndcoup + 0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 196 OF 256 ***

      // Wavefunction(s) for diagram number 196
      // (none)

      // Amplitude(s) for diagram number 196
      FFV2_4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[22], w_fp[0], w_fp[26], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 197 OF 256 ***

      // Wavefunction(s) for diagram number 197
      FFV2_4_1<W_ACCESS, CI_ACCESS>( w_fp[0], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, 0., 0., w_fp[27] );

      // Amplitude(s) for diagram number 197
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[1], w_fp[27], w_fp[25], COUPs[ndcoup + 0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 198 OF 256 ***

      // Wavefunction(s) for diagram number 198
      // (none)

      // Amplitude(s) for diagram number 198
      FFV2_4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[1], w_fp[27], w_fp[26], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 199 OF 256 ***

      // Wavefunction(s) for diagram number 199
      FFV2_4_2<W_ACCESS, CI_ACCESS>( w_fp[1], w_fp[30], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, 0., 0., w_fp[20] );

      // Amplitude(s) for diagram number 199
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[20], w_fp[0], w_fp[25], COUPs[ndcoup + 0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 200 OF 256 ***

      // Wavefunction(s) for diagram number 200
      // (none)

      // Amplitude(s) for diagram number 200
      FFV2_4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[20], w_fp[0], w_fp[26], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 201 OF 256 ***

      // Wavefunction(s) for diagram number 201
      FFV1P0_3<W_ACCESS, CI_ACCESS>( w_fp[1], w_fp[23], COUPs[ndcoup + 0], 1.0, 0., 0., w_fp[26] );

      // Amplitude(s) for diagram number 201
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[11], w_fp[26], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 202 OF 256 ***

      // Wavefunction(s) for diagram number 202
      FFV2_4_3<W_ACCESS, CI_ACCESS>( w_fp[1], w_fp[23], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, cIPD[1], cIPD[2], w_fp[25] );

      // Amplitude(s) for diagram number 202
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[11], w_fp[25], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 203 OF 256 ***

      // Wavefunction(s) for diagram number 203
      FFV1P0_3<W_ACCESS, CI_ACCESS>( w_fp[22], w_fp[0], COUPs[ndcoup + 0], 1.0, 0., 0., w_fp[30] );

      // Amplitude(s) for diagram number 203
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[11], w_fp[30], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 204 OF 256 ***

      // Wavefunction(s) for diagram number 204
      FFV2_4_3<W_ACCESS, CI_ACCESS>( w_fp[22], w_fp[0], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, cIPD[1], cIPD[2], w_fp[19] );

      // Amplitude(s) for diagram number 204
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[11], w_fp[19], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 205 OF 256 ***

      // Wavefunction(s) for diagram number 205
      FFV1P0_3<W_ACCESS, CI_ACCESS>( w_fp[1], w_fp[27], COUPs[ndcoup + 0], 1.0, 0., 0., w_fp[16] );

      // Amplitude(s) for diagram number 205
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[11], w_fp[16], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 206 OF 256 ***

      // Wavefunction(s) for diagram number 206
      FFV2_4_3<W_ACCESS, CI_ACCESS>( w_fp[1], w_fp[27], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, cIPD[1], cIPD[2], w_fp[15] );

      // Amplitude(s) for diagram number 206
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[11], w_fp[15], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 207 OF 256 ***

      // Wavefunction(s) for diagram number 207
      FFV1P0_3<W_ACCESS, CI_ACCESS>( w_fp[20], w_fp[0], COUPs[ndcoup + 0], 1.0, 0., 0., w_fp[13] );

      // Amplitude(s) for diagram number 207
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[11], w_fp[13], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 208 OF 256 ***

      // Wavefunction(s) for diagram number 208
      FFV2_4_3<W_ACCESS, CI_ACCESS>( w_fp[20], w_fp[0], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, cIPD[1], cIPD[2], w_fp[12] );

      // Amplitude(s) for diagram number 208
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[11], w_fp[12], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 209 OF 256 ***

      // Wavefunction(s) for diagram number 209
      // (none)

      // Amplitude(s) for diagram number 209
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[1], w_fp[23], w_fp[38], COUPs[ndcoup + 0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 210 OF 256 ***

      // Wavefunction(s) for diagram number 210
      // (none)

      // Amplitude(s) for diagram number 210
      FFV2_4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[1], w_fp[23], w_fp[39], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 211 OF 256 ***

      // Wavefunction(s) for diagram number 211
      // (none)

      // Amplitude(s) for diagram number 211
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[22], w_fp[0], w_fp[38], COUPs[ndcoup + 0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 212 OF 256 ***

      // Wavefunction(s) for diagram number 212
      // (none)

      // Amplitude(s) for diagram number 212
      FFV2_4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[22], w_fp[0], w_fp[39], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 213 OF 256 ***

      // Wavefunction(s) for diagram number 213
      // (none)

      // Amplitude(s) for diagram number 213
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[1], w_fp[27], w_fp[38], COUPs[ndcoup + 0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 214 OF 256 ***

      // Wavefunction(s) for diagram number 214
      // (none)

      // Amplitude(s) for diagram number 214
      FFV2_4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[1], w_fp[27], w_fp[39], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 215 OF 256 ***

      // Wavefunction(s) for diagram number 215
      // (none)

      // Amplitude(s) for diagram number 215
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[20], w_fp[0], w_fp[38], COUPs[ndcoup + 0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 216 OF 256 ***

      // Wavefunction(s) for diagram number 216
      // (none)

      // Amplitude(s) for diagram number 216
      FFV2_4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[20], w_fp[0], w_fp[39], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 3], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 217 OF 256 ***

      // Wavefunction(s) for diagram number 217
      // (none)

      // Amplitude(s) for diagram number 217
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[31], w_fp[26], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 218 OF 256 ***

      // Wavefunction(s) for diagram number 218
      // (none)

      // Amplitude(s) for diagram number 218
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[31], w_fp[25], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 219 OF 256 ***

      // Wavefunction(s) for diagram number 219
      // (none)

      // Amplitude(s) for diagram number 219
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[31], w_fp[30], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 220 OF 256 ***

      // Wavefunction(s) for diagram number 220
      // (none)

      // Amplitude(s) for diagram number 220
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[31], w_fp[19], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 221 OF 256 ***

      // Wavefunction(s) for diagram number 221
      // (none)

      // Amplitude(s) for diagram number 221
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[31], w_fp[16], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 222 OF 256 ***

      // Wavefunction(s) for diagram number 222
      // (none)

      // Amplitude(s) for diagram number 222
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[31], w_fp[15], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 223 OF 256 ***

      // Wavefunction(s) for diagram number 223
      // (none)

      // Amplitude(s) for diagram number 223
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[31], w_fp[13], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 224 OF 256 ***

      // Wavefunction(s) for diagram number 224
      // (none)

      // Amplitude(s) for diagram number 224
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[31], w_fp[12], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 225 OF 256 ***

      // Wavefunction(s) for diagram number 225
      // (none)

      // Amplitude(s) for diagram number 225
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[2], w_fp[26], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 226 OF 256 ***

      // Wavefunction(s) for diagram number 226
      // (none)

      // Amplitude(s) for diagram number 226
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[2], w_fp[25], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 227 OF 256 ***

      // Wavefunction(s) for diagram number 227
      // (none)

      // Amplitude(s) for diagram number 227
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[2], w_fp[30], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 228 OF 256 ***

      // Wavefunction(s) for diagram number 228
      // (none)

      // Amplitude(s) for diagram number 228
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[2], w_fp[19], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 229 OF 256 ***

      // Wavefunction(s) for diagram number 229
      // (none)

      // Amplitude(s) for diagram number 229
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[2], w_fp[16], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 230 OF 256 ***

      // Wavefunction(s) for diagram number 230
      // (none)

      // Amplitude(s) for diagram number 230
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[2], w_fp[15], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 231 OF 256 ***

      // Wavefunction(s) for diagram number 231
      // (none)

      // Amplitude(s) for diagram number 231
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[2], w_fp[13], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 232 OF 256 ***

      // Wavefunction(s) for diagram number 232
      // (none)

      // Amplitude(s) for diagram number 232
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[21], w_fp[2], w_fp[12], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 233 OF 256 ***

      // Wavefunction(s) for diagram number 233
      // (none)

      // Amplitude(s) for diagram number 233
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[2], w_fp[26], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 234 OF 256 ***

      // Wavefunction(s) for diagram number 234
      // (none)

      // Amplitude(s) for diagram number 234
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[2], w_fp[25], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 235 OF 256 ***

      // Wavefunction(s) for diagram number 235
      // (none)

      // Amplitude(s) for diagram number 235
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[2], w_fp[30], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 236 OF 256 ***

      // Wavefunction(s) for diagram number 236
      // (none)

      // Amplitude(s) for diagram number 236
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[2], w_fp[19], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 237 OF 256 ***

      // Wavefunction(s) for diagram number 237
      // (none)

      // Amplitude(s) for diagram number 237
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[2], w_fp[16], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 238 OF 256 ***

      // Wavefunction(s) for diagram number 238
      // (none)

      // Amplitude(s) for diagram number 238
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[2], w_fp[15], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 239 OF 256 ***

      // Wavefunction(s) for diagram number 239
      // (none)

      // Amplitude(s) for diagram number 239
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[2], w_fp[13], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 240 OF 256 ***

      // Wavefunction(s) for diagram number 240
      // (none)

      // Amplitude(s) for diagram number 240
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[5], w_fp[2], w_fp[12], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 241 OF 256 ***

      // Wavefunction(s) for diagram number 241
      // (none)

      // Amplitude(s) for diagram number 241
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[6], w_fp[26], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 242 OF 256 ***

      // Wavefunction(s) for diagram number 242
      // (none)

      // Amplitude(s) for diagram number 242
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[6], w_fp[25], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 243 OF 256 ***

      // Wavefunction(s) for diagram number 243
      // (none)

      // Amplitude(s) for diagram number 243
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[2], w_fp[26], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 244 OF 256 ***

      // Wavefunction(s) for diagram number 244
      // (none)

      // Amplitude(s) for diagram number 244
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[2], w_fp[25], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 245 OF 256 ***

      // Wavefunction(s) for diagram number 245
      // (none)

      // Amplitude(s) for diagram number 245
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[6], w_fp[30], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 246 OF 256 ***

      // Wavefunction(s) for diagram number 246
      // (none)

      // Amplitude(s) for diagram number 246
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[6], w_fp[19], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 247 OF 256 ***

      // Wavefunction(s) for diagram number 247
      // (none)

      // Amplitude(s) for diagram number 247
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[2], w_fp[30], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 248 OF 256 ***

      // Wavefunction(s) for diagram number 248
      // (none)

      // Amplitude(s) for diagram number 248
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[2], w_fp[19], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 249 OF 256 ***

      // Wavefunction(s) for diagram number 249
      // (none)

      // Amplitude(s) for diagram number 249
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[6], w_fp[16], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 250 OF 256 ***

      // Wavefunction(s) for diagram number 250
      // (none)

      // Amplitude(s) for diagram number 250
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[6], w_fp[15], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 251 OF 256 ***

      // Wavefunction(s) for diagram number 251
      // (none)

      // Amplitude(s) for diagram number 251
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[2], w_fp[16], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 252 OF 256 ***

      // Wavefunction(s) for diagram number 252
      // (none)

      // Amplitude(s) for diagram number 252
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[2], w_fp[15], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 253 OF 256 ***

      // Wavefunction(s) for diagram number 253
      // (none)

      // Amplitude(s) for diagram number 253
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[6], w_fp[13], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 254 OF 256 ***

      // Wavefunction(s) for diagram number 254
      // (none)

      // Amplitude(s) for diagram number 254
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[3], w_fp[6], w_fp[12], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 255 OF 256 ***

      // Wavefunction(s) for diagram number 255
      // (none)

      // Amplitude(s) for diagram number 255
      FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[2], w_fp[13], COUPs[ndcoup + 1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 256 OF 256 ***

      // Wavefunction(s) for diagram number 256
      // (none)

      // Amplitude(s) for diagram number 256
      FFV2_3_0<W_ACCESS, A_ACCESS, CI_ACCESS>( w_fp[7], w_fp[2], w_fp[12], COUPs[ndcoup + 2], 1.0, COUPs[ndcoup + 4], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

#/* clang-format off */

      // *** COLOR CHOICE BELOW ***

      // Store the leading color flows for choice of color
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
#ifndef MGONGPUCPP_GPUIMPL
      if( jamp2_sv ) // disable color choice if nullptr
      {
        for( int icol = 0; icol < ncolor; icol++ )
          jamp2_sv[ncolor * iParity + icol] += cxabs2( jamp_sv[icol] ); // may underflow #831
      }
#else
      assert( iParity == 0 ); // sanity check for J2_ACCESS
      using J2_ACCESS = DeviceAccessJamp2;
      if( colAllJamp2s ) // disable color choice if nullptr
      {
        for( int icol = 0; icol < ncolor; icol++ )
          // NB: atomicAdd is needed after moving to cuda streams with one helicity per stream!
          atomicAdd( &J2_ACCESS::kernelAccessIcol( colAllJamp2s, icol ), cxabs2( jamp_sv[icol] ) );
      }
#endif /* clang-format on */
#endif

      // *** PREPARE OUTPUT JAMPS ***
#ifdef MGONGPUCPP_GPUIMPL
      //printf( "calculate_jamps: dcNGoodHel=%d\n", dcNGoodHel );
      // In CUDA, copy the local jamp to the output global-memory jamp
      constexpr int ihel0 = 0; // the allJamps buffer already points to a specific helicity _within a super-buffer for dcNGoodHel helicities_
      using J_ACCESS = DeviceAccessJamp;
      for( int icol = 0; icol < ncolor; icol++ )
        J_ACCESS::kernelAccessIcolIhelNhel( allJamps, icol, ihel0, dcNGoodHel ) = jamp_sv[icol];
#else
      // In C++, copy the local jamp to the output array passed as function argument
      for( int icol = 0; icol < ncolor; icol++ )
        allJamp_sv[iParity * ncolor + icol] = jamp_sv[icol];
#endif
    }
    // END LOOP ON IPARITY

    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  CPPProcess::CPPProcess( bool verbose,
                          bool debug )
    : m_verbose( verbose )
    , m_debug( debug )
#ifndef MGONGPU_HARDCODE_PARAM
    , m_pars( 0 )
#endif
    , m_masses()
  {
    // Helicities for the process [NB do keep 'static' for this constexpr array, see issue #283]
    // *** NB There is no automatic check yet that these are in the same order as Fortran! #569 ***
    static constexpr short tHel[ncomb][npar] = {
      { -1, 1, -1, 1, -1, 1, -1, -1 },
      { -1, 1, -1, 1, -1, 1, -1, 1 },
      { -1, 1, -1, 1, -1, 1, 1, -1 },
      { -1, 1, -1, 1, -1, 1, 1, 1 },
      { -1, 1, -1, 1, -1, -1, -1, -1 },
      { -1, 1, -1, 1, -1, -1, -1, 1 },
      { -1, 1, -1, 1, -1, -1, 1, -1 },
      { -1, 1, -1, 1, -1, -1, 1, 1 },
      { -1, 1, -1, 1, 1, 1, -1, -1 },
      { -1, 1, -1, 1, 1, 1, -1, 1 },
      { -1, 1, -1, 1, 1, 1, 1, -1 },
      { -1, 1, -1, 1, 1, 1, 1, 1 },
      { -1, 1, -1, 1, 1, -1, -1, -1 },
      { -1, 1, -1, 1, 1, -1, -1, 1 },
      { -1, 1, -1, 1, 1, -1, 1, -1 },
      { -1, 1, -1, 1, 1, -1, 1, 1 },
      { -1, 1, -1, -1, -1, 1, -1, -1 },
      { -1, 1, -1, -1, -1, 1, -1, 1 },
      { -1, 1, -1, -1, -1, 1, 1, -1 },
      { -1, 1, -1, -1, -1, 1, 1, 1 },
      { -1, 1, -1, -1, -1, -1, -1, -1 },
      { -1, 1, -1, -1, -1, -1, -1, 1 },
      { -1, 1, -1, -1, -1, -1, 1, -1 },
      { -1, 1, -1, -1, -1, -1, 1, 1 },
      { -1, 1, -1, -1, 1, 1, -1, -1 },
      { -1, 1, -1, -1, 1, 1, -1, 1 },
      { -1, 1, -1, -1, 1, 1, 1, -1 },
      { -1, 1, -1, -1, 1, 1, 1, 1 },
      { -1, 1, -1, -1, 1, -1, -1, -1 },
      { -1, 1, -1, -1, 1, -1, -1, 1 },
      { -1, 1, -1, -1, 1, -1, 1, -1 },
      { -1, 1, -1, -1, 1, -1, 1, 1 },
      { -1, 1, 1, 1, -1, 1, -1, -1 },
      { -1, 1, 1, 1, -1, 1, -1, 1 },
      { -1, 1, 1, 1, -1, 1, 1, -1 },
      { -1, 1, 1, 1, -1, 1, 1, 1 },
      { -1, 1, 1, 1, -1, -1, -1, -1 },
      { -1, 1, 1, 1, -1, -1, -1, 1 },
      { -1, 1, 1, 1, -1, -1, 1, -1 },
      { -1, 1, 1, 1, -1, -1, 1, 1 },
      { -1, 1, 1, 1, 1, 1, -1, -1 },
      { -1, 1, 1, 1, 1, 1, -1, 1 },
      { -1, 1, 1, 1, 1, 1, 1, -1 },
      { -1, 1, 1, 1, 1, 1, 1, 1 },
      { -1, 1, 1, 1, 1, -1, -1, -1 },
      { -1, 1, 1, 1, 1, -1, -1, 1 },
      { -1, 1, 1, 1, 1, -1, 1, -1 },
      { -1, 1, 1, 1, 1, -1, 1, 1 },
      { -1, 1, 1, -1, -1, 1, -1, -1 },
      { -1, 1, 1, -1, -1, 1, -1, 1 },
      { -1, 1, 1, -1, -1, 1, 1, -1 },
      { -1, 1, 1, -1, -1, 1, 1, 1 },
      { -1, 1, 1, -1, -1, -1, -1, -1 },
      { -1, 1, 1, -1, -1, -1, -1, 1 },
      { -1, 1, 1, -1, -1, -1, 1, -1 },
      { -1, 1, 1, -1, -1, -1, 1, 1 },
      { -1, 1, 1, -1, 1, 1, -1, -1 },
      { -1, 1, 1, -1, 1, 1, -1, 1 },
      { -1, 1, 1, -1, 1, 1, 1, -1 },
      { -1, 1, 1, -1, 1, 1, 1, 1 },
      { -1, 1, 1, -1, 1, -1, -1, -1 },
      { -1, 1, 1, -1, 1, -1, -1, 1 },
      { -1, 1, 1, -1, 1, -1, 1, -1 },
      { -1, 1, 1, -1, 1, -1, 1, 1 },
      { -1, -1, -1, 1, -1, 1, -1, -1 },
      { -1, -1, -1, 1, -1, 1, -1, 1 },
      { -1, -1, -1, 1, -1, 1, 1, -1 },
      { -1, -1, -1, 1, -1, 1, 1, 1 },
      { -1, -1, -1, 1, -1, -1, -1, -1 },
      { -1, -1, -1, 1, -1, -1, -1, 1 },
      { -1, -1, -1, 1, -1, -1, 1, -1 },
      { -1, -1, -1, 1, -1, -1, 1, 1 },
      { -1, -1, -1, 1, 1, 1, -1, -1 },
      { -1, -1, -1, 1, 1, 1, -1, 1 },
      { -1, -1, -1, 1, 1, 1, 1, -1 },
      { -1, -1, -1, 1, 1, 1, 1, 1 },
      { -1, -1, -1, 1, 1, -1, -1, -1 },
      { -1, -1, -1, 1, 1, -1, -1, 1 },
      { -1, -1, -1, 1, 1, -1, 1, -1 },
      { -1, -1, -1, 1, 1, -1, 1, 1 },
      { -1, -1, -1, -1, -1, 1, -1, -1 },
      { -1, -1, -1, -1, -1, 1, -1, 1 },
      { -1, -1, -1, -1, -1, 1, 1, -1 },
      { -1, -1, -1, -1, -1, 1, 1, 1 },
      { -1, -1, -1, -1, -1, -1, -1, -1 },
      { -1, -1, -1, -1, -1, -1, -1, 1 },
      { -1, -1, -1, -1, -1, -1, 1, -1 },
      { -1, -1, -1, -1, -1, -1, 1, 1 },
      { -1, -1, -1, -1, 1, 1, -1, -1 },
      { -1, -1, -1, -1, 1, 1, -1, 1 },
      { -1, -1, -1, -1, 1, 1, 1, -1 },
      { -1, -1, -1, -1, 1, 1, 1, 1 },
      { -1, -1, -1, -1, 1, -1, -1, -1 },
      { -1, -1, -1, -1, 1, -1, -1, 1 },
      { -1, -1, -1, -1, 1, -1, 1, -1 },
      { -1, -1, -1, -1, 1, -1, 1, 1 },
      { -1, -1, 1, 1, -1, 1, -1, -1 },
      { -1, -1, 1, 1, -1, 1, -1, 1 },
      { -1, -1, 1, 1, -1, 1, 1, -1 },
      { -1, -1, 1, 1, -1, 1, 1, 1 },
      { -1, -1, 1, 1, -1, -1, -1, -1 },
      { -1, -1, 1, 1, -1, -1, -1, 1 },
      { -1, -1, 1, 1, -1, -1, 1, -1 },
      { -1, -1, 1, 1, -1, -1, 1, 1 },
      { -1, -1, 1, 1, 1, 1, -1, -1 },
      { -1, -1, 1, 1, 1, 1, -1, 1 },
      { -1, -1, 1, 1, 1, 1, 1, -1 },
      { -1, -1, 1, 1, 1, 1, 1, 1 },
      { -1, -1, 1, 1, 1, -1, -1, -1 },
      { -1, -1, 1, 1, 1, -1, -1, 1 },
      { -1, -1, 1, 1, 1, -1, 1, -1 },
      { -1, -1, 1, 1, 1, -1, 1, 1 },
      { -1, -1, 1, -1, -1, 1, -1, -1 },
      { -1, -1, 1, -1, -1, 1, -1, 1 },
      { -1, -1, 1, -1, -1, 1, 1, -1 },
      { -1, -1, 1, -1, -1, 1, 1, 1 },
      { -1, -1, 1, -1, -1, -1, -1, -1 },
      { -1, -1, 1, -1, -1, -1, -1, 1 },
      { -1, -1, 1, -1, -1, -1, 1, -1 },
      { -1, -1, 1, -1, -1, -1, 1, 1 },
      { -1, -1, 1, -1, 1, 1, -1, -1 },
      { -1, -1, 1, -1, 1, 1, -1, 1 },
      { -1, -1, 1, -1, 1, 1, 1, -1 },
      { -1, -1, 1, -1, 1, 1, 1, 1 },
      { -1, -1, 1, -1, 1, -1, -1, -1 },
      { -1, -1, 1, -1, 1, -1, -1, 1 },
      { -1, -1, 1, -1, 1, -1, 1, -1 },
      { -1, -1, 1, -1, 1, -1, 1, 1 },
      { 1, 1, -1, 1, -1, 1, -1, -1 },
      { 1, 1, -1, 1, -1, 1, -1, 1 },
      { 1, 1, -1, 1, -1, 1, 1, -1 },
      { 1, 1, -1, 1, -1, 1, 1, 1 },
      { 1, 1, -1, 1, -1, -1, -1, -1 },
      { 1, 1, -1, 1, -1, -1, -1, 1 },
      { 1, 1, -1, 1, -1, -1, 1, -1 },
      { 1, 1, -1, 1, -1, -1, 1, 1 },
      { 1, 1, -1, 1, 1, 1, -1, -1 },
      { 1, 1, -1, 1, 1, 1, -1, 1 },
      { 1, 1, -1, 1, 1, 1, 1, -1 },
      { 1, 1, -1, 1, 1, 1, 1, 1 },
      { 1, 1, -1, 1, 1, -1, -1, -1 },
      { 1, 1, -1, 1, 1, -1, -1, 1 },
      { 1, 1, -1, 1, 1, -1, 1, -1 },
      { 1, 1, -1, 1, 1, -1, 1, 1 },
      { 1, 1, -1, -1, -1, 1, -1, -1 },
      { 1, 1, -1, -1, -1, 1, -1, 1 },
      { 1, 1, -1, -1, -1, 1, 1, -1 },
      { 1, 1, -1, -1, -1, 1, 1, 1 },
      { 1, 1, -1, -1, -1, -1, -1, -1 },
      { 1, 1, -1, -1, -1, -1, -1, 1 },
      { 1, 1, -1, -1, -1, -1, 1, -1 },
      { 1, 1, -1, -1, -1, -1, 1, 1 },
      { 1, 1, -1, -1, 1, 1, -1, -1 },
      { 1, 1, -1, -1, 1, 1, -1, 1 },
      { 1, 1, -1, -1, 1, 1, 1, -1 },
      { 1, 1, -1, -1, 1, 1, 1, 1 },
      { 1, 1, -1, -1, 1, -1, -1, -1 },
      { 1, 1, -1, -1, 1, -1, -1, 1 },
      { 1, 1, -1, -1, 1, -1, 1, -1 },
      { 1, 1, -1, -1, 1, -1, 1, 1 },
      { 1, 1, 1, 1, -1, 1, -1, -1 },
      { 1, 1, 1, 1, -1, 1, -1, 1 },
      { 1, 1, 1, 1, -1, 1, 1, -1 },
      { 1, 1, 1, 1, -1, 1, 1, 1 },
      { 1, 1, 1, 1, -1, -1, -1, -1 },
      { 1, 1, 1, 1, -1, -1, -1, 1 },
      { 1, 1, 1, 1, -1, -1, 1, -1 },
      { 1, 1, 1, 1, -1, -1, 1, 1 },
      { 1, 1, 1, 1, 1, 1, -1, -1 },
      { 1, 1, 1, 1, 1, 1, -1, 1 },
      { 1, 1, 1, 1, 1, 1, 1, -1 },
      { 1, 1, 1, 1, 1, 1, 1, 1 },
      { 1, 1, 1, 1, 1, -1, -1, -1 },
      { 1, 1, 1, 1, 1, -1, -1, 1 },
      { 1, 1, 1, 1, 1, -1, 1, -1 },
      { 1, 1, 1, 1, 1, -1, 1, 1 },
      { 1, 1, 1, -1, -1, 1, -1, -1 },
      { 1, 1, 1, -1, -1, 1, -1, 1 },
      { 1, 1, 1, -1, -1, 1, 1, -1 },
      { 1, 1, 1, -1, -1, 1, 1, 1 },
      { 1, 1, 1, -1, -1, -1, -1, -1 },
      { 1, 1, 1, -1, -1, -1, -1, 1 },
      { 1, 1, 1, -1, -1, -1, 1, -1 },
      { 1, 1, 1, -1, -1, -1, 1, 1 },
      { 1, 1, 1, -1, 1, 1, -1, -1 },
      { 1, 1, 1, -1, 1, 1, -1, 1 },
      { 1, 1, 1, -1, 1, 1, 1, -1 },
      { 1, 1, 1, -1, 1, 1, 1, 1 },
      { 1, 1, 1, -1, 1, -1, -1, -1 },
      { 1, 1, 1, -1, 1, -1, -1, 1 },
      { 1, 1, 1, -1, 1, -1, 1, -1 },
      { 1, 1, 1, -1, 1, -1, 1, 1 },
      { 1, -1, -1, 1, -1, 1, -1, -1 },
      { 1, -1, -1, 1, -1, 1, -1, 1 },
      { 1, -1, -1, 1, -1, 1, 1, -1 },
      { 1, -1, -1, 1, -1, 1, 1, 1 },
      { 1, -1, -1, 1, -1, -1, -1, -1 },
      { 1, -1, -1, 1, -1, -1, -1, 1 },
      { 1, -1, -1, 1, -1, -1, 1, -1 },
      { 1, -1, -1, 1, -1, -1, 1, 1 },
      { 1, -1, -1, 1, 1, 1, -1, -1 },
      { 1, -1, -1, 1, 1, 1, -1, 1 },
      { 1, -1, -1, 1, 1, 1, 1, -1 },
      { 1, -1, -1, 1, 1, 1, 1, 1 },
      { 1, -1, -1, 1, 1, -1, -1, -1 },
      { 1, -1, -1, 1, 1, -1, -1, 1 },
      { 1, -1, -1, 1, 1, -1, 1, -1 },
      { 1, -1, -1, 1, 1, -1, 1, 1 },
      { 1, -1, -1, -1, -1, 1, -1, -1 },
      { 1, -1, -1, -1, -1, 1, -1, 1 },
      { 1, -1, -1, -1, -1, 1, 1, -1 },
      { 1, -1, -1, -1, -1, 1, 1, 1 },
      { 1, -1, -1, -1, -1, -1, -1, -1 },
      { 1, -1, -1, -1, -1, -1, -1, 1 },
      { 1, -1, -1, -1, -1, -1, 1, -1 },
      { 1, -1, -1, -1, -1, -1, 1, 1 },
      { 1, -1, -1, -1, 1, 1, -1, -1 },
      { 1, -1, -1, -1, 1, 1, -1, 1 },
      { 1, -1, -1, -1, 1, 1, 1, -1 },
      { 1, -1, -1, -1, 1, 1, 1, 1 },
      { 1, -1, -1, -1, 1, -1, -1, -1 },
      { 1, -1, -1, -1, 1, -1, -1, 1 },
      { 1, -1, -1, -1, 1, -1, 1, -1 },
      { 1, -1, -1, -1, 1, -1, 1, 1 },
      { 1, -1, 1, 1, -1, 1, -1, -1 },
      { 1, -1, 1, 1, -1, 1, -1, 1 },
      { 1, -1, 1, 1, -1, 1, 1, -1 },
      { 1, -1, 1, 1, -1, 1, 1, 1 },
      { 1, -1, 1, 1, -1, -1, -1, -1 },
      { 1, -1, 1, 1, -1, -1, -1, 1 },
      { 1, -1, 1, 1, -1, -1, 1, -1 },
      { 1, -1, 1, 1, -1, -1, 1, 1 },
      { 1, -1, 1, 1, 1, 1, -1, -1 },
      { 1, -1, 1, 1, 1, 1, -1, 1 },
      { 1, -1, 1, 1, 1, 1, 1, -1 },
      { 1, -1, 1, 1, 1, 1, 1, 1 },
      { 1, -1, 1, 1, 1, -1, -1, -1 },
      { 1, -1, 1, 1, 1, -1, -1, 1 },
      { 1, -1, 1, 1, 1, -1, 1, -1 },
      { 1, -1, 1, 1, 1, -1, 1, 1 },
      { 1, -1, 1, -1, -1, 1, -1, -1 },
      { 1, -1, 1, -1, -1, 1, -1, 1 },
      { 1, -1, 1, -1, -1, 1, 1, -1 },
      { 1, -1, 1, -1, -1, 1, 1, 1 },
      { 1, -1, 1, -1, -1, -1, -1, -1 },
      { 1, -1, 1, -1, -1, -1, -1, 1 },
      { 1, -1, 1, -1, -1, -1, 1, -1 },
      { 1, -1, 1, -1, -1, -1, 1, 1 },
      { 1, -1, 1, -1, 1, 1, -1, -1 },
      { 1, -1, 1, -1, 1, 1, -1, 1 },
      { 1, -1, 1, -1, 1, 1, 1, -1 },
      { 1, -1, 1, -1, 1, 1, 1, 1 },
      { 1, -1, 1, -1, 1, -1, -1, -1 },
      { 1, -1, 1, -1, 1, -1, -1, 1 },
      { 1, -1, 1, -1, 1, -1, 1, -1 },
      { 1, -1, 1, -1, 1, -1, 1, 1 } };
#ifdef MGONGPUCPP_GPUIMPL
    gpuMemcpyToSymbol( cHel, tHel, ncomb * npar * sizeof( short ) );
#else
    memcpy( cHel, tHel, ncomb * npar * sizeof( short ) );
#endif

    // Enable SIGFPE traps for Floating Point Exceptions
#ifdef MGONGPUCPP_DEBUG
    fpeEnable();
#endif
  }

  //--------------------------------------------------------------------------

  CPPProcess::~CPPProcess() {}

  //--------------------------------------------------------------------------

#ifndef MGONGPU_HARDCODE_PARAM
  // Initialize process (with parameters read from user cards)
  void
  CPPProcess::initProc( const std::string& param_card_name )
  {
    // Instantiate the model class and set parameters that stay fixed during run
    m_pars = Parameters_sm::getInstance();
    SLHAReader slha( param_card_name, m_verbose );
    m_pars->setIndependentParameters( slha );
    m_pars->setIndependentCouplings();
    //m_pars->setDependentParameters(); // now computed event-by-event (running alphas #373)
    //m_pars->setDependentCouplings(); // now computed event-by-event (running alphas #373)
    if( m_verbose )
    {
      m_pars->printIndependentParameters();
      m_pars->printIndependentCouplings();
      //m_pars->printDependentParameters(); // now computed event-by-event (running alphas #373)
      //m_pars->printDependentCouplings(); // now computed event-by-event (running alphas #373)
    }
    // Set external particle masses for this matrix element
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->mdl_MB );
    m_masses.push_back( m_pars->mdl_MB );
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->ZERO );
#ifdef MGONGPUCPP_GPUIMPL
    // Create the normalized color matrix in device memory
    createNormalizedColorMatrix();
#endif
    // Read physics parameters like masses and couplings from user configuration files (static: initialize once)
    // Then copy them to CUDA constant memory (issue #39) or its C++ emulation in file-scope static memory
    const fptype tIPD[nIPD] = { (fptype)m_pars->mdl_MB, (fptype)m_pars->mdl_MZ, (fptype)m_pars->mdl_WZ, (fptype)m_pars->mdl_MH, (fptype)m_pars->mdl_WH };
    const cxtype tIPC[nIPC] = { cxmake( m_pars->GC_3 ), cxmake( m_pars->GC_1 ), cxmake( m_pars->GC_50 ), cxmake( m_pars->GC_59 ), cxmake( m_pars->GC_58 ), cxmake( m_pars->GC_81 ), cxmake( m_pars->GC_83 ) };
#ifdef MGONGPUCPP_GPUIMPL
    gpuMemcpyToSymbol( cIPD, tIPD, nIPD * sizeof( fptype ) );
    gpuMemcpyToSymbol( cIPC, tIPC, nIPC * sizeof( cxtype ) );
#ifdef MGONGPUCPP_NBSMINDEPPARAM_GT_0
    if( Parameters_sm::nBsmIndepParam > 0 )
      gpuMemcpyToSymbol( bsmIndepParam, m_pars->mdl_bsmIndepParam, Parameters_sm::nBsmIndepParam * sizeof( double ) );
#endif
#else
    memcpy( cIPD, tIPD, nIPD * sizeof( fptype ) );
    memcpy( cIPC, tIPC, nIPC * sizeof( cxtype ) );
#ifdef MGONGPUCPP_NBSMINDEPPARAM_GT_0
    if( Parameters_sm::nBsmIndepParam > 0 )
      memcpy( bsmIndepParam, m_pars->mdl_bsmIndepParam, Parameters_sm::nBsmIndepParam * sizeof( double ) );
#endif
#endif
    //for ( int i=0; i<nIPD; i++ ) std::cout << std::setprecision(17) << "tIPD[i] = " << tIPD[i] << std::endl;
    //for ( int i=0; i<nIPC; i++ ) std::cout << std::setprecision(17) << "tIPC[i] = " << tIPC[i] << std::endl;
    //for ( int i=0; i<Parameters_sm::nBsmIndepParam; i++ ) std::cout << std::setprecision(17) << "m_pars->mdl_bsmIndepParam[i] = " << m_pars->mdl_bsmIndepParam[i] << std::endl;
  }
#else
  // Initialize process (with hardcoded parameters)
  void
  CPPProcess::initProc( const std::string& /*param_card_name*/ )
  {
    // Use hardcoded physics parameters
    if( m_verbose )
    {
      Parameters_sm::printIndependentParameters();
      Parameters_sm::printIndependentCouplings();
      //Parameters_sm::printDependentParameters(); // now computed event-by-event (running alphas #373)
      //Parameters_sm::printDependentCouplings(); // now computed event-by-event (running alphas #373)
    }
    // Set external particle masses for this matrix element
    m_masses.push_back( Parameters_sm::ZERO );
    m_masses.push_back( Parameters_sm::ZERO );
    m_masses.push_back( Parameters_sm::mdl_MB );
    m_masses.push_back( Parameters_sm::mdl_MB );
    m_masses.push_back( Parameters_sm::ZERO );
    m_masses.push_back( Parameters_sm::ZERO );
    m_masses.push_back( Parameters_sm::ZERO );
    m_masses.push_back( Parameters_sm::ZERO );
#ifdef MGONGPUCPP_GPUIMPL
    // Create the normalized color matrix in device memory
    createNormalizedColorMatrix();
#endif
  }
#endif

  //--------------------------------------------------------------------------

  // Retrieve the compiler that was used to build this module
  const std::string
  CPPProcess::getCompiler()
  {
    std::stringstream out;
    // HIP version (HIPCC)
    // [Use __HIPCC__ instead of MGONGPUCPP_GPUIMPL here!]
    // [This tests if 'hipcc' was used even to build a .cc file, even if not necessarily 'nvcc -x cu' for a .cu file]
    // [Check 'hipcc -dM -E -x hip -I ../../src CPPProcess.cc | grep HIP']
#ifdef __HIPCC__
#if defined HIP_VERSION_MAJOR && defined HIP_VERSION_MINOR && defined HIP_VERSION_PATCH
    out << "hipcc " << HIP_VERSION_MAJOR << "." << HIP_VERSION_MINOR << "." << HIP_VERSION_PATCH;
#else
    out << "hipcc UNKNOWN";
#endif
    out << " (";
#endif
    // CUDA version (NVCC)
    // [Use __NVCC__ instead of MGONGPUCPP_GPUIMPL here!]
    // [This tests if 'nvcc' was used even to build a .cc file, even if not necessarily 'nvcc -x cu' for a .cu file]
    // [Check 'nvcc --compiler-options -dM -E dummy.c | grep CUDA': see https://stackoverflow.com/a/53713712]
#ifdef __NVCC__
#if defined __CUDACC_VER_MAJOR__ && defined __CUDACC_VER_MINOR__ && defined __CUDACC_VER_BUILD__
    out << "nvcc " << __CUDACC_VER_MAJOR__ << "." << __CUDACC_VER_MINOR__ << "." << __CUDACC_VER_BUILD__;
#else
    out << "nvcc UNKNOWN";
#endif
    out << " (";
#endif
    // ICX version (either as CXX or as host compiler inside NVCC)
#if defined __INTEL_COMPILER
#error "icc is no longer supported: please use icx"
#elif defined __INTEL_LLVM_COMPILER // alternative: __INTEL_CLANG_COMPILER
    out << "icx " << __INTEL_LLVM_COMPILER;
#ifdef __NVCC__
    out << ", ";
#else
    out << " (";
#endif
#endif
    // CLANG version (either as CXX or as host compiler inside NVCC or inside ICX)
#if defined __clang__
#if defined __clang_major__ && defined __clang_minor__ && defined __clang_patchlevel__
#ifdef __APPLE__
    out << "Apple clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
#else
    out << "clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
    /*
    // === AV 26-Jan-2024 DISABLE THIS CODE (START)
    // === AV 26-Jan-2024 First, it is totally wrong to assume that the CXX environment variable is used in the build!
    // === AV 26-Jan-2024 Second and worse, here we need build time values, while CXX in this code is evaluated at runtime!
    // GCC toolchain version inside CLANG
    std::string tchainout;
    std::string tchaincmd = "readelf -p .comment $(${CXX} -print-libgcc-file-name) |& grep 'GCC: (GNU)' | grep -v Warning | sort -u | awk '{print $5}'";
    std::unique_ptr<FILE, decltype( &pclose )> tchainpipe( popen( tchaincmd.c_str(), "r" ), pclose );
    if( !tchainpipe ) throw std::runtime_error( "`readelf ...` failed?" );
    std::array<char, 128> tchainbuf;
    while( fgets( tchainbuf.data(), tchainbuf.size(), tchainpipe.get() ) != nullptr ) tchainout += tchainbuf.data();
    tchainout.pop_back(); // remove trailing newline
#if defined __NVCC__ or defined __INTEL_LLVM_COMPILER
    out << ", gcc " << tchainout;
#else
    out << " (gcc " << tchainout << ")";
#endif
    // === AV 26-Jan-2024 DISABLE THIS CODE (END)
    */
#endif
#else
    out << "clang UNKNOWKN";
#endif
#else
    // GCC version (either as CXX or as host compiler inside NVCC)
#if defined __GNUC__ && defined __GNUC_MINOR__ && defined __GNUC_PATCHLEVEL__
    out << "gcc " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__;
#else
    out << "gcc UNKNOWKN";
#endif
#endif
#if defined __HIPCC__ or defined __NVCC__ or defined __INTEL_LLVM_COMPILER
    out << ")";
#endif
    return out.str();
  }

  //--------------------------------------------------------------------------

  __global__ void /* clang-format off */
  computeDependentCouplings( const fptype* allgs, // input: Gs[nevt]
                             fptype* allcouplings // output: couplings[nevt*ndcoup*2]
#ifndef MGONGPUCPP_GPUIMPL
                             , const int nevt     // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
  ) /* clang-format on */
  {
#ifdef MGONGPUCPP_GPUIMPL
    using namespace mg5amcGpu;
    using G_ACCESS = DeviceAccessGs;
    using C_ACCESS = DeviceAccessCouplings;
    G2COUP<G_ACCESS, C_ACCESS>( allgs, allcouplings, bsmIndepParam );
#else
    using namespace mg5amcCpu;
    using G_ACCESS = HostAccessGs;
    using C_ACCESS = HostAccessCouplings;
    for( int ipagV = 0; ipagV < nevt / neppV; ++ipagV )
    {
      const int ievt0 = ipagV * neppV;
      const fptype* gs = MemoryAccessGs::ieventAccessRecordConst( allgs, ievt0 );
      fptype* couplings = MemoryAccessCouplings::ieventAccessRecord( allcouplings, ievt0 );
      G2COUP<G_ACCESS, C_ACCESS>( gs, couplings, bsmIndepParam );
    }
#endif
  }

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  void /* clang-format off */
  sigmaKin_getGoodHel( const fptype* allmomenta,   // input: momenta[nevt*npar*4]
                       const fptype* allcouplings, // input: couplings[nevt*ndcoup*2]
                       fptype* allMEs,             // output: allMEs[nevt], |M|^2 final_avg_over_helicities
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                       fptype* allNumerators,      // output: multichannel numerators[nevt], running_sum_over_helicities
                       fptype* allDenominators,    // output: multichannel denominators[nevt], running_sum_over_helicities
#endif
                       fptype_sv* allJamps,        // tmp: jamp[ncolor*2*nevt] _for one helicity_ (reused in the getGoodHel helicity loop)
                       bool* isGoodHel,            // output: isGoodHel[ncomb] - host array
                       const int nevt )            // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
  { /* clang-format on */
    const int maxtry0 = 16;
    fptype hstMEs[maxtry0];
    const int maxtry = std::min( maxtry0, nevt ); // 16, but at most nevt (avoid invalid memory access if nevt<maxtry0)
    //std::cout << "sigmaKin_getGoodHel nevt=" << nevt << " maxtry=" << maxtry << std::endl;
    for( int ihel = 0; ihel < ncomb; ihel++ )
    {
      const int gpublocks = 1;
      const int gputhreads = maxtry;
      constexpr int nOneHel = 1; // use a jamp buffer for a single helicity
      gpuMemcpyToSymbol( dcNGoodHel, &nOneHel, sizeof( int ) );
      // NEW IMPLEMENTATION OF GETGOODHEL (#630): RESET THE RUNNING SUM OVER HELICITIES TO 0 BEFORE ADDING A NEW HELICITY
      gpuMemset( allMEs, 0, maxtry * sizeof( fptype ) );
      // NB: color_sum ADDS |M|^2 for one helicity to the running sum of |M|^2 over helicities for the given event(s)
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      constexpr fptype_sv* allJamp2s = nullptr;        // no need for color selection during helicity filtering
      constexpr unsigned int* allChannelIds = nullptr; // disable multichannel single-diagram enhancement
      gpuLaunchKernel( calculate_jamps, gpublocks, gputhreads, ihel, allmomenta, allcouplings, allJamps, allChannelIds, allNumerators, allDenominators, allJamp2s, gpublocks * gputhreads );
#else
      gpuLaunchKernel( calculate_jamps, gpublocks, gputhreads, ihel, allmomenta, allcouplings, allJamps, gpublocks * gputhreads );
#endif
      gpuLaunchKernel( color_sum_kernel, gpublocks, gputhreads, allMEs, allJamps, nOneHel );
      gpuMemcpy( hstMEs, allMEs, maxtry * sizeof( fptype ), gpuMemcpyDeviceToHost );
      //std::cout << "sigmaKin_getGoodHel ihel=" << ihel << std::endl;
      for( int ievt = 0; ievt < maxtry; ++ievt )
      {
        //std::cout << "sigmaKin_getGoodHel hstMEs[ievt]=" << hstMEs[ievt] << std::endl;
        if( hstMEs[ievt] != 0 ) // NEW IMPLEMENTATION OF GETGOODHEL (#630): COMPARE EACH HELICITY CONTRIBUTION TO 0
        {
          //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
          isGoodHel[ihel] = true;
        }
      }
    }
  }
#else
  void
  sigmaKin_getGoodHel( const fptype* allmomenta,   // input: momenta[nevt*npar*4]
                       const fptype* allcouplings, // input: couplings[nevt*ndcoup*2]
                       fptype* allMEs,             // output: allMEs[nevt], |M|^2 final_avg_over_helicities
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                       fptype* allNumerators,      // output: multichannel numerators[nevt], running_sum_over_helicities
                       fptype* allDenominators,    // output: multichannel denominators[nevt], running_sum_over_helicities
#endif
                       bool* isGoodHel,            // output: isGoodHel[ncomb] - host array
                       const int nevt )            // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
  {
    //assert( (size_t)(allmomenta) % mgOnGpu::cppAlign == 0 ); // SANITY CHECK: require SIMD-friendly alignment [COMMENT OUT TO TEST MISALIGNED ACCESS]
    //assert( (size_t)(allMEs) % mgOnGpu::cppAlign == 0 ); // SANITY CHECK: require SIMD-friendly alignment [COMMENT OUT TO TEST MISALIGNED ACCESS]
    // Allocate arrays at build time to contain at least 16 events (or at least neppV events if neppV>16, e.g. in future VPUs)
    constexpr int maxtry0 = std::max( 16, neppV ); // 16, but at least neppV (otherwise the npagV loop does not even start)
    // Loop over only nevt events if nevt is < 16 (note that nevt is always >= neppV)
    assert( nevt >= neppV );
    const int maxtry = std::min( maxtry0, nevt ); // 16, but at most nevt (avoid invalid memory access if nevt<maxtry0)
    //std::cout << "sigmaKin_getGoodHel nevt=" << nevt << " maxtry=" << maxtry << std::endl;
    // HELICITY LOOP: CALCULATE WAVEFUNCTIONS
    const int npagV = maxtry / neppV;
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT /* clang-format off */
    // Mixed fptypes #537: float for color algebra and double elsewhere
    // Delay color algebra and ME updates (only on even pages)
    assert( npagV % 2 == 0 ); // SANITY CHECK for mixed fptypes: two neppV-pages are merged to one 2*neppV-page
    const int npagV2 = npagV / 2; // loop on two SIMD pages (neppV events) at a time
#else
    const int npagV2 = npagV; // loop on one SIMD page (neppV events) at a time
#endif /* clang-format on */
    for( int ipagV2 = 0; ipagV2 < npagV2; ++ipagV2 )
    {
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT /* clang-format off */
      const int ievt00 = ipagV2 * neppV * 2; // loop on two SIMD pages (neppV events) at a time
#else
      const int ievt00 = ipagV2 * neppV; // loop on one SIMD page (neppV events) at a time
#endif /* clang-format on */
      for( int ihel = 0; ihel < ncomb; ihel++ )
      {
        //std::cout << "sigmaKin_getGoodHel ihel=" << ihel << std::endl;
        // NEW IMPLEMENTATION OF GETGOODHEL (#630): RESET THE RUNNING SUM OVER HELICITIES TO 0 BEFORE ADDING A NEW HELICITY
        for( int ieppV = 0; ieppV < neppV; ++ieppV )
        {
          const int ievt = ievt00 + ieppV;
          allMEs[ievt] = 0;
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
          const int ievt2 = ievt00 + ieppV + neppV;
          allMEs[ievt2] = 0;
#endif
        }
        constexpr fptype_sv* jamp2_sv = nullptr; // no need for color selection during helicity filtering
        //std::cout << "sigmaKin_getGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
        cxtype_sv jamp_sv[2 * ncolor] = {}; // all zeros
#else
        cxtype_sv jamp_sv[ncolor] = {};  // all zeros
#endif
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL /* clang-format off */
        constexpr unsigned int channelId = 0; // disable multichannel single-diagram enhancement
        calculate_jamps( ihel, allmomenta, allcouplings, jamp_sv, channelId, allNumerators, allDenominators, jamp2_sv, ievt00 ); //maxtry?
#else
        calculate_jamps( ihel, allmomenta, allcouplings, jamp_sv, ievt00 ); //maxtry?
#endif /* clang-format on */
        color_sum_cpu( allMEs, jamp_sv, ievt00 );
        for( int ieppV = 0; ieppV < neppV; ++ieppV )
        {
          const int ievt = ievt00 + ieppV;
          //std::cout << "sigmaKin_getGoodHel allMEs[ievt]=" << allMEs[ievt] << std::endl;
          if( allMEs[ievt] != 0 ) // NEW IMPLEMENTATION OF GETGOODHEL (#630): COMPARE EACH HELICITY CONTRIBUTION TO 0
          {
            //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
            isGoodHel[ihel] = true;
          }
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
          const int ievt2 = ievt00 + ieppV + neppV;
          if( allMEs[ievt2] != 0 ) // NEW IMPLEMENTATION OF GETGOODHEL (#630): COMPARE EACH HELICITY CONTRIBUTION TO 0
          {
            //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
            isGoodHel[ihel] = true;
          }
#endif
        }
      }
    }
  }
#endif

  //--------------------------------------------------------------------------

  int                                          // output: nGoodHel (the number of good helicity combinations out of ncomb)
  sigmaKin_setGoodHel( const bool* isGoodHel ) // input: isGoodHel[ncomb] - host array (CUDA and C++)
  {
    int nGoodHel = 0;
    int goodHel[ncomb] = { 0 }; // all zeros https://en.cppreference.com/w/c/language/array_initialization#Notes
    for( int ihel = 0; ihel < ncomb; ihel++ )
    {
      //std::cout << "sigmaKin_setGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
      if( isGoodHel[ihel] )
      {
        goodHel[nGoodHel] = ihel;
        nGoodHel++;
      }
    }
#ifdef MGONGPUCPP_GPUIMPL
    gpuMemcpyToSymbol( dcNGoodHel, &nGoodHel, sizeof( int ) );
    gpuMemcpyToSymbol( dcGoodHel, goodHel, ncomb * sizeof( int ) );
#endif
    cNGoodHel = nGoodHel;
    for( int ihel = 0; ihel < ncomb; ihel++ ) cGoodHel[ihel] = goodHel[ihel];
    return nGoodHel;
  }

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL /* clang-format off */
  __global__ void
  normalise_output( fptype* allMEs,                    // output: allMEs[nevt], |M|^2 running_sum_over_helicities
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                    fptype* ghelAllNumerators,         // input/tmp: allNumerators super-buffer for nGoodHel <= ncomb individual helicities (index is ighel)
                    fptype* ghelAllDenominators,       // input/tmp: allNumerators super-buffer for nGoodHel <= ncomb individual helicities (index is ighel)
                    const unsigned int* allChannelIds, // input: multichannel channelIds[nevt] (1 to #diagrams); nullptr to disable SDE enhancement (fix #899/#911)
#endif
                    const fptype globaldenom ) /* clang-format on */
  {
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread)
    allMEs[ievt] /= globaldenom;
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    const int nevt = gridDim.x * blockDim.x;
    if( allChannelIds != nullptr ) // fix segfault #892 (not 'channelIds[0] != 0')
    {
      fptype* totAllNumerators = ghelAllNumerators;     // reuse "helicity #0" buffer to compute the total over all helicities
      fptype* totAllDenominators = ghelAllDenominators; // reuse "helicity #0" buffer to compute the total over all helicities
      for( int ighel = 1; ighel < dcNGoodHel; ighel++ ) // NB: the loop starts at ighel=1
      {
        fptype* hAllNumerators = ghelAllNumerators + ighel * nevt;
        fptype* hAllDenominators = ghelAllDenominators + ighel * nevt;
        totAllNumerators[ievt] += hAllNumerators[ievt];
        totAllDenominators[ievt] += hAllDenominators[ievt];
      }
      allMEs[ievt] *= totAllNumerators[ievt] / totAllDenominators[ievt];
    }
#endif
    return;
  }
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  __global__ void
  add_and_select_hel( int* allselhel,          // output: helicity selection[nevt]
                      const fptype* allrndhel, // input: random numbers[nevt] for helicity selection
                      fptype* ghelAllMEs,      // input/tmp: allMEs for nGoodHel <= ncomb individual/runningsum helicities (index is ighel)
                      fptype* allMEs,          // output: allMEs[nevt], final sum over helicities
                      const int nevt )         // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
  {
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread)
    // Compute the sum of MEs over all good helicities (defer this after the helicity loop to avoid breaking streams parall>
    for( int ighel = 0; ighel < dcNGoodHel; ighel++ )
    {
      allMEs[ievt] += ghelAllMEs[ighel * nevt + ievt];
      ghelAllMEs[ighel * nevt + ievt] = allMEs[ievt]; // reuse the buffer to store the running sum for helicity selection
    }
    // Event-by-event random choice of helicity #403
    //printf( "select_hel: ievt=%4d rndhel=%f\n", ievt, allrndhel[ievt] );
    for( int ighel = 0; ighel < dcNGoodHel; ighel++ )
    {
      if( allrndhel[ievt] < ( ghelAllMEs[ighel * nevt + ievt] / allMEs[ievt] ) )
      {
        const int ihelF = dcGoodHel[ighel] + 1; // NB Fortran [1,ncomb], cudacpp [0,ncomb-1]
        allselhel[ievt] = ihelF;
        //printf( "select_hel: ievt=%4d ihel=%4d\n", ievt, ihelF );
        break;
      }
    }
    return;
  }
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
  __global__ void
  select_col( int* allselcol,                    // output: color selection[nevt]
              const fptype* allrndcol,           // input: random numbers[nevt] for color selection
              const unsigned int* allChannelIds, // input: multichannel channelIds[nevt] (1 to #diagrams); nullptr to disable SDE enhancement (fix #899/#911)
              const int* allIgraph,              // input: per-event MLM graph (0 = no MLM); nullptr if no MLM
              const fptype_sv* allJamp2s,        // input: jamp2[ncolor][nevt] for color choice (nullptr if disabled)
              const int nevt )                   // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
  {
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread)
    // SCALAR channelId for the current event (CUDA)
    unsigned int channelId = gpu_channelId( allChannelIds );
    // Per-event MLM graph (0 = no MLM)
    const int igraph = ( allIgraph != nullptr ) ? allIgraph[ievt] : 0;
    // Event-by-event random choice of color #402
    if( channelId != 0 || igraph != 0 ) // no event-by-event choice of color if both channelId and igraph are 0 (fix FPE #783)
    {
      // Determine the jamp2 for this event (TEMPORARY? could do this with a dedicated memory accessor instead...)
      fptype_sv jamp2_sv[ncolor] = { 0 };
      assert( allJamp2s != nullptr ); // sanity check
      using J2_ACCESS = DeviceAccessJamp2;
      for( int icolC = 0; icolC < ncolor; icolC++ )
        jamp2_sv[icolC] = J2_ACCESS::kernelAccessIcolConst( allJamp2s, icolC );
      // Determine iconfig: use per-event MLM graph if provided, otherwise use channel2iconfig
      int iconfig;
      if( igraph != 0 )
      {
        iconfig = igraph; // use MLM-matched graph directly as iconfig (F-indexed, 1-based)
      }
      else
      {
        if( channelId > mgOnGpu::nchannels )
        {
          printf( "INTERNAL ERROR! Cannot choose an event-by-event random color for channelId=%d which is greater than nchannels=%d\n", channelId, mgOnGpu::nchannels );
          assert( channelId <= mgOnGpu::nchannels ); // SANITY CHECK #919 #910
        }
        // NB (see #877): in the array channel2iconfig, the input index uses C indexing (channelId -1), the output index uses F indexing (iconfig)
        // NB (see #917): mgOnGpu::channel2iconfig returns an int (which may be -1), not an unsigned int!
        iconfig = mgOnGpu::channel2iconfig[channelId - 1]; // map N_diagrams to N_config <= N_diagrams configs (fix LHE color mismatch #856: see also #826, #852, #853)
        if( iconfig <= 0 )
        {
          printf( "INTERNAL ERROR! Cannot choose an event-by-event random color for channelId=%d which has no associated SDE iconfig\n", channelId );
          assert( iconfig > 0 ); // SANITY CHECK #917
        }
        else if( iconfig > (int)mgOnGpu::nconfigSDE )
        {
          printf( "INTERNAL ERROR! Cannot choose an event-by-event random color for channelId=%d (invalid SDE iconfig=%d\n > nconfig=%d)", channelId, iconfig, mgOnGpu::nconfigSDE );
          assert( iconfig <= (int)mgOnGpu::nconfigSDE ); // SANITY CHECK #917
        }
      }
      fptype targetamp[ncolor] = { 0 };
      // NB (see #877): explicitly use 'icolC' rather than 'icol' to indicate that icolC uses C indexing in [0, N_colors-1]
      for( int icolC = 0; icolC < ncolor; icolC++ )
      {
        if( icolC == 0 )
          targetamp[icolC] = 0;
        else
          targetamp[icolC] = targetamp[icolC - 1];
        // NB (see #877): in the array icolamp, the input index uses C indexing (iconfig -1)
        if( mgOnGpu::icolamp[iconfig - 1][icolC] ) targetamp[icolC] += jamp2_sv[icolC];
      }
      //printf( "sigmaKin: ievt=%4d rndcol=%f\n", ievt, allrndcol[ievt] );
      for( int icolC = 0; icolC < ncolor; icolC++ )
      {
        if( allrndcol[ievt] < ( targetamp[icolC] / targetamp[ncolor - 1] ) )
        {
          allselcol[ievt] = icolC + 1; // NB Fortran [1,ncolor], cudacpp [0,ncolor-1]
          //printf( "sigmaKin: ievt=%d icol=%d\n", ievt, icolC+1 );
          break;
        }
      }
    }
    else
    {
      allselcol[ievt] = 0; // no color selected in Fortran range [1,ncolor] if both channelId and igraph are 0 (see #931)
    }
    return;
  }
#endif
#endif

  //--------------------------------------------------------------------------
  // Evaluate |M|^2, part independent of incoming flavour

  void /* clang-format off */
  sigmaKin( const fptype* allmomenta,           // input: momenta[nevt*npar*4]
            const fptype* allcouplings,         // input: couplings[nevt*ndcoup*2]
            const fptype* allrndhel,            // input: random numbers[nevt] for helicity selection
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
            const fptype* allrndcol,            // input: random numbers[nevt] for color selection
            const unsigned int* allChannelIds,  // input: multichannel channelIds[nevt] (1 to #diagrams); nullptr to disable single-diagram enhancement (fix #899/#911)
            const int* allIgraph,               // input: per-event MLM graph (0 = no MLM); nullptr if no MLM
#endif
            fptype* allMEs,                     // output: allMEs[nevt], |M|^2 final_avg_over_helicities
            int* allselhel,                     // output: helicity selection[nevt]
#ifdef MGONGPUCPP_GPUIMPL
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
            int* allselcol,                     // output: helicity selection[nevt]
            fptype* colAllJamp2s,               // tmp: allJamp2s super-buffer for ncolor individual colors, running sum over colors and helicities
            fptype* ghelAllNumerators,          // tmp: allNumerators super-buffer for nGoodHel <= ncomb individual helicities (index is ighel)
            fptype* ghelAllDenominators,        // tmp: allDenominators super-buffer for nGoodHel <= ncomb individual helicities (index is ighel)
#endif
            fptype* ghelAllMEs,                 // tmp: allMEs super-buffer for nGoodHel <= ncomb individual helicities (index is ighel)
            fptype* ghelAllJamps,               // tmp: jamp[2*ncolor*nGoodHel*nevt] super-buffer for nGoodHel <= ncomb individual helicities
            fptype2* ghelAllBlasTmp,            // tmp: allBlasTmp super-buffer for nGoodHel <= ncomb individual helicities (index is ighel)
            gpuBlasHandle_t* pBlasHandle,       // input: cuBLAS/hipBLAS handle
            gpuStream_t* ghelStreams,           // input: cuda streams (index is ighel: only the first nGoodHel <= ncomb are non-null)
            const int gpublocks,                // input: cuda gpublocks
            const int gputhreads                // input: cuda gputhreads
#else
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
            int* allselcol,                     // output: helicity selection[nevt]
            fptype* allNumerators,              // tmp: multichannel numerators[nevt], running_sum_over_helicities
            fptype* allDenominators,            // tmp: multichannel denominators[nevt], running_sum_over_helicities
#endif
            const int nevt                      // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
            ) /* clang-format on */
  {
    mgDebugInitialise();

    // SANITY CHECKS for cudacpp code generation (see issues #272 and #343 and PRs #619, #626, #360, #396 and #754)
    // These variable are not used anywhere else in the code and their scope is limited to this sanity check
    {
      // nprocesses == 2 may happen for "mirror processes" such as P0_uux_ttx within pp_tt012j (see PR #754)
      constexpr int nprocesses = 1;
      static_assert( nprocesses == 1 || nprocesses == 2, "Assume nprocesses == 1 or 2" );
      constexpr int process_id = 1; // code generation source: standalone_cudacpp
      static_assert( process_id == 1, "Assume process_id == 1" );
    }

    // Denominators: spins, colors and identical particles
    constexpr int helcolDenominators[1] = { 8 }; // assume nprocesses == 1 (#272 and #343)

#ifndef MGONGPUCPP_GPUIMPL
    //assert( (size_t)(allmomenta) % mgOnGpu::cppAlign == 0 ); // SANITY CHECK: require SIMD-friendly alignment [COMMENT OUT TO TEST MISALIGNED ACCESS]
    //assert( (size_t)(allMEs) % mgOnGpu::cppAlign == 0 ); // SANITY CHECK: require SIMD-friendly alignment [COMMENT OUT TO TEST MISALIGNED ACCESS]
    using E_ACCESS = HostAccessMatrixElements; // non-trivial access: buffer includes all events
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    using NUM_ACCESS = HostAccessNumerators;   // non-trivial access: buffer includes all events
    using DEN_ACCESS = HostAccessDenominators; // non-trivial access: buffer includes all events
    using CID_ACCESS = HostAccessChannelIds;   // non-trivial access: buffer includes all events
#endif
#endif

    // Start sigmaKin_lines
#include "GpuAbstraction.h"

    // === PART 0 - INITIALISATION (before calculate_jamps) ===
    // Reset the "matrix elements" - running sums of |M|^2 over helicities for the given event
#ifdef MGONGPUCPP_GPUIMPL
    // *** PART 0a - CUDA ***
    const int nevt = gpublocks * gputhreads;
    gpuMemset( allMEs, 0, nevt * sizeof( fptype ) );
    // all function parameter are type int, the first is casted to size_t to avoid integer overflow if the product is greater than max(int)
    gpuMemset( ghelAllJamps, 0, static_cast<size_t>( cNGoodHel ) * ncolor * mgOnGpu::nx2 * nevt * sizeof( fptype ) );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    gpuMemset( colAllJamp2s, 0, ncolor * nevt * sizeof( fptype ) );
    gpuMemset( ghelAllNumerators, 0, cNGoodHel * nevt * sizeof( fptype ) );
    gpuMemset( ghelAllDenominators, 0, cNGoodHel * nevt * sizeof( fptype ) );
#endif
    gpuMemset( ghelAllMEs, 0, cNGoodHel * nevt * sizeof( fptype ) );
#else
    // *** PART 0b - C++ ***
    const int npagV = nevt / neppV;
    for( int ipagV = 0; ipagV < npagV; ++ipagV )
    {
      const int ievt0 = ipagV * neppV;
      fptype* MEs = E_ACCESS::ieventAccessRecord( allMEs, ievt0 );
      fptype_sv& MEs_sv = E_ACCESS::kernelAccess( MEs );
      MEs_sv = fptype_sv{ 0 };
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      fptype* numerators = NUM_ACCESS::ieventAccessRecord( allNumerators, ievt0 );
      fptype* denominators = DEN_ACCESS::ieventAccessRecord( allDenominators, ievt0 );
      fptype_sv& numerators_sv = NUM_ACCESS::kernelAccess( numerators );
      fptype_sv& denominators_sv = DEN_ACCESS::kernelAccess( denominators );
      numerators_sv = fptype_sv{ 0 };
      denominators_sv = fptype_sv{ 0 };
#endif
    }
#endif

    // === PART 1 - HELICITY LOOP: CALCULATE WAVEFUNCTIONS ===
    // (in both CUDA and C++, using precomputed good helicities)

#ifdef MGONGPUCPP_GPUIMPL // CUDA OR C++

    // *** START OF PART 1a - CUDA (one event per GPU thread) ***
    // Use CUDA/HIP streams to process different helicities in parallel (one good helicity per stream)
    // (1) First, within each helicity stream, compute the QCD partial amplitudes jamp's for each helicity
    // In multichannel mode, also compute the running sums over helicities of numerators, denominators and squared jamp2s
    for( int ighel = 0; ighel < cNGoodHel; ighel++ )
    {
      const int ihel = cGoodHel[ighel];
      fptype* hAllJamps = ghelAllJamps + ighel * nevt; // HACK: bypass DeviceAccessJamp (consistent with layout defined there)
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      fptype* hAllNumerators = ghelAllNumerators + ighel * nevt;
      fptype* hAllDenominators = ghelAllDenominators + ighel * nevt;
      gpuLaunchKernelStream( calculate_jamps, gpublocks, gputhreads, ghelStreams[ighel], ihel, allmomenta, allcouplings, hAllJamps, allChannelIds, hAllNumerators, hAllDenominators, colAllJamp2s, nevt );
#else
      gpuLaunchKernelStream( calculate_jamps, gpublocks, gputhreads, ghelStreams[ighel], ihel, allmomenta, allcouplings, hAllJamps, nevt );
#endif
    }
    // (2) Then compute the ME for that helicity from the color sum of QCD partial amplitudes jamps
    color_sum_gpu( ghelAllMEs, ghelAllJamps, ghelAllBlasTmp, pBlasHandle, ghelStreams, cNGoodHel, gpublocks, gputhreads );
    checkGpu( gpuDeviceSynchronize() ); // do not start helicity/color selection until the loop over helicities has completed
    // (3) Wait for all helicity streams to complete, then finally compute the ME sum over all helicities and choose one helicity and one color
    // Event-by-event random choice of helicity #403 and ME sum over helicities (defer this after the helicity loop to avoid breaking streams parallelism)
    gpuLaunchKernel( add_and_select_hel, gpublocks, gputhreads, allselhel, allrndhel, ghelAllMEs, allMEs, gpublocks * gputhreads );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Event-by-event random choice of color #402
    gpuLaunchKernel( select_col, gpublocks, gputhreads, allselcol, allrndcol, allChannelIds, allIgraph, colAllJamp2s, gpublocks * gputhreads );
#endif
    // *** END OF PART 1a - CUDA (one event per GPU thread) ***

#else // CUDA OR C++

    // *** START OF PART 1b - C++ (loop on event pages)
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
    // Mixed fptypes #537: float for color algebra and double elsewhere
    // Delay color algebra and ME updates (only on even pages)
    assert( npagV % 2 == 0 );     // SANITY CHECK for mixed fptypes: two neppV-pages are merged to one 2*neppV-page
    const int npagV2 = npagV / 2; // loop on two SIMD pages (neppV events) at a time
#else
    const int npagV2 = npagV;            // loop on one SIMD page (neppV events) at a time
#endif
#ifdef _OPENMP
    // OMP multithreading #575 (NB: tested only with gcc11 so far)
    // See https://www.openmp.org/specifications/
    // - default(none): no variables are shared by default
    // - shared: as the name says
    // - private: give each thread its own copy, without initialising
    // - firstprivate: give each thread its own copy, and initialise with value from outside
#define _OMPLIST0 allcouplings, allMEs, allmomenta, allrndcol, allrndhel, allselcol, allselhel, cGoodHel, cNGoodHel, npagV2
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
#define _OMPLIST1 , allDenominators, allNumerators, allChannelIds, allIgraph, mgOnGpu::icolamp, mgOnGpu::channel2iconfig
#else
#define _OMPLIST1
#endif
#pragma omp parallel for default( none ) shared( _OMPLIST0 _OMPLIST1 )
#undef _OMPLIST0
#undef _OMPLIST1
#endif // _OPENMP
    for( int ipagV2 = 0; ipagV2 < npagV2; ++ipagV2 )
    {
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
      const int ievt00 = ipagV2 * neppV * 2; // loop on two SIMD pages (neppV events) at a time
#else
      const int ievt00 = ipagV2 * neppV; // loop on one SIMD page (neppV events) at a time
#endif
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // SCALAR channelId for the whole SIMD neppV2 event page (C++), i.e. one or two neppV event page(s)
      // The cudacpp implementation ASSUMES (and checks! #898) that all channelIds are the same in a neppV2 SIMD event page
      // **NB! in "mixed" precision, using SIMD, calculate_jamps computes MEs for TWO neppV pages with a single channelId! #924
      unsigned int channelId = 0; // disable multichannel single-diagram enhancement unless allChannelIds != nullptr
      if( allChannelIds != nullptr )
      {
        // First - and/or only - neppV page of channels (iParity=0 => ievt0 = ievt00 + 0 * neppV)
        const unsigned int* channelIds = CID_ACCESS::ieventAccessRecordConst( allChannelIds, ievt00 ); // fix bug #899/#911
        uint_sv channelIds_sv = CID_ACCESS::kernelAccessConst( channelIds );                           // fix #895 (compute this only once for all diagrams)
#ifndef MGONGPU_CPPSIMD
        // NB: channelIds_sv is a scalar in no-SIMD C++
        channelId = channelIds_sv;
#else
        // NB: channelIds_sv is a vector in SIMD C++
        channelId = channelIds_sv[0];    // element[0]
        for( int i = 1; i < neppV; ++i ) // elements[1...neppV-1]
        {
          assert( channelId == channelIds_sv[i] ); // SANITY CHECK #898: check that all events in a SIMD vector have the same channelId
        }
#endif
        assert( channelId > 0 ); // SANITY CHECK: scalar channelId must be > 0 if multichannel is enabled (allChannelIds != nullptr)
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
        // Second neppV page of channels (iParity=1 => ievt0 = ievt00 + 1 * neppV)
        const unsigned int* channelIds2 = CID_ACCESS::ieventAccessRecordConst( allChannelIds, ievt00 + neppV ); // fix bug #899/#911
        uint_v channelIds2_v = CID_ACCESS::kernelAccessConst( channelIds2 );                                    // fix #895 (compute this only once for all diagrams)
        // **NB! in "mixed" precision, using SIMD, calculate_jamps computes MEs for TWO neppV pages with a single channelId! #924
        for( int i = 0; i < neppV; ++i )
        {
          assert( channelId == channelIds2_v[i] ); // SANITY CHECKS #898 #924: all events in the 2nd SIMD vector have the same channelId as that of the 1st SIMD vector
        }
#endif
      }
#endif
      // Running sum of partial amplitudes squared for event by event color selection (#402)
      // (jamp2[nParity][ncolor][neppV] for the SIMD vector - or the two SIMD vectors - of events processed in calculate_jamps)
      fptype_sv jamp2_sv[nParity * ncolor] = {};
      fptype_sv MEs_ighel[ncomb] = {};  // sum of MEs for all good helicities up to ighel (for the first - and/or only - neppV page)
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
      fptype_sv MEs_ighel2[ncomb] = {}; // sum of MEs for all good helicities up to ighel (for the second neppV page)
#endif
      for( int ighel = 0; ighel < cNGoodHel; ighel++ )
      {
        const int ihel = cGoodHel[ighel];
        cxtype_sv jamp_sv[nParity * ncolor] = {}; // fixed nasty bug (omitting 'nParity' caused memory corruptions after calling calculate_jamps)
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
        // **NB! in "mixed" precision, using SIMD, calculate_jamps computes MEs for TWO neppV pages with a single channelId! #924
        calculate_jamps( ihel, allmomenta, allcouplings, jamp_sv, channelId, allNumerators, allDenominators, jamp2_sv, ievt00 );
#else
        calculate_jamps( ihel, allmomenta, allcouplings, jamp_sv, ievt00 );
#endif
        color_sum_cpu( allMEs, jamp_sv, ievt00 );
        MEs_ighel[ighel] = E_ACCESS::kernelAccess( E_ACCESS::ieventAccessRecord( allMEs, ievt00 ) );
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
        MEs_ighel2[ighel] = E_ACCESS::kernelAccess( E_ACCESS::ieventAccessRecord( allMEs, ievt00 + neppV ) );
#endif
      }
      // Event-by-event random choice of helicity #403
      for( int ieppV = 0; ieppV < neppV; ++ieppV )
      {
        const int ievt = ievt00 + ieppV;
        //printf( "sigmaKin: ievt=%4d rndhel=%f\n", ievt, allrndhel[ievt] );
        for( int ighel = 0; ighel < cNGoodHel; ighel++ )
        {
#if defined MGONGPU_CPPSIMD
          //printf( "sigmaKin: ievt=%4d ighel=%d MEs_ighel=%f\n", ievt, ighel, MEs_ighel[ighel][ieppV] );
          const bool okhel = allrndhel[ievt] < ( MEs_ighel[ighel][ieppV] / MEs_ighel[cNGoodHel - 1][ieppV] );
#else
          //printf( "sigmaKin: ievt=%4d ighel=%d MEs_ighel=%f\n", ievt, ighel, MEs_ighel[ighel] );
          const bool okhel = allrndhel[ievt] < ( MEs_ighel[ighel] / MEs_ighel[cNGoodHel - 1] );
#endif
          if( okhel )
          {
            const int ihelF = cGoodHel[ighel] + 1; // NB Fortran [1,ncomb], cudacpp [0,ncomb-1]
            allselhel[ievt] = ihelF;
            //printf( "sigmaKin: ievt=%4d ihel=%4d\n", ievt, ihelF );
            break;
          }
        }
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
        const int ievt2 = ievt00 + ieppV + neppV;
        //printf( "sigmaKin: ievt=%4d rndhel=%f\n", ievt2, allrndhel[ievt2] );
        for( int ighel = 0; ighel < cNGoodHel; ighel++ )
        {
          //printf( "sigmaKin: ievt=%4d ighel=%d MEs_ighel=%f\n", ievt2, ighel, MEs_ighel2[ighel][ieppV] );
          if( allrndhel[ievt2] < ( MEs_ighel2[ighel][ieppV] / MEs_ighel2[cNGoodHel - 1][ieppV] ) )
          {
            const int ihelF = cGoodHel[ighel] + 1; // NB Fortran [1,ncomb], cudacpp [0,ncomb-1]
            allselhel[ievt2] = ihelF;
            //printf( "sigmaKin: ievt=%4d ihel=%4d\n", ievt2, ihelF );
            break;
          }
        }
#endif
      }
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL // multichannel enabled (random color choice)
      // Event-by-event random choice of color #402
      // NB: with MLM, different events in a SIMD page may have different igraph values, so iconfig must be per-event
      for( int ieppV = 0; ieppV < neppV; ++ieppV )
      {
        const int ievt = ievt00 + ieppV;
        // Use per-event MLM graph if provided, otherwise use channel2iconfig
        const int igraph1_ievt = ( allIgraph != nullptr ) ? allIgraph[ievt] : 0; // per-event igraph (may differ across SIMD page with MLM)
        if( channelId != 0 || igraph1_ievt != 0 )                                // no event-by-event choice of color if both channelId and igraph are 0 (fix FPE #783)
        {
          // Determine iconfig: use per-event MLM graph if provided, otherwise use channel2iconfig
          int iconfig;
          if( igraph1_ievt != 0 )
          {
            iconfig = igraph1_ievt; // use MLM-matched graph directly as iconfig (F-indexed, 1-based)
          }
          else
          {
            if( channelId > mgOnGpu::nchannels )
            {
              printf( "INTERNAL ERROR! Cannot choose an event-by-event random color for channelId=%d which is greater than nchannels=%d\n", channelId, mgOnGpu::nchannels );
              assert( channelId <= mgOnGpu::nchannels ); // SANITY CHECK #919 #910
            }
            // NB (see #877): in the array channel2iconfig, the input index uses C indexing (channelId -1), the output index uses F indexing (iconfig)
            // NB (see #917): mgOnGpu::channel2iconfig returns an int (which may be -1), not an unsigned int!
            iconfig = mgOnGpu::channel2iconfig[channelId - 1]; // map N_diagrams to N_config <= N_diagrams configs (fix LHE color mismatch #856: see also #826, #852, #853)
            if( iconfig <= 0 )
            {
              printf( "INTERNAL ERROR! Cannot choose an event-by-event random color for channelId=%d which has no associated SDE iconfig\n", channelId );
              assert( iconfig > 0 ); // SANITY CHECK #917
            }
            else if( iconfig > (int)mgOnGpu::nconfigSDE )
            {
              printf( "INTERNAL ERROR! Cannot choose an event-by-event random color for channelId=%d (invalid SDE iconfig=%d\n > nconfig=%d)", channelId, iconfig, mgOnGpu::nconfigSDE );
              assert( iconfig <= (int)mgOnGpu::nconfigSDE ); // SANITY CHECK #917
            }
          }
          // NB (see #877): explicitly use 'icolC' rather than 'icol' to indicate that icolC uses C indexing in [0, N_colors-1]
          // NB: targetamp is a scalar fptype (not fptype_sv) - iconfig is per-event so we extract the scalar lane from jamp2_sv
          fptype targetamp[ncolor] = { 0 };
          for( int icolC = 0; icolC < ncolor; icolC++ )
          {
            if( icolC == 0 )
              targetamp[icolC] = fptype{ 0 };
            else
              targetamp[icolC] = targetamp[icolC - 1];
              // NB (see #877): in the array icolamp, the input index uses C indexing (iconfig -1)
#if defined MGONGPU_CPPSIMD
            if( mgOnGpu::icolamp[iconfig - 1][icolC] ) targetamp[icolC] += jamp2_sv[icolC][ieppV];
#else
            if( mgOnGpu::icolamp[iconfig - 1][icolC] ) targetamp[icolC] += jamp2_sv[icolC];
#endif
          }
          //printf( "sigmaKin: ievt=%4d rndcol=%f\n", ievt, allrndcol[ievt] );
          for( int icolC = 0; icolC < ncolor; icolC++ )
          {
            const bool okcol = allrndcol[ievt] < ( targetamp[icolC] / targetamp[ncolor - 1] );
            if( okcol )
            {
              allselcol[ievt] = icolC + 1; // NB Fortran [1,ncolor], cudacpp [0,ncolor-1]
              //printf( "sigmaKin: ievt=%d icol=%d\n", ievt, icolC+1 );
              break;
            }
          }
        }
        else
        {
          allselcol[ievt] = 0; // no color selected in Fortran range [1,ncolor] if both channelId and igraph are 0 (see #931)
        }
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
        const int ievt2 = ievt00 + ieppV + neppV;
        const int igraph1_ievt2 = ( allIgraph != nullptr ) ? allIgraph[ievt2] : 0; // per-event igraph (may differ across SIMD page with MLM)
        if( channelId != 0 || igraph1_ievt2 != 0 )                                 // no event-by-event choice of color if both channelId and igraph are 0 (fix FPE #783)
        {
          // Determine iconfig2: use per-event MLM graph if provided, otherwise use channel2iconfig
          int iconfig2;
          if( igraph1_ievt2 != 0 )
          {
            iconfig2 = igraph1_ievt2; // use MLM-matched graph directly as iconfig (F-indexed, 1-based)
          }
          else
          {
            iconfig2 = mgOnGpu::channel2iconfig[channelId - 1]; // same channelId as for ievt (sanity checks already done above)
          }
          fptype targetamp2[ncolor] = { 0 };
          for( int icolC = 0; icolC < ncolor; icolC++ )
          {
            if( icolC == 0 )
              targetamp2[icolC] = fptype{ 0 };
            else
              targetamp2[icolC] = targetamp2[icolC - 1];
            // NB (see #877): in the array icolamp, the input index uses C indexing (iconfig -1)
            if( mgOnGpu::icolamp[iconfig2 - 1][icolC] ) targetamp2[icolC] += jamp2_sv[ncolor + icolC][ieppV];
          }
          //printf( "sigmaKin: ievt=%4d rndcol=%f\n", ievt2, allrndcol[ievt2] );
          for( int icolC = 0; icolC < ncolor; icolC++ )
          {
            if( allrndcol[ievt2] < ( targetamp2[icolC] / targetamp2[ncolor - 1] ) )
            {
              allselcol[ievt2] = icolC + 1; // NB Fortran [1,ncolor], cudacpp [0,ncolor-1]
              //printf( "sigmaKin: ievt2=%d icol=%d\n", ievt2, icolC+1 );
              break;
            }
          }
        }
        else
        {
          allselcol[ievt2] = 0; // no color selected in Fortran range [1,ncolor] if channelId == 0 (see #931)
        }
#endif
      }
#endif // multichannel enabled (random color choice)
    }
    // *** END OF PART 1b - C++ (loop on event pages)

#endif // CUDA or C++

    // PART 2 - FINALISATION (after calculate_jamps)
    // Get the final |M|^2 as an average over helicities/colors of the running sum of |M|^2 over helicities for the given event
    // [NB 'sum over final spins, average over initial spins', eg see
    // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
#ifdef MGONGPUCPP_GPUIMPL
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    gpuLaunchKernel( normalise_output, gpublocks, gputhreads, allMEs, ghelAllNumerators, ghelAllDenominators, allChannelIds, helcolDenominators[0] );
#else
    gpuLaunchKernel( normalise_output, gpublocks, gputhreads, allMEs, helcolDenominators[0] );
#endif
#else
    for( int ipagV = 0; ipagV < npagV; ++ipagV )
    {
      const int ievt0 = ipagV * neppV;
      fptype* MEs = E_ACCESS::ieventAccessRecord( allMEs, ievt0 );
      fptype_sv& MEs_sv = E_ACCESS::kernelAccess( MEs );
      MEs_sv /= helcolDenominators[0];
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( allChannelIds != nullptr ) // fix segfault #892 (not 'channelIds[0] != 0')
      {
        fptype* numerators = NUM_ACCESS::ieventAccessRecord( allNumerators, ievt0 );
        fptype* denominators = DEN_ACCESS::ieventAccessRecord( allDenominators, ievt0 );
        fptype_sv& numerators_sv = NUM_ACCESS::kernelAccess( numerators );
        fptype_sv& denominators_sv = DEN_ACCESS::kernelAccess( denominators );
        MEs_sv *= numerators_sv / denominators_sv;
      }
#endif
      //for( int ieppV = 0; ieppV < neppV; ieppV++ )
      //{
      //  const unsigned int ievt = ipagV * neppV + ieppV;
      //  printf( "sigmaKin: ievt=%2d me=%f\n", ievt, allMEs[ievt] );
      //}
    }
#endif
    mgDebugFinalise();
  }

  //--------------------------------------------------------------------------

} // end namespace

//==========================================================================
