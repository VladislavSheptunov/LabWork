#ifndef _SETTINGS_H_
#define _SETTINGS_H_

// Pre-defined C/C++ Compiler Macros
#ifdef __INTEL_COMPILER
  #include "immintrin.h"

  #define ALIGNMENT
  #define VECTORIZATION
  #define PARALLELIZATION   1
#else
  #define PARALLELIZATION   1
#endif // __INTEL_COMPILER

#if PARALLELIZATION
  #define DEFAULT_SETTING   1
    #if !DEFAULT_SETTING
      #define COUNT_THREAD  7
    #endif // !DEFAULT_SETTING
#endif // PARALLELIZATION

#endif //_SETTINGS_H_
