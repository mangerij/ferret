#ifndef FCUDAGLOBAL_HPP
#define FCUDAGLOBAL_HPP

#include "../../Utils/FGlobal.hpp"

// Manage special case for nvcc
#if defined(__CUDACC__) || defined(__NVCC__)
#else
#endif

#include <cuda.h>

#endif // FCUDAGLOBAL_HPP

