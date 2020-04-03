cmake_minimum_required (VERSION 2.8.3 FATAL_ERROR)

if(CMAKE_MAJOR_VERSION GREATER 1 AND CMAKE_MINOR_VERSION GREATER 7 AND CMAKE_PATCH_VERSION GREATER 11)
  cmake_policy(SET CMP0022 OLD)
endif()

cmake_policy(SET CMP0004 OLD)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

#===========================================================================
# Project Declaration
#===========================================================================
project(SCALFMM C CXX)

# check if compiling into source directories
string(COMPARE EQUAL "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}" insource)
if(insource)
  message(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto ./Build and tapes cmake ../")
endif(insource)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/)
set(SCALFMM_CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)
#
# Adds the CMAKE_DEPENDENT_OPTION command
INCLUDE(CMakeDependentOption)
#  Add to check CPU info
include(GetCpuInfos)
GetCpuInfos()

#===========================================================================
# Version Number
#===========================================================================

# SCALFMM version number.  An even minor number corresponds to releases.
set(SCALFMM_MAJOR_VERSION 1)
set(SCALFMM_MINOR_VERSION 4)
set(SCALFMM_PATCH_VERSION 0)
set(SCALFMM_VERSION  "${SCALFMM_MAJOR_VERSION}.${SCALFMM_MINOR_VERSION}.${SCALFMM_PATCH_VERSION}" )

set( MORSE_DISTRIB_DIR "" CACHE PATH "Directory of MORSE distribution")

if (MORSE_DISTRIB_DIR OR EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/morse/")

  # Add extra cmake module path and initialize morse cmake modules
  # --------------------------------------------------------------
  if(MORSE_DISTRIB_DIR)
    list(APPEND CMAKE_MODULE_PATH ${MORSE_DISTRIB_DIR}/cmake_modules)
    list(APPEND CMAKE_MODULE_PATH "${MORSE_DISTRIB_DIR}/cmake_modules/morse")
    set(MORSE_CMAKE_MODULE_PATH ${MORSE_DISTRIB_DIR}/cmake_modules/morse )
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/morse/")
    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/morse/)
    set(MORSE_CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/morse )
  endif()
  include(MorseInit)
  #
  # Options
  #
  option( SCALFMM_USE_MPI              "Set to ON to build ScaFMM with MPI"          OFF )
  option( SCALFMM_USE_BLAS             "Set to ON to build ScaFMM with BLAS"         ON )
  option( SCALFMM_USE_FFT              "Set to ON to build ScaFMM with FFTW"         ON )
  option( SCALFMM_BUILD_EXAMPLES       "Set to ON to compile examples "              ON  )
  option( SCALFMM_BUILD_TESTS          "Set to ON to build functionalities Tests"    OFF )
  option( SCALFMM_BUILD_UTESTS         "Set to ON to build UTests"                   OFF )
  option( SCALFMM_INSTALL_DATA         "Set to ON to include the data when proceeding to the installtion " OFF )
  option( SCALFMM_BUILD_DOC            "Set to ON to build the Doxygen documentation " OFF )
  option( SCALFMM_BUILD_DEBUG          "Set to ON to build in Debug"                 OFF )
  option( SCALFMM_USE_LOG              "Set to ON to print output debug information" OFF )
  option( SCALFMM_USE_MEM_STATS        "Set to ON to profile memory"                 OFF )
  option( SCALFMM_ATTACHE_SOURCE       "Set to ON to compile with -g"                OFF )
  option( SCALFMM_USE_ADDONS           "Set to ON to compile add ons"                OFF )
  option( SCALFMM_USE_SIGNALS          "Set to ON to catch various signal an print backtrace"  OFF )
  option( SCALFMM_USE_ASSERT           "Set to ON to enable safe tests during execution" ON  )
  option( SCALFMM_USE_MIC_NATIVE       "Set to ON to compile in native mode for MIC" OFF  )
  option( SCALFMM_ONLY_DEVEL           "Set to ON to compile Development tools (only scalfmm team)" OFF )
  option( SCALFMM_USE_EZTRACE          "Set to ON to compile with eztrace framwork" OFF )
  option( SCALFMM_USE_STARPU           "Set to ON to build SCALFMM with StarPU"     OFF )
  option( SCALFMM_BUILD_UTILS          "Set to ON to build utils Tests"             OFF )

#
#  VECTORISATION
#
  if( APPLE ) # to fix problem with  GCC and avx
    CMAKE_DEPENDENT_OPTION( SCALFMM_USE_SSE              "Set to ON to compile with SSE support (and use intrinsec SSE P2P)" ON "CPUOPTION_SSE3;NOT CPUOPTION_AVX2" OFF  )
    CMAKE_DEPENDENT_OPTION( SCALFMM_USE_AVX              "Set to ON to compile with AVX support (and use intrinsec AVX P2P)" OFF "CPUOPTION_AVX; NOT CPUOPTION_AVX2" OFF  )
  else(APPLE)
    CMAKE_DEPENDENT_OPTION( SCALFMM_USE_SSE              "Set to ON to compile with SSE support (and use intrinsec SSE P2P)" ON "CPUOPTION_SSE3;NOT CPUOPTION_AVX;NOT CPUOPTION_AVX2" OFF  )
    CMAKE_DEPENDENT_OPTION( SCALFMM_USE_AVX              "Set to ON to compile with AVX support (and use intrinsec AVX P2P)" ON "CPUOPTION_AVX; NOT CPUOPTION_AVX2" OFF  )
  endif(APPLE)
  CMAKE_DEPENDENT_OPTION( SCALFMM_USE_AVX2             "Set to ON to compile with AVX support (and use intrinsec AVX2 P2P)" ON "CPUOPTION_AVX2" OFF )

  if( SCALFMM_ONLY_DEVEL )
# OPENMP 4/5 support 
    option( OPENMP_SUPPORT_COMMUTE   "Set to ON to let tasks commute (KSTAR/StarPU compiler only)" OFF )
    option( OPENMP_SUPPORT_PRIORITY  "Set to ON to enable tasks priority (KSTAR/StarPU compiler only)" OFF )
    option( SCALFMM_DISABLE_NATIVE_OMP4 "Set to ON to disable the gcc/intel omp4"    OFF )
    option( SCALFMM_TIME_OMPTASKS "Set to ON to time omp4 tasks and generate output file"    OFF )
# STARPU options
    CMAKE_DEPENDENT_OPTION(SCALFMM_STARPU_USE_COMMUTE  "Set to ON to enable commute with StarPU" ON  "SCALFMM_USE_STARPU" OFF)
    CMAKE_DEPENDENT_OPTION(SCALFMM_STARPU_USE_REDUX  "Set to ON to enable redux with StarPU"     OFF "SCALFMM_USE_STARPU" OFF)
    CMAKE_DEPENDENT_OPTION(SCALFMM_STARPU_USE_PRIO  "Set to ON to enable priority with StarPU"   ON  "SCALFMM_USE_STARPU" OFF)
    CMAKE_DEPENDENT_OPTION(SCALFMM_STARPU_FORCE_NO_SCHEDULER  "Set to ON to disable heteroprio even if supported" OFF "SCALFMM_USE_STARPU" OFF)
  endif()
  message(STATUS "AVANT  ${CMAKE_CXX_COMPILER_ID}" )
#
#  MPI
#
  if( SCALFMM_USE_MPI )
    try_compile(COMPILE_INTEL ${CMAKE_CURRENT_BINARY_DIR}
      ${SCALFMM_CMAKE_MODULE_PATH}/compileTestIntel.cpp
      COMPILE_DEFINITIONS "${CMAKE_CXX_FLAGS}")
    if (COMPILE_INTEL)
      set(CMAKE_CXX_COMPILER_ID  "Intel")
    endif()
  endif()
  message(STATUS "CXX  ${CMAKE_CXX_COMPILER_ID}" )
  # Set scalfmm to default libraries
  set(SCALFMM_LIBRARIES "")
  set(SCALFMM_CXX_FLAGS "-std=c++11 -fpic -Wall")
  #
  #
  # Test if openmp is here
  #
  find_package (OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  else(OPENMP_FOUND)
    message(WARNING "OPENMP NOT FOUND")
  endif(OPENMP_FOUND)
  message(STATUS "SCALFMM_BUILD_DEBUG = ${SCALFMM_BUILD_DEBUG}" )
  #
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(SCALFMM_CXX_FLAGS  "${SCALFMM_CXX_FLAGS} -m64")
  endif()
  ##############################################################################
  #                           Compile options                                  #
  ##############################################################################
  #  -xHost -mfpmath=sse
  # -Wall Wnosign-conversion
  #
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    # INTEL
    IF (APPLE)
      IF( CPUOPTION_SSE42 ) 
        set(SSE_FLAGS  "-msse4  -mfpmath=sse")   # -mtune=native -march=native
      ELSEIF (CPUOPTION_SSE3)
        set(SSE_FLAGS  "-msse3  -mfpmath=sse")   # -mtune=native -march=native
      ENDIF (CPUOPTION_SSE42)  
    else(APPLE)
      set(AVX_FLAGS  "-march=native -axCORE-AVX2,CORE-AVX-I,AVX") #-mavx
      set(AVX2_FLAGS  "-march=native  -axCORE-AVX2,CORE-AVX-I") #-march=core-avx2
      set(SSE_FLAGS  "-axSSE4.2  -march=native")
    endif(APPLE)
    #-Wshadow -Wpointer-arith -Wcast-qual -Wconversion  -Wall -Wnosign-conversion ")
  else()
    # NOT      INTEL
    if(NOT SCALFMM_USE_MPI)
      include(CheckCCompilerFlag)
      check_c_compiler_flag(-Wzero-as-null-pointer-constant HAS_WZERO_NULL_PTR_FLAG)
      if(HAS_WZERO_NULL_PTR_FLAG)
        set(SCALFMM_CXX_FLAGS  "${SCALFMM_CXX_FLAGS} -Wzero-as-null-pointer-constant")
      endif()
    else()
      include(CheckCCompilerFlag)
      check_c_compiler_flag(-Wno-literal-suffix HAS_NO_LITERAL_SUFFIX_FLAG)
      if(HAS_NO_LITERAL_SUFFIX_FLAG)
        set(SCALFMM_CXX_FLAGS  "${SCALFMM_CXX_FLAGS} -Wno-literal-suffix")
      endif()
    endif()
    IF (APPLE)
      #      set(SSE_FLAGS  "-msse4  -mfpmath=sse")   # -mtune=native -march=native
      IF( CPUOPTION_SSE42 ) 
        set(SSE_FLAGS  "-msse4  -mfpmath=sse")   # -mtune=native -march=native
      ELSEIF (CPUOPTION_SSE3)
        set(SSE_FLAGS  "-msse3  -mfpmath=sse")   # -mtune=native -march=native
      ENDIF (CPUOPTION_SSE42)
      set(AVX_FLAGS "-mtune=native -march=avx")
      set(AVX2_FLAGS "-mtune=native -march=native -mmic")
    else(APPLE)
      set(SSE_FLAGS  "-mtune=native -march=native")
      set(AVX_FLAGS "-mtune=native -march=native")
      set(AVX2_FLAGS "-mtune=native -march=native -mmic")
    endif(APPLE)
  endif()

  ##############################################################################
  #                           FUSE list                                        #
  ##############################################################################
  set(FUSE_LIST "")
  # then do list(APPEND FUSE_LIST "BLAS") to protect from FUSE_BLAS
  list(APPEND FUSE_LIST "MPI")

  ##############################################################################
  #                        Debug      flags                                    #
  ##############################################################################
  #
  if( SCALFMM_BUILD_DEBUG )
    #
    # Compile Debug flags
    #
    set(CMAKE_BUILD_TYPE Debug)
    ADD_DEFINITIONS(-O0)
    if( APPLE )
      set(SCALFMM_FLAGS_OPTI_DEBUG "-funroll-loops" CACHE STRING
        "Set your optimization flags for debug mode.")
    else(APPLE)
      if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        set(SCALFMM_FLAGS_OPTI_DEBUG "-fp-model strict  -funroll-loops" CACHE STRING
          "Set your optimization flags for debug mode.")
      else()
        set(SCALFMM_FLAGS_OPTI_DEBUG "-funroll-loops" CACHE STRING
          "Set your optimization flags for debug mode.")
      endif()
    endif(APPLE)
    #  ADD_DEFINITIONS(${SCALFMM_FLAGS_OPTI_DEBUG})
    set(SCALFMM_CXX_FLAGS "${SCALFMM_CXX_FLAGS} ${SCALFMM_FLAGS_OPTI_DEBUG}")
    set(CMAKE_CXX_FLAGS_TYPE "${CMAKE_CXX_FLAGS_DEBUG}")

  else()
    ##############################################################################
    #                          Release flags                                     #
    ##############################################################################

    set(CMAKE_BUILD_TYPE Release)
    set(CMAKE_CXX_FLAGS_TYPE  "${CMAKE_CXX_FLAGS_RELEASE}")
    # if compiler is intel add -ip

    if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
      set(SCALFMM_CXX_FLAGS  "${SCALFMM_CXX_FLAGS} -ipo  -fstrict-aliasing ")# -no-prec-div -no-prec-sqrt")
    endif()
    # Compile optimization
    if(APPLE)
      set(SCALFMM_FLAGS_OPTI_RELEASE "-funroll-loops" CACHE STRING
        "Set your optimization flags for release mode.")
    else(APPLE)
      # Not apple system - Check the compiler flags
      if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        #set(SCALFMM_FLAGS_OPTI_RELEASE "-fp-model precise -fp-model source -fimf-precision=low  -funroll-loops -ftree-vectorize"
        set(SCALFMM_FLAGS_OPTI_RELEASE "-funroll-loops -ftree-vectorize"
          CACHE STRING "Set your optimization flags for release mode.")
  # set(SCALFMM_FLAGS_OPTI_RELEASE "-funroll-loops -ftree-vectorize" CACHE STRING
  #   "Set your optimization flags for release mode.")
      else()
        set(SCALFMM_FLAGS_OPTI_RELEASE "-ffast-math -funroll-loops -ftree-vectorize" CACHE STRING
          "Set your optimization flags for release mode.")
      endif()
    endif(APPLE)
    #
    set(SCALFMM_CXX_FLAGS  "${SCALFMM_CXX_FLAGS} ${SCALFMM_FLAGS_OPTI_RELEASE}")
  endif()
  MESSAGE(STATUS   " %%%%%%%%%% SCALFMM_CXX_FLAGS  ${SCALFMM_CXX_FLAGS} %%%%%%%%%%%%%")
  MESSAGE(STATUS   " %%%%%%%%%% CMAKE_CXX_FLAGS_RELEASE  ${CMAKE_CXX_FLAGS_RELEASE} %%%%%%%%%%%%%")

  #
  ##############################################################################
  #                Attach source code to exec                                  #
  ##############################################################################
  message(STATUS "SCALFMM_ATTACHE_SOURCE       = ${SCALFMM_ATTACHE_SOURCE}" )
  if( SCALFMM_ATTACHE_SOURCE )
    message(STATUS "Option -g is used, the code is attached to the binary." )
    set(SCALFMM_CXX_FLAGS  "${SCALFMM_CXX_FLAGS} -g")
  endif(SCALFMM_ATTACHE_SOURCE)
  #
  ##############################################################################
  #                                 find MPI                                   #
  ##############################################################################
  #
  if( SCALFMM_USE_MPI )

    # allows to use an external mpi compilation by setting compilers with
    # -DMPI_C_COMPILER=path/to/mpicc -DMPI_CXX_COMPILER=path/to/mpicxx
    # at cmake configure
    if(NOT MPI_C_COMPILER)
      set(MPI_C_COMPILER mpicc)
    endif()
    if(NOT MPI_CXX_COMPILER)
      set(MPI_CXX_COMPILER mpicxx)
    endif()

    find_package(MPI REQUIRED)

    if (MPI_CXX_INCLUDE_PATH)
      include_directories( ${MPI_CXX_INCLUDE_PATH} )
    endif()
    if (MPI_CXX_COMPILE_FLAGS)
      set(SCALFMM_CXX_FLAGS "${SCALFMM_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    endif()
    if (MPI_CXX_INCLUDE_PATH)
      set(SCALFMM_INCLUDES  "${SCALFMM_INCLUDES}; ${MPI_CXX_INCLUDE_PATH}")
    endif()
    if (MPI_CXX_LINK_FLAGS)
      set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES};${MPI_CXX_LINK_FLAGS}")
    endif()
    if (MPI_CXX_LIBRARIES)
      set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES};${MPI_CXX_LIBRARIES}")
    endif()

  endif()
  #
  ##############################################################################
  #                             Blas option                                    #
  ##############################################################################
  #
  if( SCALFMM_USE_BLAS )
    message(STATUS "CMAKE_CXX_COMPILER_ID STREQUAL  ${CMAKE_CXX_COMPILER_ID}")

    option( SCALFMM_USE_MKL_AS_BLAS "Set to ON to use MKL CBLAS" OFF )

    if( SCALFMM_USE_MKL_AS_BLAS )
      set(BLA_VENDOR "Intel10_64lp_seq")
      find_package(BLASEXT QUIET) # not REQUIRED

      if(BLAS_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${BLAS_LIBRARY_DIRS}")
      endif()
      unset(LAPACK_LIBRARIES)
      if (BLAS_LIBRARIES)
        set(BLASLAPACK_LIBRARIES ${BLAS_LIBRARIES})
      endif()
    elseif(SCALFMM_USE_EXTERNAL_BLAS)
      message(STATUS "BLAS SET BY EXTERNAL PROGRAM = ${BLAS_LIBRARIES}")
    else()
      find_package(BLASEXT)  # not REQUIRED
      find_package(LAPACKEXT)  # not REQUIRED
      set(BLASLAPACK_LIBRARIES)
      if (LAPACK_LIBRARIES)
        set(BLASLAPACK_LIBRARIES "${LAPACK_LIBRARIES}")
      endif()
      if (BLAS_LIBRARIES)
        list(APPEND BLASLAPACK_LIBRARIES "${BLAS_LIBRARIES}")
      endif()
    endif()

    if(BLAS_FOUND)
      set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES};${BLASLAPACK_LIBRARIES}")
      #message(STATUS "SCALFMM_LIBRARIES          = ${SCALFMM_LIBRARIES}")
    else()
      message(WARNING "BLAS has not been found, SCALFMM will continue to compile but some applications will be disabled.")
      message(WARNING "If you have BLAS set BLAS_LIBDIR, BLAS_INCDIR or BLAS_DIR (CMake variables using -D or environment variables).")        
      set(SCALFMM_USE_BLAS OFF)
    endif()

  endif(SCALFMM_USE_BLAS)
  list(APPEND FUSE_LIST "BLAS")
  #
  # FFT option
  #
  CMAKE_DEPENDENT_OPTION(SCALFMM_USE_MKL_AS_FFTW "Set to ON to use MKL FFTW" ON  "SCALFMM_USE_FFT;SCALFMM_USE_MKL_AS_BLAS" OFF )
  if( SCALFMM_USE_FFT )
    message(STATUS "SCALFMM USE FFT Configure:")

    if( SCALFMM_USE_MKL_AS_FFTW )

      message(STATUS "         SCALFMM USE MKL ")

      if( SCALFMM_USE_MKL_AS_BLAS )

        unset(FFT_LIBRARIES CACHE)
        message(STATUS "         SCALFMM USE MKL already defined")
        set(FFT_INCLUDES "$ENV{MKLROOT}/include/fftw" CACHE STRING "Set your MKL flags")
        if (BLAS_FOUND)
          set(FFTW_FOUND ON)
        endif()

      else(SCALFMM_USE_MKL_AS_BLAS)

        # The package can be used with the following COMPONENTS:
        # MKL, THREADS|OMP and/or SIMPLE|DOUBLE|LONG|QUAD
        # Default is DOUBLE and without THREADS|OMP
        find_package(FFTW  COMPONENTS MKL) # not REQUIRED
        if (FFTW_LIBRARY_DIRS_DEP)
          set(FFT_LIBRARIES "-L${FFTW_LIBRARY_DIRS_DEP};" CACHE STRING "Set your MKL flags")
        endif()
        if (FFTW_LIBRARIES_DEP)
          foreach (fft_lib ${FFTW_LIBRARIES_DEP})
            set(FFT_LIBRARIES "${FFT_LIBRARIES};${fft_lib};")
          endforeach()
        endif()
        set(FFT_INCLUDES "${FFTW_INCLUDE_DIRS_DEP}" )
        if (FFT_LIBRARIES)
          set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES};${FFT_LIBRARIES}")
        endif()

      endif(SCALFMM_USE_MKL_AS_BLAS)

    else(SCALFMM_USE_MKL_AS_FFTW)

      message(STATUS "         SCALFMM USE FFTW ")
      # The package can be used with the following COMPONENTS:
      # MKL, THREADS|OMP and/or SIMPLE|DOUBLE|LONG|QUAD
      # Default is DOUBLE and without THREADS|OMP
      find_package(FFTW COMPONENTS SIMPLE) # not REQUIRED
      if (FFTW_LIBRARY_DIRS_DEP)
        set(FFT_LIBRARIES "-L${FFTW_LIBRARY_DIRS_DEP};" CACHE STRING "Set your MKL flags")
      endif()
      if (FFTW_LIBRARIES_DEP)
        foreach (fft_lib ${FFTW_LIBRARIES_DEP})
          set(FFT_LIBRARIES "${FFT_LIBRARIES};${fft_lib};")
        endforeach()
      endif()
      set(FFT_INCLUDES "${FFTW_INCLUDE_DIRS_DEP}" )
      if (FFT_LIBRARIES)
        set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES};${FFT_LIBRARIES}")
      endif()

    endif(SCALFMM_USE_MKL_AS_FFTW)

    if (FFT_INCLUDES)
      set(SCALFMM_INCLUDES "${SCALFMM_INCLUDES}; ${FFT_INCLUDES}")
    endif()

    if(FFTW_FOUND)
      message(STATUS "         SCALFMM_LIBRARIES          = ${SCALFMM_LIBRARIES}")
      message(STATUS "         SCALFMM_INCLUDES          = ${SCALFMM_INCLUDES}")
    else()
      message(WARNING "FFTW has not been found, SCALFMM will continue to compile but some applications will be disabled.")
      message(WARNING "If you have FFTW set FFTW_LIBDIR, FFTW_INCDIR or FFTW_DIR (CMake variables using -D or environment variables).")        
      set(SCALFMM_USE_FFT OFF)
    endif()
  endif(SCALFMM_USE_FFT)
  list(APPEND FUSE_LIST "FFT")
  message(STATUS "         SCALFMM_LIBRARIES          = ${SCALFMM_LIBRARIES}")
  message(STATUS "         SCALFMM_INCLUDES          = ${SCALFMM_INCLUDES}")
  
  message(STATUS "SCALFMM_USE_FFT       =  ${SCALFMM_USE_FFT}")

  #
  # C++ 2011
  message( STATUS "CMAKE_SYSTEM:   ${CMAKE_SYSTEM}  CMAKE_CXX_COMPILER_ID ${CMAKE_CXX_COMPILER_ID}  CMAKE_CXX_COMPILER ${NAME}")
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE INTEL_VERSION)
    message( STATUS "Intel:   ${INTEL_VERSION}")
    set(SCALFMM_CXX_FLAGS "${SCALFMM_CXX_FLAGS} -Wall")
  else()
    # Compile option valid on GCC/G++
    set(SCALFMM_CXX_FLAGS "${SCALFMM_CXX_FLAGS} -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wconversion -Wcast-qual -Wcast-align -Woverloaded-virtual")
  endif()
  #
  ##################################################################
  #              STARPU Section
  ##################################################################
  #
  message( STATUS "SCALFMM_USE_STARPU             = ${SCALFMM_USE_STARPU}" )
  if( SCALFMM_USE_STARPU )
    # CUDA could be used with StarPU enabled
    option( SCALFMM_USE_CUDA "Set to ON to use CUDA with StarPU" OFF )
    message( STATUS "SCALFMM_USE_CUDA             = ${SCALFMM_USE_CUDA}" )
    if(SCALFMM_USE_CUDA)
      execute_process(COMMAND nvcc --version ERROR_VARIABLE cuda_error_output OUTPUT_QUIET)
      if(cuda_error_output)
        message( FATAL_ERROR "nvcc is needed with CUDA." )
      endif()
      if(NOT DEFINED CUSTOM_CUDA_FLAGS)
        set( CUSTOM_CUDA_FLAGS "-std=c++11;-arch=sm_20" CACHE
          STRING "Set your CUDA flags, for example : -arch=sm_20;-ptxas-options=-v;-use_fast_math")
      endif()
      # This is needed to remove backslash after space in ADD_CUSTOM_COMMAND
      separate_arguments(CUSTOM_CUDA_FLAGS)
      message( STATUS "CUSTOM_CUDA_FLAGS             = ${CUSTOM_CUDA_FLAGS}" )

      find_package(CUDA REQUIRED)

      if (CUDA_INCLUDE_DIRS)
        include_directories(${CUDA_INCLUDE_DIRS})
      endif()
      if (CUDA_LIBRARIES)
        set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES};${CUDA_LIBRARIES}")
      endif()
      
      set(CUDA_NEEDED_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/Src)
    endif()

    # Find StarPU with a list of optional components
    set(SCALFMM_STARPU_VERSION "1.1" CACHE STRING "oldest STARPU version desired")
    # create list of components in order to make a single call to find_package(starpu...)
    set(STARPU_COMPONENT_LIST "HWLOC")
    if(SCALFMM_USE_MPI)
      list(APPEND STARPU_COMPONENT_LIST "MPI")
    endif()
    if(SCALFMM_USE_CUDA)
      list(APPEND STARPU_COMPONENT_LIST "CUDA")
    endif()

    find_package(STARPU ${SCALFMM_STARPU_VERSION} REQUIRED
      COMPONENTS ${STARPU_COMPONENT_LIST})

    # Append list of libraries and include dirs
    include_directories(${STARPU_INCLUDE_DIRS_DEP})
    foreach (starpu_libdir ${STARPU_LIBRARY_DIRS_DEP})
      if (${starpu_libdir} MATCHES "^ *-L")
        set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES};${starpu_libdir}")
      else()
        set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES};-L${starpu_libdir}")
      endif()
    endforeach()
    foreach (starpu_lib ${STARPU_LIBRARIES_DEP})
      if (EXISTS ${starpu_lib} OR ${starpu_lib} MATCHES "^ *-")
        set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES};${starpu_lib}")
      else()
        set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES};-l${starpu_lib}")
      endif()
    endforeach()

    # TODO: is this very useful? CUDA is already a component of find starpu
    if (CUDA_LIBRARIES)
      set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES};${CUDA_LIBRARIES}")
    endif()

    # Message
    message(STATUS "         STARPU_LIBRARIES          = ${STARPU_LIBRARIES}")
    if (STARPU_INCLUDE_DIRS)
      message(STATUS "         STARPU_INCLUDES           = ${STARPU_INCLUDE_DIRS}")
      set(SCALFMM_INCLUDES  "${SCALFMM_INCLUDES}; ${STARPU_INCLUDE_DIRS}")
    endif()

    # TODO: replace this by a component of find starpu
    OPTION( SCALFMM_USE_OPENCL      "Set to ON to use OPENCL with StarPU" OFF )
    MESSAGE( STATUS "SCALFMM_USE_OPENCL             = ${SCALFMM_USE_OPENCL}" )
    if(SCALFMM_USE_OPENCL)
      include_directories($ENV{OPENCL_INC})
      SET(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES};-L$ENV{OPENCL_LIB};-lOpenCL")
    endif()
  endif(SCALFMM_USE_STARPU)
  list(APPEND FUSE_LIST "STARPU")
  list(APPEND FUSE_LIST "CUDA")
  list(APPEND FUSE_LIST "OPENCL")

  ##################################################################
  #                         Use SSE                                #
  ##################################################################

  message( STATUS "SCALFMM_USE_SSE              = ${SCALFMM_USE_SSE}" )
  if( SCALFMM_USE_SSE )
    if(NOT EXISTS ${SCALFMM_CMAKE_MODULE_PATH}/compileTestSse.cpp)
      message(FATAL_ERROR "The CompileTestSseFile does not exist (${SCALFMM_CMAKE_MODULE_PATH}/compileTestSse.cpp)" )
    endif()
    message( STATUS "SSE_FLAGS ${SSE_FLAGS}  -- ${CMAKE_CXX_FLAGS}  ")
    try_compile(COMPILE_SSE  ${CMAKE_CURRENT_BINARY_DIR}
      ${SCALFMM_CMAKE_MODULE_PATH}/compileTestSse.cpp
      COMPILE_DEFINITIONS "${CMAKE_CXX_FLAGS} ${SSE_FLAGS}"
      OUTPUT_VARIABLE COMPILE_SSE_OUTPUT)

    if(${COMPILE_SSE})
      set(SCALFMM_CXX_FLAGS "${SCALFMM_CXX_FLAGS} ${SSE_FLAGS}")

      try_compile(COMPILE_RESULT_VAR ${CMAKE_CURRENT_BINARY_DIR}
        ${SCALFMM_CMAKE_MODULE_PATH}/checkSSEpe.cpp
        COMPILE_DEFINITIONS "${CMAKE_CXX_FLAGS} ${SSE_FLAGS}")
      if( NOT ${COMPILE_RESULT_VAR})
        set(__SSEPE_INTEL_COMPILER ON)
      endif()
      #set(SCALFMM_USE_AVX OFF)
    else(${COMPILE_SSE})
      message(FATAL_ERROR "SSE NOT SUPPORTED ; Set SCALFMM_USE_SSE  to OFF \n Output from test is : ${COMPILE_SSE_OUTPUT}")
    endif(${COMPILE_SSE})
  endif()
  list(APPEND FUSE_LIST "SSE")

  ##################################################################
  #                           Use AVX                              #
  ##################################################################

  message(STATUS "SCALFMM_USE_AVX               = ${SCALFMM_USE_AVX}")
  if(SCALFMM_USE_AVX)
    if(NOT EXISTS ${SCALFMM_CMAKE_MODULE_PATH}/compileTestAvx.cpp)
      message(FATAL_ERROR "The CompileTestSseFile does not exist (${SCALFMM_CMAKE_MODULE_PATH}/compileTestAvx.cpp)" )
    endif()

    try_compile(COMPILE_AVX ${CMAKE_CURRENT_BINARY_DIR}
      ${SCALFMM_CMAKE_MODULE_PATH}/compileTestAvx.cpp
      COMPILE_DEFINITIONS "${CMAKE_CXX_FLAGS} ${AVX_FLAGS}"
      OUTPUT_VARIABLE COMPILE_AVX_OUTPUT)
    if(${COMPILE_AVX})
      message(STATUS "%%%%%%%%%%%% COMPILE_AVX               = ${COMPILE_AVX}  %%%%<    ${AVX_FLAGS}")
      
      set(SCALFMM_CXX_FLAGS "${SCALFMM_CXX_FLAGS}   ${AVX_FLAGS}")
      message(STATUS "%%%%%%%%%%%% SCALFMM_CXX_FLAGS               = ${SCALFMM_CXX_FLAGS}")
      #set( SCALFMM_USE_SSE   OFF   FORCE) # ne marche pas
      try_compile(COMPILE_RESULT_AVSPE ${CMAKE_CURRENT_BINARY_DIR}
        ${SCALFMM_CMAKE_MODULE_PATH}/checkAVXpe.cpp
        COMPILE_DEFINITIONS "${CMAKE_CXX_FLAGS} ${AVX_FLAGS}")
      if( NOT ${COMPILE_RESULT_AVSPE})


        set(__AVXPE_INTEL_COMPILER ON)
      endif()

      message(STATUS ${CMAKE_CXX_FLAGS} )
    else(${COMPILE_AVX})
      message(FATAL_ERROR "AVX NOT SUPPORTED ; Set SCALFMM_USE_AVX  to OFF \n Output from test is : ${COMPILE_AVX_OUTPUT} ")
    endif(${COMPILE_AVX})
  endif(SCALFMM_USE_AVX)
  list(APPEND FUSE_LIST "AVX")
  #
  # Error if both SCALFMM_USE_AVX AND SCALFMM_USE_SSE are set
  #
  if( SCALFMM_USE_AVX AND SCALFMM_USE_SSE)
    message(FATAL_ERROR "Check SCALFMM_USE_SSE or SCALFMM_USE_AVX BUT NOT BOTH. ")
  endif(SCALFMM_USE_AVX AND SCALFMM_USE_SSE)
  ##################################################################
  #                           Use AVX2                             #
  ##################################################################

  message(STATUS "SCALFMM_USE_AVX2               = ${SCALFMM_USE_AVX2}")
  if(SCALFMM_USE_AVX2)
    if(NOT EXISTS ${SCALFMM_CMAKE_MODULE_PATH}/compileTestAvx2.cpp)
      message(FATAL_ERROR "The CompileTestSseFile does not exist (${SCALFMM_CMAKE_MODULE_PATH}/compileTestAvx.cpp)" )
    endif()

    try_compile(COMPILE_AVX2 ${CMAKE_CURRENT_BINARY_DIR}
      ${SCALFMM_CMAKE_MODULE_PATH}/compileTestAvx2.cpp
      COMPILE_DEFINITIONS "${CMAKE_CXX_FLAGS} ${AVX2_FLAGS}"
      OUTPUT_VARIABLE COMPILE_AVX2_OUTPUT)
    if(${COMPILE_AVX2})
      set(SCALFMM_CXX_FLAGS "${SCALFMM_CXX_FLAGS}   ${AVX2_FLAGS}")
      #set( SCALFMM_USE_SSE   OFF   FORCE) # ne marche pas
      try_compile(COMPILE_RESULT_AVSPE ${CMAKE_CURRENT_BINARY_DIR}
        ${SCALFMM_CMAKE_MODULE_PATH}/checkAVX2pe.cpp
        COMPILE_DEFINITIONS "${CMAKE_CXX_FLAGS} ${AVX2_FLAGS}")
      if( NOT ${COMPILE_RESULT_AVSPE})
        set(__AVX2PE_INTEL_COMPILER ON)
      endif()

      message(STATUS ${CMAKE_CXX_FLAGS} )
    else(${COMPILE_AVX2})
      message(FATAL_ERROR "AVX2 NOT SUPPORTED ; Set SCALFMM_USE_AVX2  to OFF \n Output from test is : ${COMPILE_AVX_OUTPUT} ")
    endif(${COMPILE_AVX2})
  endif(SCALFMM_USE_AVX2)
  list(APPEND FUSE_LIST "AVX2")
  #
  # Error if both SCALFMM_USE_AVX2 AND SCALFMM_USE_SSE are set
  #
  if( SCALFMM_USE_AVX2 AND SCALFMM_USE_SSE)
    message(FATAL_ERROR "Check SCALFMM_USE_SSE or SCALFMM_USE_AVX2 BUT NOT BOTH. ")
  endif(SCALFMM_USE_AVX2 AND SCALFMM_USE_SSE)
  ##################################################################
  #                     Use  native MIC compilation                #
  ##################################################################
  if( SCALFMM_USE_MIC_NATIVE )
    include(CheckCCompilerFlag)
    check_c_compiler_flag(-mmic HAS_MMIC_FLAG)
    if(NOT HAS_MMIC_FLAG)
      message(FATAL_ERROR "MIC NOT SUPPORTED ; Set SCALFMM_USE_MIC_NATIVE to OFF")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmic")
  else()
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xhost")
  endif()
  list(APPEND FUSE_LIST "MIC")
  ##################################################################
  #
  #   Set EZTRACE
  #
  ##################################################################
  If(SCALFMM_USE_EZTRACE)
    find_package (PkgConfig)
    if(PKG_CONFIG_FOUND)
      set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH "ON")
      pkg_search_module( EZTrace REQUIRED eztrace)
      if(EZTrace_FOUND)
        link_directories(${EZTrace_LIBRARY_DIRS})
        link_libraries( ${EZTrace_LIBRARIES})
        IF( SCALFMM_USE_MPI )
          link_libraries(-leztrace-autostart-mpi)
  ENDIF(SCALFMM_USE_MPI)     
  include_directories(${EZTrace_INCLUDE_DIRS})
  MESSAGE(STATUS "EZTRACE:  ${EZTrace_INCLUDE_DIRS}   ${EZTrace_LIBRARY_DIRS}  ${EZTrace_LIBRARIES}")
  CMAKE_DEPENDENT_OPTION(SCALFMM_TRACE_ALGO "Set to ON to trace the full algorithm (all operators)" ON  "SCALFMM_USE_EZTRACE" OFF )
  CMAKE_DEPENDENT_OPTION(SCALFMM_TRACE_P2M   "Set to ON to trace P2M operator" OFF  "SCALFMM_USE_EZTRACE" OFF )
  CMAKE_DEPENDENT_OPTION(SCALFMM_TRACE_M2M   "Set to ON to trace M2M operator" OFF  "SCALFMM_USE_EZTRACE" OFF )
  CMAKE_DEPENDENT_OPTION(SCALFMM_TRACE_M2L   "Set to ON to trace M2L operator" OFF  "SCALFMM_USE_EZTRACE" OFF )
  CMAKE_DEPENDENT_OPTION(SCALFMM_TRACE_L2L   "Set to ON to trace L2L operator" OFF  "SCALFMM_USE_EZTRACE" OFF )
  CMAKE_DEPENDENT_OPTION(SCALFMM_TRACE_P2P   "Set to ON to trace P2P operator" OFF  "SCALFMM_USE_EZTRACE" OFF )
      else(EZTrace_FOUND)
  MESSAGE(WARNING "Eztrace not found - EZTRACE Is set to OFF")
  set(SCALFMM_USE_EZTRACE OFF)  
      endif(EZTrace_FOUND)
    else(PKG_CONFIG_FOUND)
      MESSAGE(WARNING "PKG-CONFIG not found - EZTRACE Is set to OFF")
      set(SCALFMM_USE_EZTRACE OFF)
    endif(PKG_CONFIG_FOUND)

  endif(SCALFMM_USE_EZTRACE)
  
  ##################################################################
  #
  #   To catch signals
  #
  ##################################################################  
    if((WIN32 OR (${CMAKE_SYSTEM_NAME} MATCHES "Windows")) AND SCALFMM_USE_SIGNALS)
      # We are on Windows and we cannot catch signals
      message( FATAL_ERROR "Catching the signals requieres an UNIX system." )
    else()
        IF( NOT APPLE)
        SET(SCALFMM_CXX_FLAGS "${SCALFMM_CXX_FLAGS} -rdynamic")
        ENDIF()
        
    endif()
  ##################################################################
  #                                                                #
  #                     END SETTING VARIABLES                      #
  #                                                                #
  ##################################################################
  #
  # Generate  ScalFmmConfig.h
  #
  set(SCALFMM_COMPILE_FLAGS "${SCALFMM_CXX_FLAGS} ${CMAKE_CXX_FLAGS_TYPE} ${CMAKE_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS "${SCALFMM_COMPILE_FLAGS}")
  # We need the libraries without spaces (inside the config file)
  set(SCALFMM_COMPILE_LIBS   "")
  foreach(lib_var ${SCALFMM_LIBRARIES})
    string(STRIP ${lib_var} lib_var)
    LIST(APPEND SCALFMM_COMPILE_LIBS   ${lib_var})        
  endforeach()
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/Src/ScalFmmConfig.h.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/Src/ScalFmmConfig.h  )
  #
  ##################################################################
  #                      Build - lib                               #
  ##################################################################
  #
  add_subdirectory(Src)
  # Link with scalfmm lib
  set(scalfmm_lib scalfmm)
  #
  ##################################################################
  #            Remove compilation  of the drivers                  #
  ##################################################################
  if( SCALFMM_BUILD_EXAMPLES )
    # Build - Examples and drivers
    add_subdirectory(Examples)
  endif()
  ##################################################################
  #           Build tools (Internal use)                           #
  ##################################################################
  # Build - UTILs
  if( SCALFMM_ONLY_DEVEL AND SCALFMM_BUILD_UTILS)
    add_subdirectory(Utils)
  endif()

  ##################################################################
  #                Build -Tests to see examples                    #
  ##################################################################
  # Build - Tests
  message( STATUS "SCALFMM_BUILD_TESTS          = ${SCALFMM_BUILD_TESTS}" )
  if( SCALFMM_BUILD_TESTS )
    add_subdirectory(Tests)
  endif()

  ##################################################################
  #              Build - Unitary and numerical Tests               #
  ##################################################################
  #
  message( STATUS "SCALFMM_BUILD_UTESTS         = ${SCALFMM_BUILD_UTESTS}" )
  if( SCALFMM_BUILD_UTESTS )
    include(CTest)
    enable_testing()
    add_subdirectory(UTests)
  endif()

  ##################################################################
  #                           Add - doc                            #
  ##################################################################
  if(SCALFMM_BUILD_DOC)
    add_subdirectory(Doc)
  endif()
  ##################################################################
  #                       Build - Addons                           #
  ##################################################################
  #
  if( SCALFMM_USE_ADDONS )
    file(GLOB all_addons ./Addons/*)
    # Options
    option( SCALFMM_ADDON_CKERNELAPI "Set to ON to build ScaFMM FMM API interface" OFF )

    foreach(addon_name ${all_addons})
      message(STATUS "Addons += ${addon_name}" )
      add_subdirectory(${addon_name})
    endforeach()
  endif()

  # Add custom clean command
  if("${CMAKE_GENERATOR}" MATCHES "Make")
    add_custom_target(reset
      COMMAND ${CMAKE_MAKE_PROGRAM} clean && rm -r ./CMake* && cmake ..
      WORKING_DIRECTORY ${CMAKE_CURRENT_DIR}
      COMMENT "Cleaning Build dir and CMake cache, then do cmake.."
      )
  endif()
  #
  ##################################################################
  #                                                                #
  #                  Build and export cmake files                  #
  #                                                                #
  ##################################################################
  #
  message(STATUS "SCALFMM_BINARY_DIR:          " ${SCALFMM_BINARY_DIR})
  configure_file(${SCALFMM_SOURCE_DIR}/ScalFMMConfig.cmake.in
    ${SCALFMM_BINARY_DIR}/ScalFMMConfig.cmake
    @ONLY IMMEDIATE
    )
  #
  install(FILES  ${SCALFMM_BINARY_DIR}/ScalFMMConfig.cmake  DESTINATION lib/  )
  if( SCALFMM_INSTALL_DATA )
    install(FILES  ${SCALFMM_SOURCE_DIR}/Data/test20k.fma  DESTINATION Data/   )
  endif()

  ##################################################################
  #                                                                #
  #              build a CPack driven installer package            #
  #                                                                #
  ##################################################################
  #
  # Use git commit number as CPACK_PACKAGE_VERSION_PATCH
  set(CPACK_PACKAGE_VERSION_PATCH 0)
  execute_process(
    COMMAND git rev-list HEAD --count
    OUTPUT_VARIABLE CPACK_PACKAGE_VERSION_PATCH
    RESULT_VARIABLE RET
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  #string_TRIM(PATCH1 PATCH)
  set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY "ON")
  #
  set (CPACK_RESOURCE_FILE_LICENSE  "${CMAKE_CURRENT_SOURCE_DIR}/Licence.txt")
  set(CPACK_PACKAGE_VERSION_MAJOR "${SCALFMM_MAJOR_VERSION}")
  set(CPACK_PACKAGE_VERSION_MINOR "${SCALFMM_MINOR_VERSION}")
  #
  set(PACK_PACKAGE_VERSION "${SCALFMM_MAJOR_VERSION}.${SCALFMM_MINOR_VERSION}-${CPACK_PACKAGE_VERSION_PATCH}")
  set(CPACK_SOURCE_GENERATOR "TGZ")
  set(CPACK_SOURCE_PACKAGE_FILE_NAME  "SCALFMM-${SCALFMM_MAJOR_VERSION}.${SCALFMM_MINOR_VERSION}-${CPACK_PACKAGE_VERSION_PATCH}")
  set(CPACK_SOURCE_IGNORE_FILES "\\\\.git;.DS_Store;.*~;/*.aux;/*.idx;/*.log;/*.out;/*.toc;/*.ilg;scalfmm.pro*;org.eclipse.core.resources.prefs;.cproject;.project")
  list(APPEND CPACK_SOURCE_IGNORE_FILES  "${CMAKE_CURRENT_BINARY_DIR};${CMAKE_CURRENT_SOURCE_DIR}/Utils/;Notes;Deprecated;/Build*;/noDist/;/Bench;ToRemove;Olivier;Addons/HMat")
  # Uniform;GroupTree;Adaptive;testUnif*;/*Lagrange*")
  #list(APPEND CPACK_SOURCE_IGNORE_FILES  "Stages;Uniform;Adaptive;testUnif*;*Lagrange*" )
  #
  include(CPack)
  #
  #
  ##################################################################
  #                                                                #
  #                      PRINT messages                            #
  #                                                                #
  ##################################################################
  message(STATUS "CPACK_SOURCE_IGNORE_FILES = ${CPACK_SOURCE_IGNORE_FILES}")
  message(STATUS "CPACK_SOURCE_PACKAGE_FILE_NAME = ${CPACK_SOURCE_PACKAGE_FILE_NAME}")
  #
  #                    Use Mem stats
  message(STATUS "SCALFMM_USE_MEM_STATS        = ${SCALFMM_USE_MEM_STATS}" )
  #                    Use Log
  message(STATUS "SCALFMM_USE_LOG              = ${SCALFMM_USE_LOG}" )
  #                    Use Assert
  message(STATUS "SCALFMM_USE_ASSERT           = ${SCALFMM_USE_ASSERT}" )
  #
  #message(STATUS "CMAKE_CXX_FLAGS   =  ${CMAKE_CXX_FLAGS}")
  #set(CMAKE_CXX_FLAGS "${SCALFMM_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
  message(STATUS "SCALFMM_USE_BLAS     =  ${SCALFMM_USE_BLAS}")
  message(STATUS "SCALFMM_USE_FFT       =  ${SCALFMM_USE_FFT}")
  message(STATUS "SCALFMM_USE_MKL      =   ${SCALFMM_USE_MKL}")
  # 
  message(STATUS "CMAKE_CXX_FLAGS     =  ${CMAKE_CXX_FLAGS}")
  message(STATUS "SCALFMM_CXX_FLAGS   =  ${SCALFMM_CXX_FLAGS}")
  message(STATUS "SCALFMM_LIBRARIES   =  ${SCALFMM_LIBRARIES}")
  message(STATUS "SCALFMM_INCLUDES    =  ${SCALFMM_INCLUDES}")
  ##################################################################
  #                               END                              #
  ##################################################################

else(MORSE_DISTRIB_DIR OR EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/morse/")

  message(STATUS "MORSE_DISTRIB_DIR is not set")
  message(STATUS "Please indicate where is located your MORSE distribution directory."
    " This is necessary to find cmake_modules.")
  message(STATUS "If you are a MORSE developer, you normally have access to the svn repository of the MORSE_DISTRIB_DIR: "
    "svn checkout svn+ssh://yourlogin@scm.gforge.inria.fr/svnroot/morse/trunk/morse_distrib")

endif(MORSE_DISTRIB_DIR OR EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/morse/")
