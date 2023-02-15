cmake_minimum_required (VERSION 3.10.0 FATAL_ERROR)

#
# Project Declaration
#--------------------
project(scalfmm CXX)
# check if compiling into source directories
string(COMPARE EQUAL "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}" insource)
if(insource)
  message(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto ./Build and tapes cmake ../")
endif(insource)

#
# Version Number
# --------------

# SCALFMM version number.  An even minor number corresponds to releases.
set(${CMAKE_PROJECT_NAME}_MAJOR_VERSION 2)
set(${CMAKE_PROJECT_NAME}_MINOR_VERSION 0)
set(${CMAKE_PROJECT_NAME}_PATCH_VERSION rc0)
set(${CMAKE_PROJECT_NAME}_VERSION  "${${CMAKE_PROJECT_NAME}_MAJOR_VERSION}.${${CMAKE_PROJECT_NAME}_MINOR_VERSION}.${${CMAKE_PROJECT_NAME}_PATCH_VERSION}" )

#
# Set global c++ support
# ----------------------
set(CMAKE_CXX_STANDARD 14)

#
# Update cmake path
# -----------------
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/)

#
# Creating main lib
# -----------------
list(APPEND SCALFMM_SOURCES FAlgorithmTimers.cpp FLog.cpp FMemStats.cpp FSignalHandler.cpp)
list(TRANSFORM SCALFMM_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/")
add_library(${CMAKE_PROJECT_NAME} STATIC ${SCALFMM_SOURCES})
# Set library includes
# --------------------
target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${${CMAKE_PROJECT_NAME}_BINARY_DIR}/include>
  $<BUILD_INTERFACE:${${CMAKE_PROJECT_NAME}_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  )

if(MORSE_DISTRIB_DIR OR EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/morse_cmake")

  # Add extra cmake module path and initialize morse cmake modules
  # --------------------------------------------------------------
  if(MORSE_DISTRIB_DIR)
    list(APPEND CMAKE_MODULE_PATH ${MORSE_DISTRIB_DIR}/modules)
    set(MORSE_CMAKE_MODULE_PATH ${MORSE_DISTRIB_DIR}/modules )
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/morse_cmake")
    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/morse_cmake/modules)
    set(MORSE_CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/morse_cmake/modules )
  endif()

  include(MorseInit)

  # Set the RPATH config
  # These options allow all installed executables to run correctly
  # when installed to the install diretory
  # --------------------
  # use, i.e. don't skip the full RPATH for the build tree
  set(CMAKE_SKIP_BUILD_RPATH  FALSE)
  # when building, don't use the install RPATH already
  # (but later on when installing)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
  # the RPATH to be used when installing
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


  # Build ScalFmm as static or dynamic lib
  # --------------------------------------
  # option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
  #
  # Options
  #
  # Cmake options for dependencies
  option( SCALFMM_USE_MPI              "Set to ON to build ScaFMM with MPI"         OFF )
  option( SCALFMM_USE_BLAS             "Set to ON to build ScaFMM with BLAS"        ON  )
  option( SCALFMM_USE_FFT              "Set to ON to build ScaFMM with FFTW"        OFF )
  option( SCALFMM_USE_STARPU           "Set to ON to build SCALFMM with StarPU"     OFF )

  # Cmake options that impact the build tree
  option( SCALFMM_BUILD_EXAMPLES       "Set to ON to compile examples "                ON  )
  option( SCALFMM_BUILD_TESTS          "Set to ON to build functionalities Tests"      OFF )
  option( SCALFMM_BUILD_UTESTS         "Set to ON to build UTests"                     OFF )
  option( SCALFMM_BUILD_UTILS          "Set to ON to build utils Tests"                OFF )
  option( SCALFMM_BUILD_DOC            "Set to ON to build the Doxygen documentation " OFF )

  # Cmake options to install additional content
  option( SCALFMM_INSTALL_DATA         "Set to ON to include the data when proceeding to the installtion " OFF )

  # Cmake options related to trace, logging and statistics
  option( SCALFMM_USE_LOG              "Set to ON to print output debug information" OFF )
  option( SCALFMM_USE_MEM_STATS        "Set to ON to profile memory"                 OFF )
  option( SCALFMM_USE_ASSERT           "Set to ON to enable safe tests during execution" ON  )
  option( SCALFMM_ONLY_DEVEL           "Set to ON to compile Development tools (only scalfmm team)" OFF )

  # Cmake options related to the targeted architecture
  option( SCALFMM_USE_MIC_NATIVE       "Set to ON to compile in native mode for MIC" OFF  )



  #
  # Set scalfmm global lists
  # ------------------------
  list(APPEND SCALFMM_CXX_FLAGS -fpic -Wall)

  # TODO: Find new way to resolve dependencies...
  # Set FUSE LIST
  # -------------
  set(FUSE_LIST "") # Construct during configure
  set(FUSE_DEP_AVAILABLE "MPI;CUDA;STARPU;OMP4;BLAS;FFT;OPENCL") # List to compare against when compiling tests, etc...
  #
  # OpenMP
  # ------
  if(NOT OpenMP_CXX_FOUND)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
      target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
      list(APPEND SCALFMM_LIBRARIES OpenMP)
      if(OpenMP_CXX_VERSION_MAJOR GREATER_EQUAL 4)
        list(APPEND FUSE_LIST OMP4)
      endif()
    else(OpenMP_CXX_FOUND)
      message(WARNING "OPENMP NOT FOUND")
    endif(OpenMP_CXX_FOUND)
  endif(NOT OpenMP_CXX_FOUND)

  #
  # Inastemp
  # --------
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/inastemp")

    # Set option() CMP0077 to old for inastemp...
    set(CMAKE_POLICY_DEFAULT_CMP0077 OLD)
    set(INASTEMP_JUST_LIB TRUE)
    set(INASTEMP_DO_INSTALL TRUE)
    # Create target
    add_subdirectory(inastemp)
    add_library(inastemp INTERFACE)
    target_include_directories(inastemp INTERFACE
      $<BUILD_INTERFACE:${INASTEMP_SOURCE_DIR}/Src>
      $<BUILD_INTERFACE:${INASTEMP_BINARY_DIR}/Src>
      $<INSTALL_INTERFACE:include>
      )
    string(REPLACE " -" ";-" INASTEMP_FLAGS "${INASTEMP_CXX_FLAGS}")
    list(REMOVE_DUPLICATES INASTEMP_FLAGS)
    target_compile_options(inastemp INTERFACE ${INASTEMP_FLAGS})
    target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC inastemp)

    message(STATUS "Inastemp flags : ${INASTEMP_FLAGS}")
    message(WARNING "Inastemp configures compile flags according to the current architure "
            "you are compiling on. If you target a different machine, be carefull "
            "with the available SIMD instruction set supported."
            "Thus, use INASTEMP_USE_XXX to target the correct instruction set.")

  else(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/inastemp")

    message(STATUS "Inastemp is not available.")
    message(STATUS "Please, initialize all git submodule.")
    message(FATAL_ERROR "Exit : Inastemp dependency not found.")

  endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/inastemp")

  #
  # Specific Debug flags
  # --------------------
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    list(APPEND SCALFMM_CXX_FLAGS -funroll-loops)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
      list(APPEND SCALFMM_CXX_FLAGS -fp-model\ strict)
    endif()
  endif()

  #
  # Specific Release flags
  # ----------------------
  if(CMAKE_BUILD_TYPE STREQUAL "Release")
    list(APPEND SCALFMM_CXX_FLAGS -funroll-loops)
    if(NOT APPLE)
      list(APPEND -ftree-vectorize)
      if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        list(APPEND SCALFMM_CXX_FLAGS  -ipo  -fstrict-aliasing)
      endif()
      list(APPEND SCALFMM_CXX_FLAGS -ffast-math -ftree-vectorize)
    endif(NOT APPLE)
  endif()

  #
  # MPI
  # ---
  if(SCALFMM_USE_MPI)
    if(NOT MPI_FOUND)
      try_compile(COMPILE_INTEL ${CMAKE_CURRENT_BINARY_DIR}
      ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/compileTestIntel.cpp
      COMPILE_DEFINITIONS "${CMAKE_CXX_FLAGS}")
      if(COMPILE_INTEL)
        set(CMAKE_CXX_COMPILER_ID  "Intel")
      endif()
      message(STATUS "CXX  ${CMAKE_CXX_COMPILER_ID}" )

      find_package(MPI REQUIRED)
      if(MPI_CXX_FOUND)
        target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC MPI::MPI_CXX)
        list(APPEND SCALFMM_LIBRARIES MPI)
        list(APPEND FUSE_LIST MPI)
	target_sources(${CMAKE_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/FMpi.cpp)
      else(MPI_CXX_FOUND)
        message(FATAL_ERROR "MPI_CXX is required but was not found. "
                "Please provide an MPI compiler in your environment."
                "Note that you can give the path to your MPI installation "
                "by setting MPI_DIR cmake variable.")
      endif(MPI_CXX_FOUND)
    endif(NOT MPI_FOUND)
  endif(SCALFMM_USE_MPI)
  #
  # Blas options
  # ------------
  if(SCALFMM_USE_BLAS)

    option( SCALFMM_USE_MKL_AS_BLAS   "Set to ON to use MKL BLAS/LAPACK" OFF )
    option( SCALFMM_USE_ESSL_AS_BLAS  "Set to ON to use ESSL BLAS/LAPACK" OFF )
    option( SCALFMM_USE_EXTERNAL_BLAS "Set to ON to use external BLAS. Please, set BLAS_LIBRARIES variable." OFF )

    if( SCALFMM_USE_MKL_AS_BLAS )
      set(BLA_VENDOR "Intel10_64lp_seq")
      find_package(BLAS ) # not REQUIRED
    elseif(SCALFMM_USE_ESSL_AS_BLAS)
      set(BLA_VENDOR "IBMESSL")
      find_package(BLAS   QUIET) # not REQUIRED
      find_package(LAPACK QUIET) # not REQUIRED
    elseif(SCALFMM_USE_EXTERNAL_BLAS)
      message(STATUS "BLAS SET BY EXTERNAL PROGRAM = ${BLAS_LIBRARIES}")
      list(APPEND BLASLAPACK_LIBRARIES "${BLAS_LIBRARIES}")
    else()
      find_package(BLAS)  # not REQUIRED
      find_package(LAPACK)  # not REQUIRED
    endif()

    if(BLAS_FOUND)
      target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC blas::blas)
      list(APPEND SCALFMM_LIBRARIES BLAS)
      list(APPEND FUSE_LIST BLAS)

      if(LAPACK_FOUND)
        target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC lapack::lapack)
        list(APPEND SCALFMM_LIBRARIES LAPACK)
      endif()

      message (STATUS "Check BLAS Fortran mangling")
      # add options to let the user be able to force a behavior
      option( SCALFMM_BLAS_ADD_     "Set to ON to force calls to BLAS Fortran symbols with _ (ex: dgemm_)"        ON )
      option( SCALFMM_BLAS_UPCASE   "Set to ON to force calls to BLAS Fortran symbols in capital (ex: DGEMM)"     OFF )
      option( SCALFMM_BLAS_NOCHANGE "Set to ON to force calls to BLAS Fortran symbols with no change (ex: dgemm)" OFF )
      # if options not changed by user then auto-detection
      if (NOT SCALFMM_BLAS_ADD_ AND NOT SCALFMM_BLAS_UPCASE AND NOT SCALFMM_BLAS_NOCHANGE)
        # give blas libraries and check dgemm symbol
        set(CMAKE_REQUIRED_LIBRARIES "${BLAS_LIBRARIES}")
        check_function_exists(dgemv_ DGEMV_ADD_)
        check_function_exists(DGEMV  DGEMV_UPCASE)
        check_function_exists(dgemv  DGEMV_NOCHANGE)
        # here we consider that the first kind of symbol found will be the one used
        # current order is: ADD_, UPCASE, NOCHANGE
        if (DGEMV_ADD_)
          set (SCALFMM_BLAS_ADD_     ON)
          set (SCALFMM_BLAS_UPCASE   OFF)
          set (SCALFMM_BLAS_NOCHANGE OFF)
          message (STATUS "BLAS dgemv_ symbol found, SCALFMM_BLAS_ADD_ is ON")
        else (DGEMV_ADD_)
          if (DGEMV_UPCASE)
            set (SCALFMM_BLAS_ADD_     OFF)
            set (SCALFMM_BLAS_UPCASE   ON)
            set (SCALFMM_BLAS_NOCHANGE OFF)
            message (STATUS "BLAS DGEMV symbol found, SCALFMM_BLAS_UPCASE is ON")
          else (DGEMV_UPCASE)
            if (DGEMV_NOCHANGE)
              set (SCALFMM_BLAS_ADD_     OFF)
              set (SCALFMM_BLAS_UPCASE   OFF)
              set (SCALFMM_BLAS_NOCHANGE ON)
              message (STATUS "BLAS dgemv symbol found, SCALFMM_BLAS_NOCHANGE is ON")
            endif (DGEMV_NOCHANGE)
          endif (DGEMV_UPCASE)
        endif (DGEMV_ADD_)
        if ( (NOT DGEMV_ADD_) AND (NOT DGEMV_UPCASE) AND (NOT DGEMV_NOCHANGE) )
          message(FATAL_ERROR "BLAS Fortran mangling is not properly detected - please check your BLAS libraries")
        endif ()
      else  (NOT SCALFMM_BLAS_ADD_ AND NOT SCALFMM_BLAS_UPCASE AND NOT SCALFMM_BLAS_NOCHANGE)
        if (SCALFMM_BLAS_ADD_)
          message (STATUS "SCALFMM_BLAS_ADD_ is set to ON")
        endif()
        if (SCALFMM_BLAS_UPCASE)
          message (STATUS "SCALFMM_BLAS_UPCASE is set to ON")
        endif()
        if (SCALFMM_BLAS_NOCHANGE)
          message (STATUS "SCALFMM_BLAS_NOCHANGE is set to ON")
        endif()
      endif (NOT SCALFMM_BLAS_ADD_ AND NOT SCALFMM_BLAS_UPCASE AND NOT SCALFMM_BLAS_NOCHANGE)
    else()
      MESSAGE(WARNING "BLAS: ${BLAS_FOUND}")

      message(WARNING "BLAS has not been found, SCALFMM will continue to compile but some applications will be disabled.")
      message(WARNING "If you have BLAS set BLAS_LIBDIR, BLAS_INCDIR or BLAS_DIR (CMake variables using -D or environment variables).")
      set(SCALFMM_USE_BLAS OFF)
    endif()

  endif(SCALFMM_USE_BLAS)

  #
  # FFT options
  # -----------
  CMAKE_DEPENDENT_OPTION(SCALFMM_USE_MKL_AS_FFTW  "Set to ON to use MKL FFTW"  ON  "SCALFMM_USE_FFT;SCALFMM_USE_MKL_AS_BLAS"  OFF )
  CMAKE_DEPENDENT_OPTION(SCALFMM_USE_ESSL_AS_FFTW "Set to ON to use ESSL FFTW" ON  "SCALFMM_USE_FFT;SCALFMM_USE_ESSL_AS_BLAS" OFF )

  if(SCALFMM_USE_FFT)
    message(STATUS "SCALFMM USE FFT Configure:")

    # The package find_package(FFTW) can be used with the following COMPONENTS:
    # MKL, ESSL, THREADS|OMP and/or SIMPLE|LONG|QUAD
    # Default will find the real double precision fftw library version without THREADS|OMP
    if(SCALFMM_USE_MKL_AS_FFTW)
      message(STATUS "         SCALFMM USE FFT from MKL")
      set(BLA_VENDOR "Intel10_64lp_seq")
      if(NOT TARGET blas::blas)
        find_package(FFTW COMPONENTS MKL)
      else()
        set(FFTW_FOUND 1)
        find_package(FFTW COMPONENTS MKL QUIET)
        target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC ${FFTW_INCLUDE_DIRS})
      endif()
    elseif (SCALFMM_USE_ESSL_AS_FFTW)
      message(STATUS "         SCALFMM USE FFT from ESSL ")
      find_package(FFTW COMPONENTS ESSL)
    else()
      message(STATUS "         SCALFMM USE FFTW")
      find_package(FFTW COMPONENTS SIMPLE)
    endif()

    if(FFTW_FOUND)
      list(APPEND FUSE_LIST FFT)
      if(NOT SCALFMM_USE_MKL_AS_FFTW)
        target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC fftw::fftw)
        list(APPEND SCALFMM_LIBRARIES FFTW)
      endif()
    else(FFTW_FOUND)
      message(WARNING "FFTW has not been found, SCALFMM will continue to compile but some applications will be disabled.")
      message(WARNING "If you have FFTW set FFTW_LIBDIR, FFTW_INCDIR or FFTW_DIR (CMake variables using -D or environment variables).")
      set(SCALFMM_USE_FFT OFF)
    endif(FFTW_FOUND)

  endif(SCALFMM_USE_FFT)
  message(STATUS "SCALFMM_USE_FFT = ${SCALFMM_USE_FFT}")

  #
  # C++ 2014
  # --------
  target_compile_features(${CMAKE_PROJECT_NAME} INTERFACE cxx_std_14)

  #
  # Additional flags for the library
  # --------------------------------
  list(APPEND SCALFMM_CXX_FLAGS -Wshadow -Wpointer-arith )
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    list(APPEND SCALFMM_CXX_FLAGS -Wcast-qual -Wconversion -Wcast-align -Woverloaded-virtual)
  endif()

  #
  # STARPU Section
  # --------------
  message( STATUS "SCALFMM_USE_STARPU = ${SCALFMM_USE_STARPU}" )
  CMAKE_DEPENDENT_OPTION( SCALFMM_USE_OPENCL "Set to ON to use OPENCL with StarPU" OFF "SCALFMM_USE_STARPU" OFF)
  CMAKE_DEPENDENT_OPTION( SCALFMM_USE_CUDA "Set to ON to use OPENCL with StarPU" OFF "SCALFMM_USE_STARPU" OFF)
  if( SCALFMM_USE_STARPU )
    enable_language(C)
    # Find StarPU with a list of optional components
    set(SCALFMM_STARPU_VERSION "1.2.8" CACHE STRING " STARPU version desired")
    # create list of components in order to make a single call to find_package(starpu...)
    set(STARPU_COMPONENT_LIST "HWLOC")
    if(SCALFMM_USE_CUDA)
      message(WARNING "This release doest not include a full support for CUDA/OpenCl.")
      #list(APPEND STARPU_COMPONENT_LIST "CUDA")
    endif()
    if(SCALFMM_USE_MPI)
      list(APPEND STARPU_COMPONENT_LIST "MPI")
    endif()
    if(SCALFMM_USE_OPENCL)
      message(WARNING "This release doest not include a full support for CUDA/OpenCl.")
      #list(APPEND STARPU_COMPONENT_LIST "OPENCL")
    endif()

    find_package(STARPU ${SCALFMM_STARPU_VERSION} REQUIRED
      COMPONENTS ${STARPU_COMPONENT_LIST})

    if(STARPU_FOUND)
      target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC starpu::starpu_dep)
      list(APPEND FUSE_LIST "STARPU")
      list(APPEND SCALFMM_LIBRARIES "STARPU")
      if(SCALFMM_USE_CUDA)
        #list(APPEND FUSE_LIST "CUDA")
      endif()
      if(SCALFMM_USE_OPENCL)
        #list(APPEND FUSE_LIST "OPENCL")
      endif()
    else(STARPU_FOUND)
      message(FATAL_ERROR "StarPU not found.")
    endif(STARPU_FOUND)
  endif(SCALFMM_USE_STARPU)

  # Additional options for developers
  if( SCALFMM_ONLY_DEVEL )
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
      option( SCALFMM_USE_OMP4 "Set to ON to disable the gcc/intel omp4"    OFF )
    else()
      option( SCALFMM_USE_OMP4 "Set to ON to disable the gcc/intel omp4"    ON )
    endif()
    option( SCALFMM_TIME_OMPTASKS "Set to ON to time omp4 tasks and generate output file"    OFF )
    # SIMGRID and peformance models options
    option( SCALFMM_SIMGRID_NODATA "Set to ON to avoid the allocation of numerical parts in the group tree" OFF )
    option( STARPU_SIMGRID_MLR_MODELS "Set to ON to enable MLR models need for calibration and simulation" OFF )
    # OPENMP 4/5 support
    message(WARNING "OpenMP_VERSION:  ${OpenMP_VERSION}")
    option( OPENMP_SUPPORT_PRIORITY  "Set to ON to enable tasks priority (OMP4.x/KSTAR/StarPU compiler only)" OFF )
    # STARPU options
    cmake_dependent_option(OPENMP_SUPPORT_COMMUTE        "Set to ON to let tasks commute (KSTAR/StarPU compiler only)" OFF
      "SCALFMM_USE_STARPU" OFF)
    cmake_dependent_option(OPENMP_SUPPORT_TASK_NAME       "Set to ON to enable a taskname clause for tasks (KSTAR/StarPU compiler only)" OFF
      "SCALFMM_USE_STARPU" OFF)
    cmake_dependent_option(SCALFMM_STARPU_USE_COMMUTE  "Set to ON to enable commute with StarPU" ON  "SCALFMM_USE_STARPU" OFF)
    cmake_dependent_option(SCALFMM_STARPU_USE_REDUX  "Set to ON to enable redux with StarPU"     OFF "SCALFMM_USE_STARPU" OFF)
    cmake_dependent_option(SCALFMM_STARPU_USE_PRIO  "Set to ON to enable priority with StarPU"   ON  "SCALFMM_USE_STARPU" OFF)
    cmake_dependent_option(SCALFMM_STARPU_FORCE_NO_SCHEDULER  "Set to ON to disable heteroprio even if supported" OFF "SCALFMM_USE_STARPU" OFF)
    cmake_dependent_option(SCALFMM_USE_STARPU_EXTRACT  "Set to ON to enable extract with StarPU mpi implicit" ON  "SCALFMM_USE_STARPU" OFF)
  endif()

  #
  # Set ScalFMM compile flags
  # -------------------------
  string(REPLACE ";" " " SCALFMM_COMPILE_FLAGS "${SCALFMM_CXX_FLAGS}")
  set(SCALFMM_CXX_FLAGS "${SCALFMM_CXX_FLAGS}" CACHE STRING "Global compile flags for ScalFMM")
  message(STATUS "ScalFMM final flags : ${SCALFMM_CXX_FLAGS}")
  # And adding flags to target
  # --------------------------
  target_compile_options(${CMAKE_PROJECT_NAME} PUBLIC "${SCALFMM_CXX_FLAGS}")

  #
  # Set ScalFMM dependencies
  # ------------------------
  string(REPLACE ";" " " SCALFMM_COMPILE_LIBS "${SCALFMM_LIBRARIES}")
  set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES}" CACHE STRING "ScalFMM libraries dependencies")
  message(STATUS "ScalFMM dependencies : ${SCALFMM_LIBRARIES}")

  #  #
  #  # Set ScalFMM link flags
  #  # ----------------------
  #  string(REPLACE ";" " " SCALFMM_LINK_FLAGS "${SCALFMM_LINK_FLAGS}")
  #  set(SCALFMM_LINK_FLAGS "${SCALFMM_LINK_FLAGS}" CACHE STRING "Global link flags for ScalFMM")
  #  message(STATUS "ScalFMM link flags : ${SCALFMM_LINK_FLAGS}")
  #  # And adding lin flags to the target
  #  # ----------------------------------
  #  target_link_options(${CMAKE_PROJECT_NAME} INTERFACE "${SCALFMM_LIBRARIES}")

  #
  # Set ScalFMM compile definitions
  # -------------------------------
  set(SCALFMM_COMPILE_DEFINITIONS "${SCALFMM_COMPILE_DEFINITIONS}" CACHE STRING "Global compile definitions for ScalFMM")
  message(STATUS "ScalFMM compile definitions : ${SCALFMM_COMPILE_DEFINITIONS}")
  # And adding definitions to target
  # --------------------------------
  target_compile_definitions(${CMAKE_PROJECT_NAME} PUBLIC "${SCALFMM_COMPILE_DEFINITIONS}")


  #
  # Generate ScalFmmConfig.h
  # ------------------------
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/ScalFmmConfig.h.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/include/ScalFmmConfig.h)

  #
  # Static analysis during build
  # ----------------------------
  option(SCALFMM_USE_STATIC_ANALYSIS "Set to ON to enable static analysis tools." OFF)
  get_target_property(SCALFMM_INCLUDE_DIRS ${CMAKE_PROJECT_NAME} INCLUDE_DIRECTORIES)

  if(SCALFMM_USE_STATIC_ANALYSIS)
    find_program(CPPCHECK "cppcheck")
    if (CPPCHECK)
      set(CMAKE_CXX_CPPCHECK "${CPPCHECK}"
        "--language=c++"
        "$<GENEX_EVAL:-I$<JOIN:$<TARGET_PROPERTY:INCLUDE_DIRECTORIES>, -I>>"
        "-v"
        "--enable=all"
        "--force"
        "--inline-suppr"
        )
      message(STATUS "CPPCHECK analysis is ON.")
    endif()

    find_program(CLANGTIDY "clang-tidy")
    if (CLANGTIDY)
      set(CMAKE_CXX_CLANG_TIDY "${CLANGTIDY};-checks=-*,boost-*,cppcoreguidelines-*,clang-analyser-cplusplus*,modernize-*,mpi-*,performance-*,portability-*,readability-*")
      message(STATUS "Clang Tidy analysis is ON.")
    endif()

    find_program(IWYU "include-what-you-use")
    if(IWYU)
      set(CMAKE_CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
      message(STATUS "Include What You Use analysis is ON.")
    endif()
  endif()

  #
  # Build Examples
  # --------------
  message(STATUS "SCALFMM_BUILD_EXAMPLES = ${SCALFMM_BUILD_EXAMPLES}" )
  if( SCALFMM_BUILD_EXAMPLES )
    # Build - Examples and drivers
    add_subdirectory(Examples)
  endif()

  #
  # Build tools (Internal use)
  # --------------------------
  # Build - Utils
  if( SCALFMM_ONLY_DEVEL AND SCALFMM_BUILD_UTILS)
    add_subdirectory(Utils)
  endif()

  #
  # Build Tests to see examples
  # ----------------------------
  # Build - Tests
  message( STATUS "SCALFMM_BUILD_TESTS = ${SCALFMM_BUILD_TESTS}" )
  if( SCALFMM_BUILD_TESTS )
    add_subdirectory(Tests)
  endif()

  #
  # Build - Unitary and numerical Tests
  # -----------------------------------
  message( STATUS "SCALFMM_BUILD_UTESTS = ${SCALFMM_BUILD_UTESTS}" )
  if( SCALFMM_BUILD_UTESTS )
    enable_testing()
    include(CTest)
    add_subdirectory(UTests)
  endif()

  #
  # Add - doc
  # ---------
  message(STATUS "SCALFMM_BUILD_DOC = ${SCALFMM_BUILD_DOC}" )
  if(SCALFMM_BUILD_DOC)
    add_subdirectory(Doc)
  endif()

  #
  # Export Library
  # --------------
  install(TARGETS ${CMAKE_PROJECT_NAME} inastemp
    EXPORT ${CMAKE_PROJECT_NAME}-targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin)

  install(DIRECTORY ${${CMAKE_PROJECT_NAME}_SOURCE_DIR}/include
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    PATTERN "ScalFmmConfig.h.cmake" EXCLUDE)

  install(FILES ${${CMAKE_PROJECT_NAME}_BINARY_DIR}/include/ScalFmmConfig.h
    DESTINATION include)

  install(EXPORT ${CMAKE_PROJECT_NAME}-targets
    FILE ${CMAKE_PROJECT_NAME}-targets.cmake
    NAMESPACE ${CMAKE_PROJECT_NAME}::
    DESTINATION lib/cmake/${CMAKE_PROJECT_NAME})

  include(CMakePackageConfigHelpers)
  write_basic_package_version_file("${CMAKE_PROJECT_NAME}ConfigVersion.cmake"
    VERSION ${${CMAKE_PROJECT_NAME}_VERSION}
    COMPATIBILITY SameMajorVersion)

  if( SCALFMM_INSTALL_DATA )
    install(FILES  ${SCALFMM_SOURCE_DIR}/Data/test20k.fma  DESTINATION Data/   )
  endif()

  configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/${CMAKE_PROJECT_NAME}Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/${CMAKE_PROJECT_NAME}/)

  install(FILES "${CMAKE_BINARY_DIR}/${CMAKE_PROJECT_NAME}ConfigVersion.cmake"
    "${PROJECT_BINARY_DIR}/${CMAKE_PROJECT_NAME}Config.cmake"
    DESTINATION lib/cmake/${CMAKE_PROJECT_NAME})

  #
  # build a CPack driven installer package
  # --------------------------------------
  # See CPackConfig.cmake for the details.
  include(CPack)

  #
  # Print configure status
  # ----------------------
  #                    Use Mem stats
  message(STATUS "SCALFMM_USE_MEM_STATS = ${SCALFMM_USE_MEM_STATS}" )
  #                    Use Log
  message(STATUS "SCALFMM_USE_LOG = ${SCALFMM_USE_LOG}" )
  #                    Use Assert
  message(STATUS "SCALFMM_USE_ASSERT = ${SCALFMM_USE_ASSERT}" )
  #
  message(STATUS "SCALFMM_USE_MPI = ${SCALFMM_USE_MPI}" )

  message(STATUS "SCALFMM_USE_BLAS = ${SCALFMM_USE_BLAS}")
  message(STATUS "SCALFMM_USE_FFT = ${SCALFMM_USE_FFT}")
  message(STATUS "SCALFMM_USE_MKL = ${SCALFMM_USE_MKL}")
  #
  if(SCALFMM_USE_STARPU)
    message(STATUS "SCALFMM_USE_STARPU = ${SCALFMM_USE_STARPU}" )
    message(STATUS "STARPU_LIBRARIES = ${STARPU_LIBRARIES}")
    message(STATUS "SCALFMM_USE_CUDA = ${SCALFMM_USE_CUDA}" )
    MESSAGE(STATUS "SCALFMM_USE_OPENCL = ${SCALFMM_USE_OPENCL}" )
  endif(SCALFMM_USE_STARPU)
  message(STATUS "SCALFMM_CXX_FLAGS = ${SCALFMM_CXX_FLAGS}")
  message(STATUS "SCALFMM_LIBRARIES = ${SCALFMM_LIBRARIES}")
  message(STATUS "FUSE_LIST = ${FUSE_LIST}")

else(MORSE_DISTRIB_DIR OR EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/morse_cmake")

  message(STATUS "MORSE_DISTRIB_DIR is not set")
  message(STATUS "Please indicate where is located your MORSE distribution directory."
    " This is necessary to find cmake_modules.")
  message(STATUS "If you are a MORSE developer, you normally have access to the svn repository of the MORSE_DISTRIB_DIR: "
    "svn checkout svn+ssh://yourlogin@scm.gforge.inria.fr/svnroot/morse/trunk/morse_distrib")

endif(MORSE_DISTRIB_DIR OR EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/morse_cmake")
