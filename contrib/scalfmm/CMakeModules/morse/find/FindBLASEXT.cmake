###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2014 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
# - Find BLAS EXTENDED for MORSE projects: find include dirs and libraries
#
# This module allows to find BLAS libraries by calling the official FindBLAS module
# and handles the creation of different library lists whether the user wishes to link
# with a sequential BLAS or a multihreaded (BLAS_SEQ_LIBRARIES and BLAS_PAR_LIBRARIES).
# BLAS is detected with a FindBLAS call then if the BLAS vendor is Intel10_64lp or ACML
# then the module tries to find the corresponding multithreaded libraries.
#
# The following variables have been added to manage links with sequential or multithreaded
# versions:
#  BLAS_INCLUDE_DIRS  - BLAS include directories
#  BLAS_LIBRARY_DIRS  - Link directories for BLAS libraries
#  BLAS_SEQ_LIBRARIES - BLAS component libraries to be linked (sequential)
#  BLAS_PAR_LIBRARIES - BLAS component libraries to be linked (multithreaded)

#=============================================================================
# Copyright 2012-2013 Inria
# Copyright 2012-2013 Emmanuel Agullo
# Copyright 2012-2013 Mathieu Faverge
# Copyright 2012      Cedric Castagnede
# Copyright 2013      Florent Pruvost
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file MORSE-Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of Morse, substitute the full
#  License text for the above reference.)


# add a cache variable to let the user specify the BLAS vendor
set(BLA_VENDOR "" CACHE STRING "list of possible BLAS vendor:
    Open, Eigen, Goto, ATLAS PhiPACK, CXML, DXML, SunPerf, SCSL, SGIMATH, IBMESSL,
    Intel10_32 (intel mkl v10 32 bit),
    Intel10_64lp (intel mkl v10 64 bit, lp thread model, lp64 model),
    Intel10_64lp_seq (intel mkl v10 64 bit, sequential code, lp64 model),
    Intel( older versions of mkl 32 and 64 bit),
    ACML, ACML_MP, ACML_GPU, Apple, NAS, Generic")


if (NOT BLAS_FOUND)
    # First try to detect two cases:
    # 1: only SEQ libs are handled
    # 2: both SEQ and PAR libs are handled
    if(BLASEXT_FIND_REQUIRED)
        find_package(BLAS REQUIRED)
    else()
        find_package(BLAS)
    endif()
endif ()

# detect the cases where SEQ and PAR libs are handled
if(BLA_VENDOR STREQUAL "All" AND
   (BLAS_mkl_core_LIBRARY OR BLAS_mkl_core_dll_LIBRARY)
   )
    set(BLA_VENDOR "Intel")
    if(BLAS_mkl_intel_LIBRARY)
        set(BLA_VENDOR "Intel10_32")
    endif()
    if(BLAS_mkl_intel_lp64_LIBRARY)
        set(BLA_VENDOR "Intel10_64lp")
    endif()
    if(NOT BLASEXT_FIND_QUIETLY)
        message(STATUS "A BLAS library has been found (${BLAS_LIBRARIES}) but we"
            "have also potentially detected some BLAS libraries from the MKL."
            "We try to use this one.")
        message(STATUS "If you want to force the use of one specific library, "
            "please specify the BLAS vendor by setting -DBLA_VENDOR=blas_vendor_name"
            "at cmake configure.")
        message(STATUS "List of possible BLAS vendor: Goto, ATLAS PhiPACK, CXML, "
            "DXML, SunPerf, SCSL, SGIMATH, IBMESSL, Intel10_32 (intel mkl v10 32 bit),"
            "Intel10_64lp (intel mkl v10 64 bit, lp thread model, lp64 model),"
            "Intel10_64lp_seq (intel mkl v10 64 bit, sequential code, lp64 model),"
            "Intel( older versions of mkl 32 and 64 bit),"
            "ACML, ACML_MP, ACML_GPU, Apple, NAS, Generic")
    endif()
    set(BLAS_FOUND "")
elseif(BLA_VENDOR STREQUAL "All" AND BLAS_acml_LIBRARY)
    set(BLA_VENDOR "ACML")
    if(NOT BLASEXT_FIND_QUIETLY)
        message(STATUS "A BLAS library has been found (${BLAS_LIBRARIES}) but we"
            "have also potentially detected some BLAS libraries from the ACML."
            "We try to use this one.")
        message(STATUS "If you want to force the use of one specific library, "
            "please specify the BLAS vendor by setting -DBLA_VENDOR=blas_vendor_name"
            "at cmake configure.")
        message(STATUS "List of possible BLAS vendor: Goto, ATLAS PhiPACK, CXML, "
            "DXML, SunPerf, SCSL, SGIMATH, IBMESSL, Intel10_32 (intel mkl v10 32 bit),"
            "Intel10_64lp (intel mkl v10 64 bit, lp thread model, lp64 model),"
            "Intel10_64lp_seq (intel mkl v10 64 bit, sequential code, lp64 model),"
            "Intel( older versions of mkl 32 and 64 bit),"
            "ACML, ACML_MP, ACML_GPU, Apple, NAS, Generic")
    endif()
    set(BLAS_FOUND "")
endif()

# Intel case
if(BLA_VENDOR MATCHES "Intel*")

    ###
    # look for include path if the BLAS vendor is Intel
    ###

    # gather system include paths
    unset(_inc_env)
    if(WIN32)
        string(REPLACE ":" ";" _inc_env "$ENV{INCLUDE}")
    else()
        string(REPLACE ":" ";" _path_env "$ENV{INCLUDE}")
        list(APPEND _inc_env "${_path_env}")
        string(REPLACE ":" ";" _path_env "$ENV{C_INCLUDE_PATH}")
        list(APPEND _inc_env "${_path_env}")
        string(REPLACE ":" ";" _path_env "$ENV{CPATH}")
        list(APPEND _inc_env "${_path_env}")
        string(REPLACE ":" ";" _path_env "$ENV{INCLUDE_PATH}")
        list(APPEND _inc_env "${_path_env}")
    endif()
    list(APPEND _inc_env "${CMAKE_PLATFORM_IMPLICIT_INCLUDE_DIRECTORIES}")
    list(APPEND _inc_env "${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES}")
    list(REMOVE_DUPLICATES _inc_env)

    # find mkl.h inside known include paths
    set(BLAS_mkl.h_INCLUDE_DIRS "BLAS_mkl.h_INCLUDE_DIRS-NOTFOUND")
    if(BLAS_INCDIR)
        set(BLAS_mkl.h_INCLUDE_DIRS "BLAS_mkl.h_INCLUDE_DIRS-NOTFOUND")
        find_path(BLAS_mkl.h_INCLUDE_DIRS
                NAMES mkl.h
                HINTS ${BLAS_INCDIR})
    else()
        if(BLAS_DIR)
            set(BLAS_mkl.h_INCLUDE_DIRS "BLAS_mkl.h_INCLUDE_DIRS-NOTFOUND")
            find_path(BLAS_mkl.h_INCLUDE_DIRS
                    NAMES mkl.h
                    HINTS ${BLAS_DIR}
                    PATH_SUFFIXES include)
        else()
            set(BLAS_mkl.h_INCLUDE_DIRS "BLAS_mkl.h_INCLUDE_DIRS-NOTFOUND")
            find_path(BLAS_mkl.h_INCLUDE_DIRS
                    NAMES mkl.h
                    HINTS ${_inc_env})
        endif()
    endif()
    mark_as_advanced(BLAS_mkl.h_INCLUDE_DIRS)
    ## Print status if not found
    ## -------------------------
    #if (NOT BLAS_mkl.h_INCLUDE_DIRS AND MORSE_VERBOSE)
    #    Print_Find_Header_Status(blas mkl.h)
    #endif ()
    set(BLAS_INCLUDE_DIRS "")
    if(BLAS_mkl.h_INCLUDE_DIRS)
        list(APPEND BLAS_INCLUDE_DIRS "${BLAS_mkl.h_INCLUDE_DIRS}" )
    endif()

    ###
    # look for libs
    ###
    # if Intel 10 64 bit -> look for sequential and multithreaded versions
    if(BLA_VENDOR MATCHES "Intel10_64lp*")

        ## look for the sequential version
        set(BLA_VENDOR "Intel10_64lp_seq")
        if(NOT BLASEXT_FIND_QUIETLY)
            message(STATUS "Look for the sequential version Intel10_64lp_seq")
        endif()
#         if(NOT BLAS_FOUND AND BLASEXT_FIND_REQUIRED)
        if(BLASEXT_FIND_REQUIRED)
            find_package(BLAS REQUIRED)
        else()
            find_package(BLAS)
        endif()
        if(BLAS_FOUND)
            set(BLAS_SEQ_LIBRARIES "${BLAS_LIBRARIES}")
        else()
            set(BLAS_SEQ_LIBRARIES "${BLAS_SEQ_LIBRARIES-NOTFOUND}")
        endif()

        ## look for the multithreaded version
        set(BLA_VENDOR "Intel10_64lp")
        if(NOT BLASEXT_FIND_QUIETLY)
            message(STATUS "Look for the multithreaded version Intel10_64lp")
        endif()
        find_package(BLAS)
        if(BLAS_FOUND)
            set(BLAS_PAR_LIBRARIES "${BLAS_LIBRARIES}")
        else()
            set(BLAS_PAR_LIBRARIES "${BLAS_PAR_LIBRARIES-NOTFOUND}")
        endif()

    else()

        if(BLAS_FOUND)
            set(BLAS_SEQ_LIBRARIES "${BLAS_LIBRARIES}")
        else()
            set(BLAS_SEQ_LIBRARIES "${BLAS_SEQ_LIBRARIES-NOTFOUND}")
        endif()

    endif()

# ACML case
elseif(BLA_VENDOR MATCHES "ACML*")

    ## look for the sequential version
    set(BLA_VENDOR "ACML")
    if(BLASEXT_FIND_REQUIRED)
        find_package(BLAS REQUIRED)
    else()
        find_package(BLAS)
    endif()
    if(BLAS_FOUND)
        set(BLAS_SEQ_LIBRARIES "${BLAS_LIBRARIES}")
    else()
        set(BLAS_SEQ_LIBRARIES "${BLAS_SEQ_LIBRARIES-NOTFOUND}")
    endif()

    ## look for the multithreaded version
    set(BLA_VENDOR "ACML_MP")
    if(BLASEXT_FIND_REQUIRED)
        find_package(BLAS REQUIRED)
    else(BLAS_FOUND)
        find_package(BLAS)
    endif()
    if(BLAS_FOUND)
        set(BLAS_PAR_LIBRARIES "${BLAS_LIBRARIES}")
    else()
        set(BLAS_PAR_LIBRARIES "${BLAS_PAR_LIBRARIES-NOTFOUND}")
    endif()

else()

    if(BLAS_FOUND)
        # define the SEQ libs as the BLAS_LIBRARIES
        set(BLAS_SEQ_LIBRARIES "${BLAS_LIBRARIES}")
    else()
        set(BLAS_SEQ_LIBRARIES "${BLAS_SEQ_LIBRARIES-NOTFOUND}")
    endif()
    set(BLAS_PAR_LIBRARIES "${BLAS_PAR_LIBRARIES-NOTFOUND}")

endif()


if(BLAS_SEQ_LIBRARIES)
    set(BLAS_LIBRARIES "${BLAS_SEQ_LIBRARIES}")
endif()
# message("DEBUG BLAS: ${BLAS_SEQ_LIBRARIES} ${BLAS_LIBRARIES}")

# extract libs paths
# remark: because it is not given by find_package(BLAS)
set(BLAS_LIBRARY_DIRS "")
foreach(blas_lib ${BLAS_LIBRARIES})
    get_filename_component(a_blas_lib_dir "${blas_lib}" PATH)
    list(APPEND BLAS_LIBRARY_DIRS "${a_blas_lib_dir}" )
endforeach()
if (BLAS_LIBRARY_DIRS)
    list(REMOVE_DUPLICATES BLAS_LIBRARY_DIRS)
endif ()


# message(STATUS "BLAS_FOUND: ${BLAS_FOUND}")
# message(STATUS "BLA_VENDOR: ${BLA_VENDOR}")
# message(STATUS "BLAS_LIBRARIES: ${BLAS_LIBRARIES}")
# message(STATUS "BLAS_SEQ_LIBRARIES: ${BLAS_SEQ_LIBRARIES}")
# message(STATUS "BLAS_PAR_LIBRARIES: ${BLAS_PAR_LIBRARIES}")
# message(STATUS "BLAS_INCLUDE_DIRS: ${BLAS_INCLUDE_DIRS}")
# message(STATUS "BLAS_LIBRARY_DIRS: ${BLAS_LIBRARY_DIRS}")

# check that BLAS has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
if(BLA_VENDOR MATCHES "Intel*")
    if(BLA_VENDOR MATCHES "Intel10_64lp*")
        if(NOT BLASEXT_FIND_QUIETLY)
            message(STATUS "BLAS found is Intel MKL:"
                           "we manage two lists of libs,"
                           " one sequential and one parallel if found (see BLAS_SEQ_LIBRARIES and BLAS_PAR_LIBRARIES)")
            message(STATUS "BLAS sequential libraries stored in BLAS_SEQ_LIBRARIES")
        endif()
        find_package_handle_standard_args(BLAS DEFAULT_MSG
                                          BLAS_SEQ_LIBRARIES
                                          BLAS_LIBRARY_DIRS
                                          BLAS_INCLUDE_DIRS)
        if(BLAS_PAR_LIBRARIES)
            if(NOT BLASEXT_FIND_QUIETLY)
                message(STATUS "BLAS parallel libraries stored in BLAS_PAR_LIBRARIES")
            endif()
            find_package_handle_standard_args(BLAS DEFAULT_MSG
                                              BLAS_PAR_LIBRARIES)
        endif()
    else()
        if(NOT BLASEXT_FIND_QUIETLY)
            message(STATUS "BLAS sequential libraries stored in BLAS_SEQ_LIBRARIES")
        endif()
        find_package_handle_standard_args(BLAS DEFAULT_MSG
                                          BLAS_SEQ_LIBRARIES
                                          BLAS_LIBRARY_DIRS
                                          BLAS_INCLUDE_DIRS)
    endif()
elseif(BLA_VENDOR MATCHES "ACML*")
    if(NOT BLASEXT_FIND_QUIETLY)
        message(STATUS "BLAS found is ACML:"
                        "we manage two lists of libs,"
                        " one sequential and one parallel if found (see BLAS_SEQ_LIBRARIES and BLAS_PAR_LIBRARIES)")
        message(STATUS "BLAS sequential libraries stored in BLAS_SEQ_LIBRARIES")
    endif()
    find_package_handle_standard_args(BLAS DEFAULT_MSG
                                      BLAS_SEQ_LIBRARIES
                                      BLAS_LIBRARY_DIRS
                                      BLAS_INCLUDE_DIRS)
    if(BLAS_PAR_LIBRARIES)
        if(NOT BLASEXT_FIND_QUIETLY)
            message(STATUS "BLAS parallel libraries stored in BLAS_PAR_LIBRARIES")
        endif()
        find_package_handle_standard_args(BLAS DEFAULT_MSG
                                          BLAS_PAR_LIBRARIES)
    endif()
else()
    if(NOT BLASEXT_FIND_QUIETLY)
        message(STATUS "BLAS sequential libraries stored in BLAS_SEQ_LIBRARIES")
    endif()
    find_package_handle_standard_args(BLAS DEFAULT_MSG
                                      BLAS_SEQ_LIBRARIES
                                      BLAS_LIBRARY_DIRS)
endif()
