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
#  @file AuxilaryFlags.cmake
#
#  @project MORSE
#  MORSE is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 0.9.0
#  @author Xavier Lacoste
#  @date 30-01-2015
#
# Define auxilary variables:
#  - CMAKE_Fortran_PREPROCESS_FLAGS : force C preprocessor.
#  - CMAKE_Fortran_FREEFORM_FLAG : Force free format.
#  - CMAKE_Fortran
###


IF(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  list(APPEND CMAKE_Fortran_PREPROCESS_FLAGS "-cpp")
  list(APPEND CMAKE_Fortran_FREEFORM_FLAG "-ffree-form")

ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  list(APPEND CMAKE_Fortran_PREPROCESS_FLAG "-fpp")
  list(APPEND CMAKE_Fortran_FREEFORM_FLAG "")
ENDIF()
