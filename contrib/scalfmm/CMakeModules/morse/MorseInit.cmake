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
#  @file MorseInit.cmake
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
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 13-07-2012
#
###

# This include is required to check symbols of libs in the main CMakeLists.txt
include(CheckFunctionExists)

# To colorize messages
include(ColorizeMessage)

# To find headers and libs
include(FindHeadersAndLibs)

# Some macros to print status when search for headers and libs
# PrintFindStatus.cmake is in cmake_modules/morse/find directory
include(PrintFindStatus)

# Define some auxilary flags
include(AuxilaryFlags)

# Define some variables to et info about ressources
include(Ressources)

# Add the path where we handle our FindFOO.cmake to seek for liraries
list(APPEND CMAKE_MODULE_PATH ${MORSE_CMAKE_MODULE_PATH}/find)

option(MORSE_VERBOSE_FIND_PACKAGE "Add additional messages concerning no found packages" OFF)
#message("-- ${BoldGreen}MORSE_VERBOSE_FIND_PACKAGE is set to OFF, turn it ON to get"
#        "information about packages not found${ColourReset}")

##
## @end file MorseInit.cmake
##
