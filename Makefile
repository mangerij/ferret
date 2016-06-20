###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project
# HERD_TRUNK_DIR   - Location of the HERD repository
# FRAMEWORK_DIR    - Location of the MOOSE framework
#
###############################################################################
MOOSE_DIR          ?= $(shell dirname `pwd`)/moose
FRAMEWORK_DIR      ?= $(MOOSE_DIR)/framework
ifndef APPLICATION_DIR
  FERRET_DIR := $(shell pwd)
else
  FERRET_DIR := $(APPLICATION_DIR)
endif
###############################################################################

# framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk

################################## MODULES ####################################
TENSOR_MECHANICS := yes
PHASE_FIELD := yes
MISC := yes
include           $(MOOSE_DIR)/modules/modules.mk
###############################################################################

# dep apps
# Observe that APPLICATION_DIR is defined, but incorrectly -- left over from the
# modules.mk definitions. Must define it explicitly, that's why we have a stashed
# FERRET_DIR around.
APPLICATION_DIR    := $(FERRET_DIR)
APPLICATION_NAME   := ferret
BUILD_EXEC         := yes
DEP_APPS           := $(shell $(FRAMEWORK_DIR)/scripts/find_dep_apps.py $(APPLICATION_NAME))
include            $(FRAMEWORK_DIR)/app.mk

###############################################################################
# Additional special case targets should be added here
