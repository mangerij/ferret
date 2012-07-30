###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Required Environment variables
# LIBMESH_DIR	- location of the libMesh library
#
# Required Make variables
# APP_NAME	- the name of this application (all lower case)
# MOOSE_DIR	- location of the MOOSE framework
# ELK_DIR	- location of ELK (if enabled)
# BISON_DIR	- location of BISON
# MARMOT_DIR	- location of MARMOT
#
# Optional Environment variables
# CURR_DIR	- current directory (DO NOT MODIFY THIS VARIABLE)
#
#
# Note: Make sure that there is no whitespace after the word 'yes' if enabling
# an application
###############################################################################
CURR_DIR        ?= $(shell pwd)
ROOT_DIR        ?= $(shell dirname `pwd`)

ifeq ($(MOOSE_DEV),true)
	MOOSE_DIR ?= $(ROOT_DIR)/devel/moose
else
	MOOSE_DIR ?= $(ROOT_DIR)/moose
endif

LIBMESH_DIR     ?= $(ROOT_DIR)/libmesh
ELK_DIR         ?= $(ROOT_DIR)/elk
FERRET_DIR     ?= $(ROOT_DIR)/ferret

APPLICATION_NAME := ferret

DEP_APPS    ?= $(shell $(MOOSE_DIR)/scripts/find_dep_apps.py $(APPLICATION_NAME))

################################## ELK MODULES ################################
ALL_ELK_MODULES := yes
###############################################################################


include $(MOOSE_DIR)/build.mk

include $(MOOSE_DIR)/moose.mk
include $(ELK_DIR)/elk.mk
include $(FERRET_DIR)/ferret.mk

###############################################################################
# Additional special case targets should be added here
