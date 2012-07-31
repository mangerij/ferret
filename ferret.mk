ferret_SRC_DIRS := $(FERRET_DIR)/src/*/*

ferret_INC_DIRS := $(shell find $(FERRET_DIR)/include -type d -not -path "*/.svn*")
ferret_INCLUDE  := $(foreach i, $(ferret_INC_DIRS), -I$(i))

libmesh_INCLUDE := $(ferret_INCLUDE) $(libmesh_INCLUDE)

ferret_LIB := $(FERRET_DIR)/libferret-$(METHOD)$(libext)

ferret_APP := $(FERRET_DIR)/ferret-$(METHOD)

# source files
ferret_srcfiles    := $(shell find $(ferret_SRC_DIRS) -name *.C)
ferret_csrcfiles   := $(shell find $(ferret_SRC_DIRS) -name *.c)
ferret_fsrcfiles   := $(shell find $(ferret_SRC_DIRS) -name *.f)
ferret_f90srcfiles := $(shell find $(ferret_SRC_DIRS) -name *.f90)
# object files
ferret_objects := $(patsubst %.C, %.$(obj-suffix), $(ferret_srcfiles))
ferret_objects += $(patsubst %.c, %.$(obj-suffix), $(ferret_csrcfiles))
ferret_objects += $(patsubst %.f, %.$(obj-suffix), $(ferret_fsrcfiles))
ferret_objects += $(patsubst %.f90, %.$(obj-suffix), $(ferret_f90srcfiles))

ferret_app_objects := $(patsubst %.C, %.$(obj-suffix), $(FERRET_DIR)/src/main.C)

# plugin files
ferret_plugfiles   := $(shell find $(FERRET_DIR)/plugins/ -name *.C 2>/dev/null)
ferret_cplugfiles  := $(shell find $(FERRET_DIR)/plugins/ -name *.c 2>/dev/null)
ferret_fplugfiles  := $(shell find $(FERRET_DIR)/plugins/ -name *.f 2>/dev/null)
ferret_f90plugfiles:= $(shell find $(FERRET_DIR)/plugins/ -name *.f90 2>/dev/null)

# plugins
ferret_plugins     := $(patsubst %.C, %-$(METHOD).plugin, $(ferret_plugfiles))
ferret_plugins     += $(patsubst %.c, %-$(METHOD).plugin, $(ferret_cplugfiles))
ferret_plugins     += $(patsubst %.f, %-$(METHOD).plugin, $(ferret_fplugfiles))
ferret_plugins     += $(patsubst %.f90, %-$(METHOD).plugin, $(ferret_f90plugfiles))

all:: $(ferret_LIB)

# build rule for lib FERRET
ifeq ($(enable-shared),yes)
# Build dynamic library
$(ferret_LIB): $(ferret_objects) $(ferret_plugins)
	@echo "Linking "$@"..."
	@$(libmesh_CC) $(libmesh_CXXSHAREDFLAG) -o $@ $(ferret_objects) $(libmesh_LDFLAGS)
else
# Build static library
ifeq ($(findstring darwin,$(hostos)),darwin)
$(ferret_LIB): $(ferret_objects)
	@echo "Linking "$@"..."
	@libtool -static -o $@ $(ferret_objects)
else
$(ferret_LIB): $(ferret_objects)
	@echo "Linking "$@"..."
	@$(AR) rv $@ $(ferret_objects)
endif
endif

# include FERRET dep files
-include $(FERRET_DIR)/src/*/*.d


# how to build FERRET application
ifeq ($(APPLICATION_NAME),ferret)
all:: ferret

ferret: $(ferret_APP)

$(ferret_APP): $(moose_LIB) $(elk_MODULES) $(ferret_LIB) $(ferret_app_objects)
	@echo "Linking "$@"..."
	@$(libmesh_CXX) $(libmesh_CXXFLAGS) $(ferret_app_objects) -o $@ $(ferret_LIB) $(elk_MODULES) $(moose_LIB) $(libmesh_LIBS) $(libmesh_LDFLAGS) $(ADDITIONAL_LIBS)

-include $(FERRET_DIR)/src/*.d
endif

#
# Maintenance
#
delete_list := $(ferret_APP) $(ferret_LIB)

###############################################################################
# Additional special case targets should be added here
