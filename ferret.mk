ferret_SRC_DIRS := $(FERRET_DIR)/src/*/*

ferret_INC_DIRS := $(shell find $(FERRET_DIR)/include -type d -not -path "*/.svn*")
ferret_INCLUDE  := $(foreach i, $(ferret_INC_DIRS), -I$(i))

libmesh_INCLUDE := $(ferret_INCLUDE) $(libmesh_INCLUDE)

ferret_LIB := $(FERRET_DIR)/libferret-$(METHOD).la

ferret_APP := $(FERRET_DIR)/ferret-$(METHOD)

# source files
ferret_srcfiles    := $(shell find $(ferret_SRC_DIRS) -name "*.C")
ferret_csrcfiles   := $(shell find $(ferret_SRC_DIRS) -name "*.c")
ferret_fsrcfiles   := $(shell find $(ferret_SRC_DIRS) -name "*.f")
ferret_f90srcfiles := $(shell find $(ferret_SRC_DIRS) -name "*.f90")
# object files
ferret_objects := $(patsubst %.C, %.$(obj-suffix), $(ferret_srcfiles))
ferret_objects += $(patsubst %.c, %.$(obj-suffix), $(ferret_csrcfiles))
ferret_objects += $(patsubst %.f, %.$(obj-suffix), $(ferret_fsrcfiles))
ferret_objects += $(patsubst %.f90, %.$(obj-suffix), $(ferret_f90srcfiles))

# plugin files
ferret_plugfiles   := $(shell find $(FERRET_DIR)/plugins/ -name "*.C" 2>/dev/null)
ferret_cplugfiles  := $(shell find $(FERRET_DIR)/plugins/ -name "*.c" 2>/dev/null)
ferret_fplugfiles  := $(shell find $(FERRET_DIR)/plugins/ -name "*.f" 2>/dev/null)
ferret_f90plugfiles:= $(shell find $(FERRET_DIR)/plugins/ -name "*.f90" 2>/dev/null)

# plugins
ferret_plugins     := $(patsubst %.C, %-$(METHOD).plugin, $(ferret_plugfiles))
ferret_plugins     += $(patsubst %.c, %-$(METHOD).plugin, $(ferret_cplugfiles))
ferret_plugins     += $(patsubst %.f, %-$(METHOD).plugin, $(ferret_fplugfiles))
ferret_plugins     += $(patsubst %.f90, %-$(METHOD).plugin, $(ferret_f90plugfiles))

# ferret main
ferret_main_src    := $(FERRET_DIR)/src/main.C
ferret_app_objects := $(patsubst %.C, %.$(obj-suffix), $(ferret_main_src))

# dependency files
ferret_deps := $(patsubst %.C, %.$(obj-suffix).d, $(ferret_srcfiles)) \
               $(patsubst %.c, %.$(obj-suffix).d, $(ferret_csrcfiles)) \
               $(patsubst %.C, %.$(obj-suffix).d, $(ferret_main_src))

# If building shared libs, make the plugins a dependency, otherwise don't.
ifeq ($(libmesh_shared),yes)
  ferret_plugin_deps := $(ferret_plugins)
else
  ferret_plugin_deps :=
endif

all:: $(ferret_LIB)

$(ferret_LIB): $(ferret_objects) $(ferret_plugin_deps)
	@echo "Linking "$@"..."
	@$(libmesh_LIBTOOL) --tag=CXX $(LIBTOOLFLAGS) --mode=link --quiet \
	  $(libmesh_CXX) $(libmesh_CXXFLAGS) -o $@ $(ferret_objects) $(libmesh_LIBS) $(libmesh_LDFLAGS) $(EXTERNAL_FLAGS) -rpath $(FERRET_DIR)
	@$(libmesh_LIBTOOL) --mode=install --quiet install -c $(ferret_LIB) $(FERRET_DIR)

# include FERRET dep files
-include $(ferret_deps)

# how to build FERRET application
ifeq ($(APPLICATION_NAME),ferret)
all:: ferret

ferret: $(ferret_APP)

$(ferret_APP): $(moose_LIB) $(elk_MODULES) $(ferret_LIB) $(ferret_app_objects)
	@echo "Linking "$@"..."
	@$(libmesh_LIBTOOL) --tag=CXX $(LIBTOOLFLAGS) --mode=link --quiet \
          $(libmesh_CXX) $(libmesh_CXXFLAGS) -o $@ $(ferret_app_objects) $(ferret_LIB) $(elk_MODULES) $(moose_LIB) $(libmesh_LIBS) $(libmesh_LDFLAGS) $(ADDITIONAL_LIBS)

endif

#
# Maintenance
#
delete_list := $(ferret_APP) $(ferret_LIB) $(FERRET_DIR)/libferret-$(METHOD).*

###############################################################################
# Additional special case targets should be added here
