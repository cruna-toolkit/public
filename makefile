####################################################################################################
# CRUNA MAKEFILE FOR GMAKE       ###################################################################
# please do not make any changes ###################################################################
####################################################################################################

# disable change folder statement during run (GNUmake only)
MAKEFLAGS += --no-print-directory

# system environment
export HOST  = $(shell hostname)
export IAM   = $(shell whoami)
export CTIME = $(shell date)

# cruna root path
export CRUNA_ROOT = $(PWD)

# directories and names
export BIN = $(CRUNA_ROOT)/bin
export EXT = $(CRUNA_ROOT)/ext
export LIB = $(CRUNA_ROOT)/lib
export OBJ = $(CRUNA_ROOT)/obj
export SRC = $(CRUNA_ROOT)/src

include $(CRUNA_ROOT)/cruna.conf
# exports USER, CASE, GENERIC_MAKEFILE and LOCAL_MAKEFILE names
# exports user specific compiler/linker variables to:
# COPT		compiler flags (aka FOPTS)						>>> used by makefile_cruna_objects
# LOPT		linker flags (aka LDOPTS)						>>> used by this makefile

include $(CRUNA_ROOT)/makefiles/$(GENERIC_MAKEFILE)
# exports generic compiler/linker variables to:
# COPT		compiler flags (aka FOPTS)						>>> used by makefile_cruna_objects
# LOPT		linker flags (aka LDOPTS)						>>> used by this makefile

include $(CRUNA_ROOT)/makefiles/local/$(LOCAL_MAKEFILE)
# exports mashine specific compiler/linker variables to:
# LCDIR       	local compiler library directories to include, eg -I/libpath/to/include	>>> used by makefile_cruna_objects
# LCOPT       	local compiler flags and library names (aka FFLAGS)			>>> used by makefile_cruna_objects
# LLDIR		local linker library directories (aka LDFLAGS), eg -L/libpath/to/lib 	>>> used by this makefile
# LLOPT       	local linker library names (aka LDLIBS), eg -lfoo 			>>> used by this makefile

export USR = $(strip $(CRUNA_ROOT)/usr/$(USER))
export CAS = $(strip $(CRUNA_ROOT)/usr/$(USER)/$(CASE))
export TAR = $(MAKECMDGOALS)

# complog
.NOTPARALLEL: complog
.PHONY: complog
complog:
	@echo ""
	@echo "============= "
	@echo "MAKE COMPLOG: "
	@echo "============= "
	@$(MAKE) -f makefiles/makefile_complog

# reporting
.NOTPARALLEL: reporting
.PHONY: reporting
reporting: complog
	@echo ""
	@echo "=============== "
	@echo "MAKE REPORTING: "
	@echo "=============== "
	@$(MAKE) -f makefiles/makefile_reporting

#.NOTPARALLEL: dep
.PHONY: dep
dep: complog reporting
	@echo ""
	@echo "========= "
	@echo "MAKE DEP: "
	@echo "========= "
	@$(MAKE) -j -f makefiles/makefile_dep

.PHONY: cruna_objects
cruna_objects: complog reporting dep
	@echo ""
	@echo "================ "
	@echo "COMPILE SRC/CAS: "
	@echo "================ "
#       compile everything in SRC and CAS
	@cd $(SRC) && $(MAKE) -j -f ../makefiles/makefile_cruna_objects $(MAKECMDGOALS)

.PHONY: folders cruna
cruna: complog reporting dep cruna_objects
	@echo ""
	@echo "=========== "
	@echo "LINK CRUNA: "
	@echo "=========== "
#       link objects (*.o) and libs (see local makefile)
	@echo         $(FC) *.o $(LLDIR) $(LLOPT) $(LOPT) -o $(BIN)/cruna
	@cd $(OBJ) && $(FC) *.o $(LLDIR) $(LLOPT) $(LOPT) -o $(BIN)/cruna

.PHONY: folders cruna_hpc
cruna_hpc: complog reporting dep cruna_objects
	@echo ""
	@echo "=============== "
	@echo "LINK CRUNA_HPC: "
	@echo "=============== "
#       link objects (*.o) and libs (see local makefile)
	@echo         $(FC) *.o $(LLDIR) $(LLOPT) $(LOPT) -o $(BIN)/cruna_hpc
	@cd $(OBJ) && $(FC) *.o $(LLDIR) $(LLOPT) $(LOPT) -o $(BIN)/cruna_hpc

.PHONY: folders cruna_cdps
cruna_cdps: complog reporting dep cruna_objects
	@echo ""
	@echo "================ "
	@echo "LINK CRUNA_CDPS: "
	@echo "================ "
#       link objects (*.o) and libs (see local makefile)
	@echo         $(FC) *.o $(LLDIR) $(LLOPT) $(LOPT) -o $(BIN)/cruna_cdps
	@cd $(OBJ) && $(FC) *.o $(LLDIR) $(LLOPT) $(LOPT) -o $(BIN)/cruna_cdps	

####################################################################################################
# HELPER
####################################################################################################

list:
	@echo ""
	@echo "================== "
	@echo "AVAILABLE TARGETS: "
	@echo "================== "
	@LC_ALL=C $(MAKE) -pRrq -f $(firstword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/(^|\n)# Files(\n|$$)/,/(^|\n)# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | grep -E -v -e '^[^[:alnum:]]' -e '^$@$$'

# check
check: runtime_env
#	@echo "user:" $(USER)
#	@echo "case:" $(CASE)

# check and build env. file
.PHONY: runtime_env
runtime_env:
	@echo ""
	@echo "export USER=$(USER)" | tee runtime_env
	@echo "export CASE=$(CASE)" | tee -a runtime_env
	@echo "export CAS=$(CAS)" | tee -a runtime_env
	@echo "export CC=$(CC)" | tee -a runtime_env
	@echo "export FC=$(FC)" | tee -a runtime_env
	@echo "export GMAKEFILE=makefiles/$(GENERIC_MAKEFILE)" | tee -a runtime_env
	@echo "export LMAKEFILE=makefiles/local/$(LOCAL_MAKEFILE)" | tee -a runtime_env
	@echo ""
	@echo "to export these runtime variables call"
	@echo "source runtime_env"

# clean
.PHONY: clean
clean:
	@rm -f $(BIN)/cruna*
	@rm -f $(OBJ)/*.d
	@rm -f $(OBJ)/*.mod
	@rm -f $(OBJ)/*.o
	@rm -f runtime_env

# clean_all
.PHONY: clean_all
clean_all: clean
	@rm -rf $(LIB)/*
	@rm -rf $(OBJ)/*

# emacs tags
.NOTPARALLEL: tags
.PHONY: tags
tags: 
	@rm -f etags	
	@FILES=`find $$CAS $$SRC -name "*.f90"`; \
	etags -o ./etags $$FILES

####################################################################################################
# external libaries
####################################################################################################
lib: folders hdf5 fftw3 openBLAS
	@echo ""
	@echo "========== "
	@echo "MAKE LIB : "
	@echo "========== "

hdf5:
	@$(MAKE) -f makefiles/makefile_ext_hdf5

fftw3:
	@$(MAKE) -f makefiles/makefile_ext_fftw3

openBLAS:
	@$(MAKE) -f makefiles/makefile_ext_openBLAS

folders:
	@mkdir -p bin
	@mkdir -p lib
	@mkdir -p obj
	@ln -sf lib ./include

