# Makefile example -- modified by Leah

# A simple hand-made makefile for a package including applications
# built from Fortran 90 sources, taking into account the usual
# dependency cases.

# This makefile works with the GNU make command found on
# GNU/Linux systems and often called gmake on non-GNU systems
# ======================================================================
# declarations -- modify to suit your computer, compiler, paths, etc.
# ======================================================================

# FC: The compiler
FC = pgf90
# FCFLAGS: flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -Mvect=cachesize:524288 -Munroll -Mnoframe -O2 -pc 64 -Mfree # optimization flags
FCFLAGS += -Ktrap=fp -traceback -Mbounds # array bounds checking is slow!
# flags forall (e.g. look for system .mod files, required in gfortran)
#FCFLAGS += -I/usr/include

# libraries needed for linking, unused but kept for examples
#LDFLAGS = -li_need_this_lib

# specify paths
# MODELDIR: core model directory
MODELDIR = /home/nfalk/HFCACPG/prod21
# SRCDIR: directory where src code (.f90 files) are located
SRCDIR = $(MODELDIR)/src
# RUNDIR: directory from which to run the model
RUNDIR = $(MODELDIR)/run
# BINDIR: directory containing this Makefile where .o and .mod files are compiled
BINDIR = $(MODELDIR)/bin

# Executables
# PROGRAM: name of model program file (e.g. driver.f90) without the .f90 extension 
# the executable compiled in this bin directory will be named after the .f90 program file
PROGRAM = MAC
# EXE: executable name to link into run directory; can be the same as PROGRAM or not
EXE = MAC.exe

# type "make" to compile the code (compiles all)
# after compiling, link the executable (named after the PROGRAM.f90 file) into the run 
# directory under the specified EXE name above
all: $(PROGRAM)
	ln -sf $(BINDIR)/$(PROGRAM) $(RUNDIR)/$(EXE)

# =============================================================================
# rules -- modify for your file and module .f90 names, include all dependencies
# =============================================================================

# The process of building an executable is typically done in two steps:

# 1.Compilation: every source file required for our program (typically
# x.f90 or x.F90 in Fortran) is compiled into an object file (usually x.o)
# 2.Linking: the final executable file is built by "linking" together
# all the object files compiled in the previous step; in this step,
# additional pre-compiled libraries of can be added.

save_vars_no_mphys.o: constants.o constants_b.o model_arrays.o looping_stuff.o
damping.o: constants.o constants_b.o model_arrays.o
save_base_state.o: constants.o constants_b.o model_arrays.o looping_stuff.o
save_saved_timeseries.o: looping_stuff.o saved_timeseries.o
saved_timeseries.o: constants_b.o
store_saved_timeseries.o: model_arrays.o constants_b.o saved_timeseries.o
clear_model_arrays.o: model_arrays.o
looping_stuff.o:
model_arrays.o: constants.o constants_b.o looping_stuff.o
base_state_dry_isentropic.o: constants.o constants_b.o model_arrays.o looping_stuff.o
base_pertubation.o: constants.o constants_b.o model_arrays.o looping_stuff.o 
save_vars.o: constants.o constants_b.o model_arrays.o looping_stuff.o
get_k.o: constants.o constants_b.o model_arrays.o looping_stuff.o 
solve_prog.o: constants.o constants_b.o model_arrays.o looping_stuff.o 
model_integration.o: solve_prog.o constants.o constants_b.o\
 model_arrays.o save_vars_no_mphys.o get_k.o looping_stuff.o\
 store_saved_timeseries.o save_saved_timeseries.o save_base_state.o damping.o
$(PROGRAM):  base_pertubation.o model_integration.o constants.o constants_b.o model_arrays.o solve_prog.o\
get_k.o base_state_dry_isentropic.o looping_stuff.o clear_model_arrays.o store_saved_timeseries.o save_saved_timeseries.o\
save_base_state.o damping.o save_vars_no_mphys.o saved_timeseries.o

# ======================================================================
# general rules, these should not require modification
# ======================================================================

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files.
# Added $(SRCDIR) here in order to compile .mod and .o files in a separate 
# directory from where src code files are located
%.o: $(SRCDIR)/%.f90
	$(FC) $(FCFLAGS) -c $<

%.o: $(SRCDIR)/%.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

# type "make clean" to remove old .mod and .o files
clean:
	rm -f *.o *.mod *.MOD

# type "make veryclean" to remove .mod and .o files as well as remove old executables
veryclean: clean
	rm -f *~ $(PROGRAM) $(RUNDIR)/$(EXE)

