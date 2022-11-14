
include ../environment_variables

OBJ := Obj
SRC := .
#F90 := gfortran #ifort
#AR := ar

LIB_OBJECTS := \
	    $(OBJ)/units.o \
            $(OBJ)/artn_params_mod.o \
            $(OBJ)/center.o \
            $(OBJ)/diag.o \
	    $(OBJ)/lanczos.o \
            $(OBJ)/check_force_convergence.o \
            $(OBJ)/displacement_validation.o \
	    $(OBJ)/pbc.o \
	    $(OBJ)/invmat3x3.o \
            $(OBJ)/push_init.o \
	    $(OBJ)/sum_force.o \
            $(OBJ)/perpforce.o \
            $(OBJ)/smooth_interpol.o \
            $(OBJ)/push_over_procedure.o \
	    $(OBJ)/move_mode.o \
	    $(OBJ)/write_report.o \
	    $(OBJ)/write_struct.o \
	    $(OBJ)/restart.o \
            $(OBJ)/nextmin.o \
            $(OBJ)/get_functions.o \
            $(OBJ)/make_filename.o \
            $(OBJ)/read_guess.o \
            $(OBJ)/start_guess.o \
            $(OBJ)/nperp_limitation.o \
            $(OBJ)/fperp_min_alignment.o \
            $(OBJ)/clean_artn.o \
            $(OBJ)/artn.o \
            $(OBJ)/artn_QE.o
	    #$(OBJ)/DEBUG.o \


FDEBUG := -g -fcheck=all -Wall #-check all  #-fcheck=all

ifeq ($(F90),ifort)
  MOD   := -module $(OBJ)
  FLAGS := $(FFLAGS) $(FDEBUG) $(LAPACK_LIBS) -fPIC
else #($(F90),gfortran)
  MOD   := -I$(OBJ) -J$(OBJ) 
  FLAGS := $(FFLAGS) $(FDEBUG) $(LAPACK_LIBS) -fPIC
endif


# ---------------------------------------------------------------------------

default: lib


lib: create-folder libartn.a

# --- create folder for .o files
create-folder:
	@if [ ! -d $(OBJ) ]; then mkdir $(OBJ) ; fi
# ----------------------------------------------------------------------------


libartn.a: $(LIB_OBJECTS)
	$(AR) $(ARFLAGS) $@ $?

# ...............................................................................

$(OBJ)/%.o: $(SRC)/%.f90
	$(F90) $(MOD) -c $(FLAGS)  $^ -o $@


# --- CLEANUP
clean:
	rm -rf *.o *.a $(OBJ)


