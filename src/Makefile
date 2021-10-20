
#include ${QE_PATH}/make.inc

OBJ  := Obj
SRC  := .
F90 := gfortran
#AR := ar

LIB_OBJECTS := \
            $(OBJ)/units.o \
            $(OBJ)/artn_params_mod.o \
            $(OBJ)/center.o \
            $(OBJ)/diag.o \
	    $(OBJ)/lanczos.o \
            $(OBJ)/displacement_validation.o \
	    $(OBJ)/pbc.o \
	    $(OBJ)/invmat3x3.o \
            $(OBJ)/push_init.o \
	    $(OBJ)/sum_force.o \
            $(OBJ)/perpforce.o \
	    $(OBJ)/move_mode.o \
      	    $(OBJ)/write_report.o \
	    $(OBJ)/write_struct.o \
            $(OBJ)/artn.o \
            $(OBJ)/artn_QE.o


FDEBUG := #-check all  #-fcheck=all

ifeq ($(F90),ifort)
  MOD   := -module $(OBJ)
  FLAGS := $(FFLAGS) $(FDEBUG) $(LAPACK_LIBS) -fPIC
else #($(F90),gfortran)
  MOD   := -I$(OBJ) -J$(OBJ) # -module $(OBJ)
  FLAGS := $(FFLAGS) $(FDEBUG) $(LAPACK_LIBS)
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


