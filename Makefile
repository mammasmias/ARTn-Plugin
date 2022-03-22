
include environment_variables

default : help

lib :
	( cd src; $(MAKE); cd - )

clean :
	( cd src; $(MAKE) clean; cd - )


# -------------------------------------------------------------------------- LAMMPS
patch-lammps:
	@$(call check_defined, LAMMPS_PATH)
	@cp Files_LAMMPS/*artn.*  ${LAMMPS_PATH}/src/
	@echo " ***** Fix/artn files copied"
	@echo " ***** You have to compile LAMMPS before useing it"

unpatch-lammps:
	@$(call check_defined, LAMMPS_PATH)
	@rm ${LAMMPS_PATH}/src/fix_artn.* ${LAMMPS_PATH}/src/artn.h
	@echo " ***** Fix/artn files removed from $LAMMPS_PATH/src"
	@echo " ***** You have to compile LAMMPS to have the effect"


# -------------------------------------------------------------------------- Quantum ESPRESSO
patch-qe: patch-QE
patch-QE :
	@$(call check_defined, QE_PATH)
	echo "LIBOBJS += ${ART_PATH}/src/libartn.a" >> ${QE_PATH}/make.inc
	#sh patch-ARTn.sh
	cp ${QE_PATH}/PW/src/plugin_ext_forces.f90 Files_QE/.
	cat Files_QE/PW-src-modified/plugin_ext_forces.f90 > ${QE_PATH}/PW/src/plugin_ext_forces.f90
	#cp ${QE_PATH}/Modules/plugin_flags.f90 Files_QE/.
	#cat Files_QE/Modules-modified/plugin_flags.f90 > ${QE_PATH}/Modules/plugin_flags.f90
	#cp ${QE_PATH}/Modules/plugin_arguments.f90 Files_QE/.
	#cat Files_QE/Modules-modified/plugin_arguments.f90 > ${QE_PATH}/Modules/plugin_arguments.f90

	( cd ${QE_PATH}; $(MAKE) pw; cd - )

unpatch-qe: unpatch-QE
unpatch-QE :
	@$(call check_defined, QE_PATH)
	#sh unpatch-ARTn.sh
	cp Files_QE/plugin_ext_forces.f90 ${QE_PATH}/PW/src/plugin_ext_forces.f90
	#cp Files_QE/plugin_flags.f90 ${QE_PATH}/Modules/plugin_flags.f90
	#cp Files_QE/plugin_arguments.f90 ${QE_PATH}/Modules/plugin_arguments.f90
	head -n -1 ${QE_PATH}/make.inc > ${QE_PATH}/make.inc_tmp
	mv ${QE_PATH}/make.inc_tmp ${QE_PATH}/make.inc
	rm Files_QE/plugin_ext_forces.f90
	#rm Files_QE/plugin_flags.f90
	#rm Files_QE/plugin_arguments.f90

	( cd ${QE_PATH}; $(MAKE) pw; cd - )



# ---------------------------------------------
check_defined = \
    $(strip $(foreach 1,$1, \
        $(call __check_defined,$1,$(strip $(value 2)))))
__check_defined = \
    $(if $(value $1),, \
      $(error Undefined $1$(if $2, ($2)): Define it in the file environment_variables))

verif_defined = \
    $(strip $(foreach 1,$1, \
        $(call __verif_defined,$1,$(strip $(value 2)))))
__verif_defined = \
    $(if $(value $1),, \
      echo "*WARNING** Undefined $1$(if $2, ($2)): Define it in the file environment_variables")




# ---------------------------------------------
help:
	@echo "*******************************************************************************"
	@echo "*                    Plugin-ARTn Library "
	@echo "*******************************************************************************"
	@echo ""
	@echo "* COMPILATION:"
	@$(call verif_defined, F90)
	@echo "make lib		compile the libartn.a library in src/ folder with ${F90} compiler"
	@echo "make clean		delete the object files and libartn.a from src/ "
	@echo ""
	@echo ""
	@echo "* LAMMPS Interface:"
	@$(call verif_defined, LAMMPS_PATH)
	@echo "make patch-lammps	copy the Files_LAMMPS/fix_artn.* and Files_LAMMPS/artn.h to LAMMPS_PATH/src"
	@echo "make unpatch-lammps	delete the fix_artn.* and artn.h files from  LAMMPS_PATH/src"
	@echo ""
	@echo "* QE Interface:"
	@$(call verif_defined, QE_PATH)
	@echo "make patch-qe		copy Files_QE/plugin_ext_forces.f90 to QE_PATH/src"
	@echo "make unpatch-qe		delete the changes in plugin_ext_forces.f90 from QE_PATH/src"
	@echo ""
	@echo "*******************************************************************************"
	@echo ""
	






