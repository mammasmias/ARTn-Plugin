
include environment_variables

default : help

lib :
	@$(call check_defined, F90)
	( cd src; $(MAKE); cd - )

clean :
	( cd src; $(MAKE) clean; cd - )
	@rm Files_LAMMPS/*.o 

lammps:
	@$(call check_defined, CXX)
	@$(call check_defined, LAMMPS_PATH)
	@echo cxx is: "${CXX}"
	if  echo "${CXX}" | grep -q "mpi" ; then \
	$(CXX) -fPIC -c Files_LAMMPS/fix_artn.cpp -o Files_LAMMPS/fix_artn.o -I${LAMMPS_PATH}/src; \
	$(CXX) -fPIC -c Files_LAMMPS/artnplugin.cpp -o Files_LAMMPS/artnplugin.o -I${LAMMPS_PATH}/src; \
	else \
	$(CXX) -fPIC -c Files_LAMMPS/fix_artn.cpp -o Files_LAMMPS/fix_artn.o -I${LAMMPS_PATH}/src -I${LAMMPS_PATH}/src/STUBS; \
	$(CXX) -fPIC -c Files_LAMMPS/artnplugin.cpp -o Files_LAMMPS/artnplugin.o -I${LAMMPS_PATH}/src -I${LAMMPS_PATH}/src/STUBS; \
	fi
	
	#$(CXX) -fPIC -c Files_LAMMPS/fix_artn.cpp -o Files_LAMMPS/fix_artn.o -I${LAMMPS_PATH}/src -I${LAMMPS_PATH}/src/STUBS
	#$(CXX) -fPIC -c Files_LAMMPS/artnplugin.cpp -o Files_LAMMPS/artnplugin.o -I${LAMMPS_PATH}/src -I${LAMMPS_PATH}/src/STUBS 

sharelib: lib lammps
	@echo "";echo ">>>> environment_variable verification"
	@$(call check_defined, CC)
	@$(call check_defined, FORT_LIB)
	@$(call check_defined, BLAS_LIB)
	@echo "<<<< OK "; echo ""
	${CC} -shared -rdynamic -o libartn.so src/Obj/*.o Files_LAMMPS/*.o $(FORT_LIB) $(BLAS_LIB)
	@echo ">>>> Shared library done" ; echo ""
	@echo " 1) In LAMMPS Package PLUGIN must be loaded"
	@echo " 2) LAMMPS must be compiled with mode=shared"
	@echo " 3) To launch LAMMPS, lammps/src path should be loaded in LD_LIBRARY_PATH"
	@echo " 4) Enjoy ;) "
	@echo ""

# -------------------------------------------------------------------------- LAMMPS
patch-lammps:
	@$(call check_defined, LAMMPS_PATH)
	@cp Files_LAMMPS/*artn.*  ${LAMMPS_PATH}/src/
	@echo " ***** File artnplugin.cpp copied"
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
	echo "QELIBS += ${ART_PATH}/src/libartn.a" >> ${QE_PATH}/make.inc
	cp ${QE_PATH}/PW/src/plugin_ext_forces.f90 Files_QE/.
	cat Files_QE/PW-src-modified/plugin_ext_forces.f90 > ${QE_PATH}/PW/src/plugin_ext_forces.f90

	( cd ${QE_PATH}; $(MAKE) pw; cd - )

unpatch-qe: unpatch-QE
unpatch-QE :
	@$(call check_defined, QE_PATH)
	#sh unpatch-ARTn.sh
	cp Files_QE/plugin_ext_forces.f90 ${QE_PATH}/PW/src/plugin_ext_forces.f90
	head -n -1 ${QE_PATH}/make.inc > ${QE_PATH}/make.inc_tmp
	mv ${QE_PATH}/make.inc_tmp ${QE_PATH}/make.inc
	rm Files_QE/plugin_ext_forces.f90

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
	@echo ""
	@echo "*******************************************************************************"
	@echo "*                    Plugin-ARTn Library "
	@echo "*******************************************************************************"
	@echo ""
	@echo "* WARNNG: The user should take care to fill correctly the file environment_variable"
	@echo "          before to compile ARTn"
	@echo ""
	@echo "* COMPILATION:"
	@$(call verif_defined, F90)
	@echo "make lib		compile the libartn.a library in src/ folder with ${F90} compiler"
	@echo "make clean		delete the object files and libartn.a from src/ "
	@echo ""
	@echo ""
	@echo "* LAMMPS Interface:"
	@$(call verif_defined, LAMMPS_PATH)
	@echo "make sharelib		compile dynamic library libartn.so with plugin interfaces for LAMMPS"
	@#echo "make patch-lammps	copy the Files_LAMMPS/fix_artn.* and Files_LAMMPS/artn.h to LAMMPS_PATH/src"
	@#echo "make unpatch-lammps	delete the fix_artn.* and artn.h files from  LAMMPS_PATH/src"
	@echo ""
	@echo "* QE Interface:"
	@$(call verif_defined, QE_PATH)
	@echo "make patch-qe		copy Files_QE/plugin_ext_forces.f90 to QE_PATH/src"
	@echo "make unpatch-qe		delete the changes in plugin_ext_forces.f90 from QE_PATH/src"
	@echo ""
	@echo "*******************************************************************************"
	@echo ""
	






