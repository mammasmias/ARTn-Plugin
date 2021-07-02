include environment_variables

default : source

source :
	( cd src; $(MAKE); cd - )

clean :
	( cd src; $(MAKE) clean; cd - )

patch :
	echo "LIBOBJS += ${ART_PATH}/src/libartn.a" >> ${QE_PATH}/make.inc
#	sh patch-ARTn.sh
	cp ${QE_PATH}/PW/src/plugin_ext_forces.f90 .
	cat PW-src-modified/plugin_ext_forces.f90 > ${QE_PATH}/PW/src/plugin_ext_forces.f90
	cp ${QE_PATH}/Modules/plugin_flags.f90 .
	cat Modules-modified/plugin_flags.f90 > ${QE_PATH}/Modules/plugin_flags.f90
	cp ${QE_PATH}/Modules/plugin_arguments.f90 .
	cat Modules-modified/plugin_arguments.f90 > ${QE_PATH}/Modules/plugin_arguments.f90

	( cd ${QE_PATH}; $(MAKE) pw; cd - )

unpatch :
	#sh unpatch-ARTn.sh
	cp plugin_ext_forces.f90 ${QE_PATH}/PW/src/plugin_ext_forces.f90
	cp plugin_flags.f90 ${QE_PATH}/Modules/plugin_flags.f90
	cp plugin_arguments.f90 ${QE_PATH}/Modules/plugin_arguments.f90
	head -n -1 ${QE_PATH}/make.inc > ${QE_PATH}/make.inc_tmp
	mv ${QE_PATH}/make.inc_tmp ${QE_PATH}/make.inc
	rm plugin_ext_forces.f90
	rm plugin_flags.f90
	rm plugin_arguments.f90

	( cd ${QE_PATH}; $(MAKE) pw; cd - )


