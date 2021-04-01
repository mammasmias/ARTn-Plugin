
QE_PATH = /home/mgunde/qe/artn-plugin-qe/qe-6.6
export QE_PATH

include ${QE_PATH}/make.inc

default : source

source :
	( cd src; $(MAKE); cd - )

clean :
	( cd src; $(MAKE) clean; cd - )

patch :
	sh patch-ARTn.sh

unpatch :
	sh unpatch-ARTn.sh
