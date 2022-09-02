/*
!> @author 
!!   Nicolas Salles,
!!   Matic Poberznic,
!!   Miha Gubde

!> @brief
!!   Plugin interface for the class FixARTn

!> @note 
!!   Should be compiled with lammps/src/
*/

#include "lammpsplugin.h"
#include "version.h"
#include <cstring>
#include "fix_artn.h"

using namespace LAMMPS_NS;


static Fix *artn2creator( LAMMPS *lmp, int argc, char **argv ){

  return new FixARTn( lmp, argc, argv );

}


extern "C" void lammpsplugin_init( void *lmp, void *handle, void *regfunc ){

  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  plugin.version = LAMMPS_VERSION;
  plugin.style = "fix";
  plugin.name = "artn";
  plugin.info = "ARTn plugin fix style v1.0";
  plugin.author = "Nicoolas Salles (nsalles33@gmail.com)";
  plugin.creator.v2 = (lammpsplugin_factory2 *) &artn2creator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);

}

