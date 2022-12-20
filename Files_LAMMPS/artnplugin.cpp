
#include "lammpsplugin.h"
#include "version.h"
#include <cstring>
#include "fix_artn.h"

using namespace LAMMPS_NS;


/**
 * @author
 *   Matic Poberznic
 *   Miha Gunde
 *   Nicolas Salles
 *
 * @brief 
 *   Function return the class FixARTn
 * 
 * @param[in]  lmp        PTR, on the Class lammps
 * @param[in]  argc       INT, number words are in arguments of the executable
 * @param[in]  argv       CHAR, Array of string of the argument
 * @return  pointer on Class FixARTn
 *
 * @ingroup Interface
 * @snippet artnplugin.cpp creator 
 */
static Fix *artn2creator( LAMMPS *lmp, int argc, char **argv ){
  //! [creator]
  return new FixARTn( lmp, argc, argv );
  //! [creator]
}



/**
 * @author
 *   Matic Poberznic
 *   Miha Gunde
 *   Nicolas Salles
 *
 * @brief 
 *   Interface between PLUGIN of lammps amd the shared library of ARTn
 * 
 * @param[in]  lmp        VOID PTR, on the Class lammps
 * @param[in]  handle     VOID PTR, on the shared library (dynamic)
 * @param[out] regfunc    VOID PTR, on the constructor of FixARTn
 *
 * @ingroup Interface
 * @snippet artnplugin.cpp init 
 */
extern "C" void lammpsplugin_init( void *lmp, void *handle, void *regfunc ){
//! [init]
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  plugin.version = LAMMPS_VERSION;
  plugin.style = "fix";
  plugin.name = "artn";
  plugin.info = "ARTn plugin fix style v1.0";
  plugin.author = "Nicolas Salles (nsalles33@gmail.com)";
  plugin.creator.v2 = (lammpsplugin_factory2 *) &artn2creator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);
//! [init]
}

