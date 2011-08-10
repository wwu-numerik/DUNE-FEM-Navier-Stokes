#ifndef DUNE_NAVIERSTOKES_GLOBAL_DEFINES_HH
#define DUNE_NAVIERSTOKES_GLOBAL_DEFINES_HH

#ifdef HAVE_CMAKE_CONFIG
	#include "cmake_config.h"
#endif

//gcc header cyclic dependency error workaround
namespace std { class type_info; }

#define TESTCASE_NAME "TESTCASE"

#if ! defined(TESTCASE)
	#define TESTCASE TestCase3D
#endif

#if ( ( defined(SGRID) || defined(ALUGRID_SIMPLEX) ||  defined(ALUGRID_CUBE) ) && ( GRIDDIM == 3 ) ) || defined(UGGRID) || defined(YASPGRID)
	//this is no mistake, ALU is indeed only incompatible in 3d
	#define OLD_DUNE_GRID_VERSION
#endif

#if (GRIDDIM==3) //why ??
	#define MODEL_PROVIDES_LOCALFUNCTION 1
#else
	#define MODEL_PROVIDES_LOCALFUNCTION 0
#endif



#ifndef NAVIER_DATA_NAMESPACE
	#error "no data namspeace given"
#endif

//the adaption manager might be troublesome with certain gridparts/spaces, so we needed a easy way to disable it
#ifndef ENABLE_ADAPTIVE
	#define ENABLE_ADAPTIVE 1
#endif

#ifndef COMMIT
	#define COMMIT "undefined"
#endif

#include <string>
static const std::string commit_string (COMMIT);


#endif // DUNE_NAVIERSTOKES_GLOBAL_DEFINES_HH
