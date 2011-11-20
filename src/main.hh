#ifndef MAIN_H
#define MAIN_H

/**
 *  \file   dune_stokes.cc
 *
 *  \brief  brief
 **/
#ifdef HAVE_CMAKE_CONFIG
    #include "cmake_config.h"
#endif
#include <dune/navier/global_defines.hh>

#include <cstdio>
#if defined(USE_PARDG_ODE_SOLVER) && defined(USE_BFG_CG_SCHEME)
    #warning ("USE_PARDG_ODE_SOLVER enabled, might conflict with custom solvers")
#endif

#if defined(UGGRID) && defined(DEBUG)
    #warning ("UGGRID in debug mode is likely to produce a segfault")
#endif

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include <dune/fem/misc/mpimanager.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/grid/common/capabilities.hh>

#include <dune/grid/utility/gridtype.hh>
typedef Dune::GridSelector::GridType
    GridType;

#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/pass/pass.hh>
#include <dune/fem/function/adaptivefunction.hh> // for AdaptiveDiscreteFunction
#include <dune/fem/misc/gridwidth.hh>

#include <dune/stokes/discretestokesfunctionspacewrapper.hh>
#include <dune/stokes/discretestokesmodelinterface.hh>
#include <dune/stokes/stokespass.hh>
#include <dune/stokes/boundarydata.hh>

#include <dune/stuff/printing.hh>
#include <dune/stuff/femeoc.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/logging.hh>
#include <dune/stuff/parametercontainer.hh>
#include <dune/stuff/postprocessing.hh>
#include <dune/stuff/profiler.hh>
#include <dune/stuff/timeseries.hh>
#include <dune/stuff/signals.hh>
#include <dune/stuff/runinfo.hh>

#include <dune/navier/thetascheme_runner.hh>

#if ENABLE_MPI
        typedef Dune::CollectiveCommunication< MPI_Comm > CollectiveCommunication;
#else
        typedef Dune::CollectiveCommunication< double > CollectiveCommunication;
#endif

//! the strings used for column headers in tex output
typedef std::vector<std::string>
    ColumnHeaders;

void eocCheck( const Stuff::RunInfoVector& runInfos )
{
    bool ups = false;
    Stuff::RunInfoVector::const_iterator it = runInfos.begin();
    Stuff::RunInfo last = *it;
    ++it;
    for ( ; it != runInfos.end(); ++it ) {
        ups = ( last.L2Errors[0] < it->L2Errors[0]
            || last.L2Errors[1] < it->L2Errors[1] );
        last = *it;
    }
    if ( ups ) {
        Logger().Err() 	<< 	"----------------------------------------------------------\n"
                        <<	"-                                                        -\n"
                        <<	"-                  negative EOC                          -\n"
                        <<	"-                                                        -\n"
                        << 	"----------------------------------------------------------\n"
                        << std::endl;
    }
}

CollectiveCommunication init( int argc, char** argv )
{
    Dune::MPIManager::initialize(argc, argv);
    assert( Dune::Capabilities::isParallel< GridType >::v );


    /* ********************************************************************** *
     * initialize all the stuff we need                                       *
     * ********************************************************************** */
    if ( argc < 2 ) {
        std::cerr << "\nUsage: " << argv[0] << " parameterfile \n" << "\n\t --- OR --- \n";
        std::cerr << "\nUsage: " << argv[0] << " paramfile:"<<"file" << " more-opts:val ..." << std::endl;
        std::cerr << "\nUsage: " << argv[0] << " -d paramfile "<< "\n\t(for displaying solutions in grape) "<< std::endl;
        Parameters().PrintParameterSpecs( std::cerr );
        std::cerr << std::endl;
        return 2;
    }

    if ( !(  Parameters().ReadCommandLine( argc, argv ) ) )
        DUNE_THROW( Dune::IOError, "Error parsing command line arguments" );

    // LOG_NONE = 1, LOG_ERR = 2, LOG_INFO = 4,LOG_DEBUG = 8,LOG_CONSOLE = 16,LOG_FILE = 32
    //--> LOG_ERR | LOG_INFO | LOG_DEBUG | LOG_CONSOLE | LOG_FILE = 62
    const bool useLogger = false;
    Logger().Create( Parameters().getParam( "loglevel",         62,                         useLogger ),
                     Parameters().getParam( "logfile",          std::string("dune_stokes"), useLogger ),
                     Parameters().getParam( "fem.io.datadir",   std::string("data"),        useLogger ),
                     Parameters().getParam( "fem.io.logdir",    std::string(),              useLogger )
                    );

    return CollectiveCommunication();//( Dune::MPIManager::helper().getCommunicator() );
}

#endif // MAIN_H