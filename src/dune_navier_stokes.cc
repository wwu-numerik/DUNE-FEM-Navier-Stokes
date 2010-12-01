/**
 *  \file   dune_stokes.cc
 *
 *  \brief  brief
 **/

#include "cmake_config.h"

#include <cstdio>
#if defined(USE_PARDG_ODE_SOLVER) && defined(USE_BFG_CG_SCHEME)
	#warning ("USE_PARDG_ODE_SOLVER enabled, might conflict with custom solvers")
#endif

//the adaption manager might be troublesome with certain gridparts/spaces, so we needed a easy way to disable it
#ifndef ENABLE_ADAPTIVE
	#define ENABLE_ADAPTIVE 1
#endif

#if defined(UGGRID) && defined(DEBUG)
	#warning ("UGGRID in debug mode is likely to produce a segfault")
#endif

#if ! defined(TESTCASE)
	#define TESTCASE TestCase3D
#endif

#define TESTCASE_NAME "TESTCASE"

#if ( ( defined(SGRID) || defined(ALUGRID_SIMPLEX) ||  defined(ALUGRID_CUBE) ) && ( GRIDDIM == 3 ) ) || defined(UGGRID) || defined(YASPGRID)
	//this is no mistake, ALU is indeed only incompatible in 3d
	#define OLD_DUNE_GRID_VERSION
#endif

#define USE_GRPAE_VISUALISATION (HAVE_GRAPE && !defined( AORTA_PROBLEM ))

#if (GRIDDIM==3)
	#define MODEL_PROVIDES_LOCALFUNCTION 1
#endif

//#define NS Dune::NavierStokes::TESTCASE
//#define NS Testing::AdapterFunctionsVisco
#define NS Testing::AdapterFunctionsVectorial
#define TESTING_NS NS
#include "testing.hh"

#include <vector>
#include <string>

#include <iostream>
#include <cmath>
#include <dune/fem/misc/mpimanager.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/grid/common/capabilities.hh>

//!ATTENTION: undef's GRIDDIM
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh> // for the grid

#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/space/dgspace/dgadaptiveleafgridpart.hh>
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

#include <dune/navier/thetascheme.hh>
#include <dune/navier/testdata.hh>
#include "testing.hh"

#ifndef COMMIT
	#define COMMIT "undefined"
#endif

#define MODEL_PROVIDES_LOCALFUNCTION 1

static const std::string commit_string (COMMIT);

#if ENABLE_MPI
		typedef Dune::CollectiveCommunication< MPI_Comm > CollectiveCommunication;
#else
		typedef Dune::CollectiveCommunication< double > CollectiveCommunication;
#endif

//! the strings used for column headers in tex output
typedef std::vector<std::string>
	ColumnHeaders;


/** \brief one single application of the discretisation and solver

	\param  mpicomm
			mostly useless atm, but mandatory
	\param  refine_level_factor
			integer to be multiplied by Dune::DGFGridInfo< GridType >::refineStepsForHalf()
			to get the used refine level for the constructed grid
	\param  stabil_coeff
			the set of coefficients to be used in the run. Default is used in all run types but StabRun().

**/
RunInfoVector singleRun(  CollectiveCommunication& mpicomm,
					int refine_level_factor  );

//! output alert for neg. EOC
void eocCheck( const RunInfoVector& runInfos );

/**
 *  \brief  main function
 *
 *  ParameterContainer and Logger setup, select run type from parameter file and go
 *
 *  \param  argc
 *          number of arguments from command line
 *  \param  argv
 *          array of arguments from command line
 **/
int main( int argc, char** argv )
{
	Stuff::Signals::installSignalHandler();
#ifdef NDEBUG
	try
#endif
	{

		Dune::MPIManager::initialize(argc, argv);
		//assert( Dune::Capabilities::isParallel< GridType >::v );
		CollectiveCommunication mpicomm( Dune::MPIManager::helper().getCommunicator() );

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
		#if USE_GRPAE_VISUALISATION
		if ( !strcmp( argv[1], "-d" ) || !strcmp( argv[1], "-r" ) ) {
			return display( argc, argv );
		}
		#endif
		if ( !(  Parameters().ReadCommandLine( argc, argv ) ) ) {
			return 1;
		}

		// LOG_NONE = 1, LOG_ERR = 2, LOG_INFO = 4,LOG_DEBUG = 8,LOG_CONSOLE = 16,LOG_FILE = 32
		//--> LOG_ERR | LOG_INFO | LOG_DEBUG | LOG_CONSOLE | LOG_FILE = 62
		const bool useLogger = false;
		Logger().Create( Parameters().getParam( "loglevel",         62,                         useLogger ),
						 Parameters().getParam( "logfile",          std::string("dune_stokes"), useLogger ),
						 Parameters().getParam( "fem.io.logdir",    std::string(),              useLogger )
						);

		int err = 0;
		const unsigned int minref = Parameters().getParam( "minref", 0 );
		RunInfoVectorMap rf;
		if ( Parameters().getParam( "runtype", 5 ) == 6 ) {
			Logger().Info() << "Time refine runs\n";
			const int dt_steps = Parameters().getParam( "dt_steps", 3 );
			profiler().Reset( dt_steps - 1 );
			int current_step = 0;
			for ( double dt = Parameters().getParam( "fem.timeprovider.dt", 0.1 );
				  dt_steps > current_step;
				  ++current_step )
			{
				rf[current_step] = singleRun( mpicomm, minref );
				assert( rf.size() );
				rf[current_step].at(0).refine_level = minref;//just in case the key changes from ref to sth else
				profiler().NextRun();
				dt /= 2.0;
				Parameters().setParam( "fem.timeprovider.dt", dt );
			}

		}
		else {
			// ensures maxref>=minref
			const unsigned int maxref = Stuff::clamp( Parameters().getParam( "maxref", (unsigned int)(0) ), minref, Parameters().getParam( "maxref", (unsigned int)(0) ) );
			profiler().Reset( maxref - minref + 1 );
			Logger().Info() << "Grid refine runs\n";
			for ( unsigned int ref = minref;
				  ref <= maxref;
				  ++ref )
			{
				rf[ref] = singleRun( mpicomm, ref );
				rf[ref].at(0).refine_level = ref;//just in case the key changes from ref to sth else
				profiler().NextRun();
			}
		}
//			profiler().OutputMap( mpicomm, rf );//! \TODO find out why this ain't working in refine runs

		Stuff::TimeSeriesOutput out( rf );
		out.writeTex( Parameters().getParam("fem.io.datadir", std::string(".") ) + std::string("/timeseries") );

		Logger().Dbg() << "\nRun from: " << commit_string << std::endl;
		return err;
	}


#ifdef NDEBUG
  catch (Dune::Exception &e){
	std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch ( std::bad_alloc& b ) {
	  std::cerr << "Memory allocation failed: " << b.what() ;
	  Logger().Info().Resume();
	  Stuff::meminfo( Logger().Info() );
  }
  catch ( assert_exception& a ) {
	  std::cerr << "Exception thrown at:\n" << a.what() << std::endl ;
  }
  catch (...){
	std::cerr << "Unknown exception thrown!" << std::endl;
  }
#endif
}

RunInfoVector singleRun(  CollectiveCommunication& mpicomm,
					int refine_level_factor )
{
	profiler().StartTiming( "SingleRun" );
	Logging::LogStream& infoStream = Logger().Info();
	Logging::LogStream& debugStream = Logger().Dbg();

	infoStream << "\n- initialising grid" << std::endl;
	const int gridDim = GridType::dimensionworld;
	Dune::GridPtr< GridType > gridPtr( Parameters().DgfFilename( gridDim ) );
	int refine_level = ( refine_level_factor  ) * Dune::DGFGridInfo< GridType >::refineStepsForHalf();
	gridPtr->globalRefine( refine_level );
	typedef Dune::AdaptiveLeafGridPart< GridType >
		GridPartType;
	GridPartType gridPart( *gridPtr );

	const int polOrder = POLORDER;
	debugStream << "  - polOrder: " << polOrder << std::endl;
	const double grid_width = Dune::GridWidth::calcGridWidth( gridPart );
	infoStream << "  - max grid width: " << grid_width << std::endl;

	// model traits
	typedef Dune::NavierStokes::ThetaSchemeTraits<
					CollectiveCommunication,
					GridPartType,
					NS::Force,
					NS::DirichletData,
					NS::Pressure,
					NS::Velocity,
					gridDim,
					polOrder,
					VELOCITY_POLORDER,
					PRESSURE_POLORDER >
		ThetaSchemeTraitsType;

	Dune::NavierStokes::ThetaScheme<ThetaSchemeTraitsType>
			thetaScheme( gridPart );/*,
						 1 - std::pow( 2.0, -1/2.0 ),
						 mpicomm,
						 1,
						 0);*/
	return thetaScheme.run();
}

void eocCheck( const RunInfoVector& runInfos )
{
	bool ups = false;
	RunInfoVector::const_iterator it = runInfos.begin();
	RunInfo last = *it;
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
