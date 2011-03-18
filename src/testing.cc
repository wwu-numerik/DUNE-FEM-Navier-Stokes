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
//#define TESTING_NS Testing::AdapterFunctionsVisco
#define TESTING_NS Testing::AdapterFunctionsVectorial
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

#include <dune/navier/thetascheme.hh>
#include <dune/navier/testdata.hh>

#include <dune/stokes/discretestokesmodelinterface.hh>
#include <dune/stokes/stokespass.hh>
#include <dune/navier/fractionaltimeprovider.hh>
#include <dune/navier/stokestraits.hh>
#include <dune/navier/exactsolution.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/stuff/datawriter.hh>
#include <dune/stuff/customprojection.hh>
#include <dune/common/collectivecommunication.hh>
#include <dune/stuff/error.hh>
#include <dune/stuff/functionadapter.hh>
#include <boost/format.hpp>

#ifndef COMMIT
	#define COMMIT "undefined"
#endif

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

		const int minref = Parameters().getParam( "minref", 0 );
		// ensures maxref>=minref
		const int maxref = Stuff::clamp( Parameters().getParam( "maxref", 0 ), minref, Parameters().getParam( "maxref", 0 ) );
		profiler().Reset( maxref - minref + 1 );
		for ( int ref = minref;
			  ref <= maxref;
			  ++ref )
		{
			singleRun( mpicomm, ref );
			profiler().NextRun();
			std::cout << "\n-------------------------------------------------------------\n";
		}

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
#include "testing.hh"
RunInfoVector singleRun(  CollectiveCommunication& mpicomm,
					int refine_level_factor )
{
	profiler().StartTiming( "SingleRun" );
	Logging::LogStream& infoStream = Logger().Info();
	Logging::LogStream& debugStream = Logger().Dbg();
	RunInfoVector runInfoVector;


	/* ********************************************************************** *
	 * initialize the grid                                                    *
	 * ********************************************************************** */
	infoStream << "\n- initialising grid" << std::endl;
	const int gridDim = GridType::dimensionworld;
	Dune::GridPtr< GridType > gridPtr( Parameters().DgfFilename( gridDim ) );
	int refine_level = ( refine_level_factor  ) * Dune::DGFGridInfo< GridType >::refineStepsForHalf();
	gridPtr->globalRefine( refine_level );

	typedef Dune::AdaptiveLeafGridPart< GridType >
		GridPartType;
	GridPartType gridPart_( *gridPtr );

	/* ********************************************************************** *
	 * initialize problem                                                     *
	 * ********************************************************************** */
	infoStream << "\n- initialising problem" << std::endl;

	const int polOrder = POLORDER;
	debugStream << "  - polOrder: " << polOrder << std::endl;
	typedef Stuff::L2Error<GridPartType>
			L2ErrorType;
	L2ErrorType l2Error( gridPart_ );
	// model traits
	typedef Dune::NavierStokes::ThetaSchemeTraits<
					CollectiveCommunication,
					GridPartType,
					TESTING_NS::Force,
					TESTING_NS::DirichletData,
					TESTING_NS::Pressure,
					TESTING_NS::Velocity,
					gridDim,
					polOrder,
					VELOCITY_POLORDER,
					PRESSURE_POLORDER >
		Traits;
	typedef Dune::NavierStokes::ThetaScheme<Traits>
		ThetaSchemeType;
	Traits::CommunicatorType communicator_= Dune::MPIManager::helper().getCommunicator();

	Traits::DiscreteStokesFunctionSpaceWrapperType functionSpaceWrapper_( gridPart_ );

	Traits::StokesModelTraits::VelocityFunctionSpaceType
			continousVelocitySpace_;

		typedef Traits::DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType
			DiscreteVelocityFunctionType;
//	DiscreteVelocityFunctionType velocity_convection_discrete("velocity_convection_discrete", functionSpaceWrapper_.discreteVelocitySpace() );
//	DiscreteVelocityFunctionType velocity_laplace_discrete("velocity_laplace_discrete", exactSolution_.discreteVelocity().space() );
//	DiscreteVelocityFunctionType pressure_gradient_discrete("pressure_gradient_discrete", exactSolution_.discreteVelocity().space() );

//	Dune::BetterL2Projection
//		::project(timeprovider_,velocity_laplace, velocity_laplace_discrete );

//	L2ErrorType::Errors errors_pressure_gradient = l2Error.get(	thetaScheme.rhsDatacontainer().pressure_gradient,
//														pressure_gradient_discrete );
//	L2ErrorType::Errors errors_velocity_gradient = l2Error.get(	thetaScheme.rhsDatacontainer().velocity_gradient,
//														velocity_gradient_discrete );
//	L2ErrorType::Errors errors_velocity_laplace = l2Error.get(	thetaScheme.rhsDatacontainer().velocity_laplace,
//														velocity_laplace_discrete );
//	std::cout << errors_pressure_gradient.str()
//			  << errors_velocity_laplace.str()
//			  << errors_velocity_gradient.str();

//	typedef Dune::Tuple<	const DiscreteVelocityFunctionType*,
//							const DiscreteVelocityFunctionType*,
//							const DiscreteVelocityFunctionType*,
//							const DiscreteVelocityFunctionType*,
//							const DiscreteVelocityFunctionType*,
//							const DiscreteVelocityFunctionType*,
//							const DiscreteVelocityFunctionType*,
//							const DiscreteVelocityFunctionType*,
//							const DiscreteVelocityFunctionType*
//						>
//		OutputTupleType;
//	typedef Dune::TimeAwareDataWriter<	Traits::TimeProviderType,
//										GridPartType::GridType,
//										OutputTupleType >
//		DataWriterType;

//	{
//		GradientSplitterFunctionType velocity_gradient_splitter(	exactSolution_.discreteVelocity().space(),
//										thetaScheme.rhsDatacontainer().velocity_gradient );
//		OutputTupleType out(
//							 &pressure_gradient_discrete,
//							 &(thetaScheme.rhsDatacontainer().velocity_laplace),
//							 &velocity_laplace_discrete,
//							&(thetaScheme.rhsDatacontainer().pressure_gradient),
//							 &pressure_gradient_discrete,
//							velocity_gradient_splitter[0].get(),
//							velocity_gradient_splitter[1].get(),
//							exact_gradient_splitter[0].get(),
//							exact_gradient_splitter[1].get()
//							);
//		DataWriterType( timeprovider_,
//						   gridPart_.grid(),
//						   out ).write();
//	}

//	RunInfo info_dummy;
//	thetaScheme.nextStep(1,info_dummy);

	return runInfoVector;
}


