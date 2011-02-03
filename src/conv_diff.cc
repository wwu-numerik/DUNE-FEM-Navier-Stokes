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

#if ! defined(POLORDER)
	#define POLORDER 0
	#warning ("using default polorder 0 for all spaces")
#endif

#if ! defined(PRESSURE_POLORDER)
	#define PRESSURE_POLORDER POLORDER
#endif

#if ! defined(VELOCITY_POLORDER)
	#define VELOCITY_POLORDER POLORDER
#endif

#if ! defined(TESTCASE)
	#define TESTCASE TestCase3D
#endif

#define TESTCASE_NAME "TESTCASE"

#if ( ( defined(SGRID) || defined(ALUGRID_SIMPLEX) ||  defined(ALUGRID_CUBE) ) && ( GRIDDIM == 3 ) ) || defined(UGGRID) || defined(YASPGRID)
	//this is no mistake, ALU is indeed only incompatible in 3d
	#define OLD_DUNE_GRID_VERSION
#endif
#define TESTING_NS Testing::AdapterFunctionsVectorial
#include "testing.hh"
#define USE_GRPAE_VISUALISATION (HAVE_GRAPE && !defined( AORTA_PROBLEM ))

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
#include <dune/fem/misc/l2error.hh>

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
#include <dune/stuff/tuple.hh>
#include <dune/stuff/error.hh>
#include <dune/stuff/functionadapter.hh>


#include <dune/navier/thetascheme_traits.hh>

#include "conv_diff.hh"

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
	Stuff::Signals::installSignalHandler(SIGINT);
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
		RunInfoVectorMap rf;
		for ( unsigned int ref = minref;
			  ref <= maxref;
			  ++ref )
		{
			singleRun( mpicomm, ref );
			profiler().NextRun();
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
	catch ( std::exception& e ) {
		std::cerr << "Exception thrown at:\n" << e.what() << std::endl ;
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
	GridPartType gridPart( *gridPtr );

	/* ********************************************************************** *
	 * initialize problem                                                     *
	 * ********************************************************************** */
	infoStream << "\n- initialising problem" << std::endl;

	const int polOrder = POLORDER;
	debugStream << "  - polOrder: " << polOrder << std::endl;


	Parameters().setParam( "reduced_oseen_solver", true );

	const double theta_ = 1.0;
	const double alpha = Parameters().getParam( "alpha", 1.0 );
	const double viscosity = Parameters().getParam( "viscosity", 1.0 );

	Dune::StabilizationCoefficients stab_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();
	stab_coeff.FactorFromParams( "D12", 0 );
	stab_coeff.FactorFromParams( "C12", 0 );
	stab_coeff.Add( "E12", 0.5 );

	typedef Dune::ConvDiff::Traits<
			CollectiveCommunication,
			GridPartType,
			gridDim,
			polOrder,
			VELOCITY_POLORDER,
			PRESSURE_POLORDER >
		ConvDiffTraits;

	CollectiveCommunication comm = Dune::MPIManager::helper().getCommunicator();

	ConvDiffTraits::TimeProviderType timeprovider_( ConvDiffTraits::SchemeDescriptionType::crank_nicholson( 0.5 ), comm );
	ConvDiffTraits::OseenModelTraits::DiscreteStokesFunctionSpaceWrapperType functionSpaceWrapper ( gridPart );

	typedef ConvDiffTraits::OseenModelTraits::DiscreteStokesFunctionWrapperType
		DiscreteStokesFunctionWrapperType;
	DiscreteStokesFunctionWrapperType currentFunctions(  "current_",
						functionSpaceWrapper,
						gridPart );
	DiscreteStokesFunctionWrapperType nextFunctions(  "next_",
					functionSpaceWrapper,
					gridPart );
	DiscreteStokesFunctionWrapperType errorFunctions(  "error_",
					functionSpaceWrapper,
					gridPart );
	ConvDiffTraits::ExactSolutionType exactSolution( timeprovider_,
					gridPart,
					functionSpaceWrapper );
	exactSolution.project();

	Logging::MatlabLogStream& matlabLogStream = Logger().Matlab();
//	Stuff::printDiscreteFunctionMatlabStyle( exactSolution.discretePressure(), "p_exakt", matlabLogStream );
//	Stuff::printDiscreteFunctionMatlabStyle( exactSolution.discreteVelocity(), "u_exakt", matlabLogStream );

	ConvDiffTraits::StartPassType startPass;
	ConvDiffTraits::OseenModelTraits::AnalyticalDirichletDataType stokesDirichletData =
			ConvDiffTraits::OseenModelTraits ::AnalyticalDirichletDataTraitsImplementation
							::getInstance( timeprovider_,
										   functionSpaceWrapper );

	ConvDiffTraits::OseenModelTraits::VelocityFunctionSpaceType
			continousVelocitySpace;

	ConvDiffTraits::OseenModelTraits::AnalyticalForceFunctionType force( viscosity, continousVelocitySpace, alpha );
	ConvDiffTraits::OseenModelType
			stokesModel( stab_coeff ,
						force,
						stokesDirichletData,
						viscosity,
						alpha );
	currentFunctions.assign( exactSolution );

	ConvDiffTraits::OseenModelTraits::DiscreteSigmaFunctionSpaceType sigma_space ( gridPart );
	ConvDiffTraits::OseenModelTraits::DiscreteSigmaFunctionType discrete_velocityGradient( "velocityGradient", sigma_space );
	ConvDiffTraits::OseenModelTraits::SigmaFunctionSpaceType cont_sigma_space;
	ConvDiffTraits::VelocityGradientType velocityGradient( timeprovider_, cont_sigma_space );

	ConvDiffTraits::ConvectionType convection( timeprovider_, continousVelocitySpace );
	DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType discrete_convection( "convection", currentFunctions.discreteVelocity().space() );
	DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType discrete_exactConvection( "exact_convection", currentFunctions.discreteVelocity().space() );
	DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType convection_diff( "convection_diff", currentFunctions.discreteVelocity().space() );
	Dune::L2Projection< double,
						double,
						ConvDiffTraits::ConvectionType,
						DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType >
		()(convection, discrete_convection);

	ConvDiffTraits::OseenPassType::RhsDatacontainer rhs_container ( currentFunctions.discreteVelocity().space(),
																 sigma_space );
	ConvDiffTraits::OseenPassType oseenPass( startPass,
							stokesModel,
							gridPart,
							functionSpaceWrapper,
							discrete_convection,
							true );
	oseenPass.apply( currentFunctions, nextFunctions, &rhs_container, &velocityGradient );

	ConvDiffTraits::ExactConvectionType exactConvection( timeprovider_, continousVelocitySpace );
	Dune::L2Projection< double,
						double,
						ConvDiffTraits::ExactConvectionType,
						DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType >
		()(exactConvection, discrete_exactConvection);

	typedef Stuff::GradientSplitterFunction<	DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType,
												ConvDiffTraits::OseenModelTraits::DiscreteSigmaFunctionType >
			GradientSplitterFunctionType;
	GradientSplitterFunctionType gradient_splitter(	functionSpaceWrapper.discreteVelocitySpace(),
									rhs_container.velocity_gradient );

	ConvDiffTraits::VelocityGradientYType velocityGradientY( timeprovider_, continousVelocitySpace );
	DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType discrete_velocityGradientY( "velocityGradientY", currentFunctions.discreteVelocity().space() );
	Dune::L2Projection< double,
						double,
						ConvDiffTraits::VelocityGradientYType,
						DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType >
		()(velocityGradientY, discrete_velocityGradientY);

	ConvDiffTraits::VelocityGradientXType velocityGradientX( timeprovider_, continousVelocitySpace );
	DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType discrete_velocityGradientX( "velocityGradientX", currentFunctions.discreteVelocity().space() );
	Dune::L2Projection< double,
						double,
						ConvDiffTraits::VelocityGradientXType,
						DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType >
		()(velocityGradientX, discrete_velocityGradientX);

	Dune::BetterL2Projection
		::project(0.0,velocityGradient, discrete_velocityGradient);

	GradientSplitterFunctionType exact_gradient_splitter(	functionSpaceWrapper.discreteVelocitySpace(),
									discrete_velocityGradient );

	ConvDiffTraits::VelocityLaplaceType velocityLaplace( timeprovider_, continousVelocitySpace );
	DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType discrete_velocityLaplace( "exact_laplace", currentFunctions.discreteVelocity().space() );
	Dune::L2Projection< double,
						double,
						ConvDiffTraits::VelocityLaplaceType,
						DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType >
		()(velocityLaplace, discrete_velocityLaplace);

	typedef Stuff::L2Error<GridPartType>
			L2ErrorType;
	L2ErrorType l2Error( gridPart );
	L2ErrorType::Errors errors_convection = l2Error.get( rhs_container.convection,
															discrete_exactConvection,
															convection_diff );
	L2ErrorType::Errors errors_velocity	 = l2Error.get( nextFunctions.discreteVelocity(),
															exactSolution.discreteVelocity(),
															errorFunctions.discreteVelocity() );

	L2ErrorType::Errors errors_gradient = l2Error.get(	rhs_container.velocity_gradient,
														discrete_velocityGradient );
	L2ErrorType::Errors errors_laplace = l2Error.get(	rhs_container.velocity_laplace,
														discrete_velocityLaplace );
	double GD = Stuff::boundaryIntegral( stokesDirichletData, nextFunctions.discreteVelocity().space() );

	Logger().Info().Resume();
	Logger().Info()
//			<< "L2-Error Pressure (abs|rel): " << std::setw(8) << l2_error_pressure << " | " << relative_l2_error_pressure << "\n"
					<< errors_convection.str()
					<< errors_velocity.str()
					<< errors_gradient.str()
					<< errors_laplace.str()
//					<< "Mean pressure (exact|discrete): " << meanPressure_exact << " | " << meanPressure_discrete << std::endl
					<< "GD: " << GD << std::endl;

	typedef Stuff::FullTuple<	const DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType* >
		OutputTupleType;
	typedef Dune::TimeAwareDataWriter<	ConvDiffTraits::TimeProviderType,
										GridPartType::GridType,
										OutputTupleType >
		DataWriterType;
	OutputTupleType out( &nextFunctions.discreteVelocity(),
						 &exactSolution.discreteVelocity(),
						 &errorFunctions.discreteVelocity(),
						 &rhs_container.convection,
						 &discrete_exactConvection,
						gradient_splitter[0].get(),
	                    gradient_splitter[1].get(),
						exact_gradient_splitter[0].get(),
	                    exact_gradient_splitter[1].get()
						 );

	DataWriterType dt( timeprovider_,
					   gridPart.grid(),
					   out );
	dt.write();

	OutputTupleType out2( &discrete_velocityLaplace,
						 &rhs_container.velocity_laplace,
						 0,
						 0,
						 0,
						0,
						0,
						0,
						0
						 );

	DataWriterType dt2( timeprovider_,
					   gridPart.grid(),
					   out2 );
	dt2.write();

	return runInfoVector;
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

