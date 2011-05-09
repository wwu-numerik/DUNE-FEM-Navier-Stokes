/**
 *  \file   dune_stokes.cc
 *
 *  \brief  brief
 **/

#include <dune/navier/global_defines.hh>

#include <cstdio>
#if defined(USE_PARDG_ODE_SOLVER) && defined(USE_BFG_CG_SCHEME)
	#warning ("USE_PARDG_ODE_SOLVER enabled, might conflict with custom solvers")
#endif

#if defined(UGGRID) && defined(DEBUG)
	#warning ("UGGRID in debug mode is likely to produce a segfault")
#endif

#define USE_GRPAE_VISUALISATION (HAVE_GRAPE && !defined( AORTA_PROBLEM ))

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
#include <dune/stuff/signals.hh>

#include "oseen.hh"
#include <dune/stuff/runinfo.hh>

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
Stuff::RunInfoVector singleRun(  CollectiveCommunication& mpicomm,
					int refine_level_factor  );

//! output alert for neg. EOC
void eocCheck( const Stuff::RunInfoVector& runInfos );

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
						 Parameters().getParam( "fem.io.datadir",	std::string("data"),		useLogger ),
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

Stuff::RunInfoVector singleRun(  CollectiveCommunication& /*mpicomm*/,
					int refine_level_factor )
{
	profiler().StartTiming( "SingleRun" );
	Logging::LogStream& infoStream = Logger().Info();
	Logging::LogStream& debugStream = Logger().Dbg();
	Stuff::RunInfoVector runInfoVector;


	/* ********************************************************************** *
	 * initialize the grid                                                    *
	 * ********************************************************************** */
	infoStream << "\n- initialising grid" << std::endl;
	const int gridDim = GridType::dimensionworld;
	Dune::GridPtr< GridType > gridPtr( Parameters().DgfFilename( gridDim ) );
	static bool firstRun = true;
	int refine_level = ( refine_level_factor  ) * Dune::DGFGridInfo< GridType >::refineStepsForHalf();
	if ( firstRun && refine_level_factor > 0 ) {
		refine_level = ( refine_level_factor ) * Dune::DGFGridInfo< GridType >::refineStepsForHalf();
		gridPtr->globalRefine( refine_level );
	}

	typedef Dune::AdaptiveLeafGridPart< GridType >
		GridPartType;
	GridPartType gridPart( *gridPtr );

	infoStream << "\n- initialising problem" << std::endl;

	const int polOrder = POLORDER;
	debugStream << "  - polOrder: " << polOrder << std::endl;

	const double grid_width = Dune::GridWidth::calcGridWidth( gridPart );
	infoStream << "  - max grid width: " << grid_width << std::endl;

//	Dune::CompileTimeChecker< ( VELOCITY_POLORDER >= 2 ) > RHS_ADAPTER_CRAPS_OUT_WITH_VELOCITY_POLORDER_LESS_THAN_2;

//	const double reynolds = Parameters().getParam( "reynolds", 1.0 );
//	const double theta_ = 1.0;
//	const double d_t = 1.0;
//	const double operator_weight_beta_ = 1.0;
//	const double operator_weight_alpha_ = 1.0;
	const double oseen_alpha = Parameters().getParam( "alpha", 1.0 );
	const double oseen_viscosity = Parameters().getParam( "viscosity", 1.0 );
//	const double lambda = ( reynolds * 0.5 )
//						  - std::sqrt(
//								  ( std::pow( reynolds, 2 ) * 0.25 )
//								  + ( 4 * std::pow( M_PI, 2 ) )
//									  ) ;
//	const double pressure_C = ( std::exp( 3 * lambda ) - std::exp(-1  * lambda ) ) / ( - 8 * lambda );

//	const double lambda = - 8 *M_PI * M_PI / ( reynolds + std::sqrt(reynolds*reynolds + 64 * M_PI * M_PI));

//	Parameters().setParam( "lambda", lambda );
	Parameters().setParam( "viscosity", oseen_viscosity );
	Dune::StabilizationCoefficients stab_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();
//	stab_coeff.FactorFromParams( "D12", 0 );
//	stab_coeff.FactorFromParams( "C12", 0 );
//	stab_coeff.Add( "E12", 0.5 );
//	stab_coeff.FactorFromParams( "E12", 0.5 );

	typedef Dune::Oseen::Traits<
			CollectiveCommunication,
			GridPartType,
			gridDim,
			polOrder,
			VELOCITY_POLORDER,
			PRESSURE_POLORDER >
		OseenTraits;

	CollectiveCommunication comm = Dune::MPIManager::helper().getCommunicator();
	OseenTraits::TimeProviderType timeprovider_( OseenTraits::SchemeDescriptionType::crank_nicholson( 0.5 ), comm );
	OseenTraits::OseenModelTraits::DiscreteStokesFunctionSpaceWrapperType functionSpaceWrapper ( gridPart );

	typedef OseenTraits::OseenModelTraits::DiscreteStokesFunctionWrapperType
		DiscreteStokesFunctionWrapperType;
	DiscreteStokesFunctionWrapperType currentFunctions(  "current_",
						functionSpaceWrapper,
						gridPart );
	DiscreteStokesFunctionWrapperType nextFunctions(  "next_",
					functionSpaceWrapper,
					gridPart );
	DiscreteStokesFunctionWrapperType tmpFunctions(  "tmp_",
					functionSpaceWrapper,
					gridPart );
	DiscreteStokesFunctionWrapperType errorFunctions(  "error_",
					functionSpaceWrapper,
					gridPart );
	OseenTraits::ExactSolutionType exactSolution( timeprovider_,
					gridPart,
					functionSpaceWrapper );
	exactSolution.project();
//	exactSolution.exactPressure().setShift( pressure_C );
	Logging::MatlabLogStream& matlabLogStream = Logger().Matlab();
	Stuff::printDiscreteFunctionMatlabStyle( exactSolution.discretePressure(), "p_exakt", matlabLogStream );
	Stuff::printDiscreteFunctionMatlabStyle( exactSolution.discreteVelocity(), "u_exakt", matlabLogStream );

	OseenTraits::StartPassType startPass;
	OseenTraits::OseenModelTraits::AnalyticalDirichletDataType stokesDirichletData( timeprovider_,
										   functionSpaceWrapper );

//	OseenTraits::OseenModelTraits::PressureFunctionSpaceType
//			continousPressureSpace;
	OseenTraits::OseenModelTraits::VelocityFunctionSpaceType
			continousVelocitySpace;

	OseenTraits::OseenModelTraits::AnalyticalForceFunctionType force( timeprovider_, continousVelocitySpace, oseen_viscosity, oseen_alpha );
	OseenTraits::OseenModelType
			stokesModel( stab_coeff ,
						force,
						stokesDirichletData,
						oseen_viscosity, /*viscosity*/
						oseen_alpha, /*alpha*/
						1,//Parameters().getParam( "cscale", 1.0 ),/*convection_scale_factor*/
						1/*pressure_gradient_scale_factor*/);

//	currentFunctions.assign( exactSolution );
	currentFunctions.clear();
	nextFunctions.clear();

	OseenTraits::ConvectionType convection( timeprovider_, continousVelocitySpace );
	DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType discrete_convection( "convetion", currentFunctions.discreteVelocity().space() );
	Dune::L2Projection< double,
						double,
						OseenTraits::ConvectionType,
						DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType >
		()(convection, discrete_convection);

	OseenTraits::OseenModelTraits::DiscreteSigmaFunctionSpaceType sigma_space ( gridPart );
	OseenTraits::OseenPassType::RhsDatacontainer rhs_container ( currentFunctions.discreteVelocity().space(),
																 sigma_space );
	OseenTraits::OseenPassType oseenPass( startPass,
							stokesModel,
							gridPart,
							functionSpaceWrapper,
							exactSolution.discreteVelocity(),
							true );
	nextFunctions.clear();
	oseenPass.printInfo();
	oseenPass.apply( currentFunctions, nextFunctions, &rhs_container, &convection );

	errorFunctions.discretePressure().assign( exactSolution.discretePressure() );
	errorFunctions.discretePressure() -= nextFunctions.discretePressure();
	errorFunctions.discreteVelocity().assign( exactSolution.discreteVelocity() );
	errorFunctions.discreteVelocity() -= nextFunctions.discreteVelocity();

	double meanPressure_exact = Stuff::integralAndVolume( exactSolution.exactPressure(), nextFunctions.discretePressure().space() ).first;
	double meanPressure_discrete = Stuff::meanValue( nextFunctions.discretePressure(), nextFunctions.discretePressure().space() );
//	typedef OseenTraits::OseenModelTraits::PressureFunctionSpaceType
//			PressureFunctionSpaceType;
//	PressureFunctionSpaceType pressureFunctionSpace;
//	Stuff::ConstantFunction<PressureFunctionSpaceType> vol(pressureFunctionSpace, meanPressure_discrete );
//	Dune::BetterL2Projection
//		::project( 0.0, vol, tmpFunctions.discretePressure() );
//	nextFunctions.discretePressure() -= tmpFunctions.discretePressure();
	double meanPressure_discrete_after = Stuff::meanValue( nextFunctions.discretePressure(), nextFunctions.discretePressure().space() );

	double GD = Stuff::boundaryIntegral( stokesDirichletData, nextFunctions.discreteVelocity().space() );

	Dune::L2Norm< GridPartType > l2_Error( gridPart );

	const double l2_error_pressure				= l2_Error.norm( errorFunctions.discretePressure() );
	const double l2_error_velocity				= l2_Error.norm( errorFunctions.discreteVelocity() );
	const double relative_l2_error_pressure		= l2_error_pressure / l2_Error.norm( exactSolution.discretePressure() );
	const double relative_l2_error_velocity		= l2_error_velocity / l2_Error.norm( exactSolution.discreteVelocity() );

//	errorFunctions.discretePressure().assign( exactSolution.discretePressure() );
//	errorFunctions.discretePressure() -= nextFunctions.discretePressure();

	const double l2_error_pressure_after				= l2_Error.norm( errorFunctions.discretePressure() );
	const double relative_l2_error_pressure_after		= l2_error_pressure_after / l2_Error.norm( exactSolution.discretePressure() );


	Logger().Info().Resume();
	Logger().Info() << "L2-Error Pressure (abs|rel): " << std::setw(8) << l2_error_pressure << " | " << relative_l2_error_pressure << "\n"
					<< "L2-Error Pressure after (abs|rel): " << std::setw(8) << l2_error_pressure_after << " | " << relative_l2_error_pressure_after << "\n"
					<< "L2-Error Velocity (abs|rel): " << std::setw(8) << l2_error_velocity << " | " << relative_l2_error_velocity << "\n"
					<< "Mean pressure (exact|discrete|after): " << meanPressure_exact << " | "
																<< meanPressure_discrete << " | "
																<< meanPressure_discrete_after << std::endl
					<< "GD: " << GD << "\n"
//					<< "lambda: " << lambda
					<< "current time: " << timeprovider_.time()
					<< std::endl;

	typedef Dune::Tuple<	const DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType*,
							const DiscreteStokesFunctionWrapperType::DiscretePressureFunctionType*,
							const DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType*,
							const DiscreteStokesFunctionWrapperType::DiscretePressureFunctionType*,
							const DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType*,
							const DiscreteStokesFunctionWrapperType::DiscretePressureFunctionType*,
							const DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType*
						>
		OutputTupleType;
	typedef Dune::TimeAwareDataWriter<	OseenTraits::TimeProviderType,
										GridPartType::GridType,
										OutputTupleType >
		DataWriterType;
	OutputTupleType out( &nextFunctions.discreteVelocity(),
						 &nextFunctions.discretePressure(),
						 &exactSolution.discreteVelocity(),
						 &exactSolution.discretePressure(),
						 &errorFunctions.discreteVelocity(),
						 &errorFunctions.discretePressure(),
						 &discrete_convection);

	DataWriterType dt( timeprovider_,
					   gridPart.grid(),
					   out );
	dt.write();

	return runInfoVector;
}

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

