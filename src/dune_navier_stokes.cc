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

#if (GRIDDIM==3)
	#define MODEL_PROVIDES_LOCALFUNCTION 1
#endif

#define MODEL_PROVIDES_LOCALFUNCTION 1

#define NS Dune::NavierStokes::TESTCASE
//#define NS Testing::AdapterFunctionsVisco
//#define NS Testing::AdapterFunctionsVectorial
//#define NS Testing::AdapterFunctionsScalar
#define TESTING_NS NS
#include "testing.hh"
#include <dune/navier/testdata.hh>

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
RunInfoTimeMap singleRun(	CollectiveCommunication& mpicomm,
							const int refine_level_factor,
							const int scheme_type );
//! output alert for neg. EOC
//void eocCheck( const RunInfoVector& runInfos );

bool setSchemeTypeFromString();

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

	if ( !(  Parameters().ReadCommandLine( argc, argv ) ) ) {
		return 1;
	}

	// LOG_NONE = 1, LOG_ERR = 2, LOG_INFO = 4,LOG_DEBUG = 8,LOG_CONSOLE = 16,LOG_FILE = 32
	//--> LOG_ERR | LOG_INFO | LOG_DEBUG | LOG_CONSOLE | LOG_FILE = 62
	const bool useLogger = false;
	Logger().Create( Parameters().getParam( "loglevel",         62,                         useLogger ),
					 Parameters().getParam( "logfile",          std::string("dune_stokes"), useLogger ),
					 Parameters().getParam( "fem.io.datadir",   std::string(),              useLogger )
					);
	if ( setSchemeTypeFromString() )
		Logger().Info() << "overrode scheme id from string" << std::endl;

	int err = 0;
	const unsigned int minref = Parameters().getParam( "minref", 0, Dune::ValidateNotLess<int>(0) );
	RunInfoTimeMapMap rf;
	const int runtype = Parameters().getParam( "runtype", 5 );
	switch( runtype ) {
		case 8: {
			Logger().Info() << "Reynolds runs\n";
			const int dt_steps = Parameters().getParam( "dt_steps", 3, Dune::ValidateNotLess<int>(2) );
			profiler().Reset( dt_steps - 1 );
			int current_step = 0;
			Stuff::LoopTimer<int,Logging::LogStream> loop_timer( current_step, dt_steps, Logger().Info() );
			for ( double viscosity = Parameters().getParam( "viscosity", 0.1, Dune::ValidateNotLess<double>(0.0) );
				  dt_steps > current_step;
				  ++loop_timer )
			{
				rf[current_step] = singleRun( mpicomm, minref, Parameters().getParam( "scheme_type", 1, true ) );
				assert( rf.size() );
				rf[current_step].begin()->second.refine_level = minref;//just in case the key changes from ref to sth else
				profiler().NextRun();
				viscosity /= 10.0f;
				Parameters().setParam( "viscosity", viscosity );
			}
			break;
		}
		case 6: {
			Logger().Info() << "Time refine runs\n";
			const int dt_steps = Parameters().getParam( "dt_steps", 3, Dune::ValidateNotLess<int>(2) );
			profiler().Reset( dt_steps - 1 );
			int current_step = 0;
			Stuff::LoopTimer<int,Logging::LogStream,Stuff::QuadraticWeights> loop_timer( current_step, dt_steps, Logger().Info() );
			for ( double dt = Parameters().getParam( "fem.timeprovider.dt", 0.1, Dune::ValidateNotLess<double>(0.0) );
				  dt_steps > current_step;
				  ++loop_timer )
			{
				rf[current_step] = singleRun( mpicomm, minref, Parameters().getParam( "scheme_type", 1, true ) );
				assert( rf.size() );
				rf[current_step].begin()->second.refine_level = minref;//just in case the key changes from ref to sth else
				profiler().NextRun();
				dt /= 2.0f;
				Parameters().setParam( "fem.timeprovider.dt", dt );
			}
			break;
		}
		case 7: {
			Logger().Info() << "Scheme runs\n";
			profiler().Reset( 4 );
			int current_scheme = 2;
			Stuff::LoopTimer<int,Logging::LogStream> loop_timer( current_scheme, 5, Logger().Info() );
			for ( ;
				  current_scheme < 6;
				  ++loop_timer )
			{
				rf[current_scheme] = singleRun( mpicomm, minref, current_scheme );
				assert( rf.size() );
				rf[current_scheme].begin()->second.refine_level = minref;//just in case the key changes from ref to sth else
				profiler().NextRun();
			}
			break;
		}
		case 5:
			Parameters().setParam( "maxref", minref );//only one run with ref=minref
		case 0:
		default: {
			// ensures maxref>=minref
			const unsigned int maxref = Stuff::clamp( Parameters().getParam( "maxref", (unsigned int)(0) ), minref, Parameters().getParam( "maxref", (unsigned int)(0) ) );
			profiler().Reset( maxref - minref + 1 );
			Logger().Info() << "Grid refine runs\n";
			unsigned int ref = minref;
			Stuff::LoopTimer<unsigned int,Logging::LogStream,Stuff::LinearWeights> loop_timer( ref, maxref - minref + 1, Logger().Info() );
			for ( ;
				  ref <= maxref;
				  ++loop_timer )
			{
				rf[ref] = singleRun( mpicomm, ref, Parameters().getParam( "scheme_type", 1, true ) );
				rf[ref].begin()->second.refine_level = ref;//just in case the key changes from ref to sth else
				profiler().NextRun();
			}
			break;
		}
	}
	profiler().OutputMap( mpicomm, rf );

	Stuff::TimeSeriesOutput out( rf );
	out.writeTex( Parameters().getParam("fem.io.datadir", std::string(".") ) + std::string("/timeseries") );

	Logger().Dbg() << "\nRun from: " << commit_string << std::endl;
	return err;
}

template < class GridPartType, class CollectiveCommunicationType >
class ThetaschemeRunner {
	private:
		typedef Dune::NavierStokes::ThetaSchemeTraits<
						CollectiveCommunicationType,
						GridPartType,
						NS::Force,
						NS::DirichletData,
						NS::Pressure,
						NS::Velocity,
						1,//number of substeps
						GridType::dimensionworld,
						POLORDER,
						VELOCITY_POLORDER,
						PRESSURE_POLORDER >
			OneStepThetaSchemeTraitsType;
		typedef Dune::NavierStokes::ThetaScheme<OneStepThetaSchemeTraitsType>
			OneStepThetaSchemeType;
		typedef typename OneStepThetaSchemeTraitsType::ThetaSchemeDescriptionType
			OneStepThetaSchemeDescriptionType;
		typedef Dune::NavierStokes::ThetaSchemeTraits<
						CollectiveCommunicationType,
						GridPartType,
						NS::Force,
						NS::DirichletData,
						NS::Pressure,
						NS::Velocity,
						3,//number of substeps
						GridType::dimensionworld,
						POLORDER,
						VELOCITY_POLORDER,
						PRESSURE_POLORDER >
			ThreeStepThetaSchemeTraitsType;
		typedef Dune::NavierStokes::ThetaScheme<ThreeStepThetaSchemeTraitsType>
			ThreeStepThetaSchemeType;
		typedef typename ThreeStepThetaSchemeTraitsType::ThetaSchemeDescriptionType
			ThreeStepThetaSchemeDescriptionType;
	public:
		ThetaschemeRunner( const GridPartType& grid_part, CollectiveCommunicationType& comm )
			:grid_part_(grid_part),
			comm_(comm)
		{}

		RunInfoTimeMap run(const int scheme_type)
		{
			const double dt_ = Parameters().getParam( "fem.timeprovider.dt", double(0.1) );
			switch ( scheme_type ) {
				case 1: return OneStepThetaSchemeType(grid_part_,
													  OneStepThetaSchemeDescriptionType::forward_euler( dt_ ) )
									.run();
				case 2: return OneStepThetaSchemeType(grid_part_,
													  OneStepThetaSchemeDescriptionType::backward_euler( dt_ ) )
									.run();
				case 3: return OneStepThetaSchemeType(grid_part_,
													  OneStepThetaSchemeDescriptionType::crank_nicholson( dt_ ) )
									.run();
				case 4: return ThreeStepThetaSchemeType(grid_part_,
													  ThreeStepThetaSchemeDescriptionType::fs0( dt_ ) )
									.run();
				default: Logger().Info() << "Using default value for theta scheme type\n";
				case 5: return ThreeStepThetaSchemeType(grid_part_,
													  ThreeStepThetaSchemeDescriptionType::fs1( dt_ ) )
									.run();
			}
		}

	private:
		const GridPartType& grid_part_;
		CollectiveCommunicationType& comm_;
};

RunInfoTimeMap singleRun(  CollectiveCommunication& mpicomm,
					const int refine_level_factor, const int scheme_type )
{
	Profiler::ScopedTiming pf_t( "SingleRun" );
	Logging::LogStream& infoStream = Logger().Info();
	Logging::LogStream& debugStream = Logger().Dbg();

	infoStream << "\n- initialising grid" << std::endl;
	const int gridDim = GridType::dimensionworld;
	Dune::GridPtr< GridType > gridPtr( Parameters().DgfFilename( gridDim ) );
	const int refine_level = ( refine_level_factor  ) * Dune::DGFGridInfo< GridType >::refineStepsForHalf();
	gridPtr->globalRefine( refine_level );
	typedef Dune::AdaptiveLeafGridPart< GridType >
		GridPartType;
	GridPartType gridPart( *gridPtr );

	const int polOrder = POLORDER;
	debugStream << "  - polOrder: " << polOrder << std::endl;
	const double grid_width = Dune::GridWidth::calcGridWidth( gridPart );
	infoStream << "  - max grid width: " << grid_width << std::endl;

#ifdef NDEBUG
	try {
#endif
		return ThetaschemeRunner<GridPartType,CollectiveCommunication>(gridPart,mpicomm).run( scheme_type );
#ifdef NDEBUG
	}
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
	return RunInfoTimeMap();
#endif
}

bool setSchemeTypeFromString()
{
	bool changed = false;
	if ( Dune::Parameter::exists("scheme_type_string") )
	{
		const std::vector<std::string>& scheme_names = Dune::NavierStokes::ThetaSchemeDescription<0>::scheme_names;
		Stuff::ValidateInList<std::string> validator( scheme_names );
		const std::string scheme_string = Parameters().getParam("scheme_type_string", scheme_names.back(), validator );
		const int scheme_id = Stuff::getIdx( scheme_names, scheme_string );
		Parameters().setParam( "scheme_type", scheme_id );
		changed = true;
	}
	return changed;
}

//void eocCheck( const RunInfoVector& runInfos )
//{
//	bool ups = false;
//	RunInfoVector::const_iterator it = runInfos.begin();
//	RunInfo last = *it;
//	++it;
//	for ( ; it != runInfos.end(); ++it ) {
//		ups = ( last.L2Errors[0] < it->L2Errors[0]
//			|| last.L2Errors[1] < it->L2Errors[1] );
//		last = *it;
//	}
//	if ( ups ) {
//		Logger().Err() 	<< 	"----------------------------------------------------------\n"
//						<<	"-                                                        -\n"
//						<<	"-                  negative EOC                          -\n"
//						<<	"-                                                        -\n"
//						<< 	"----------------------------------------------------------\n"
//						<< std::endl;
//	}
//}
