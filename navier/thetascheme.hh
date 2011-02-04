#ifndef METADATA_HH
#define METADATA_HH

#include <dune/navier/exactsolution.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/stuff/datawriter.hh>
#include <dune/stuff/tuple.hh>
#include <dune/stuff/customprojection.hh>
#include <dune/stuff/error.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/profiler.hh>
#include <dune/common/collectivecommunication.hh>
#include <cmath>
#include <boost/scoped_ptr.hpp>
#include <algorithm>

#include <dune/navier/thetascheme_traits.hh>

namespace Dune {
	namespace NavierStokes {
		template < class TraitsImp >
		class ThetaScheme {
			protected:
				typedef TraitsImp
					Traits;
				typedef typename Traits::CommunicatorType
					CommunicatorType;
				typedef typename Traits::ExactSolutionType
					ExactSolutionType;
				typedef Stuff::TupleSerializer<	typename Traits::DiscreteStokesFunctionWrapperType,
											typename Traits::DiscreteStokesFunctionWrapperType,
											ExactSolutionType,
											typename Traits::DiscreteStokesFunctionWrapperType>
					TupleSerializerType1;
				typedef typename TupleSerializerType1::TupleType
					OutputTupleType1;
				typedef TimeAwareDataWriter<	typename Traits::TimeProviderType,
												typename Traits::GridPartType::GridType,
												OutputTupleType1 >
					DataWriterType1;
				typedef CheckPointer< typename Traits::GridPartType::GridType,
									  OutputTupleType1 >
					CheckPointerType;
				typedef Stuff::TupleSerializer<	typename Traits::DiscreteStokesFunctionWrapperType >
					TupleSerializerType2;
				typedef typename TupleSerializerType2::TupleType
					OutputTupleType2;
				typedef TimeAwareDataWriter<	typename Traits::TimeProviderType,
												typename Traits::GridPartType::GridType,
												OutputTupleType2 >
					DataWriterType2;
				typedef typename Traits::DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType
					DiscreteVelocityFunctionType;
				typedef typename Traits::DiscreteStokesFunctionWrapperType::DiscretePressureFunctionType
					DiscretePressureFunctionType;

				mutable typename Traits::GridPartType gridPart_;
				const typename Traits::ThetaSchemeDescriptionType& scheme_params_;

		protected:
				CommunicatorType& communicator_;
				typename Traits::TimeProviderType timeprovider_;
				typename Traits::DiscreteStokesFunctionSpaceWrapperType functionSpaceWrapper_;
				mutable typename Traits::DiscreteStokesFunctionWrapperType currentFunctions_;
				mutable typename Traits::DiscreteStokesFunctionWrapperType nextFunctions_;
				typename Traits::DiscreteStokesFunctionWrapperType errorFunctions_;
				ExactSolutionType exactSolution_;
				mutable typename Traits::DiscreteStokesFunctionWrapperType dummyFunctions_;
				mutable typename Traits::DiscreteStokesFunctionWrapperType updateFunctions_;
				mutable typename Traits::DiscreteStokesFunctionWrapperType rhsFunctions_;
				OutputTupleType1& data_tuple_1;
				DataWriterType1 dataWriter1_;
				CheckPointerType check_pointer_;
				DataWriterType2 dataWriter2_;
				const typename Traits::OseenPassType::Traits::DiscreteSigmaFunctionSpaceType sigma_space_;
				mutable typename Traits::OseenPassType::RhsDatacontainer rhsDatacontainer_;
				mutable typename Traits::DiscreteStokesFunctionWrapperType lastFunctions_;

				typedef Stuff::L2Error< typename Traits::GridPartType >
					L2ErrorType;
				L2ErrorType l2Error_;


			public:
				const double viscosity_;
				const double d_t_;
				const double reynolds_;
				double current_max_gridwidth_;

			public:
				ThetaScheme( typename Traits::GridPartType gridPart,
							 const typename Traits::ThetaSchemeDescriptionType& scheme_params,
							 CommunicatorType comm			= Dune::MPIManager::helper().getCommunicator()
						)
					: gridPart_( gridPart ),
					scheme_params_( scheme_params ),
					communicator_( comm ),
					timeprovider_( scheme_params_, communicator_ ),
					functionSpaceWrapper_( gridPart_ ),
					currentFunctions_(  "current_",
										functionSpaceWrapper_,
										gridPart_ ),
					nextFunctions_(  "next_",
									functionSpaceWrapper_,
									gridPart_ ),
					errorFunctions_(  "error_",
									functionSpaceWrapper_,
									gridPart_ ),
					exactSolution_( timeprovider_,
									gridPart_,
									functionSpaceWrapper_ ),
					dummyFunctions_("dummy",
									functionSpaceWrapper_,
									gridPart_ ),
					updateFunctions_("updates",
									  functionSpaceWrapper_,
									  gridPart_ ),
					rhsFunctions_("rhs-adapter",
								 functionSpaceWrapper_,
								 gridPart_ ),
					data_tuple_1( TupleSerializerType1::getTuple(
							  currentFunctions_,
							  errorFunctions_,
							  exactSolution_,
							  dummyFunctions_) ),
					dataWriter1_( timeprovider_,
								 gridPart_.grid(),
								 data_tuple_1
								),
					check_pointer_(	gridPart_.grid(),
									"myGridName",
									data_tuple_1,
									timeprovider_.startTime(),
									timeprovider_.endTime(),
									static_cast<const LoadBalancerInterface*>( 0 )
								  ),
					dataWriter2_( timeprovider_,
								 gridPart_.grid(),
								 TupleSerializerType2::getTuple(
										 updateFunctions_,
										 rhsFunctions_)
								),
					sigma_space_( gridPart_ ),
					rhsDatacontainer_( currentFunctions_.discreteVelocity().space(), sigma_space_ ),
					  lastFunctions_("last",
										functionSpaceWrapper_,
										gridPart_ ),
					l2Error_( gridPart ),
					viscosity_( Parameters().getParam( "viscosity", 1.0, Dune::ValidateNotLess<double>(0.0) ) ),
					d_t_( timeprovider_.deltaT() ),
					reynolds_( 1.0 / viscosity_ ),
					current_max_gridwidth_( Dune::GridWidth::calcGridWidth( gridPart_ ) )
				{
					Logger().Info() << scheme_params_;
				}

				void nextStep( const int step, RunInfo& info )
				{
					current_max_gridwidth_ = Dune::GridWidth::calcGridWidth( gridPart_ );
					lastFunctions_.assign( currentFunctions_ );
					currentFunctions_.assign( nextFunctions_ );
					exactSolution_.project();
					const bool last_substep = ( step == ( Traits::ThetaSchemeDescriptionType::numberOfSteps_ -1) );

					//error calc
					if ( Parameters().getParam( "calculate_errors", true ) ) {
						Profiler::ScopedTiming error_time("error_calc");

						errorFunctions_.discretePressure().assign( exactSolution_.discretePressure() );
						errorFunctions_.discretePressure() -= currentFunctions_.discretePressure();
						errorFunctions_.discreteVelocity().assign( exactSolution_.discreteVelocity() );
						errorFunctions_.discreteVelocity() -= currentFunctions_.discreteVelocity();

						double meanPressure_exact = Stuff::integralAndVolume( exactSolution_.exactPressure(), currentFunctions_.discretePressure().space() ).first;
						double meanPressure_discrete = Stuff::integralAndVolume( currentFunctions_.discretePressure(), currentFunctions_.discretePressure().space() ).first;

						Dune::L2Norm< typename Traits::GridPartType > l2_Error( gridPart_ );
						Dune::H1Norm< typename Traits::GridPartType > h1_Error( gridPart_ );

//						if ( Parameters().getParam( "error_scaling", false ) ) {
//								const double scale		= 1 / std::sqrt( viscosity_ );
//								errorFunctions_.discretePressure() *= scale;
//								errorFunctions_.discreteVelocity() *= scale;
//						}

						const double l2_error_pressure_				= l2_Error.norm( errorFunctions_.discretePressure() );
						const double l2_error_velocity_				= l2_Error.norm( errorFunctions_.discreteVelocity() );
						const double h1_error_pressure_				= h1_Error.norm( errorFunctions_.discretePressure() );
						const double h1_error_velocity_				= h1_Error.norm( errorFunctions_.discreteVelocity() );
						const double relative_l2_error_pressure_	= l2_error_pressure_ / l2_Error.norm( exactSolution_.discretePressure() );
						const double relative_l2_error_velocity_	= l2_error_velocity_ / l2_Error.norm( exactSolution_.discreteVelocity() );
						const double relative_h1_error_velocity_	= h1_error_velocity_ / h1_Error.norm( exactSolution_.discreteVelocity() );
						std::vector<double> error_vector;
						error_vector.push_back( l2_error_velocity_ );
						error_vector.push_back( l2_error_pressure_ );
						std::vector<double> h1_error_vector;
						h1_error_vector.push_back( h1_error_velocity_ );
						h1_error_vector.push_back( h1_error_pressure_ );

					#ifdef NDEBUG
						if ( last_substep ) //no need to be so verbose otherwise
					#endif
						{
							Logger().Info().Resume();
							Logger().Info() << boost::format ("L2-Error Pressure (abs|rel): %e | %e \t Velocity (abs|rel): %e | %e")
												% l2_error_pressure_ % relative_l2_error_pressure_
												% l2_error_velocity_ % relative_l2_error_velocity_
										#ifndef NDEBUG
											<< boost::format ("\nH1-Error Velocity (abs|rel): %e | %e \t Mean pressure (exact|discrete) %e | %e")
												% h1_error_velocity_ % relative_h1_error_velocity_
												% meanPressure_exact % meanPressure_discrete
										#endif
											<< std::endl;
						}
						const double max_l2_error = Parameters().getParam( "max_error", 1e2, Dune::ValidateGreater<double>(0.0) );
						info.L2Errors		= error_vector;
						info.H1Errors		= h1_error_vector;
						if ( l2_error_velocity_ > max_l2_error || ( !Parameters().getParam( "parabolic", false ) && l2_error_pressure_ > max_l2_error ) )
							throw Stuff::singlerun_abort_exception( "Aborted, L2 error above " + Stuff::toString(max_l2_error) );
						if ( !Parameters().getParam( "parabolic", false ) && (std::isnan( l2_error_velocity_ ) || std::isnan( l2_error_pressure_ ) )  )
							throw Stuff::singlerun_abort_exception("L2 error is Nan");
					}
					//end error calc

					if ( last_substep ) {
						typedef Dune::StabilizationCoefficients::ValueType
							Pair;
						Dune::StabilizationCoefficients stabil_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();

						info.codim0			= gridPart_.grid().size( 0 );
						info.grid_width		= current_max_gridwidth_;
						info.run_time		= profiler().GetTiming( "full_step" );
						info.delta_t		= timeprovider_.deltaT();
						info.current_time	= timeprovider_.subTime();
						info.viscosity		= viscosity_;
						info.reynolds		= reynolds_;

						info.c11			= Pair( stabil_coeff.Power( "C11" ), stabil_coeff.Factor( "C11" ) );
						info.c12			= Pair( stabil_coeff.Power( "C12" ), stabil_coeff.Factor( "C12" ) );
						info.d11			= Pair( stabil_coeff.Power( "D11" ), stabil_coeff.Factor( "D11" ) );
						info.d12			= Pair( stabil_coeff.Power( "D12" ), stabil_coeff.Factor( "D12" ) );
						info.bfg			= Parameters().getParam( "do-bfg", true );
						info.gridname		= gridPart_.grid().name();
						info.refine_level	= Parameters().getParam( "minref", 0, Dune::ValidateNotLess<int>(0) );

						info.polorder_pressure	= Traits::OseenModelTraits::pressureSpaceOrder;
						info.polorder_sigma		= Traits::OseenModelTraits::sigmaSpaceOrder;
						info.polorder_velocity	= Traits::OseenModelTraits::velocitySpaceOrder;

						info.solver_accuracy		= Parameters().getParam( "absLimit", 1e-4 );
						info.inner_solver_accuracy	= Parameters().getParam( "inner_absLimit", 1e-4 );
						info.bfg_tau				= Parameters().getParam( "bfg-tau", 0.1 );

						info.problemIdentifier	= TESTCASE_NAME;
						info.algo_id			= scheme_params_.algo_id;
						info.extra_info			= (boost::format("%s on %s") % COMMIT % std::getenv("HOSTNAME") ).str();

						Logger().Info() << boost::format("current time (substep %d ): %f (%f)\n")
												% step
												% timeprovider_.subTime()
												% timeprovider_.previousSubTime();
					}

					if ( last_substep || !Parameters().getParam( "write_fulltimestep_only", false ) )
						writeData();
					timeprovider_.nextFractional();
				}


				void Init()
				{
					typename Traits::TimeProviderType::StepZeroGuard
						step0( timeprovider_.stepZeroGuard( d_t_ ) );
					//initial flow field at t = 0
					exactSolution_.project();
					currentFunctions_.assign( exactSolution_ );
					nextFunctions_.assign( exactSolution_ );
					writeData();
					//the guard dtor sets current time to t_0 + dt_k
				}

				RunInfoTimeMap run()
				{
					RunInfoTimeMap runInfoMap;
					Init();

					for( ;timeprovider_.time() <= timeprovider_.endTime(); )
					{
						assert( timeprovider_.time() > 0.0 );
						RunInfo info;
						if ( Parameters().getParam("old_timestep", false) )
							info = operator_split_fullstep();
						else
							info = full_timestep();
						const double real_time = timeprovider_.subTime();
						try {
							nextStep( Traits::substep_count -1 , info );
						}
						catch ( Stuff::singlerun_abort_exception& e ) {
							Logger().Err() << e.what() << std::endl;
							//fill up the map with dummy data so it can still be used in output
							runInfoMap[real_time] = info;
							for( ;timeprovider_.time() <= timeprovider_.endTime(); ) {
								timeprovider_.nextFractional();
								runInfoMap[timeprovider_.subTime()] = RunInfo::dummy();
							}
							return runInfoMap;
						}
						timeprovider_.printRemainderEstimate( Logger().Info() );
						runInfoMap[real_time] = info;
					}
					assert( runInfoMap.size() > 0 );
					return runInfoMap;
				}

				RunInfo full_timestep()
				{
					Profiler::ScopedTiming fullstep_time("full_step");
					RunInfo info;
					for ( int i=0; i < Traits::substep_count; ++i )
					{
						const double dt_k = scheme_params_.step_sizes_[i];
						if ( Parameters().getParam( "alternate_fixpoint", false ) )
							alternative_substep( dt_k, scheme_params_.thetas_[i] );
						else
							substep( dt_k, scheme_params_.thetas_[i] );
						if ( i != Traits::substep_count - 1 )
							//the last step increase is done after one call level up
							nextStep( i, info );
					}
					return info;
				}

				void substep( const double dt_k, const typename Traits::ThetaSchemeDescriptionType::ThetaValueArray& theta_values )
				{
					//build rhs
					const bool first_step = timeprovider_.timeStep() <= 1;
					const typename Traits::AnalyticalForceType force ( viscosity_,
																 currentFunctions_.discreteVelocity().space() );
					if ( Parameters().getParam( "rhs_cheat", false ) )
						cheatRHS();


					if ( scheme_params_.algo_id == Traits::ThetaSchemeDescriptionType::scheme_names[3] /*CN*/)
					{
						DiscreteVelocityFunctionType beta = currentFunctions_.discreteVelocity();
						beta *= 3.0;
						beta -= lastFunctions_.discreteVelocity();
						beta *= 0.5;
						Dune::BruteForceReconstruction< typename Traits::OseenPassType::RhsDatacontainer, typename Traits::OseenModelType >
															::getConvection( beta, rhsDatacontainer_.velocity_gradient, rhsDatacontainer_.convection );
					}
					boost::scoped_ptr< typename Traits::OseenForceAdapterFunctionType >
							ptr_oseenForce( first_step //in our very first step no previous computed data is avail. in rhs_container
												? new typename Traits::OseenForceAdapterFunctionType (	timeprovider_,
																										currentFunctions_.discreteVelocity(),
																										force,
																										reynolds_,
																										theta_values )
												: new typename Traits::OseenForceAdapterFunctionType (	timeprovider_,
																										currentFunctions_.discreteVelocity(),
																										force,
																										reynolds_,
																										theta_values,
																										rhsDatacontainer_ )
											);
					typename Traits::DiscreteStokesFunctionWrapperType
							exactSolution_at_next_time ( "reoh", exactSolution_.space(), gridPart_ );
					exactSolution_.atTime( timeprovider_.subTime(), exactSolution_at_next_time  );
					dummyFunctions_.assign( exactSolution_at_next_time );
					rhsFunctions_.discreteVelocity().assign( *ptr_oseenForce );
					typename Traits::StokesStartPassType stokesStartPass;
					typename Traits::AnalyticalDirichletDataType oseenDirichletData =
							Traits::OseenModelTraits::AnalyticalDirichletDataTraitsImplementation
											::getInstance( timeprovider_,
														   functionSpaceWrapper_ );

					unsigned int oseen_iterations = Parameters().getParam( "oseen_iterations", (unsigned int)(1), ValidateGreater<unsigned int>( 0 ) );
					const double dt_n = timeprovider_.deltaT();
					const typename L2ErrorType::Errors old_error_velocity
							= l2Error_.get( currentFunctions().discreteVelocity(), exactSolution_at_next_time.discreteVelocity() );
					const typename L2ErrorType::Errors old_error_pressure
							= l2Error_.get( currentFunctions().discretePressure(), exactSolution_at_next_time.discretePressure() );
					double velocity_error_reduction = 1.0;
					double pressure_error_reduction = 1.0;
					unsigned int i = 0;
					do
					{
						bool abort_loop = false;
						DiscreteVelocityFunctionType beta = currentFunctions_.discreteVelocity();
						if ( scheme_params_.algo_id == Traits::ThetaSchemeDescriptionType::scheme_names[3] /*CN*/)
						{
							beta *= 3.0;
							beta -= lastFunctions_.discreteVelocity();
							beta *= 0.5;
							abort_loop = true; // linCN only needs a single "iteration"
						}
						typename Traits::OseenModelType
								oseenModel( Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients(),
											*ptr_oseenForce,
											oseenDirichletData,
											theta_values[0] * dt_n / reynolds_, /*viscosity*/
											1.0f, /*alpha*/
											dt_k,/*pressure_gradient_scale_factor*/
											theta_values[0] * dt_n /*convection_scale_factor*/
						                   );
						typename Traits::OseenPassType oseenPass( stokesStartPass,
												oseenModel,
												gridPart_,
												functionSpaceWrapper_,
												beta /*beta*/,
												!Parameters().getParam( "parabolic", false ) /*do_oseen_disc*/ );
						if ( timeprovider_.timeStep() <= 1 && i < 1)
							oseenPass.printInfo();
						if ( Parameters().getParam( "silent_stokes", true ) )
							Logger().Info().Suspend( Logging::LogStream::default_suspend_priority + 10 );
						oseenPass.apply( currentFunctions_, nextFunctions_, &rhsDatacontainer_ );
						Logger().Info().Resume( Logging::LogStream::default_suspend_priority + 10 );

						{
							Profiler::ScopedTiming error_time("error_calc");
							const typename L2ErrorType::Errors new_error_velocity
									= l2Error_.get( nextFunctions_.discreteVelocity(), exactSolution_at_next_time.discreteVelocity() );
							const typename L2ErrorType::Errors new_error_pressure
									= l2Error_.get( nextFunctions_.discretePressure(), exactSolution_at_next_time.discretePressure() );
							velocity_error_reduction = old_error_velocity.absolute() / new_error_velocity.absolute();
							pressure_error_reduction = old_error_pressure.absolute() / new_error_pressure.absolute() ;
							Logger().Dbg() << boost::format(" abs diff velo %1.20e \tpress %1.20e\nabs new velo %1.20e \tpress %1.20e")
											  % ( old_error_velocity.absolute() - new_error_velocity.absolute() )
											  % ( old_error_pressure.absolute() - new_error_pressure.absolute() )
											  % new_error_velocity.absolute()
											  % new_error_pressure.absolute()
										   << std::endl;
						}

						setUpdateFunctions();
						currentFunctions_.assign( nextFunctions_ );

//						if ( ( ( pressure_error_reduction < 1.0 )
//							  && ( velocity_error_reduction < 1.0 ) ) )
//						{
//							Logger().Info() << "Oseen iteration increased error, aborting.. -- ";
//							abort_loop = true;
//						}

//						else if ( ( pressure_error_reduction > 10.0 )
//								|| ( velocity_error_reduction > 10.0 ) )
//						{
//							Logger().Info() << "Oseen iteration reduced error by factor 10, aborting.. -- ";
//							abort_loop = true;
//						}
//						else if (  ( ! ( ( last_pressure_error_reduction != pressure_error_reduction )
//									|| ( last_velocity_error_reduction != velocity_error_reduction ) ) )
//								|| ( pressure_error_reduction < Parameters().getParam( "min_error_reduction", 1.05 ) )
//								|| ( velocity_error_reduction < Parameters().getParam( "min_error_reduction", 1.05 ) ) )
//						{
//							Logger().Info() << "Oseen iteration reduced no error, aborting.. -- ";
//							abort_loop = true;
//						}
						if ( abort_loop || i++ >= oseen_iterations )
						{
							break;
						}
						else
						{
							Logger().Dbg() << boost::format(" iteration %d, error reduction: pressure  %e | velocity %e")
																		   % i % pressure_error_reduction % velocity_error_reduction
																		<< std::endl;
						}
					} while ( true ) ;
					Logger().Info() << boost::format(" iteration %d, error reduction: pressure  %e | velocity %e")
																   % i % pressure_error_reduction % velocity_error_reduction
																<< std::endl;
				}

				void alternative_substep( const double dt_k, const typename Traits::ThetaSchemeDescriptionType::ThetaValueArray& theta_values )
				{
					//build rhs
					const bool first_step = timeprovider_.timeStep() <= 1;
					const typename Traits::AnalyticalForceType force ( viscosity_,
																 currentFunctions_.discreteVelocity().space() );
					if ( Parameters().getParam( "rhs_cheat", false ) )
						cheatRHS();


//					rhsFunctions_.discreteVelocity().assign( *ptr_oseenForce );
					typename Traits::StokesStartPassType stokesStartPass;
					typename Traits::AnalyticalDirichletDataType oseenDirichletData =
							Traits::OseenModelTraits::AnalyticalDirichletDataTraitsImplementation
											::getInstance( timeprovider_,
														   functionSpaceWrapper_ );

					unsigned int oseen_iterations = Parameters().getParam( "oseen_iterations", (unsigned int)(1), ValidateGreater<unsigned int>( 0 ) );
					const double dt_n = timeprovider_.deltaT();
					const typename L2ErrorType::Errors old_error_velocity
							= l2Error_.get( currentFunctions().discreteVelocity(), exactSolution_.discreteVelocity() );
					const	typename L2ErrorType::Errors old_error_pressure
							= l2Error_.get( currentFunctions().discretePressure(), exactSolution_.discretePressure() );
					double velocity_error_reduction = 1.0;
					double pressure_error_reduction = 1.0;
					unsigned int oseen_iteration_number = 0;
					while ( true )
					{
						boost::scoped_ptr< typename Traits::OseenForceAdapterFunctionType >
								ptr_oseenForce( first_step //in our very first step no previous computed data is avail. in rhs_container
													? new typename Traits::OseenForceAdapterFunctionType (	timeprovider_,
																											currentFunctions_.discreteVelocity(),
																											force,
																											reynolds_,
																											theta_values )
													: new typename Traits::OseenForceAdapterFunctionType (	timeprovider_,
																											currentFunctions_.discreteVelocity(),
																											force,
																											reynolds_,
																											theta_values,
																											rhsDatacontainer_ )
												);
						typename Traits::OseenAltRhsForceAdapterFunctionType
							real_rhs( rhsDatacontainer_.convection );
						real_rhs *= -theta_values[0]*dt_n;
						real_rhs += *ptr_oseenForce;


						const double last_velocity_error_reduction = velocity_error_reduction;
						const double last_pressure_error_reduction = pressure_error_reduction;
						typename Traits::OseenModelAltRhsType
								oseenModel( Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients(),
											real_rhs,
											oseenDirichletData,
											theta_values[0] * dt_n / reynolds_, /*viscosity*/
											1.0f, /*alpha*/
											dt_k,/*pressure_gradient_scale_factor*/
											theta_values[0] * dt_n /*convection_scale_factor*/
										   );
						typename Traits::OseenPassAltRhsType oseenPass( stokesStartPass,
												oseenModel,
												gridPart_,
												functionSpaceWrapper_,
												currentFunctions_.discreteVelocity() /*beta*/,
												false /*do_oseen_disc*/ );
						if ( timeprovider_.timeStep() <= 1 )
							oseenPass.printInfo();
						if ( Parameters().getParam( "silent_stokes", true ) )
						{
							Logger().Info().Suspend( Logging::LogStream::default_suspend_priority + 10 );
							Logger().Dbg().Suspend( Logging::LogStream::default_suspend_priority + 10 );
						}
						oseenPass.apply( currentFunctions_, nextFunctions_, &rhsDatacontainer_ );
						Logger().Info().Resume( Logging::LogStream::default_suspend_priority + 10 );
						Logger().Dbg().Resume( Logging::LogStream::default_suspend_priority + 10 );
						Logger().Flush();
						bool velocity_error_reduced;
						bool pressure_error_reduced;
						{
							Profiler::ScopedTiming error_time("error_calc");
							typename L2ErrorType::Errors new_error_velocity
									= l2Error_.get( nextFunctions_.discreteVelocity(), exactSolution_.discreteVelocity() );
							typename L2ErrorType::Errors new_error_pressure
									= l2Error_.get( nextFunctions_.discretePressure(), exactSolution_.discretePressure() );
							velocity_error_reduction = old_error_velocity.absolute() / new_error_velocity.absolute();
							pressure_error_reduction = old_error_pressure.absolute() / new_error_pressure.absolute() ;
							const double v_diff = new_error_velocity.absolute() - old_error_velocity.absolute();
							const double p_diff = new_error_pressure.absolute() - old_error_pressure.absolute();
							Logger().Dbg()	<< boost::format(" iteration %d, new error: pressure  %e | velocity %e\ndiff %e \t%e")
													% oseen_iteration_number % new_error_pressure.absolute() % new_error_velocity.absolute()
											   % p_diff % v_diff
											<< std::endl;

							velocity_error_reduced = ( v_diff ) < 0;
							pressure_error_reduced = ( p_diff ) < 0;
						}

						currentFunctions_.assign( nextFunctions_ );

						bool abort_loop = false;
//						if ( ! ( velocity_error_reduced || pressure_error_reduced ) )
//						{
//							Logger().Info() << "Oseen iteration increased error, aborting.. -- ";
//							abort_loop = true;
//						}

						/*else*/ if ( ( pressure_error_reduction > 10.0 )
								|| ( velocity_error_reduction > 10.0 ) )
						{
							Logger().Info() << "Oseen iteration reduced error by factor 10, aborting.. -- ";
							abort_loop = true;
						}
//						else if (  ( ! ( ( last_pressure_error_reduction != pressure_error_reduction )
//									|| ( last_velocity_error_reduction != velocity_error_reduction ) ) )
//								|| ( pressure_error_reduction < Parameters().getParam( "min_error_reduction", 1.05 ) )
//								|| ( velocity_error_reduction < Parameters().getParam( "min_error_reduction", 1.05 ) ) )
//						{
//							Logger().Info() << "Oseen iteration reduced no error, aborting.. -- ";
//							abort_loop = true;
//						}
						if ( oseen_iteration_number++ >= oseen_iterations || abort_loop  )
						{
							Logger().Info() << boost::format(" iteration %d, error reduction: pressure  %e | velocity %e")
																		   % (oseen_iteration_number-1) % pressure_error_reduction % velocity_error_reduction
											<< std::endl;
							break;
						}
						Logger().Dbg()	<< boost::format(" iteration %d, error reduction: pressure  %e | velocity %e")
												% (oseen_iteration_number-1) % pressure_error_reduction % velocity_error_reduction

										<< std::endl;

					}
				}

				void setUpdateFunctions() const
				{
					updateFunctions_.assign( nextFunctions_);
					updateFunctions_ -= currentFunctions_ ;
				}

				void writeData()
				{
					Profiler::ScopedTiming io_time("IO");
					dataWriter1_.write();
					dataWriter2_.write();
//					check_pointer_.write( timeprovider_.time(), timeprovider_.timeStep() );
				}

				typename Traits::OseenPassType::RhsDatacontainer& rhsDatacontainer()
				{
					return rhsDatacontainer_;
				}

				const ExactSolutionType& exactSolution() const
				{
					return exactSolution_;
				}

				const typename Traits::DiscreteStokesFunctionWrapperType& currentFunctions() const
				{
					return currentFunctions_;
				}

				const typename Traits::TimeProviderType& timeprovider() const
				{
					return timeprovider_;
				}

				void cheatRHS()
				{
					typedef typename DiscreteVelocityFunctionType::FunctionSpaceType::FunctionSpaceType
						VelocityFunctionSpaceType;
					VelocityFunctionSpaceType continousVelocitySpace_;

					// ----
					typedef TESTING_NS::VelocityLaplace<	VelocityFunctionSpaceType,
																				typename Traits::TimeProviderType >
							VelocityLaplace;
					VelocityLaplace velocity_laplace( timeprovider_, continousVelocitySpace_ );
					Dune::BetterL2Projection
						::project( timeprovider_.previousSubTime(), velocity_laplace, rhsDatacontainer_.velocity_laplace );
					currentFunctions_.discreteVelocity().assign( exactSolution_.discreteVelocity() );
				}

				RunInfo operator_split_fullstep()
				{
					RunInfo info;
					{
						Profiler::ScopedTiming fullstep_time("full_step");
						RunInfo info_dummy;
						//stokes step A
						DiscreteVelocityFunctionType u_n( "u_n", dummyFunctions_.discreteVelocity().space() );
						u_n.assign( currentFunctions_.discreteVelocity() );
						stokesStep( scheme_params_.step_sizes_[0], scheme_params_.thetas_[0] );
						nextStep( 0, info_dummy );

						Parameters().setParam( "reduced_oseen_solver", true );
						//Nonlinear step
						nonlinearStep( scheme_params_.step_sizes_[1], scheme_params_.thetas_[1], u_n );
						nextStep( 1, info_dummy );
						Parameters().setParam( "reduced_oseen_solver", false );

						//stokes step B
						info = stokesStep( scheme_params_.step_sizes_[2], scheme_params_.thetas_[2] );
					}
					nextStep( 2, info );

					return info;
				}

				struct DiscretizationWeights {
					const double theta,alpha,beta,theta_times_delta_t,viscosity,one_neg_two_theta_dt;
					DiscretizationWeights( const double d_t, const double visc ):
						theta ( 1.0 - (std::sqrt(2)/2.0f) ),
						alpha ( ( 1.0-2*theta ) / ( 1.0-theta ) ),
						beta ( 1.0 - alpha ),
						theta_times_delta_t(theta * d_t),
						viscosity( visc ),
						one_neg_two_theta_dt( ( 1. - 2. * theta ) * d_t )
					{}
				};

				RunInfo stokesStep( const double /*dt_k*/,const typename Traits::ThetaSchemeDescriptionType::ThetaValueArray& /*theta_values*/ ) const
				{
					DiscretizationWeights discretization_weights(d_t_, viscosity_);

					if ( Parameters().getParam( "silent_stokes", true ) )
						Logger().Suspend( Logging::LogStream::default_suspend_priority + 1 );

					const bool first_stokes_step = timeprovider_.timeStep() <= 1;
					const typename Traits::AnalyticalForceType force ( viscosity_,
																 currentFunctions_.discreteVelocity().space() );

					boost::scoped_ptr< typename Traits::StokesForceAdapterType >
							ptr_stokesForce_vanilla ( first_stokes_step
												? new typename Traits::StokesForceAdapterType ( timeprovider_,
																										  currentFunctions_.discreteVelocity(),
																										  force,
																										  discretization_weights )
												: new typename Traits::StokesForceAdapterType ( timeprovider_,
																										  currentFunctions_.discreteVelocity(),
																										  force,
																										  discretization_weights,
																										  rhsDatacontainer_ )
											);

					typedef Stuff::L2Error<typename Traits::GridPartType>
							L2ErrorType;
					L2ErrorType l2Error( gridPart_ );

					// CHEAT (projecting the anaylitcal evals into the container filled by last pass
					const bool do_cheat = Parameters().getParam( "rhs_cheat", false ) && !first_stokes_step ;
					dummyFunctions_.discreteVelocity().assign( currentFunctions_.discreteVelocity() );
//					if ( do_cheat ) //do cheat rhs assembly unconditionally, below we'll choose according to do_cheat which rhs to put into the model
					{
						typedef typename DiscreteVelocityFunctionType::FunctionSpaceType::FunctionSpaceType
							VelocityFunctionSpaceType;
						VelocityFunctionSpaceType continousVelocitySpace_;
						typedef TESTING_NS::VelocityConvection<	VelocityFunctionSpaceType,
																typename Traits::TimeProviderType >
							VelocityConvection;
						VelocityConvection velocity_convection( timeprovider_, continousVelocitySpace_ );
						Dune::BetterL2Projection //we need evals from the _previous_ (t_0) step
							::project( timeprovider_.previousSubTime(), velocity_convection, rhsDatacontainer_.convection );
//						// ----
						typedef TESTING_NS::VelocityLaplace<	VelocityFunctionSpaceType,
																					typename Traits::TimeProviderType >
								VelocityLaplace;
						VelocityLaplace velocity_laplace( timeprovider_, continousVelocitySpace_ );
						Dune::BetterL2Projection //this seems currently inconsequential to the produced error
							::project( timeprovider_.previousSubTime(), velocity_laplace, rhsDatacontainer_.velocity_laplace );

//						typename L2ErrorType::Errors errors_convection = l2Error.get(	exactSolution_.discreteVelocity() ,
//																			currentFunctions_.discreteVelocity(),
//																			dummyFunctions_.discreteVelocity() );
//						std::cerr << "BLAH " << errors_convection.str();

						currentFunctions_.discreteVelocity().assign( exactSolution_.discreteVelocity() );
					}// END CHEAT

					boost::scoped_ptr< typename Traits::StokesForceAdapterType >
							ptr_stokesForce ( first_stokes_step
												? new typename Traits::StokesForceAdapterType ( timeprovider_,
																										  currentFunctions_.discreteVelocity(),
																										  force,
																										  discretization_weights )
												: new typename Traits::StokesForceAdapterType ( timeprovider_,
																										  currentFunctions_.discreteVelocity(),
																										  force,
																										  discretization_weights,
																										  rhsDatacontainer_ )
											);
					typename L2ErrorType::Errors errors_rhs = l2Error.get(	static_cast<typename Traits::StokesForceAdapterType::BaseType>(*ptr_stokesForce),
																		static_cast<typename Traits::StokesForceAdapterType::BaseType>(*ptr_stokesForce_vanilla),
																		dummyFunctions_.discreteVelocity() );
					std::cerr << "RHS " << errors_rhs.str();

					Dune::StabilizationCoefficients stab_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();

					{
//					if ( Parameters().getParam( "stab_coeff_visc_scale", true ) ) {
//						stab_coeff.Factor( "D11", ( 1 / stokes_viscosity ) );
//						stab_coeff.Factor( "C11", stokes_viscosity );
//					}
//					else {
//						stab_coeff.FactorFromParams("D11");
//						stab_coeff.FactorFromParams("C11");
//					}
//					stab_coeff.FactorFromParams("D12");
//					stab_coeff.FactorFromParams("C12");
//					stab_coeff.Add( "E12", 0.5 );
//					stab_coeff.print( Logger().Info() );
					}

					typename Traits::AnalyticalDirichletDataType stokesDirichletData =
							Traits::StokesModelTraits::AnalyticalDirichletDataTraitsImplementation
											::getInstance( timeprovider_,
														   functionSpaceWrapper_ );

					typename Traits::StokesModelType
							stokesModel(stab_coeff,
										do_cheat ? *ptr_stokesForce : *ptr_stokesForce_vanilla,
										stokesDirichletData,
										discretization_weights.alpha * discretization_weights.theta_times_delta_t, /*viscosity*/
										1.0, /*alpha*/
										0.0,/*convection_scale_factor*/
										discretization_weights.theta_times_delta_t /*pressure_gradient_scale_factor*/);
					typename Traits::StokesStartPassType stokesStartPass;
					typename Traits::StokesPassType stokesPass( stokesStartPass,
											stokesModel,
											gridPart_,
											functionSpaceWrapper_,
											dummyFunctions_.discreteVelocity(),
											false );

					stokesPass.apply( currentFunctions_, nextFunctions_, &rhsDatacontainer_ );
					setUpdateFunctions();
					RunInfo info;
					stokesPass.getRuninfo( info );
					if ( Parameters().getParam( "silent_stokes", true ) )
						Logger().Resume( Logging::LogStream::default_suspend_priority + 1 );
					return info;
				}

				void nonlinearStep( const double /*dt_k*/,const typename Traits::ThetaSchemeDescriptionType::ThetaValueArray& /*theta_values*/, const DiscreteVelocityFunctionType& u_n )
				{
					DiscretizationWeights discretization_weights(d_t_, viscosity_);


					const typename Traits::AnalyticalForceType force ( viscosity_,
																	  currentFunctions_.discreteVelocity().space() );

					// CHEAT (projecting the anaylitcal evals into the container filled by last pass
					if ( Parameters().getParam( "rhs_cheat", false ) ) {
						typedef typename DiscreteVelocityFunctionType::FunctionSpaceType::FunctionSpaceType
								VelocityFunctionSpaceType;
						VelocityFunctionSpaceType continousVelocitySpace_;

						typedef TESTING_NS::PressureGradient<	VelocityFunctionSpaceType,
								typename Traits::TimeProviderType >
								PressureGradient;
						PressureGradient pressure_gradient( timeprovider_, continousVelocitySpace_ );
						Dune::BetterL2Projection //we need evals from the _previous_ (t_0) step
								::project( timeprovider_.previousSubTime(), pressure_gradient, rhsDatacontainer_.pressure_gradient );
						// ----
						typedef TESTING_NS::VelocityLaplace<	VelocityFunctionSpaceType,
								typename Traits::TimeProviderType >
								VelocityLaplace;
						VelocityLaplace velocity_laplace( timeprovider_, continousVelocitySpace_ );
						Dune::BetterL2Projection
								::project( timeprovider_.previousSubTime(), velocity_laplace, rhsDatacontainer_.velocity_laplace );
						currentFunctions_.discreteVelocity().assign( exactSolution_.discreteVelocity() );
					}// END CHEAT

					typename Traits::NonlinearForceAdapterType nonlinearForce( timeprovider_,
																					  currentFunctions_.discreteVelocity(),
																					  force,
																					  discretization_weights,
																					  rhsDatacontainer_ );

					rhsFunctions_.discreteVelocity().assign( nonlinearForce );
					unsigned int oseen_iterations = Parameters().getParam( "oseen_iterations", (unsigned int)(1) );
					assert( oseen_iterations > 0 );
					nonlinearStepSingle( nonlinearForce, discretization_weights, u_n );
				}

				template < class T >
				void nonlinearStepSingle(	const T& nonlinearForce,
											const DiscretizationWeights& discretization_weights,
											const DiscreteVelocityFunctionType& u_n)
				{
					typename Traits::StokesStartPassType stokesStartPass;

					typename Traits::AnalyticalDirichletDataType stokesDirichletData =
							Traits::StokesModelTraits::AnalyticalDirichletDataTraitsImplementation
							::getInstance( timeprovider_,
										  functionSpaceWrapper_ );
					Dune::StabilizationCoefficients stab_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();
//					if ( Parameters().getParam( "stab_coeff_visc_scale", true ) ) {
//						stab_coeff.Factor( "D11", ( 1 / oseen_viscosity )  );
//						stab_coeff.Factor( "C11", oseen_viscosity );
//					}
//					else {
//						stab_coeff.FactorFromParams("D11");
//						stab_coeff.FactorFromParams("C11");
//					}
//					stab_coeff.FactorFromParams("D12");
//					stab_coeff.FactorFromParams("C12");
//					stab_coeff.Add( "E12", 0.5 );

//					stab_coeff.print( Logger().Info() );

					DiscreteVelocityFunctionType beta( "beta", dummyFunctions_.discreteVelocity().space() );
					DiscreteVelocityFunctionType tmp( "tmp", dummyFunctions_.discreteVelocity().space() );
					tmp.assign( u_n );
					const double theta = discretization_weights.theta;
					tmp *= (2.0*theta) / (1.0 - theta);
					beta.assign( currentFunctions_.discreteVelocity() );
					beta *= theta / (1.0-theta);
					beta += tmp;


					typename Traits::NonlinearModelType
							stokesModel(stab_coeff,
										nonlinearForce,
										stokesDirichletData,
										discretization_weights.beta * discretization_weights.one_neg_two_theta_dt, /*viscosity*/
										1.0, /*alpha*/
										discretization_weights.one_neg_two_theta_dt, /*convection_scale_factor*/
										0.0 /*pressure_gradient_scale_factor*/ );
					typename Traits::NonlinearPassType oseenPass( stokesStartPass,
															 stokesModel,
															 gridPart_,
															 functionSpaceWrapper_,
															 beta,
															 true );
					oseenPass.apply( currentFunctions_, nextFunctions_, &rhsDatacontainer_ );


				}
		};
	}//end namespace NavierStokes
}//end namespace Dune

#endif // METADATA_HH
