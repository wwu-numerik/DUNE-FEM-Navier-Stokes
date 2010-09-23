#ifndef METADATA_HH
#define METADATA_HH

#include <dune/stokes/discretestokesmodelinterface.hh>
#include <dune/stokes/stokespass.hh>
#include <dune/navier/fractionaltimeprovider.hh>
#include <dune/navier/stokestraits.hh>
#include <dune/navier/testdata.hh>
#include <dune/navier/exactsolution.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/stuff/datawriter.hh>
#include <dune/stuff/functions.hh>
#include <dune/stuff/tuple.hh>
#include <dune/stuff/customprojection.hh>
#include <dune/common/collectivecommunication.hh>
#include <cmath>

namespace Dune {
	namespace NavierStokes {
		template <	class CommunicatorImp,
					class GridPartImp,
					template < class > class AnalyticalForceImp,
					template < class > class AnalyticalDirichletDataImp,
					template < class,class > class ExactPressureImp,
					template < class,class > class ExactVelocityImp,
					int gridDim, int sigmaOrder, int velocityOrder = sigmaOrder, int pressureOrder = sigmaOrder >
		struct ThetaSchemeTraits {
			typedef ThetaSchemeTraits<CommunicatorImp,
										GridPartImp,
										AnalyticalForceImp,
										AnalyticalDirichletDataImp,
										ExactPressureImp,
										ExactVelocityImp,
										gridDim, sigmaOrder, velocityOrder , pressureOrder >
				ThisType;
			typedef GridPartImp
				GridPartType;
			typedef FractionalTimeProvider<CommunicatorImp>
				TimeProviderType;

			typedef StokesStep::DiscreteStokesModelTraits<
						TimeProviderType,
						GridPartType,
						AnalyticalForceImp,
						AnalyticalDirichletDataImp,
						gridDim,
						sigmaOrder,
						velocityOrder,
						pressureOrder >
					StokesModelTraits;
			typedef typename StokesModelTraits::PressureFunctionSpaceType
				PressureFunctionSpaceType;
			typedef typename StokesModelTraits::VelocityFunctionSpaceType
				VelocityFunctionSpaceType;

			typedef Dune::DiscreteStokesModelDefault< StokesModelTraits >
				StokesModelType;
			typedef typename StokesModelTraits::DiscreteStokesFunctionSpaceWrapperType
				DiscreteStokesFunctionSpaceWrapperType;

			typedef typename StokesModelTraits::DiscreteStokesFunctionWrapperType
				DiscreteStokesFunctionWrapperType;
			typedef typename StokesModelTraits::AnalyticalForceFunctionType
				AnalyticalForceType;
			typedef typename StokesModelTraits::AnalyticalForceAdapterType
				StokesAnalyticalForceAdapterType;
			typedef typename StokesModelTraits::AnalyticalDirichletDataType
				AnalyticalDirichletDataType;

			typedef Dune::StartPass< DiscreteStokesFunctionWrapperType, -1 >
				StokesStartPassType;
			typedef Dune::StokesPass< StokesModelType, StokesStartPassType, 0 >
				StokesPassType;

			typedef CommunicatorImp
				CommunicatorType;

			typedef ExactPressureImp< typename StokesModelTraits::PressureFunctionSpaceType,
									  TimeProviderType >
				ExactPressureType;
			typedef ExactVelocityImp< typename StokesModelTraits::VelocityFunctionSpaceType,
									  TimeProviderType >
				ExactVelocityType;

			typedef ExactSolution<ThisType>
				ExactSolutionType;

			typedef NonlinearStep::ForceAdapterFunction<	TimeProviderType,
															AnalyticalForceType,
															typename DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType,
															typename DiscreteStokesFunctionWrapperType::DiscretePressureFunctionType,
															typename StokesModelTraits::DiscreteSigmaFunctionType>
				NonlinearForceAdapterFunctionType;
			typedef NonlinearStep::DiscreteStokesModelTraits<
						TimeProviderType,
						GridPartType,
						NonlinearForceAdapterFunctionType,
						typename ExactSolutionType::DiscreteVelocityFunctionType,
						AnalyticalDirichletDataImp,
						gridDim,
						sigmaOrder,
						velocityOrder,
						pressureOrder >
				OseenModelTraits;
			typedef Dune::DiscreteStokesModelDefault< OseenModelTraits >
				OseenModelType;
			typedef Dune::StokesPass< OseenModelType,StokesStartPassType, 0 >
				OseenpassType;
		};

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

				typename Traits::GridPartType gridPart_;
				const double theta_;
				const double operator_weight_alpha_;
				const double operator_weight_beta_;
				CommunicatorType& communicator_;
				typename Traits::TimeProviderType timeprovider_;
				typename Traits::DiscreteStokesFunctionSpaceWrapperType functionSpaceWrapper_;
				typename Traits::DiscreteStokesFunctionWrapperType currentFunctions_;
				typename Traits::DiscreteStokesFunctionWrapperType nextFunctions_;
				typename Traits::DiscreteStokesFunctionWrapperType errorFunctions_;
				typename Traits::DiscreteStokesFunctionWrapperType dummyFunctions_;
				typename Traits::DiscreteStokesFunctionWrapperType updateFunctions_;
				ExactSolutionType exactSolution_;
				DataWriterType1 dataWriter1_;
				DataWriterType2 dataWriter2_;
				const typename Traits::StokesPassType::DiscreteSigmaFunctionSpaceType sigma_space_;
				typename Traits::StokesPassType::RhsDatacontainer rhsDatacontainer_;

				//constants
				const double viscosity_;
				const double d_t_;
				const double reynolds_;
				const double beta_qout_re_;
				double current_max_gridwidth_;

			public:
				ThetaScheme( typename Traits::GridPartType gridPart,
							 const double theta = 1 - std::pow( 2.0, -1/2.0 ),
							 CommunicatorType comm = Dune::MPIManager::helper().getCommunicator()
						)
					: gridPart_( gridPart ),
					theta_(theta),
					operator_weight_alpha_( ( 1-2*theta_ ) / ( 1-theta_ ) ),
					operator_weight_beta_( 1 - operator_weight_alpha_ ),
					communicator_( comm ),
					timeprovider_( theta_,operator_weight_alpha_,operator_weight_beta_, communicator_ ),
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
					dummyFunctions_("force",
									functionSpaceWrapper_,
									gridPart_ ),
					updateFunctions_("updates",
									  functionSpaceWrapper_,
									  gridPart_ ),
					dataWriter1_( timeprovider_,
								 gridPart_.grid(),
								 TupleSerializerType1::getTuple(
										 currentFunctions_,
										 errorFunctions_,
										 exactSolution_,
										 dummyFunctions_)
								),
					dataWriter2_( timeprovider_,
								 gridPart_.grid(),
								 TupleSerializerType2::getTuple(
										 updateFunctions_)
								),
					sigma_space_( gridPart_ ),
					rhsDatacontainer_( currentFunctions_.discreteVelocity().space(), sigma_space_ ),
					viscosity_( Parameters().getParam( "viscosity", 1.0 ) ),
					d_t_( timeprovider_.deltaT() ),
					reynolds_( 1.0 / viscosity_ ),
					beta_qout_re_( operator_weight_beta_ / reynolds_ ),
					current_max_gridwidth_( Dune::GridWidth::calcGridWidth( gridPart_ ) )
				{}

				void nextStep( const int step, RunInfo& info )
				{
					current_max_gridwidth_ = Dune::GridWidth::calcGridWidth( gridPart_ );
					currentFunctions_.assign( nextFunctions_ );

					//error calc
					if ( Parameters().getParam( "calculate_errors", true ) ) {
						exactSolution_.project();
						errorFunctions_.discretePressure().assign( exactSolution_.discretePressure() );
						errorFunctions_.discretePressure() -= currentFunctions_.discretePressure();
						errorFunctions_.discreteVelocity().assign( exactSolution_.discreteVelocity() );
						errorFunctions_.discreteVelocity() -= currentFunctions_.discreteVelocity();

						double meanPressure_exact = Stuff::integralAndVolume( exactSolution_.exactPressure(), currentFunctions_.discretePressure().space() ).first;
						double meanPressure_discrete = Stuff::integralAndVolume( currentFunctions_.discretePressure(), currentFunctions_.discretePressure().space() ).first;

						Dune::L2Norm< typename Traits::GridPartType > l2_Error( gridPart_ );

						if ( Parameters().getParam( "error_scaling", false ) ) {
								const double scale		= 1 / std::sqrt( viscosity_ );
								errorFunctions_.discretePressure() *= scale;
								errorFunctions_.discreteVelocity() *= scale;
						}

						const double l2_error_pressure_				= l2_Error.norm( errorFunctions_.discretePressure() );
						const double l2_error_velocity_				= l2_Error.norm( errorFunctions_.discreteVelocity() );
						const double relative_l2_error_pressure_	= l2_error_pressure_ / l2_Error.norm( exactSolution_.discretePressure() );
						const double relative_l2_error_velocity_	= l2_error_velocity_ / l2_Error.norm( exactSolution_.discreteVelocity() );
						std::vector<double> error_vector;
						error_vector.push_back( l2_error_velocity_ );
						error_vector.push_back( l2_error_pressure_ );


						Logger().Info().Resume();
						Logger().Info() << "L2-Error Pressure (abs|rel): " << std::setw(8) << l2_error_pressure_ << " | " << relative_l2_error_pressure_ << "\n"
										<< "L2-Error Velocity (abs|rel): " << std::setw(8) << l2_error_velocity_ << " | " << relative_l2_error_velocity_ << "\n"
										<< "Mean pressure (exact|discrete): " << meanPressure_exact << " | " << meanPressure_discrete << std::endl;
						const double max_l2_error = 1e4;
						if ( l2_error_velocity_ > max_l2_error )
							DUNE_THROW(MathError, "Aborted, L2 error above " << max_l2_error );
						info.L2Errors		= error_vector;
					}
					//end error calc

					if (step == 3) {
						typedef Dune::StabilizationCoefficients::ValueType
							Pair;
						Dune::StabilizationCoefficients stabil_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();

						info.codim0			= gridPart_.grid().size( 0 );
						info.grid_width		= current_max_gridwidth_;
						info.run_time		= profiler().GetTiming( "Timestep" );
						info.delta_t		= timeprovider_.deltaT();
						info.current_time	= timeprovider_.time();
						info.viscosity		= viscosity_;
						info.reynolds		= reynolds_;

						info.c11			= Pair( stabil_coeff.Power( "C11" ), stabil_coeff.Factor( "C11" ) );
						info.c12			= Pair( stabil_coeff.Power( "C12" ), stabil_coeff.Factor( "C12" ) );
						info.d11			= Pair( stabil_coeff.Power( "D11" ), stabil_coeff.Factor( "D11" ) );
						info.d12			= Pair( stabil_coeff.Power( "D12" ), stabil_coeff.Factor( "D12" ) );
						info.bfg			= Parameters().getParam( "do-bfg", true );
						info.gridname		= gridPart_.grid().name();
						info.refine_level	= Parameters().getParam( "minref", 0 );

						info.polorder_pressure	= Traits::StokesModelTraits::pressureSpaceOrder;
						info.polorder_sigma		= Traits::StokesModelTraits::sigmaSpaceOrder;
						info.polorder_velocity	= Traits::StokesModelTraits::velocitySpaceOrder;

						info.solver_accuracy		= Parameters().getParam( "absLimit", 1e-4 );
						info.inner_solver_accuracy	= Parameters().getParam( "inner_absLimit", 1e-4 );
						info.bfg_tau				= Parameters().getParam( "bfg-tau", 0.1 );

						info.problemIdentifier = TESTCASE_NAME;
					}

					std::cout << "current time (substep " << step << "): " << timeprovider_.subTime() << std::endl;

					if ( step == 3 || !Parameters().getParam( "write_fulltimestep_only", false ) )
						writeData();
					timeprovider_.nextFractional();
				}

				RunInfo stokesStep()
				{
					const bool scale_equations = Parameters().getParam( "scale_equations", false );
					const double delta_t_factor = theta_ * d_t_;
					double stokes_alpha,scale_factor,stokes_viscosity;
					const double stokes_alpha_unscaled = 1 / delta_t_factor;
					if ( scale_equations ) {
						stokes_alpha = 1;
						scale_factor = delta_t_factor;
						stokes_viscosity = operator_weight_alpha_ * delta_t_factor / reynolds_;
					}
					else {
						stokes_alpha = stokes_alpha_unscaled;
						scale_factor = 1;
						stokes_viscosity = operator_weight_alpha_ / reynolds_;
					}


					if ( Parameters().getParam( "silent_stokes", true ) )
						Logger().Suspend( Logging::LogStream::default_suspend_priority + 1 );
					const bool first_stokes_step = timeprovider_.timeStep() <= 1;
					const bool add_extra_terms = Parameters().getParam( "add_extra_terms", true ) || first_stokes_step;
					const double stokes_T = theta_ * d_t_;
					const typename Traits::AnalyticalForceType force ( viscosity_,
																 currentFunctions_.discreteVelocity().space() );
					typename Traits::StokesAnalyticalForceAdapterType stokesForce( timeprovider_,
																				   currentFunctions_.discreteVelocity(),
																				   force,
																				   beta_qout_re_,
																				   stokes_alpha_unscaled,
																				   add_extra_terms );
					typename Traits::AnalyticalDirichletDataType stokesDirichletData =
							Traits::StokesModelTraits::AnalyticalDirichletDataTraitsImplementation
											::getInstance( timeprovider_,
														   functionSpaceWrapper_ );
					double meanGD
							= Stuff::boundaryIntegral( stokesDirichletData, currentFunctions_.discreteVelocity().space() );
					if ( !add_extra_terms ) {
						// F = f + \alpha / \Re * laplace u + ( 1/(theta * tau) ) u - ( u * nable ) u
						rhsDatacontainer_.velocity_laplace *= ( operator_weight_beta_ / reynolds_ );
						dummyFunctions_.discreteVelocity().assign( currentFunctions_.discreteVelocity() );
						dummyFunctions_.discreteVelocity() *= stokes_alpha_unscaled;
						stokesForce += dummyFunctions_.discreteVelocity();
						stokesForce -= rhsDatacontainer_.convection;
						stokesForce += rhsDatacontainer_.velocity_laplace;

						typename Traits::StokesAnalyticalForceAdapterType stokesForce_full( timeprovider_,
																					   currentFunctions_.discreteVelocity(),
																					   force,
																					   beta_qout_re_,
																					   stokes_alpha_unscaled,
																					   true );

						// L2 error
						Dune::L2Norm< typename Traits::GridPartType > l2_Error( gridPart_ );
						const double force_abs = l2_Error.norm( stokesForce_full ) ;
						stokesForce_full -= stokesForce;
						const double force_error_abs = l2_Error.norm( stokesForce_full ) ;
						Logger().Info() << boost::format( "rhs error %f\n" ) % force_error_abs;
					}

					Dune::StabilizationCoefficients stab_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();

					if ( Parameters().getParam( "stab_coeff_visc_scale", true ) ) {
						stab_coeff.Factor( "D11", ( 1 / stokes_viscosity ) );
						stab_coeff.Factor( "C11", stokes_viscosity );
					}
					else {
						stab_coeff.FactorFromParams("D11");
						stab_coeff.FactorFromParams("C11");
					}
					stab_coeff.FactorFromParams("D12");
					stab_coeff.FactorFromParams("C12");
					stab_coeff.Add( "E12", 0.5 );

					dummyFunctions_.discreteVelocity().assign( stokesForce );
//					Logger().Info() << "stokes a/RE|b/RE|y " << stokes_viscosity_ << " | "
//														<< beta_qout_re_ << " | "
//														<< quasi_stokes_alpha_ << "\n";
					stab_coeff.print( Logger().Info() );

					stokesForce *= scale_factor;
					typename Traits::StokesModelType
							stokesModel(stab_coeff,
										stokesForce,
										stokesDirichletData,
										stokes_viscosity ,
										stokes_alpha, scale_factor, scale_factor );
					typename Traits::StokesStartPassType stokesStartPass;
					typename Traits::StokesPassType stokesPass( stokesStartPass,
											stokesModel,
											gridPart_,
											functionSpaceWrapper_ );

					stokesPass.apply( currentFunctions_, nextFunctions_, &rhsDatacontainer_ );
					setUpdateFunctions();
					RunInfo info;
					stokesPass.getRuninfo( info );
					if ( Parameters().getParam( "silent_stokes", true ) )
						Logger().Resume( Logging::LogStream::default_suspend_priority + 1 );
					Logger().Info() << "Dirichlet boundary integral " << meanGD << std::endl;
					return info;
				}

				RunInfoVector run()
				{
					RunInfoVector runInfoVector;

					typename Traits::DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType::RangeType meanVelocity
							= Stuff::meanValue( currentFunctions_.discreteVelocity(), currentFunctions_.discreteVelocity().space() );
					const double typicalVelocity = meanVelocity.two_norm();
					Stuff::printFieldVector( meanVelocity, "meanVelocity", Logger().Info() );
					Logger().Info() << "typicalVelocity " << typicalVelocity << std::endl;
					const double typicalLength = 1.0;

					const double grid_width = Dune::GridWidth::calcGridWidth( gridPart_ );
					double dt_new = grid_width * grid_width * 2;

					timeprovider_.init( d_t_ );
					//initial flow field at t = 0
					exactSolution_.project();
					currentFunctions_.assign( exactSolution_ );
					nextFunctions_.assign( exactSolution_ );
					writeData();
					//set current time to t_0 + theta
					timeprovider_.nextFractional();

					for( ;timeprovider_.time() < timeprovider_.endTime(); )
					{
						RunInfo info_dummy;
						profiler().StartTiming( "Timestep" );
						//stokes step A
						if( Parameters().getParam( "enable_stokesA", true ) )
							stokesStep();
						else {
							exactSolution_.project();
							nextFunctions_.assign( exactSolution_ );
						}
						nextStep( 1, info_dummy );
						if( Parameters().getParam( "pressure_cheat", false ) )
							nextFunctions_.discretePressure().assign( exactSolution_.discretePressure() );
						//Nonlinear step
						if( Parameters().getParam( "enable_oseen", true ) )
							oseenStep();
						else {
							exactSolution_.project();
							nextFunctions_.assign( exactSolution_ );
						}

						nextStep( 2, info_dummy );
						if( Parameters().getParam( "pressure_cheat", false ) )
							nextFunctions_.discretePressure().assign( exactSolution_.discretePressure() );
						//stokes step B
						RunInfo info;
						if( Parameters().getParam( "enable_stokesB", true ) )
							stokesStep();
						else {
							exactSolution_.project();
							nextFunctions_.assign( exactSolution_ );
						}

						profiler().StopTiming( "Timestep" );
						nextStep( 3, info );
						if( Parameters().getParam( "pressure_cheat", false ) )
							nextFunctions_.discretePressure().assign( exactSolution_.discretePressure() );
						runInfoVector.push_back( info );
					}
					return runInfoVector;
				}

				void oseenStep()
				{
					const bool scale_equations = Parameters().getParam( "scale_equations", false );
					const double delta_t_factor = ( 1 - 2 * theta_ ) * d_t_;
					double oseen_alpha, oseen_viscosity,scale_factor;
					const double oseen_alpha_unscaled = 1 / delta_t_factor;
					if ( scale_equations ) {
						oseen_alpha = 1;
						scale_factor = delta_t_factor;
						oseen_viscosity = beta_qout_re_ * scale_factor;
					}
					else {
						oseen_alpha = oseen_alpha_unscaled;
						oseen_viscosity = beta_qout_re_;
						scale_factor = 1;
					}

					const typename Traits::AnalyticalForceType force ( viscosity_,
																 currentFunctions_.discreteVelocity().space() );
					const bool add_extra_terms = Parameters().getParam( "add_extra_terms", true );

					typename Traits::NonlinearForceAdapterFunctionType nonlinearForce( timeprovider_,
																	  currentFunctions_.discreteVelocity(),
																	  currentFunctions_.discretePressure(),
																	  force,
																	  operator_weight_alpha_ / reynolds_,
																	  oseen_alpha_unscaled,
																	  add_extra_terms );
					typename Traits::NonlinearForceAdapterFunctionType nonlinearForce_full( timeprovider_,
																	  currentFunctions_.discreteVelocity(),
																	  currentFunctions_.discretePressure(),
																	  force,
																	  operator_weight_alpha_ / reynolds_,
																	  oseen_alpha_unscaled,
																	  true );


					// F = f + \alpha \Re \delta u - \nabla p + ( 1/(1-2 \theta) ) * u
					if ( !add_extra_terms ) {
						nonlinearForce -= rhsDatacontainer_.pressure_gradient;
						rhsDatacontainer_.velocity_laplace *= operator_weight_alpha_ / reynolds_;
						nonlinearForce += rhsDatacontainer_.velocity_laplace;
						dummyFunctions_.discreteVelocity().assign( currentFunctions_.discreteVelocity() );
						dummyFunctions_.discreteVelocity() *= ( 1 / delta_t_factor );
						nonlinearForce += dummyFunctions_.discreteVelocity();
					}
					nonlinearForce *= scale_factor;
					for( unsigned int i = 0; i<Parameters().getParam("oseen_iterations",int(1)); ++i )
					{
						oseenStepSingle( nonlinearForce, oseen_viscosity, oseen_alpha, scale_factor );
						setUpdateFunctions();
						currentFunctions_.assign( nextFunctions_ );
					}
				}

				void oseenStepSingle(	const typename Traits::NonlinearForceAdapterFunctionType& nonlinearForce,
										const double oseen_viscosity,
										const double oseen_alpha,
										const double scale_factor )
				{

					//
					typename Traits::StokesStartPassType stokesStartPass;

					typename Traits::AnalyticalDirichletDataType stokesDirichletData =
							Traits::StokesModelTraits::AnalyticalDirichletDataTraitsImplementation
											::getInstance( timeprovider_,
														   functionSpaceWrapper_ );
					Dune::StabilizationCoefficients stab_coeff = Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();
					if ( Parameters().getParam( "stab_coeff_visc_scale", true ) ) {
						stab_coeff.Factor( "D11", ( 1 / oseen_viscosity )  );
						stab_coeff.Factor( "C11", oseen_viscosity );
					}
					else {
						stab_coeff.FactorFromParams("D11");
						stab_coeff.FactorFromParams("C11");
					}
					stab_coeff.FactorFromParams("D12");
					stab_coeff.FactorFromParams("C12");
					stab_coeff.Add( "E12", 0.5 );

//					Logger().Info() << "oseen a/RE|b/RE|y " << operator_weight_alpha_ / reynolds_ << " | "
//														<< beta_qout_re_ << " | "
//														<< oseen_alpha << "\n";
					stab_coeff.print( Logger().Info() );

					typename Traits::OseenModelType
							stokesModel(stab_coeff,
										nonlinearForce,
										stokesDirichletData,
										oseen_viscosity,
										oseen_alpha,
										scale_factor, scale_factor);
					typename Traits::OseenpassType oseenPass( stokesStartPass,
											stokesModel,
											gridPart_,
											functionSpaceWrapper_,
											&currentFunctions_.discreteVelocity() );
					oseenPass.apply( currentFunctions_, nextFunctions_ );

				}

				void setUpdateFunctions()
				{
					updateFunctions_.assign( currentFunctions_ );
					updateFunctions_ -= nextFunctions_;
				}

				void writeData()
				{
					dataWriter1_.write();
					dataWriter2_.write();
				}

		};
	}//end namespace NavierStokes
}//end namespace Dune

#endif // METADATA_HH
