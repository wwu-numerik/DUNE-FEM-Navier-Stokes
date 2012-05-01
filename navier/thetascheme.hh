#ifndef METADATA_HH
#define METADATA_HH

#include <dune/navier/thetascheme_base.hh>

namespace Dune {
	namespace NavierStokes {
		template < class Traits >
		class ThetaScheme : public ThetaSchemeBase< Traits > {
			protected:
				typedef ThetaSchemeBase< Traits >
					BaseType;

				using BaseType::gridPart_;
				using BaseType::scheme_params_;

		protected:
				using BaseType::timeprovider_;
				using BaseType::functionSpaceWrapper_;
				using BaseType::currentFunctions_;
				using BaseType::nextFunctions_;
				using BaseType::exactSolution_;
				using BaseType::dummyFunctions_;
				using BaseType::rhsFunctions_;
				using BaseType::rhsDatacontainer_;
				using BaseType::lastFunctions_;
				using BaseType::l2Error_;


			public:
				using BaseType::viscosity_;
				using BaseType::reynolds_;

			public:
				ThetaScheme( typename Traits::GridPartType gridPart,
							 const typename Traits::ThetaSchemeDescriptionType& scheme_params,
							 typename BaseType::CommunicatorType comm			= typename BaseType::CommunicatorType()
							)
					: BaseType( gridPart, scheme_params, comm )
				{}


				virtual Stuff::RunInfo full_timestep()
				{
					Stuff::Profiler::ScopedTiming fullstep_time("full_step");
					Stuff::RunInfo info;
					for ( int i=0; i < Traits::substep_count; ++i )
					{
						const double dt_k = scheme_params_.step_sizes_[i];
						substep( dt_k, scheme_params_.thetas_[i] );
						if ( i != Traits::substep_count - 1 )
							//the last step increase is done after one call level up
							BaseType::nextStep( i, info );
					}
					return info;
				}

				void substep( const double dt_k, const typename Traits::ThetaSchemeDescriptionType::ThetaValueArray& theta_values )
				{
					//build rhs
					const bool first_step = timeprovider_.timeStep() <= 2;
					const typename Traits::AnalyticalForceType force ( timeprovider_,
																	   currentFunctions_.discreteVelocity().space() ,
																	   viscosity_,
																	   0.0 /*stokes alpha*/ );
					const bool do_cheat = Parameters().getParam( "rhs_cheat", false );

					if ( !Parameters().getParam( "parabolic", false )
							&& ( scheme_params_.algo_id == Traits::ThetaSchemeDescriptionType::scheme_names[3] /*CN*/) )
					{
						typename BaseType::DiscreteVelocityFunctionType beta = currentFunctions_.discreteVelocity();
						beta *= 3.0;
						beta -= lastFunctions_.discreteVelocity();
						beta *= 0.5;
                        Dune::BruteForceReconstruction< typename Traits::OseenModelType >
															::getConvection( beta, rhsDatacontainer_.velocity_gradient, rhsDatacontainer_.convection );
					}
					auto null_f ( exactSolution_.discreteVelocity() );
					null_f *= 0.0;
					boost::scoped_ptr< typename Traits::OseenForceAdapterFunctionType >
							ptr_oseenForceVanilla( first_step //in our very first step no previous computed data is avail. in rhs_container
												? new typename Traits::OseenForceAdapterFunctionType (	timeprovider_,
																										exactSolution_.discreteVelocity(),
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
//					if ( do_cheat )
						BaseType::cheatRHS();
					boost::scoped_ptr< typename Traits::OseenForceAdapterFunctionType >
							ptr_oseenForce( first_step //in our very first step no previous computed data is avail. in rhs_container
												? new typename Traits::OseenForceAdapterFunctionType (	timeprovider_,
																										exactSolution_.discreteVelocity(),
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
					typename BaseType::L2ErrorType::Errors errors_rhs = l2Error_.get(	static_cast<typename Traits::StokesForceAdapterType::BaseType>(*ptr_oseenForce),
																		static_cast<typename Traits::StokesForceAdapterType::BaseType>(*ptr_oseenForceVanilla) );
					Logger().Dbg().Resume(9000);
					Logger().Dbg() << "RHS " << errors_rhs.str() << std::endl;
					typename Traits::DiscreteOseenFunctionWrapperType
							exactSolution_at_next_time ( "reoh", exactSolution_.space(), gridPart_ );
					exactSolution_.atTime( timeprovider_.subTime(), exactSolution_at_next_time  );

					rhsFunctions_.discreteVelocity().assign( *ptr_oseenForce );
					typename Traits::AnalyticalDirichletDataType oseenDirichletData ( timeprovider_,functionSpaceWrapper_ );
					Dune::BetterL2Projection
						::project( timeprovider_, oseenDirichletData, dummyFunctions_.discreteVelocity());
					const double boundaryInt = Stuff::boundaryIntegral( oseenDirichletData, BaseType::currentFunctions().discreteVelocity().space() );
					Logger().Dbg() << boost::format("discrete Boundary integral: %e\n") % boundaryInt;

					unsigned int oseen_iterations = Parameters().getParam( "oseen_iterations", (unsigned int)(1), ValidateGreater<unsigned int>( 0 ) );
					const double dt_n = timeprovider_.deltaT();
					const typename BaseType::L2ErrorType::Errors old_error_velocity
							= l2Error_.get( BaseType::currentFunctions().discreteVelocity(), exactSolution_at_next_time.discreteVelocity() );
					const typename BaseType::L2ErrorType::Errors old_error_pressure
							= l2Error_.get( BaseType::currentFunctions().discretePressure(), exactSolution_at_next_time.discretePressure() );
					double velocity_error_reduction = 1.0;
					double pressure_error_reduction = 1.0;
					unsigned int i = 0;
					do
					{
						const bool do_convection_disc = ! ( Parameters().getParam( "navier_no_convection", false )
															|| Parameters().getParam( "parabolic", false ) );
						bool abort_loop = !do_convection_disc;
						typename BaseType::DiscreteVelocityFunctionType beta = currentFunctions_.discreteVelocity();
						if ( do_convection_disc
								&& ( scheme_params_.algo_id == Traits::ThetaSchemeDescriptionType::scheme_names[3] /*CN*/) )
						{
							beta *= 3.0;
                            auto dummy = lastFunctions_.discreteVelocity();
                            dummy *= 0.5;
                            beta -= dummy;
							abort_loop = true; // linCN only needs a single "iteration"
						}
						else if ( !do_convection_disc )
							beta.clear();

						Dune::StabilizationCoefficients stab_coeff  =
								Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients();
						stab_coeff.FactorFromParams( "C11" );
						stab_coeff.FactorFromParams( "C12" );
						stab_coeff.FactorFromParams( "D11" );
						stab_coeff.FactorFromParams( "D12" );
						typename Traits::OseenModelType
								oseenModel( stab_coeff,
											do_cheat ? *ptr_oseenForce : *ptr_oseenForceVanilla,
											oseenDirichletData,
											theta_values[0] / reynolds_, /*viscosity*/
											1.0f/ dt_n, /*alpha*/
//											do_convection_disc ? theta_values[0] * dt_n : 0.0, /*convection_scale_factor*/
											theta_values[0] , /*convection_scale_factor*/
											theta_values[0] /*pressure_gradient_scale_factor*/
						                   );

                        typename Traits::OseenPassType oseenPass( oseenModel,
												gridPart_,
												functionSpaceWrapper_,
												beta /*beta*/,
												do_convection_disc /*do_oseen_disc*/ );
						if ( timeprovider_.timeStep() <= 2 && i < 1)
							oseenPass.printInfo();
						if ( Parameters().getParam( "silent_stokes", true ) )
							Logger().Info().Suspend( Stuff::Logging::LogStream::default_suspend_priority + 10 );
						//currentFunctions_.clear();
						if ( Parameters().getParam( "clear_u" , true ) )
							nextFunctions_.discreteVelocity().clear();
						if ( Parameters().getParam( "clear_p" , true ) )
							nextFunctions_.discretePressure().clear();
						oseenPass.apply( currentFunctions_, nextFunctions_, &rhsDatacontainer_ );
						Logger().Info().Resume( Stuff::Logging::LogStream::default_suspend_priority + 10 );

						{
							Stuff::Profiler::ScopedTiming error_time("error_calc");
							const typename BaseType::L2ErrorType::Errors new_error_velocity
									= l2Error_.get( nextFunctions_.discreteVelocity(), exactSolution_at_next_time.discreteVelocity() );
							const typename BaseType::L2ErrorType::Errors new_error_pressure
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

						BaseType::setUpdateFunctions();
						currentFunctions_.assign( nextFunctions_ );
						dummyFunctions_.discreteVelocity().assign( rhsDatacontainer_.convection );

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
		};
	}//end namespace NavierStokes
}//end namespace Dune

#endif // METADATA_HH
