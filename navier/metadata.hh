#ifndef METADATA_HH
#define METADATA_HH

#include <dune/stokes/discretestokesmodelinterface.hh>
#include <dune/stokes/stokespass.hh>
#include <dune/navier/fractionaltimeprovider.hh>
#include <dune/navier/stokestraits.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/stuff/datawriter.hh>
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
			typedef Dune::DiscreteStokesModelDefault< StokesModelTraits >
				StokesModelType;
			typedef typename StokesModelTraits::DiscreteStokesFunctionSpaceWrapperType
				DiscreteStokesFunctionSpaceWrapperType;

			typedef typename StokesModelTraits::DiscreteStokesFunctionWrapperType
				DiscreteStokesFunctionWrapperType;
			typedef typename StokesModelTraits::AnalyticalForceType
				AnalyticalForceType;
			typedef typename StokesModelTraits::AnalyticalDirichletDataType
				AnalyticalDirichletDataType;

			typedef Dune::StartPass< DiscreteStokesFunctionWrapperType, -1 >
				StokesStartPassType;
			typedef Dune::StokesPass< StokesModelType, StokesStartPassType, 0 >
				StokesPassType;

			typedef CommunicatorImp
				CommunicatorType;

			typedef TimeAwareDataWriter<	TimeProviderType,
											typename GridPartType::GridType,
											typename DiscreteStokesFunctionWrapperType::Traits::FunctionTupleType >
				DataWriterType;
			typedef ExactPressureImp< typename StokesModelTraits::PressureFunctionSpaceType,
									  TimeProviderType >
				ExactPressureType;
			typedef ExactVelocityImp< typename StokesModelTraits::VelocityFunctionSpaceType,
									  TimeProviderType >
				ExactVelocityType;
		};

		template < class ThetaSchemeTraitsType >
		class ExactSolution : public ThetaSchemeTraitsType ::DiscreteStokesFunctionWrapperType {
				typedef typename ThetaSchemeTraitsType ::DiscreteStokesFunctionWrapperType
					BaseType;
				const typename ThetaSchemeTraitsType::TimeProviderType&
						timeprovider_;
				typename ThetaSchemeTraitsType::StokesModelTraits::PressureFunctionSpaceType
						continousPressureSpace_;
				typename ThetaSchemeTraitsType::StokesModelTraits::VelocityFunctionSpaceType
						continousVelocitySpace_;
				const typename ThetaSchemeTraitsType::ExactVelocityType
						velocity_;
				const typename ThetaSchemeTraitsType::ExactPressureType
						pressure_;
			public:
				ExactSolution(	const typename ThetaSchemeTraitsType::TimeProviderType& timeprovider,
								typename ThetaSchemeTraitsType::GridPartType& gridPart,
							  typename ThetaSchemeTraitsType::DiscreteStokesFunctionSpaceWrapperType& space_wrapper)
					: BaseType( "exact",
								space_wrapper,
								gridPart ),
					timeprovider_( timeprovider ),
					velocity_( timeprovider_, continousVelocitySpace_ ),
					pressure_( timeprovider_, continousPressureSpace_ )
				{

				}

		};

		template < class TraitsImp >
		class ThetaScheme {
			protected:
				typedef TraitsImp
						Traits;
				typedef typename Traits::CommunicatorType
						CommunicatorType;
				typedef ExactSolution<Traits>
						ExactSolutionType;
				CommunicatorType& communicator_;
				const double theta_;
				const double operator_weight_alpha_;
				const double operator_weight_beta_;
				typename Traits::GridPartType gridPart_;
				typename Traits::TimeProviderType timeprovider_;
				typename Traits::DiscreteStokesFunctionSpaceWrapperType functionSpaceWrapper_;
				typename Traits::DiscreteStokesFunctionWrapperType currentFunctions_;
				typename Traits::DiscreteStokesFunctionWrapperType nextFunctions_;
				ExactSolutionType exactSolution_;
				typename Traits::DataWriterType dataWriter_;

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
					exactSolution_( timeprovider_,
									gridPart_,
									functionSpaceWrapper_ ),
					dataWriter_( timeprovider_,
								 gridPart_.grid(),
								 nextFunctions_.functionTuple() )
				{}

				void run()
				{
					currentFunctions_.clear(); //actually load initial velocity/pressure
					const double viscosity = 48102.;
					const double alpha = 48102.;
					typename Traits::AnalyticalForceType stokesForce( timeprovider_, currentFunctions_.discreteVelocity() );
					typename Traits::AnalyticalDirichletDataType stokesDirichletData =
							Traits::StokesModelTraits::AnalyticalDirichletDataTraitsImplementation
											::getInstance( timeprovider_,functionSpaceWrapper_ );
					typename Traits::StokesModelType
							stokesModel( Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients() ,
										stokesForce,
										stokesDirichletData,
										viscosity,
										alpha );
					typename Traits::StokesStartPassType stokesStartPass;
					typename Traits::StokesPassType stokesPass( stokesStartPass,
											stokesModel,
											gridPart_,
											functionSpaceWrapper_ );

					for( timeprovider_.init( timeprovider_.deltaT() ); timeprovider_.time() < timeprovider_.endTime(); )
					{
						for ( unsigned int i =0 ; i< 3 ; ++i, timeprovider_.nextFractional() )
							std::cout << "current time (substep " << i << "): " << timeprovider_.subTime() << std::endl;
						dataWriter_.write();
					}
//					stokesPass.apply(currentFunctions_,nextFunctions_);
				}
		};
	}//end namespace NavierStokes
}//end namespace Dune

#endif // METADATA_HH
