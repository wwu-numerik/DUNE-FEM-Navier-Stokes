#ifndef METADATA_HH
#define METADATA_HH

#include <dune/stokes/discretestokesmodelinterface.hh>
#include <dune/stokes/stokespass.hh>
#include <dune/navier/fractionaltimeprovider.hh>
#include <dune/navier/stokestraits.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/common/collectivecommunication.hh>
#include <cmath>

namespace Dune {
	namespace NavierStokes {
		template <	class CommunicatorImp,
					class GridPartImp,
					template <class > class AnalyticalForceImp,
					template <class > class AnalyticalDirichletDataImp,
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
		};

		template < class TraitsImp >
		class ThetaScheme {
			public:
				typedef TraitsImp
						Traits;
				typedef typename Traits::CommunicatorType
						CommunicatorType;
				CommunicatorType& communicator_;
				const double theta_;
				const double operator_weight_alpha_;
				const double operator_weight_beta_;
				const double deltaTime_;
				typename Traits::GridPartType gridPart_;
				typename Traits::DiscreteStokesFunctionSpaceWrapperType functionSpaceWrapper_;

				typename Traits::TimeProviderType timeprovider_;
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
					deltaTime_( Parameters().getParam( "deltaTime", 1e-2 ) ),
		//			const double startTime	= Parameters().getParam( "startTime", 0.0 ),
					timeprovider_( deltaTime_, theta_, communicator_ ),
					functionSpaceWrapper_( gridPart_ )
				{

				}

				void run()
				{
					typename Traits::DiscreteStokesFunctionWrapperType
						currentFunctions(  "current_",
											functionSpaceWrapper_,
											gridPart_ );
					typename Traits::DiscreteStokesFunctionWrapperType
						nextFunctions(  "next_",
											functionSpaceWrapper_,
											gridPart_ );
					const double viscosity = 48102.;
					const double alpha = 48102.;
					typename Traits::AnalyticalForceType stokesForce( timeprovider_, currentFunctions.discreteVelocity() );
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
											stokesModel,//reference goes out-of-scope??
											gridPart_,
											functionSpaceWrapper_ );


					timeprovider_.provideCflEstimate( 1 );
					//not manually setting the delta in tp.nexxt() results in assertions cause TimepRoiver claims dt isn't valid ie unset..
					const double endTime	= Parameters().getParam( "endTime", 1.0 );
					for( timeprovider_.init( deltaTime_ ); timeprovider_.time() < endTime; )
					{
						for ( unsigned int i =0 ; i< 3 ; ++i, timeprovider_.nextFractional() )
							std::cout << "current time (substep " << i << "): " << timeprovider_.subTime() << std::endl;
					}
					stokesPass.apply(currentFunctions,nextFunctions);
				}
		};
	}//end namespace NavierStokes
}//end namespace Dune

#endif // METADATA_HH
