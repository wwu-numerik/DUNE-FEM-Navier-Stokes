#ifndef METADATA_HH
#define METADATA_HH

#include <dune/stokes/discretestokesmodelinterface.hh>
#include <dune/stokes/stokespass.hh>
#include <dune/navier/fractionaltimeprovider.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/common/collectivecommunication.hh>
#include <cmath>

namespace Dune {
	namespace NavierStokes {
		template <	class CommunicatorImp,
					class GridPartImp, template <class > class AnalyticalForceImp, class AnalyticalDirichletDataTraits,
					int gridDim, int sigmaOrder, int velocityOrder = sigmaOrder, int pressureOrder = sigmaOrder >
		class ThetaSchemeTraits {
			typedef Dune::DiscreteStokesModelDefaultTraits<
						GridPartImp,
						AnalyticalForceImp,
						AnalyticalDirichletDataTraits,
						gridDim,
						sigmaOrder,
						velocityOrder,
						pressureOrder >
					StokesModelTraits;
			typedef Dune::DiscreteStokesModelDefault< StokesModelTraits >
					StokesModelType;
			typedef StokesModelTraits::DiscreteStokesFunctionSpaceWrapperType
				DiscreteStokesFunctionSpaceWrapperType;

			typedef StokesModelTraits::DiscreteStokesFunctionWrapperType
				DiscreteStokesFunctionWrapperType;
			typedef StokesModelTraits::AnalyticalForceType
				AnalyticalForceType;
			typedef StokesModelTraits::AnalyticalDirichletDataType
				AnalyticalDirichletDataType;

			typedef Dune::StartPass< DiscreteStokesFunctionWrapperType, -1 >
				StartPassType;


			typedef Dune::StokesPass< StokesModelType, StartPassType, 0 >
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
				CommunicatorType communicator_;
				const double theta_;
				const double operator_weight_alpha_;
				const double operator_weight_beta_;
				const double deltaTime_;
				typedef FractionalTimeProvider<Communicator>
						TimeProviderType;
				TimeProviderType timeprovider_;
			public:
				ThetaScheme(
							 const double theta = 1 - std::pow( 2.0, -1/2.0 ),
							 CommunicatorType comm = Dune::MPIManager::helper().getCommunicator()
						)
				:theta_(theta),
				operator_weight_alpha_( ( 1-2*theta_ ) / ( 1-theta_ ) ),
				operator_weight_beta_( 1 - operator_weight_alpha_ ),
				communicator_( comm ),
				deltaTime_( Parameters().getParam( "deltaTime", 1e-2 ) ),
				timeprovider_( deltaTime_, theta_, communicator_ )
	//			const double startTime	= Parameters().getParam( "startTime", 0.0 );

				{
					typename Traits::AnalyticalForceType analyticalForce( viscosity , discreteStokesFunctionSpaceWrapper.discreteVelocitySpace(), alpha );
					typename Traits::AnalyticalDirichletDataType analyticalDirichletData =
							typename Traits::StokesModelTraits::AnalyticalDirichletDataTraitsImplementation::getInstance( discreteStokesFunctionSpaceWrapper );

				}

				void dummy ()
				{
					// function wrapper for the solutions
					typename Traits::DiscreteStokesFunctionSpaceWrapperType
						discreteStokesFunctionSpaceWrapper( gridPart );

					typename Traits::DiscreteStokesFunctionWrapperType
						computedSolutions(  "computed_",
											discreteStokesFunctionSpaceWrapper,
											gridPart );

					timeprovider_.provideCflEstimate( 1 );
					//not manually setting the delta in tp.nexxt() results in assertions cause TimepRoiver claims dt isn't valid ie unset..
					const double endTime	= Parameters().getParam( "endTime", 1.0 );
					for( timeprovider_.init( deltaTime_ ); timeprovider_.time() < endTime; )
					{
						for ( unsigned int i =0 ; i< 3 ; ++i, timeprovider_.nextFractional() )
							std::cout << "current time (substep " << i << "): " << timeprovider_.subTime() << std::endl;
					}
					getStokesPass().apply(computedSolutions,computedSolutions);
				}

				typename Traits::StokesPassType getStokesPass()
				{
					typename Traits::StokesModelType
							stokesModel( Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients() ,
										analyticalForce,
										analyticalDirichletData,
										viscosity,
										alpha );
					StartPassType startPass;
					StokesPassType stokesPass(  startPass,
												stokesModel,
												gridPart,
												discreteStokesFunctionSpaceWrapper );
					return stokesPass;
				}
		};
	}//end namespace NavierStokes
}//end namespace Dune

#endif // METADATA_HH
