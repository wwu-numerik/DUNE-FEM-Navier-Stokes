#ifndef METADATA_HH
#define METADATA_HH

#include <dune/navier/fractionaltimeprovider.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/common/collectivecommunication.hh>
#include <cmath>

namespace Dune {
	namespace NavierStokes {

		enum StepType {
			StokesStepA		= 0,
			NonLinearStep	= 1,
			StokesStepB		= 2,
		};

		template < class Communicator >
		class ThetaScheme {

				Communicator communicator_;
				const double theta_;
				const double operator_weight_alpha_;
				const double operator_weight_beta_;
				const double deltaTime_;
				typedef Dune::NavierStokes::FractionalTimeProvider<Communicator> TPP;
				TPP timeprovider_;
			public:
				ThetaScheme(
							 const double theta = 1 - std::pow( 2.0, -1/2.0 ),
							 Communicator comm = Dune::MPIManager::helper().getCommunicator()
						)
				:theta_(theta),
				operator_weight_alpha_( ( 1-2*theta_ ) / ( 1-theta_ ) ),
				operator_weight_beta_( 1 - operator_weight_alpha_ ),
				communicator_( comm ),
				deltaTime_( Parameters().getParam( "deltaTime", 1e-2 ) ),
				timeprovider_( deltaTime_, theta_, communicator_ )
	//			const double startTime	= Parameters().getParam( "startTime", 0.0 );

				{}


				void dummy () {
					timeprovider_.provideCflEstimate( 1 );
					//not manually setting the delta in tp.nexxt() results in assertions cause TimepRoiver claims dt isn't valid ie unset..
					const double endTime	= Parameters().getParam( "endTime", 1.0 );
					for( timeprovider_.init( deltaTime_ ); timeprovider_.time() < endTime; timeprovider_.next( deltaTime_ ) )
					{
						for ( unsigned int i =0 ; i< 3 ; ++i )
							std::cout << "current time (substep " << i << "): " << timeprovider_.subTime(i) << std::endl;
					}
				}
		};
	}//end namespace NavierStokes
}//end namespace Dune

#endif // METADATA_HH
