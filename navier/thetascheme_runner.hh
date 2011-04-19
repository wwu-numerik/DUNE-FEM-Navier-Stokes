#ifndef THETASCHEME_RUNNER_HH
#define THETASCHEME_RUNNER_HH

#include <dune/navier/thetascheme_traits.hh>
#include <dune/navier/thetascheme.hh>
#include <dune/navier/thetascheme_alt_split.hh>
#include <dune/navier/global_defines.hh>

template < class GridPartType, class CollectiveCommunicationType >
class ThetaschemeRunner {
	private:
		typedef Dune::NavierStokes::ThetaSchemeTraits<
						CollectiveCommunicationType,
						GridPartType,
						NAVIER_DATA_NAMESPACE::Force,
						NAVIER_DATA_NAMESPACE::DirichletData,
						NAVIER_DATA_NAMESPACE::Pressure,
						NAVIER_DATA_NAMESPACE::Velocity,
						1,//number of substeps
						GridPartType::GridType::dimensionworld,
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
						NAVIER_DATA_NAMESPACE::Force,
						NAVIER_DATA_NAMESPACE::DirichletData,
						NAVIER_DATA_NAMESPACE::Pressure,
						NAVIER_DATA_NAMESPACE::Velocity,
						3,//number of substeps
						GridPartType::GridType::dimensionworld,
						POLORDER,
						VELOCITY_POLORDER,
						PRESSURE_POLORDER >
			ThreeStepThetaSchemeTraitsType;
		typedef Dune::NavierStokes::ThetaScheme<ThreeStepThetaSchemeTraitsType>
			ThreeStepThetaSchemeType;
		typedef Dune::NavierStokes::ThetaSchemeAltSplitting<ThreeStepThetaSchemeTraitsType>
			ThreeStepThetaSchemeAltSplittingType;
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
				case 4: if ( Parameters().getParam("old_timestep", false) )
							return ThreeStepThetaSchemeAltSplittingType(grid_part_,
													  ThreeStepThetaSchemeDescriptionType::fs0( dt_ ) )
									.run();
						else
							return ThreeStepThetaSchemeType(grid_part_,
													  ThreeStepThetaSchemeDescriptionType::fs0( dt_ ) )
									.run();
				default: Logger().Info() << "Using default value for theta scheme type\n";
				case 5: if ( Parameters().getParam("old_timestep", false) )
						return ThreeStepThetaSchemeAltSplittingType(grid_part_,
												  ThreeStepThetaSchemeDescriptionType::fs1( dt_ ) )
								.run();
					else
						return ThreeStepThetaSchemeType(grid_part_,
												  ThreeStepThetaSchemeDescriptionType::fs1( dt_ ) )
								.run();
			}
		}

	private:
		const GridPartType& grid_part_;
		CollectiveCommunicationType& comm_;
};


#endif // THETASCHEME_RUNNER_HH
