#ifndef EXACTSOLUTION_HH
#define EXACTSOLUTION_HH

#include <dune/stuff/customprojection.hh>

namespace Dune {
	namespace NavierStokes {

		template < class TraitsImp >
		class ExactSolution : public TraitsImp::DiscreteStokesFunctionWrapperType {

				typedef TraitsImp
					TraitsType;
				typedef typename TraitsType::DiscreteStokesFunctionWrapperType
					BaseType;

				const typename TraitsType::TimeProviderType&
						timeprovider_;
				typename TraitsType::PressureFunctionSpaceType
						continousPressureSpace_;
				typename TraitsType::VelocityFunctionSpaceType
						continousVelocitySpace_;
				typename TraitsType::ExactVelocityType
						velocity_;
				typename TraitsType::ExactPressureType
						pressure_;
			public:
				ExactSolution(	const typename TraitsType::TimeProviderType& timeprovider,
								typename TraitsType::GridPartType& gridPart,
								typename TraitsType::DiscreteStokesFunctionSpaceWrapperType& space_wrapper)
					: BaseType( "exact",
								space_wrapper,
								gridPart ),
					timeprovider_( timeprovider ),
					velocity_( timeprovider_, continousVelocitySpace_ ),
					pressure_( timeprovider_, continousPressureSpace_ )
				{
					project();
				}

				void project() {
					BaseType::projectInto( velocity_, pressure_ );
				}

				const typename TraitsType::ExactVelocityType& exactVelocity() const
				{
					return velocity_;
				}

				const typename TraitsType::ExactPressureType& exactPressure() const
				{
					return pressure_;
				}

				typename TraitsType::ExactVelocityType& exactVelocity()
				{
					return velocity_;
				}

				typename TraitsType::ExactPressureType& exactPressure()
				{
					return pressure_;
				}

				void atTime( const double time, BaseType& dest ) const
				{
					Dune::BetterL2Projection
						::project( time, pressure_, dest.discretePressure() );

					Dune::BetterL2Projection
						::project( time, velocity_, dest.discreteVelocity() );
				}

			public:
				typedef typename BaseType::Traits::FunctionTupleType
						FunctionTupleType;

		};

	}//end namespace NavierStokes
}//end namespace Dune


#endif // EXACTSOLUTION_HH
