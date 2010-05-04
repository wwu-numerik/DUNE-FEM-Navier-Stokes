#ifndef EXACTSOLUTION_HH
#define EXACTSOLUTION_HH

namespace Dune {
	namespace NavierStokes {

		template < class ThetaSchemeTraitsImp >
		class ExactSolution : public ThetaSchemeTraitsImp ::DiscreteStokesFunctionWrapperType {

				typedef ThetaSchemeTraitsImp
					ThetaSchemeTraitsType;
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
					project();
				}

				void project() {
					projectInto( velocity_, pressure_ );
				}

				const typename ThetaSchemeTraitsType::ExactVelocityType& exactVelocity() const
				{
					return velocity_;
				}

				const typename ThetaSchemeTraitsType::ExactPressureType& exactPressure() const
				{
					return pressure_;
				}

			public:
				typedef typename BaseType::Traits::FunctionTupleType
						FunctionTupleType;

		};

	}//end namespace NavierStokes
}//end namespace Dune


#endif // EXACTSOLUTION_HH
