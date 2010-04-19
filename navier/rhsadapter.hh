#ifndef RHSADAPTER_HH
#define RHSADAPTER_HH

#include <dune/navier/fractionaltimeprovider.hh>

namespace Dune {
	namespace NavierStokes {
		namespace StokesStep {
			//! take previous step solution and analytical RHS to form function to be passed in either StokesStep
			template < class TimeProviderType, class AnalyticalForceType, class VelocityDiscreteFunctionType >
			class ForceAdapterFunction :
					public Function< typename AnalyticalForceType::FunctionSpaceType,
									ForceAdapterFunction<AnalyticalForceType, VelocityDiscreteFunctionType> >
			{
				protected:
					typedef ForceAdapterFunction<AnalyticalForceType, VelocityDiscreteFunctionType>
							ThisType;
					typedef Function< typename AnalyticalForceType::FunctionSpaceType, ThisType >
							BaseType;
					const AnalyticalForceType force_;
					const VelocityDiscreteFunctionType& velocity_;
					const TimeProviderType& timeProvider_;
				public:
					ForceAdapterFunction( const TimeProviderType& timeProvider,
										  const VelocityDiscreteFunctionType& velocity )
					: BaseType( velocity.functionSpace() ),
					timeProvider_( timeProvider ),
					force_(sometuff),
					velocity_( velocity )
					{}

			};
		} //namespace StokesStep
	}//end namespace NavierStokes
} //end namespace Dune


#endif // RHSADAPTER_HH
