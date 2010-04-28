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
									ForceAdapterFunction<TimeProviderType, AnalyticalForceType, VelocityDiscreteFunctionType> >
			{
				protected:
					typedef ForceAdapterFunction<TimeProviderType,AnalyticalForceType, VelocityDiscreteFunctionType>
							ThisType;
					typedef Function< typename AnalyticalForceType::FunctionSpaceType, ThisType >
							BaseType;
					const AnalyticalForceType force_;
					const VelocityDiscreteFunctionType& velocity_;
					const TimeProviderType& timeProvider_;
				public:
					ForceAdapterFunction( const TimeProviderType& timeProvider,
										  const VelocityDiscreteFunctionType& velocity )
						: BaseType( velocity.space() ),
						timeProvider_( timeProvider ),
						force_(0.0 /*visc*/, velocity.space()),
						velocity_( velocity )
					{}

					inline void evaluate( const typename AnalyticalForceType::DomainType& arg,
										  typename AnalyticalForceType::RangeType& ret ) const
					{
						const double time = timeProvider_.time();
						force_.evaluate( time, arg, ret );
						typename VelocityDiscreteFunctionType::JacobianRangeType
								velocity_jacobian;
//						velocity_.jacobian( arg, velocity_jacobian );
						typename VelocityDiscreteFunctionType::RangeType
								velocity_eval;
						velocity_.evaluate( arg, velocity_eval );

						FieldVector<deriType,2> a;
						typename VelocityDiscreteFunctionType::RangeType d;
//						velocity_.evaluate( a, arg, d );
						ret = d;
					}
			};

			template < class TimeProviderType, class AnalyticalDirichletType >
			class DirichletAdapterFunction :
					public Function< typename AnalyticalDirichletType::FunctionSpaceType,
									DirichletAdapterFunction<TimeProviderType, AnalyticalDirichletType > >
			{
				protected:
					typedef DirichletAdapterFunction<TimeProviderType, AnalyticalDirichletType >
							ThisType;
					typedef Function< typename AnalyticalDirichletType::FunctionSpaceType, ThisType >
							BaseType;
					const AnalyticalDirichletType gd_;
					const TimeProviderType& timeProvider_;
				public:
					DirichletAdapterFunction( const TimeProviderType& timeProvider,
										  const typename AnalyticalDirichletType::FunctionSpaceType space )
						: BaseType( space ),
						timeProvider_( timeProvider ),
						gd_(space)
					{}

					template < class IntersectionIteratorType >
					inline void evaluate( const typename AnalyticalDirichletType::DomainType& arg,
										  typename AnalyticalDirichletType::RangeType& ret,
										  const IntersectionIteratorType intIt ) const
					{
						const double time = timeProvider_.time();
						gd_.evaluate( time, arg, ret, intIt );
					}

					inline void evaluate(	const typename AnalyticalDirichletType::DomainType& arg,
											typename AnalyticalDirichletType::RangeType& ret ) const
					{assert(false);}

			};
			template < template < class > class DiricheltDataImp,
						class TimeProviderType >
			struct DirichletAdapterFunctionTraits {

				template < class FunctionSpaceImp, class GridPartImp >
				struct Implementation {
					typedef DirichletAdapterFunction< TimeProviderType, DiricheltDataImp< FunctionSpaceImp > >
						AnalyticalDirichletDataType;

					template <class DiscreteStokesFunctionWrapper >
					static AnalyticalDirichletDataType getInstance( TimeProviderType& timeProvider, const DiscreteStokesFunctionWrapper& wrapper ) {
						return 	AnalyticalDirichletDataType( timeProvider, wrapper.discreteVelocitySpace() );
					}
				};
			};


		} //namespace StokesStep
	}//end namespace NavierStokes
} //end namespace Dune


#endif // RHSADAPTER_HH
