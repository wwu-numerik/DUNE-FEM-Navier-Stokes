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
					const TimeProviderType& timeProvider_;
					const AnalyticalForceType& force_;
					const VelocityDiscreteFunctionType& velocity_;
				public:
					ForceAdapterFunction( const TimeProviderType& timeProvider,
										  const VelocityDiscreteFunctionType& velocity,
										  const AnalyticalForceType& force )
						: BaseType( velocity.space() ),
						timeProvider_( timeProvider ),
						force_( force ),
						velocity_( velocity )
					{

					}

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
					const TimeProviderType& timeProvider_;
					const AnalyticalDirichletType gd_;
				public:
					DirichletAdapterFunction( const TimeProviderType& timeProvider,
										  const typename AnalyticalDirichletType::FunctionSpaceType space )
						: BaseType( space ),
						timeProvider_( timeProvider ),
						gd_( space )
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
					{
						NEEDS_IMPLEMENTATION
					}

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
		namespace NonlinearStep {
			//! take previous step solution and analytical RHS to form function to be passed in either StokesStep
			template <	class TimeProviderType,
						class AnalyticalForceType,
						class DiscreteVelocityFunctionType,
						class DiscretePressureFunctionType >
			class ForceAdapterFunction :
					public DiscreteVelocityFunctionType
			{
				protected:
					typedef ForceAdapterFunction<	TimeProviderType,
													AnalyticalForceType,
													DiscreteVelocityFunctionType,
													DiscretePressureFunctionType >
						ThisType;
					typedef DiscreteVelocityFunctionType
						BaseType;
					const TimeProviderType& timeProvider_;
//					const AnalyticalForceType force_;
//					const DiscreteVelocityFunctionType& velocity_;
				public:
					ForceAdapterFunction( const TimeProviderType& timeProvider,
										  const DiscreteVelocityFunctionType& velocity,
										  const DiscretePressureFunctionType& pressure,
										  const AnalyticalForceType& force,
										  int polOrd = -1)
						: BaseType( "rhsdapater" , velocity.space()),
						timeProvider_( timeProvider )
//						force_( 0.0 /*visc*/, velocity.space()),
//						velocity_( velocity )
					{
						typedef typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType
							DiscreteFunctionSpaceType;
						typedef typename DiscreteVelocityFunctionType::LocalFunctionType
							LocalFuncType;
						typedef typename DiscreteFunctionSpaceType::Traits::GridPartType
							GridPartType;
						typedef typename DiscreteFunctionSpaceType::Traits::IteratorType
							Iterator;
						typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
							BaseFunctionSetType ;
						typedef typename GridPartType::GridType
							GridType;
						typedef typename DiscreteVelocityFunctionType::LocalFunctionType
							LocalFType;

						typename DiscreteFunctionSpaceType::RangeType ret (0.0);
						typename DiscreteFunctionSpaceType::RangeType phi (0.0);
						const DiscreteFunctionSpaceType& space =  velocity.space();

						// type of quadrature
						typedef CachingQuadrature<GridPartType,0> QuadratureType;
						// type of local mass matrix
						typedef LocalDGMassMatrix< DiscreteFunctionSpaceType, QuadratureType > LocalMassMatrixType;

						const int quadOrd = (polOrd == -1) ? (2 * space.order()) : polOrd;

						// create local mass matrix object
						LocalMassMatrixType massMatrix( space, quadOrd );

						// check whether geometry mappings are affine or not
						const bool affineMapping = massMatrix.affine();

						// clear destination
						BaseType::clear();

						const Iterator endit = space.end();
						for(Iterator it = space.begin(); it != endit ; ++it)
						{
						  // get entity
						  const typename GridType::template Codim<0>::Entity& entity = *it;
						  // get geometry
						  const typename GridType::template Codim<0>::Geometry& geo = entity.geometry();

						  // get quadrature
						  QuadratureType quad(entity, quadOrd);

						  // get local function of destination
						  LocalFuncType self_local = localFunction(entity);
						  // get local function of argument
						  const LocalFType velocity_local = velocity.localFunction(entity);

						  // get base function set
						  const BaseFunctionSetType & baseset = self_local.baseFunctionSet();

						  const int quadNop = quad.nop();
						  const int numDofs = self_local.numDofs();

						  for(int qP = 0; qP < quadNop ; ++qP)
						  {
							const double intel = (affineMapping) ?
								 quad.weight(qP) : // affine case
								 quad.weight(qP) * geo.integrationElement( quad.point(qP) ); // general case

							// evaluate function
							typename DiscreteFunctionSpaceType::RangeType
									dummy;
							velocity_local.evaluate(quad.point(qP), dummy);
							ret =dummy;

							typename DiscreteFunctionSpaceType::DomainType
								xWorld = geo.global( quad.point(qP) );
							typename AnalyticalForceType::RangeType force_eval;
							force.evaluate( xWorld, force_eval );

							// do projection
							for(int i=0; i<numDofs; ++i)
							{
							  baseset.evaluate(i, quad[qP], phi);
							  self_local[i] += intel * (ret * phi) ;
							}
						  }

						  // in case of non-linear mapping apply inverse
						  if ( ! affineMapping )
						  {
							massMatrix.applyInverse( entity, self_local );
						  }
						}
					}

			};
		}//end namespace NonlinearStep
	}//end namespace NavierStokes
} //end namespace Dune


#endif // RHSADAPTER_HH
