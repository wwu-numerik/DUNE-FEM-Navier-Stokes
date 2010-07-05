#ifndef RHSADAPTER_HH
#define RHSADAPTER_HH

#include <dune/navier/fractionaltimeprovider.hh>

namespace Dune {
	namespace NavierStokes {
		namespace StokesStep {
			/** \brief take previous step solution and analytical RHS to form function to be passed to either StokesStep
			  * given analytical force \f$f_{ana}\f$ and discrete function \f$u\f$ representing previous time step's velocity solution,
			  *	this calculates new right hand side \f$f := f_{ana} + \frac{\beta}{Re} \Delta u - \left( u \cdot \nabla \right) u \f$
			  *
			  */
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
//						typename VelocityDiscreteFunctionType::JacobianRangeType
//								velocity_jacobian;
////						velocity_.jacobian( arg, velocity_jacobian );
//						typename VelocityDiscreteFunctionType::RangeType
//								velocity_eval;
//						velocity_.evaluate( arg, velocity_eval );
//
//						FieldVector<deriType,2> a;
//						typename VelocityDiscreteFunctionType::RangeType d;
////						velocity_.evaluate( a, arg, d );
//						ret = d;
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
			/** \brief take previous step solution and analytical RHS to form function to be passed to localdg code
			  * given analytical force \f$f_{ana}\f$ and discrete functions \f$u,p\f$ representing previous time step's velocity and pressure solution,
			  *	this calculates new right hand side \f$f := f_{ana} + \frac{\alpha}{Re} \Delta u - \nabla p\f$
			  *
			  */
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
				public:
					ForceAdapterFunction( const TimeProviderType& timeProvider,
										  const DiscreteVelocityFunctionType& velocity,
										  const DiscretePressureFunctionType& pressure,
										  const AnalyticalForceType& force,
										  const double alpha_re_qoutient,
										  int polOrd = -1)
						: BaseType( "rhsdapater" , velocity.space()),
						timeProvider_( timeProvider )
					{
						typedef typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType
							DiscreteVelocityFunctionSpaceType;
						typedef typename DiscreteVelocityFunctionType::LocalFunctionType
							LocalFuncType;
						typedef typename DiscreteVelocityFunctionSpaceType::Traits::GridPartType
							GridPartType;
						typedef typename DiscreteVelocityFunctionSpaceType::Traits::IteratorType
							Iterator;
						typedef typename DiscreteVelocityFunctionSpaceType::BaseFunctionSetType
							BaseFunctionSetType ;
						typedef typename GridPartType::IntersectionIteratorType
							IntersectionIteratorType;
						typedef typename DiscreteVelocityFunctionType::LocalFunctionType
							LocalFType;

						typename DiscreteVelocityFunctionSpaceType::RangeType ret (0.0);

						const DiscreteVelocityFunctionSpaceType& space =  velocity.space();
						const GridPartType& gridPart = space.gridPart();
						// type of quadrature
						typedef CachingQuadrature<GridPartType,0> VolumeQuadratureType;
						typedef CachingQuadrature<GridPartType,1> FaceQuadratureType;
						// type of local mass matrix
						typedef LocalDGMassMatrix< DiscreteVelocityFunctionSpaceType, VolumeQuadratureType> LocalMassMatrixType;

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
							typedef  typename GridType::template Codim<0>::Geometry
							Geometry;
							const Geometry& geo = entity.geometry();

							// get quadrature
							VolumeQuadratureType quad(entity, quadOrd);

							// get local function of destination
							LocalFuncType self_local = BaseType::localFunction(entity);
							// get local function of argument
							const LocalFType velocity_local = velocity.localFunction(entity);

							// get base function set
							const BaseFunctionSetType & baseset = self_local.baseFunctionSet();

							const int quadNop = quad.nop();
							const int numDofs = self_local.numDofs();

							//volume part
							for(int qP = 0; qP < quadNop ; ++qP)
							{
								const typename DiscreteVelocityFunctionSpaceType::DomainType xLocal = quad.point(qP);

								const double intel = (affineMapping) ?
								quad.weight(qP): // affine case
								quad.weight(qP)* geo.integrationElement( xLocal ); // general case

								typename DiscreteVelocityFunctionSpaceType::DomainType
								xWorld = geo.global( xLocal );

								// evaluate function
								typename DiscreteVelocityFunctionSpaceType::RangeType
								velocity_eval;
								velocity_local.evaluate( xLocal, velocity_eval );

								typename DiscreteVelocityFunctionSpaceType::JacobianRangeType
								velocity_jacobian_eval;
								velocity_local.jacobian( xLocal, velocity_jacobian_eval );

								typename AnalyticalForceType::RangeType force_eval;
								force.evaluate( timeProvider_.subTime(), xWorld, force_eval );

								typename DiscretePressureFunctionType::JacobianRangeType pressure_jacobian_eval;
								pressure.localFunction( entity ).jacobian( xLocal, pressure_jacobian_eval );

								// do projection
								for(int i=0; i<numDofs; ++i)
								{
									typename DiscreteVelocityFunctionSpaceType::RangeType phi (0.0);
									typename DiscreteVelocityFunctionSpaceType::JacobianRangeType phi_jacobian (0.0);
									baseset.evaluate(i, quad[qP], phi);
									baseset.jacobian(i, quad[qP], phi_jacobian );
									const double velocity_jacobian_eval_times_phi_jacobian = alpha_re_qoutient * Stuff::colonProduct( velocity_jacobian_eval, phi_jacobian  );
									const double force_eval_times_phi = force_eval * phi;
									double pressure_jacobian_eval_times_phi = pressure_jacobian_eval[0] *  phi;
									self_local[i] += intel * ( velocity_jacobian_eval_times_phi_jacobian  + force_eval_times_phi - pressure_jacobian_eval_times_phi ) ;
								}
							}

							//surface part
							IntersectionIteratorType intItEnd = gridPart.iend( *it );
							for (   IntersectionIteratorType intIt = gridPart.ibegin( *it );
									intIt != intItEnd;
									++intIt )
							{

								FaceQuadratureType faceQuadrature( gridPart,
														 *intIt,
														 ( 4 * space.order() ) + 1,
														 FaceQuadratureType::INSIDE );
								for ( size_t qP = 0; qP < faceQuadrature.nop(); ++qP )
								{
									const typename DiscreteVelocityFunctionSpaceType::DomainType x = faceQuadrature.point(qP);

									const typename FaceQuadratureType::LocalCoordinateType xLocal = faceQuadrature.localPoint( qP );
									const double intel =
									  faceQuadrature.weight(qP) * intIt->intersectionGlobal().integrationElement( xLocal ); // general case
									const typename DiscreteVelocityFunctionSpaceType::RangeType outerNormal = intIt->unitOuterNormal( xLocal );
									typename DiscreteVelocityFunctionSpaceType::JacobianRangeType extra_u_jacobian;
									velocity_local.jacobian( x, extra_u_jacobian );
									typename DiscreteVelocityFunctionSpaceType::RangeType extra_u_jacobian_times_normal;
									extra_u_jacobian.mv( outerNormal, extra_u_jacobian_times_normal );
									for(int i=0; i<baseset.numBaseFunctions(); ++i)
									{
										typename DiscreteVelocityFunctionSpaceType::RangeType phi (0.0);
										baseset.evaluate(i, faceQuadrature[qP], phi);
										const double extra_u_jacobian_times_normal_times_phi = alpha_re_qoutient * ( extra_u_jacobian_times_normal * phi );
										self_local[i] += intel * extra_u_jacobian_times_normal_times_phi ;
									}
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
