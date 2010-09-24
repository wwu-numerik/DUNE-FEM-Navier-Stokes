#ifndef RHSADAPTER_HH
#define RHSADAPTER_HH

#include <dune/navier/fractionaltimeprovider.hh>
#include <dune/stuff/printing.hh>

namespace Dune {

	template < int dim, class RangeType, class JacobianRangeType >
	struct GradientJacobianToLaplacian : public RangeType
	{
		GradientJacobianToLaplacian( const JacobianRangeType& jacobian )
		{
			Dune::CompileTimeChecker< ( dim == 1 || dim > 3 ) > NotImplemented;
		}
	};

	template < class RangeType, class JacobianRangeType >
	struct GradientJacobianToLaplacian< 2, RangeType, JacobianRangeType>  : public RangeType
	{
		GradientJacobianToLaplacian( const JacobianRangeType& jacobian )
		{
			(*this)[0] = jacobian[0][0];
			(*this)[1] = jacobian[3][1];
		}
	};

	template < class RangeType, class JacobianRangeType >
	struct GradientJacobianToLaplacian< 3, RangeType, JacobianRangeType>  : public RangeType
	{
		GradientJacobianToLaplacian( const JacobianRangeType& jacobian )
		{
			(*this)[0] = jacobian[0][0];
			(*this)[1] = jacobian[4][1];
			(*this)[2] = jacobian[8][2];
		}
	};

	template <	class TimeProviderType,
				class DiscreteVelocityFunctionType,
				class SigmaFunctionType >
	class GradientAdapterFunction :
			public SigmaFunctionType
	{
		protected:
			typedef GradientAdapterFunction <	TimeProviderType,
												DiscreteVelocityFunctionType,
												SigmaFunctionType >
				ThisType;
			typedef SigmaFunctionType
				BaseType;
			const TimeProviderType& timeProvider_;

		public:
			GradientAdapterFunction ( const TimeProviderType& timeProvider,
								  const DiscreteVelocityFunctionType& velocity,
								  SigmaFunctionType& dummy,
								  int polOrd = -1 )
				: BaseType( "grad" , dummy.space()),
				timeProvider_( timeProvider )
			{
				typedef SigmaFunctionType
					DiscreteFunctionType;
				typedef typename SigmaFunctionType::DiscreteFunctionSpaceType
					DiscreteFunctionSpaceType;
				typedef typename DiscreteFunctionType::LocalFunctionType
					LocalFuncType;
				typedef typename DiscreteFunctionSpaceType::Traits::GridPartType
					GridPartType;
				typedef typename DiscreteFunctionSpaceType::Traits::IteratorType
					Iterator;
				typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
					BaseFunctionSetType ;
				typedef typename GridPartType::IntersectionIteratorType
					IntersectionIteratorType;
				typedef typename DiscreteVelocityFunctionType::LocalFunctionType
					LocalFType;

				const DiscreteFunctionSpaceType space ( velocity.space().gridPart() );
				const GridPartType& gridPart = space.gridPart();
				// type of quadrature
				typedef Dune::CachingQuadrature<GridPartType,0> VolumeQuadratureType;
				typedef Dune::CachingQuadrature<GridPartType,1> FaceQuadratureType;
				// type of local mass matrix
				typedef Dune::LocalDGMassMatrix< DiscreteFunctionSpaceType, VolumeQuadratureType> LocalMassMatrixType;

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
						const typename DiscreteFunctionSpaceType::DomainType xLocal = quad.point(qP);

						const double intel = (affineMapping) ?
							quad.weight(qP): // affine case
							quad.weight(qP)* geo.integrationElement( xLocal ); // general case

						typename DiscreteFunctionSpaceType::DomainType
							xWorld = geo.global( xLocal );

						// evaluate function
						typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType::RangeType
							velocity_eval;
						velocity_local.evaluate( quad[qP], velocity_eval );

						typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType::JacobianRangeType
							velocity_jacobian_eval;
						velocity_local.jacobian( quad[qP], velocity_jacobian_eval );

						// do projection
						for(int i=0; i<numDofs; ++i)
						{
							typename DiscreteFunctionType::DiscreteFunctionSpaceType::RangeType phi (0.0);
							baseset.evaluate(i, quad[qP], phi);
							self_local[i] += intel * ( Stuff::colonProduct(velocity_jacobian_eval, phi) );
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

	namespace NavierStokes {
		namespace StokesStep {
			/** \brief take previous step solution and analytical RHS to form function to be passed to either StokesStep
			  * given analytical force \f$f_{ana}\f$ and discrete function \f$u\f$ representing previous time step's velocity solution,
			  *	this calculates new right hand side \f$f := f_{ana} + frac{1}{\theta \tau}u + \frac{\beta}{Re} \Delta u - \left( u \cdot \nabla \right) u \f$
			  *
			  */
			template <	class TimeProviderType,
						class AnalyticalForceType,
						class DiscreteVelocityFunctionType,
						class DiscreteVelocityJacobianFunctionType >
			class ForceAdapterFunction :
					public DiscreteVelocityFunctionType
			{
				protected:
					typedef ForceAdapterFunction<	TimeProviderType,
													AnalyticalForceType,
													DiscreteVelocityFunctionType,
													DiscreteVelocityJacobianFunctionType >
						ThisType;
					typedef DiscreteVelocityFunctionType
						BaseType;
					const TimeProviderType& timeProvider_;
					const bool add_extra_terms_;

				public:
					ForceAdapterFunction( const TimeProviderType& timeProvider,
										  const DiscreteVelocityFunctionType& velocity,
										  const AnalyticalForceType& force,
										  const double beta_re_qoutient,
										  const double quasi_stokes_alpha,
										  const bool add_extra_terms,
										  int polOrd = -1 )
						: BaseType( "stokes-rhsdapater" , velocity.space()),
						timeProvider_( timeProvider ),
						add_extra_terms_( add_extra_terms )
					{
						typedef DiscreteVelocityJacobianFunctionType
							DiscreteSigmaFunctionType;

						typedef typename DiscreteSigmaFunctionType::DiscreteFunctionSpaceType
							DiscreteSigmaFunctionSpaceType;

						typedef GradientAdapterFunction<	TimeProviderType,
															DiscreteVelocityFunctionType,
															DiscreteSigmaFunctionType >
							GradientType;
						DiscreteSigmaFunctionSpaceType discrete_sigma_space( velocity.space().gridPart() );
						DiscreteSigmaFunctionType dummy("d", discrete_sigma_space );

						GradientType velocity_jacobian_function ( timeProvider_,
												 velocity,
												 dummy );


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
						typedef CachingQuadrature<GridPartType,0>
							VolumeQuadratureType;
						typedef CachingQuadrature<GridPartType,1>
							FaceQuadratureType;
						// type of local mass matrix
						typedef LocalDGMassMatrix<	DiscreteVelocityFunctionSpaceType,
													VolumeQuadratureType>
							LocalMassMatrixType;

						const int quadOrd = (polOrd == -1) ? (2 * space.order()) : polOrd;

						// create local mass matrix object
						LocalMassMatrixType massMatrix( space, quadOrd );

						// check whether geometry mappings are affine or not
						const bool affineMapping = massMatrix.affine();
						const int out_i = 1;

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

							const typename GradientType::LocalFunctionType& velocity_jacobian_function_local
									= velocity_jacobian_function.localFunction( entity );

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
								velocity_local.evaluate( quad[qP], velocity_eval );

								typename DiscreteVelocityFunctionSpaceType::JacobianRangeType velocity_jacobian_eval;
								velocity_local.jacobian( quad[qP], velocity_jacobian_eval );

								typename DiscreteSigmaFunctionSpaceType::JacobianRangeType velocity_jacobian_function_eval;
								velocity_jacobian_function_local.jacobian( quad[qP], velocity_jacobian_function_eval );
//								Stuff::printFieldMatrix( velocity_jacobian_function_eval, "sjwo", std::cout, "GRA");

								typename AnalyticalForceType::RangeType force_eval;
								force.evaluate( timeProvider_.subTime(), xWorld, force_eval );

								typename DiscreteVelocityFunctionSpaceType::RangeType nonlin;
								for ( unsigned int d = 0; d < nonlin.dim(); ++d ) {
									nonlin[d] = velocity_eval * velocity_jacobian_eval[d];
								}

								GradientJacobianToLaplacian<	DiscreteVelocityFunctionSpaceType::RangeType::size,
																typename DiscreteVelocityFunctionSpaceType::RangeType,
																typename DiscreteSigmaFunctionSpaceType::JacobianRangeType >
									velocity_real_laplacian ( velocity_jacobian_function_eval );


								// do projection
								for(int i=0; i<numDofs; ++i)
								{
									typename DiscreteVelocityFunctionSpaceType::RangeType phi (0.0);
									typename DiscreteVelocityFunctionSpaceType::JacobianRangeType phi_jacobian (0.0);
									baseset.evaluate(i, quad[qP], phi);
									baseset.jacobian( i, quad[qP], phi_jacobian );

									const double force_eval_times_phi = force_eval * phi;
									const double velocity_times_phi = quasi_stokes_alpha * ( velocity_eval * phi );
									const double nonlin_times_phi = nonlin * phi;
									const double velocity_real_laplacian_times_phi = beta_re_qoutient * ( velocity_real_laplacian * phi );

									if ( add_extra_terms_ ) {
										self_local[i] += intel *
														(
															velocity_times_phi
															+ velocity_real_laplacian_times_phi
															+ force_eval_times_phi
															- nonlin_times_phi
														);
									}
									else {
										self_local[i] += intel * force_eval_times_phi;
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
			  *	this calculates new right hand side \f$f := f_{ana} + \frac{\alpha}{Re} \Delta u - \nabla p + \frac{1}{(1-2\theta)dt} u\f$
			  *
			  */
			template <	class TimeProviderType,
						class AnalyticalForceType,
						class DiscreteVelocityFunctionType,
						class DiscretePressureFunctionType,
						class DiscreteVelocityJacobianFunctionType >
			class ForceAdapterFunction :
					public DiscreteVelocityFunctionType
			{
				protected:
					typedef ForceAdapterFunction<	TimeProviderType,
													AnalyticalForceType,
													DiscreteVelocityFunctionType,
													DiscretePressureFunctionType,
													DiscreteVelocityJacobianFunctionType >
						ThisType;
					typedef DiscreteVelocityFunctionType
						BaseType;
					const TimeProviderType& timeProvider_;
					const bool add_extra_terms_;

				public:
					ForceAdapterFunction( const TimeProviderType& timeProvider,
										  const DiscreteVelocityFunctionType& velocity,
										  const DiscretePressureFunctionType& pressure,
										  const AnalyticalForceType& force,
										  const double alpha_re_qoutient,
										  const double u_factor,
										  const bool add_extra_terms,
										  int polOrd = -1)
						: BaseType( "nonlinear-rhsdapater" , velocity.space()),
						timeProvider_( timeProvider ),
						add_extra_terms_( add_extra_terms )
					{
						typedef DiscreteVelocityJacobianFunctionType
							DiscreteSigmaFunctionType;

						typedef typename DiscreteSigmaFunctionType::DiscreteFunctionSpaceType
							DiscreteSigmaFunctionSpaceType;

						typedef GradientAdapterFunction<	TimeProviderType,
															DiscreteVelocityFunctionType,
															DiscreteSigmaFunctionType >
							GradientType;
						DiscreteSigmaFunctionSpaceType discrete_sigma_space( velocity.space().gridPart() );
						DiscreteSigmaFunctionType dummy("d", discrete_sigma_space );

						GradientType velocity_jacobian_function ( timeProvider_,
												 velocity,
												 dummy );

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

							const typename GradientType::LocalFunctionType& velocity_jacobian_function_local
									= velocity_jacobian_function.localFunction( entity );

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

								typename DiscreteSigmaFunctionSpaceType::JacobianRangeType velocity_jacobian_function_eval;
								velocity_jacobian_function_local.jacobian( quad[qP], velocity_jacobian_function_eval );

								GradientJacobianToLaplacian<	DiscreteVelocityFunctionSpaceType::RangeType::size,
																typename DiscreteVelocityFunctionSpaceType::RangeType,
																typename DiscreteSigmaFunctionSpaceType::JacobianRangeType >
									velocity_real_laplacian ( velocity_jacobian_function_eval );

								// do projection
								for(int i=0; i<numDofs; ++i)
								{
									typename DiscreteVelocityFunctionSpaceType::RangeType phi (0.0);
									typename DiscreteVelocityFunctionSpaceType::JacobianRangeType phi_jacobian (0.0);
									baseset.evaluate(i, quad[qP], phi);
									baseset.jacobian(i, quad[qP], phi_jacobian );

									const double force_eval_times_phi = force_eval * phi;
									const double scaled_velocity_times_phi = ( velocity_eval * phi ) * u_factor;
									const double pressure_jacobian_eval_times_phi = pressure_jacobian_eval[0] *  phi;
									const double velocity_real_laplacian_times_phi = alpha_re_qoutient * ( velocity_real_laplacian * phi );

									if ( add_extra_terms_ ) {
										self_local[i] += intel * (
																	velocity_real_laplacian_times_phi
//																   + force_eval_times_phi
//																   - pressure_jacobian_eval_times_phi
//																   + scaled_velocity_times_phi
																) ;
									}
									else
										self_local[i] += intel * force_eval_times_phi;
								}
							}

//							// in case of non-linear mapping apply inverse
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
