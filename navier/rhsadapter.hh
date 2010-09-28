#ifndef RHSADAPTER_HH
#define RHSADAPTER_HH

#include <dune/navier/fractionaltimeprovider.hh>
#include <dune/stuff/printing.hh>
#include <dune/stuff/customprojection.hh>

namespace Dune {
	namespace NavierStokes {
		namespace StokesStep {
			/** \brief take previous step solution and analytical RHS to form function to be passed to either StokesStep
			  * given analytical force \f$f_{ana}\f$ and discrete function \f$u\f$ representing previous time step's velocity solution,
			  *	this calculates new right hand side \f$f := f_{ana} + frac{1}{\theta \tau}u + \frac{\beta}{Re} \Delta u - \left( u \cdot \nabla \right) u \f$
			  *
			  */
			template <	class TimeProviderType,
						class AnalyticalForceType,
						class DiscreteVelocityFunctionType >
			class ForceAdapterFunction :
					public DiscreteVelocityFunctionType
			{
				protected:
					typedef ForceAdapterFunction<	TimeProviderType,
															AnalyticalForceType,
															DiscreteVelocityFunctionType >
						ThisType;
					typedef DiscreteVelocityFunctionType
						BaseType;
					const TimeProviderType& timeProvider_;
					const AnalyticalForceType& force_;
					const double beta_re_qoutient_;
					const double quasi_stokes_alpha_;

				public:
					//! this sginature is used in the first stokes where we have analytical data to derive from
					ForceAdapterFunction(const TimeProviderType& timeProvider,
												   const DiscreteVelocityFunctionType& velocity,
												   const AnalyticalForceType& force,
												   const double beta_re_qoutient,
												   const double quasi_stokes_alpha,
												   int polOrd = -1 )
								 : BaseType( "stokes-ana-rhsdapater" , velocity.space()),
								 timeProvider_( timeProvider ),
								 force_( force ),
								 beta_re_qoutient_( beta_re_qoutient ),
								 quasi_stokes_alpha_( quasi_stokes_alpha )
					{

						//						TESTING_NS
						typedef typename DiscreteVelocityFunctionType::FunctionSpaceType::FunctionSpaceType
							VelocityFunctionSpaceType;
						VelocityFunctionSpaceType continousVelocitySpace_;

						typedef TESTING_NS::VelocityLaplace<	VelocityFunctionSpaceType,
																TimeProviderType >
								VelocityLaplace;
						VelocityLaplace velocity_laplace( timeProvider_, continousVelocitySpace_ );
						typedef TESTING_NS::VelocityConvection<	VelocityFunctionSpaceType,
																TimeProviderType >
							VelocityConvection;
						VelocityConvection velocity_convection( timeProvider_, continousVelocitySpace_ );

						DiscreteVelocityFunctionType velocity_convection_discrete("velocity_convection_discrete", velocity.space() );
						DiscreteVelocityFunctionType velocity_laplace_discrete("velocity_laplace_discrete", velocity.space() );
						DiscreteVelocityFunctionType velocity_copy( "velocity", velocity.space() );
						velocity_copy.assign( velocity );

						Dune::BetterL2Projection
							::project( velocity_convection, velocity_convection_discrete );
						Dune::BetterL2Projection
							::project( velocity_laplace, velocity_laplace_discrete );

						AddCommon( velocity, velocity_convection_discrete, velocity_laplace_discrete );
					}

					//! this signature is used in all other stokes steps where we get the data from the previous step's discretisation
					template < class RhsContainerType >
					ForceAdapterFunction(const TimeProviderType& timeProvider,
												   const DiscreteVelocityFunctionType& velocity,
												   const AnalyticalForceType& force,
												   const double beta_re_qoutient,
												   const double quasi_stokes_alpha,
												   const RhsContainerType& rhs_container,
												   int polOrd = -1 )
								 : BaseType( "stokes-ana-rhsdapater" , velocity.space()),
								 timeProvider_( timeProvider ),
								 force_( force ),
								 beta_re_qoutient_( beta_re_qoutient ),
								 quasi_stokes_alpha_( quasi_stokes_alpha )
					{
						AddCommon( velocity, rhs_container.convection, rhs_container.velocity_laplace );
					}

				protected:
					//! F = f + \alpha / \Re * laplace u + ( 1/(theta * tau) ) u - ( u * nable ) u
					void AddCommon( const DiscreteVelocityFunctionType& velocity,
									const DiscreteVelocityFunctionType& convection,
									const DiscreteVelocityFunctionType& velocity_laplace )
					{
						Dune::BetterL2Projection
							::project( timeProvider_, force_, *this );//this = f

						DiscreteVelocityFunctionType tmp("rhs-ana-tmp", velocity.space() );

						tmp.assign( velocity_laplace );
						tmp *= beta_re_qoutient_;
						*this += tmp;// this = f + beta_re_qoutient * laplace

						tmp.assign( convection );
						*this -= tmp;

						tmp.assign( velocity );
						tmp *= quasi_stokes_alpha_;
						*this += tmp;
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
						class DiscreteVelocityFunctionType >
			class ForceAdapterFunction :
					public DiscreteVelocityFunctionType
			{
				protected:
					typedef ForceAdapterFunction<	TimeProviderType,
													AnalyticalForceType,
													DiscreteVelocityFunctionType >
						ThisType;
					typedef DiscreteVelocityFunctionType
						BaseType;
					const TimeProviderType& timeProvider_;

				public:
					template < class RhsContainerType >
					ForceAdapterFunction( const TimeProviderType& timeProvider,
										  const DiscreteVelocityFunctionType& velocity,
										  const AnalyticalForceType& force,
										  const double alpha_re_qoutient,
										  const double oseen_alpha_unscaled,
										  const RhsContainerType& rhs_container,
										  int polOrd = -1)
						: BaseType( "nonlinear-rhsdapater" , velocity.space()),
						timeProvider_( timeProvider )
					{
						// F = f + \alpha \Re \delta u - \nabla p + ( 1/(1-2 \theta) ) * u
						Dune::BetterL2Projection
							::project( timeProvider_, force, *this );//this = f
						*this -= rhs_container.pressure_gradient;

						DiscreteVelocityFunctionType tmp("nonlinear-rhsdapater-tmp", velocity.space() );
						tmp.assign( rhs_container.velocity_laplace );
						tmp *= alpha_re_qoutient;
						*this += tmp;

						tmp.assign( velocity );
						tmp *= ( oseen_alpha_unscaled );
						*this += tmp;
					}
			};
		}//end namespace NonlinearStep
	}//end namespace NavierStokes
} //end namespace Dune


#endif // RHSADAPTER_HH
