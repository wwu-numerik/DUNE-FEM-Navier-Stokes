#ifndef OSEEN_HH
#define OSEEN_HH

#include <dune/stokes/discretestokesmodelinterface.hh>
#include <dune/stokes/stokespass.hh>
#include <dune/navier/fractionaltimeprovider.hh>
#include <dune/navier/stokestraits.hh>
#include <dune/navier/exactsolution.hh>
#include <dune/navier/oseen/oseenpass.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/stuff/datawriter.hh>
#include <dune/stuff/functions.hh>
#include <dune/stuff/timefunction.hh>
#include <dune/stuff/customprojection.hh>
#include <dune/common/collectivecommunication.hh>
#include <cmath>

namespace Dune {
namespace Oseen {

		template <	class TimeProviderType,
					class GridPartImp,
					template < class > class ForceFuntionType,
					template < class > class AnalyticalDirichletDataImp,
					int gridDim, int sigmaOrder, int velocityOrder = sigmaOrder, int pressureOrder = sigmaOrder >
		class DiscreteModelTraits
		{
			public:

				//! for CRTP trick
				typedef DiscreteStokesModelDefault < DiscreteModelTraits >
					DiscreteModelType;

				//! we use caching quadratures for the entities
				typedef Dune::CachingQuadrature< GridPartImp, 0 >
					VolumeQuadratureType;

				//! we use caching quadratures for the faces
				typedef Dune::CachingQuadrature< GridPartImp, 1 >
					FaceQuadratureType;

				//! polynomial order for the discrete sigma function space
				static const int sigmaSpaceOrder = sigmaOrder;
				//! polynomial order for the discrete velocity function space
				static const int velocitySpaceOrder = velocityOrder;
				//! polynomial order for the discrete pressure function space
				static const int pressureSpaceOrder = pressureOrder;

		//    private:

				//! function space type for the velocity
				typedef Dune::FunctionSpace< double, double, gridDim, gridDim >
					VelocityFunctionSpaceType;

				//! discrete function space type for the velocity
				typedef Dune::DiscontinuousGalerkinSpace<   VelocityFunctionSpaceType,
															GridPartImp,
															velocitySpaceOrder >
					DiscreteVelocityFunctionSpaceType;

				//! function space type for the pressure
				typedef Dune::FunctionSpace< double, double, gridDim, 1 >
					PressureFunctionSpaceType;

				//! discrete function space type for the pressure
				typedef Dune::DiscontinuousGalerkinSpace<   PressureFunctionSpaceType,
															GridPartImp,
															pressureSpaceOrder >
					DiscretePressureFunctionSpaceType;

			public:

				//! discrete function space wrapper type
				typedef Dune::DiscreteStokesFunctionSpaceWrapper< Dune::DiscreteStokesFunctionSpaceWrapperTraits<
							DiscreteVelocityFunctionSpaceType,
							DiscretePressureFunctionSpaceType > >
					DiscreteStokesFunctionSpaceWrapperType;

			private:

				//! discrete function type for the velocity
				typedef Dune::AdaptiveDiscreteFunction< typename DiscreteStokesFunctionSpaceWrapperType::DiscreteVelocityFunctionSpaceType >
					DiscreteVelocityFunctionType;

				//! discrete function type for the pressure
				typedef Dune::AdaptiveDiscreteFunction< typename DiscreteStokesFunctionSpaceWrapperType::DiscretePressureFunctionSpaceType >
					DiscretePressureFunctionType;

			public:

				//! discrete function wrapper type
				typedef Dune::DiscreteStokesFunctionWrapper< Dune::DiscreteStokesFunctionWrapperTraits<
							DiscreteStokesFunctionSpaceWrapperType,
							DiscreteVelocityFunctionType,
							DiscretePressureFunctionType > >
					DiscreteStokesFunctionWrapperType;

			private:

				//! function space type for sigma
				typedef Dune::MatrixFunctionSpace<  double,
													double,
													gridDim,
													gridDim,
													gridDim >
					SigmaFunctionSpaceType;

				//! discrete function space type for sigma
				typedef Dune::DiscontinuousGalerkinSpace<   SigmaFunctionSpaceType,
															GridPartImp,
															sigmaSpaceOrder >
					DiscreteSigmaFunctionSpaceType;

			public:

				//! discrete function type for sigma
				typedef Dune::AdaptiveDiscreteFunction< DiscreteSigmaFunctionSpaceType >
					DiscreteSigmaFunctionType;

				//! function type for the analytical force
				typedef ForceFuntionType < VelocityFunctionSpaceType >
					AnalyticalForceFunctionType;

				typedef AnalyticalForceFunctionType
					AnalyticalForceType;

				//! function type for the analytical dirichlet data
				typedef typename Dune::NavierStokes::StokesStep::DirichletAdapterFunctionTraits< AnalyticalDirichletDataImp, TimeProviderType >
									::template Implementation<VelocityFunctionSpaceType,GridPartImp >
						AnalyticalDirichletDataTraitsImplementation;
				typedef typename AnalyticalDirichletDataTraitsImplementation::AnalyticalDirichletDataType
					AnalyticalDirichletDataType;

				typedef DiscreteVelocityFunctionType
					ExtraDataFunctionType;
				/**
				 *  \name   types needed for the pass
				 *  \{
				 **/
				//! return type of the pass
				typedef DiscreteStokesFunctionWrapperType
					DestinationType;
				/**
				 *  \}
				 **/

		};

	namespace TestCase2D {
		template < class FunctionSpaceImp >
		class Force : public Function < FunctionSpaceImp , Force < FunctionSpaceImp > >
		{
			  public:
				  typedef Force< FunctionSpaceImp >
					  ThisType;
				  typedef Function < FunctionSpaceImp ,ThisType >
					  BaseType;
				  typedef typename BaseType::DomainType
					  DomainType;
				  typedef typename BaseType::RangeType
					  RangeType;

				  /**
				   *  \brief  constructor
				   *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
				   **/
				  Force( const double viscosity, const FunctionSpaceImp& space, const double alpha = 0.0 )
					  : BaseType ( space ),
						viscosity_( viscosity ),
						alpha_( alpha ),
						lambda_( Parameters().getParam( "lambda", 0.0 ) ),
						gamma_( Parameters().getParam( "alpha", 0.0 ) )
				  {}

				  /**
				   *  \brief  destructor
				   *  doing nothing
				   **/
				  ~Force()
				  {}

				  /**
				   *  \brief  evaluates the force
				   *  \param  arg
				   *          point to evaluate at
				   *  \param  ret
				   *          value of force at given point
				   **/
				  inline void evaluate( const double /*time*/, const DomainType& arg, RangeType& ret ) const
				  {
					  const double lambda			= lambda_;
					  const double x				= arg[0];
					  const double y				= arg[1];
					  const double e_lambda_x		= std::exp( lambda * x );
					  const double cos				= std::cos( 2 * M_PI * y );
					  const double sin				= std::sin( 2 * M_PI * y );
					  const double lambda_square	= lambda * lambda;
					  ret[0] = gamma_ * ( 1 - e_lambda_x * cos ) + e_lambda_x * (
							  ( 1 - e_lambda_x * cos ) * (-1 * lambda ) * cos
							  + ( std::pow(lambda,3)/ (2 * std::pow(M_PI,2) )  ) * e_lambda_x * std::pow( sin, 2 )
							  - viscosity_ * 2 * lambda * M_PI * sin
							  - e_lambda_x
							  );
					  ret[1] = gamma_ * ( lambda / ( 2 * M_PI ) ) * e_lambda_x * cos
							   + e_lambda_x * sin * (
									( lambda_square / ( 2 * M_PI ) ) * e_lambda_x * cos
									+ ( 1 - e_lambda_x * cos ) * 2 * M_PI
									+ viscosity_ * 2 * lambda * M_PI
									   );
//					  ret = RangeType(0);
				  }
				  inline void evaluate( const DomainType& arg, RangeType& ret ) const {
					  evaluate(0, arg, ret);
				  }

			  private:
				  const double viscosity_;
				  const double alpha_;
				  const double lambda_;
				  const double gamma_;
				  static const int dim_ = FunctionSpaceImp::dimDomain;
		};

		template < class DomainType, class RangeType >
		void VelocityEvaluate( const double lambda, const double time, const DomainType& arg, RangeType& ret)
		{
			const double x				= arg[0];
			const double y				= arg[1];
			const double e_lambda_x		= std::exp( lambda * x );

			ret[0] = 1 - e_lambda_x * 	std::cos( 2 * M_PI * y );
			ret[1] = (lambda/(2*M_PI)) * e_lambda_x * 	std::sin( 2 * M_PI * y );
		}

		/**
		*  \brief  describes the dirichlet boundary data
		*
		*  \tparam DirichletTraitsImp
		*          types like functionspace, range type, etc
		*
		*  \todo   extensive docu with latex
		**/
		template < class FunctionSpaceImp >
		class DirichletData : public Function < FunctionSpaceImp , DirichletData < FunctionSpaceImp > >
		{
			public:
				typedef DirichletData< FunctionSpaceImp >
					ThisType;
				typedef Function< FunctionSpaceImp, ThisType >
					BaseType;
				typedef typename BaseType::DomainType
					DomainType;
				typedef typename BaseType::RangeType
					RangeType;

				/**
				*  \brief  constructor
				*
				*  doing nothing besides Base init
				**/
				DirichletData( const FunctionSpaceImp& space,
							 const double parameter_a = M_PI /2.0 ,
							 const double parameter_d = M_PI /4.0)
					: BaseType( space ),
					lambda_( Parameters().getParam( "lambda", 0.0 ) )
				{}

				/**
				*  \brief  destructor
				*
				*  doing nothing
				**/
				~DirichletData()
				{}

				template < class IntersectionType >
				void evaluate( const double time, const DomainType& arg, RangeType& ret, const IntersectionType& /*intersection */) const
				{
					Dune::CompileTimeChecker< ( dim_ == 2 ) > DirichletData_Unsuitable_WorldDim;
					VelocityEvaluate( lambda_, time, arg, ret);
				}

				/**
				* \brief  evaluates the dirichlet data
				* \param  arg
				*         point to evaluate at
				* \param  ret
				*         value of dirichlet boundary data at given point
				**/
				inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}

			private:
				  static const int dim_ = FunctionSpaceImp::dimDomain ;
				  const double lambda_;
		};

		template < class FunctionSpaceImp, class TimeProviderImp >
		class Velocity : public TimeFunction < FunctionSpaceImp , Velocity< FunctionSpaceImp,TimeProviderImp >, TimeProviderImp >
		{
			public:
				typedef Velocity< FunctionSpaceImp, TimeProviderImp >
					ThisType;
				typedef TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
					BaseType;
				typedef typename BaseType::DomainType
					DomainType;
				typedef typename BaseType::RangeType
					RangeType;

				/**
				*  \brief  constructor
				*
				*  doing nothing besides Base init
				**/
				Velocity(	const TimeProviderImp& timeprovider,
							const FunctionSpaceImp& space,
							const double parameter_a = M_PI /2.0 ,
							const double parameter_d = M_PI /4.0)
					: BaseType( timeprovider, space ),
					lambda_( Parameters().getParam( "lambda", 0.0 ) )
				{}

				/**
				*  \brief  destructor
				*
				*  doing nothing
				**/
				~Velocity()
				{}

				void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
				{
					Dune::CompileTimeChecker< ( dim_ == 2 ) > DirichletData_Unsuitable_WorldDim;
					VelocityEvaluate( lambda_, time, arg, ret);
				}

			   /**
				* \brief  evaluates the dirichlet data
				* \param  arg
				*         point to evaluate at
				* \param  ret
				*         value of dirichlet boundary data at given point
				**/
//					inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}

			private:
				static const int dim_ = FunctionSpaceImp::dimDomain ;
				const double lambda_;
		};

		template <	class FunctionSpaceImp,
					class TimeProviderImp >
		class Pressure : public TimeFunction <	FunctionSpaceImp ,
												Pressure < FunctionSpaceImp,TimeProviderImp >,
												TimeProviderImp >
		{
			public:
				typedef Pressure< FunctionSpaceImp, TimeProviderImp >
					ThisType;
				typedef TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
					BaseType;
				typedef typename BaseType::DomainType
					DomainType;
				typedef typename BaseType::RangeType
					RangeType;

			  /**
			   *  \brief  constructor
			   *
			   *  doing nothing besides Base init
			   **/
			  Pressure( const TimeProviderImp& timeprovider,
						const FunctionSpaceImp& space,
						const double parameter_a = M_PI /2.0 ,
						const double parameter_d = M_PI /4.0)
				  : BaseType( timeprovider, space ),
				  lambda_( Parameters().getParam( "lambda", 0.0 ) ),
				  shift_(0.0)
			  {}

			  /**
			   *  \brief  destructor
			   *
			   *  doing nothing
			   **/
			   ~Pressure()
			   {}

				void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
				{
					Dune::CompileTimeChecker< ( dim_ == 2 ) > Pressure_Unsuitable_WorldDim;
					const double x			= arg[0];
					const double y			= arg[1];
					const double e_2lambda_x= std::exp( 2 * lambda_ * x );

					ret[0] = 0.5 * e_2lambda_x + shift_;
				}

				template < class DiscreteFunctionSpace >
				void setShift( const DiscreteFunctionSpace& space )
				{
					shift_ = Stuff::meanValue( *this, space );
					Logger().Info() <<  "mean pressure pre-shift: " << shift_ << std::endl;
					shift_ = 0;//std::sqrt( shift_ / 10) ;
				}

				/**
				* \brief  evaluates the dirichlet data
				* \param  arg
				*         point to evaluate at
				* \param  ret
				*         value of dirichlet boundary data at given point
				**/
//					inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}

			private:
				static const int dim_ = FunctionSpaceImp::dimDomain ;
				const double lambda_;
				double shift_;
		};

	} //end namespace Oseem::TestCase2D

	namespace TestCaseTaylor2D {
		template < class FunctionSpaceImp >
		class Force : public Function < FunctionSpaceImp , Force < FunctionSpaceImp > >
		{
			  public:
				  typedef Force< FunctionSpaceImp >
					  ThisType;
				  typedef Function < FunctionSpaceImp ,ThisType >
					  BaseType;
				  typedef typename BaseType::DomainType
					  DomainType;
				  typedef typename BaseType::RangeType
					  RangeType;

				  /**
				   *  \brief  constructor
				   *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
				   **/
				  Force( const double viscosity, const FunctionSpaceImp& space, const double alpha = 0.0 )
					  : BaseType ( space ),
						viscosity_( viscosity ),
						alpha_( alpha )
				  {}

				  /**
				   *  \brief  destructor
				   *  doing nothing
				   **/
				  ~Force()
				  {}

				  /**
				   *  \brief  evaluates the force
				   *  \param  arg
				   *          point to evaluate at
				   *  \param  ret
				   *          value of force at given point
				   **/
//				  inline void evaluate( const double /*time*/, const DomainType& arg, RangeType& ret ) const
//				  {
//						const double x				= arg[0];
//						const double y				= arg[1];
//						ret[0] = M_PI * std::cos( M_PI * x ) * std::sin( M_PI * y );
//						ret[1] = M_PI * std::sin( M_PI * x ) * std::cos( M_PI * y );
//						ret[0]= 4230;
//						ret[1]=4213;
//				  }
				  inline void evaluate( const DomainType& arg, RangeType& ret ) const
				  {
					  const double x			= arg[0];
					  const double y			= arg[1];
					  const double v			= viscosity_;
					  const double P			= M_PI;//pi_factor;
					  const double E			= 1;//std::exp( -2 * std::pow( P, 2 ) * viscosity_ * time );
					  const double F			= 1;//std::exp( -4 * std::pow( P, 2 ) * viscosity_ * time );
					  const double S_x			= std::sin( P * x );
					  const double S_y			= std::sin( P * y );
					  const double S_2x			= std::sin( 2 * P * x );
					  const double S_2y			= std::sin( 2 * P * y );
					  const double C_x			= std::cos( P * x );
					  const double C_y			= std::cos( P * y );
					  ret[0] = - C_x * E * P * ( S_x * E + v * S_y * P )	+ 0.5 * P * F * S_2x;
					  ret[1] = - C_y * E * P * ( S_y * E - v * S_x * P )	+ 0.5 * P * F * S_2y;
				  }

			  private:
				  const double viscosity_;
				  const double alpha_;
				  static const int dim_ = FunctionSpaceImp::dimDomain;
		};

		template < class DomainType, class RangeType >
		void VelocityEvaluate( const double lambda, const double time, const DomainType& arg, RangeType& ret)
		{
			const double x				= arg[0];
			const double y				= arg[1];
			const double pi_factor		= M_PI;
			const double v				= Parameters().getParam( "viscosity", 1.0 );
			const double e_x			= std::exp( -2 * std::pow(M_PI,2) * v * time );

			const double e_minus_2_t	= std::exp( -2 * std::pow( pi_factor, 2 ) * v * time );

			ret[0] = -1 *	std::cos( pi_factor * x ) * std::sin( pi_factor * y ) * e_minus_2_t;
			ret[1] =		std::sin( pi_factor * x ) * std::cos( pi_factor * y ) * e_minus_2_t;		}

		/**
		*  \brief  describes the dirichlet boundary data
		*
		*  \tparam DirichletTraitsImp
		*          types like functionspace, range type, etc
		*
		*  \todo   extensive docu with latex
		**/
		template < class FunctionSpaceImp >
		class DirichletData : public Function < FunctionSpaceImp , DirichletData < FunctionSpaceImp > >
		{
			public:
				typedef DirichletData< FunctionSpaceImp >
					ThisType;
				typedef Function< FunctionSpaceImp, ThisType >
					BaseType;
				typedef typename BaseType::DomainType
					DomainType;
				typedef typename BaseType::RangeType
					RangeType;

				/**
				*  \brief  constructor
				*
				*  doing nothing besides Base init
				**/
				DirichletData( const FunctionSpaceImp& space,
							 const double parameter_a = M_PI /2.0 ,
							 const double parameter_d = M_PI /4.0)
					: BaseType( space ),
					lambda_( Parameters().getParam( "lambda", 0.0 ) )
				{}

				/**
				*  \brief  destructor
				*
				*  doing nothing
				**/
				~DirichletData()
				{}

				template < class IntersectionType >
				void evaluate( const double time, const DomainType& arg, RangeType& ret, const IntersectionType& /*intersection */) const
				{
					Dune::CompileTimeChecker< ( dim_ == 2 ) > DirichletData_Unsuitable_WorldDim;
					VelocityEvaluate( lambda_, time, arg, ret);
				}

				/**
				* \brief  evaluates the dirichlet data
				* \param  arg
				*         point to evaluate at
				* \param  ret
				*         value of dirichlet boundary data at given point
				**/
				inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}

			private:
				  static const int dim_ = FunctionSpaceImp::dimDomain ;
				  const double lambda_;
		};

		template < class FunctionSpaceImp, class TimeProviderImp >
		class Velocity : public TimeFunction < FunctionSpaceImp , Velocity< FunctionSpaceImp,TimeProviderImp >, TimeProviderImp >
		{
			public:
				typedef Velocity< FunctionSpaceImp, TimeProviderImp >
					ThisType;
				typedef TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
					BaseType;
				typedef typename BaseType::DomainType
					DomainType;
				typedef typename BaseType::RangeType
					RangeType;

				/**
				*  \brief  constructor
				*
				*  doing nothing besides Base init
				**/
				Velocity(	const TimeProviderImp& timeprovider,
							const FunctionSpaceImp& space,
							const double parameter_a = M_PI /2.0 ,
							const double parameter_d = M_PI /4.0)
					: BaseType( timeprovider, space ),
					lambda_( Parameters().getParam( "lambda", 0.0 ) )
				{}

				/**
				*  \brief  destructor
				*
				*  doing nothing
				**/
				~Velocity()
				{}

				void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
				{
					Dune::CompileTimeChecker< ( dim_ == 2 ) > DirichletData_Unsuitable_WorldDim;
					VelocityEvaluate( lambda_, time, arg, ret);
				}

			   /**
				* \brief  evaluates the dirichlet data
				* \param  arg
				*         point to evaluate at
				* \param  ret
				*         value of dirichlet boundary data at given point
				**/
//					inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}

			private:
				static const int dim_ = FunctionSpaceImp::dimDomain ;
				const double lambda_;
		};

		template <	class FunctionSpaceImp,
					class TimeProviderImp >
		class Pressure : public TimeFunction <	FunctionSpaceImp ,
												Pressure < FunctionSpaceImp,TimeProviderImp >,
												TimeProviderImp >
		{
			public:
				typedef Pressure< FunctionSpaceImp, TimeProviderImp >
					ThisType;
				typedef TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
					BaseType;
				typedef typename BaseType::DomainType
					DomainType;
				typedef typename BaseType::RangeType
					RangeType;

			  /**
			   *  \brief  constructor
			   *
			   *  doing nothing besides Base init
			   **/
			  Pressure( const TimeProviderImp& timeprovider,
						const FunctionSpaceImp& space,
						const double parameter_a = M_PI /2.0 ,
						const double parameter_d = M_PI /4.0)
				  : BaseType( timeprovider, space ),
				  lambda_( Parameters().getParam( "lambda", 0.0 ) ),
				  shift_(0.0)
			  {}

			  /**
			   *  \brief  destructor
			   *
			   *  doing nothing
			   **/
			   ~Pressure()
			   {}

				void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
				{
					Dune::CompileTimeChecker< ( dim_ == 2 ) > Pressure_Unsuitable_WorldDim;
					const double x				= arg[0];
					const double y				= arg[1];
					const double v				= Parameters().getParam( "viscosity", 1.0 );
					const double pi_factor		= M_PI;
					const double e_minus_4_t	= std::exp( -4 * std::pow( pi_factor, 2 ) * time * v );

					ret[0] = -0.25 * (
										std::cos( 2 * pi_factor * x ) + std::cos( 2 * pi_factor * y )
									) * e_minus_4_t;
				}

				template < class DiscreteFunctionSpace >
				void setShift( const DiscreteFunctionSpace& space )
				{
//					shift_ = -1 * Stuff::meanValue( *this, space );
				}

				/**
				* \brief  evaluates the dirichlet data
				* \param  arg
				*         point to evaluate at
				* \param  ret
				*         value of dirichlet boundary data at given point
				**/
//					inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}

			private:
				static const int dim_ = FunctionSpaceImp::dimDomain ;
				const double lambda_;
				double shift_;
		};
	}//end namespace TestCaseTaylor2D

#ifndef OSEEN_DATA_NAMESPACE
	#define OSEEN_DATA_NAMESPACE Oseen::TestCaseTaylor2D
#endif

	template <	class CommunicatorImp,
				class GridPartImp,
				int gridDim, int sigmaOrder, int velocityOrder = sigmaOrder, int pressureOrder = sigmaOrder >
	struct Traits {
		typedef Traits<	CommunicatorImp,
						GridPartImp,
						gridDim, sigmaOrder,
						velocityOrder, pressureOrder >
			ThisType;
		typedef GridPartImp
			GridPartType;
		typedef Dune::NavierStokes::FractionalTimeProvider<CommunicatorImp>
			TimeProviderType;

		typedef DiscreteModelTraits<
					TimeProviderType,
					GridPartType,
					OSEEN_DATA_NAMESPACE::Force,
					OSEEN_DATA_NAMESPACE::DirichletData,
					gridDim,
					sigmaOrder,
					velocityOrder,
					pressureOrder >
			OseenModelTraits;
		typedef Dune::DiscreteStokesModelDefault< OseenModelTraits >
			OseenModelType;
		typedef Dune::StartPass< typename OseenModelTraits::DiscreteStokesFunctionWrapperType, -1 >
			StartPassType;
		typedef Dune::NavierStokes::NonlinearStep::OseenPass< OseenModelType,StartPassType >
			OseenPassType;

		typedef OSEEN_DATA_NAMESPACE::Pressure< typename OseenModelTraits::PressureFunctionSpaceType,
								  TimeProviderType >
			ExactPressureType;
		typedef OSEEN_DATA_NAMESPACE::Velocity< typename OseenModelTraits::VelocityFunctionSpaceType,
								  TimeProviderType >
			ExactVelocityType;
		typedef Dune::NavierStokes::ExactSolution<ThisType>
			ExactSolutionType;
		typedef typename OseenModelType::DiscreteStokesFunctionWrapperType
			DiscreteStokesFunctionWrapperType;
		typedef typename OseenModelType::DiscreteStokesFunctionSpaceWrapperType
			DiscreteStokesFunctionSpaceWrapperType;
		typedef typename DiscreteStokesFunctionSpaceWrapperType::DiscretePressureFunctionSpaceType::FunctionSpaceType
			PressureFunctionSpaceType;
		typedef typename DiscreteStokesFunctionSpaceWrapperType::DiscreteVelocityFunctionSpaceType::FunctionSpaceType
			VelocityFunctionSpaceType;


	};

} //end namespace Oseen
}//end namespace Dune

#endif // OSEEN_HH
