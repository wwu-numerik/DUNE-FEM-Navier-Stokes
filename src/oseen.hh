#ifndef OSEEN_HH
#define OSEEN_HH

#include <dune/stokes/discretestokesmodelinterface.hh>
#include <dune/stokes/stokespass.hh>
#include <dune/navier/fractionaltimeprovider.hh>
#include <dune/navier/stokestraits.hh>
#include <dune/navier/exactsolution.hh>
#include <dune/navier/nonlinear/models.hh>
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
				  inline void evaluate( const double /*time*/, const DomainType& /*arg*/, RangeType& ret ) const
				  {
					  ret = RangeType(0);
				  }
				  inline void evaluate( const DomainType& /*arg*/, RangeType& ret ) const {ret = RangeType(0);}

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
					shift_ = -1 * Stuff::meanValue( *this, space );
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

#ifndef OSEEN_DATA_NAMESPACE
	#define OSEEN_DATA_NAMESPACE Oseen::TestCase2D
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
