#ifndef NAVIER_PROBLEMS_TRIVIAL_HH
#define NAVIER_PROBLEMS_TRIVIAL_HH

#include <dune/stuff/functions.hh>
#include <dune/stuff/timefunction.hh>
#include <dune/stuff/parametercontainer.hh>

namespace NavierProblems {
namespace Trivial {

static const std::string identifier = "Trvial";
static const bool hasExactSolution	= true;

template < class FunctionSpaceImp, class TimeProviderImp >
class Force : public Dune::TimeFunction < FunctionSpaceImp , Force< FunctionSpaceImp,TimeProviderImp >, TimeProviderImp >
{
	public:
		typedef Force< FunctionSpaceImp, TimeProviderImp >
			ThisType;
		typedef Dune::TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
			BaseType;
		  typedef typename BaseType::DomainType
			  DomainType;
		  typedef typename BaseType::RangeType
			  RangeType;

		  /**
		   *  \brief  constructor
		   *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
		   **/
		  Force( const TimeProviderImp& timeprovider, const FunctionSpaceImp& space, const double  viscosity = 0.0, const double alpha = 0.0 )
			  : BaseType ( timeprovider, space ),
				viscosity_( viscosity ),
				alpha_( alpha )
		  {}

		  ~Force()
		  {}


		  void evaluateTime( const double /*time*/, const DomainType& arg, RangeType& ret ) const
		  {
			  const double x			= arg[0];
			  const double y			= arg[1];
			  ret = RangeType( 0 );
			  //pressure gradient
			  ret[0] = 2*x;
			  ret[1] = -2*y;
			  //conv
			  ret[0] += -x;
			  ret[1] += y;
		  }

	  private:
		  const double viscosity_;
		  const double alpha_;
		  static const int dim_ = FunctionSpaceImp::dimDomain;
};

template < class DomainType, class RangeType >
void VelocityEvaluate( const double /*lambda*/, const double /*time*/, const DomainType& arg, RangeType& ret)
{
	const double x				= arg[0];
	const double y				= arg[1];
	ret[0] = -y;
	ret[1] = x;
}

template < class FunctionSpaceImp, class TimeProviderImp >
class VelocityConvection : public Dune::TimeFunction < FunctionSpaceImp , VelocityConvection< FunctionSpaceImp,TimeProviderImp >, TimeProviderImp >
{
	public:
		typedef VelocityConvection< FunctionSpaceImp, TimeProviderImp >
			ThisType;
		typedef Dune::TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
			BaseType;
		  typedef typename BaseType::DomainType
			  DomainType;
		  typedef typename BaseType::RangeType
			  RangeType;

		  /**
		   *  \brief  constructor
		   *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
		   **/
		  VelocityConvection( const TimeProviderImp& timeprovider, const FunctionSpaceImp& space, const double  viscosity = 0.0, const double alpha = 0.0 )
			  : BaseType ( timeprovider, space ),
				viscosity_( viscosity ),
				alpha_( alpha )
		  {}

		  ~VelocityConvection()
		  {}


		  void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
		  {
			  VelocityEvaluate( 0.0f, time, arg, ret );
		  }

	  private:
		  const double viscosity_;
		  const double alpha_;
		  static const int dim_ = FunctionSpaceImp::dimDomain;
};

template < class FunctionSpaceImp >
class DirichletData : public Dune::Function < FunctionSpaceImp , DirichletData < FunctionSpaceImp > >
{
	public:
		typedef DirichletData< FunctionSpaceImp >
			ThisType;
		typedef Dune::Function< FunctionSpaceImp, ThisType >
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
class Velocity : public Dune::TimeFunction < FunctionSpaceImp , Velocity< FunctionSpaceImp,TimeProviderImp >, TimeProviderImp >
{
	public:
		typedef Velocity< FunctionSpaceImp, TimeProviderImp >
			ThisType;
		typedef Dune::TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
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

template <	class FunctionSpaceImp, class TimeProviderImp >
class Pressure : public Dune::TimeFunction < FunctionSpaceImp , Pressure < FunctionSpaceImp,TimeProviderImp >,TimeProviderImp >
{
	public:
		typedef Pressure< FunctionSpaceImp, TimeProviderImp >
			ThisType;
		typedef Dune::TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
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
			const double F				= std::exp( -16 * std::pow( M_PI, 2 ) * time );
			const double C1				= std::cos(4*M_PI* ( x + 0.25 ) );
			const double C2				= std::cos(4*M_PI* ( y + 0.5 ) );

			ret = x*x - y*y;
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

NULLFUNCTION_TP(VelocityLaplace)
NULLFUNCTION_TP(PressureGradient)

}//end ns
}//end ns

#endif //NAVIER_PROBLEMS_TRIVIAL_HH
