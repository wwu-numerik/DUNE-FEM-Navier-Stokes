#ifndef NAVIER_PROBLEMS_TESTCASE_TWOD_HH
#define NAVIER_PROBLEMS_TESTCASE_TWOD_HH

#include <dune/stuff/functions.hh>
#include <dune/stuff/timefunction.hh>
#include <dune/stuff/parametercontainer.hh>

namespace NavierProblems {
namespace TwoDee {

static const std::string identifier = "Testcase2D";
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
	inline void evaluateTime( const double /*time*/, const DomainType& arg, RangeType& ret ) const
	{
		const double lambda			= lambda_;
		const double x				= arg[0];
		const double y				= arg[1];
		const double e_lambda_x		= std::exp( lambda * x );
		const double cos				= std::cos( 2 * M_PI * y );
		const double sin				= std::sin( 2 * M_PI * y );
		const double lambda_square	= lambda * lambda;
		RangeType u;
		VelocityEvaluate( lambda_, 0, arg, u);
		//convection
		//					  assert( false ); //M_2_PI == 2 / PI
		ret[0] = u[0] * ( -lambda * e_lambda_x * cos ) + u[1] * ( -2 * M_PI * e_lambda_x * sin );
		ret[1] = u[0] * ( ( lambda_square / (2 * M_PI)) * e_lambda_x * sin ) + u[1] * lambda * e_lambda_x * cos;
		//laplace
		ret[0] -= viscosity_ * ( - lambda_square * e_lambda_x * cos );
		ret[1] -= viscosity_ * ( - lambda * e_lambda_x * sin );
		//pressure grad
		ret[0] += lambda * std::exp( lambda * 2 * x );
		ret[1] += 0;
		//u
		ret[0] += gamma_ * u[0];
		ret[1] += gamma_ * u[1];
		//					  ret *= 0 ;
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
template < class FunctionSpaceImp , class TimeProviderImp >
class VelocityConvection :  public Dune::TimeFunction < FunctionSpaceImp , VelocityConvection< FunctionSpaceImp,TimeProviderImp >, TimeProviderImp >
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
   *
   *  doing nothing besides Base init
   **/
	VelocityConvection(	const TimeProviderImp& timeprovider,
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
	~VelocityConvection()
	{}

	template < class IntersectionType >
	void evaluate( const double time, const DomainType& arg, RangeType& ret, const IntersectionType& /*intersection */) const
	{
		dune_static_assert( FunctionSpaceImp::dimDomain == 2, "__CLASS__ evaluate not implemented for world dimension");
		VelocityEvaluate( lambda_, time, arg, ret);
	}

	void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
	{VelocityEvaluate( lambda_, 0, arg, ret);}

	/**
   * \brief  evaluates the dirichlet data
   * \param  arg
   *         point to evaluate at
   * \param  ret
   *         value of dirichlet boundary data at given point
   **/
	inline void evaluate( const DomainType& arg, RangeType& ret ) const {VelocityEvaluate( lambda_, 0, arg, ret);}

private:
	static const int dim_ = FunctionSpaceImp::dimDomain ;
	const double lambda_;
};

/**
 *  \brief  describes the dirichlet boundary data
 *
 *  \tparam DirichletTraitsImp
 *          types like functionspace, range type, etc
 *
 *  \todo   extensive docu with latex
 **/
template < class FunctionSpaceImp, class TimeProviderImp >
class DirichletData : public Dune::IntersectionTimeFunction < FunctionSpaceImp , DirichletData< FunctionSpaceImp,TimeProviderImp >, TimeProviderImp >
{
public:
	typedef DirichletData< FunctionSpaceImp, TimeProviderImp >
		ThisType;
	typedef Dune::IntersectionTimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
		BaseType;
	typedef typename BaseType::DomainType
		DomainType;
	typedef typename BaseType::RangeType
		RangeType;

	/**
	  *  \brief  constructor
	  *  \param  viscosity,alpha   dummies
	  **/
	DirichletData( const TimeProviderImp& timeprovider,
				   const FunctionSpaceImp& space,
				   const double /*viscosity*/ = 0.0,
				   const double /*alpha*/ = 0.0 )
		: BaseType ( timeprovider, space )
	{}

	/**
   *  \brief  destructor
   *
   *  doing nothing
   **/
	~DirichletData()
	{}

	template < class IntersectionType >
	void evaluateTime( const double time, const DomainType& arg, RangeType& ret, const IntersectionType& /*intersection */) const
	{
		dune_static_assert( FunctionSpaceImp::dimDomain == 2, "__CLASS__ evaluate not implemented for world dimension");
		VelocityEvaluate( 0.0, time, arg, ret);
	}
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
		dune_static_assert( FunctionSpaceImp::dimDomain == 2, "__CLASS__ evaluate not implemented for world dimension");
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
class Pressure : public Dune::TimeFunction <	FunctionSpaceImp ,
		Pressure < FunctionSpaceImp,TimeProviderImp >,
		TimeProviderImp >
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
		dune_static_assert( FunctionSpaceImp::dimDomain == 2, "__CLASS__ evaluate not implemented for world dimension");
		const double x			= arg[0];
		const double y			= arg[1];
		const double e_2lambda_x= std::exp( 2 * lambda_ * x );

		ret[0] = 0.5 * e_2lambda_x + shift_;
	}

	void setShift( const double shift )
	{
		shift_ = shift;
		Logger().Info() <<  "Set pressure shift to: " << shift_ << std::endl;
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

}//end ns
}//end ns

#endif //NAVIER_PROBLEMS_TESTCASE_TWOD_HH
