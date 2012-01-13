#ifndef NAVIER_PROBLEMS_TIMEDISC_HH
#define NAVIER_PROBLEMS_TIMEDISC_HH

#include <dune/stuff/functions.hh>
#include <dune/stuff/timefunction.hh>
#include <dune/stuff/parametercontainer.hh>
#include "common.hh"

namespace NavierProblems {
namespace TimeDisc {

static const std::string identifier = "TimeDisc";
static const bool hasExactSolution	= true;

ALLGOOD_SETUPCHECK;

template < class DomainType, class RangeType >
static void evaluateTimeVelocity( const double time, const DomainType& arg, RangeType& ret )
{
	ret[0] = std::pow(time,3.0)* arg[1]*arg[1];
	ret[1] = std::pow(time,2.0)* arg[0];
}

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
	void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
	{
		const double x			= arg[0];
		const double y			= arg[1];
		const double v			= viscosity_;
		RangeType u;
		evaluateTimeVelocity( time, arg, u );
		//					  ret[0] = std::pow(time,3.0)* arg[1] * arg[1];// * Parameters().getParam( "alpha", 1.0 ) ;
		//					  ret[1] = std::pow(time,2.0)* arg[0];// * Parameters().getParam( "alpha", 1.0 ) ;
		//					  ret *= alpha;
		//laplce
		ret[0] = -2*std::pow(time,3.0)*v;// * Parameters().getParam( "viscosity", 1.0 );
		ret[1] = 0;//-2*std::pow(time,2.0)*v;// * Parameters().getParam( "viscosity", 1.0 );
		//grad p
		ret[0] += time;
		ret[1] += 1;
		//conv
		if ( !Parameters().getParam( "navier_no_convection", false ) ) {
			ret[0] += 2* std::pow(time,5.0)*x*y;
			ret[1] +=  std::pow(time,5.0)*y*y;
		}
		//dt u
		ret[0] += std::pow(time,2.0)*3*y*y;
		ret[1] += 2*time*x;

		//					  ret *=Parameters().getParam( "fscale", 1.0 );
		//					  ret *= 0;


	}

private:
	const double viscosity_;
	const double alpha_;
	static const int dim_ = FunctionSpaceImp::dimDomain;
};

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

	~DirichletData()
	{}
    void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
    {
        evaluateTimeVelocity( time, arg, ret );
    }
	template < class IntersectionType >
	void evaluateTime( const double time, const DomainType& arg, RangeType& ret, const IntersectionType& /*intersection*/ ) const
	{
		evaluateTimeVelocity( time, arg, ret );
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

	Velocity(	const TimeProviderImp& timeprovider,
				const FunctionSpaceImp& space,
				const double parameter_a = M_PI /2.0 ,
				const double parameter_d = M_PI /4.0)
		: BaseType( timeprovider, space ),
		  parameter_a_( parameter_a ),
		  parameter_d_( parameter_d )
	{}

	~Velocity()
	{}

	void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
	{
		dune_static_assert( FunctionSpaceImp::dimDomain == 2  , "Wrong world dim");
		evaluateTimeVelocity( time, arg, ret );
	}

private:
	static const int dim_ = FunctionSpaceImp::dimDomain ;
	const double parameter_a_;
	const double parameter_d_;
};

template < class FunctionSpaceImp, class TimeProviderImp >
class PressureGradient : public Dune::TimeFunction < FunctionSpaceImp , PressureGradient< FunctionSpaceImp,TimeProviderImp >, TimeProviderImp >
{
public:
	typedef PressureGradient< FunctionSpaceImp, TimeProviderImp >
		ThisType;
	typedef Dune::TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
		BaseType;
	typedef typename BaseType::DomainType
		DomainType;
	typedef typename BaseType::RangeType
		RangeType;

	PressureGradient(	const TimeProviderImp& timeprovider,
						const FunctionSpaceImp& space,
						const double parameter_a = M_PI /2.0 ,
						const double parameter_d = M_PI /4.0)
		: BaseType( timeprovider, space ),
		  parameter_a_( parameter_a ),
		  parameter_d_( parameter_d )
	{}

	~PressureGradient()
	{}

	void evaluateTime( const double time, const DomainType& /*arg*/, RangeType& ret ) const
	{
		ret[0] = 1;
		ret[1] = time;
	}

private:
	const double parameter_a_;
	const double parameter_d_;
};

template <	class FunctionSpaceImp, class TimeProviderImp >
class Pressure : public Dune::TimeFunction < FunctionSpaceImp , Pressure < FunctionSpaceImp,TimeProviderImp >, TimeProviderImp >
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

	Pressure( const TimeProviderImp& timeprovider,
			  const FunctionSpaceImp& space,
			  const double parameter_a = M_PI /2.0 ,
			  const double parameter_d = M_PI /4.0)
		: BaseType( timeprovider, space ),
		  parameter_a_( parameter_a ),
		  parameter_d_( parameter_d )
	{}

	~Pressure()
	{}

	void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
	{

		ret = time * arg[0] + arg[1] - ( ( time + 1 ) / 2.0 );
	}

private:
	static const int dim_ = FunctionSpaceImp::dimDomain ;
	const double parameter_a_;
	const double parameter_d_;
};

template < class FunctionSpaceImp, class TimeProviderImp >
class VelocityLaplace : public Dune::TimeFunction < FunctionSpaceImp , VelocityLaplace< FunctionSpaceImp,TimeProviderImp >, TimeProviderImp >
{
public:
	typedef VelocityLaplace< FunctionSpaceImp, TimeProviderImp >
		ThisType;
	typedef Dune::TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
		BaseType;
	typedef typename BaseType::DomainType
		DomainType;
	typedef typename BaseType::RangeType
		RangeType;

	VelocityLaplace(	const TimeProviderImp& timeprovider,
						const FunctionSpaceImp& space,
						const double parameter_a = M_PI /2.0 ,
						const double parameter_d = M_PI /4.0)
		: BaseType( timeprovider, space ),
		  parameter_a_( parameter_a ),
		  parameter_d_( parameter_d )
	{}

	~VelocityLaplace() {}

	void evaluateTime( const double time, const DomainType& /*arg*/, RangeType& ret ) const
	{
		ret[0] = 2*std::pow(time,3.0);
		ret[1] = 0;
	}

private:
	static const int dim_ = FunctionSpaceImp::dimDomain ;
	const double parameter_a_;
	const double parameter_d_;
};

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

	VelocityConvection(	const TimeProviderImp& timeprovider,
						const FunctionSpaceImp& space,
						const double parameter_a = M_PI /2.0 ,
						const double parameter_d = M_PI /4.0)
		: BaseType( timeprovider, space ),
		  parameter_a_( parameter_a ),
		  parameter_d_( parameter_d )
	{}

	~VelocityConvection()
	{}

	void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
	{
		evaluateTimeVelocity( time, arg, ret );
	}

    inline void jacobianTime (	const double time,	const DomainType& /*arg*/,
								typename BaseType::BaseType::JacobianRangeType& ret ) const
	{
		ret[0][0] = 0;
		ret[0][1] = std::pow(time,3.0);
		ret[1][0] = std::pow(time,2.0);
		ret[1][1] = 0;

	}

private:
	static const int dim_ = FunctionSpaceImp::dimDomain ;
	const double parameter_a_;
	const double parameter_d_;
};

}//end ns
}//end ns

#endif //NAVIER_PROBLEMS_TIMEDISC_HH
