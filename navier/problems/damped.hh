#ifndef DAMPED_HH
#define DAMPED_HH

#include <dune/stuff/functions.hh>
#include <dune/stuff/timefunction.hh>
#include <dune/stuff/parametercontainer.hh>
#include <dune/oseen/boundarydata.hh>
#include "common.hh"

namespace NavierProblems {

namespace DampedParallel {

static const std::string identifier = "TwoDeeDampedParallel";
static const bool hasExactSolution	= true;
ALLGOOD_SETUPCHECK;

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
    inline void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
    {
        const double x				= arg[0];
        const double y				= arg[1];
        RangeType u;
        VelocityEvaluate( lambda_, 0, arg, u);
        //convection
        //					  assert( false ); //M_2_PI == 2 / PI
        ret[0] = 0;
        ret[1] = 0;
        //laplace
        ret[0] -= 0;
        ret[1] -= 0;
        //pressure grad
        ret[0] += (BaseType::timeProvider_.endTime() - time);
        ret[1] += 0;
        //dtu
        ret[0] -= 1;
        ret[1] += 0;
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
void VelocityEvaluate( const double endtime, const double time, const DomainType& arg, RangeType& ret)
{
    const double x				= arg[0];
    const double y				= arg[1];
    ret[0] = (endtime - time);
    ret[1] = 0;
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
                        const double /*parameter_a*/ = M_PI /2.0 ,
                        const double /*parameter_d*/ = M_PI /4.0)
        : BaseType( timeprovider, space )
    {}

    /**
   *  \brief  destructor
   *
   *  doing nothing
   **/
    ~VelocityConvection()
    {}

    template < class IntersectionType >
    void evaluateTime( const double /*time*/, const DomainType& /*arg*/, RangeType& ret, const IntersectionType& /*intersection */) const
    {
        dune_static_assert( FunctionSpaceImp::dimDomain == 2, "__CLASS__ evaluate not implemented for world dimension");
        ret = RangeType(0);
    }

    void evaluateTime( const double /*time*/, const DomainType& /*arg*/, RangeType& ret ) const
    {ret = RangeType(0);}

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
        VelocityEvaluate( BaseType::timeProvider_.endTime(), time, arg, ret);
    }
    void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
    {
        VelocityEvaluate( BaseType::timeProvider_.endTime(), time, arg, ret);
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
                const double /*parameter_a*/ = M_PI /2.0 ,
                const double /*parameter_d*/ = M_PI /4.0)
        : BaseType( timeprovider, space )
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
        VelocityEvaluate( BaseType::timeProvider_.endTime(), time, arg, ret);
    }

    /**
   * \brief  evaluates the dirichlet data
   * \param  arg
   *         point to evaluate at
   * \param  ret
   *         value of dirichlet boundary data at given point
   **/
    //					inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}
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
              const double /*parameter_a*/ = M_PI /2.0 ,
              const double /*parameter_d*/ = M_PI /4.0)
        : BaseType( timeprovider, space )
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
        ret = (BaseType::timeProvider_.endTime() - time) * arg[0];
    }

    void setShift( const double /*shift*/ )
    {

    }

    /**
   * \brief  evaluates the dirichlet data
   * \param  arg
   *         point to evaluate at
   * \param  ret
   *         value of dirichlet boundary data at given point
   **/
    //					inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}
};

NULLFUNCTION_TP(PressureGradient)
NULLFUNCTION_TP(VelocityLaplace)
}//end namespace DampedParallel
}//end namespace NavierProblems

#endif // DAMPED_HH
