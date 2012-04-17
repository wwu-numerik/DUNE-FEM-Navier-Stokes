#ifndef NAVIER_PROBLEMS_RUNTIME_HH
#define NAVIER_PROBLEMS_RUNTIME_HH

#include <dune/stuff/functions.hh>
#include <dune/stuff/timefunction.hh>
#include <dune/stuff/runtimefunction.hh>
#include <dune/stuff/parametercontainer.hh>
#include "common.hh"

namespace NavierProblems {
namespace Runtime {

static const std::string identifier = "Runtime";
static const bool hasExactSolution	= true;
ALLGOOD_SETUPCHECK;

NV_RUNTIME_FUNC(Force);
NV_RUNTIME_FUNC(VelocityConvection);
NV_RUNTIME_FUNC(Velocity);
NV_RUNTIME_FUNC(Pressure);
NV_RUNTIME_FUNC(VelocityLaplace);
NV_RUNTIME_FUNC(PressureGradient);
NV_RUNTIME_FUNC(Beta);

template < class FunctionSpaceImp, class TimeProviderImp >
class DirichletData : public Dune::IntersectionTimeFunction < FunctionSpaceImp ,
                                                              DirichletData< FunctionSpaceImp,TimeProviderImp >,
                                                              TimeProviderImp >
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
            typedef Velocity< FunctionSpaceImp, TimeProviderImp >
                VelocityType;
        /**
        *  \brief  constructor
        *
        *  doing nothing besides Base init
        **/
        DirichletData( const TimeProviderImp& timeprovider,
                   const FunctionSpaceImp& space,
                   const double  /*viscosity*/ = 0.0,
                   const double /*alpha*/ = 0.0 )
            : BaseType ( timeprovider, space ),
              velocity_( timeprovider, space )
        {}

        ~DirichletData()
        {}

        template < class IntersectionType >
        void evaluateTime( const double time, const DomainType& arg,
                           RangeType& ret, const IntersectionType& /*intersection */) const
        {
            velocity_.evaluateTime( time, arg, ret );
        }

        void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
        {
            velocity_.evaluateTime( time, arg, ret );
        }
    private:
          const VelocityType velocity_;
};

}//end ns

}//end ns

#endif //NAVIER_PROBLEMS_RUNTIME_HH