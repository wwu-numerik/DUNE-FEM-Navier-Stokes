#ifndef WEIGHED_FUNCTION_HH
#define WEIGHED_FUNCTION_HH


#include <dune/stuff/functions.hh>
#include <dune/stuff/timefunction.hh>

namespace Dune {
namespace NavierStokes {

/** simple wrapper that allows us to create functions
 *
 */
template < class FunctionSpaceImp, class TimeProviderImp, class FunctionImp >
class WeighedIntersectionFunction : public Dune::IntersectionTimeFunction < FunctionSpaceImp ,
        WeighedIntersectionFunction< FunctionSpaceImp, TimeProviderImp, FunctionImp >, TimeProviderImp >
{
        public:
            typedef WeighedIntersectionFunction< FunctionSpaceImp, TimeProviderImp, FunctionImp >
                ThisType;
            typedef Dune::IntersectionTimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
                BaseType;
            typedef typename BaseType::DomainType
                DomainType;
            typedef typename BaseType::RangeType
                RangeType;
            typedef FunctionImp
                FunctionType;
        /**
        *  \brief  constructor
        *
        *  doing nothing besides Base init
        **/
        WeighedIntersectionFunction( const TimeProviderImp& timeprovider,
                   const FunctionSpaceImp& space,
                   const double weight_a = 1.0,
                   const double weight_b = 0.0 )
            : BaseType ( timeprovider, space ),
              weight_a_(weight_a),
              weight_b_(weight_b),
              function_(timeprovider, space)
        {}

        /**
        *  \brief  destructor
        *
        *  doing nothing
        **/
        ~WeighedIntersectionFunction()
        {}
        void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
        {
            RangeType b;
            function_.evaluateTime(time, arg, ret);
            function_.evaluateTime(time - BaseType::timeProvider_.deltaT(), arg, b);
            ret *= weight_a_;
            b *= weight_b_;
            ret += b;
        }

        template < class IntersectionType >
        void evaluateTime( const double time, const DomainType& arg, RangeType& ret, const IntersectionType& intersection ) const
        {
            RangeType b;
            function_.evaluateTime(time, arg, ret, intersection);
            function_.evaluateTime(time - BaseType::timeProvider_.deltaT(), arg, b, intersection);
            ret *= weight_a_;
            b *= weight_b_;
            ret += b;
        }

    private:
        const double weight_a_;
        const double weight_b_;
        FunctionType function_;
};

} //namespace Dune
} //namespace NavierStokes
#endif // WEIGHED_FUNCTION_HH
