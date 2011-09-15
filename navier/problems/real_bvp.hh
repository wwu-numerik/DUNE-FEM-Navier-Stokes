#ifndef REAL_BVP_HH
#define REAL_BVP_HH

#include <dune/stuff/functions.hh>
#include <dune/stuff/timefunction.hh>
#include <dune/stuff/parametercontainer.hh>
#include <dune/stokes/boundarydata.hh>

namespace NavierProblems {

namespace Null {
static const std::string identifier = "NULL";
static const bool hasExactSolution	= true;
    NULLFUNCTION_TP(PressureGradient)
    NULLFUNCTION_TP(Pressure)
    NULLFUNCTION_TP(VelocityConvection)
    NULLFUNCTION_TP(VelocityLaplace)
    NULLFUNCTION_TP(Velocity)
    NULLFUNCTION_TP(Force)
    NULLFUNCTION_TP_BOUNDARY(DirichletData)
}

namespace BVP {

static const std::string identifier = "BVP";
static const bool hasExactSolution	= false;
//currently not possible, see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=45114
//    template < class FunctionSpaceImp, class TimeProviderImp >
//    using DirichletData = Stuff::InstationaryBoundaryFluxFunction<FunctionSpaceImp, TimeProviderImp>;
//    template < class FunctionSpaceImp, class TimeProviderImp >
//    class DirichletData : public Stuff::InstationaryBoundaryFluxFunction<FunctionSpaceImp, TimeProviderImp>
//    {
//	    typedef Stuff::InstationaryBoundaryFluxFunction<FunctionSpaceImp, TimeProviderImp>
//		BaseType;
//	public:
//	    DirichletData( const TimeProviderImp& timeprovider,
//			   const FunctionSpaceImp& space,
//			   const double /*constant*/ = 0.0 )
//	    : BaseType ( timeprovider, space )
//	    {}
//    };

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
		    *
		    *  doing nothing besides Base init
		    **/
		    DirichletData( const TimeProviderImp& timeprovider,
				       const FunctionSpaceImp& space,
				       const double  /*viscosity*/ = 0.0,
				       const double /*alpha*/ = 0.0 )
			    : BaseType ( timeprovider, space ),
			    z_max( Parameters().getParam( "z_max", 3.0 ) )
		    {}

		    /**
		    *  \brief  destructor
		    *
		    *  doing nothing
		    **/
		    ~DirichletData()
		    {}
		    void evaluateTime( const double /*time*/, const DomainType& /*arg*/, RangeType& ret ) const
		    {
			ret = RangeType( 0 );
		    }

		    template < class IntersectionType >
		    void evaluateTime( const double time, const DomainType& arg, RangeType& ret, const IntersectionType& /*intersection */) const
		    {
			dune_static_assert( dim_ == 3 ,"__CLASS__ evaluate not implemented for world dimension");
			DomainType normal(0);
			normal[2] =  1 ;
			const double gd_factor = time * Parameters().getParam( "gd_factor", 1.0 );
			if ( arg[2] > 0.0 && arg[2] < z_max )
			{
			    normal *= 0.0;
			}
			ret = normal;
			ret *= gd_factor;
		    }

	    private:
		      static const int dim_ = FunctionSpaceImp::dimDomain ;
		      const double z_max;
    };

    NULLFUNCTION_TP(PressureGradient)
    NULLFUNCTION_TP(Pressure)
    NULLFUNCTION_TP(VelocityConvection)
    NULLFUNCTION_TP(VelocityLaplace)
    NULLFUNCTION_TP(Velocity)
    NULLFUNCTION_TP(Force)
}//end ns
}//end ns


#endif // REAL_BVP_HH
