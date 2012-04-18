#ifndef TWODTUBE_HH
#define TWODTUBE_HH

#include <dune/stuff/functions.hh>
#include <dune/stuff/timefunction.hh>
#include <dune/stuff/parametercontainer.hh>
#include <dune/oseen/boundarydata.hh>
#include "common.hh"

namespace NavierProblems {

namespace TwoDeeTube {

static const std::string identifier = "TwoDeeTube";
static const bool hasExactSolution	= false;
ALLGOOD_SETUPCHECK;

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
		void evaluateTime( const double time, const DomainType& arg, RangeType& ret, const IntersectionType& intersection ) const
		{
		    dune_static_assert( dim_ == 2 ,"__CLASS__ evaluate not implemented for world dimension");
		    const int id = intersection.boundaryId();
		    const double y = arg[1];
		    const double f = -4 * ( y - 0.5) * ( y + 0.5) ;
		    ret = RangeType( 0.0 );
		    if ( id == 2 ) { // bottom
			ret[ 0 ] = 0.0;
			ret[ 1 ] = 0.0;
		    }
		    else if ( id == 3 ) { // right
			ret[ 0 ] = std::abs( std::sin( 2.0 * M_PI * time ) ) * f;
			ret[ 1 ] = 0.0;
		    }
		    else if ( id == 4 ) { // top
			ret[ 0 ] = 0.0;
			ret[ 1 ] = 0.0;
		    }
		    else if ( id == 5 ) { // left
			ret[ 0 ] = std::abs( std::sin( 2.0 * M_PI * time ) ) * f;
			ret[ 1 ] = 0.0;
		    }
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
}//end namespace TwoDeeTube


}//end namespace NavierProblems

#endif // TWODTUBE_HH
