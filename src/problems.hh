#ifndef NAVIER_PROBLEMS_HH
#define NAVIER_PROBLEMS_HH

//#include "problems/timedisc.hh"
#include "problems/trivial.hh"
//#include "problems/taylor.hh"
//#include "problems/testcase2d.hh"
//#include "problems/testcase3d.hh"
//#include "problems/cockburn.hh"

//#include <dune/fem/space/

#if defined(NAVIER_DATA_NAMESPACE) || defined(OSEEN_DATA_NAMESPACE)
	#if defined(NAVIER_DATA_NAMESPACE)
		#define CHECKING_NS NAVIER_DATA_NAMESPACE
	#else
		#define CHECKING_NS OSEEN_DATA_NAMESPACE
	#endif
namespace NavierProblems {
namespace {
	template < class FunctionSpaceImp, class TimeProviderImp >
	struct Interface {
		typedef CHECKING_NS::PressureGradient<FunctionSpaceImp, TimeProviderImp>
			PressureGradient;
		typedef CHECKING_NS::VelocityConvection<FunctionSpaceImp, TimeProviderImp>
			VelocityConvection;
		typedef CHECKING_NS::VelocityLaplace<FunctionSpaceImp, TimeProviderImp>
			VelocityLaplace;
		typedef CHECKING_NS::Velocity<FunctionSpaceImp, TimeProviderImp>
			Velocity;
		typedef CHECKING_NS::Force<FunctionSpaceImp, TimeProviderImp>
			Force;
		typedef CHECKING_NS::DirichletData<FunctionSpaceImp>
			DirichletData;
	};
}
}
#undef CHECKING_NS
#endif //defined(NAVIER_DATA_NAMESPACE) || defined(OSEEN_DATA_NAMESPACE)

#endif // NAVIER_PROBLEMS_HH
