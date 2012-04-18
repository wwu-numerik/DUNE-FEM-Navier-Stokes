#ifndef NAVIER_PROBLEMS_HH
#define NAVIER_PROBLEMS_HH

#ifndef NAVIER_DATA_NAMESPACE
	#define NAVIER_DATA_NAMESPACE NavierProblems::Trivial
#endif
#ifndef OSEEN_DATA_NAMESPACE
	#define OSEEN_DATA_NAMESPACE NavierProblems::Trivial
#endif

#include "problems/timedisc.hh"
#include "problems/trivial.hh"
#include "problems/taylor.hh"
#include "problems/heat.hh"
#include "problems/testcase2d.hh"
#include "problems/testcase3d.hh"
#include "problems/cockburn.hh"
#include "problems/real_bvp.hh"
#include "problems/2dtube.hh"
#include "problems/runtime.hh"
#include "problems/damped.hh"

namespace NavierProblems {
namespace {
	//! a dummy interface class that throws errors if we've forgotten to implement data funcs in given namespace
	template < class FunctionSpaceImp, class TimeProviderImp >
	struct Interface {
		typedef OSEEN_DATA_NAMESPACE::PressureGradient<FunctionSpaceImp, TimeProviderImp>
			PressureGradient;
		typedef OSEEN_DATA_NAMESPACE::VelocityConvection<FunctionSpaceImp, TimeProviderImp>
			VelocityConvection;
		typedef OSEEN_DATA_NAMESPACE::VelocityLaplace<FunctionSpaceImp, TimeProviderImp>
			VelocityLaplace;
		typedef OSEEN_DATA_NAMESPACE::Velocity<FunctionSpaceImp, TimeProviderImp>
			Velocity;
		typedef OSEEN_DATA_NAMESPACE::Force<FunctionSpaceImp, TimeProviderImp>
			Force;
		typedef OSEEN_DATA_NAMESPACE::DirichletData<FunctionSpaceImp, TimeProviderImp>
			DirichletData;
	};
}
}

#endif // NAVIER_PROBLEMS_HH
