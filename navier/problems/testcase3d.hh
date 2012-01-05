#ifndef NAVIER_PROBLEMS_THREEDEE_HH
#define NAVIER_PROBLEMS_THREEDEE_HH

#include <dune/stuff/functions.hh>
#include <dune/stuff/timefunction.hh>
#include <dune/stuff/parametercontainer.hh>
#include "common.hh"

namespace NavierProblems {
namespace ThreeDee {
ALLGOOD_SETUPCHECK;
static const std::string identifier = "Testcase3D";
static const bool hasExactSolution	= true;

}//end ns
}//end ns

#endif // NAVIER_PROBLEMS_THREEDEE_HH
