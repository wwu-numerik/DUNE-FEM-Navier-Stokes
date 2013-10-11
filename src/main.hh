#ifndef MAIN_H
#define MAIN_H

/**
 *  \file   dune_stokes.cc
 *
 *  \brief  brief
 **/
#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#endif
#include <dune/grid/utility/gridtype.hh>
#include <dune/fem/nvs/global_defines.hh>

#include <cstdio>
#if defined(USE_PARDG_ODE_SOLVER)
#warning("USE_PARDG_ODE_SOLVER enabled, might conflict with custom solvers")
#endif

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include <dune/fem/misc/mpimanager.hh> // An initializer of MPI
#include <dune/common/exceptions.hh>   // We use exceptions
#include <dune/grid/common/capabilities.hh>

typedef Dune::GridSelector::GridType GridType;

#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/pass/pass.hh>
#include <dune/fem/function/adaptivefunction.hh> // for AdaptiveDiscreteFunction
#include <dune/fem/misc/gridwidth.hh>

#include <dune/fem/oseen/functionspacewrapper.hh>
#include <dune/fem/oseen/modelinterface.hh>
#include <dune/fem/oseen/ldg_method.hh>
#include <dune/fem/oseen/boundarydata.hh>

#include <dune/stuff/fem/femeoc.hh>
#include <dune/stuff/common/misc.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/common/profiler.hh>
//#include <dune/stuff/fem/timeseries.hh>
#include <dune/stuff/common/signals.hh>
#include <dune/fem/oseen/runinfo.hh>

#include <dune/fem/nvs/thetascheme_runner.hh>
#include <dune/fem/nvs/fractionaldatawriter.hh>

#if ENABLE_MPI
typedef Dune::CollectiveCommunication<MPI_Comm> CollectiveCommunication;
#else
typedef Dune::CollectiveCommunication<double> CollectiveCommunication;
#endif

//! the strings used for column headers in tex output
typedef std::vector<std::string> ColumnHeaders;

void eocCheck(const DSC::RunInfoVector& runInfos) {
  bool ups = false;
  DSC::RunInfoVector::const_iterator it = runInfos.begin();
  DSC::RunInfo last = *it;
  ++it;
  for (; it != runInfos.end(); ++it) {
    ups = (last.L2Errors[0] < it->L2Errors[0] || last.L2Errors[1] < it->L2Errors[1]);
    last = *it;
  }
  if (ups) {
    DSC_LOG_ERROR << "----------------------------------------------------------\n"
                  << "-                                                        -\n"
                  << "-                  negative EOC                          -\n"
                  << "-                                                        -\n"
                  << "----------------------------------------------------------\n" << std::endl;
  }
}

CollectiveCommunication init(int argc, char** argv) {
  Dune::MPIManager::initialize(argc, argv);
  //    assert( Dune::Capabilities::isParallel< GridType >::v );

  /* ********************************************************************** *
   * initialize all the stuff we need                                       *
   * ********************************************************************** */
  if (argc < 2) {
    std::cerr << "\nUsage: " << argv[0] << " parameterfile \n"
              << "\n\t --- OR --- \n";
    std::cerr << "\nUsage: " << argv[0] << " paramfile:"
              << "file"
              << " more-opts:val ..." << std::endl;
    std::cerr << "\nUsage: " << argv[0] << " -d paramfile "
              << "\n\t(for displaying solutions in grape) " << std::endl;
    std::cerr << std::endl;
    return 2;
  }

  DSC_CONFIG.readCommandLine(argc, argv);

  // LOG_NONE = 1, LOG_ERR = 2, LOG_INFO = 4,LOG_DEBUG = 8,LOG_CONSOLE = 16,LOG_FILE = 32
  //--> LOG_ERR | LOG_INFO | LOG_DEBUG | LOG_CONSOLE | LOG_FILE = 62
  const bool useLogger = false;
  DSC_LOG.create(DSC_CONFIG_GETB("loglevel", 62, useLogger),
                 DSC_CONFIG_GETB("logfile", std::string("dune_stokes"), useLogger),
                 DSC_CONFIG_GETB("fem.io.datadir", std::string("data"), useLogger),
                 DSC_CONFIG_GETB("fem.io.logdir", std::string(), useLogger));

  return CollectiveCommunication(); //( Dune::MPIManager::helper().getCommunicator() );
}

#endif // MAIN_H

/** Copyright (c) 2012, Rene Milk
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies,
 * either expressed or implied, of the FreeBSD Project.
**/
