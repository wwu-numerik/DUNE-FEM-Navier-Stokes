
#include "main.hh"
#include <dune/fem/nvs/problems.hh>
#include <dune/stuff/common/loop_timer.hh>

/** \brief one single application of the discretisation and solver

	\param  mpicomm
			mostly useless atm, but mandatory
	\param  refine_level_factor
			integer to be multiplied by Dune::DGFGridInfo< GridType >::refineStepsForHalf()
			to get the used refine level for the constructed grid
	\param  stabil_coeff
			the set of coefficients to be used in the run. Default is used in all run types but StabRun().

**/
DSC::RunInfoTimeMap singleRun(	CollectiveCommunication& mpicomm,
							const int refine_level_factor,
							const int scheme_type );
//! output alert for neg. EOC
//void eocCheck( const RunInfoVector& runInfos );

bool setSchemeTypeFromString();

/**
 *  \brief  main function
 *
 *  ParameterContainer and Logger setup, select run type from parameter file and go
 *
 *  \param  argc
 *          number of arguments from command line
 *  \param  argv
 *          array of arguments from command line
 **/
int main( int argc, char** argv )
{
    CollectiveCommunication mpicomm( init(argc,argv) );//( Dune::MPIManager::helper().getCommunicator() );

	if ( setSchemeTypeFromString() )
        DSC_LOG_INFO << "overrode scheme id from string" << std::endl;

	int err = 0;
    const unsigned int minref = DSC_CONFIG_GETV( "minref", 0, DSC::ValidateNotLess<int>(0) );
    DSC::RunInfoTimeMapMap rf;
	const int runtype = DSC_CONFIG_GET( "runtype", 5 );
	switch( runtype ) {
		case 8: {
            DSC_LOG_INFO << "Reynolds runs\n";
            const int dt_steps = DSC_CONFIG_GETV( "dt_steps", 3, DSC::ValidateNotLess<int>(2) );
            DSC_PROFILER.reset( dt_steps - 1 );
			int current_step = 0;
            DSC::LoopTimer<int> loop_timer( current_step, dt_steps, DSC_LOG_INFO );
            for ( double viscosity = DSC_CONFIG_GETV( "viscosity", 0.1, DSC::ValidateNotLess<double>(0.0) );
				  dt_steps > current_step;
				  ++loop_timer )
			{
                rf[current_step] = singleRun( mpicomm, minref, DSC_CONFIG_GETB( "scheme_type", 1, true ) );
				assert( rf.size() );
				rf[current_step].begin()->second.refine_level = minref;//just in case the key changes from ref to sth else
                DSC_PROFILER.nextRun();
				viscosity /= 10.0f;
                DSC_CONFIG.set( "viscosity", viscosity );
			}
			break;
		}
		case 6: {
            DSC_LOG_INFO << "Time refine runs\n";
            const int dt_steps = DSC_CONFIG_GETV( "dt_steps", 3, DSC::ValidateNotLess<int>(2) );
            DSC_PROFILER.reset( dt_steps - 1 );
			int current_step = 0;
            DSC::LoopTimer<int,DSC::QuadraticWeights> loop_timer( current_step, dt_steps, DSC_LOG_INFO );
            for ( double dt = DSC_CONFIG_GETV( "fem.timeprovider.dt", 0.1, DSC::ValidateNotLess<double>(0.0) );
				  dt_steps > current_step;
				  ++loop_timer )
			{
                rf[current_step] = singleRun( mpicomm, minref, DSC_CONFIG_GETB( "scheme_type", 1, true ) );
				assert( rf.size() );
				rf[current_step].begin()->second.refine_level = minref;//just in case the key changes from ref to sth else
                DSC_PROFILER.nextRun();
				dt /= 2.0f;
                DSC_CONFIG.set( "fem.timeprovider.dt", dt );
			}
			break;
		}
		case 7: {
            DSC_LOG_INFO << "Scheme runs\n";
            DSC_PROFILER.reset( 4 );
			int current_scheme = 2;
            DSC::LoopTimer<int> loop_timer( current_scheme, 5, DSC_LOG_INFO );
			for ( ;
				  current_scheme < 6;
				  ++loop_timer )
			{
				rf[current_scheme] = singleRun( mpicomm, minref, current_scheme );
				assert( rf.size() );
				rf[current_scheme].begin()->second.refine_level = minref;//just in case the key changes from ref to sth else
                DSC_PROFILER.nextRun();
			}
			break;
		}
		case 5:
            DSC_CONFIG.set( "maxref", minref );//only one run with ref=minref
		case 0:
		default: {
			// ensures maxref>=minref
            const unsigned int maxref = DSC::clamp( DSC_CONFIG_GET( "maxref", (unsigned int)(0) ),
                                                    minref,
                                                    DSC_CONFIG_GET( "maxref", (unsigned int)(0) ) );
            DSC_PROFILER.reset( maxref - minref + 1 );
            DSC_LOG_INFO << "Grid refine runs\n";
			unsigned int ref = minref;
            DSC::LoopTimer<unsigned int,DSC::LinearWeights> loop_timer( ref, maxref - minref + 1, DSC_LOG_INFO );
			for ( ;
				  ref <= maxref;
				  ++loop_timer )
			{
                rf[ref] = singleRun( mpicomm, ref, DSC_CONFIG_GETB( "scheme_type", 1, true ) );
				rf[ref].begin()->second.refine_level = ref;//just in case the key changes from ref to sth else
                DSC_PROFILER.nextRun();
			}
			break;
		}
	}
//    DSC_PROFILER.outputMap( mpicomm, rf );

//	if ( NAVIER_DATA_NAMESPACE::hasExactSolution && DSC_CONFIG_GET( "calculate_errors", true ) )
//	{
//        DSFe::TimeSeriesOutput out( rf );
//	    out.writeTex( DSC_CONFIG_GET("fem.io.datadir", std::string(".") ) + std::string("/timeseries") );
//	}

	DSC_LOG_DEBUG << "\nRun from: " << commit_string << std::endl;
	return err;
}

DSC::RunInfoTimeMap singleRun(  CollectiveCommunication& mpicomm,
					const int refine_level_factor, const int scheme_type )
{
    DSC::Profiler::ScopedTiming pf_t( "SingleRun" );
    DSC::LogStream& infoStream = DSC_LOG_INFO;
    DSC::LogStream& debugStream = DSC_LOG_DEBUG;

	infoStream << "\n- initialising grid" << std::endl;
    Dune::GridPtr< GridType > gridPtr( DSC_CONFIG_GET( "dgf_file", "grid_2d.dgf") );
	const int refine_level = ( refine_level_factor  ) * Dune::DGFGridInfo< GridType >::refineStepsForHalf();
	gridPtr->globalRefine( refine_level );

	const int polOrder = POLORDER;
	debugStream << "  - polOrder: " << polOrder << std::endl;
//    const double grid_width = Dune::GridWidth::calcGridWidth( gridPart );
//    infoStream << (boost::format("  - max grid width: %f\n") % grid_width) << std::endl;

	try {
        return ThetaschemeRunner<GridType,CollectiveCommunication>(*gridPtr,mpicomm).run( scheme_type );
	}
	catch (Dune::Exception &e){
		std::cerr << "Dune reported error: " << e.what() << std::endl;
	}
	catch ( std::bad_alloc& b ) {
		std::cerr << "Memory allocation failed: " << b.what() ;
        DSC_LOG_INFO.resume();
        DSC::meminfo( DSC_LOG_INFO );
	}
	catch ( assert_exception& a ) {
		std::cerr << "Exception thrown at:\n" << a.what() << std::endl ;
	}
	catch (...){
		std::cerr << "Unknown exception thrown!" << std::endl;
	}
    return DSC::RunInfoTimeMap();
}

bool setSchemeTypeFromString()
{
	bool changed = false;
	if ( Dune::Parameter::exists("scheme_type_string") )
	{
		const std::vector<std::string>& scheme_names = Dune::NavierStokes::ThetaSchemeDescription<0>::scheme_names;
        DSC::ValidateInList<std::string> validator( scheme_names );
        const std::string scheme_string = DSC_CONFIG_GETV("scheme_type_string", scheme_names.back(), validator );
        const int scheme_id = DSC::getIdx( scheme_names, scheme_string );
        DSC_CONFIG.set( "scheme_type", scheme_id );
		changed = true;
	}
	return changed;
}

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

