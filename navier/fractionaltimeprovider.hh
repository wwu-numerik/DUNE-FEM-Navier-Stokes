#ifndef FRACTIONALTIMEPROVIDER_HH
#define FRACTIONALTIMEPROVIDER_HH

#include <dune/common/deprecated.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/misc/femtimer.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/math.hh>
#include <boost/format.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace Dune {
	namespace NavierStokes {


		template< class SchemeParameterType, class CommProvider = DefaultCollectiveCommunicationType >
		class FractionalTimeProvider : public TimeProvider < CommProvider > {
					typedef FractionalTimeProvider< SchemeParameterType, CommProvider >
						ThisType;
					typedef TimeProvider< CommProvider >
						BaseType;

				public:
					using typename BaseType :: CollectiveCommunicationType;
					using BaseType :: time;
					using BaseType :: timeStep;
					using BaseType :: deltaT;

					class StepZeroGuard {
						ThisType& timeprovider_;
						public:
							StepZeroGuard( const double initial_dt_, ThisType& timeprovider )
								:timeprovider_(timeprovider)
							{
								timeprovider_.init( initial_dt_ );
							}
							~StepZeroGuard()
							{
								timeprovider_.nextFractional();
								timeprovider_.nextFractional();
							}
					};

				protected:
					using BaseType :: comm_;
					using BaseType :: cfl_;
					using BaseType :: dt_;
					using BaseType :: dtEstimate_;
					using BaseType :: dtUpperBound_;
					using BaseType :: valid_;
					using BaseType :: timeStep_;
					static const int substep_count_  = SchemeParameterType::numberOfSteps_;
					const double startTime_;
					const double endTime_;
					const SchemeParameterType& theta_scheme_parameter_;
					int current_substep_;
					ExecutionTimer step_timer_;
					Stuff::MovingAverage avg_time_per_step_;
					long total_stepcount_estimate_;

				public:
					FractionalTimeProvider (
								   const SchemeParameterType& theta_scheme_parameter,
								   const CommProvider &comm )
						: BaseType( comm ),
						startTime_ ( Parameter :: getValue( "fem.timeprovider.starttime", //this is somewhat duplicated in empty basetype ctor
															   (double)0.0 ) ),
						endTime_ ( Parameter :: getValidValue( "fem.timeprovider.endtime",
															   (double)1.0 ,
																ValidateGreater<double>(startTime_) ) ),
						theta_scheme_parameter_( theta_scheme_parameter ),
						current_substep_( -1 ),
						total_stepcount_estimate_( -1 )
					{
						dt_ = Parameter :: getValidValue( "fem.timeprovider.dt",
														 (double)0.1,
														 //assure  dt is in (0,endTime_ - startTime_]
														 ValidateInterval<double,false,true>( 0.0, endTime_ - startTime_) );
						init( dt_ );
						total_stepcount_estimate_ = long( std::ceil( ( endTime_ - startTime_ ) / dt_ ) );
						step_timer_.start();
					}

					//! equivalent of t_{k+1}
					double subTime( ) const
					{
						double current = BaseType::time();
						for ( int i = 0; i < current_substep_; ++i )
							current += theta_scheme_parameter_.step_sizes_[i];
						return current;
					}

					//! equivalent of t_{k}
					double previousSubTime( ) const
					{
						const double t = subTime() - theta_scheme_parameter_.step_sizes_[current_substep_-1];
						return Stuff::clamp( t, double(0.0), t);
					}
					//! equivalent of t_{k+2}??
					double nextSubTime( ) const
					{
						assert( false ); //currently produces wring results on the interval bounds
						const double t = subTime() + theta_scheme_parameter_.step_sizes_[current_substep_];
						return Stuff::clamp( t, double(0.0), t);
					}

					double time () const
					{
						return subTime();
					}

					void nextFractional()
					{
						if ( current_substep_ == -1 ) {
							current_substep_ = 0;
						}
						else if ( current_substep_ == substep_count_  -1 ) {
							next( deltaT() );
						}
						else
							++current_substep_;
					}

					//! return t_{n+1} - t_{n}
					double deltaT() const {
						return dt_;
					}

					//! return t_{k+1} - t_{k}
					double sub_deltaT() const {
						return theta_scheme_parameter_.step_sizes_[current_substep_];
					}

					double startTime()	const { return startTime_;	}
					double endTime()	const { return endTime_;	}

					int timeStep () const
					{
						const int ret = timeStep_ * substep_count_ + current_substep_ + 1;
						assert( ret >= 0 );
						return ret;
					}

					template < class Stream >
					void printRemainderEstimate( Stream& stream )
					{
						long remaining_steps = total_stepcount_estimate_ - ( timeStep_ );
						double remaining_seconds = remaining_steps * double(avg_time_per_step_);
						boost::posix_time::time_duration diff(0,0,remaining_seconds,0);
						boost::posix_time::ptime target = boost::posix_time::second_clock::local_time();
						target += diff;
						stream << boost::format("\n---\nTime remaining: %s -- %s(%f %%)\n---\n")
									% boost::posix_time::to_simple_string(diff)
									% boost::posix_time::to_simple_string(target)
									% (100 * ( remaining_steps/double(total_stepcount_estimate_) ) );
					}

					StepZeroGuard stepZeroGuard( const double dt )
					{
						return StepZeroGuard( dt, *this );
					}

				protected:
					void next ( const double timeStep )
					{
						assert( timeStep > 0 );
						current_substep_ = 0;
						// timer
						step_timer_.end();
						avg_time_per_step_ += std::abs( step_timer_.read() );
						BaseType::next( timeStep );
						step_timer_.start();
					}

					//! hidden since outside calling is nonsensical
					void next (  )
					{
						assert( false );//make whatever triggers this call nextFractional instead
					}

		};
	}//end namespace NavierStokes
} //end namespace Dune

#endif // FRACTIONALTIMEPROVIDER_HH
