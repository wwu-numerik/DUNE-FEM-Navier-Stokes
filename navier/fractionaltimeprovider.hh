#ifndef FRACTIONALTIMEPROVIDER_HH
#define FRACTIONALTIMEPROVIDER_HH

#include <dune/fem/solver/timeprovider.hh>


namespace Dune {
	namespace NavierStokes {


		template< class SchemeParameterType, class CommProvider = DefaultCollectiveCommunicationType >
		class FractionalTimeProvider : public TimeProvider < CommProvider > {
					typedef FractionalTimeProvider< CommProvider >
						ThisType;
					typedef TimeProvider< CommProvider >
						BaseType;

				public:
					using typename BaseType :: CollectiveCommunicationType;
					using BaseType :: time;
					using BaseType :: timeStep;
					using BaseType :: deltaT;

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

				public:
					FractionalTimeProvider (
								   const SchemeParameterType& theta_scheme_parameter,
								   const CommProvider &comm )
						: BaseType( comm ),
						startTime_ ( Parameter :: getValue( "fem.timeprovider.starttime", //this is somewhat duplicated in empty basetype ctor
															   (double)0.0 ) ),
						endTime_ ( Parameter :: getValue( "fem.timeprovider.endtime",
															   (double)1.0 ) ),
						theta_scheme_parameter_( theta_scheme_parameter ),
						current_substep_( -1 )
					{
						dt_ = Parameter :: getValidValue( "fem.timeprovider.dt",
														 (double)0.1,
														 //assure  dt is in (0,endTime_ - startTime_]
														 ValidateInterval<double,false,true>( 0.0, endTime_ - startTime_) );
						init( dt_ );
					}

					//! equivalent of t_{k+1}
					const double subTime( ) const
					{
						assert( current_substep_ > -1 );
						double current = time();
						for ( int i = 0; i < current_substep_; ++i )
							current += theta_scheme_parameter_.step_sizes_[i];
						return current;
					}

					//! equivalent of t_{k}
					const double previousSubTime( ) const
					{
						return subTime() - theta_scheme_parameter_.step_sizes_[current_substep_-1];
					}

					const double time () const
					{
						assert( false );
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

					const double startTime()	const { return startTime_;	}
					const double endTime()		const { return endTime_;	}

					int timeStep () const
					{
						return -1;
					}

				protected:
					void next ( const double timeStep )
					{
						current_substep_ = 0;
						BaseType::next( timeStep );
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
