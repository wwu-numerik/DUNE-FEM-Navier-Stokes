#ifndef FRACTIONALTIMEPROVIDER_HH
#define FRACTIONALTIMEPROVIDER_HH

#include <dune/fem/solver/timeprovider.hh>


namespace Dune {
	namespace NavierStokes {

		typedef unsigned int
				StepType;
		const StepType StokesStepA		= 0;
		const StepType NonlinearStepID	= 1;
		const StepType StokesStepB		= 2;

		template< class CommProvider = DefaultCollectiveCommunicationType >
		class FractionalTimeProvider : public TimeProvider < CommProvider > {
					typedef FractionalTimeProvider< CommProvider >
						ThisType;
					typedef TimeProvider< CommProvider >
						BaseType;

				public:
					using BaseType ::  CollectiveCommunicationType;
					using BaseType :: time;
					using BaseType :: timeStep;
					using BaseType :: deltaT;
					typedef BaseType
						SubStepTimeproviderType;

				protected:
					using BaseType :: comm_;
					using BaseType :: cfl_;
					using BaseType :: dt_;
					using BaseType :: dtEstimate_;
					using BaseType :: dtUpperBound_;
					using BaseType :: valid_;
					using BaseType :: timeStep_;
					const double startTime_;
					const double endTime_;
					const double theta_;
					const double theta_alpha_;
					const double theta_beta_;
					StepType currentStepType_;
					SubStepTimeproviderType* subProvider_;

				public:
					FractionalTimeProvider (
								   const double theta,
								   const double theta_alpha,
								   const double theta_beta,
								   const CommProvider &comm )
						: BaseType( comm ),
						startTime_ ( Parameter :: getValue( "fem.timeprovider.starttime", //this is somewhat duplicated in empty basetype ctor
															   (double)0.0 ) ),
						endTime_ ( Parameter :: getValue( "fem.timeprovider.endtime",
															   (double)1.0 ) ),
						theta_( theta ),
						theta_alpha_( theta_alpha ),
						theta_beta_( theta_beta ),
						currentStepType_( StokesStepA ),
						subProvider_( 0 )
					{
						dt_ = Parameter :: getValidValue( "fem.timeprovider.dt",
														 (double)0.1,
														 ValidateGreater<double>(0.0) );
						init( dt_ );
					}

					const double subTime( ) const
					{
						switch ( currentStepType_ ) {
							case StokesStepA:
							default: return BaseType::time();
							case NonlinearStepID: return BaseType::time()+ deltaT() * theta_ ;
							case StokesStepB: return BaseType::time()+ deltaT() * (1 - theta_);
						}
					}

					const double time () const
					{
						return subTime();
					}

					void nextFractional()
					{
						if ( currentStepType_ == StokesStepB ) {
							next( deltaT() );
						}
						else
							++currentStepType_;
					}

					const double alpha ()		const { return theta_alpha_;}
					const double beta ()		const { return theta_beta_; }
					const double startTime()	const { return startTime_;	}
					const double endTime()		const { return endTime_;	}

					SubStepTimeproviderType& subStepTimeprovider( const double =0.0 )
					{
						delete subProvider_;
						subProvider_ = new SubStepTimeproviderType(subTime(),BaseType::comm_);
						const double stepsize = ( 1 - 2 * theta_ ) * deltaT();
						subProvider_->init( stepsize );
						return *subProvider_;
					}

					int timeStep () const
					{
						return timeStep_ * 3 + currentStepType_;
					}

				protected:
					void next ( const double timeStep )
					{
						currentStepType_ = StokesStepA;
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
