#ifndef FRACTIONALTIMEPROVIDER_HH
#define FRACTIONALTIMEPROVIDER_HH

#include <dune/fem/solver/timeprovider.hh>


namespace Dune {
	namespace NavierStokes {

		typedef unsigned int
				StepType;
		const StepType StokesStepA		= 0;
		const StepType NonlinearStep	= 1;
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

				protected:
					using BaseType :: comm_;
					using BaseType :: cfl_;
					using BaseType :: dt_;
					using BaseType :: dtEstimate_;
					using BaseType :: dtUpperBound_;
					using BaseType :: valid_;
					using BaseType :: timeStep_;
					const double theta_;
					StepType currentStepType_;

				public:
					FractionalTimeProvider ( const double startTime,
								   const double theta,
								   const CommProvider &comm )
					: BaseType( startTime,
								Parameter :: getValidValue( "fem.timeprovider.factor", (double)1.0,
																	   ValidateGreater< double >( 0.0 ) ),
								comm ),
					  theta_( theta ),
					  currentStepType_( StokesStepA )
					{}

					const double subTime( ) const
					{
						switch ( currentStepType_ ) {
							case StokesStepA:
							default: return time();
							case NonlinearStep: return time()+ deltaT() * theta_ ;
							case StokesStepB: return time()+ deltaT() * (1 - theta_);
						}
					}

					void nextFractional()
					{
						if ( currentStepType_ == StokesStepB ) {
							next( deltaT() );
						}
						else
							++currentStepType_;
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
