#ifndef FRACTIONALTIMEPROVIDER_HH
#define FRACTIONALTIMEPROVIDER_HH

#include <dune/fem/solver/timeprovider.hh>


namespace Dune {
	namespace NavierStokes {
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

				public:
					FractionalTimeProvider ( const double startTime,
								   const double theta,
								   const CommProvider &comm )
					: BaseType( startTime,
								Parameter :: getValidValue( "fem.timeprovider.factor", (double)1.0,
																	   ValidateGreater< double >( 0.0 ) ),
								comm ),
					  theta_( theta )
					{}

					const double subTime( const unsigned int fractionStep = 0 ) const
					{
						assert( fractionStep >= 0 && fractionStep <= 2 );
						switch ( fractionStep ) {
							default: return time();
							case 1: return time()+ deltaT() * theta_ ;
							case 2: return time()+ deltaT() * (1 - theta_);
						}
					}

		};
	}//end namespace NavierStokes
} //end namespace Dune

#endif // FRACTIONALTIMEPROVIDER_HH
