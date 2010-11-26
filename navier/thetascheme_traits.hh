#ifndef THETASCHEME_TRAITS_HH
#define THETASCHEME_TRAITS_HH

#include <dune/common/fixedarray.hh>
#include <dune/stokes/discretestokesmodelinterface.hh>
#include <dune/stokes/stokespass.hh>
#include <dune/navier/fractionaltimeprovider.hh>
#include <dune/navier/stokestraits.hh>
#include <dune/navier/testdata.hh>
#include <dune/stuff/functions.hh>
#include <dune/stuff/misc.hh>

namespace Dune {
	namespace NavierStokes {
		//! for each step keep a set of theta values and one value for dt
		// thetas_[stepnumber][theta_subscript_index]
		template < int numberOfSteps >
		struct ThetaSchemeDescription {
			typedef ThetaSchemeDescription< numberOfSteps >
				ThisType;
			static const int numberOfSteps_ = numberOfSteps;
			typedef Dune::array< double, 4 >
				ThetaValueArray;
			typedef Dune::array< ThetaValueArray, numberOfSteps >
				ThetaArray;
			typedef Stuff::wraparound_array< double, numberOfSteps >
				TimestepArray;
			ThetaArray thetas_;
			TimestepArray step_sizes_;

			ThetaSchemeDescription(){}

			ThetaSchemeDescription( const ThetaArray& a, const TimestepArray& s )
				:thetas_( a ), step_sizes_( s )
			{}

			ThetaSchemeDescription( const ThetaArray& a, double dt )
				: thetas_( a )
			{
				Stuff::fill_entirely( step_sizes_, dt );
			}



			static ThetaSchemeDescription<1> crank_nicholson( double delta_t )
			{
				ThetaValueArray c;
				Stuff::fill_entirely( c, 0.5f );
				ThetaArray a;
				Stuff::fill_entirely( a, c );
				return ThisType ( a, delta_t );
			}
			static ThetaSchemeDescription<1> forward_euler( double delta_t )
			{
				ThetaValueArray c = { 1.0 ,  0.0 ,  0.0 ,  1.0f  }  ;
				ThetaArray a;
				Stuff::fill_entirely( a, c );
				return ThetaSchemeDescription<1> ( a, delta_t );
			}
			static ThetaSchemeDescription<3> fs0( double delta_t )
			{
				const double theta			= 1 - (std::sqrt(2)/2.0f);
				const double theta_squigly	= 1 - ( 2 * theta );
				const double tau			= theta_squigly / ( 1 - theta );
				const double eta			= 1 - tau;
				typedef ThetaSchemeDescription<3>
					ReturnType;
				ReturnType::ThetaValueArray step_one	= { tau * theta, eta * theta, eta * theta, tau * theta };
				ReturnType::ThetaValueArray step_two	= { eta * theta_squigly, tau * theta_squigly, tau * theta_squigly, eta * theta_squigly };
				ReturnType::ThetaValueArray step_three	= { tau * theta, eta * theta, eta * theta, tau * theta };
				ReturnType::ThetaArray a;
				a[0] = step_one;
				a[1] = step_two;
				a[2] = step_three;
				ReturnType::TimestepArray c;
				c[0] = theta * delta_t;
				c[1] = theta_squigly * delta_t;
				c[2] = theta * delta_t;
				return ThisType ( a, c );
			}
		};

		template <class Stream, int N>
	    inline Stream& operator<< (Stream& s, ThetaSchemeDescription<N> desc )
	    {
			s << boost::format("%d-step theta-scheme description:\n") % N;
			for ( int i = 0; i < N; ++i ) {
				s << boost::format ( "dt_k: %e\ttheta1: %e\ttheta2: %e\ttheta3: %e\ttheta4: %e\n")//ideally this could be re-used, but that result in exception 'too-few-args"..
							% desc.step_sizes_[i]
							% desc.thetas_[i][0]
							% desc.thetas_[i][1]
							% desc.thetas_[i][2]
							% desc.thetas_[i][3];
			}
	        return s;
	    }

		template <	class CommunicatorImp,
					class GridPartImp,
					template < class > class AnalyticalForceImp,
					template < class > class AnalyticalDirichletDataImp,
					template < class,class > class ExactPressureImp,
					template < class,class > class ExactVelocityImp,
					int subStepCount,
					int gridDim, int sigmaOrder, int velocityOrder = sigmaOrder, int pressureOrder = sigmaOrder >
		struct ThetaSchemeTraits {
			typedef ThetaSchemeTraits<	CommunicatorImp,
										GridPartImp,
										AnalyticalForceImp,
										AnalyticalDirichletDataImp,
										ExactPressureImp,
										ExactVelocityImp,
										subStepCount,
										gridDim, sigmaOrder, velocityOrder , pressureOrder >
				ThisType;

			typedef GridPartImp
				GridPartType;

			static const int substep_count = subStepCount;
			typedef ThetaSchemeDescription< subStepCount >
				ThetaSchemeDescriptionType;
			typedef FractionalTimeProvider< ThetaSchemeDescriptionType, CommunicatorImp>
				TimeProviderType;

			typedef NonlinearStep::DiscreteStokesModelTraits<
						TimeProviderType,
						GridPartType,
						AnalyticalForceImp,
						OseenStep::ForceAdapterFunction,
						AnalyticalDirichletDataImp,
						typename ThetaSchemeDescriptionType::ThetaValueArray,
						gridDim,
						sigmaOrder,
						velocityOrder,
						pressureOrder >
				OseenModelTraits;
			typedef typename OseenModelTraits::ForceAdatperType
				OseenForceAdapterFunctionType;

			typedef ExactPressureImp< typename OseenModelTraits::PressureFunctionSpaceType,
									  TimeProviderType >
				ExactPressureType;
			typedef ExactVelocityImp< typename OseenModelTraits::VelocityFunctionSpaceType,
									  TimeProviderType >
				ExactVelocityType;

			typedef ExactSolution<ThisType>
				ExactSolutionType;

			typedef typename OseenModelTraits::PressureFunctionSpaceType
				PressureFunctionSpaceType;
			typedef typename OseenModelTraits::VelocityFunctionSpaceType
				VelocityFunctionSpaceType;

			typedef typename OseenModelTraits::DiscreteStokesFunctionSpaceWrapperType
				DiscreteStokesFunctionSpaceWrapperType;

			typedef typename OseenModelTraits::DiscreteStokesFunctionWrapperType
				DiscreteStokesFunctionWrapperType;
			typedef typename OseenModelTraits::RealAnalyticalForceType
				AnalyticalForceType;
			typedef typename OseenModelTraits::AnalyticalDirichletDataType
				AnalyticalDirichletDataType;

			typedef Dune::StartPass< DiscreteStokesFunctionWrapperType, -1 >
				StokesStartPassType;

			typedef CommunicatorImp
				CommunicatorType;

			typedef Dune::DiscreteStokesModelDefault< OseenModelTraits >
				OseenModelType;
			typedef Dune::StokesPass< OseenModelType,StokesStartPassType, 0 >
				OseenPassType;
		};
	}
}

#endif // THETASCHEME_TRAITS_HH
