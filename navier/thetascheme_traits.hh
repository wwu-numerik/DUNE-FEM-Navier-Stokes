#ifndef THETASCHEME_TRAITS_HH
#define THETASCHEME_TRAITS_HH

#include <dune/common/fixedarray.hh>
#include <dune/stokes/discretestokesmodelinterface.hh>
#include <dune/stokes/stokespass.hh>
#include <dune/navier/fractionaltimeprovider.hh>
#include <dune/navier/stokestraits.hh>
#include <dune/navier/testdata.hh>
#include <dune/stuff/functions.hh>

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
			typedef Dune::array< double, numberOfSteps >
				TimestepArray;
			ThetaArray thetas_;
			TimestepArray step_sizes_;

			ThetaSchemeDescription(){}

			ThetaSchemeDescription( ThetaArray a, double dt )
				: thetas_( a )
			{
				Stuff::fill_entirely( step_sizes_, dt );
			}

			static ThisType crank_nicholson( double delta_t )
			{
				dune_static_assert( numberOfSteps_ == 1, "Crank Nicholson is a one step scheme" );
				ThetaValueArray c;
				Stuff::fill_entirely( c, 0.5f );
				ThetaArray a;
				Stuff::fill_entirely( a, c );
				return ThisType ( a, delta_t );
			}
			static ThisType forward_euler( double delta_t )
			{
				dune_static_assert( numberOfSteps_ == 1, "Forward Euler is a one step scheme" );
				ThetaValueArray c = { 1.0 ,  0.0 ,  0.0 ,  1.0f  }  ;
				ThetaArray a;
				Stuff::fill_entirely( a, c );
				return ThisType ( a, delta_t );
			}
		};

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
