#ifndef OSEEN_HH
#define OSEEN_HH

#include <dune/navier/global_defines.hh>
#include <dune/navier/problems.hh>

#include <dune/oseen/discretestokesmodelinterface.hh>
#include <dune/oseen/stokespass.hh>
#include <dune/navier/fractionaltimeprovider.hh>
#include <dune/navier/stokestraits.hh>
#include <dune/navier/exactsolution.hh>
#include <dune/navier/thetascheme_traits.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/navier/fractionaldatawriter.hh>
#include <dune/stuff/customprojection.hh>
#include <dune/common/collectivecommunication.hh>
#include <cmath>

namespace Dune {
namespace Oseen {

		template <	class TimeProviderType,
					class GridPartImp,
					template < class, class > class ForceFuntionType,
					template < class, class > class AnalyticalDirichletDataImp,
					int gridDim, int sigmaOrder, int velocityOrder = sigmaOrder, int pressureOrder = sigmaOrder >
		class DiscreteModelTraits
		{
			public:

				//! for CRTP trick
				typedef DiscreteStokesModelDefault < DiscreteModelTraits >
					DiscreteModelType;

				//! we use caching quadratures for the entities
				typedef Dune::CachingQuadrature< GridPartImp, 0 >
					VolumeQuadratureType;

				//! we use caching quadratures for the faces
				typedef Dune::CachingQuadrature< GridPartImp, 1 >
					FaceQuadratureType;

				//! polynomial order for the discrete sigma function space
				static const int sigmaSpaceOrder = sigmaOrder;
				//! polynomial order for the discrete velocity function space
				static const int velocitySpaceOrder = velocityOrder;
				//! polynomial order for the discrete pressure function space
				static const int pressureSpaceOrder = pressureOrder;

		//    private:

				//! function space type for the velocity
				typedef Dune::FunctionSpace< double, double, gridDim, gridDim >
					VelocityFunctionSpaceType;

				//! discrete function space type for the velocity
				typedef Dune::DiscontinuousGalerkinSpace<   VelocityFunctionSpaceType,
															GridPartImp,
															velocitySpaceOrder >
					DiscreteVelocityFunctionSpaceType;

				//! function space type for the pressure
				typedef Dune::FunctionSpace< double, double, gridDim, 1 >
					PressureFunctionSpaceType;

				//! discrete function space type for the pressure
				typedef Dune::DiscontinuousGalerkinSpace<   PressureFunctionSpaceType,
															GridPartImp,
															pressureSpaceOrder >
					DiscretePressureFunctionSpaceType;

			public:

				//! discrete function space wrapper type
				typedef Dune::DiscreteStokesFunctionSpaceWrapper< Dune::DiscreteStokesFunctionSpaceWrapperTraits<
							DiscreteVelocityFunctionSpaceType,
							DiscretePressureFunctionSpaceType > >
					DiscreteStokesFunctionSpaceWrapperType;

			private:

				//! discrete function type for the velocity
				typedef Dune::AdaptiveDiscreteFunction< typename DiscreteStokesFunctionSpaceWrapperType::DiscreteVelocityFunctionSpaceType >
					DiscreteVelocityFunctionType;

				//! discrete function type for the pressure
				typedef Dune::AdaptiveDiscreteFunction< typename DiscreteStokesFunctionSpaceWrapperType::DiscretePressureFunctionSpaceType >
					DiscretePressureFunctionType;

			public:

				//! discrete function wrapper type
				typedef Dune::DiscreteStokesFunctionWrapper< Dune::DiscreteStokesFunctionWrapperTraits<
							DiscreteStokesFunctionSpaceWrapperType,
							DiscreteVelocityFunctionType,
							DiscretePressureFunctionType > >
					DiscreteStokesFunctionWrapperType;

				//! function space type for sigma
				typedef Dune::MatrixFunctionSpace<  double,
													double,
													gridDim,
													gridDim,
													gridDim >
					SigmaFunctionSpaceType;

				//! discrete function space type for sigma
				typedef Dune::DiscontinuousGalerkinSpace<   SigmaFunctionSpaceType,
															GridPartImp,
															sigmaSpaceOrder >
					DiscreteSigmaFunctionSpaceType;

			public:

				//! discrete function type for sigma
				typedef Dune::AdaptiveDiscreteFunction< DiscreteSigmaFunctionSpaceType >
					DiscreteSigmaFunctionType;

				//! function type for the analytical force
				typedef ForceFuntionType < VelocityFunctionSpaceType,TimeProviderType >
					AnalyticalForceFunctionType;

				typedef AnalyticalForceFunctionType
					AnalyticalForceType;

				typedef AnalyticalDirichletDataImp<VelocityFunctionSpaceType,TimeProviderType>
					AnalyticalDirichletDataType;

				typedef DiscreteVelocityFunctionType
					ExtraDataFunctionType;
				/**
				 *  \name   types needed for the pass
				 *  \{
				 **/
				//! return type of the pass
				typedef DiscreteStokesFunctionWrapperType
					DestinationType;
				/**
				 *  \}
				 **/

		};

	template <	class CommunicatorImp,
				class GridPartImp,
				int gridDim, int sigmaOrder, int velocityOrder = sigmaOrder, int pressureOrder = sigmaOrder >
	struct Traits {
		typedef Traits<	CommunicatorImp,
						GridPartImp,
						gridDim, sigmaOrder,
						velocityOrder, pressureOrder >
			ThisType;
		typedef GridPartImp
			GridPartType;
		typedef Dune::NavierStokes::ThetaSchemeDescription<1>
			SchemeDescriptionType;
		typedef Dune::NavierStokes::FractionalTimeProvider<SchemeDescriptionType,CommunicatorImp>
			TimeProviderType;

		typedef DiscreteModelTraits<
					TimeProviderType,
					GridPartType,
					OSEEN_DATA_NAMESPACE::Force,
					OSEEN_DATA_NAMESPACE::DirichletData,
					gridDim,
					sigmaOrder,
					velocityOrder,
					pressureOrder >
			OseenModelTraits;
		typedef Dune::DiscreteStokesModelDefault< OseenModelTraits >
			OseenModelType;
        typedef Dune::StokesPass< OseenModelType >
			OseenPassType;

		typedef OSEEN_DATA_NAMESPACE::Pressure< typename OseenModelTraits::PressureFunctionSpaceType,
								  TimeProviderType >
			ExactPressureType;
		typedef OSEEN_DATA_NAMESPACE::Velocity< typename OseenModelTraits::VelocityFunctionSpaceType,
								  TimeProviderType >
			ExactVelocityType;
		typedef OSEEN_DATA_NAMESPACE::VelocityConvection< typename OseenModelTraits::VelocityFunctionSpaceType, TimeProviderType >
			ConvectionType;
		typedef Dune::NavierStokes::ExactSolution<ThisType>
			ExactSolutionType;
		typedef typename OseenModelType::DiscreteStokesFunctionWrapperType
			DiscreteStokesFunctionWrapperType;
		typedef typename OseenModelType::DiscreteStokesFunctionSpaceWrapperType
			DiscreteStokesFunctionSpaceWrapperType;
		typedef typename DiscreteStokesFunctionSpaceWrapperType::DiscretePressureFunctionSpaceType::FunctionSpaceType
			PressureFunctionSpaceType;
		typedef typename DiscreteStokesFunctionSpaceWrapperType::DiscreteVelocityFunctionSpaceType::FunctionSpaceType
			VelocityFunctionSpaceType;


	};

} //end namespace Oseen
}//end namespace Dune

#endif // OSEEN_HH
