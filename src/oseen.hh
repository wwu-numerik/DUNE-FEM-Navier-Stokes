#ifndef OSEEN_HH
#define OSEEN_HH

#include <dune/navier/global_defines.hh>
#include <dune/navier/problems.hh>

#include <dune/oseen/discreteoseenmodelinterface.hh>
#include <dune/oseen/oseenpass.hh>
#include <dune/navier/fractionaltimeprovider.hh>
#include <dune/navier/stokestraits.hh>
#include <dune/navier/exactsolution.hh>
#include <dune/navier/thetascheme_traits.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/navier/fractionaldatawriter.hh>
#include <dune/stuff/customprojection.hh>
#include <dune/common/collectivecommunication.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <cmath>

namespace Dune {
namespace Oseen {

		template <	class TimeProviderType,
                    class GridImp,
					template < class, class > class ForceFuntionType,
					template < class, class > class AnalyticalDirichletDataImp,
					int gridDim, int sigmaOrder, int velocityOrder = sigmaOrder, int pressureOrder = sigmaOrder >
		class DiscreteModelTraits
		{
			public:
                //! using DGAdaptiveLeafGridPart is mandated by DUNE-FEM, but not in any way checked...
                typedef Dune::DGAdaptiveLeafGridPart< GridImp >
                    GridPartType;
				//! for CRTP trick
				typedef DiscreteOseenModelDefault < DiscreteModelTraits >
					DiscreteModelType;

				//! we use caching quadratures for the entities
                typedef Dune::CachingQuadrature< GridPartType, 0 >
					VolumeQuadratureType;

				//! we use caching quadratures for the faces
                typedef Dune::CachingQuadrature< GridPartType, 1 >
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
                                                            GridPartType,
															velocitySpaceOrder >
					DiscreteVelocityFunctionSpaceType;

				//! function space type for the pressure
				typedef Dune::FunctionSpace< double, double, gridDim, 1 >
					PressureFunctionSpaceType;

				//! discrete function space type for the pressure
				typedef Dune::DiscontinuousGalerkinSpace<   PressureFunctionSpaceType,
                                                            GridPartType,
															pressureSpaceOrder >
					DiscretePressureFunctionSpaceType;

			public:

				//! discrete function space wrapper type
				typedef Dune::DiscreteOseenFunctionSpaceWrapper< Dune::DiscreteOseenFunctionSpaceWrapperTraits<
							DiscreteVelocityFunctionSpaceType,
							DiscretePressureFunctionSpaceType > >
					DiscreteOseenFunctionSpaceWrapperType;

				//! discrete function type for the velocity
				typedef Dune::AdaptiveDiscreteFunction< typename DiscreteOseenFunctionSpaceWrapperType::DiscreteVelocityFunctionSpaceType >
					DiscreteVelocityFunctionType;

				//! discrete function type for the pressure
				typedef Dune::AdaptiveDiscreteFunction< typename DiscreteOseenFunctionSpaceWrapperType::DiscretePressureFunctionSpaceType >
					DiscretePressureFunctionType;

			public:

				//! discrete function wrapper type
				typedef Dune::DiscreteOseenFunctionWrapper< Dune::DiscreteOseenFunctionWrapperTraits<
							DiscreteOseenFunctionSpaceWrapperType,
							DiscreteVelocityFunctionType,
							DiscretePressureFunctionType > >
					DiscreteOseenFunctionWrapperType;

				//! function space type for sigma
				typedef Dune::MatrixFunctionSpace<  double,
													double,
													gridDim,
													gridDim,
													gridDim >
					SigmaFunctionSpaceType;

				//! discrete function space type for sigma
				typedef Dune::DiscontinuousGalerkinSpace<   SigmaFunctionSpaceType,
                                                            GridPartType,
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
				typedef DiscreteOseenFunctionWrapperType
					DestinationType;
				/**
				 *  \}
				 **/

		};

	template <	class CommunicatorImp,
                class GridImp,
				int gridDim, int sigmaOrder, int velocityOrder = sigmaOrder, int pressureOrder = sigmaOrder >
	struct Traits {
		typedef Traits<	CommunicatorImp,
                        GridImp,
						gridDim, sigmaOrder,
						velocityOrder, pressureOrder >
			ThisType;
        //! using DGAdaptiveLeafGridPart is mandated by DUNE-FEM, but not in any way checked...
        typedef Dune::DGAdaptiveLeafGridPart< GridImp >
            GridPartType;
		typedef Dune::NavierStokes::ThetaSchemeDescription<1>
			SchemeDescriptionType;
		typedef Dune::NavierStokes::FractionalTimeProvider<SchemeDescriptionType,CommunicatorImp>
			TimeProviderType;

		typedef DiscreteModelTraits<
					TimeProviderType,
                    GridType,
					OSEEN_DATA_NAMESPACE::Force,
					OSEEN_DATA_NAMESPACE::DirichletData,
					gridDim,
					sigmaOrder,
					velocityOrder,
					pressureOrder >
			OseenModelTraits;
		typedef Dune::DiscreteOseenModelDefault< OseenModelTraits >
			OseenModelType;
        typedef Dune::OseenPass< OseenModelType >
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
		typedef typename OseenModelType::DiscreteOseenFunctionWrapperType
			DiscreteOseenFunctionWrapperType;
		typedef typename OseenModelType::DiscreteOseenFunctionSpaceWrapperType
			DiscreteOseenFunctionSpaceWrapperType;
		typedef typename DiscreteOseenFunctionSpaceWrapperType::DiscretePressureFunctionSpaceType::FunctionSpaceType
			PressureFunctionSpaceType;
		typedef typename DiscreteOseenFunctionSpaceWrapperType::DiscreteVelocityFunctionSpaceType::FunctionSpaceType
			VelocityFunctionSpaceType;


	};

} //end namespace Oseen
}//end namespace Dune

#endif // OSEEN_HH
