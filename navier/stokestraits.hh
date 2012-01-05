#ifndef STOKESTRAITS_HH
#define STOKESTRAITS_HH

#include <dune/stokes/discretestokesmodelinterface.hh>
#include <dune/navier/rhsadapter.hh>

namespace Dune {
	namespace NavierStokes {
		namespace NonlinearStep {
			template <	class TimeProviderType,
						class GridPartImp,
						template < class, class > class AnalyticalForceFunctionType,
						template < class, class, class, class> class ForceAdatperTemplateType,
						template < class, class > class AnalyticalDirichletDataImp,
						class ThetaValueArrayType,
						int gridDim, int sigmaOrder, int velocityOrder = sigmaOrder, int pressureOrder = sigmaOrder >
			class DiscreteStokesModelTraits
			{
				public:

					//! for CRTP trick
					typedef DiscreteStokesModelDefault < DiscreteStokesModelTraits >
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


                #if STOKES_USE_ISTL
                        //! discrete function type for the velocity
                        typedef Dune::BlockVectorDiscreteFunction< typename DiscreteStokesFunctionSpaceWrapperType::DiscreteVelocityFunctionSpaceType >
                            DiscreteVelocityFunctionType;

                        //! discrete function type for the pressure
                        typedef Dune::BlockVectorDiscreteFunction< typename DiscreteStokesFunctionSpaceWrapperType::DiscretePressureFunctionSpaceType >
                            DiscretePressureFunctionType;

                    public:
                        //! discrete function type for sigma
                        typedef Dune::BlockVectorDiscreteFunction< DiscreteSigmaFunctionSpaceType >
                            DiscreteSigmaFunctionType;
                #else
                        //! discrete function type for the velocity
                        typedef Dune::AdaptiveDiscreteFunction< typename DiscreteStokesFunctionSpaceWrapperType::DiscreteVelocityFunctionSpaceType >
                            DiscreteVelocityFunctionType;

                        //! discrete function type for the pressure
                        typedef Dune::AdaptiveDiscreteFunction< typename DiscreteStokesFunctionSpaceWrapperType::DiscretePressureFunctionSpaceType >
                            DiscretePressureFunctionType;

                    public:
                        //! discrete function type for sigma
                        typedef Dune::AdaptiveDiscreteFunction< DiscreteSigmaFunctionSpaceType >
                            DiscreteSigmaFunctionType;
                #endif

				public:

					//! discrete function wrapper type
					typedef Dune::DiscreteStokesFunctionWrapper< Dune::DiscreteStokesFunctionWrapperTraits<
								DiscreteStokesFunctionSpaceWrapperType,
								DiscreteVelocityFunctionType,
								DiscretePressureFunctionType > >
						DiscreteStokesFunctionWrapperType;

				public:
					//! function type for the analytical force
					typedef AnalyticalForceFunctionType< VelocityFunctionSpaceType,TimeProviderType >
						RealAnalyticalForceType;

					typedef ForceAdatperTemplateType<	TimeProviderType,
														RealAnalyticalForceType,
														DiscreteVelocityFunctionType,
														ThetaValueArrayType >
						ForceAdatperType;

					typedef ForceAdatperType
						AnalyticalForceType;

					//! function type for the analytical dirichlet data
//					typedef typename Dune::NavierStokes::DirichletAdapterFunctionTraits< AnalyticalDirichletDataImp, TimeProviderType >
//										::template Implementation<VelocityFunctionSpaceType,GridPartImp >
//							AnalyticalDirichletDataTraitsImplementation;
					typedef AnalyticalDirichletDataImp< VelocityFunctionSpaceType, TimeProviderType >
						AnalyticalDirichletDataType;

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
		} //namespace Nonlinear
	}//end namespace NavierStokes
} //end namespace Dune

#endif // STOKESTRAITS_HH
