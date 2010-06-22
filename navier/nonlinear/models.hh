#ifndef DUNE_NAVIERSTOKES_NONLINEAR_MODELS_HH
#define DUNE_NAVIERSTOKES_NONLINEAR_MODELS_HH

#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/solver/rungekutta.hh>
#include <dune/fem/io/parameter.hh>

namespace Dune {
	namespace NavierStokes {
		namespace NonlinearStep {
			// approximation order
			const int order   = POLORDER;
			const int rkSteps = POLORDER + 1;
		}//end namespace NonlinearStep
	}//end namespace NavierStokes
}//end namespace Dune

#include "dgoperator.hh"

#include "problem.hh"


namespace Dune {
	namespace NavierStokes {
		namespace NonlinearStep {
			/**
			 * @brief Traits class for AdvectionDiffusionModel
			 */
			template <	class GridPart,
						int dimRange2,
						int dimRange1=dimRange2* GridPart::GridType::dimensionworld >
			class AdvDiffModelTraits
			{
				public:
					typedef GridPart
						GridPartType;
					typedef typename GridPartType :: GridType
						GridType;

					enum { dimDomain = GridType::dimensionworld };
					enum { dimRange = dimRange2, dimGradRange = dimRange1 };

					// Definition of domain and range types
					typedef FieldVector<double, dimDomain>
						DomainType;
					typedef FieldVector<double, dimDomain-1>
						FaceDomainType;
					typedef FieldVector<double,dimRange>
						RangeType;
					typedef FieldVector<double,dimGradRange>
						GradientType;
					// ATTENTION: These are matrices (c.f. AdvectionDiffusionModel)
					typedef FieldMatrix<double,dimRange,dimDomain>
						FluxRangeType;
					typedef FieldMatrix<double,dimGradRange,dimDomain>
						DiffusionRangeType;
					typedef typename GridPartType::IntersectionIteratorType
						IntersectionIterator;
					typedef typename GridType::template Codim<0>::Entity
						EntityType;
			};

			/**
			 * @brief describes the analytical model
			 *
			 * This is an description class for the problem
			 * \f{eqnarray*}{ V + \nabla a(U)      & = & 0 \\
			 * \partial_t U + \nabla (F(U)+A(U,V)) & = & 0 \\
			 *                          U          & = & g_D \\
			 *                   \nabla U \cdot n  & = & g_N \f}
			 *
			 * where each class methods describes an analytical function.
			 * <ul>
			 * <li> \f$F\f$:   advection() </li>
			 * <li> \f$a\f$:   diffusion1() </li>
			 * <li> \f$A\f$:   diffusion2() </li>
			 * <li> \f$g_D\f$  boundaryValue() </li>
			 * <li> \f$g_N\f$  boundaryFlux1(), boundaryFlux2() </li>
			 * </ul>
			 *
			 * \attention \f$F(U)\f$ and \f$A(U,V)\f$ are matrix valued, and therefore the divergence is defined as
			 *
			 * \f[ \Delta M = \nabla \cdot (\nabla \cdot (M_{i\cdot})^t)_{i\in 1\dots n} \f]
			 *
			 * for a matrix \f$M\in \mathbf{M}^{n\times m}\f$.
			 *
			 * @param GridPart GridPart for extraction of dimension
			 * @param ProblemType Class describing the initial(t=0) and exact solution
			 */
			template <	class GridPartType,
						class ProblemType,
						class DiscreteStokesFunctionWrapperType,
						class AnalyticalForceType >
			class AdvectionDiffusionModel
			{
				public:
					enum { dimDomain = GridType::dimensionworld };
					enum { dimRange = GridType::dimensionworld };
					enum { ConstantVelocity = ProblemType :: ConstantVelocity };
					typedef typename GridPartType :: GridType
						GridType;
					typedef AdvDiffModelTraits<GridPartType,dimRange,dimRange*dimDomain>
						Traits;
					typedef typename Traits::DomainType
						DomainType;
					typedef typename Traits::RangeType
						RangeType;
					typedef typename Traits::GradientType
						GradientType;
					typedef typename Traits::FluxRangeType
						FluxRangeType;
					typedef typename Traits::DiffusionRangeType
						DiffusionRangeType;
				public:
				  /**
				   * @brief Constructor
				   *
				   * initializes model parameter
				   *
				   * @param problem Class describing the initial(t=0) and exact solution
				   */
					AdvectionDiffusionModel(const ProblemType& problem,
											const DiscreteStokesFunctionWrapperType& extraSource,
											const AnalyticalForceType& force)
						: problem_(problem),
						extraSource_(extraSource),
						force_(force),
						velocity_(0),
						epsilon(problem.epsilon)
					{
						// if diffusionTimeStep is set to non-zero in the parameterfile, the
						// deltaT in the timeprovider is updated according to the diffusion
						// parameter epsilon.
						bool diff_tstep;
						Parameter::get("femhowto.diffusionTimeStep", diff_tstep);
						tstep_eps = diff_tstep ? epsilon : 0;

						if(ConstantVelocity) {
							problem_.velocity(velocity_,velocity_);
						}
					}

					/**
					* @brief advection term \f$F\f$
					*
					* @param en entity on which to evaluate the advection term
					* @param time current time of TimeProvider
					* @param x coordinate local to entity
					* @param u \f$U\f$
					* @param f \f$f(U)\f$
					*/
					inline  void advection(const typename Traits::EntityType& en,
										 double time,
										 const DomainType& x,
										 const RangeType& u,
										 FluxRangeType & f) const
					{
						// evaluate velocity
						problem_.velocity(en.geometry().global(x),f[0]);
						// multiply with u
						//				f *= u;
						f[0] = u;
						NEEDS_IMPLEMENTATION
					}

					/**
					* @brief velocity calculation, is called by advection()
					*/
					inline  void velocity(const typename Traits::EntityType& en,
										double time,
										const DomainType& x,
										DomainType& v) const
					{
						problem_.velocity(en.geometry().global(x),v);
					}

					//! called by upwind flux calc
					inline  void velocity(const DomainType& x,
										double time,
										RangeType& v) const
					{
						problem_.evaluate(x,time,v);
					}
					/**
					* @brief diffusion term \f$a\f$
					*/
					inline  void diffusion1(typename Traits::EntityType& en,
										 double time,
										 const DomainType& x,
										 const RangeType& u,
										 DiffusionRangeType& a) const
					{
						a = 0;
						RangeType d = u;
						d *= std::sqrt(epsilon);
						for (int i=0;i<dimDomain;i++)
						//				  a[i][i]=d;
							a[i]=d;
						NEEDS_IMPLEMENTATION
					}

					/**
					* @brief diffusion term \f$A\f$
					*/
					inline double diffusion2(typename Traits::EntityType& en,
										  double time,
										  const DomainType& x,
										  const RangeType& u,
										  const GradientType& v,
										  FluxRangeType& A) const
					{
						//				A = v;
						//v ist FieldVector<9>, A Fieldmatrix <3,3>
						  NEEDS_IMPLEMENTATION
						A *= std::sqrt(epsilon);
						// QUESTION: Beieinflusst diese Größe das deltaT? Und wenn ja, wie?
						return tstep_eps;
					}

					/**
					* @brief checks for existence of dirichlet boundary values
					*/
					inline bool hasBoundaryValue(typename Traits::IntersectionIterator& it,
											   double time,
											   const typename Traits::FaceDomainType& x) const
					{
						return true;
					}

					/**
					* @brief neuman boundary values \f$g_N\f$ for pass2
					*/
					inline double boundaryFlux2(typename Traits::IntersectionIterator& it,
											 double time,
											 const typename Traits::FaceDomainType& x,
											 const RangeType& uLeft,
											 const GradientType& vLeft,
											 RangeType& gLeft) const
					{
						gLeft = 0.;
						return 0.;
					}

					/**
					* @brief neuman boundary values \f$g_N\f$ for pass1
					*/
					inline double boundaryFlux1(typename Traits::IntersectionIterator& it,
											 double time,
											 const typename Traits::FaceDomainType& x,
											 const RangeType& uLeft,
											 RangeType& gLeft) const
					{
						gLeft = 0.;
						return 0.;
					}

					/**
					* @brief dirichlet boundary values
					*/
					inline  void boundaryValue(typename Traits::IntersectionIterator& it,
											 double time,
											 const typename Traits::FaceDomainType& x,
											 const RangeType& uLeft,
											 RangeType& uRight) const
					{
						DomainType xgl=it.intersectionGlobal().global(x);
						problem_.evaluate(xgl,time,uRight);
					}


					/**
					* @brief return an estimate for deltaT that can be passed to the TimeProvider
					*/
					inline double diffusionTimeStep() const
					{
						return 2.*tstep_eps;
					}

					const DiscreteStokesFunctionWrapperType& extraSource() const
					{
					  return extraSource_;
					}

					const AnalyticalForceType& force() const
					{
					  return force_;
					}

				protected:
					const ProblemType& problem_;
					const DiscreteStokesFunctionWrapperType& extraSource_;
					const AnalyticalForceType& force_;
					public:
					mutable DomainType velocity_;
					protected:
					double epsilon;
					double tstep_eps;
			};

				/**
				 * @brief defines the advective flux
				 */
				template <class ModelType>
				class UpwindFlux {

				public:
				  typedef ModelType Model;
				  typedef typename Model::Traits Traits;
				  enum { dimRange = Model::dimRange };
				  typedef typename Model::DomainType DomainType;
				  typedef typename Model::RangeType RangeType;
				  typedef typename Model::FluxRangeType FluxRangeType;
				  typedef typename Model::DiffusionRangeType DiffusionRangeType;
				protected:
				  template <class Model, bool constVelo>
					struct Velocity
					{
					  /**
					   * @brief computes and returns the wind direction
					   */
					  static inline double upwind(const Model& model,
												  typename Traits::IntersectionIterator& it,
												  double time,
												  const typename Traits::FaceDomainType& x,
												  const RangeType& uLeft)
					  {
						const typename Traits::DomainType normal = it.integrationOuterNormal(x);
						RangeType velocity;
						const DomainType global_x = it.intersectionSelfLocal().global(x);
						model.velocity( global_x,
									   time,velocity);
						return normal*velocity;
					  }
					};

				  template <class Model>
					struct Velocity<Model,true>
					{
					  /**
					   * @brief computes and returns the wind direction for models with
					   * constant velocity
					   */
					  static inline double upwind(const Model& model,
												  typename Traits::IntersectionIterator& it,
												  double time,
												  const typename Traits::FaceDomainType& x,
												  const RangeType& uLeft)
					  {
						const typename Traits::DomainType normal = it.integrationOuterNormal(x);
						return normal * model.velocity_;
					  }
					};

				public:
				  /**
				   * @brief Constructor
				   */
				  UpwindFlux(const Model& mod) : model_(mod) {}

				  const Model& model() const {return model_;}

				  /**
				   * @brief evaluates the flux \f$g(u,v)\f$
				   *
				   * @return maximum wavespeed * normal
				   */
				  inline  double numericalFlux(typename Traits::IntersectionIterator& it,
											   double time,
											   const typename Traits::FaceDomainType& x,
											   const RangeType& uLeft,
											   const RangeType& uRight,
											   RangeType& gLeft,
											   RangeType& gRight) const
				  {
					const double upwind = Velocity<Model,Model::ConstantVelocity>::
					  upwind(model_,it,time,x,uLeft);

					if (upwind>0)
					  gLeft = uLeft;
					else
					  gLeft = uRight;
					gLeft *= upwind;
					gRight = gLeft;
					return std::abs(upwind);
				  }
				private:
				  mutable DomainType velocity_;
				  const Model& model_;
				};


				/************************************************/  /*@LST0@*/
				/* Definition of model and solver                *
				 ************************************************/
				template <	class GridPartImp,
							class DiscreteStokesFunctionWrapperImp,
							class AnalyticalForceType >
				struct Traits {

					// The initial function u_0 and the exact solution
					typedef ProblemAdapter<typename DiscreteStokesFunctionWrapperImp::DiscreteVelocityFunctionType> InitialDataType;
					// An analytical version of our model
					typedef AdvectionDiffusionModel<GridPartImp, InitialDataType, DiscreteStokesFunctionWrapperImp, AnalyticalForceType> ModelType;
					// The flux for the discretization of advection terms
					typedef UpwindFlux<ModelType> FluxType;
					// The DG Operator (using 2 Passes)
					typedef DGAdvectionDiffusionOperator<ModelType,UpwindFlux,order> DgType;
					// The ODE Solver
					typedef DuneODE::ExplicitRungeKuttaSolver<typename DgType::DestinationType> ODEType;

					// This is needed for the dataWriter that can write solutions to harddisk.
					typedef Tuple< typename DgType :: DestinationType * > IOTupleType;  /*@LST1@*/
				};

		}//end namespace NonlinearStep
	}//end namespace NavierStokes
}//end namespace Dune

#endif // MODELS_HH
