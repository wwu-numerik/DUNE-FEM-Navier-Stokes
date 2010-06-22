#ifndef DUNE_NAVIERSTOKES_NONLINEAR_UPWIND_HH
#define DUNE_NAVIERSTOKES_NONLINEAR_UPWIND_HH

namespace Dune {
	namespace NavierStokes {
		namespace NonlinearStep {

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



		}//end namespace NonlinearStep
	}//end namespace NavierStokes
}//end namespace Dune


#endif // UPWIND_HH
