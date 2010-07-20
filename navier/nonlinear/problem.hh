#ifndef  DUNE_PROBLEM_HH__
#define  DUNE_PROBLEM_HH__


namespace Dune {

namespace NavierStokes {
	namespace NonlinearStep {
		template < class DiscreteVelocityFunctionImp >
		class ProblemAdapter
		{
			public:
				typedef DiscreteVelocityFunctionImp
					DiscreteVelocityFunctionType;
				enum { ConstantVelocity = false };
				enum { dimDomain = GridType::dimensionworld };
				typedef typename DiscreteVelocityFunctionType::DomainType
					DomainType;
				typedef typename DiscreteVelocityFunctionType::RangeType
					RangeType;
				typedef typename DomainType::field_type
					DomainFieldType;
				typedef typename RangeType::field_type
					RangeFieldType;


			/**
			* @brief define problem parameters
		   */
			ProblemAdapter( const DiscreteVelocityFunctionType initialVelocity )
				: initialVelocity_( initialVelocity ),
				startTime_(Parameter::getValue<double>("fem.timeprovider.starttime",0.0)),
				epsilon(Parameter::getValue<double>("femhowto.epsilon",0.1))
			{

			}


			/**
		   * @brief getter for the velocity
		   */
			void velocity(const DomainType& x, RangeType& v) const {
				NEEDS_IMPLEMENTATION
			}

			/**
		   * @brief evaluates \f$ u_0(x) \f$
		   */
			void evaluate(const DomainType& arg, RangeType& res) const {             /*@LST0@@LST1@*/
				evaluate(arg, startTime_, res);
			}

			/**
			   * @brief old version of the exact solution
			   *
			   * old version of evaluate(const DomainType& arg, double t, RangeType& res),
			   * which is still needed by the DataWriter
			   */
			inline void evaluate(double t,  const DomainType& arg, RangeType& res) const {
				evaluate(arg, t, res);
			}

			/**
			   * @brief evaluate exact solution
			   */
			void evaluate(const DomainType& arg, double t, RangeType& res) const
			{
				NEEDS_IMPLEMENTATION
			}

			/**
			   * @brief latex output for EocOutput
			   */
			std::string description()
			{
				NEEDS_IMPLEMENTATION
				return std::string();
			}

		private:
			const DiscreteVelocityFunctionType& initialVelocity_;
			double startTime_;
		public:
			double epsilon;
			std::string myName;
		};
	}//end namespace NonlinearStep
}//end namespace NavierStokes

}//end namespace dune
#endif  /*DUNE_PROBLEM_HH__*/
