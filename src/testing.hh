#ifndef NAVIER_TESTING_HH
#define NAVIER_TESTING_HH

#include <dune/stuff/timefunction.hh>

namespace Testing {
namespace AdapterFunctions {

	template < class FunctionSpaceImp >
	class Force : public Dune::Function < FunctionSpaceImp , Force < FunctionSpaceImp > >
	{
		  public:
			  typedef Force< FunctionSpaceImp >
				  ThisType;
			  typedef Dune::Function < FunctionSpaceImp ,ThisType >
				  BaseType;
			  typedef typename BaseType::DomainType
				  DomainType;
			  typedef typename BaseType::RangeType
				  RangeType;

			  /**
			   *  \brief  constructor
			   *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
			   **/
			  Force( const double viscosity, const FunctionSpaceImp& space, const double alpha = 0.0 )
				  : BaseType ( space ),
					viscosity_( viscosity ),
					alpha_( alpha )
			  {}

			  /**
			   *  \brief  destructor
			   *  doing nothing
			   **/
			  ~Force()
			  {}

			  /**
			   *  \brief  evaluates the force
			   *  \param  arg
			   *          point to evaluate at
			   *  \param  ret
			   *          value of force at given point
			   **/
			  inline void evaluate( const double /*time*/, const DomainType& /*arg*/, RangeType& ret ) const
			  {
				  ret = RangeType(0);
			  }
			  inline void evaluate( const DomainType& /*arg*/, RangeType& ret ) const {ret = RangeType(0);}

		  private:
			  const double viscosity_;
			  const double alpha_;
			  static const int dim_ = FunctionSpaceImp::dimDomain;
	};

	template < class DomainType, class RangeType >
	void VelocityEvaluate( const double /*parameter_a*/, const double /*parameter_d*/, const double time, const DomainType& arg, RangeType& ret)
	{
		const double x				= arg[0];
		const double y				= arg[1];
		const double e_minus_2_t	= std::exp( -2 * time );

		ret[0] = -1 *	std::cos( x ) * std::sin( y ) * e_minus_2_t;
		ret[1] =		std::sin( x ) * std::cos( y ) * e_minus_2_t;
	}

	/**
	*  \brief  describes the dirichlet boundary data
	*
	*  \tparam DirichletTraitsImp
	*          types like functionspace, range type, etc
	*
	*  \todo   extensive docu with latex
	**/
	template < class FunctionSpaceImp >
	class DirichletData : public Dune::Function < FunctionSpaceImp , DirichletData < FunctionSpaceImp > >
	{
		public:
			typedef DirichletData< FunctionSpaceImp >
				ThisType;
			typedef Dune::Function< FunctionSpaceImp, ThisType >
				BaseType;
			typedef typename BaseType::DomainType
				DomainType;
			typedef typename BaseType::RangeType
				RangeType;

			/**
			*  \brief  constructor
			*
			*  doing nothing besides Base init
			**/
			DirichletData( const FunctionSpaceImp& space,
						 const double parameter_a = M_PI /2.0 ,
						 const double parameter_d = M_PI /4.0)
				: BaseType( space ),
				parameter_a_( parameter_a ),
				parameter_d_( parameter_d )
			{}

			/**
			*  \brief  destructor
			*
			*  doing nothing
			**/
			~DirichletData()
			{}

			template < class IntersectionType >
			void evaluate( const double time, const DomainType& arg, RangeType& ret, const IntersectionType& /*intersection */) const
			{
				Dune::CompileTimeChecker< ( dim_ == 2 ) > DirichletData_Unsuitable_WorldDim;
				VelocityEvaluate( parameter_a_, parameter_d_, time, arg, ret);
			}

			/**
			* \brief  evaluates the dirichlet data
			* \param  arg
			*         point to evaluate at
			* \param  ret
			*         value of dirichlet boundary data at given point
			**/
			inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}

		private:
			  static const int dim_ = FunctionSpaceImp::dimDomain ;
			  const double parameter_a_;
			  const double parameter_d_;
	};

	template < class FunctionSpaceImp, class TimeProviderImp >
	class Velocity : public Dune::TimeFunction < FunctionSpaceImp , Velocity< FunctionSpaceImp,TimeProviderImp >, TimeProviderImp >
	{
		public:
			typedef Velocity< FunctionSpaceImp, TimeProviderImp >
				ThisType;
			typedef Dune::TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
				BaseType;
			typedef typename BaseType::DomainType
				DomainType;
			typedef typename BaseType::RangeType
				RangeType;

			/**
			*  \brief  constructor
			*
			*  doing nothing besides Base init
			**/
			Velocity(	const TimeProviderImp& timeprovider,
						const FunctionSpaceImp& space,
						const double parameter_a = M_PI /2.0 ,
						const double parameter_d = M_PI /4.0)
				: BaseType( timeprovider, space ),
				parameter_a_( parameter_a ),
				parameter_d_( parameter_d )
			{}

			/**
			*  \brief  destructor
			*
			*  doing nothing
			**/
			~Velocity()
			{}

			void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
			{
				Dune::CompileTimeChecker< ( dim_ == 2 ) > DirichletData_Unsuitable_WorldDim;
				VelocityEvaluate( parameter_a_, parameter_d_, time, arg, ret);
			}

		   /**
			* \brief  evaluates the dirichlet data
			* \param  arg
			*         point to evaluate at
			* \param  ret
			*         value of dirichlet boundary data at given point
			**/
					inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}

		private:
			static const int dim_ = FunctionSpaceImp::dimDomain ;
			const double parameter_a_;
			const double parameter_d_;
	};

	template <	class FunctionSpaceImp,
				class TimeProviderImp >
	class Pressure : public Dune::TimeFunction <	FunctionSpaceImp ,
											Pressure < FunctionSpaceImp,TimeProviderImp >,
											TimeProviderImp >
	{
		public:
			typedef Pressure< FunctionSpaceImp, TimeProviderImp >
				ThisType;
			typedef Dune::TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
				BaseType;
			typedef typename BaseType::DomainType
				DomainType;
			typedef typename BaseType::RangeType
				RangeType;

		  /**
		   *  \brief  constructor
		   *
		   *  doing nothing besides Base init
		   **/
		  Pressure( const TimeProviderImp& timeprovider,
					const FunctionSpaceImp& space,
					const double parameter_a = M_PI /2.0 ,
					const double parameter_d = M_PI /4.0)
			  : BaseType( timeprovider, space ),
			  parameter_a_( parameter_a ),
			  parameter_d_( parameter_d )
		  {}

		  /**
		   *  \brief  destructor
		   *
		   *  doing nothing
		   **/
		   ~Pressure()
		   {}

			void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
			{
//				Dune::CompileTimeChecker< ( dim_ == 2 ) > Pressure_Unsuitable_WorldDim;
				const double x			= arg[0];
				const double y			= arg[1];

				ret[0] = std::sin( x );
			}

			inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}

		private:
			static const int dim_ = FunctionSpaceImp::dimDomain ;
			const double parameter_a_;
			const double parameter_d_;
	};

	template <	class FunctionSpaceImp,
				class TimeProviderImp >
	class PressureGradient : public Dune::TimeFunction <	FunctionSpaceImp ,
											PressureGradient < FunctionSpaceImp,TimeProviderImp >,
											TimeProviderImp >
	{
		public:
			typedef PressureGradient< FunctionSpaceImp, TimeProviderImp >
				ThisType;
			typedef Dune::TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
				BaseType;
			typedef typename BaseType::DomainType
				DomainType;
			typedef typename BaseType::RangeType
				RangeType;

		  /**
		   *  \brief  constructor
		   *
		   *  doing nothing besides Base init
		   **/
		  PressureGradient( const TimeProviderImp& timeprovider,
					const FunctionSpaceImp& space,
					const double parameter_a = M_PI /2.0 ,
					const double parameter_d = M_PI /4.0)
			  : BaseType( timeprovider, space ),
			  parameter_a_( parameter_a ),
			  parameter_d_( parameter_d )
		  {}

		  /**
		   *  \brief  destructor
		   *
		   *  doing nothing
		   **/
		   ~PressureGradient()
		   {}

			void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
			{
//				Dune::CompileTimeChecker< ( dim_ == 2 ) > Pressure_Unsuitable_WorldDim;
				const double x			= arg[0];
				const double y			= arg[1];

				ret[0] = std::cos( x );
				ret[1] = 0;
			}

			inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}

		private:
			static const int dim_ = FunctionSpaceImp::dimDomain ;
			const double parameter_a_;
			const double parameter_d_;
	};

	template < class FunctionSpaceImp, class TimeProviderImp >
	class VelocityLaplace : public Dune::TimeFunction < FunctionSpaceImp , VelocityLaplace< FunctionSpaceImp,TimeProviderImp >, TimeProviderImp >
	{
		public:
			typedef VelocityLaplace< FunctionSpaceImp, TimeProviderImp >
				ThisType;
			typedef Dune::TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
				BaseType;
			typedef typename BaseType::DomainType
				DomainType;
			typedef typename BaseType::RangeType
				RangeType;

			/**
			*  \brief  constructor
			*
			*  doing nothing besides Base init
			**/
			VelocityLaplace(	const TimeProviderImp& timeprovider,
						const FunctionSpaceImp& space,
						const double parameter_a = M_PI /2.0 ,
						const double parameter_d = M_PI /4.0)
				: BaseType( timeprovider, space ),
				parameter_a_( parameter_a ),
				parameter_d_( parameter_d )
			{}

			/**
			*  \brief  destructor
			*
			*  doing nothing
			**/
			~VelocityLaplace()
			{}

			void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
			{
				Dune::CompileTimeChecker< ( dim_ == 2 ) > DirichletData_Unsuitable_WorldDim;

			}

		   /**
			* \brief  evaluates the dirichlet data
			* \param  arg
			*         point to evaluate at
			* \param  ret
			*         value of dirichlet boundary data at given point
			**/
			inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}

		private:
			static const int dim_ = FunctionSpaceImp::dimDomain ;
			const double parameter_a_;
			const double parameter_d_;
	};

	template < class FunctionSpaceImp, class TimeProviderImp >
	class VelocityConvection : public Dune::TimeFunction < FunctionSpaceImp , VelocityConvection< FunctionSpaceImp,TimeProviderImp >, TimeProviderImp >
	{
		public:
			typedef VelocityConvection< FunctionSpaceImp, TimeProviderImp >
				ThisType;
			typedef Dune::TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
				BaseType;
			typedef typename BaseType::DomainType
				DomainType;
			typedef typename BaseType::RangeType
				RangeType;

			/**
			*  \brief  constructor
			*
			*  doing nothing besides Base init
			**/
			VelocityConvection(	const TimeProviderImp& timeprovider,
						const FunctionSpaceImp& space,
						const double parameter_a = M_PI /2.0 ,
						const double parameter_d = M_PI /4.0)
				: BaseType( timeprovider, space ),
				parameter_a_( parameter_a ),
				parameter_d_( parameter_d )
			{}

			/**
			*  \brief  destructor
			*
			*  doing nothing
			**/
			~VelocityConvection()
			{}

			void evaluateTime( const double time, const DomainType& arg, RangeType& ret ) const
			{
				Dune::CompileTimeChecker< ( dim_ == 2 ) > DirichletData_Unsuitable_WorldDim;

			}

		   /**
			* \brief  evaluates the dirichlet data
			* \param  arg
			*         point to evaluate at
			* \param  ret
			*         value of dirichlet boundary data at given point
			**/
			inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}

		private:
			static const int dim_ = FunctionSpaceImp::dimDomain ;
			const double parameter_a_;
			const double parameter_d_;
	};
}//end namespace TestCase2D
}//end namespace AdapterFunctions
#endif // TESTING_HH
