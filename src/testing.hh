#ifndef NAVIER_TESTING_HH
#define NAVIER_TESTING_HH

#include <dune/stuff/misc.hh>
#include <dune/stuff/timefunction.hh>
#include <dune/stuff/parametercontainer.hh>
#include <dune/common/tuples.hh>

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
				const double x				= arg[0];
				const double y				= arg[1];

				ret[0] = x*x + y;
				ret[1] = y*y + x;
			}

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
				const double x				= arg[0];
				const double y				= arg[1];

				ret[0] = 2;
				ret[1] = 2;
			}

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
				Dune::CompileTimeChecker< ( dim_ == 2 ) > VelocityConvection_Unsuitable_WorldDim;
				const double x				= arg[0];
				const double y				= arg[1];
				const double u_1 = x*x + y;
				const double u_2 = y*y + x;
				ret[0] = u_1 * 2 * x	+ u_2;
				ret[1] = u_1			+ u_2 * 2 * y;
			}

		private:
			static const int dim_ = FunctionSpaceImp::dimDomain ;
			const double parameter_a_;
			const double parameter_d_;
	};
}//end namespace TestCase2D

namespace AdapterFunctionsScalar {

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
				Dune::CompileTimeChecker< ( dim_ == 1 ) > DirichletData_Unsuitable_WorldDim;
				const double x				= arg[0];

				ret[0] = std::sin( x );
			}

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
				Dune::CompileTimeChecker< ( dim_ == 1 ) > Pressure_Unsuitable_WorldDim;
				const double x			= arg[0];
				const double y			= arg[1];

				ret[0] = std::sin( x );
			}

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
//				const double x			= arg[0];
//				const double y			= arg[1];

//				ret[0] = std::cos( x );
				ret[1] = 0;
			}

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
				Dune::CompileTimeChecker< ( dim_ == 1 ) > DirichletData_Unsuitable_WorldDim;
				const double x				= arg[0];
//				const double y				= arg[1];
				ret[0] =  - std::sin( x );
//				ret[1] = 2;
			}

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
				Dune::CompileTimeChecker< ( dim_ == 1 ) > DirichletData_Unsuitable_WorldDim;
				const double x				= arg[0];
				const double y				= arg[1];
				const double u_1 = x*x + y;
				const double u_2 = y*y + x;
				ret[0] = u_1 * 2 * x	+ u_2;
				ret[1] = u_1			+ u_2 * 2 * y;
			}

		private:
			static const int dim_ = FunctionSpaceImp::dimDomain ;
			const double parameter_a_;
			const double parameter_d_;
	};
}//end namespace TestCase2D

namespace AdapterFunctionsVectorial {

	static const double pi_factor = M_PI;//controls number of vortices
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
			  inline void evaluate( const double time, const DomainType& arg, RangeType& ret ) const
			  {
				  const double x			= arg[0];
				  const double y			= arg[1];
				  const double v			= viscosity_;
				  const double P			= pi_factor;
				  const double E			= std::exp( -2 * std::pow( P, 2 ) * viscosity_ * time );
				  const double F			= std::exp( -4 * std::pow( P, 2 ) * viscosity_ * time );
				  const double S_x			= std::sin( P * x );
				  const double S_y			= std::sin( P * y );
				  const double S_2x			= std::sin( 2 * P * x );
				  const double S_2y			= std::sin( 2 * P * y );
				  const double C_x			= std::cos( P * x );
				  const double C_y			= std::cos( P * y );
//						  ret[0] = - C_x * E * P * ( S_x * E + v * S_y * P )	+ 0.5 * P * F * S_2x;
//						  ret[1] = - C_y * E * P * ( S_y * E - v * S_x * P )	+ 0.5 * P * F * S_2y;

				  //diff
				  ret[0] = - 2 * v * C_x * S_y * E * P * P;
				  ret[1] =   2 * v * C_y * S_x * E * P * P;

				  //druck
				  ret[0] += 0.5 * P * F * S_2x;
				  ret[1] += 0.5 * P * F * S_2y;

				  //conv
				  RangeType conv;
				  VelocityConvectionEvaluateTime( time, arg, conv );
				  ret += conv;


				  //zeitableitung
				  RangeType u;
				  VelocityEvaluate( 0, 0, time, arg, u);
				  ret[0] += ( -2 * M_PI * M_PI * v ) * u[0];
				  ret[1] += ( -2 * M_PI * M_PI * v ) * u[1];

				  ret *=  Parameters().getParam( "rhs_factor", 1.0 );
			  }
			  inline void evaluate( const DomainType& /*arg*/, RangeType& ret ) const {assert(false);}

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
		const double v				= Parameters().getParam( "viscosity", 1.0 );
		const double e_minus_2_t	= std::exp( -2 * std::pow( pi_factor, 2 ) * v * time );

		ret[0] = -1 *	std::cos( pi_factor * x ) * std::sin( pi_factor * y ) * e_minus_2_t;
		ret[1] =		std::sin( pi_factor * x ) * std::cos( pi_factor * y ) * e_minus_2_t;
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
//					inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}

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
				Dune::CompileTimeChecker< ( dim_ == 2 ) > Pressure_Unsuitable_WorldDim;
				const double x				= arg[0];
				const double y				= arg[1];
				const double v				= Parameters().getParam( "viscosity", 1.0 );
				const double e_minus_4_t	= std::exp( -4 * std::pow( pi_factor, 2 ) * time * v );

				ret[0] = -0.25 * (
									std::cos( 2 * pi_factor * x ) + std::cos( 2 * pi_factor * y )
								) * e_minus_4_t;
			}

			/**
			* \brief  evaluates the dirichlet data
			* \param  arg
			*         point to evaluate at
			* \param  ret
			*         value of dirichlet boundary data at given point
			**/
//					inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}

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
				const double x			= arg[0];
				const double y			= arg[1];
				const double v			= Parameters().getParam( "viscosity", 1.0 );
				const double P			= pi_factor;
				const double E			= std::exp( -2 * std::pow( P, 2 ) * v * time );
				const double F			= std::exp( -4 * std::pow( P, 2 ) * v * time );
				const double S_x			= std::sin( P * x );
				const double S_y			= std::sin( P * y );
				const double S_2x			= std::sin( 2 * P * x );
				const double S_2y			= std::sin( 2 * P * y );
				const double C_x			= std::cos( P * x );
				const double C_y			= std::cos( P * y );
				ret[0] = 0.5 * P * F * S_2x;
				ret[1] = 0.5 * P * F * S_2y;
			}

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
//				Dune::CompileTimeChecker< ( dim_ == 2 ) > DirichletData_Unsuitable_WorldDim;
				const double x			= arg[0];
				const double y			= arg[1];
				const double v			= Parameters().getParam( "viscosity", 1.0 );
				const double P			= pi_factor;
				const double E			= std::exp( -2 * std::pow( P, 2 ) * v * time );
				const double F			= std::exp( -4 * std::pow( P, 2 ) * v * time );
				const double S_x			= std::sin( P * x );
				const double S_y			= std::sin( P * y );
				const double S_2x			= std::sin( 2 * P * x );
				const double S_2y			= std::sin( 2 * P * y );
				const double C_x			= std::cos( P * x );
				const double C_y			= std::cos( P * y );
				ret[0] = - 2 * C_x * E * P * (  v * S_y * P )	;
				ret[1] = - 2 * C_y * E * P * ( - v * S_x * P );

			}

		private:
			static const int dim_ = FunctionSpaceImp::dimDomain ;
			const double parameter_a_;
			const double parameter_d_;
	};
	template < class DomainType, class RangeType >
	void VelocityConvectionEvaluateTime( const double time, const DomainType& arg, RangeType& ret )
	{
//				Dune::CompileTimeChecker< ( dim_ == 2 ) > DirichletData_Unsuitable_WorldDim;

		const double x			= arg[0];
		const double y			= arg[1];
		const double v			= Parameters().getParam( "viscosity", 1.0 );;
		const double P			= pi_factor;
		const double E			= std::exp( -2 * std::pow( P, 2 ) * v * time );
		const double F			= std::exp( -4 * std::pow( P, 2 ) * v * time );
		const double S_x			= std::sin( P * x );
		const double S_y			= std::sin( P * y );
		const double S_2x			= std::sin( 2 * P * x );
		const double S_2y			= std::sin( 2 * P * y );
		const double C_x			= std::cos( P * x );
		const double C_y			= std::cos( P * y );
		ret[0] = - E * E *P * C_x * S_x ;//eigentlich richtig
		ret[1] = - E * E *P * S_y * C_y;

//		ret[0] =  E * E *P * C_x * S_x * ( C_y * C_y - S_y * S_y );//eigentlich falsch
//		ret[1] = - E * E *P * S_y * C_y * ( S_x * S_x - C_x * C_x );

	}
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
				VelocityConvectionEvaluateTime( time, arg, ret );
			}

		private:
			static const int dim_ = FunctionSpaceImp::dimDomain ;
			const double parameter_a_;
			const double parameter_d_;
	};
}//end namespace AdapterFunctionsVectorial

}//end namespace AdapterFunctions

template < class T1, class T2, class T3 >
struct TupleSerializer {
	typedef Dune::Tuple<	const typename T1::DiscreteVelocityFunctionType*,
							const typename T1::DiscretePressureFunctionType*,
							const typename T2::DiscreteVelocityFunctionType*,
							const typename T2::DiscretePressureFunctionType*,
							const typename T3::DiscreteVelocityFunctionType*,
							const typename T3::DiscretePressureFunctionType* >
		TupleType;

	static TupleType& getTuple( T1& t1,
								T2& t2,
								T3& t3 )
	{
		static TupleType t( &(t1.discreteVelocity()),
							&(t1.discretePressure()),
							&(t2.discreteVelocity()),
							&(t2.discretePressure()),
							&(t3.discreteVelocity()),
							&(t3.discretePressure()));
		return t;
	}
};

#endif // TESTING_HH
