#ifndef NAVIER_TESTING_HH
#define NAVIER_TESTING_HH

#include <dune/stuff/misc.hh>
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
//				Dune::CompileTimeChecker< ( dim_ == 2 ) > DirichletData_Unsuitable_WorldDim;
				const double x				= arg[0];
				const double y				= arg[1];
				const double z				= arg[2];

//				ret[0] = std::sin(x);// * std::cos( x );
//				ret[1] = std::sin(y);// * std::cos( x );;


				ret[0] = x*x ;
				ret[1] = 2*y*y;
				ret[2] = 4*z*z;
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
				const double z			= arg[2];

				ret[0] = std::sin(z) + std::sin(x) + std::sin(y);
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
				const double z			= arg[2];

				ret[0] = 2*x;
				ret[1] = 2*y;
				ret[2] = 2*z;



				ret[0] = std::cos(x);
				ret[1] = std::cos(y);
				ret[2] = std::cos(z);
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
				const double x				= arg[0];
				const double y				= arg[1];
				ret[0] =  2;
				ret[1] =  4;
				ret[2] =  8;
//				ret[0] =  - std::sin( x );
//				ret[1] =  - std::sin( y );

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
//				Dune::CompileTimeChecker< ( dim_ == 2 ) > DirichletData_Unsuitable_WorldDim;
				const double x				= arg[0];
				const double y				= arg[1];
				const double z				= arg[2];
				const double u_1 = x*x;
				const double u_2 = 2*y*y;
				const double u_3 = 4*z*z;
				ret[0] = u_1 * 2 * x	;
				ret[1] = u_2 * 4 * y	;
				ret[2] = u_3 * 8 * z	;
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
