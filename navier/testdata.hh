#ifndef TESTDATA_HH
#define TESTDATA_HH

#include <dune/stuff/timefunction.hh>

namespace Dune {
	namespace NavierStokes {
		namespace TestCase3D {
			template < class FunctionSpaceImp >
			class Force : public Function < FunctionSpaceImp , Force < FunctionSpaceImp > >
			{
				  public:
					  typedef Force< FunctionSpaceImp >
						  ThisType;
					  typedef Function < FunctionSpaceImp ,ThisType >
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
					  inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}

				  private:
					  const double viscosity_;
					  const double alpha_;
					  static const int dim_ = FunctionSpaceImp::dimDomain;
			};

			template < class DomainType, class RangeType >
			void VelocityEvaluate( const double parameter_a, const double parameter_d, const double time, const DomainType& arg, RangeType& ret)
			{
				const double x		= arg[0];
				const double y		= arg[1];
				const double z		= arg[2];
				const double e_d_d_t	= std::exp( -1 * parameter_d * parameter_d *  time );
			  #define PM +
				ret[0] = - parameter_a * (
						( std::exp( parameter_a * x ) * std::sin( parameter_a * y PM/*+-*/ parameter_d * z ) ) +
						( std::exp( parameter_a * z ) * std::cos( parameter_a * x PM/*+-*/ parameter_d * y ) )
					  ) * e_d_d_t;
				ret[1] = - parameter_a * (
						( std::exp( parameter_a * y ) * std::sin( parameter_a * z PM/*+-*/ parameter_d * x ) ) +
						( std::exp( parameter_a * x ) * std::cos( parameter_a * y PM/*+-*/ parameter_d * z ) )
					  ) * e_d_d_t;
				ret[2] = - parameter_a * (
						( std::exp( parameter_a * z ) * std::sin( parameter_a * x PM/*+-*/ parameter_d * y ) ) +
						( std::exp( parameter_a * y ) * std::cos( parameter_a * z PM/*+-*/ parameter_d * x ) )
					  ) * e_d_d_t;
			  #undef PM
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
			class DirichletData : public Function < FunctionSpaceImp , DirichletData < FunctionSpaceImp > >
			{
				public:
					typedef DirichletData< FunctionSpaceImp >
						ThisType;
					typedef Function< FunctionSpaceImp, ThisType >
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
						Dune::CompileTimeChecker< ( dim_ == 3 ) > DirichletData_Unsuitable_WorldDim;
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
			class Velocity : public TimeFunction < FunctionSpaceImp , Velocity< FunctionSpaceImp,TimeProviderImp >, TimeProviderImp >
			{
				public:
					typedef Velocity< FunctionSpaceImp, TimeProviderImp >
						ThisType;
					typedef TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
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
						Dune::CompileTimeChecker< ( dim_ == 3 ) > DirichletData_Unsuitable_WorldDim;
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
			class Pressure : public TimeFunction <	FunctionSpaceImp ,
													Pressure < FunctionSpaceImp,TimeProviderImp >,
													TimeProviderImp >
			{
				public:
					typedef Pressure< FunctionSpaceImp, TimeProviderImp >
						ThisType;
					typedef TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
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
						Dune::CompileTimeChecker< ( dim_ == 3 ) > Pressure_Unsuitable_WorldDim;
						const double x			= arg[0];
						const double y			= arg[1];
						const double z			= arg[2];
						const double e_d_d_t	= std::exp( -2 * parameter_d_*parameter_d_ *  time );
					#define PM +
						ret = ( - parameter_a_ * parameter_a_/ 2.0 ) * (
							  std::exp( 2.0 * parameter_a_ * x ) +
							  std::exp( 2.0 * parameter_a_ * y ) +
							  std::exp( 2.0 * parameter_a_ * z ) +
							  ( 2 * std::sin( parameter_a_ * x PM parameter_d_ * y ) *
									std::cos( parameter_a_ * z PM parameter_d_ * x ) *
									std::exp( parameter_a_ * ( y + z ) ) ) +
							  ( 2 * std::sin( parameter_a_ * y PM parameter_d_ * z ) *
									std::cos( parameter_a_ * x PM parameter_d_ * y ) *
									std::exp( parameter_a_ * ( z + x ) ) ) +
							  ( 2 * std::sin( parameter_a_ * z PM parameter_d_ * x ) *
									std::cos( parameter_a_ * y PM parameter_d_ * z ) *
									std::exp( parameter_a_ * ( x + y ) ) )
							) * e_d_d_t;
					#undef PM
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

		}//end namespace TestCase3D

		namespace TestCase2D {
			template < class FunctionSpaceImp >
			class Force : public Function < FunctionSpaceImp , Force < FunctionSpaceImp > >
			{
				  public:
					  typedef Force< FunctionSpaceImp >
						  ThisType;
					  typedef Function < FunctionSpaceImp ,ThisType >
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
			class DirichletData : public Function < FunctionSpaceImp , DirichletData < FunctionSpaceImp > >
			{
				public:
					typedef DirichletData< FunctionSpaceImp >
						ThisType;
					typedef Function< FunctionSpaceImp, ThisType >
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
			class Velocity : public TimeFunction < FunctionSpaceImp , Velocity< FunctionSpaceImp,TimeProviderImp >, TimeProviderImp >
			{
				public:
					typedef Velocity< FunctionSpaceImp, TimeProviderImp >
						ThisType;
					typedef TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
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
			class Pressure : public TimeFunction <	FunctionSpaceImp ,
													Pressure < FunctionSpaceImp,TimeProviderImp >,
													TimeProviderImp >
			{
				public:
					typedef Pressure< FunctionSpaceImp, TimeProviderImp >
						ThisType;
					typedef TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
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
						const double x			= arg[0];
						const double y			= arg[1];
						const double e_minus_4_t	= std::exp( -4 * time );

						ret[0] = -0.25 * (
											std::cos( 2 * x ) + std::cos( 2 * y )
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

		}//end namespace TestCase2D

		namespace TrivialTestCase {
			template < class FunctionSpaceImp >
			class Force : public Function < FunctionSpaceImp , Force < FunctionSpaceImp > >
			{
				  public:
					  typedef Force< FunctionSpaceImp >
						  ThisType;
					  typedef Function < FunctionSpaceImp ,ThisType >
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
			void VelocityEvaluate( const double /*parameter_a*/, const double /*parameter_d*/, const double /*time*/, const DomainType& /*arg*/, RangeType& ret)
			{
				ret = RangeType( 0 );
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
			class DirichletData : public Function < FunctionSpaceImp , DirichletData < FunctionSpaceImp > >
			{
				public:
					typedef DirichletData< FunctionSpaceImp >
						ThisType;
					typedef Function< FunctionSpaceImp, ThisType >
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
			class Velocity : public TimeFunction < FunctionSpaceImp , Velocity< FunctionSpaceImp,TimeProviderImp >, TimeProviderImp >
			{
				public:
					typedef Velocity< FunctionSpaceImp, TimeProviderImp >
						ThisType;
					typedef TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
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
			class Pressure : public TimeFunction <	FunctionSpaceImp ,
													Pressure < FunctionSpaceImp,TimeProviderImp >,
													TimeProviderImp >
			{
				public:
					typedef Pressure< FunctionSpaceImp, TimeProviderImp >
						ThisType;
					typedef TimeFunction< FunctionSpaceImp, ThisType, TimeProviderImp >
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

					void evaluateTime( const double /*time*/, const DomainType& /*arg*/, RangeType& ret ) const
					{
						ret = RangeType( 0 );
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

		}//end namespace TrivialTestCase
	}//end namespace NavierStokes
}//end namespace Dune
#endif // TESTDATA_HH
