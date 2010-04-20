#ifndef TESTDATA_HH
#define TESTDATA_HH



namespace Dune {
	namespace NavierStokes {
		namespace TestCase {
			template < class FunctionSpaceImp >
					class Force : public Dune::Function < FunctionSpaceImp , Force < FunctionSpaceImp > > {
				  public:
					  typedef Force< FunctionSpaceImp >
						  ThisType;
					  typedef Dune::Function< FunctionSpaceImp, ThisType >
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

				/**
				*  \brief  describes the dirichlet boundary data
				*
				*  \tparam DirichletTraitsImp
				*          types like functionspace, range type, etc
				*
				*  \todo   extensive docu with latex
				**/
				template < class FunctionSpaceImp >
				class DirichletData : public Dune::Function < FunctionSpaceImp, DirichletData < FunctionSpaceImp > >
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
						  Dune::CompileTimeChecker< ( dim_ == 3 ) > DirichletData_Unsuitable_WorldDim;
						  const double x		= arg[0];
						  const double y		= arg[1];
						  const double z		= arg[2];
						  const double e_d_d_t	= std::exp( -1 * parameter_d_*parameter_d_ *  time );
						#define PM +
						  ret[0] = - parameter_a_ * (
								  ( std::exp( parameter_a_ * x ) * std::sin( parameter_a_ * y PM/*+-*/ parameter_d_ * z ) ) +
								  ( std::exp( parameter_a_ * z ) * std::cos( parameter_a_ * x PM/*+-*/ parameter_d_ * y ) )
								) * e_d_d_t;
						  ret[1] = - parameter_a_ * (
								  ( std::exp( parameter_a_ * y ) * std::sin( parameter_a_ * z PM/*+-*/ parameter_d_ * x ) ) +
								  ( std::exp( parameter_a_ * x ) * std::cos( parameter_a_ * y PM/*+-*/ parameter_d_ * z ) )
								) * e_d_d_t;
						  ret[2] = - parameter_a_ * (
								  ( std::exp( parameter_a_ * z ) * std::sin( parameter_a_ * x PM/*+-*/ parameter_d_ * y ) ) +
								  ( std::exp( parameter_a_ * y ) * std::cos( parameter_a_ * z PM/*+-*/ parameter_d_ * x ) )
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
					  inline void evaluate( const DomainType& arg, RangeType& ret ) const {assert(false);}

				  private:
					  static const int dim_ = FunctionSpaceImp::dimDomain ;
					  const double parameter_a_;
					  const double parameter_d_;
				};

		}//end namespace TestCase
	}//end namespace NavierStokes
}//end namespace Dune
#endif // TESTDATA_HH
