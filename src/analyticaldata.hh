/**
 *  \file   analyticaldata.hh
 *  \brief  contains classes representing analytical data
 **/

#ifndef ANALYTICALDATA_HH
#define ANALYTICALDATA_HH

#include <cmath>

//#include <dune/common/fvector.hh>
//#include <dune/common/fmatrix.hh>


/**
 *  \todo   texdoc
 **/
template < class FunctionSpaceImp >
class Force : public Dune::Function < FunctionSpaceImp , Force < FunctionSpaceImp > >
{
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
        inline void evaluate( const DomainType& arg, RangeType& ret ) const
        {
			Dune::CompileTimeChecker< ( dim_ > 1 && dim_ < 4 ) > Force_Unsuitable_WorldDim;

			if ( dim_ == 2 ) {
#ifdef SIMPLE_PROBLEM
                ret[0] = 0.0;//arg[1];
                ret[1] = 0.0;//-1.0;//arg[0];
#elif defined(CONSTANT_PROBLEM)
                ret[0] = 0.0;//arg[1];
                ret[1] = -1.0;//arg[0];
#elif defined(ROTATE_PROBLEM)
                ret[0] = arg[1];
                ret[1] = -1.0 * arg[0];
#elif defined(GENRALIZED_STOKES_PROBLEM)
                const double x = arg[0];
                const double y = arg[1];
                const double tmp = alpha_ * std::cos( M_PI_2 * (x+y) ) + M_PI_2 * M_PI * std::cos( M_PI_2 * ( x + y ) ) + M_PI_2 * std::cos( M_PI_2 * ( x - y ) ) ;
                ret[0]  =      tmp;
                ret[1] =    -  tmp;
#elif defined(COCKBURN_PROBLEM)
                ret[0] = 0.0;//arg[1];
                ret[1] = 0.0;//arg[0];
#else
				LOGIC_ERROR
#endif
            }
            else if ( dim_ == 3 ) {
#ifdef SIMPLE_PROBLEM
                ret[0] = 0.0;//arg[1];
                ret[1] = 0.0;//-1.0;//arg[0];
                ret[2] = 0.0;
#elif defined(CONSTANT_PROBLEM)
                ret[0] = 0.0;//arg[1];
                ret[1] = 0.0;//arg[0];
                ret[2] = -1.0;//arg[0];
#elif defined(AORTA_PROBLEM)
                ret[0] = 0.0;//arg[1];
                ret[1] = 0.0;//-1.0;//arg[0];
                ret[2] = 0.0;
#else
				LOGIC_ERROR
#endif
            }
        }

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
        DirichletData( const FunctionSpaceImp& space )
			: BaseType( space )
        {}

        /**
         *  \brief  destructor
         *
         *  doing nothing
         **/
         ~DirichletData()
         {}

		template < class IntersectionType >
		void evaluate( const DomainType& arg, RangeType& ret, const IntersectionType& intersection ) const
        {
			const int id = intersection.boundaryId();
			Dune::CompileTimeChecker< ( dim_ > 1 && dim_ < 4 ) > DirichletData_Unsuitable_WorldDim;
			if ( dim_ == 2 ) {
                // some computations
                const double x1 = arg[0];
                const double x2 = arg[1];
#ifdef SIMPLE_PROBLEM
                ret[0] = 1.0;
                ret[1] = 0.0;
#elif defined(CONSTANT_PROBLEM)
                ret[0] = 0.0;
                ret[1] = 0.0;
#elif defined(ROTATE_PROBLEM)
                ret[0] = 0.0;
                ret[1] = 0.0;
#elif defined(POROSITY_PROBLEM)
                if ( id == 2 ) { // faces on inner hole
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( id == 3 ) { // bottom faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 1.0;
                }
                else if ( id == 4 ) { // right faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( id == 5 ) { // top faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 1.0;
                }
                else if ( id == 6 ) { // left faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
#elif defined(POROSITY_PROBLEM_WOIDS)
                const double x1 = arg[0];
                const double x2 = arg[1];
                if ( !( x2 > 0.0 ) ) { // bottom faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( !( x1 < 1.0 ) ) { // right faces
                    ret[ 0 ] = 1.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( !( x2 < 1.0 ) ) { // top faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( !( x1 > 0.0 ) ) { // left faces
                    ret[ 0 ] = 1.0;
                    ret[ 1 ] = 0.0;
                }
                else {
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
#elif defined(GENRALIZED_STOKES_PROBLEM)
                const double tmp = std::cos( ( M_PI_2 ) * ( x1 + x2 ) );
                ret[0] = tmp;
                ret[1] = -1.0 * tmp;
#elif defined(DARCY_PROBLEM)
                ret[0] = 0.0;
                ret[1] = 0.0;
#elif defined(MICRO_PROBLEM_X)
                if ( id == 2 ) { // faces on inner hole
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( id == 3 ) { // bottom faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( id == 4 ) { // right faces
                    ret[ 0 ] = 1.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( id == 5 ) { // top faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( id == 6 ) { // left faces
                    ret[ 0 ] = 1.0;
                    ret[ 1 ] = 0.0;
                }
#elif defined(MICRO_PROBLEM_Y)
                if ( id == 2 ) { // faces on inner hole
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( id == 3 ) { // bottom faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 1.0;
                }
                else if ( id == 4 ) { // right faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( id == 5 ) { // top faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 1.0;
                }
                else if ( id == 6 ) { // left faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
#elif defined(COCKBURN_PROBLEM)
                if ( id == 2 ) { // faces on inner hole
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else {
                    double exp_of_x1 = std::exp( x1 );
                    double sin_of_x2 = std::sin( x2 );
                    double cos_of_x2 = std::cos( x2 );
                    //return
                    ret[0] = x2 * cos_of_x2;
                    ret[0] += sin_of_x2;
                    ret[0] *= -1.0 * exp_of_x1;
                    ret[1] = exp_of_x1 * x2 * sin_of_x2;
                }
#else
				LOGIC_ERROR
#endif
            }
            else if ( dim_ == 3 ) {
#ifdef SIMPLE_PROBLEM
                ret[0] = 1.0;
                ret[1] = 0.0;
                ret[2] = 0.0;
#elif defined(CONSTANT_PROBLEM)
                ret[0] = 0.0;
                ret[1] = 0.0;
                ret[2] = 0.0;
#elif defined(ROTATE_PROBLEM)
                assert( !"ROTATE_PROBLEM not implemented in 3D!" );
#elif defined(POROSITY_PROBLEM)
                assert( !"POROSITY_PROBLEM not implemented in 3D!" );
#elif defined(POROSITY_PROBLEM_WOIDS)
                assert( !"POROSITY_PROBLEM_WOIDS not implemented in 3D!" );
#elif defined(GENRALIZED_STOKES_PROBLEM)
                assert( !"GENRALIZED_STOKES_PROBLEM not implemented in 3D!" );
#elif defined(DARCY_PROBLEM)
                assert( !"DARCY_PROBLEM not implemented in 3D!" );
#elif defined(AORTA_PROBLEM)

                switch ( id ) {
                    case 1: {
                        ret[0] = 0.0;//arg[1];
                        ret[1] = 0.0;//-1.0;//arg[0];
                        ret[2] = 0.0;
                        return;
                    }
                    case 2: {
                        ret[0] = 1000.0;//arg[1];
                        ret[1] = 1000.0;//-1.0;//arg[0];
                        ret[2] = 1000.0;
                        return;
                    }
                    case 6:
                    case 5:
                    case 4:
                    case 3: {
                        ret[0] = 1000.0;//arg[1];
                        ret[1] = 1000.0;//-1.0;//arg[0];
                        ret[2] = 1000.0;
                        return;
                    }
                    default:
                        assert( false );
                        return;
                }
#else
				LOGIC_ERROR
#endif
            }
        }

         /**
          * \brief  evaluates the dirichlet data
          * \param  arg
          *         point to evaluate at
          * \param  ret
          *         value of dirichlet boundary data at given point
          **/
        inline void evaluate( const DomainType& arg, RangeType& ret ) const {}

    private:
		static const int dim_ = FunctionSpaceImp::dimDomain ;
};

#endif // end of analyticaldata.hh
