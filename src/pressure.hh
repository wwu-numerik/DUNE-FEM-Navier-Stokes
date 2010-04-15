/**
 *  \file   pressure.hh
 *
 *  \brief  contains a class Pressure with traitsclass PressureTraits
 **/

#ifndef NAVIER_PRESSURE_HH
#define NAVIER_PRESSURE_HH

#include <cmath>

#include <dune/common/fvector.hh>

#include <dune/stuff/logging.hh>

/**
 *  \brief  containing typedefs needed by Pressure
 *
 *  \tparam gridDim (unused)
 *          dimension of the grid
 *
 *  \tparam PressureFunctionSpaceImp
 *          (continuous) FunctionSpace
 **/
template < int gridDim, class PressureFunctionSpaceImp >
class PressureTraits
{
    public:
        typedef PressureFunctionSpaceImp
            FunctionSpaceType;
        typedef typename FunctionSpaceType::DomainType
            DomainType;
        typedef typename FunctionSpaceType::RangeType
            RangeType;
        typedef typename FunctionSpaceType::JacobianRangeType
        //typedef Dune::FieldVector< double, gridDim >
            GradientRangeType;

};

/**
 *  \brief  describes the presure
 *
 *  \tparam PressureTraitsImp
 *          types like functionspace, range type, etc
 *
 *  \todo   extensive docu with latex
 **/
template < class PressureTraitsImp >
class Pressure : public Dune::Function < typename PressureTraitsImp::FunctionSpaceType, Pressure < PressureTraitsImp > >
{
    public:
        typedef PressureTraitsImp
            Traits;
        typedef typename Traits::DomainType
            DomainType;
        typedef typename Traits::RangeType
            RangeType;
        typedef typename Traits::GradientRangeType
            GradientRangeType;
        typedef typename Traits::FunctionSpaceType
            PressureFunctionSpaceType;
        typedef Dune::Function < typename PressureTraitsImp::FunctionSpaceType , Pressure < PressureTraitsImp > >
            BaseType;
        /**
         *  \brief  constructor
         *
         *  doing nothing besides Base init
         **/
        Pressure( const PressureFunctionSpaceType& press_space )
            : BaseType( press_space ),
            dim_( PressureTraitsImp::FunctionSpaceType::dimDomain )
        {}

        /**
         *  \brief  destructor
         *
         *  doing nothing
         **/
        ~Pressure()
        {}

        /**
         *  \brief  evaluates the pressure
         *
         *  \param  arg
         *          point to evaluate at
         *  \param  ret
         *          value of pressure at given point
         **/
        inline void evaluate( const DomainType& arg, RangeType& ret ) const
        {
            if ( dim_ == 1 ) {
                assert( !"pressure not implemented in 1D" );
            }
            else if ( dim_ == 2 ) {
                double x1 = arg[0];
                double x2 = arg[1];
#ifdef SIMPLE_PROBLEM
                ret[0] = -x2;
#elif defined(CONSTANT_PROBLEM)
                ret[0] = -x2;
#elif defined(GENRALIZED_STOKES_PROBLEM)
                ret[0] = std::sin( ( M_PI_2 ) * ( x1 - x2 ) );
#else
                ret[0] = 2.0 * std::exp( x1 ) * std::sin( x2 );
#endif
            }
            else if ( dim_ == 3 ) {
//                double x1 = arg[0];
//                double x2 = arg[1];
//                double x3 = arg[2];
#ifdef SIMPLE_PROBLEM
                assert( !"SIMPLE_PROBLEM not implemented in 3D" );
#elif defined(CONSTANT_PROBLEM)
                ret[0] = 0;
                ret[1] = 0;
                ret[2] = 0;
#elif defined(GENRALIZED_STOKES_PROBLEM)
                assert( !"GENRALIZED_STOKES_PROBLEM not implemented in 3D" );
#elif defined(AORTA_PROBLEM)
                ret[0] = 0.0;//arg[1];
                ret[1] = 0.0;//-1.0;//arg[0];
                ret[2] = 0.0;
#else
                assert( !"pressure not implemented in 3D" );
#endif
            }
            else {
                assert( !"pressure not implemented for more than 3 dimensions" );
            }
        }

        /**
         *  \brief  evaluates the pressure
         *
         *  \param  arg
         *          point to evaluate at
         *
         *  \return value of pressure at given point
         **/
        RangeType operator () ( const DomainType& arg)
        {
            RangeType ret;
            evaluate( arg, ret );
            return ret;
        }

        /**
         *  \brief  evaluates the gradient of the pressure
         *
         *  \param  arg
         *          point to evaluate at
         *  \param  ret
         *          value of gradient of the pressure at given point
         **/
//        inline void gradient(const DomainType& arg, GradientRangeType& ret ) const;

        /**
         *  \brief  a simple test of all class' functionalities
         **/
        void testMe() const;

    private:
        const int dim_;
 };

#endif // end of pressure.hh
