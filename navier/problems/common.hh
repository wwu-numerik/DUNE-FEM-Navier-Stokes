#ifndef NAVIER_COMMON_HH
#define NAVIER_COMMON_HH

#include <string>
#include <sstream>

#ifndef ALLGOOD_SETUPCHECK
#define ALLGOOD_SETUPCHECK struct SetupCheck { \
    template < typename ...Types > \
    bool operator()( const Types&... /*args*/ ) { return true; } \
    std::string error() { return "";} }
#endif

#define NV_RUNTIME_FUNC(name) \
    template < class FunctionSpaceImp, class TimeProviderImp >\
    struct name : public Stuff::RuntimeFunction < FunctionSpaceImp, TimeProviderImp >\
    {\
            typedef Stuff::RuntimeFunction < FunctionSpaceImp, TimeProviderImp >\
                BaseType;\
            name(	const TimeProviderImp& timeprovider,\
                        const FunctionSpaceImp& /*space*/,\
                        const double /*parameter_a*/ = M_PI /2.0 ,\
                        const double /*parameter_d */= M_PI /4.0)\
                : BaseType( #name, timeprovider)\
            {}\
    }

#endif // NAVIER_COMMON_HH
