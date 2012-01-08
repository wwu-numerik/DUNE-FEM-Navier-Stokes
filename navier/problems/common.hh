#ifndef NAVIER_COMMON_HH
#define NAVIER_COMMON_HH

#include <string>
#include <sstream>

#ifndef ALLGOOD_SETUPCHECK
#define ALLGOOD_SETUPCHECK struct SetupCheck { \
    template < typename ...Types > \
    bool operator()( const Types&... args ) { return true; } \
    std::string error() { return "";} }
#endif
#endif // NAVIER_COMMON_HH
