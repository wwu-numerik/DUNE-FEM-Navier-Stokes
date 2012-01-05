#ifndef NAVIER_COMMON_HH
#define NAVIER_COMMON_HH

#include <string>
#include <sstream>

#define ALLGOOD_SETUPCHECK struct SetupCheck { \
    template < typename ...Types > \
    bool check( const Types&... args ) { return true; } \
    std::string error() { return "";} }

#endif // NAVIER_COMMON_HH
