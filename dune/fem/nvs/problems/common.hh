#ifndef NAVIER_COMMON_HH
#define NAVIER_COMMON_HH

#include <string>
#include <sstream>
#include <dune/stuff/aliases.hh>
#include "runtimefunction.hh"

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

/** Copyright (c) 2012, Rene Milk 
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer. 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution. 
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies, 
 * either expressed or implied, of the FreeBSD Project.
**/

