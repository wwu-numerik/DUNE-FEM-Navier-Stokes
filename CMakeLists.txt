PROJECT(dune-navier-stokes)

#set minimum cmake version
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/../dune-stuff/cmake )
include(DuneUtils)
include(UseDoxygen)

FIND_PACKAGE(Boost 1.35.0 REQUIRED)
INCLUDE_DIR(${Boost_INCLUDE_DIR})

#----------------------------------------------------------------------------------------------------
# General Settings
#----------------------------------------------------------------------------------------------------
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
SET( CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".lib" ".la" ".a")
SET( CMAKE_COLOR_MAKEFILE ON)

EXECUTE_PROCESS(
	COMMAND ${CMAKE_SOURCE_DIR}/getRevision.sh
	OUTPUT_VARIABLE COMMIT
	ERROR_VARIABLE shell_error
	OUTPUT_STRIP_TRAILING_WHITESPACE
)

#----------------------------------------------------------------------------------------------------
# Options, that can be changed be the User in order to customise dune-navier-stokes
#----------------------------------------------------------------------------------------------------
SET( BLAS_LIB
	"blas" CACHE STRING
	"" )

SET( GRIDTYPE
	"ALUGRID_SIMPLEX" CACHE STRING
	"GRIDTYPE" )

SET( GRIDDIM
	"2" CACHE STRING
	"GRIDDIM" )

SET( TESTCASE
	"GreenTaylor" CACHE STRING
	"navier test data define" )

SET( POLORDER
	"2" CACHE STRING
	"Sigma POLORDER" )

SET( PRESSUREPOLORDER
	${POLORDER} CACHE STRING
	"Pressure POLORDER" )

SET( VELOCITYPOLORDER
	${POLORDER} CACHE STRING
	"Velocity POLORDER" )

SET( INNER_CG_SOLVERTYPE
	"OEMCGOp" CACHE STRING
	"INNER_CG_SOLVERTYPE" )

SET( OUTER_CG_SOLVERTYPE
	"OEMCGOp" CACHE STRING
	"OUTER_CG_SOLVERTYPE" )

SET( PROBLEM
	"COCKBURN_PROBLEM" CACHE STRING
	"PROBLEM" )

SET( OSEEN_DATA_NAMESPACE
	"TestCase2D" CACHE STRING
	"OSEEN_DATA_NAMESPACE" )

SET( DIRICHLET_DATA 
	"DirichletData" CACHE STRING
	"DIRICHLET_DATA implementation for default dirichlet traits" )

SET( ALUGRID_BASE_PATH
	"/share/dune/Modules/modules_x86_64/ALUGrid-1.22" CACHE STRING
	"ALUGRID_BASE_PATH" )

SET( UGGRID_BASE_PATH
	"/share/dune/Modules/modules_x86_64/ug" CACHE STRING
	"UGGRID_BASE_PATH" )

SET( GRAPE_PATH
	"/share/dune/Modules/modules_x86_64/grape" CACHE STRING
	"GRAPE_PATH" )

SET( CXX_FLAGS
	-Wall -O3 -fomit-frame-pointer -funroll-loops -DDNDEBUG -g -ggdb -DENABLE_ADAPTIVE=0 -DADAPTIVE_SOLVER -DUSE_BFG_CG_SCHEME CACHE STRING
	"")

SET( ENABLE_PARALELL
	0 CACHE BOOL
	"Enable parallel features" )

SET ( METIS_DIR
	"/share/dune/Modules/modules_x86_64/ParMetis-3.1.1" CACHE STRING
	"metis toplevel directory" )

#needs to come first to include correct config.h
INCLUDE_DIR( ${dune-navier-stokes_SOURCE_DIR} )

IF( ${GRIDTYPE} STREQUAL "UGGRID" )
	ADD_CXX_FLAGS( -I${UGGRID_BASE_PATH}/include -DENABLE_UG )
	LINK_DIRECTORIES( ${UGGRID_BASE_PATH}/lib/ )
	SET( UGLIBS "ugS3" "ugS2" "devS")
ENDIF( ${GRIDTYPE} STREQUAL "UGGRID" )

ADD_CXX_FLAGS(-DHAVE_CONFIG_H -DTESTCASE=${TESTCASE} -DOSEEN_DATA_NAMESPACE=${OSEEN_DATA_NAMESPACE} -fno-strict-aliasing )
ADD_CXX_FLAGS(-D${GRIDTYPE} -DGRIDDIM=${GRIDDIM} -D${PROBLEM} -DCOMMIT="${COMMIT}" -DDIRICHLET_DATA=${DIRICHLET_DATA} )
ADD_CXX_FLAGS(-DPOLORDER=${POLORDER} -DPRESSURE_POLORDER=${PRESSUREPOLORDER} -DVELOCITY_POLORDER=${VELOCITYPOLORDER} )
ADD_CXX_FLAGS(-DINNER_CG_SOLVERTYPE=${INNER_CG_SOLVERTYPE} -DOUTER_CG_SOLVERTYPE=${OUTER_CG_SOLVERTYPE} )
ADD_CXX_FLAGS( -I${GRAPE_PATH}
	-I${ALUGRID_BASE_PATH}/include
	-I${ALUGRID_BASE_PATH}/include/serial 
	-I${ALUGRID_BASE_PATH}/include/duneinterface -DENABLE_ALUGRID )

IF ( ENABLE_PARALELL )
	LINK_DIRECTORIES( ${METIS_DIR} )
	Set ( ParaLIBS "-lparmetis -lmetis" )
	ADD_CXX_FLAGS( -I${ALUGRID_BASE_PATH}/include/parallel
					-DMPIPP_H -DENABLE_MPI=1
					)
ENDIF ( ENABLE_PARALELL )

ADD_DUNE_MODULES( common istl grid fem stuff stokes )

LINK_DIRECTORIES( "/usr/lib" "${ALUGRID_BASE_PATH}/lib" "${GRAPE_PATH}" )
ADD_CXX_FLAGS( ${CXX_FLAGS} )

#----------------------------------------------------------------------------------------------------
# Source listing
#----------------------------------------------------------------------------------------------------
SET(dune-navier-stokesSrc
	src/dune_navier_stokes.cc
)

FILE( GLOB header "src/*.hh" )
FILE( GLOB_RECURSE navier "navier/*.hh" )
FILE( GLOB_RECURSE oseen "oseen/*.hh" )
FILE( GLOB_RECURSE stokes "../dune-stokes/stokes/*.hh" )
FILE( GLOB stokes_src1 "../dune-stokes/src/*.hh" )
FILE( GLOB stokes_src2 "../dune-stokes/src/*.cc" )
LIST( APPEND stokes ${stokes_src1} ${stokes_src2} )
FILE( GLOB stuff "../dune-stuff/stuff/*.hh" )


set_source_files_properties( ${DUNE_HEADERS} ${stuff} ${navier} ${oseen} ${stokes} ${header} PROPERTIES HEADER_FILE_ONLY 1 )

set( COMMON_LIBS "fem" "grid" "common" ${BLAS_LIB} "alugrid" ${UGLIBS}  "GL" ${ParaLIBS} )
set( COMMON_HEADER ${header} ${stokes} ${stuff} ${navier} ${DUNE_HEADERS} )

ADD_EXECUTABLE(testing src/testing.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(testing ${COMMON_LIBS} )

ADD_EXECUTABLE(navier ${dune-navier-stokesSrc} ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(navier ${COMMON_LIBS} )

ADD_EXECUTABLE(oseen src/oseen.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(oseen ${COMMON_LIBS} )

ADD_EXECUTABLE(conv_diff src/conv_diff.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(conv_diff ${COMMON_LIBS} )

HEADERCHECK( ${header} ${oseen} ${navier} )
