PROJECT(dune-navier-stokes)

#set minimum cmake version
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/../dune-stuff/cmake )
include(DuneUtils)
include(UseDoxygen)

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
	"YASPGRID" CACHE STRING
	"GRIDTYPE" )

SET( GRIDDIM
	"2" CACHE STRING
	"GRIDDIM" )

SET( NAVIER_DATA_NAMESPACE
	"NavierProblems::Trivial" CACHE STRING
	"navier test data define" )

SET( POLORDER
	"1" CACHE STRING
	"Sigma POLORDER" )

SET( PRESSUREPOLORDER
	"1" CACHE STRING
	"Pressure POLORDER" )

SET( VELOCITYPOLORDER
	"2" CACHE STRING
	"Velocity POLORDER" )

SET( INNER_CG_SOLVERTYPE
	"OEMBICGSTABOp" CACHE STRING
	"INNER_CG_SOLVERTYPE" )

SET( OUTER_CG_SOLVERTYPE
	"OEMBICGSTABOp" CACHE STRING
	"OUTER_CG_SOLVERTYPE" )

SET( OSEEN_DATA_NAMESPACE
	"NavierProblems::Trivial" CACHE STRING
	"OSEEN_DATA_NAMESPACE" )

SET( DIRICHLET_DATA 
	"DirichletData" CACHE STRING
	"DIRICHLET_DATA implementation for default dirichlet traits" )

SET( ALUGRID_BASE_PATH
	"/opt/dune/modules/alugrid/alugrid-1.22" CACHE STRING
	"ALUGRID_BASE_PATH" )

SET( UGGRID_BASE_PATH
	"/share/dune/Modules/modules_x86_64/ug" CACHE STRING
	"UGGRID_BASE_PATH" )

SET( ENABLE_PARALELL
	0 CACHE BOOL
	"Enable parallel features" )

SET ( METIS_DIR
	"/share/dune/Modules/modules_x86_64/ParMetis-3.1.1" CACHE STRING
	"metis toplevel directory" )

SET_CONFIGHEADER_VARS()
CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_config.h )
#needs to come first to include correct config.h
INCLUDE_DIR( ${dune-navier-stokes_BINARY_DIR} )
INCLUDE_DIR( ${dune-navier-stokes_SOURCE_DIR} )
INCLUDE_DIR( ${dune-navier-stokes_SOURCE_DIR}/.. )

IF( ${GRIDTYPE} STREQUAL "UGGRID" )
	ADD_CXX_FLAGS( -isystem$${UGGRID_BASE_PATH}/include -DENABLE_UG )
	LINK_DIRECTORIES( ${UGGRID_BASE_PATH}/lib/ )
	SET( UGLIBS "ugS3" "ugS2" "devS")
ENDIF( ${GRIDTYPE} STREQUAL "UGGRID" )

IF( ${ALUGRID_FOUND} )
	INCLUDE_SYS_DIR( ${ALUGRID_BASE_PATH}/include )
	INCLUDE_SYS_DIR( ${ALUGRID_BASE_PATH}/include/serial  )
	INCLUDE_SYS_DIR( ${ALUGRID_BASE_PATH}/include/duneinterface  )
	SET( ALUGRID_LIBS "alugrid" )
ENDIF( ${ALUGRID_FOUND} )	

IF ( ENABLE_PARALELL )
	LINK_DIRECTORIES( ${METIS_DIR} )
	Set ( ParaLIBS "-lparmetis -lmetis" )
	INCLUDE_SYS_DIR( ${ALUGRID_BASE_PATH}/include/parallel )
	ADD_CXX_FLAGS( -DMPIPP_H -DENABLE_MPI=1 )
ENDIF ( ENABLE_PARALELL )

ADD_DUNE_MODULES( common istl grid fem )
ADD_MY_MODULES( stuff stokes )
INCLUDE_SYS_DIR( ${dune-navier-stokes_SOURCE_DIR}/../dune-stokes/stokes/oemsolver/ )
ADD_CXX_FLAGS( ${CUSTOM_FLAGS} )

LINK_DIRECTORIES( "/usr/lib" "${ALUGRID_BASE_PATH}/lib" ${CCGNU_LIBRARY_DIRS} )

FILE( GLOB_RECURSE header "src/*.hh" )
FILE( GLOB_RECURSE navier "navier/*.hh" )
FILE( GLOB_RECURSE oseen "oseen/*.hh" )
FILE( GLOB_RECURSE stokes "../dune-stokes/stokes/*.hh" )
FILE( GLOB stokes_src1 "../dune-stokes/src/*.hh" )
FILE( GLOB stokes_src2 "../dune-stokes/src/*.cc" )
LIST( APPEND stokes ${stokes_src1} ${stokes_src2} )

set( COMMON_LIBS "dunefem" "dunegrid" "dunecommon" ${BLAS_LIB} ${ALUGRID_LIBS} ${UGLIBS} ${ParaLIBS} ${CCGNU_LIBRARIES} "boost_date_time" "boost_filesystem" "boost_system" )
set( COMMON_HEADER ${header} ${stokes} ${stuff} ${navier} ${DUNE_HEADERS} )
set_source_files_properties( ${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1 )

ADD_EXECUTABLE(testing src/testing.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(testing ${COMMON_LIBS} )

ADD_EXECUTABLE(navier src/dune_navier_stokes.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(navier ${COMMON_LIBS} )

ADD_EXECUTABLE(oseen src/oseen.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(oseen ${COMMON_LIBS} )

ADD_EXECUTABLE(conv_diff src/conv_diff.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(conv_diff ${COMMON_LIBS} )

HEADERCHECK( ${header} ${oseen} ${navier} )
ADD_DEFINITIONS( ${MY_CXX_FLAGS} )
ADD_CPPCHECK( ${CMAKE_CURRENT_SOURCE_DIR}/src/dune_navier_stokes.cc ${COMMON_HEADER} )
