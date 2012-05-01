PROJECT(dune_navier)

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
# Options, that can be changed be the User in order to customise dune_navier
#----------------------------------------------------------------------------------------------------
SET( BLAS_LIB
	"blas" CACHE STRING
	"" )

SET( USE_ISTL
	0 CACHE STRING
	"Use fem-wrapped istl matrix" )
		
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

SET( INNER_SOLVER
	"CG" CACHE STRING
	"INNER_CG_SOLVERTYPE" )

SET( OUTER_SOLVER
	"CG" CACHE STRING
	"OUTER_CG_SOLVERTYP" )

SET( OSEEN_DATA_NAMESPACE
	"NavierProblems::Trivial" CACHE STRING
	"OSEEN_DATA_NAMESPACE" )

SET( DIRICHLET_DATA 
	"DirichletData" CACHE STRING
	"DIRICHLET_DATA implementation for default dirichlet traits" )

SET( ALUGRID_BASE_PATH
	"/opt/dune/modules/alugrid/alugrid-1.50" CACHE STRING
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
INCLUDE_DIR( ${dune_navier_BINARY_DIR} )
INCLUDE_DIR( ${dune_navier_SOURCE_DIR} )
INCLUDE_DIR( ${dune_navier_SOURCE_DIR}/.. )

IF( ${GRIDTYPE} STREQUAL "UGGRID" )
	ADD_CXX_FLAGS( -DENABLE_UG )
	INCLUDE_SYS_DIR( ${UGGRID_BASE_PATH}/include )
	LINK_DIRECTORIES( ${UGGRID_BASE_PATH}/lib/ )
	SET( UGLIBS "ugS3" "ugS2" "devS")
ENDIF( ${GRIDTYPE} STREQUAL "UGGRID" )

IF( ${ALUGRID_FOUND} )
	INCLUDE_SYS_DIR( ${ALUGRID_BASE_PATH}/include )
	INCLUDE_SYS_DIR( ${ALUGRID_BASE_PATH}/include/serial  )
	INCLUDE_SYS_DIR( ${ALUGRID_BASE_PATH}/include/duneinterface  )
	SET( ALUGRID_LIBS "alugrid" )
ENDIF( ${ALUGRID_FOUND} )	

SET( ENABLE_ALUGRID OFF )
IF( ${GRIDTYPE} MATCHES "ALUGRID" )
#	IF( ${ALUGRID_FOUND} )
#		MESSAGE( FATAL ) 
#	ENDIF( ${ALUGRID_FOUND} )	
	INCLUDE_SYS_DIR( 	${GRAPE_PATH} 
						${ALUGRID_BASE_PATH}/include 
						${ALUGRID_BASE_PATH}/include/serial 
						${ALUGRID_BASE_PATH}/include/duneinterface )
	LINK_DIRECTORIES( "${ALUGRID_BASE_PATH}/lib"  )
	SET( ALULIBS "alugrid" )
	SET( ENABLE_ALUGRID ON )
	IF ( ENABLE_PARALELL )
		INCLUDE_SYS_DIR( ${ALUGRID_BASE_PATH}/include/parallel )
	ENDIF ( ENABLE_PARALELL )
	ADD_CXX_FLAGS( "-DENABLE_ALUGRID" )
ENDIF( ${GRIDTYPE} MATCHES "ALUGRID" )

IF ( ENABLE_PARALELL )
	LINK_DIRECTORIES( ${METIS_DIR} )
	Set ( ParaLIBS "-lparmetis -lmetis" )
	INCLUDE_SYS_DIR( ${ALUGRID_BASE_PATH}/include/parallel )
	ADD_CXX_FLAGS( -DMPIPP_H -DENABLE_MPI=1 )
ENDIF ( ENABLE_PARALELL )

SET( USE_OMP "0" )
IF( ENABLE_OMP )
	ADD_CXX_FLAGS( -fopenmp -DUSE_OMP=1)
	LIST( APPEND ParaLIBS "-fopenmp" )
ELSE( ENABLE_OMP )
	ADD_CXX_FLAGS( -DUSE_OMP=0)
ENDIF( ENABLE_OMP )

IF( USE_ISTL )
	SET( OUTER_CG_SOLVERTYPE "ISTL${OUTER_SOLVER}Op" )
	SET( INNER_CG_SOLVERTYPE "ISTL${INNER_SOLVER}Op" )
ELSE( USE_ISTL )
	SET( OUTER_CG_SOLVERTYPE "OEM${OUTER_SOLVER}Op" )
	SET( INNER_CG_SOLVERTYPE "OEM${INNER_SOLVER}Op" )
ENDIF( USE_ISTL )

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_config.h )

ADD_DUNE_MODULES( common istl grid fem geometry )
ADD_MY_MODULES( stuff oseen )
INCLUDE_SYS_DIR( ${dune_navier_SOURCE_DIR}/../dune-oseen/oseen/oemsolver/ )
ADD_CXX_FLAGS( ${CUSTOM_FLAGS} )

LINK_DIRECTORIES( "/usr/lib" "${ALUGRID_BASE_PATH}/lib" ${CCGNU_LIBRARY_DIRS} )

FILE( GLOB_RECURSE header "src/*.hh" )
FILE( GLOB_RECURSE navier "navier/*.hh" )
FILE( GLOB_RECURSE oseen "oseen/*.hh" )
FILE( GLOB_RECURSE stokes "../dune-oseen/stokes/*.hh" )
FILE( GLOB stokes_src1 "../dune-oseen/src/*.hh" )
FILE( GLOB stokes_src2 "../dune-oseen/src/*.cc" )
LIST( APPEND stokes ${stokes_src1} ${stokes_src2} )

set( COMMON_LIBS "dunefem" "dunegrid" "dunecommon" "dunegeometry" ${BLAS_LIB} ${ALUGRID_LIBS} ${UGLIBS} ${ParaLIBS} ${CCGNU_LIBRARIES} "boost_date_time" "boost_filesystem" "boost_system" )
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

SET( PROBLEMS TwoDeeTube Cockburn Null Taylor TwoDee TimeDisc Trivial Runtime DampedParallel)
SET( INCOMPLETEPROBLEMS ThreeDee )

FOREACH( problem ${PROBLEMS} )
	SET( targetName nv_${problem} )
	ADD_EXECUTABLE(${targetName} ${dune_navier_SOURCE_DIR}/src/dune_navier_stokes.cc ${COMMON_HEADER} )
	TARGET_LINK_LIBRARIES( ${targetName} ${COMMON_LIBS} )
	set_target_properties(${targetName} PROPERTIES COMPILE_FLAGS -DNAVIER_DATA_NAMESPACE=NavierProblems::${problem} )
ENDFOREACH( problem )
