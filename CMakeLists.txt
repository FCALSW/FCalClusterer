CMAKE_MINIMUM_REQUIRED( VERSION 3.3 )

PROJECT( FCalClusterer CXX C )

SET( PROJECT_VERSION_MAJOR 0 )
SET( PROJECT_VERSION_MINOR 6 )
SET( CMAKE_BUILD_TYPE RelWithDebInfo )

SET( FCAL_USE_ROOT True ) ##not optional
SET( FCAL_USE_LCIO True ) ##not optional
SET( FCAL_USE_GEAR True ) ##not optional
SET( FCAL_USE_streamlog True ) ##not optional
OPTION( FCAL_USE_Marlin  " Build Marlin Processors" True )

SET( CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR} )

FIND_PACKAGE( ILCUTIL 1.3.0 REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )
# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )

LIST( APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" )

# FIND_PACKAGE( Cppcheck  )
# IF ( CPPCHECK_FOUND )
#   INCLUDE( CppcheckTargets )
#   MESSAGE (STATUS "Cppcheck ${CPPCHECK_FOUND}" )
# ENDIF()

FOREACH( pkg LCIO GEAR streamlog Marlin ) #DD4hep
  IF ( ${FCAL_USE_${pkg}} )
    MESSAGE ( STATUS "Building with ${pkg}" )
    FIND_PACKAGE (${pkg} REQUIRED)
    INCLUDE_DIRECTORIES( SYSTEM ${${pkg}_INCLUDE_DIRS} )
    ADD_DEFINITIONS( ${${pkg}_DEFINITIONS} )
  ENDIF()
ENDFOREACH()



FIND_PACKAGE( DD4hep COMPONENTS DDRec )
IF( DD4hep_FOUND )
  INCLUDE_DIRECTORIES( SYSTEM ${DD4hep_INCLUDE_DIRS} )
  ADD_DEFINITIONS( ${DD4hep_DEFINITIONS} )
  ADD_DEFINITIONS( "-DFCAL_WITH_DD4HEP" )
ENDIF()


IF ( ${FCAL_USE_ROOT} )
  MESSAGE ( STATUS "Building with ROOT" )
  FIND_PACKAGE ( ROOT REQUIRED COMPONENTS Minuit2 Unuran MathMore GenVector)
  INCLUDE_DIRECTORIES( SYSTEM ${ROOT_INCLUDE_DIRS} )
  ADD_DEFINITIONS( ${ROOT_DEFINITIONS} )
ENDIF()

FIND_PROGRAM(CLANG_TIDY "clang-tidy")
IF(CLANG_TIDY)
  MESSAGE(STATUS "FOUND clang-tidy: ${CLANG_TIDY}")
  SET(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY})
ENDIF()

ADD_SUBDIRECTORY(source/)

MESSAGE (STATUS "CMAKE_FLAGS ${CMAKE_CXX_FLAGS}" )
MESSAGE (STATUS "CMAKE_FLAGS DEBUG ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" )

option(BUILD_TESTING "Enable and build tests" ON)
if(BUILD_TESTING)
  include(CTest)
  enable_testing()
  add_subdirectory(Tests)
endif(BUILD_TESTING)

### FORMATTING ##############################################################
# Set the source files to clang-format (FIXME: determine this better
FILE(GLOB_RECURSE
     CHECK_CXX_SOURCE_FILES
     /source/ *.cc *.cpp *.h *.hh
     )
INCLUDE("cmake/clang-cpp-checks.cmake")
