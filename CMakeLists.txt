CMAKE_MINIMUM_REQUIRED( VERSION 3.3 )

PROJECT( FCalClusterer CXX C )

SET( PROJECT_VERSION_MAJOR 1 )
SET( PROJECT_VERSION_MINOR 0 )
SET( CMAKE_BUILD_TYPE RelWithDebInfo )

SET( FCAL_USE_ROOT True ) ##not optional
SET( FCAL_USE_LCIO True ) ##not optional
SET( FCAL_USE_GEAR True ) ##not optional
SET( FCAL_USE_streamlog True ) ##not optional
OPTION( FCAL_USE_Marlin  " Build Marlin Processors" True )
OPTION( FCAL_USE_DD4hep  " Build With DD4hep support" True )
OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

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

FOREACH(pkg LCIO GEAR streamlog Marlin ROOT DD4hep)
  IF ( ${FCAL_USE_${pkg}} )
    MESSAGE ( STATUS "Building with ${pkg}" )
    FIND_PACKAGE (${pkg} REQUIRED)
  ENDIF()
ENDFOREACH()

IF(${FCAL_USE_DD4hep})
  FIND_PACKAGE(DD4hep COMPONENTS DDRec)
ENDIF()

FIND_PROGRAM(CLANG_TIDY "clang-tidy")
IF(CLANG_TIDY)
  MESSAGE(STATUS "FOUND clang-tidy: ${CLANG_TIDY}")
  SET(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY})
ENDIF()

ADD_SUBDIRECTORY(source/)
IF(INSTALL_DOC)
  ADD_SUBDIRECTORY(doc/)
ENDIF()

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
