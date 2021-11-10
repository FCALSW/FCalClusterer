SET( BeamCalReco_SOURCES 
  src/BeamCal.cpp
  src/BCUtilities.cpp
  src/BCRootUtilities.cpp
  src/BeamCalGeo.cpp
  src/BeamCalGeoGear.cpp
  src/BeamCalGeoCached.cpp
  src/BeamCalBkg.cpp
  src/BeamCalBkgPregen.cpp
  src/BeamCalBkgParam.cpp
  src/BeamCalBkgGauss.cpp
  src/BeamCalBkgAverage.cpp
  src/BeamCalBkgEmpty.cpp
  src/BeamCalBkgFactory.cpp
  src/BeamCalFitShower.cpp
  src/BeamCalPadGeometry.cpp
  src/BCPadEnergies.cpp
  src/BeamCalCluster.cpp
  src/BCPCuts.cpp
  src/BCRecoObject.cpp
 )

IF( DD4hep_FOUND )
  LIST(APPEND BeamCalReco_SOURCES src/BeamCalGeoDD.cpp )
ENDIF()

ADD_LIBRARY(BeamCalReco SHARED ${BeamCalReco_SOURCES})
TARGET_INCLUDE_DIRECTORIES(BeamCalReco PUBLIC include)

FIND_PACKAGE(ROOT REQUIRED COMPONENTS Minuit2 Unuran MathMore GenVector)
TARGET_INCLUDE_DIRECTORIES(BeamCalReco SYSTEM PUBLIC ${ROOT_INCLUDE_DIRS} )
TARGET_LINK_LIBRARIES(BeamCalReco PUBLIC FCalUtils)
TARGET_LINK_LIBRARIES(BeamCalReco PRIVATE ${ROOT_LIBRARIES} )
FIND_PACKAGE(ROOT REQUIRED) # reset ROOT_LIBRARIES

IF(DD4hep_FOUND)
  TARGET_INCLUDE_DIRECTORIES(BeamCalReco SYSTEM PUBLIC ${DD4hep_INCLUDE_DIRS})
  TARGET_COMPILE_DEFINITIONS(BeamCalReco PUBLIC ${DD4hep_DEFINITIONS} FCAL_WITH_DD4HEP)
  TARGET_LINK_LIBRARIES(BeamCalReco PUBLIC ${DD4hep_LIBRARIES} ${DD4hep_COMPONENT_LIBRARIES} )
ENDIF()

FOREACH( pkg GEAR streamlog LCIO)
  IF(FCAL_USE_${pkg})
    TARGET_INCLUDE_DIRECTORIES(BeamCalReco SYSTEM BEFORE PUBLIC ${${pkg}_INCLUDE_DIRS} )
    TARGET_COMPILE_DEFINITIONS(BeamCalReco PUBLIC ${${pkg}_DEFINITIONS} )
    TARGET_LINK_LIBRARIES(BeamCalReco PUBLIC ${${pkg}_LIBRARIES})
  ENDIF()
ENDFOREACH()

IF( CPPCHECK_FOUND )
  MESSAGE( STATUS "Will cppcheck beamCalReco" )
  ADD_CPPCHECK( BeamCalReco UNUSED_FUNCTIONS STYLE POSSIBLE_ERRORS FORCE FAIL_ON_WARNINGS )
ENDIF()

INSTALL( TARGETS
  BeamCalReco
  LIBRARY DESTINATION lib)

INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/source/include/BeamCalReco
    DESTINATION include
    FILES_MATCHING PATTERN "*.hh"
    PATTERN ".svn" EXCLUDE )
