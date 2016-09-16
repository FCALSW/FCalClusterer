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
  src/BeamCalFitShower.cpp
  src/BeamCalPadGeometry.cpp
  src/BCPadEnergies.cpp
  src/BeamCalCluster.cpp
  src/BCPCuts.cpp
  src/BCRecoObject.cpp
 )

INCLUDE_DIRECTORIES ( ./include )
INCLUDE_DIRECTORIES ( ../Utilities/include )

IF( DD4hep_FOUND )
  LIST(APPEND BeamCalReco_SOURCES src/BeamCalGeoDD.cpp )
ENDIF()

ADD_LIBRARY( BeamCalReco SHARED ${BeamCalReco_SOURCES} )

IF( CPPCHECK_FOUND )
  MESSAGE( STATUS "Will cppcheck beamCalReco" )
  ADD_CPPCHECK( BeamCalReco UNUSED_FUNCTIONS STYLE POSSIBLE_ERRORS FORCE FAIL_ON_WARNINGS )
ENDIF()

TARGET_LINK_LIBRARIES( BeamCalReco FCalUtils ${ROOT_LIBRARIES}
${ROOT_COMPONENT_LIBRARIES} ${LCIO_LIBRARIES} ${GEAR_LIBRARIES} ${Marlin_LIBRARIES} )

INSTALL( TARGETS
  BeamCalReco
  LIBRARY DESTINATION lib)

INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/source/include/BeamCalReco
    DESTINATION include
    FILES_MATCHING PATTERN "*.hh"
    PATTERN ".svn" EXCLUDE )
