SET( BeamCalReco_SOURCES 
  src/BeamCal.cpp
  src/BCUtilities.cpp
  src/BCRootUtilities.cpp
  src/BeamCalGeo.cpp
  src/BeamCalGeoGear.cpp
  src/BeamCalGeoCached.cpp
  src/BCPadEnergies.cpp
  src/BeamCalCluster.cpp
  src/BCPCuts.cpp
  src/BCRecoObject.cpp
 )

INCLUDE_DIRECTORIES ( ./include )

IF( DD4hep_FOUND )
  LIST(APPEND BeamCalReco_SOURCES BeamCalGeoDD.cpp )
ENDIF()

ADD_LIBRARY( BeamCalReco SHARED ${BeamCalReco_SOURCES} )
TARGET_LINK_LIBRARIES( BeamCalReco ${ROOT_LIBRARIES} ${LCIO_LIBRARIES} ${GEAR_LIBRARIES} ${Marlin_LIBRARIES} )

INSTALL( TARGETS
  BeamCalReco
  LIBRARY DESTINATION lib)

INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/source/include/BeamCalReco
    DESTINATION include
    FILES_MATCHING PATTERN "*.hh"
    PATTERN ".svn" EXCLUDE )