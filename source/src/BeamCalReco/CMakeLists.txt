SET( BeamCalReco_SOURCES BeamCal.cpp BCUtilities.cpp BCRootUtilities.cpp BeamCalGeo.cpp BeamCalGeoGear.cpp BeamCalGeoCached.cpp  BCPadEnergies.cpp BeamCalCluster.cpp BCPCuts.cpp BCRecoObject.cpp  )

INCLUDE_DIRECTORIES ( ${CMAKE_SOURCE_DIR}/source/include/BeamCalReco )

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
    FILES_MATCHING PATTERN "*.hh")