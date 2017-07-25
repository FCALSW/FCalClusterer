#ADD_LIBRARY( BeamCal SHARED BeamCal.cpp BCUtilities.cpp BeamCalGeoGear.cpp )
#ADD_LIBRARY( GearBeamCal SHARED GearBeamCal.cpp )
#ADD_LIBRARY( BCPadEnergies SHARED BCPadEnergies.cpp BeamCalCluster.cpp BCPCuts.cpp BCRecoObject.cpp )
#ADD_LIBRARY( BeamCalReco SHARED BeamCal GearBeamCal BCPadEnergies )


#TARGET_LINK_LIBRARIES( GearBeamCal ${ROOT_LIBRARIES} ${GEAR_LIBRARIES} ${DD4hep_LIBRARIES} )
#TARGET_LINK_LIBRARIES( BeamCal BCPadEnergies ${ROOT_LIBRARIES} )

INCLUDE_DIRECTORIES ( ${CMAKE_SOURCE_DIR}/source/LumiCalReco/include )
INCLUDE_DIRECTORIES ( ${CMAKE_SOURCE_DIR}/source/BeamCalReco/include )

ADD_EXECUTABLE (ReconstructBecas ReconstructBecas.cpp)
TARGET_LINK_LIBRARIES ( ReconstructBecas BeamCalReco )

ADD_EXECUTABLE ( DrawBeamCals DrawBeamCals.cpp)
TARGET_LINK_LIBRARIES ( DrawBeamCals BeamCalReco )

IF( DD4hep_FOUND )
  ADD_EXECUTABLE (TestBeamCalReco TestBeamCalReco.cpp)
  TARGET_LINK_LIBRARIES ( TestBeamCalReco BeamCalReco )
  INSTALL( TARGETS
    TestBeamCalReco
    RUNTIME DESTINATION bin)
ENDIF()
