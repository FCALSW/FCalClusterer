#ADD_LIBRARY( BeamCal SHARED BeamCal.cpp BCUtilities.cpp BeamCalGeoGear.cpp )
#ADD_LIBRARY( GearBeamCal SHARED GearBeamCal.cpp )
#ADD_LIBRARY( BCPadEnergies SHARED BCPadEnergies.cpp BeamCalCluster.cpp BCPCuts.cpp BCRecoObject.cpp )
#ADD_LIBRARY( BeamCalReco SHARED BeamCal GearBeamCal BCPadEnergies )

ADD_SUBDIRECTORY(BeamCalReco)

#TARGET_LINK_LIBRARIES( GearBeamCal ${ROOT_LIBRARIES} ${GEAR_LIBRARIES} ${DD4hep_LIBRARIES} )
#TARGET_LINK_LIBRARIES( BeamCal BCPadEnergies ${ROOT_LIBRARIES} )

INCLUDE_DIRECTORIES ( ${CMAKE_SOURCE_DIR}/source/include )

ADD_EXECUTABLE (ReconstructBecas ReconstructBecas.cpp)
TARGET_LINK_LIBRARIES ( ReconstructBecas BeamCalReco )
