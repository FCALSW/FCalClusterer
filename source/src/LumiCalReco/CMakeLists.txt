SET( LumiCalReco_SOURCES 
  ClusterClass.cpp
  GlobalMethodsClass.cpp
  LCCluster.cpp
  LumiCalClusterer.cpp
  LumiCalClusterer_auxiliary.cpp
  LumiCalClusterer_buildClusters.cpp
  LumiCalClusterer_buildClusters_auxiliary.cpp
  LumiCalClusterer_clusterMerger.cpp
  LumiCalClusterer_energyCorrections.cpp
  LumiCalClusterer_fiducialVolumeCuts.cpp
  LumiCalClusterer_getCalHits.cpp
  SuperTrueClusterWeights.cpp
  VirtualCluster.cpp
)

INCLUDE_DIRECTORIES ( ${CMAKE_SOURCE_DIR}/source/include/LumiCalReco )
INCLUDE_DIRECTORIES ( ./include )
INCLUDE_DIRECTORIES ( ./include/LumiCalClusterer )

ADD_LIBRARY( LumiCalReco SHARED ${LumiCalReco_SOURCES} )
TARGET_LINK_LIBRARIES( LumiCalReco ${ROOT_LIBRARIES} ${LCIO_LIBRARIES} ${GEAR_LIBRARIES} ${Marlin_LIBRARIES} )

INSTALL( TARGETS
  LumiCalReco
  LIBRARY DESTINATION lib)

INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/source/include/LumiCalReco
    DESTINATION include
    FILES_MATCHING PATTERN "*.hh")