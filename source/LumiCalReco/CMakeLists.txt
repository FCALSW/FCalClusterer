SET( LumiCalReco_SOURCES 
  src/ClusterClass.cpp
  src/GlobalMethodsClass.cpp
  src/LCCluster.cpp
  src/ProjectionInfo.cpp
  src/LumiCalClusterer.cpp
  src/LumiCalClusterer_auxiliary.cpp
  src/LumiCalClusterer_buildClusters.cpp
  src/LumiCalClusterer_buildClusters_auxiliary.cpp
  src/LumiCalClusterer_clusterMerger.cpp
  src/LumiCalClusterer_energyCorrections.cpp
  src/LumiCalClusterer_fiducialVolumeCuts.cpp
  src/LumiCalClusterer_getCalHits.cpp
  src/SuperTrueClusterWeights.cpp
  src/VirtualCluster.cpp
  src/MarlinLumiCalClusterer.cpp
  src/OutputManagerClass.cpp
)

INCLUDE_DIRECTORIES ( ./include )


ADD_LIBRARY( LumiCalReco SHARED ${LumiCalReco_SOURCES} )
TARGET_LINK_LIBRARIES( LumiCalReco 
  ${ROOT_LIBRARIES} ${LCIO_LIBRARIES} ${GEAR_LIBRARIES} ${Marlin_LIBRARIES}
  ${DD4hep_LIBRARIES} ${DD4hep_COMPONENT_LIBRARIES}
)

IF( CPPCHECK_FOUND )
  ADD_CPPCHECK( LumiCalReco UNUSED_FUNCTIONS STYLE POSSIBLE_ERRORS FORCE FAIL_ON_WARNINGS )
ENDIF()

INSTALL( TARGETS
  LumiCalReco
  LIBRARY DESTINATION lib)

INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/source/include/LumiCalReco
    DESTINATION include
    FILES_MATCHING PATTERN "*.hh"
    PATTERN ".svn" EXCLUDE )