SET( LumiCalReco_SOURCES 
  src/ClusterClass.cpp
  src/MCInfo.cpp
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

ADD_LIBRARY( LumiCalReco SHARED ${LumiCalReco_SOURCES} )

TARGET_INCLUDE_DIRECTORIES(LumiCalReco PUBLIC include)
TARGET_LINK_LIBRARIES(LumiCalReco ${ROOT_LIBRARIES} ${LCIO_LIBRARIES} ${GEAR_LIBRARIES} ${Marlin_LIBRARIES})

TARGET_INCLUDE_DIRECTORIES(LumiCalReco SYSTEM PUBLIC ${ROOT_INCLUDE_DIR} ${ROOT_INCLUDE_DIRS} )
TARGET_COMPILE_DEFINITIONS(LumiCalReco PUBLIC ${ROOT_DEFINITIONS} )

FOREACH( pkg LCIO GEAR streamlog Marlin)
  IF ( ${FCAL_USE_${pkg}} )
    TARGET_INCLUDE_DIRECTORIES(LumiCalReco SYSTEM PUBLIC ${${pkg}_INCLUDE_DIRS} )
    TARGET_COMPILE_DEFINITIONS(LumiCalReco PUBLIC ${${pkg}_DEFINITIONS} )
  ENDIF()
ENDFOREACH()

IF( DD4hep_FOUND )
  TARGET_INCLUDE_DIRECTORIES(LumiCalReco SYSTEM PUBLIC ${DD4hep_INCLUDE_DIRS} )
  TARGET_COMPILE_DEFINITIONS(LumiCalReco PUBLIC ${DD4hep_DEFINITIONS} FCAL_WITH_DD4HEP )
  TARGET_LINK_LIBRARIES(LumiCalReco ${DD4hep_LIBRARIES} ${DD4hep_COMPONENT_LIBRARIES})
ENDIF()


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