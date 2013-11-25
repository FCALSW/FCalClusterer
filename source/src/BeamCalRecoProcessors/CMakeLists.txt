## BeamCalReco Processors
SET( BeamCalProcessor_SOURCE ReadBeamCal.cpp BeamCalClusterReco)
INCLUDE_DIRECTORIES ( ${CMAKE_SOURCE_DIR}/source/include/BeamCalRecoProcessors )

#LumiCalReco Processors
SET( LumiCalProcessor_SOURCE MarlinLumiCalClusterer.cc)
INCLUDE_DIRECTORIES ( ${CMAKE_SOURCE_DIR}/source/include/LumiCalReco )
INCLUDE_DIRECTORIES ( ${CMAKE_SOURCE_DIR}/source/src/LumiCalReco/include )
INCLUDE_DIRECTORIES ( ${CMAKE_SOURCE_DIR}/source/src/LumiCalReco/include/LumiCalClusterer )

ADD_LIBRARY ( FCalRecoProcessors SHARED ${BeamCalProcessor_SOURCE} ${LumiCalProcessor_SOURCE} )

TARGET_LINK_LIBRARIES ( FCalRecoProcessors BeamCalReco LumiCalReco )


INSTALL( TARGETS
  FCalRecoProcessors
  LIBRARY DESTINATION lib)
