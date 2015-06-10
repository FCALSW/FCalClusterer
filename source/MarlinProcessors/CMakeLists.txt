INCLUDE_DIRECTORIES ( ./include )
INCLUDE_DIRECTORIES ( ${CMAKE_SOURCE_DIR}/source/BeamCalReco/include )
INCLUDE_DIRECTORIES ( ${CMAKE_SOURCE_DIR}/source/LumiCalReco/include )

## BeamCalReco Processors
SET( BeamCalProcessor_SOURCE src/ReadBeamCal.cpp src/BeamCalClusterReco.cpp)

#LumiCalReco Processors
SET( LumiCalProcessor_SOURCE src/MarlinLumiCalClusterer.cc)


ADD_LIBRARY ( FCalClusterer SHARED ${BeamCalProcessor_SOURCE} ${LumiCalProcessor_SOURCE} )

TARGET_LINK_LIBRARIES ( FCalClusterer BeamCalReco LumiCalReco )


INSTALL( TARGETS
  FCalClusterer
  LIBRARY DESTINATION lib)
