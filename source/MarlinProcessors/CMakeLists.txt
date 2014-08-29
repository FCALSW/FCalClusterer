INCLUDE_DIRECTORIES ( ./include )
INCLUDE_DIRECTORIES ( ${CMAKE_SOURCE_DIR}/source/BeamCalReco/include )
INCLUDE_DIRECTORIES ( ${CMAKE_SOURCE_DIR}/source/LumiCalReco/include )
INCLUDE_DIRECTORIES ( ${CMAKE_SOURCE_DIR}/source/ShowerTemplates/include )

## BeamCalReco Processors
SET( BeamCalProcessor_SOURCE src/ReadBeamCal.cpp src/ReadBC_SepEvt.cpp src/BeamCalClusterReco.cpp)

#LumiCalReco Processors
SET( LumiCalProcessor_SOURCE src/MarlinLumiCalClusterer.cc)


ADD_LIBRARY ( FCalRecoProcessors SHARED ${BeamCalProcessor_SOURCE} ${LumiCalProcessor_SOURCE} )

TARGET_LINK_LIBRARIES ( FCalRecoProcessors BeamCalReco LumiCalReco showerTemplates)


INSTALL( TARGETS
  FCalRecoProcessors
  LIBRARY DESTINATION lib)
