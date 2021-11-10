## BeamCalReco Processors
SET( BeamCalProcessor_SOURCE src/ReadBeamCal.cpp src/BeamCalClusterReco.cpp)

#LumiCalReco Processors
SET( LumiCalProcessor_SOURCE src/MarlinLumiCalClusterer.cc)


IF( DD4hep_FOUND )
  SET( DD4hepProcessor_SOURCE src/DrawBeamCalFromDD4hep.cpp )
ENDIF()

ADD_LIBRARY ( FCalClusterer SHARED ${BeamCalProcessor_SOURCE} ${LumiCalProcessor_SOURCE} ${DD4hepProcessor_SOURCE} )
TARGET_INCLUDE_DIRECTORIES(FCalClusterer PUBLIC ./include)

#LibrariesToLink
SET( FCAL_LIBRARIES LumiCalReco BeamCalReco )
IF( DD4hep_FOUND )
  SET( FCAL_LIBRARIES ${FCAL_LIBRARIES} ${DD4hep_LIBRARIES} )
ENDIF()

TARGET_LINK_LIBRARIES ( FCalClusterer ${FCAL_LIBRARIES} )
TARGET_INCLUDE_DIRECTORIES(FCalClusterer SYSTEM BEFORE PRIVATE ${LCIO_INCLUDE_DIRS})

IF( CPPCHECK_FOUND )
  ADD_CPPCHECK( FCalClusterer UNUSED_FUNCTIONS STYLE POSSIBLE_ERRORS FORCE FAIL_ON_WARNINGS )
ENDIF()


INSTALL( TARGETS
  FCalClusterer
  LIBRARY DESTINATION lib)
