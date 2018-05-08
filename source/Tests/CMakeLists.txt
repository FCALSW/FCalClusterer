IF(DD4hep_FOUND)

  INCLUDE_DIRECTORIES(./include)

  ADD_EXECUTABLE(TestBeamCalReco src/TestBeamCalReco.cpp)
  TARGET_LINK_LIBRARIES(TestBeamCalReco BeamCalReco LumiCalReco)

  ADD_EXECUTABLE(TestLumiCalReco src/TestLumiCalReco.cpp)
  TARGET_LINK_LIBRARIES(TestLumiCalReco LumiCalReco BeamCalReco)

  ADD_EXECUTABLE(TestLumi2Clu src/TestLumi2Clu.cpp)
  TARGET_LINK_LIBRARIES(TestLumi2Clu LumiCalReco BeamCalReco)

  INSTALL(TARGETS TestLumiCalReco TestBeamCalReco TestLumi2Clu
    RUNTIME DESTINATION bin)
ENDIF()
