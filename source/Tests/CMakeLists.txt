INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/source/BeamCalReco/include)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/source/LumiCalReco/include)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/source/Tests/include)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/source/Utilities/include)

IF(DD4hep_FOUND)
  ADD_EXECUTABLE(TestBeamCalReco src/TestBeamCalReco.cpp)
  TARGET_LINK_LIBRARIES(TestBeamCalReco BeamCalReco LumiCalReco)

  ADD_EXECUTABLE(TestLumiCalReco src/TestLumiCalReco.cpp)
  TARGET_LINK_LIBRARIES(TestLumiCalReco LumiCalReco)

  ADD_EXECUTABLE(TestLumi2Clu src/TestLumi2Clu.cpp)
  TARGET_LINK_LIBRARIES(TestLumi2Clu LumiCalReco)

  INSTALL(TARGETS TestLumiCalReco TestBeamCalReco TestLumi2Clu
    RUNTIME DESTINATION bin)
ENDIF()
