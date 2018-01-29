INCLUDE_DIRECTORIES ( ${CMAKE_SOURCE_DIR}/source/BeamCalReco/include )

ADD_EXECUTABLE (ReconstructBecas ReconstructBecas.cpp)
TARGET_LINK_LIBRARIES ( ReconstructBecas BeamCalReco )

ADD_EXECUTABLE ( DrawBeamCals DrawBeamCals.cpp)
TARGET_LINK_LIBRARIES ( DrawBeamCals BeamCalReco )

