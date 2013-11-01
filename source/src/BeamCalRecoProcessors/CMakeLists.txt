ADD_LIBRARY ( BeamCalRecoProcessors SHARED ReadBeamCal.cpp )

INCLUDE_DIRECTORIES ( ${CMAKE_SOURCE_DIR}/source/include/BeamCalRecoProcessors )

TARGET_LINK_LIBRARIES ( BeamCalRecoProcessors BeamCalReco )

INSTALL( TARGETS
  BeamCalRecoProcessors
  LIBRARY DESTINATION lib)
