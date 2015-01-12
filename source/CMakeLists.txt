INCLUDE_DIRECTORIES ( ${CMAKE_SOURCE_DIR}/source/include )

ADD_SUBDIRECTORY ( Utilities )
ADD_SUBDIRECTORY ( BeamCalReco )
ADD_SUBDIRECTORY ( LumiCalReco )

if(Marlin_FOUND)
  ADD_SUBDIRECTORY ( MarlinProcessors )
endif(Marlin_FOUND)


ADD_SUBDIRECTORY ( src/ )

