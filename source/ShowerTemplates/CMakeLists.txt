SET( ShowerTemplates_SOURCES 
  src/showerTemplates.cpp
 )

INCLUDE_DIRECTORIES ( ./include )

ADD_LIBRARY( showerTemplates SHARED ${ShowerTemplates_SOURCES} )
TARGET_LINK_LIBRARIES( showerTemplates ${ROOT_LIBRARIES} )

INSTALL( TARGETS
  showerTemplates
  LIBRARY DESTINATION lib)

INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/source/include/showerTemplates
    DESTINATION include
    FILES_MATCHING PATTERN "*.hh"
    PATTERN ".svn" EXCLUDE )