SET( Utilities_SOURCES 
  src/RootUtils.cpp
  src/BCBackgroundPar.cpp
  src/BackgroundFitter.cpp
 )

INCLUDE_DIRECTORIES ( ./include )

ADD_LIBRARY( FCalUtils SHARED ${Utilities_SOURCES} )
TARGET_LINK_LIBRARIES( FCalUtils ${ROOT_LIBRARIES}  )

ADD_EXECUTABLE( BCBackgroundPar src/BCBackgroundPar.cpp src/BackgroundFitter.cpp )

TARGET_LINK_LIBRARIES( BCBackgroundPar ${ROOT_LIBRARIES}  )

INSTALL( TARGETS
  FCalUtils
  LIBRARY DESTINATION lib)

INSTALL( TARGETS
  BCBackgroundPar
  DESTINATION bin)
