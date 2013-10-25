###############################################
# cmake configuration file for MYSTUFF
# @author Andre Sailer, CERN
###############################################

SET( BeamCalReco_FOUND FALSE )
MARK_AS_ADVANCED( BeamCalReco_FOUND )

# do not store find results in cache
SET( BeamCalReco_INCLUDE_DIR BeamCalReco_INCLUDE_DIR-NOTFOUND )

FIND_PATH( BeamCalReco_INCLUDE_DIR
	NAMES BeamCalReco/BeamCal.hh
	PATHS ${BeamCalReco_DIR}
	PATH_SUFFIXES include
	NO_DEFAULT_PATH
)
IF( NOT BeamCalReco_INCLUDE_DIR )
    MESSAGE( STATUS "Check for BeamCalReco: ${BeamCalReco_DIR}"
					" -- failed to find BeamCalReco include directory!!" )
ELSE( NOT BeamCalReco_INCLUDE_DIR )
    MARK_AS_ADVANCED( BeamCalReco_INCLUDE_DIR )
ENDIF( NOT BeamCalReco_INCLUDE_DIR )


# do not store find results in cache
SET( BeamCalReco_LIB BeamCalReco_LIB-NOTFOUND )

FIND_LIBRARY( BeamCalReco_LIB
	NAMES BeamCalReco
	PATHS ${BeamCalReco_DIR}
	PATH_SUFFIXES lib
	NO_DEFAULT_PATH
)

IF( NOT BeamCalReco_LIB )
    MESSAGE( STATUS "Check for BeamCalReco: ${BeamCalReco_DIR}"
					" -- failed to find BeamCalReco library!!" )
ENDIF( NOT BeamCalReco_LIB )


# set variables and display results
IF( BeamCalReco_INCLUDE_DIR AND BeamCalReco_LIB )
    SET( BeamCalReco_FOUND TRUE )
    SET( BeamCalReco_INCLUDE_DIRS ${BeamCalReco_INCLUDE_DIR} )
	SET( BeamCalReco_LIBRARIES ${BeamCalReco_LIB} )
	MARK_AS_ADVANCED( BeamCalReco_INCLUDE_DIRS BeamCalReco_LIBRARIES )
	MESSAGE( STATUS "Check for BeamCalReco: ${BeamCalReco_DIR} -- works" )
ELSE( BeamCalReco_INCLUDE_DIR AND BeamCalReco_LIB )
	IF( BeamCalReco_FIND_REQUIRED )
		MESSAGE( FATAL_ERROR "Check for BeamCalReco: ${BeamCalReco_DIR} -- failed!!" )
    ELSE( BeamCalReco_FIND_REQUIRED )
        MESSAGE( STATUS "Check for BeamCalReco: ${BeamCalReco_DIR}"
						" -- failed!! will skip this package..." )
    ENDIF( BeamCalReco_FIND_REQUIRED )
ENDIF( BeamCalReco_INCLUDE_DIR AND BeamCalReco_LIB )
