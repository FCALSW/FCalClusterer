IF( DD4hep_FOUND )

  SET( test_name "BeamCalReco" )
  ADD_TEST( NAME t_${test_name}
    COMMAND
    ${CMAKE_SOURCE_DIR}/bin/TestBeamCalReco $ENV{lcgeo_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml
    )
  SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
    FAIL_REGULAR_EXPRESSION  "phi not correctly reconstructed"
    REQUIRED_FILES "$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml"
    )

ENDIF()
