IF(NOT DD4hep_FOUND)
  RETURN()
ENDIF()

  SET( test_name "BeamCalReco" )
  ADD_TEST( NAME t_${test_name}
    COMMAND
    TestBeamCalReco $ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v14/CLIC_o3_v14.xml BeamCal
    )
  SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
    FAIL_REGULAR_EXPRESSION  "phi not correctly reconstructed"
    REQUIRED_FILES "$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v14/CLIC_o3_v14.xml"
    )

  SET( test_name "LumiCalReco" )
  ADD_TEST( NAME t_${test_name}
    COMMAND
    TestLumiCalReco
    )
  SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
    FAIL_REGULAR_EXPRESSION  "Phi not correctly reconstructed|Theta not correctly reconstructed|Wrong number of"
    )

  SET( test_name "TestLumi2Clu" )
  ADD_TEST( NAME t_${test_name}
    COMMAND
    TestLumi2Clu
    )
  SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
    FAIL_REGULAR_EXPRESSION  "Phi not correctly reconstructed|Theta not correctly reconstructed|Wrong number of"
    )

  SET( test_name "NewLumimCalReco" )
  ADD_TEST( NAME t_${test_name}
    COMMAND
    TestBeamCalReco $ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v14/CLIC_o3_v14.xml LumiCal
    )
  SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
    FAIL_REGULAR_EXPRESSION  "phi not correctly reconstructed"
    REQUIRED_FILES "$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v14/CLIC_o3_v14.xml"
    )

ADD_TEST( NAME test_Sim_BeamCal
  COMMAND ddsim --compactFile=${CMAKE_SOURCE_DIR}/Tests/CLIC_o3_v14/CLIC_o3_v14.xml --runType=run -G --outputFile=BeamCal_sim.slcio --crossingAngleBoost=0.010 --macroFile beamcal.mac
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/Tests/ CONFIGURATIONS SIM
  )

# extract libMarlinDD4hep from existing MARLIN_DLL
string(REGEX MATCH ":\?[-a-z0-9A-Z_/.]*libMarlinDD4hep.so[.0-9]*" MY_MARLIN_DLL $ENV{MARLIN_DLL})

ADD_TEST( NAME test_Reco_BeamCal
  COMMAND Marlin
  --InitDD4hep.DD4hepXMLFile=${CMAKE_SOURCE_DIR}/Tests/CLIC_o3_v14/CLIC_o3_v14.xml
  --global.LCIOInputFiles=BeamCal_sim.slcio
  --global.MaxRecordNumber=-1
  --LCIOoutput.LCIOOutputFile=BeamCal_rec.slcio
  TestReconstruction.xml
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/Tests/ CONFIGURATIONS SIM
  )

ADD_TEST( NAME test_Sim_LumiCal
  COMMAND ddsim --compactFile=${CMAKE_SOURCE_DIR}/Tests/CLIC_o3_v14/CLIC_o3_v14.xml --runType=run -G --outputFile=LumiCal_sim.slcio --crossingAngleBoost=0.010 --macroFile lumical.mac
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/Tests/ CONFIGURATIONS SIM
  )

ADD_TEST( NAME test_Reco_LumiCal
  COMMAND Marlin
  --InitDD4hep.DD4hepXMLFile=${CMAKE_SOURCE_DIR}/Tests/CLIC_o3_v14/CLIC_o3_v14.xml
  --global.LCIOInputFiles=LumiCal_sim.slcio
  --global.MaxRecordNumber=-1
  --LCIOoutput.LCIOOutputFile=LumiCal_rec.slcio
  TestReconstruction.xml
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/Tests/ CONFIGURATIONS SIM
  )

SET_TESTS_PROPERTIES(
  test_Sim_LumiCal
  test_Sim_BeamCal
  test_Reco_LumiCal
  test_Reco_BeamCal
  PROPERTIES
  FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error"
  # does not matter for Sim tests, does not harm
  ENVIRONMENT MARLIN_DLL=$<TARGET_FILE:FCalClusterer>:${MY_MARLIN_DLL}
  )
