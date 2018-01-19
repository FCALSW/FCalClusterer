IF( DD4hep_FOUND )

  SET( test_name "BeamCalReco" )
  ADD_TEST( NAME t_${test_name}
    COMMAND
    ${CMAKE_SOURCE_DIR}/bin/TestBeamCalReco $ENV{lcgeo_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml BeamCal
    )
  SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
    FAIL_REGULAR_EXPRESSION  "phi not correctly reconstructed"
    REQUIRED_FILES "$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml"
    )

  SET( test_name "LumiCalReco" )
  ADD_TEST( NAME t_${test_name}
    COMMAND
    ${CMAKE_SOURCE_DIR}/bin/TestLumiCalReco
    )
  SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
    FAIL_REGULAR_EXPRESSION  "Phi not correctly reconstructed|Theta not correctly reconstructed|Wrong number of"
    )



  SET( test_name "TestLumi2Clu" )
  ADD_TEST( NAME t_${test_name}
    COMMAND
    ${CMAKE_SOURCE_DIR}/bin/TestLumi2Clu
    )
  SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
    FAIL_REGULAR_EXPRESSION  "Phi not correctly reconstructed|Theta not correctly reconstructed|Wrong number of"
    )

  SET( test_name "NewLumimCalReco" )
  ADD_TEST( NAME t_${test_name}
    COMMAND
    ${CMAKE_SOURCE_DIR}/bin/TestBeamCalReco $ENV{lcgeo_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml LumiCal
    )
  SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
    FAIL_REGULAR_EXPRESSION  "phi not correctly reconstructed"
    REQUIRED_FILES "$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml"
    )

SET( test_name "test_Sim_BeamCal" )
ADD_TEST( NAME t_${test_name}
  COMMAND "$ENV{lcgeo_DIR}/bin/run_test_lcgeo.sh"
  ddsim --compactFile=${CMAKE_SOURCE_DIR}/Tests/CLIC_o3_v14/CLIC_o3_v14.xml --runType=run -G --outputFile=BeamCal_sim.slcio --crossingAngleBoost=0.010 --macroFile beamcal.mac
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/Tests/ CONFIGURATIONS SIM
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_Reco_BeamCal" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_SOURCE_DIR}/Tests/bin/MarlinTests.sh
  ${CMAKE_SOURCE_DIR}/lib/libFCalClusterer.so  Marlin
  --InitDD4hep.DD4hepXMLFile=${CMAKE_SOURCE_DIR}/Tests/CLIC_o3_v14/CLIC_o3_v14.xml
  --global.LCIOInputFiles=BeamCal_sim.slcio
  --global.MaxRecordNumber=-1
  --LCIOoutput.LCIOOutputFile=BeamCal_rec.slcio
  TestReconstruction.xml
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/Tests/ CONFIGURATIONS SIM
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )

SET( test_name "test_Sim_LumiCal" )
ADD_TEST( NAME t_${test_name}
  COMMAND "$ENV{lcgeo_DIR}/bin/run_test_lcgeo.sh"
  ddsim --compactFile=${CMAKE_SOURCE_DIR}/Tests/CLIC_o3_v14/CLIC_o3_v14.xml --runType=run -G --outputFile=LumiCal_sim.slcio --crossingAngleBoost=0.010 --macroFile lumical.mac
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/Tests/ CONFIGURATIONS SIM
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_Reco_LumiCal" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_SOURCE_DIR}/Tests/bin/MarlinTests.sh
  ${CMAKE_SOURCE_DIR}/lib/libFCalClusterer.so  Marlin
  --InitDD4hep.DD4hepXMLFile=${CMAKE_SOURCE_DIR}/Tests/CLIC_o3_v14/CLIC_o3_v14.xml
  --global.LCIOInputFiles=LumiCal_sim.slcio
  --global.MaxRecordNumber=-1
  --LCIOoutput.LCIOOutputFile=LumiCal_rec.slcio
  TestReconstruction.xml
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/Tests/ CONFIGURATIONS SIM
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )


ENDIF()
