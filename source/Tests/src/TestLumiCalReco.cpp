#include "TestMain.hh"

#include "GlobalMethodsClass.h"
#include "LCCluster.hh"
#include "LumiCalClusterer.h"

#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>

#include <cmath>
#include <iostream>
#include <vector>

int runTest(int, char**) {
  bool               failed = false;
  GlobalMethodsClass gmc;
  initializeGMC(gmc);
  LumiCalClustererClass lcc("lumiCollection");
  try {
    lcc.init(gmc);
  } catch (std::out_of_range& e) {
    std::cerr << "Failed to setup LCC: " << e.what() << std::endl;
    return 1;
  }

  std::vector<int>    testCases    = {0, 8, 16, 24, 32, 40, 48, 56};  //pads for eight directions
  std::vector<double> expectedPhis = {
      0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0,  //backward
      0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0,  //forward
  };

  int counter = 0;
  for (int directionID = -1; directionID <= 1; directionID += 2) {
    // //double expectedTheta = (100.00 + 200.0/40*11.5)/3500.0;
    // if(directionID < 0) expectedTheta = M_PI - expectedTheta;
    for (auto padIDToFill : testCases) {
      streamlog_out(MESSAGE) << "**********************************************************************" << std::endl;
      auto              expectedPhi   = expectedPhis[counter++];
      double            expectedTheta = 0.0;  //filled later
      IMPL::LCEventImpl myEvt;
      auto*             col = new IMPL::LCCollectionVec(LCIO::CALORIMETERHIT);
      col->parameters().setValue(LCIO::CellIDEncoding, "system:8,barrel:3,layer:8,slice:8,r:32:-16,phi:-16");
      fillLumiCal(gmc, col, padIDToFill, directionID == -1 ? 2 : directionID, expectedTheta);
      myEvt.addCollection(col, "lumiCollection");

      lcc.processEvent(&myEvt);

      for (int armNow = -1; armNow < 2; armNow += 2) {
        streamlog_out(MESSAGE) << " Arm  " << std::setw(4) << armNow
                               << "\t Number of clusters: " << lcc._superClusterIdToCellId[armNow].size() << std::endl;
      }

      if (lcc._superClusterIdToCellId[directionID].size() != 1) {
        streamlog_out(ERROR) << "ERROR: Wrong number of reconstructed clusters" << std::endl;
        failed = true;
      }
      for (auto const& pairIDCells : lcc._superClusterIdToCellId[directionID]) {
        const int  clusterId       = pairIDCells.first;
        LCCluster& thisClusterInfo = lcc._superClusterIdClusterInfo[directionID][clusterId];
        thisClusterInfo.recalculatePositionFromHits(gmc);
        streamlog_out(MESSAGE) << thisClusterInfo << std::endl;
        failed = not checkReconstructedObject(gmc, thisClusterInfo, expectedPhi, expectedTheta);
      }
      streamlog_out(MESSAGE) << "**********************************************************************" << std::endl;
    }
  }  // forwardAndBackward

  return failed ? 1 : 0;
}
