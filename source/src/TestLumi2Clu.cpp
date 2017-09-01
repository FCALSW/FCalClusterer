#include "GlobalMethodsClass.h"
#include "LCCluster.hh"
#include "LumiCalClusterer.h"
#include "TestUtilities.hh"

#include <streamlog/streamlog.h>

#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

int testLumiCal() {
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

  const int                 deltaPad(5);
  const double              deltaPhi(double(deltaPad) / double(N_PHI_CELLS) * 360.0);
  const std::vector<int>    testCases    = {0, 8, 16, 24, 32, 40, 48, 56};  //pads for eight directions
  const std::vector<double> expectedPhis = {
      0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0,  //backward
      0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0,  //forward
  };
  int counter = 0;
  for (int directionID = -1; directionID <= 1; directionID += 2) {
    // //double expectedTheta = (100.00 + 200.0/40*11.5)/3500.0;
    // if(directionID < 0) expectedTheta = M_PI - expectedTheta;
    for (auto padIDToFill : testCases) {
      streamlog_out(MESSAGE) << "**********************************************************************" << std::endl;
      const double        expectedPhi_1(expectedPhis[counter++]);
      const double        expectedPhi_2(expectedPhi_1 + deltaPhi);
      std::vector<double> expectedPhi{expectedPhi_1, expectedPhi_2};
      std::vector<double> expectedTheta(2, 0.0);  //filled later
      IMPL::LCEventImpl   myEvt;
      auto*               col = new IMPL::LCCollectionVec(LCIO::CALORIMETERHIT);
      col->parameters().setValue(LCIO::CellIDEncoding, "system:8,barrel:3,layer:8,slice:8,r:32:-16,phi:-16");
      HitMap hits;
      fillLumiCal(gmc, hits, padIDToFill, directionID == -1 ? 2 : directionID, expectedTheta[0]);
      fillLumiCal(gmc, hits, padIDToFill + deltaPad, directionID == -1 ? 2 : directionID, expectedTheta[1], 20);
      fillCollection(hits, col);

      myEvt.addCollection(col, "lumiCollection");

      lcc.processEvent(&myEvt);

      for (int armNow = -1; armNow < 2; armNow += 2) {
        streamlog_out(MESSAGE) << " Arm  " << std::setw(4) << armNow
                               << "\t Number of clusters: " << lcc._superClusterIdToCellId[armNow].size() << std::endl;
      }

      if (lcc._superClusterIdToCellId[directionID].size() != 2) {
        streamlog_out(ERROR) << "ERROR: Wrong number of reconstructed clusters" << std::endl;
        failed = true;
      }

      for (auto const& pairIDCells : lcc._superClusterIdToCellId[directionID]) {
        const int  clusterId       = pairIDCells.first;
        LCCluster& thisClusterInfo = lcc._superClusterIdClusterInfo[directionID][clusterId];
        thisClusterInfo.recalculatePositionFromHits(gmc);
        streamlog_out(MESSAGE) << thisClusterInfo << std::endl;
        bool thisSucceeds = checkReconstructedObject(gmc, thisClusterInfo, expectedPhi[0], expectedTheta[0]) or
                            checkReconstructedObject(gmc, thisClusterInfo, expectedPhi[1], expectedTheta[1]);
        if (not thisSucceeds) {
          streamlog_out(ERROR) << "ERROR: something is wrong here" << std::endl;

          failed = true;
        }
      }
      streamlog_out(MESSAGE) << "**********************************************************************" << std::endl;
    }
  }  // forwardAndBackward

  return failed ? 1 : 0;
}

int main() {
  streamlog::out.init(std::cout, "Testing");
  initLogLevels();

  streamlog::logscope scope(streamlog::out);
  scope.setLevel("DEBUG6");

  try {
    return testLumiCal();
  } catch (std::out_of_range& e) {
    std::cerr << "Out of Range error:" << e.what() << std::endl;
    return 1;
  } catch (std::runtime_error& e) {
    std::cerr << "Runtime Error: " << e.what() << std::endl;
    return 1;
  } catch (std::invalid_argument& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  return 0;
}
