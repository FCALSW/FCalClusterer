#include "TestMain.hh"

#include "BCPCuts.hh"
#include "BCPadEnergies.hh"
#include "BCUtilities.hh"
#include "BeamCalCluster.hh"
#include "BeamCalGeoDD.hh"

#include <DD4hep/Detector.h>

#include <iomanip>
#include <iostream>
#include <vector>

/// Inject Energy into the beamcal vector
void fillPadWithEnergy(BCPadEnergies& pads, int ring, int padID, int maxPadID) {
  const double energy = 1.1;
  int          padM   = padID - 1;
  int          padP   = padID + 1;
  if (padM < 0)
    padM += maxPadID;
  if (padP >= maxPadID)
    padP -= maxPadID;

  for (int layer = 1; layer < 11; ++layer) {
    pads.addEnergy(layer, ring, padM, energy * 0.5);
    pads.addEnergy(layer, ring, padID, energy);
    pads.addEnergy(layer, ring, padP, energy * 0.5);
  }
}

int runTest(int argn, char** argc) {
  if (argn < 3) {
    throw std::invalid_argument("Not enough parameters\nTestBeamCalReco compactFile DetectorName");
  }

  std::string compactFile(argc[1]);
  std::string detectorName(argc[2]);
  std::string colName(argc[2]);
  colName += "Collection";

  BCPCuts cuts;
  cuts.setStartLayer(1).setSigmaCut(0.01).setMinimumTowerSize(4);

  auto& theDetector = dd4hep::Detector::getInstance();
  theDetector.fromCompact(compactFile);
  BeamCalGeo* geo = new BeamCalGeoDD(theDetector, detectorName, colName);

  //empty background
  std::vector<BCPadEnergies> backgroundBeamCals(2, geo);

  const int           ring         = 11;
  std::vector<int>    testCases    = {0, 12, 24, 36, 48, 60, 72, 84};  //ring 11, pads for eight directions
  std::vector<double> expectedPhis = {
      180.0, 225.0, 270.0, 315.0, 0.0, 45.0,  90.0,  135.0,  //forward, "Left"
      180.0, 135.0, 90.0,  45.0,  0.0, 315.0, 270.0, 225.0,  //backward "Right"
  };

  if (detectorName == "LumiCal") {
    testCases    = {0, 6, 12, 18, 24, 30, 36, 42};  //48 pads in phi
    expectedPhis = {
        0.0, 45.0,  90.0,  135.0,  180.0,  225.0,  270.0,  315.0,   // forward
        0.0, -45.0, -90.0, -135.0, -180.0, -225.0, -270.0, -315.0,  // backward
    };
  }
  streamlog_out(MESSAGE) << "Geometry\n" << *geo << std::endl;

  int counter = 0;
  for (int directionID = 0; directionID <= 1; ++directionID) {
    for (auto padIDToFill : testCases) {
      std::cout << "**********************************************************************" << std::endl;

      auto                       expectedPhi = expectedPhis[counter++];
      std::vector<BCPadEnergies> signalBeamCals(2, geo);
      fillPadWithEnergy(signalBeamCals[directionID], ring, padIDToFill, geo->getPadsInRing(ring));

      signalBeamCals[0].setSide(BCPadEnergies::BeamCalSide_t::kLeft);
      signalBeamCals[1].setSide(BCPadEnergies::BeamCalSide_t::kRight);

      bool                               detailedPrintout;
      const std::vector<BeamCalCluster>& clusters0 =
          signalBeamCals[0].lookForNeighbouringClustersOverSigma(backgroundBeamCals[0], cuts, detailedPrintout = true);
      const std::vector<BeamCalCluster>& clusters1 =
          signalBeamCals[1].lookForNeighbouringClustersOverSigma(backgroundBeamCals[1], cuts, detailedPrintout = true);

      double         reconstructedPhi = -1000;
      BeamCalCluster reco;
      ///Output the reconstructed objects
      for (std::vector<BeamCalCluster>::const_iterator it = clusters0.begin(); it != clusters0.end(); ++it) {
        std::cout << "ReconstructedCluster:\n" << *it << std::endl;
        reconstructedPhi = it->getPhi();
      }

      for (std::vector<BeamCalCluster>::const_iterator it = clusters1.begin(); it != clusters1.end(); ++it) {
        std::cout << "ReconstructedCluster:\n" << *it << std::endl;
        reconstructedPhi = it->getPhi();
      }

      std::cout << "Expected: " << std::setw(10) << expectedPhi                //
                << "\t\tReconstructed: " << std::setw(10) << reconstructedPhi  //
                << std::endl;

      if (not BCUtil::areCloseTogether(0.1, expectedPhi, 0.1, reconstructedPhi)) {
        std::cout << "ERROR: phi not correctly reconstructed:" << std::endl;
      }

      std::cout << "**********************************************************************" << std::endl;
    }
  }  // forwardAndBackward

  return 0;
}
