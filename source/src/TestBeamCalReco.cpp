#include "BCPCuts.hh"
#include "BCPadEnergies.hh"
#include "BCUtilities.hh"
#include "BeamCal.hh"
#include "BeamCalCluster.hh"
#include "BeamCalGeoDD.hh"

#include <iomanip>
#include <iostream>
#include <vector>

/// Inject Energy into the beamcal vector
void fillPadWithEnergy(BCPadEnergies& pads, int padID) {
  const int    ring   = 11;  // is a full ring
  const double energy = 1.1;
  for (int layer = 1; layer < 11; ++layer) {
    pads.addEnergy(layer, ring, padID, energy);
  }
}

int testBeamCal(int argn, char** argc) {
  if (argn < 2) {
    throw std::invalid_argument("Not enough parameters\nTestBeamCalReco compactFile");
  }

  std::string compactFile(argc[1]);

  BCPCuts cuts;
  cuts.setStartLayer(1).setSigmaCut(0.01).setMinimumTowerSize(4);

  auto& theDetector = dd4hep::Detector::getInstance();
  theDetector.fromCompact(compactFile);
  BeamCalGeo* geo = new BeamCalGeoDD(theDetector);

  //empty background
  std::vector<BCPadEnergies> backgroundBeamCals(2, geo);

  std::vector<int>    testCases    = {0, 12, 24, 36, 48, 60, 72, 84};  //ring 11, pads for eight directions
  std::vector<double> expectedPhis = {
      180.0, 225.0, 270.0, 315.0, 0.0, 45.0,  90.0,  135.0,  //forward, "Left"
      180.0, 135.0, 90.0,  45.0,  0.0, 315.0, 270.0, 225.0,  //backward "Right"
  };

  int counter = 0;
  for (int directionID = 0; directionID <= 1; ++directionID) {
    for (auto padIDToFill : testCases) {
      std::cout << "**********************************************************************" << std::endl;

      auto                       expectedPhi = expectedPhis[counter++];
      std::vector<BCPadEnergies> signalBeamCals(2, geo);
      fillPadWithEnergy(signalBeamCals[directionID], padIDToFill);

      signalBeamCals[0].setSide(BCPadEnergies::BeamCalSide_t::kLeft);
      signalBeamCals[1].setSide(BCPadEnergies::BeamCalSide_t::kRight);

      bool                               detailedPrintout;
      const std::vector<BeamCalCluster>& clusters0 =
          signalBeamCals[0].lookForNeighbouringClustersOverSigma(backgroundBeamCals[0], cuts, detailedPrintout = true);
      const std::vector<BeamCalCluster>& clusters1 =
          signalBeamCals[1].lookForNeighbouringClustersOverSigma(backgroundBeamCals[1], cuts, detailedPrintout = true);

      double         reconstructedPhi = -100;
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

int main(int argn, char** argc) {
  try {
    return testBeamCal(argn, argc);
  } catch (std::out_of_range& e) {
    std::cerr << "Geometry does not agree with energy in the trees:" << e.what() << std::endl;
    return 1;
  } catch (std::invalid_argument& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  return 0;
}
