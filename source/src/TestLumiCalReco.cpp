#include "BCUtilities.hh"
#include "GlobalMethodsClass.h"
#include "LCCluster.hh"
#include "LumiCalClusterer.h"

#include <streamlog/streamlog.h>

#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>

#include <EVENT/ReconstructedParticle.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

constexpr int N_PHI_CELLS = 64;

template <class t> inline void RotateToLumiCal(const t* vector, t* rotated, double angle) {
  if (vector[2] < 0) {
    angle *= -1;
  }
  const double CosAngle = cos(double(angle) / 1000.0);
  const double SinAngle = sin(double(angle) / 1000.0);
  rotated[0]            = CosAngle * vector[0] - SinAngle * vector[2];
  rotated[1]            = vector[1];
  rotated[2]            = SinAngle * vector[0] + CosAngle * vector[2];
}

template <class T> double getTheta(T& p);
template <class T> double getTheta(T* p);

template <> double getTheta(const double* momGlob) {
  double mom[3] = {0.0, 0.0, 0.0};
  RotateToLumiCal(momGlob, mom, 10.0);
  const double r     = sqrt(mom[0] * mom[0] + mom[1] * mom[1]);
  double       theta = std::atan(r / mom[2]);
  if (theta < 0)
    theta += M_PI;
  return theta;
}

template <class T> double getTheta(T* p) {
  const double* momGlob = p->getMomentum();
  return getTheta(momGlob);
}

template <class T> double getPhi(T* p) {
  const double* momGlob = p->getMomentum();
  double        mom[3]  = {0.0, 0.0, 0.0};
  RotateToLumiCal(momGlob, mom, 10.0);
  return std::atan2(mom[1], mom[0]);
}

/// Inject Energy into the beamcal vector
void fillPadWithEnergy(GlobalMethodsClass& gmc, IMPL::LCCollectionVec& collection, int phiID, int direction,
                       double& expectedTheta) {
  const double energy       = 1.0;
  const double centralLayer = 11;
  const int    nAddHits     = 2;
  const double scaleFactor  = 0.05;
  expectedTheta             = 0.0;
  double weightSum          = 0.0;

  int sign = direction == 2 ? -1 : 1;

  for (int layer = 0; layer < 40; ++layer) {
    for (int thisPad = phiID - nAddHits; thisPad <= phiID + nAddHits; ++thisPad) {
      for (int thisR = 20 - nAddHits; thisR <= 20 + nAddHits; ++thisR) {
        auto* hit = new IMPL::CalorimeterHitImpl();

        int thisPadID = thisPad;
        if (thisPadID < 0)
          thisPadID += N_PHI_CELLS;
        if (thisPadID >= N_PHI_CELLS)
          thisPadID -= N_PHI_CELLS;

        long cellID0 = (3 << 0) + (direction << 8) + (layer << (8 + 3));
        long cellID1 = thisR + (thisPad << 16);

        hit->setCellID0(cellID0);
        hit->setCellID1(cellID1);
        hit->setEnergy(energy * exp(-(layer - centralLayer) * (layer - centralLayer) * scaleFactor) *
                       exp(-(phiID - thisPad) * (phiID - thisPad) * scaleFactor) *
                       exp(-(20 - thisR) * (20 - thisR) * scaleFactor));

        std::map<GlobalMethodsClass::Coordinate_t, double> thetaPhi;

        float rCell = 100.00 + 200.0 / 40 * (thisR + 0.5);
        float zCell = 3500.0 + 500.0 / 40 * layer;

        float phiCell = thisPadID * 2.0 * M_PI / N_PHI_CELLS;

        double locPos[3]  = {rCell * cos(phiCell), rCell * sin(phiCell), sign * zCell};
        double globPos[3] = {0.0, 0.0, 0.0};
        gmc.rotateToGlobal(locPos, globPos);

        float floatPos[3] = {float(globPos[0]), float(globPos[1]), float(globPos[2])};

        hit->setPosition(floatPos);

        weightSum += energy;
        expectedTheta += energy * getTheta<const double>(globPos);

        streamlog_out(DEBUG2) << "Creating hit "
          
                              << std::setw(13) << layer
                              << std::setw(13) << thisPadID
                              << std::setw(13) << thisR
                              << std::setw(13) << locPos[0]
                              << std::setw(13) << locPos[1]
                              << std::setw(13) << locPos[2]
                              << std::setw(13) << phiCell*180.0/M_PI
                              << std::setw(13) << hit->getEnergy()
                              << std::endl;
        
        collection.addElement(hit);
      }
    }
  }

  expectedTheta /= weightSum;
}

std::tuple<ClusterImpl*, ReconstructedParticleImpl*> getLCIOObjects(GlobalMethodsClass const& gmc,
                                                                    LCCluster const& thisClusterInfo) {
  const double clusterEnergy = thisClusterInfo.getE();

  ClusterImpl* cluster = new ClusterImpl;
  cluster->setEnergy(clusterEnergy);

  ReconstructedParticleImpl* particle = new ReconstructedParticleImpl;
  const float                mass     = 0.0;
  const float                charge   = 1e+19;
  particle->setMass(mass);
  particle->setCharge(charge);
  particle->setEnergy(clusterEnergy);
  particle->addCluster(cluster);

  const float locPos[3] = {float(thisClusterInfo.getX()), float(thisClusterInfo.getY()), float(thisClusterInfo.getZ())};
  float       gP[3]     = {0.0, 0.0, 0.0};
  gmc.rotateToGlobal(locPos, gP);
  cluster->setPosition(gP);

  const float norm               = clusterEnergy / sqrt(gP[0] * gP[0] + gP[1] * gP[1] + gP[2] * gP[2]);
  const float clusterMomentum[3] = {float(gP[0] * norm), float(gP[1] * norm), float(gP[2] * norm)};
  particle->setMomentum(clusterMomentum);
  for (auto& lumicalHit : thisClusterInfo.getCaloHits()) {
    for (auto* hit : lumicalHit->getHits()) {
      cluster->addHit(hit, 1.0);
    }
  }
  cluster->subdetectorEnergies().resize(6);
  //LCAL_INDEX=3 in DDPFOCreator.hh
  cluster->subdetectorEnergies()[3] = clusterEnergy;

  return std::make_tuple(cluster, particle);
}

int testLumiCal() {
  bool               failed = false;
  GlobalMethodsClass gmc;

  gmc.GlobalParamD[GlobalMethodsClass::RMin]              = 100;
  gmc.GlobalParamD[GlobalMethodsClass::RMax]              = 300;
  gmc.GlobalParamD[GlobalMethodsClass::ZStart]            = 2500;
  gmc.GlobalParamD[GlobalMethodsClass::ZEnd]              = 3000;
  gmc.GlobalParamD[GlobalMethodsClass::RCellLength]       = 0.5;
  gmc.GlobalParamD[GlobalMethodsClass::RCellOffset]       = 0.0;
  gmc.GlobalParamD[GlobalMethodsClass::PhiCellLength]     = 2 * M_PI / double(N_PHI_CELLS);
  gmc.GlobalParamD[GlobalMethodsClass::PhiCellOffset]     = 0.0;  //no offset
  gmc.GlobalParamI[GlobalMethodsClass::NumCellsR]         = 40;
  gmc.GlobalParamI[GlobalMethodsClass::NumCellsPhi]       = N_PHI_CELLS;
  gmc.GlobalParamI[GlobalMethodsClass::NumCellsZ]         = 40;
  gmc.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle] = 0.020;

  gmc.GlobalParamD[GlobalMethodsClass::ZLayerThickness] = 500.0 / 40.0;
  gmc.GlobalParamD[GlobalMethodsClass::ZLayerZOffset]   = 0.0;

  gmc.GlobalParamS[GlobalMethodsClass::WeightingMethod] = "LogMethod";

  gmc.GlobalParamI[GlobalMethodsClass::ClusterMinNumHits]                = 3;
  gmc.GlobalParamD[GlobalMethodsClass::MinHitEnergy]                     = 5e-6;
  gmc.GlobalParamD[GlobalMethodsClass::ZLayerThickness]                  = 4.5;
  gmc.GlobalParamD[GlobalMethodsClass::ZLayerPhiOffset]                  = 0.0;
  gmc.GlobalParamD[GlobalMethodsClass::ElementsPercentInShowerPeakLayer] = 0.03;
  gmc.GlobalParamI[GlobalMethodsClass::NumOfNearNeighbor]                = 6;
  gmc.GlobalParamD[GlobalMethodsClass::MiddleEnergyHitBoundFrac]         = 0.01;
  gmc.GlobalParamD[GlobalMethodsClass::LogWeightConstant]                = 6.1;
  gmc.GlobalParamD[GlobalMethodsClass::MoliereRadius]                    = 20;
  gmc.GlobalParamD[GlobalMethodsClass::MinSeparationDist]                = 20;
  gmc.GlobalParamD[GlobalMethodsClass::MinClusterEngyGeV]                = 0.1;
  gmc.GlobalParamD[GlobalMethodsClass::ThetaMin]                         = 40;
  gmc.GlobalParamD[GlobalMethodsClass::ThetaMax]                         = 100;

  gmc.initializeAdditionalParameters();

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
      fillPadWithEnergy(gmc, *col, padIDToFill, directionID == -1 ? 2 : directionID, expectedTheta);
      myEvt.addCollection(col, "lumiCollection");

      lcc.processEvent(&myEvt);

      for (int armNow = -1; armNow < 2; armNow += 2) {
        streamlog_out(MESSAGE) << " Arm  " << std::setw(4) << armNow
                               << "\t Number of clusters: " << lcc._superClusterIdToCellId[armNow].size() << std::endl;
      }

      for (auto const& pairIDCells : lcc._superClusterIdToCellId[directionID]) {
        streamlog_out(MESSAGE) << "Got some clusterInfo" << std::endl;
        const int  clusterId       = pairIDCells.first;
        LCCluster& thisClusterInfo = lcc._superClusterIdClusterInfo[directionID][clusterId];
        thisClusterInfo.recalculatePositionFromHits(gmc);

        streamlog_out(MESSAGE) << thisClusterInfo << std::endl;

        auto objectTuple(getLCIOObjects(gmc, thisClusterInfo));
        if (std::get<0>(objectTuple) == nullptr) {
          streamlog_out(ERROR) << "FAILED To reconstruct" << std::endl;
          failed = true;
        }

        auto* cluster = std::get<0>(objectTuple);
        auto* rp      = std::get<1>(objectTuple);

        double reconstructedPhi   = getPhi(rp) * 180 / M_PI;
        double reconstructedTheta = getTheta(rp);

        if (reconstructedPhi < 0)
          reconstructedPhi += 360;

        streamlog_out(MESSAGE) << "Expected and Reconstructed Angles: " << std::endl
                               << std::setw(13) << expectedPhi << std::setw(13) << reconstructedPhi << std::endl
                               << std::setw(13) << expectedTheta << std::setw(13) << reconstructedTheta << std::endl;

        if (not BCUtil::areCloseTogether(0.1, expectedPhi, 0.1, reconstructedPhi)) {
          streamlog_out(ERROR) << "ERROR: Phi not correctly reconstructed: (exp, rec) " << expectedPhi << "   "
                               << reconstructedPhi << std::endl;
          failed = true;
        }

        if (fabs(expectedTheta - reconstructedTheta) > 0.1) {
          streamlog_out(ERROR) << "ERROR: Theta not correctly reconstructed: (exp, rec) " << expectedTheta << "   "
                               << reconstructedTheta << std::endl;
          failed = true;
        }

        delete cluster;
        delete rp;
      }
      streamlog_out(MESSAGE) << "**********************************************************************" << std::endl;
    }
  }  // forwardAndBackward

  return failed ? 1 : 0;
}

int main() {
  streamlog::out.init(std::cout, "Testing");
  streamlog::out.addLevelName<marlin::DEBUG>();
  streamlog::out.addLevelName<marlin::DEBUG0>();
  streamlog::out.addLevelName<marlin::DEBUG1>();
  streamlog::out.addLevelName<marlin::DEBUG2>();
  streamlog::out.addLevelName<marlin::DEBUG3>();
  streamlog::out.addLevelName<marlin::DEBUG4>();
  streamlog::out.addLevelName<marlin::DEBUG5>();
  streamlog::out.addLevelName<marlin::DEBUG6>();
  streamlog::out.addLevelName<marlin::DEBUG7>();
  streamlog::out.addLevelName<marlin::DEBUG8>();
  streamlog::out.addLevelName<marlin::DEBUG9>();
  streamlog::out.addLevelName<marlin::MESSAGE>();
  streamlog::out.addLevelName<marlin::MESSAGE0>();
  streamlog::out.addLevelName<marlin::MESSAGE1>();
  streamlog::out.addLevelName<marlin::MESSAGE2>();
  streamlog::out.addLevelName<marlin::MESSAGE3>();
  streamlog::out.addLevelName<marlin::MESSAGE4>();
  streamlog::out.addLevelName<marlin::MESSAGE5>();
  streamlog::out.addLevelName<marlin::MESSAGE6>();
  streamlog::out.addLevelName<marlin::MESSAGE7>();
  streamlog::out.addLevelName<marlin::MESSAGE8>();
  streamlog::out.addLevelName<marlin::MESSAGE9>();
  streamlog::out.addLevelName<marlin::WARNING>();
  streamlog::out.addLevelName<marlin::WARNING0>();
  streamlog::out.addLevelName<marlin::WARNING1>();
  streamlog::out.addLevelName<marlin::WARNING2>();
  streamlog::out.addLevelName<marlin::WARNING3>();
  streamlog::out.addLevelName<marlin::WARNING4>();
  streamlog::out.addLevelName<marlin::WARNING5>();
  streamlog::out.addLevelName<marlin::WARNING6>();
  streamlog::out.addLevelName<marlin::WARNING7>();
  streamlog::out.addLevelName<marlin::WARNING8>();
  streamlog::out.addLevelName<marlin::WARNING9>();
  streamlog::out.addLevelName<marlin::ERROR>();
  streamlog::out.addLevelName<marlin::ERROR0>();
  streamlog::out.addLevelName<marlin::ERROR1>();
  streamlog::out.addLevelName<marlin::ERROR2>();
  streamlog::out.addLevelName<marlin::ERROR3>();
  streamlog::out.addLevelName<marlin::ERROR4>();
  streamlog::out.addLevelName<marlin::ERROR5>();
  streamlog::out.addLevelName<marlin::ERROR6>();
  streamlog::out.addLevelName<marlin::ERROR7>();
  streamlog::out.addLevelName<marlin::ERROR8>();
  streamlog::out.addLevelName<marlin::ERROR9>();
  streamlog::out.addLevelName<marlin::SILENT>();

  streamlog::logscope scope(streamlog::out);
  scope.setLevel("MESSAGE");

  try {
    testLumiCal();
  } catch (std::out_of_range& e) {
    std::cerr << "Out of Range error:" << e.what() << std::endl;
    return 1;
  } catch (std::invalid_argument& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  return 0;
}
