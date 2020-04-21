#ifndef TESTUTILITIES_HH
#define TESTUTILITIES_HH 1

#include "BCUtilities.hh"
#include "GlobalMethodsClass.h"

#include <EVENT/ReconstructedParticle.h>

#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>

#include <streamlog/streamlog.h>
#include <marlin/VerbosityLevels.h>

#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <iostream>


using HitMap = std::map<long long, IMPL::CalorimeterHitImpl*>;

constexpr int N_PHI_CELLS = 64;


template <class T> double getTheta(GlobalMethodsClass const& gmc, T* p);

template <> double getTheta(GlobalMethodsClass const& gmc, const double* momGlob) {
  double mom[3] = {0.0, 0.0, 0.0};
  gmc.rotateToLumiCal(momGlob, mom);
  const double r     = sqrt(mom[0] * mom[0] + mom[1] * mom[1]);
  double       theta = std::atan(r / mom[2]);
  if (theta < 0)
    theta += M_PI;
  return theta;
}

template <class T> double getTheta(GlobalMethodsClass const& gmc, T* p) {
  const double* momGlob = p->getMomentum();
  return getTheta(gmc, momGlob);
}

template <class T> double getPhi(GlobalMethodsClass const& gmc, T* p) {
  const double* momGlob = p->getMomentum();
  double        mom[3]  = {0.0, 0.0, 0.0};
  gmc.rotateToLumiCal(momGlob, mom);
  return std::atan2(mom[1], mom[0]);
}

/// Inject Energy into the LumiCal collection
void fillLumiCal(GlobalMethodsClass& gmc, HitMap& hits, int phiID, int direction,
                 double& expectedTheta, const int centralR=20) {
  const double energy       = 1.0;
  const double centralLayer = 11;
  const int    nAddHits     = 2;
  const double scaleFactor  = 0.05;
  expectedTheta             = 0.0;
  double weightSum          = 0.0;

  int sign = direction == 2 ? -1 : 1;

  for (int layer = 0; layer < 40; ++layer) {
    for (int thisPad = phiID - nAddHits; thisPad <= phiID + nAddHits; ++thisPad) {
      for (int thisR = centralR - nAddHits; thisR <= centralR + nAddHits; ++thisR) {
        auto hit = std::make_unique<IMPL::CalorimeterHitImpl>();

        int thisPadID = thisPad;
        // counting sense inverted in backward direction (because of LumiCalConstruction)
        if (sign < 0) {
          thisPadID = N_PHI_CELLS - thisPadID;
        }
        //phi goes from -NPhi/2 to NPhi/2-1
        if (thisPadID >= N_PHI_CELLS / 2)
          thisPadID -= N_PHI_CELLS;
        if (thisPadID < -N_PHI_CELLS / 2)
          thisPadID += N_PHI_CELLS;

        double hitEnergy = energy * exp(-(layer - centralLayer) * (layer - centralLayer) * scaleFactor) *
            exp(-(phiID - thisPad) * (phiID - thisPad) * scaleFactor) *
            exp(-(centralR - thisR) * (centralR - thisR) * scaleFactor);
        
        long cellID0 = (3 << 0) + (direction << 8) + (layer << (8 + 3));
        long cellID1 = thisR + (thisPad << 16);

        unsigned long long llCell = (unsigned long long)(cellID0) + (unsigned long long)(cellID1 << 32);

        auto existingHitIt = hits.find(llCell);
        if(existingHitIt != hits.end()){
          auto* existingHit = existingHitIt->second;
          streamlog_out(DEBUG2) << "Found existing hit "
                                << std::setw(15) << llCell
                                << std::setw(13) << layer
                                << std::setw(13) << thisPadID
                                << std::setw(13) << thisR
                                << std::setw(13) << hitEnergy
                                << std::setw(13) << existingHit->getEnergy()
                                << std::endl;
          existingHit->setEnergy(existingHit->getEnergy() + hitEnergy);
          continue;
        }

        hit->setCellID0(cellID0);
        hit->setCellID1(cellID1);
        hit->setEnergy(hitEnergy);

        std::map<GlobalMethodsClass::Coordinate_t, double> thetaPhi;

        float rCell = 100.00 + 200.0 / 40 * (thisR + 0.5);
        float zCell = 3500.0 + 500.0 / 40 * layer;

        float phiCell = thisPadID * 2.0 * M_PI / N_PHI_CELLS;
        if (sign < 0) {
          phiCell *= -1.0;
        }

        double locPos[3]  = {rCell * cos(phiCell), rCell * sin(phiCell), sign * zCell};
        double globPos[3] = {0.0, 0.0, 0.0};
        gmc.rotateToGlobal(locPos, globPos);

        float floatPos[3] = {float(globPos[0]), float(globPos[1]), float(globPos[2])};

        hit->setPosition(floatPos);

        weightSum += energy;
        expectedTheta += energy * getTheta<const double>(gmc, globPos);

        streamlog_out(DEBUG2) << "Creating hit "
                              << std::setw(15) << llCell
                              << std::setw(13) << layer
                              << std::setw(13) << thisPadID
                              << std::setw(13) << thisR
                              << std::setw(13) << locPos[0]
                              << std::setw(13) << locPos[1]
                              << std::setw(13) << locPos[2]
                              << std::setw(13) << phiCell*180.0/M_PI
                              << std::setw(13) << hit->getEnergy()
                              << std::endl;
        
        hits[llCell] = hit.release();
      }
    }
  }

  expectedTheta /= weightSum;
}

void fillCollection(HitMap const& hits, IMPL::LCCollectionVec* col){
  for (auto const& hitPair : hits ) {
    col->addElement(hitPair.second);
  }
}

void fillLumiCal(GlobalMethodsClass& gmc, IMPL::LCCollectionVec* col, int phiID, int direction,
                       double& expectedTheta) {
  HitMap hits;
  fillLumiCal(gmc, hits, phiID, direction, expectedTheta);
  fillCollection(hits, col);
}



void initializeGMC(GlobalMethodsClass& gmc){
  
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

}

bool checkReconstructedObject(GlobalMethodsClass const& gmc, LCCluster& thisClusterInfo,
                              double expectedPhi, double expectedTheta) {

  bool success = true;

  auto objectTuple(gmc.getLCIOObjects(thisClusterInfo, 0.0, false));
  if (std::get<0>(objectTuple) == nullptr) {
    streamlog_out(ERROR) << "FAILED To reconstruct" << std::endl;
    success = false;
  }

  auto* cluster = std::get<0>(objectTuple);
  auto* rp      = std::get<1>(objectTuple);

  double reconstructedPhi   = getPhi(gmc, rp) * 180 / M_PI;
  double reconstructedTheta = getTheta(gmc, rp);

  if (reconstructedPhi < 0)
    reconstructedPhi += 360;

  streamlog_out(MESSAGE) << "Expected and Reconstructed Angles: "
                         << std::setw(13) << expectedPhi << std::setw(13) << reconstructedPhi
                         << std::setw(13) << expectedTheta << std::setw(13) << reconstructedTheta << std::endl;

  if (not BCUtil::areCloseTogether(0.1, expectedPhi, 0.1, reconstructedPhi)) {
    success = false;
  }

  if (fabs(expectedTheta - reconstructedTheta) > 0.1) {
    success = false;
  }

  delete cluster;
  delete rp;

  return success;
}

void initLogLevels(){
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
}

#endif // TESTUTILITIES_HH
