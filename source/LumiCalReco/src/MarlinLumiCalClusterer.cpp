
#include "MarlinLumiCalClusterer.h"

#include "ClusterClass.h"
#include "Global.hh"
#include "LCCluster.hh"
#include "MCInfo.h"

#include <streamlog/loglevels.h>  // for DEBUG6, DEBUG7, DEBUG2
#include <streamlog/streamlog.h>  // for streamlog_out

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCIO.h>
#include <EVENT/MCParticle.h>
#include <Exceptions.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>

#include <TH1.h>
#include <TH2.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <map>
#include <memory>
#include <ostream>
#include <utility>
#include <vector>

// IWYU pragma: no_include <bits/shared_ptr.h>
// IWYU pragma: no_include <ext/alloc_traits.h>

/* >> */

/// use std::bind to set first variable and use in std::sort
bool ThetaCompAsc(const double sortTheta, SMCInfo const& a, SMCInfo const& b) {
  return fabs(sortTheta - a->theta) < fabs(sortTheta - b->theta);
}

bool PositionCompAsc(const double sortAgainstX, const double sortAgainstY, SMCInfo const& a, SMCInfo const& b) {
  double dxa = (sortAgainstX - a->x());
  double dya = (sortAgainstY - a->y());
  double dxb = (sortAgainstX - b->x());
  double dyb = (sortAgainstY - b->y());
  return (dxa*dxa+dya*dya) < (dxb*dxb+dyb*dyb)  ;
}

/*------------------------------------------------------------------------------------------*/
void MarlinLumiCalClusterer::TryMarlinLumiCalClusterer(EVENT::LCEvent* evt) {
  try {
    /* --------------------------------------------------------------------------
       create clusters using: LumiCalClustererClass
       -------------------------------------------------------------------------- */

    if (!LumiCalClusterer.processEvent(evt))
      return;

    LCCollectionVec* LCalClusterCol = new LCCollectionVec(LCIO::CLUSTER);
    IMPL::LCFlagImpl lcFlagImpl;
    lcFlagImpl.setBit(LCIO::CLBIT_HITS);
    LCalClusterCol->setFlag(lcFlagImpl.getFlag());

    LCCollectionVec* LCalRPCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

    streamlog_out(DEBUG6) << " Transfering reco results to LCalClusterCollection....." << std::endl;

    for (int armNow = -1; armNow < 2; armNow += 2) {
      streamlog_out(DEBUG6) << " Arm  " << std::setw(4) << armNow
                            << "\t Number of clusters: " << LumiCalClusterer._superClusterIdToCellId[armNow].size()
                            << std::endl;

      for (auto const& pairIDCells : LumiCalClusterer._superClusterIdToCellId[armNow]) {
        const int  clusterId       = pairIDCells.first;
        LCCluster& thisClusterInfo = LumiCalClusterer._superClusterIdClusterInfo[armNow][clusterId];
        thisClusterInfo.recalculatePositionFromHits(gmc);
        auto objectTuple(gmc.getLCIOObjects(thisClusterInfo, _minClusterEngy, _cutOnFiducialVolume));
        if (std::get<0>(objectTuple) == nullptr)
          continue;

        LCalClusterCol->addElement(std::get<0>(objectTuple));
        LCalRPCol->addElement(std::get<1>(objectTuple));
      }
    }

    //Add collections to the event if there are clusters
    if (LCalClusterCol->getNumberOfElements() != 0) {
      evt->addCollection(LCalClusterCol, LumiClusterColName);
      evt->addCollection(LCalRPCol, LumiRecoParticleColName);
    } else {
      delete LCalClusterCol;
      delete LCalRPCol;
    }

    //Optionally write information to ROOT tree
    writeRootInfo(evt);

  }  // try

  // if an !E!9exception has been thrown (no *col for this event) than do....
  catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG7) << "Event " << NumEvt << " has an exception" << std::endl;
  }

  return;
}

void MarlinLumiCalClusterer::writeRootInfo(LCEvent* evt) {
  if (OutRootFileName == "")
    return;

  // instantiate a clusterClass object for each mcParticle which was
  // created infront of LumiCal and was destroyed after lumical.
  std::map<int, MapIntPClusterClass> clusterClassMap;

  CreateClusters(clusterClassMap, evt);

  /* --------------------------------------------------------------------------
     histograming
     -------------------------------------------------------------------------- */

  for (int armNow = -1; armNow < 2; armNow += 2) {
    double totEngyIn  = 0.;
    double totEngyOut = 0.;

    for (MapIntPClusterClass::iterator mapIntClusterClassIt = clusterClassMap[armNow].begin();
         mapIntClusterClassIt != clusterClassMap[armNow].end(); ++mapIntClusterClassIt) {
      ClusterClass* thisCluster = mapIntClusterClassIt->second;
      const double  engyNow     = thisCluster->getEnergy();
      const double  thetaNow    = thisCluster->getTheta();
      // highest energy RecoParticle
      if (thisCluster->HighestEnergyFlag == 1) {
        OutputManager.HisMap1D["higestEngyParticle_Engy"]->Fill(engyNow);
        OutputManager.HisMap1D["higestEngyParticle_Theta"]->Fill(thetaNow);
      }
      if (thisCluster->OutsideFlag == 1) {
        // beyond acceptance ( energy and/or fid. vol.
        bool reason = (thisCluster->OutsideReason == "Reconstructed outside the fiducial volume");
        reason      = (reason || (thisCluster->OutsideReason == "Cluster energy below minimum"));

        if (reason) {
          OutputManager.HisMap2D["thetaEnergyOut_DepositedEngy"]->Fill(engyNow, thetaNow);
        }
        totEngyOut += engyNow;
      } else {
        // accepted
        totEngyIn += engyNow;
      }
    }
    // total energy inside LumiCal
    if (totEngyIn > 0) {
      OutputManager.HisMap1D["totEnergyIn"]->Fill(totEngyIn);
    }
    if (totEngyOut > 0) {
      OutputManager.HisMap1D["totEnergyOut"]->Fill(totEngyOut);
    }
  }

  /* --------------------------------------------------------------------------
     write criteria for selection cuts for Bhabha events into a tree
     -------------------------------------------------------------------------- */
  int clusterInFlag = 0;
  for (int armNow = -1; armNow < 2; armNow += 2) {
    for (MapIntPClusterClass::iterator mapIntClusterClassIt = clusterClassMap[armNow].begin();
         mapIntClusterClassIt != clusterClassMap[armNow].end(); ++mapIntClusterClassIt) {
      ClusterClass* thisCluster = mapIntClusterClassIt->second;

      if (thisCluster->OutsideFlag == 1)
        continue;
      // only take into account the highest-energy particle in the arm
      // (by default this means that this particle's shower is contained)
      //	  if(thisCluster->HighestEnergyFlag == 0)	continue;
      clusterInFlag++;

      //(BP) compute x,y,z in global reference system
      double globPos[] = {0.0, 0.0, 0.0};
      gmc.rotateToGlobal(thisCluster->getPosition(), globPos);
      //-------------
      OutputManager.TreeIntV["nEvt"]         = NumEvt;
      OutputManager.TreeIntV["outFlag"]      = thisCluster->OutsideFlag;
      OutputManager.TreeIntV["mFlag"]        = thisCluster->getPDG();
      OutputManager.TreeIntV["highE"]        = thisCluster->HighestEnergyFlag;
      OutputManager.TreeIntV["sign"]         = armNow;
      OutputManager.TreeIntV["nHits"]        = thisCluster->getNHits();
      OutputManager.TreeDoubleV["distTheta"] = thisCluster->DiffTheta;
      OutputManager.TreeDoubleV["distXY"]    = thisCluster->DiffPosXY;
      OutputManager.TreeDoubleV["engy"]      = thisCluster->getEnergy();
      OutputManager.TreeDoubleV["theta"]     = thisCluster->getTheta();
      OutputManager.TreeDoubleV["engyMC"]    = thisCluster->getEnergyMC();
      OutputManager.TreeDoubleV["thetaMC"]   = thisCluster->getThetaMC();
      OutputManager.TreeDoubleV["phi"]       = thisCluster->getPhi();
      OutputManager.TreeDoubleV["rzStart"]   = thisCluster->getRZStart();
      OutputManager.TreeDoubleV["Xglob"]     = globPos[0];
      OutputManager.TreeDoubleV["Yglob"]     = globPos[1];
      OutputManager.TreeDoubleV["Zglob"]     = globPos[2];
      //--
      OutputManager.FillRootTree("LumiRecoParticleTree");
    }
  }

  storeMCParticleInfo(evt, clusterInFlag);

#if _GLOBAL_COUNTERS_UPDATE_DEBUG == 1
  // write out the counter map
  int numCounters = OutputManager.Counter.size();

  if (numCounters > 0)
    streamlog_out(DEBUG1) << std::endl << "Global counters:" << std::endl;

  OutputManager.CounterIterator = OutputManager.Counter.begin();
  for (int hisNow = 0; hisNow < numCounters; hisNow++, OutputManager.CounterIterator++) {
    std::string counterName = (std::string)(*OutputManager.CounterIterator).first;
    streamlog_out(DEBUG1) << "\t" << OutputManager.Counter[counterName] << "  \t <->  " << counterName << std::endl;
  }
#endif

  /* --------------------------------------------------------------------------
     write to the root tree
     -------------------------------------------------------------------------- */
  OutputManager.WriteToRootTree("", NumEvt);

  /* --------------------------------------------------------------------------
     clean ClusterClassMap
     -------------------------------------------------------------------------- */
  for (int armNow = -1; armNow < 2; armNow += 2) {
    for (MapIntPClusterClass::iterator mapIntClusterClassIt = clusterClassMap[armNow].begin();
         mapIntClusterClassIt != clusterClassMap[armNow].end(); ++mapIntClusterClassIt) {
      delete mapIntClusterClassIt->second;
      mapIntClusterClassIt->second = nullptr;
    }
  }
}  // WriteRootInfo

void MarlinLumiCalClusterer::CreateClusters(std::map<int, MapIntPClusterClass>& clusterClassMap, EVENT::LCEvent* evt) {
  std::vector<SMCInfo> mcParticlesVecPos;
  std::vector<SMCInfo> mcParticlesVecNeg;

  streamlog_out(DEBUG7) << "CreateClusters: event = " << evt->getEventNumber() << std::endl;
  LCCollection* particles = evt->getCollection("MCParticle");
  if (particles) {
    int numMCParticles = particles->getNumberOfElements();
    if (numMCParticles) {
      streamlog_out(DEBUG6) << "CreateClusters: numMCParticles = " << numMCParticles << std::endl;
    } else {
      streamlog_out(DEBUG6) << "CreateClusters: No MCparticles in this event !" << std::endl;
    }
    for (int jparticle = 0; jparticle < numMCParticles; jparticle++) {
      MCParticle* particle = static_cast<MCParticle*>(particles->getElementAt(jparticle));
      auto        mcInfo   = MCInfo::getMCParticleInfo(particle, gmc);
      if (mcInfo->mcp == nullptr)
        continue;
      streamlog_out(DEBUG2) << mcInfo.get() << std::endl;
      if (mcInfo->sign > 0)
        mcParticlesVecPos.push_back(std::move(mcInfo));
      else
        mcParticlesVecNeg.push_back(std::move(mcInfo));
    }
  }
  int numOfClustersNeg = LumiCalClusterer._superClusterIdClusterInfo[-1].size();
  int numOfClustersPos = LumiCalClusterer._superClusterIdClusterInfo[1].size();
  if (numOfClustersNeg || numOfClustersPos) {
    streamlog_out(DEBUG6) << "Initial Set Stats for LumiCal......." << std::endl;
    streamlog_out(DEBUG6) << "     numOfClusters arm[-1] = " << numOfClustersNeg << "\t arm[1] = " << numOfClustersPos
                          << std::endl;
  } else {
    streamlog_out(DEBUG6) << "Initial Set Stats for LumiCal: No clusters found in this event !" << std::endl;
    return;
  }

  for (int armNow = -1; armNow < 2; armNow += 2) {
    double EngyMax   = 0.;
    int    EngyMaxID = -1;
    for (auto const& cluIDLCCluster : LumiCalClusterer._superClusterIdClusterInfo.at(armNow)) {
      const int   clusterId       = cluIDLCCluster.first;
      auto const& thisClusterInfo = cluIDLCCluster.second;
      // create a new cluster and put it on map
      ClusterClass* thisCluster = clusterClassMap[armNow][clusterId] = new ClusterClass(gmc, clusterId, &thisClusterInfo);
      if (thisCluster->getEnergy() > EngyMax) {
        EngyMax   = thisCluster->getEnergy();
        EngyMaxID = clusterId;
      }
      streamlog_out(DEBUG6) << "arm =   " << armNow << "\t cluster " << clusterId << "  ...... " << std::endl
                            << *thisCluster << std::endl;
    }
    // set flag for highest energy cluster found in this event
    if (EngyMax > 0.0)
      clusterClassMap[armNow][EngyMaxID]->HighestEnergyFlag = 1;
  }

  if (mcParticlesVecNeg.size() || mcParticlesVecPos.size()) {
    streamlog_out(DEBUG6) << "Transfering MC information into ClusterClass objects ......  " << std::endl;
    streamlog_out(DEBUG6) << "MC particles arm z-: " << mcParticlesVecNeg.size()
                          << "\t arm z+: " << mcParticlesVecPos.size() << std::endl;
  } else {
    streamlog_out(DEBUG6) << "Transfering MC information into ClusterClass objects: No primary MCParticles"
                          << " entering LumiCal found!" << std::endl;
    return;
  }

  for (int armNow = -1; armNow < 2; armNow += 2) {
    std::vector<SMCInfo>& mcParticlesVec = (armNow == -1) ? mcParticlesVecNeg : mcParticlesVecPos;

    for (MapIntPClusterClass::iterator mapIntClusterClassIt = clusterClassMap[armNow].begin();
         mapIntClusterClassIt != clusterClassMap[armNow].end(); ++mapIntClusterClassIt) {
      const int     clusterId   = mapIntClusterClassIt->first;
      ClusterClass* thisCluster = mapIntClusterClassIt->second;

      double xs = thisCluster->getPositionAtFront()[0];
      double ys = thisCluster->getPositionAtFront()[1];

      using namespace std::placeholders;  //_1, _2
      if (mcParticlesVec.size()) {
        // try to match MC true particle with cluster by comparing positions at Lumical entry
        // sort(mcParticlesVec.begin(), mcParticlesVec.end(),
        //      std::bind(BindThetaCompAsc, thisCluster->Theta, _1, _2));
        sort(mcParticlesVec.begin(), mcParticlesVec.end(), std::bind(PositionCompAsc, xs, ys, _1, _2));
        streamlog_out(DEBUG6) << "Trying to match particle: " << mcParticlesVec[0]->engy << std::endl;
        double dTheta = fabs(thisCluster->getTheta() - mcParticlesVec[0]->theta);
        double dPos0  = sqrt((sqr(xs - mcParticlesVec[0]->x()) + sqr(ys - mcParticlesVec[0]->y())));
        streamlog_out(DEBUG6) << "RZStart " << thisCluster->getRZStart() << std::endl;
        streamlog_out(DEBUG6) << "  xs, ys " << std::setw(13) << xs << std::setw(13) << ys << std::endl
                              << "MCP x, y " << std::setw(13) << mcParticlesVec[0]->x() << std::setw(13)
                              << mcParticlesVec[0]->y() << std::endl;

        double dEne0 = fabs(mcParticlesVec[0]->engy - thisCluster->getEnergy());
        streamlog_out(DEBUG6) << " dTheta " << std::setw(13) << dTheta << " dPos0 " << std::setw(13) << dPos0 << " dEne0 "
                              << std::setw(13) << dEne0 << std::endl;

        thisCluster->DiffTheta = dTheta;
        if (mcParticlesVec.size() > 1) {
          double dPos1 = sqrt((sqr(xs - mcParticlesVec[1]->x()) + sqr(ys - mcParticlesVec[1]->y())));
          double dEne1 = fabs(mcParticlesVec[1]->engy - thisCluster->getEnergy());
          if (dPos1 < gmc.GlobalParamD[GlobalMethodsClass::MinSeparationDist] && dEne1 < dEne0) {
            thisCluster->SetStatsMC(mcParticlesVec[1]);
            thisCluster->DiffPosXY = dPos1;
            mcParticlesVec.erase(mcParticlesVec.begin() + 1);
          } else {
            thisCluster->SetStatsMC(mcParticlesVec[0]);
            thisCluster->DiffPosXY = dPos0;
            mcParticlesVec.erase(mcParticlesVec.begin());
          }
        } else {
          thisCluster->DiffPosXY = dPos0;
          if (dPos0 < gmc.GlobalParamD[GlobalMethodsClass::MinSeparationDist]) {
            thisCluster->SetStatsMC(mcParticlesVec[0]);
            mcParticlesVec.erase(mcParticlesVec.begin());
          }
        }
      }

      if (thisCluster->getPDG() != 0) {
        double Ereco = thisCluster->getEnergy();
        // clang-format off
        streamlog_out( DEBUG6 ) << "\tParticle Out ("
                                << thisCluster -> OutsideReason << "):   " << clusterId << std::endl
                                << "\t\t side(arm), pdg, parentId , NumMCDaughters = "
                                << "\t" << thisCluster->getSign() <<"("<<armNow<<")"
                                << "\t" << thisCluster->getPDG()
                                << "\t" << thisCluster->getParentId()
                                << "\t" << thisCluster->getNumMCDaughters() << std::endl
                                << "\t\t MCPos    X, Y, Z: "
                                << "\t" << std::setw(13) << thisCluster->getMCPosition()[0]
                                << "\t" << std::setw(13) << thisCluster->getMCPosition()[1]
                                << "\t" << std::setw(13) << thisCluster->getMCPosition()[2]<< std::endl
                                << "\t\t Cluster  X, Y, Z: "
                                << "\t" << std::setw(13) << thisCluster->getPositionAtFront()[0]
                                << "\t" << std::setw(13) << thisCluster->getPositionAtFront()[1]
                                << "\t" << std::setw(13) << thisCluster->getPositionAtFront()[2] << std::endl
                                << "\t\t engy, theta, phi (mc): "
                                <<  std::setw(13)<< thisCluster->getEnergyMC()
                                << "\t"<< std::setw(13)<<thisCluster->getThetaMC()
                                << "\t"<< std::setw(13)<<thisCluster->getPhiMC() << std::endl
                                << "\t\t engy, theta, phi (rec):"
                                << std::setw(13)<<Ereco
                                << "\t"<< std::setw(13)<<thisCluster->getTheta()
                                << "\t"<< std::setw(13)<<thisCluster->getPhi()
                                << std::endl << std::endl;
        // clang-format on

      } else {
        streamlog_out(DEBUG6) << "Arm: " << armNow << "\t Cluster: " << clusterId << " does not match MC particle !\n";
        thisCluster->PrintInfo();
      }

    }  // loop over clusters
  }    // loop over arms
  return;
}  // CreateClusters

  void MarlinLumiCalClusterer::storeMCParticleInfo(LCEvent* evt, int clusterInFlag) {
    int mcparticleInFlag = 0;
    LCCollection* particles(nullptr);
    try {
      particles = evt->getCollection( "MCParticle" );
    } catch ( DataNotAvailableException &e ){
      streamlog_out(WARNING)<< "MCParticle data not available for event "<< NumEvt << std::endl;
      return;
    }

    const double LcalZstart = gmc.GlobalParamD[GlobalMethodsClass::ZStart];

    const int numMCParticles = particles->getNumberOfElements();
    for ( int jparticle=0; jparticle<numMCParticles; jparticle++ ){
      MCParticle *particle = static_cast<MCParticle*>( particles->getElementAt(jparticle) );
      auto        mcInfo   = MCInfo::getMCParticleInfo(particle, gmc);
      if (mcInfo->mcp == nullptr)
        continue;
      const double *endPoint = particle->getEndpoint();
      const double *vx = particle->getVertex();

      OutputManager.TreeIntV["nEvt"]	= NumEvt;
      OutputManager.TreeIntV["sign"]     = mcInfo->sign;
      OutputManager.TreeIntV["pdg"]      = mcInfo->pdg;
      OutputManager.TreeDoubleV["engy"]  = mcInfo->engy;
      OutputManager.TreeDoubleV["theta"] = mcInfo->theta;
      OutputManager.TreeDoubleV["phi"]   = mcInfo->phi;
      // position at the face of LCal (global coordinates)
      OutputManager.TreeDoubleV["begX"] = mcInfo->x();
      OutputManager.TreeDoubleV["begY"] = mcInfo->y();
      OutputManager.TreeDoubleV["begZ"] = double(mcInfo->sign) * LcalZstart;
      OutputManager.TreeDoubleV["vtxX"]	= vx[0];
      OutputManager.TreeDoubleV["vtxY"]	= vx[1];
      OutputManager.TreeDoubleV["vtxZ"]	= vx[2];
      OutputManager.TreeDoubleV["endX"] = endPoint[0];
      OutputManager.TreeDoubleV["endY"] = endPoint[1];
      OutputManager.TreeDoubleV["endZ"] = endPoint[2];
      // fill tree
      OutputManager.FillRootTree("LumiMCParticleTree" );

    }// for jparticle
    OutputManager.HisMap2D["NumClustersIn_numMCParticlesIn"] -> Fill ( mcparticleInFlag, clusterInFlag );

}
