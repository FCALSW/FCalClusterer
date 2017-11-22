
#include "MarlinLumiCalClusterer.h"

#include "ClusterClass.h"
#include "MCInfo.h"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>

#include <algorithm>
#include <map>
#include <vector>
#include <iomanip>

/* >> */

/// use std::bind to set first variable and use in std::sort
bool ThetaCompAsc(const double sortTheta, const MCInfo& a, const MCInfo& b) {
  return fabs(sortTheta - a.theta) < fabs(sortTheta - b.theta);
}

bool PositionCompAsc(const double sortAgainstX, const double sortAgainstY, const MCInfo& a, const MCInfo& b) {
  double dxa = (sortAgainstX - a.x);
  double dya = (sortAgainstY - a.y);
  double dxb = (sortAgainstX - b.x);
  double dyb = (sortAgainstY - b.y);
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
        auto objectTuple(getLCIOObjects(thisClusterInfo));
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
 

// (BP) this is optional
// empty string OutRootFileName in steering file

      if ( OutRootFileName != "" ) {
        // instantiate a clusterClass object for each mcParticle which was
        // created infront of LumiCal and was destroyed after lumical.
        std::map<int, MapIntPClusterClass> clusterClassMap;

	CreateClusters( LumiCalClusterer._superClusterIdToCellId,
			LumiCalClusterer._superClusterIdToCellEngy,
			clusterClassMap,
			evt );
      
      /* --------------------------------------------------------------------------
	 histograming
	 -------------------------------------------------------------------------- */

     for(int armNow = -1; armNow < 2; armNow += 2) {
	double totEngyIn = 0.;
	double totEngyOut= 0.;

	for (MapIntPClusterClass::iterator mapIntClusterClassIt = clusterClassMap[armNow].begin();
	     mapIntClusterClassIt != clusterClassMap[armNow].end();
	     ++mapIntClusterClassIt) {

	  ClusterClass* thisCluster = mapIntClusterClassIt->second;
	  const double engyNow  = thisCluster->Engy;
	  const double thetaNow = thisCluster -> Theta;
	  // highest energy RecoParticle
	  if(thisCluster->HighestEnergyFlag == 1){
	    OutputManager.HisMap1D["higestEngyParticle_Engy"] -> Fill (engyNow);
	    OutputManager.HisMap1D["higestEngyParticle_Theta"] -> Fill (thetaNow);
	  }
	  if(thisCluster->OutsideFlag == 1){
	    // beyond acceptance ( energy and/or fid. vol.
	    bool reason = (thisCluster->OutsideReason == "Reconstructed outside the fiducial volume" );
	         reason = ( reason || ( thisCluster->OutsideReason == "Cluster energy below minimum" ) );

	    if( reason ) {
	      OutputManager.HisMap2D["thetaEnergyOut_DepositedEngy"] -> Fill (engyNow , thetaNow);
	    }
	    totEngyOut += engyNow;
	  }else{ 
	    // accepted
	    totEngyIn += engyNow;
	  }
	}
	// total energy inside LumiCal
	 if(totEngyIn > 0) {
	   OutputManager.HisMap1D["totEnergyIn"] -> Fill (totEngyIn);
	 }
	 if(totEngyOut > 0) {
	   OutputManager.HisMap1D["totEnergyOut"] -> Fill (totEngyOut);
	 }
      }


      /* --------------------------------------------------------------------------
	 write criteria for selection cuts for Bhabha events into a tree
	 -------------------------------------------------------------------------- */
      int	clusterInFlag = 0;
      for(int armNow = -1; armNow < 2; armNow += 2) {

	for (MapIntPClusterClass::iterator mapIntClusterClassIt = clusterClassMap[armNow].begin();
	     mapIntClusterClassIt != clusterClassMap[armNow].end();
	     ++mapIntClusterClassIt) {

	  ClusterClass* thisCluster = mapIntClusterClassIt->second;

	  if(thisCluster->OutsideFlag == 1)	continue;
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
      OutputManager.TreeIntV["mFlag"]        = thisCluster->Pdg;
      OutputManager.TreeIntV["highE"]        = thisCluster->HighestEnergyFlag;
      OutputManager.TreeIntV["sign"]         = armNow;
      OutputManager.TreeIntV["nHits"]        = thisCluster->NumHits;
      OutputManager.TreeDoubleV["distTheta"] = thisCluster->DiffTheta;
      OutputManager.TreeDoubleV["distXY"]    = thisCluster->DiffPosXY;
      OutputManager.TreeDoubleV["engy"]      = thisCluster->Engy;
      OutputManager.TreeDoubleV["theta"]     = thisCluster->Theta;
      OutputManager.TreeDoubleV["engyMC"]    = thisCluster->EngyMC;
      OutputManager.TreeDoubleV["thetaMC"]   = thisCluster->ThetaMC;
      OutputManager.TreeDoubleV["phi"]       = thisCluster->Phi;
      OutputManager.TreeDoubleV["rzStart"]   = thisCluster->RZStart;
      OutputManager.TreeDoubleV["Xglob"]     = globPos[0];
      OutputManager.TreeDoubleV["Yglob"]     = globPos[1];
      OutputManager.TreeDoubleV["Zglob"]     = globPos[2];
      //--
      OutputManager.FillRootTree("LumiRecoParticleTree");
    }
      }

      storeMCParticleInfo( evt, clusterInFlag );

#if _GLOBAL_COUNTERS_UPDATE_DEBUG == 1
      // write out the counter map
      int numCounters = OutputManager.Counter.size();

      if(numCounters > 0)
        streamlog_out(DEBUG1) << std::endl << "Global counters:" << std::endl;

      OutputManager.CounterIterator = OutputManager.Counter.begin();
      for(int hisNow = 0; hisNow < numCounters; hisNow++ , OutputManager.CounterIterator++) {
	std::string counterName = (std::string)(*OutputManager.CounterIterator).first;
    streamlog_out(DEBUG1) << "\t" << OutputManager.Counter[counterName] << "  \t <->  " << counterName << std::endl;
      }
#endif

     /* --------------------------------------------------------------------------
	 write to the root tree
	 -------------------------------------------------------------------------- */
      OutputManager.WriteToRootTree("" , NumEvt);

      /* --------------------------------------------------------------------------
	 clean ClusterClassMap
	 -------------------------------------------------------------------------- */
      for(int armNow = -1; armNow < 2; armNow += 2) {
	for (MapIntPClusterClass::iterator mapIntClusterClassIt = clusterClassMap[armNow].begin();
	     mapIntClusterClassIt != clusterClassMap[armNow].end();
	     ++mapIntClusterClassIt) {
	  delete mapIntClusterClassIt->second;
	  mapIntClusterClassIt->second = NULL;
	}
      }
 } // if ( _Lumi_Control_Out )

    } // try

    // if an !E!9exception has been thrown (no *col for this event) than do....
    catch( DataNotAvailableException &e ){
      streamlog_out(DEBUG7) << "Event " << NumEvt << " has an exception" << std::endl;
    }

    return;
  }

void MarlinLumiCalClusterer::CreateClusters(MapIntMapIntVInt const&    clusterIdToCellId,
                                            MapIntMapIntVDouble const& cellIdToCellEngy,
                                            std::map<int, MapIntPClusterClass>& clusterClassMap, EVENT::LCEvent* evt) {
  std::vector<MCInfo> mcParticlesVecPos;
  std::vector<MCInfo> mcParticlesVecNeg;

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
      MCInfo      p        = MCInfo::getMCParticleInfo(particle, gmc);
      if (p.mcp == NULL)
        continue;
      streamlog_out(DEBUG2) << p << std::endl;
      if (p.sign > 0)
        mcParticlesVecPos.push_back(p);
      else
        mcParticlesVecNeg.push_back(p);
    }
  }
  int numOfClustersNeg = clusterIdToCellId.at(-1).size();
  int numOfClustersPos = clusterIdToCellId.at(1).size();
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
    for (MapIntVInt::const_iterator clusterIdToCellIdIterator = clusterIdToCellId.at(armNow).begin();
         clusterIdToCellIdIterator != clusterIdToCellId.at(armNow).end(); clusterIdToCellIdIterator++) {
      const int clusterId = clusterIdToCellIdIterator->first;
      // create a new cluster and put it on map
      ClusterClass* thisCluster = clusterClassMap[armNow][clusterId] = new ClusterClass(
          gmc, clusterId, armNow, clusterIdToCellIdIterator->second, cellIdToCellEngy.at(armNow).at(clusterId));
      if (thisCluster->Engy > EngyMax) {
        EngyMax   = thisCluster->Engy;
        EngyMaxID = clusterId;
      }

      LCCluster const& thisClusterInfo = LumiCalClusterer._superClusterIdClusterInfo[armNow][clusterId];
      // clang-format off
      streamlog_out(DEBUG6) << "arm =   " << armNow << "\t cluster " << clusterId << "  ...... " << std::endl
                            << std::setw(20) << "X, Y, Z:" << std::endl
                            << std::setw(20) << "ClusterClass" << std::fixed << std::setprecision(3)
                            << std::setw(13) << thisCluster->clusterPosition[0]
                            << std::setw(13) << thisCluster->clusterPosition[1]
                            << std::setw(13) << thisCluster->clusterPosition[2] << std::endl
                            << std::setw(20) << "ClusterInfo"
                            << std::setw(13) << thisClusterInfo.getX()
                            << std::setw(13) << thisClusterInfo.getY()
                            << std::setw(13) << thisClusterInfo.getZ() << std::endl
                            << std::setw(20) << "Energy, Theta, Phi: " << std::endl
                            << std::setw(20) << "ClusterClass"
                            << std::setw(13) << thisCluster->Engy
                            << std::setw(13) << thisCluster->Theta
                            << std::setw(13) << thisCluster->Phi << std::endl
                            << std::setw(20) << "ClusterInfo"
                            << std::setw(13) << thisClusterInfo.getEnergy()
                            << std::setw(13) << thisClusterInfo.getTheta()
                            << std::setw(13) << thisClusterInfo.getPhi()
                            << std::endl
                            << std::endl;
      // clang-format on
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
    std::vector<MCInfo>& mcParticlesVec = (armNow == -1) ? mcParticlesVecNeg : mcParticlesVecPos;

    for (MapIntPClusterClass::iterator mapIntClusterClassIt = clusterClassMap[armNow].begin();
         mapIntClusterClassIt != clusterClassMap[armNow].end(); ++mapIntClusterClassIt) {
      const int     clusterId   = mapIntClusterClassIt->first;
      ClusterClass* thisCluster = mapIntClusterClassIt->second;

      double RZStart = thisCluster->RZStart;
      double xs      = RZStart * cos(thisCluster->Phi);
      double ys      = RZStart * sin(thisCluster->Phi);

      using namespace std::placeholders;  //_1, _2
      if (mcParticlesVec.size()) {
        // try to match MC true particle with cluster by comparing positions at Lumical entry
        // sort(mcParticlesVec.begin(), mcParticlesVec.end(),
        //      std::bind(BindThetaCompAsc, thisCluster->Theta, _1, _2));
        sort(mcParticlesVec.begin(), mcParticlesVec.end(), std::bind(PositionCompAsc, xs, ys, _1, _2));
        streamlog_out(DEBUG4) << "Trying to match particle: " << mcParticlesVec[0].engy << std::endl;
        double dTheta = fabs(thisCluster->Theta - mcParticlesVec[0].theta);
        double dPos0  = sqrt((sqr(xs - mcParticlesVec[0].x) + sqr(ys - mcParticlesVec[0].y)));
        streamlog_out(DEBUG4) << "RZStart " << RZStart << std::endl;
        streamlog_out(DEBUG4) << "  xs, ys " << std::setw(13) << xs << std::setw(13) << ys << std::endl
                              << "MCP x, y " << std::setw(13) << mcParticlesVec[0].x << std::setw(13)
                              << mcParticlesVec[0].y << std::endl;

        double dEne0 = fabs(mcParticlesVec[0].engy - thisCluster->Engy);
        streamlog_out(DEBUG4) << " dTheta " << std::setw(13) << dTheta << " dPos0 " << std::setw(13) << dPos0 << " dEne0 "
                              << std::setw(13) << dEne0 << std::endl;

        thisCluster->DiffTheta = dTheta;
        if (mcParticlesVec.size() > 1) {
          double dPos1 = sqrt((sqr(xs - mcParticlesVec[1].x) + sqr(ys - mcParticlesVec[1].y)));
          double dEne1 = fabs(mcParticlesVec[1].engy - thisCluster->Engy);
          if (dPos1 < gmc.GlobalParamD[GlobalMethodsClass::MinSeparationDist] && dEne1 < dEne0) {
            thisCluster->SetStatsMC(mcParticlesVec[1].mcp);
            thisCluster->DiffPosXY = dPos1;
            mcParticlesVec.erase(mcParticlesVec.begin() + 1);
          } else {
            thisCluster->SetStatsMC(mcParticlesVec[0].mcp);
            thisCluster->DiffPosXY = dPos0;
            mcParticlesVec.erase(mcParticlesVec.begin());
          }
        } else {
          thisCluster->DiffPosXY = dPos0;
          if (dPos0 < gmc.GlobalParamD[GlobalMethodsClass::MinSeparationDist]) {
            thisCluster->SetStatsMC(mcParticlesVec[0].mcp);
            mcParticlesVec.erase(mcParticlesVec.begin());
          }
        }
      }

      if (thisCluster->Pdg != 0) {
        double Ereco = thisCluster->Engy;
        // clang-format off
        streamlog_out( DEBUG6 ) << "\tParticle Out ("
                                << thisCluster -> OutsideReason << "):   " << clusterId << std::endl
                                << "\t\t side(arm), pdg, parentId , NumMCDaughters = "
                                << "\t" << thisCluster -> SignMC<<"("<<armNow<<")"
                                << "\t" << thisCluster -> Pdg
                                << "\t" << thisCluster -> ParentId
                                << "\t" << thisCluster -> NumMCDaughters << std::endl
                                << "\t\t MCPos    X, Y, Z: "
                                << "\t" << std::setw(13) << thisCluster -> mcpPosition[0]
                                << "\t" << std::setw(13) << thisCluster -> mcpPosition[1]
                                << "\t" << std::setw(13) << thisCluster -> mcpPosition[2]<< std::endl
                                << "\t\t Cluster  X, Y, Z: "
                                << "\t" << std::setw(13) << thisCluster -> clusterPosition[0]
                                << "\t" << std::setw(13) << thisCluster -> clusterPosition[1]
                                << "\t" << std::setw(13) << thisCluster -> clusterPosition[2] << std::endl
                                << "\t\t engy, theta, phi (mc): "
                                <<  std::setw(13)<< thisCluster -> EngyMC
                                << "\t"<< std::setw(13)<<thisCluster -> ThetaMC
                                << "\t"<< std::setw(13)<<thisCluster -> PhiMC << std::endl
                                << "\t\t engy, theta, phi (rec):"
                                << std::setw(13)<<Ereco
                                << "\t"<< std::setw(13)<<thisCluster -> Theta
                                << "\t"<< std::setw(13)<<thisCluster -> Phi
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
    LCCollection * particles(NULL);
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
      MCInfo p = MCInfo::getMCParticleInfo( particle, gmc );
      if( p.mcp == NULL ) continue;
      const double *endPoint = particle->getEndpoint();
      const double *vx = particle->getVertex();

      OutputManager.TreeIntV["nEvt"]	= NumEvt;
      OutputManager.TreeIntV["sign"]	= p.sign;
      OutputManager.TreeIntV["pdg"]	= p.pdg;
      OutputManager.TreeDoubleV["engy"]	= p.engy;
      OutputManager.TreeDoubleV["theta"]= p.theta;
      OutputManager.TreeDoubleV["phi"]	= p.phi;
      // position at the face of LCal (global coordinates)
      OutputManager.TreeDoubleV["begX"]	= p.x;
      OutputManager.TreeDoubleV["begY"]	= p.y;
      OutputManager.TreeDoubleV["begZ"]	= double(p.sign)*LcalZstart;
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

std::tuple<ClusterImpl*, ReconstructedParticleImpl*> MarlinLumiCalClusterer::getLCIOObjects(
    LCCluster const& thisClusterInfo) const {
  double ThetaMid =
      (gmc.GlobalParamD.at(GlobalMethodsClass::ThetaMin) + gmc.GlobalParamD.at(GlobalMethodsClass::ThetaMax)) / 2.;
  double ThetaTol =
      (gmc.GlobalParamD.at(GlobalMethodsClass::ThetaMax) - gmc.GlobalParamD.at(GlobalMethodsClass::ThetaMin)) / 2.;

  const double clusterEnergy = thisClusterInfo.getE();
  if (clusterEnergy < _minClusterEngy)
    return std::make_tuple(nullptr, nullptr);

  if (_cutOnFiducialVolume) {
    const double clusterTheta = thisClusterInfo.getTheta();
    if (fabs(clusterTheta - ThetaMid) > ThetaTol)
      return std::make_tuple(nullptr, nullptr);
  }

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
