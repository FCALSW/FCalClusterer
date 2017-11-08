
#include "MarlinLumiCalClusterer.h"

#include "ClusterClass.h"
#include "MCInfo.h"

#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/ClusterImpl.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>

#include <algorithm>
#include <map>
#include <vector>
#include <iomanip>


/* >> */ 

/* -----------------------------------------------------------------------------------------*/
// (BP) defs for sorting

double _sortAgainstThisTheta;
double _sortAgainstThisX;
double _sortAgainstThisY;

// variable _sortAgainstThisTheta must be set before sorting 
bool ThetaCompAsc( MCInfo a, MCInfo b ){
  return fabs( _sortAgainstThisTheta - a.theta) < fabs( _sortAgainstThisTheta -b.theta) ;
}

bool PositionCompAsc( MCInfo a, MCInfo b ){
  double dxa = ( _sortAgainstThisX - a.x);
  double dya = ( _sortAgainstThisY - a.y);
  double dxb = ( _sortAgainstThisX - b.x);
  double dyb = ( _sortAgainstThisY - b.y);
  return (dxa*dxa+dya*dya) < (dxb*dxb+dyb*dyb)  ;
}
// LumiCal rotations angles ( local->Global )
std::map < int , double > csbx;
std::map < int , double > snbx;
/*------------------------------------------------------------------------------------------*/
  void MarlinLumiCalClusterer::TryMarlinLumiCalClusterer(  EVENT::LCEvent * evt  ){

    try{


      double ThetaMid = (gmc.GlobalParamD[GlobalMethodsClass::ThetaMin] + gmc.GlobalParamD[GlobalMethodsClass::ThetaMax])/2.;
      double ThetaTol = (gmc.GlobalParamD[GlobalMethodsClass::ThetaMax] - gmc.GlobalParamD[GlobalMethodsClass::ThetaMin])/2.;


      csbx[-1] = cos( - gmc.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle]/2.);
      csbx[ 1] = cos( gmc.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle]/2.);
      snbx[-1] = sin( - gmc.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle]/2.);
      snbx[ 1] = sin( gmc.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle]/2.);
      /* --------------------------------------------------------------------------
	 create clusters using: LumiCalClustererClass
	 -------------------------------------------------------------------------- */

      if ( !LumiCalClusterer.processEvent(evt) ) return;

      LCCollectionVec* LCalClusterCol = new LCCollectionVec(LCIO::CLUSTER);
      LCCollectionVec* LCalRPCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

#if _CREATE_CLUSTERS_DEBUG == 1
      streamlog_out(DEBUG2) << " Transfering reco results to LCalClusterCollection....."<<std::endl;
#endif 
  
      for(int armNow = -1; armNow < 2; armNow += 2) {

        streamlog_out(DEBUG2)<<" Arm  "<< std::setw(4)<< armNow
			     << "\t Number of clusters: "<< LumiCalClusterer._superClusterIdToCellId[armNow].size()
			     <<std::endl;

	for( MapIntVInt::const_iterator clusterIdToCellIdIterator = LumiCalClusterer._superClusterIdToCellId[armNow].begin();
	     clusterIdToCellIdIterator != LumiCalClusterer._superClusterIdToCellId[armNow].end();
	     clusterIdToCellIdIterator++) {
	  const int clusterId = clusterIdToCellIdIterator->first;

	  LCCluster const& thisClusterInfo = LumiCalClusterer._superClusterIdClusterInfo[armNow][clusterId];

	  const double clusterEnergy = gmc.SignalGevConversion(GlobalMethodsClass::Signal_to_GeV , thisClusterInfo.getE());
	  if( clusterEnergy < _minClusterEngy ) continue;

          if( _cutOnFiducialVolume ) {
            const double clusterTheta = thisClusterInfo.getTheta();
            if( fabs ( clusterTheta - ThetaMid ) >  ThetaTol ) continue;
          }
 
	  const float  xloc =  float(thisClusterInfo.getX());
          const float  yloc =  float(thisClusterInfo.getY());
          const float  zloc =  float(thisClusterInfo.getZ());

	  ClusterImpl* cluster = new ClusterImpl;
	  cluster->setEnergy( clusterEnergy );

	  ReconstructedParticleImpl* particle = new ReconstructedParticleImpl;
	  const float mass = 0.0;
	  const float charge = 1e+19;
	  particle->setMass( mass ) ;
	  particle->setCharge( charge ) ;
	  particle->setEnergy ( clusterEnergy ) ;
	  particle->addCluster( cluster ) ;

	  const float gP[] = { float(  csbx[ armNow ]*xloc + snbx[ armNow ]*zloc ),
			       float(  yloc ),
			       float( -snbx[ armNow ]*xloc + csbx[ armNow ]*zloc )};
	  cluster->setPosition( gP );

	  const float norm = sqrt( gP[0]*gP[0] + gP[1]*gP[1] + gP[2]*gP[2] );
	  const float clusterMomentum[3] =  { float(gP[0]/norm * clusterEnergy),
					      float(gP[1]/norm * clusterEnergy),
					      float(gP[2]/norm * clusterEnergy) };
	  particle->setMomentum ( clusterMomentum) ;

#pragma message ("FIXME: Link Calohits to the Cluster")
	  //Get the cellID from the container Hit of the cluster and find the
	  //matching LumiCal SimCalorimeterHit in the event collection, probably
	  //have to loop over all of them until it is found

	  LCalClusterCol->addElement(cluster);
	  LCalRPCol->addElement(particle);

	}

     }	

      //Add collections to the event if there are clusters
      if ( LCalClusterCol->getNumberOfElements() != 0 ) {
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
	  const double engyNow  = gmc.SignalGevConversion(GlobalMethodsClass::Signal_to_GeV,
							  thisCluster->Engy);
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
	      //	    } else {
	      // EngyMC/ThetaMC is nowhere set (?) (BP)
	      // engyNow  = thisCluster -> EngyMC;
	      // thetaNow = thisCluster -> ThetaMC;
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

	  const double engyNow = gmc.SignalGevConversion(GlobalMethodsClass::Signal_to_GeV,
							 thisCluster -> Engy);
	  const double thetaNow = thisCluster -> Theta;
	  const double phiNow   = thisCluster -> Phi;
	  const double rzstartNow = thisCluster -> RZStart;

	  //(BP) compute x,y,z in global reference system
	  double xloc = rzstartNow * cos(phiNow);
	  double yloc = rzstartNow * sin(phiNow);
	  double zloc = rzstartNow / tan ( thetaNow );
	  double xglob = xloc*csbx[armNow] + zloc*snbx[armNow];
	  double yglob = yloc;
	  double zglob =-xloc*snbx[armNow] + zloc*csbx[armNow];
	  //	  thetaNow = atan( sqrt( sqr(xglob) + sqr(yglob) )/fabs(zglob) );
	  int mFlag = thisCluster->Pdg;
	  int nHits = thisCluster->NumHits;
	  int highE = thisCluster->HighestEnergyFlag;
	  //-------------
	  OutputManager.TreeIntV["nEvt"]	= NumEvt;
	  OutputManager.TreeIntV["outFlag"]	= thisCluster->OutsideFlag;
	  OutputManager.TreeIntV["mFlag"]	= mFlag;
	  OutputManager.TreeIntV["highE"]	= highE;
	  OutputManager.TreeIntV["sign"]	= armNow;
	  OutputManager.TreeIntV["nHits"]	= nHits;
	  OutputManager.TreeDoubleV["distTheta"]= thisCluster->DiffTheta;
	  OutputManager.TreeDoubleV["distXY"]   = thisCluster->DiffPosXY;
	  OutputManager.TreeDoubleV["engy"]	= engyNow;
	  OutputManager.TreeDoubleV["theta"]	= thetaNow;
	  OutputManager.TreeDoubleV["engyMC"]	= thisCluster->EngyMC;
	  OutputManager.TreeDoubleV["thetaMC"]	= thisCluster->ThetaMC;
	  OutputManager.TreeDoubleV["phi"]	= phiNow;
	  OutputManager.TreeDoubleV["rzStart"]	= rzstartNow;
	  // position at Zstart
	  OutputManager.TreeDoubleV["Xglob"]    = xglob;
	  OutputManager.TreeDoubleV["Yglob"]    = yglob;
	  OutputManager.TreeDoubleV["Zglob"]    = zglob;
	  //--
	  OutputManager.FillRootTree("LumiRecoParticleTree");
	}
      }

      storeMCParticleInfo( evt, clusterInFlag );

#if _GLOBAL_COUNTERS_UPDATE_DEBUG == 1
      // write out the counter map
      int numCounters = OutputManager.Counter.size();

      if(numCounters > 0)
	streamlog_out(DEBUG) << std::endl << "Global counters:"  << std::endl;

      OutputManager.CounterIterator = OutputManager.Counter.begin();
      for(int hisNow = 0; hisNow < numCounters; hisNow++ , OutputManager.CounterIterator++) {
	std::string counterName = (std::string)(*OutputManager.CounterIterator).first;
	streamlog_out(DEBUG) << "\t" << OutputManager.Counter[counterName] << "  \t <->  " << counterName << std::endl;
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
#ifdef _LC_DEBUG
      streamlog_out(DEBUG) << "Event " << NumEvt << " has an exception"<< std::endl;
#endif
    }

    return;
  }




  void MarlinLumiCalClusterer::CreateClusters(	std::map < int , std::map < int , std::vector<int> > > const& clusterIdToCellId,
						std::map < int , std::map < int , std::vector<double> > > const& cellIdToCellEngy,
						std::map < int, MapIntPClusterClass > & clusterClassMap,
						EVENT::LCEvent * evt ) {


    std::vector < MCInfo > mcParticlesVecPos;
    std::vector < MCInfo > mcParticlesVecNeg;

#if _CREATE_CLUSTERS_DEBUG == 1
	streamlog_out(MESSAGE4) << "CreateClusters: event = " << evt-> getEventNumber() << std::endl;
#endif
    LCCollection * particles = evt->getCollection( "MCParticle" );
    if ( particles ){
      int numMCParticles = particles->getNumberOfElements();
#if _CREATE_CLUSTERS_DEBUG == 1
      if( numMCParticles ){
	streamlog_out(MESSAGE4) << "CreateClusters: numMCParticles = " << numMCParticles  << std::endl;
      }else{
	streamlog_out(MESSAGE4) << "CreateClusters: No MCparticles in this event !" << std::endl;
      }
#endif
      for ( int jparticle=0; jparticle<numMCParticles; jparticle++ ){
	MCParticle * particle = static_cast<MCParticle*>( particles->getElementAt(jparticle) );
	MCInfo p = MCInfo::getMCParticleInfo( particle, gmc );
	if( p.mcp == NULL ) continue;
	streamlog_out(DEBUG2) << p << std::endl;
	if ( p.sign > 0 ) mcParticlesVecPos.push_back( p );
        else mcParticlesVecNeg.push_back( p );

      }
    }
    int numOfClustersNeg =  clusterIdToCellId.at(-1).size();
    int numOfClustersPos =  clusterIdToCellId.at(1).size();
#if _CREATE_CLUSTERS_DEBUG == 1
    if( numOfClustersNeg || numOfClustersPos ){
      streamlog_out(MESSAGE) << "Initial Set Stats for LumiCal......."<< std::endl;
      streamlog_out(MESSAGE) << "     numOfClusters arm[-1] = "<< numOfClustersNeg << "\t arm[1] = "<< numOfClustersPos << std::endl;
    }else{
      streamlog_out(MESSAGE) << "Initial Set Stats for LumiCal: No clusters found in this event !"<< std::endl;
    }
#endif
    if( numOfClustersNeg || numOfClustersPos ){ 
      for(int armNow = -1; armNow < 2; armNow += 2) {
	double EngyMax =  0.;
	int  EngyMaxID = -1 ;
	for( MapIntVInt::const_iterator clusterIdToCellIdIterator = clusterIdToCellId.at(armNow).begin();
	     clusterIdToCellIdIterator != clusterIdToCellId.at(armNow).end();
	       clusterIdToCellIdIterator++){
	  const int clusterId = clusterIdToCellIdIterator->first;
	  // create a new cluster and put it on map
	  clusterClassMap[armNow][clusterId] = new ClusterClass(clusterId, gmc);
	  ClusterClass* thisCluster = clusterClassMap[armNow][clusterId];
	  thisCluster->SignMC = armNow;

	  double engySum = 0.;
	  int cellNow=0;
	  for(std::vector<int>::const_iterator cellIt = clusterIdToCellIdIterator->second.begin();
	      cellIt != clusterIdToCellIdIterator->second.end();
	      ++cellIt, ++cellNow) {
	    const int cellId = *cellIt;
	    double engyHit = cellIdToCellEngy.at(armNow).at(clusterId).at(cellNow);
	    thisCluster->FillHit(cellId , engyHit);
	    engySum += engyHit;
	  }
	  if( engySum > EngyMax ) {
	    EngyMax   = engySum;
	    EngyMaxID = clusterId;
	  }
	
	  thisCluster->ResetStats(); // calculate energy, position for the cluster
	
	  LCCluster const& thisClusterInfo = LumiCalClusterer._superClusterIdClusterInfo[armNow][clusterId];
	  streamlog_out(DEBUG3) << "arm =   " << armNow <<"\t cluster "<< clusterId<< "  ...... " << std::endl
				<< std::setw(20) << "X, Y, Z:" << std::endl
				<< std::setw(20) << "ClusterClass"
				<< std::fixed << std::setprecision(3)
				<< std::setw(13) << thisCluster-> clusterPosition[0]
				<< std::setw(13) << thisCluster-> clusterPosition[1]
				<< std::setw(13) << thisCluster-> clusterPosition[2] << std::endl
				<< std::setw(20) << "ClusterInfo"
				<< std::setw(13) << thisClusterInfo.getX()
				<< std::setw(13) << thisClusterInfo.getY()
				<< std::setw(13) << thisClusterInfo.getZ() << std::endl
				<< std::setw(20) << "Energy, Theta, Phi: " << std::endl
				<< std::setw(20) << "ClusterClass"
				<< std::setw(13) << thisCluster-> Engy
				<< std::setw(13) << thisCluster-> Theta
				<< std::setw(13) << thisCluster-> Phi  << std::endl
				<< std::setw(20) << "ClusterInfo"
				<< std::setw(13) << thisClusterInfo.getEnergy()
				<< std::setw(13) << thisClusterInfo.getTheta()
				<< std::setw(13) << thisClusterInfo.getPhi()  << std::endl
				<< std::endl;

	}
	// set flag for highest energy cluster found in this event
	if ( EngyMax > 0.0 ) clusterClassMap[armNow][EngyMaxID]->HighestEnergyFlag = 1;
      }



#if _CREATE_CLUSTERS_DEBUG == 1
      if( mcParticlesVecNeg.size() || mcParticlesVecPos.size() ) {
	streamlog_out(MESSAGE4) << "Transfering MC information into ClusterClass objects ......  "  << std::endl;
	streamlog_out(MESSAGE4) << "MC particles arm z-: "<< mcParticlesVecNeg.size() << "\t arm z+: "<< mcParticlesVecPos.size() << std::endl;
      }else{
	streamlog_out(MESSAGE4) << "Transfering MC information into ClusterClass objects: No primary MCParticles"
				<< " entering LumiCal found!"  << std::endl;
      }
#endif
    
      if( mcParticlesVecNeg.size() || mcParticlesVecPos.size() ) {
	for(int armNow = -1; armNow < 2; armNow += 2) {

	  std::vector < MCInfo > mcParticlesVec;
	  if ( armNow == -1 ) mcParticlesVec = mcParticlesVecNeg ;
	  else                mcParticlesVec = mcParticlesVecPos ;

	  for (MapIntPClusterClass::iterator mapIntClusterClassIt = clusterClassMap[armNow].begin();
	       mapIntClusterClassIt != clusterClassMap[armNow].end();
	       ++mapIntClusterClassIt) {
	    const int clusterId = mapIntClusterClassIt->first;
	    ClusterClass* thisCluster = mapIntClusterClassIt->second;

	    double eneCL = gmc.SignalGevConversion(GlobalMethodsClass::Signal_to_GeV ,thisCluster -> Engy);
	    double phiCL = thisCluster -> Phi;
	    double RZStart = thisCluster -> RZStart;
	    double xs = RZStart*cos(phiCL);
	    double ys = RZStart*sin(phiCL);
	    _sortAgainstThisX = xs;
	    _sortAgainstThisY = ys;
	    _sortAgainstThisTheta = thisCluster -> Theta;

	    if( mcParticlesVec.size() ){
	      // try to match MC true particle with cluster by comparing positions at Lumical entry
	      //	  sort( mcParticlesVec.begin(), mcParticlesVec.end(), ThetaCompAsc );
	      sort( mcParticlesVec.begin(), mcParticlesVec.end(), PositionCompAsc );
	      streamlog_out(DEBUG2) << "Trying to match particle: " << mcParticlesVec[0].engy << std::endl;
	      double dTheta  = fabs( _sortAgainstThisTheta - mcParticlesVec[0].theta );
	      double dPos0    = sqrt( ( sqr(xs - mcParticlesVec[0].x) + sqr( ys - mcParticlesVec[0].y)));
	      streamlog_out(DEBUG4) << "RZStart " << RZStart  << std::endl;
	      streamlog_out(DEBUG4) << "  xs, ys "
				    << std::setw(13) << xs
				    << std::setw(13) << ys  << std::endl
				    << "MCP x, y "
				    << std::setw(13) << mcParticlesVec[0].x
				    << std::setw(13) << mcParticlesVec[0].y
				    << std::endl;

	      double dEne0    = fabs( mcParticlesVec[0].engy - eneCL);
	      streamlog_out(DEBUG4) << " dTheta " << std::setw(13) << dTheta
				    << " dPos0 "  << std::setw(13) << dPos0
				    << " dEne0 "  << std::setw(13) << dEne0
				    << std::endl;

	      thisCluster -> DiffTheta = dTheta;
	      if( mcParticlesVec.size() > 1 ) {
		double dPos1    = sqrt( ( sqr(xs - mcParticlesVec[1].x) + sqr( ys - mcParticlesVec[1].y)));
                double dEne1    = fabs( mcParticlesVec[1].engy - eneCL);
		if( dPos1 <  gmc.GlobalParamD[GlobalMethodsClass::MinSeparationDist] && dEne1 < dEne0 ) {
		  thisCluster->SetStatsMC( mcParticlesVec[1].mcp );
		  thisCluster -> DiffPosXY = dPos1;
		  mcParticlesVec.erase( mcParticlesVec.begin()+1 ); 
		}else{
		  thisCluster->SetStatsMC( mcParticlesVec[0].mcp );
		  thisCluster -> DiffPosXY = dPos0;
		  mcParticlesVec.erase( mcParticlesVec.begin() ); 
		}
	      }else{
		thisCluster -> DiffPosXY = dPos0;
		if( dPos0 <  gmc.GlobalParamD[GlobalMethodsClass::MinSeparationDist] ){
		  thisCluster->SetStatsMC( mcParticlesVec[0].mcp );
		  mcParticlesVec.erase( mcParticlesVec.begin() ); 
		}
	      }
	    }

	    if( thisCluster -> Pdg != 0) {
	      double Ereco = gmc.SignalGevConversion(GlobalMethodsClass::Signal_to_GeV ,thisCluster -> Engy);
	      streamlog_out( MESSAGE4 ) << "\tParticle Out ("
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
	    }else{
	      streamlog_out( MESSAGE4 ) << "Arm: "<< armNow << "\t Cluster: " << clusterId << " does not match MC particle !\n";
	      thisCluster->PrintInfo();
	    }

	  } // loop over clusters
	} // loop over arms
      } // if there are mcparticles
    } // if there are clusters
    return;

  }

void MarlinLumiCalClusterer::storeMCParticleInfo( LCEvent *evt, int clusterInFlag ) {
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
