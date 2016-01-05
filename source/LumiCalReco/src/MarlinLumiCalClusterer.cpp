
#include "MarlinLumiCalClusterer.h"

#include "SortingClass.h"
#include "ClusterClass.h"

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


/* >> */ 

  /* --------------------------------------------------------------------------
     -------------------------------------------------------------------------- */
// (BP) quick sorting defs
double _sortAgainstThisTheta;
struct myMCP{
  MCParticle * mcp;
  double theta;
};

bool ThetaCompAsc( myMCP a, myMCP b ){
  return fabs( _sortAgainstThisTheta - a.theta) < fabs( _sortAgainstThisTheta -b.theta) ;
}

  void MarlinLumiCalClusterer::TryMarlinLumiCalClusterer(  EVENT::LCEvent * evt  ){

    try{

      std::string hisName, optionName;
      int numMCParticles, particleId;
      double	engyNow, thetaNow, phiNow, rzstartNow;

      std::map < int , std::map < int , ClusterClass * > >	clusterClassMap;
      std::map < int , ClusterClass * > :: iterator	clusterClassMapIterator;
      std::map < int , double > csbx;
      std::map < int , double > snbx;

      int	numClusters, clusterId;
      double	weightNow;

      csbx[-1] = cos( M_PI - GlobalMethods.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle]/2.);
      csbx[ 1] = cos( GlobalMethods.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle]/2.);
      snbx[-1] = sin( M_PI - GlobalMethods.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle]/2.);
      snbx[ 1] = sin( GlobalMethods.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle]/2.);
      /* --------------------------------------------------------------------------
	 create clusters using: LumiCalClustererClass
	 -------------------------------------------------------------------------- */
      // instantiate a clusterClass object for each mcParticle which was created
      // infront of LumiCal and was destroyed after lumical.

      if ( !LumiCalClusterer.processEvent(evt) ) return;


      CreateClusters( LumiCalClusterer._superClusterIdToCellId,
		      LumiCalClusterer._superClusterIdToCellEngy,
		      clusterClassMap,
		      evt );


      /* --------------------------------------------------------------------------
	 flag the highest-energy cluster in each arm
	 -------------------------------------------------------------------------- */
      
      LCCollectionVec* LCalClusterCol = new LCCollectionVec(LCIO::CLUSTER);
      LCCollectionVec* LCalRPCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

      for(int armNow = -1; armNow < 2; armNow += 2) {
	//	std::vector < SortingClass > sortingClassV;

	numClusters             = clusterClassMap[armNow].size();
	clusterClassMapIterator = clusterClassMap[armNow].begin();
	int highestEnergyClusterId = -1;
	double theHighestEnergy = 0.;
	for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterClassMapIterator++) {
	  clusterId = (int)(*clusterClassMapIterator).first;

	  ClusterClass const* thisCluster = clusterClassMap[armNow][clusterId];
	  LCCluster const& thisClusterInfo = LumiCalClusterer._superClusterIdClusterInfo[armNow][clusterId];

	  ClusterImpl* cluster = new ClusterImpl;
	  const double thetaCluster = thisCluster->Theta;
	  const double energyCluster = GlobalMethodsClass::SignalGevConversion(GlobalMethodsClass::Signal_to_GeV , thisCluster->Engy);
	  const double phiCluster = thisCluster->Phi;
	  cluster->setEnergy( energyCluster );

	  const float clusterPosition[3] = { float(thisClusterInfo.getPosition()[0]),
					     float(thisClusterInfo.getPosition()[1]),
					     float(thisClusterInfo.getPosition()[2])};

	  cluster->setPosition( clusterPosition );
#pragma message ("FIXME: Link Calohits to the Cluster")
	  //Get the cellID from the container Hit of the cluster and find the
	  //matching LumiCal SimCalorimeterHit in the event collection, probably
	  //have to loop over all of them until it is found
	  float momentumCluster[3] = { float(energyCluster * sin ( thetaCluster ) * cos ( phiCluster ) * float(armNow) ),
				       float(energyCluster * sin ( thetaCluster ) * sin ( phiCluster )),
				       float(energyCluster * cos ( thetaCluster ) * float(armNow))  };

	  const float mass = 0.0;
	  const float charge = 1e+19;

	  ReconstructedParticleImpl* particle = new ReconstructedParticleImpl;
	  particle->setMass( mass ) ;
	  particle->setCharge( charge ) ;
	  particle->setMomentum ( momentumCluster ) ;
	  particle->setEnergy ( energyCluster ) ;
	  particle->addCluster( cluster ) ;

	  LCalClusterCol->addElement(cluster);
	  LCalRPCol->addElement(particle);

	  if(clusterClassMap[armNow][clusterId]->OutsideFlag == 1)	continue;

	  weightNow   = clusterClassMap[armNow][clusterId]->Engy;
	  /* (BP) make it simpler
	  weightNow   = 1./weightNow;
	  // instantiate the sorting class with the weight and insert to the vector
	  sortingClassV.push_back(SortingClass(clusterId , weightNow));
	  */
	  if( weightNow > theHighestEnergy ){
	    theHighestEnergy = weightNow;
	    highestEnergyClusterId = clusterId;
	  }
	}

	// flag the highest-energy cluster
	if ( highestEnergyClusterId+1 > 0 )  clusterClassMap[armNow][highestEnergyClusterId]->HighestEnergyFlag = 1;

	/* (BP) avoid wasting time
        arrange the vector according to inverse of the energy (lowest first)
	sort(sortingClassV.begin(),sortingClassV.end(),cmpRuleDesc);
	
	-- flag the highest-energy cluster
	numClusters = sortingClassV.size();
	for(int clusterNow = 0; clusterNow < numClusters; clusterNow++) {
	  clusterId = sortingClassV[clusterNow].Id;

	  if(clusterNow == 0)
	    clusterClassMap[armNow][clusterId]->HighestEnergyFlag = 1;
	  else
	    clusterClassMap[armNow][clusterId]->HighestEnergyFlag = 0;
	}
	*/

      }
	

      //Add collections to the event if there are clusters
      if ( LCalClusterCol->getNumberOfElements() != 0 ) {
	evt->addCollection(LCalClusterCol, LumiClusterColName);
	evt->addCollection(LCalRPCol, LumiRecoParticleColName);
      } else {
	delete LCalClusterCol;
	delete LCalRPCol;
      }




      /* --------------------------------------------------------------------------
	 histograming
	 -------------------------------------------------------------------------- */
      for(int armNow = -1; armNow < 2; armNow += 2) {
	double totEngyIn = 0.;
	double totEngyOut= 0.;

	clusterClassMapIterator = clusterClassMap[armNow].begin();
	numMCParticles          = clusterClassMap[armNow].size();
	for (int MCParticleNow = 0; MCParticleNow < numMCParticles; MCParticleNow++, clusterClassMapIterator++) {
	  particleId = (int)(*clusterClassMapIterator).first;

	  engyNow  = clusterClassMap[armNow][particleId] -> Engy;
	  engyNow  = GlobalMethods.SignalGevConversion(GlobalMethodsClass::Signal_to_GeV , engyNow);
	  thetaNow = clusterClassMap[armNow][particleId] -> Theta;
	  // highest energy RecoParticle
	  if(clusterClassMap[armNow][particleId]->HighestEnergyFlag == 1){
	    hisName = "higestEngyParticle_Engy";	OutputManager.HisMap1D[hisName] -> Fill (engyNow);
	    hisName = "higestEngyParticle_Theta";	OutputManager.HisMap1D[hisName] -> Fill (thetaNow);
	  }
	  if(clusterClassMap[armNow][particleId]->OutsideFlag == 1){
	    // beyond acceptance ( energy and/or fid. vol.
	    bool reason = (clusterClassMap[armNow][particleId]->OutsideReason == "Reconstructed outside the fiducial volume" );
	         reason = ( reason || ( clusterClassMap[armNow][particleId]->OutsideReason == "Cluster energy below minimum" ) );

	    if( reason ) {
	      hisName = "thetaEnergyOut_DepositedEngy";
	      OutputManager.HisMap2D[hisName] -> Fill (engyNow , thetaNow);
	      //	    } else {
	      // EngyMC/ThetaMC is nowhere set (?) (BP)
	      // engyNow  = clusterClassMap[armNow][particleId] -> EngyMC;
	      // thetaNow = clusterClassMap[armNow][particleId] -> ThetaMC;
	    }
	    totEngyOut += engyNow;
	  }else{ 
	    // accepted
	    totEngyIn += engyNow;
	  }
	}
	// total energy inside LumiCal
	 if(totEngyIn > 0) {
	  hisName = "totEnergyIn";	OutputManager.HisMap1D[hisName] -> Fill (totEngyIn);
	 }
	 if(totEngyOut > 0) {
	  hisName = "totEnergyOut";	OutputManager.HisMap1D[hisName] -> Fill (totEngyOut);
	 }
      }


      /* --------------------------------------------------------------------------
	 write criteria for selection cuts for Bhabha events into a tree
	 -------------------------------------------------------------------------- */
      int	clusterInFlag = 0;
      for(int armNow = -1; armNow < 2; armNow += 2) {
	clusterClassMapIterator = clusterClassMap[armNow].begin();
	numMCParticles          = clusterClassMap[armNow].size();
	for (int MCParticleNow = 0; MCParticleNow < numMCParticles; MCParticleNow++, clusterClassMapIterator++) {
	  particleId = (int)(*clusterClassMapIterator).first;

       	  if(clusterClassMap[armNow][particleId]->OutsideFlag == 1)	continue;
	  // only take into account the highest-energy particle in the arm
	  // (by default this means that this particle's shower is contained)
	  //	  if(clusterClassMap[armNow][particleId]->HighestEnergyFlag == 0)	continue;
	  clusterInFlag++;

	  engyNow  = clusterClassMap[armNow][particleId] -> Engy;
	  engyNow  = GlobalMethods.SignalGevConversion(GlobalMethodsClass::Signal_to_GeV , engyNow);
	  thetaNow = clusterClassMap[armNow][particleId] -> Theta;
	  phiNow   = clusterClassMap[armNow][particleId] -> Phi;
	  rzstartNow = clusterClassMap[armNow][particleId] -> RZStart;

	  //(BP) compute x,y,z in global reference system
	  double xloc = rzstartNow * cos(phiNow);
	  double yloc = rzstartNow * sin(phiNow);
	  double zloc = rzstartNow / tan ( thetaNow );
	  double xglob = xloc*csbx[armNow] + zloc*snbx[armNow];
	  double yglob = yloc;
	  double zglob =-xloc*snbx[armNow] + zloc*csbx[armNow];
	  //-------------
	  OutputManager.TreeIntV["nEvt"]	= NumEvt;
	  OutputManager.TreeIntV["outFlag"]	= clusterClassMap[armNow][particleId]->OutsideFlag;
	  OutputManager.TreeIntV["mcId"]	= clusterClassMap[armNow][particleId]->Pdg;
	  OutputManager.TreeIntV["sign"]	= armNow;
	  OutputManager.TreeDoubleV["engy"]	= engyNow;
	  OutputManager.TreeDoubleV["theta"]	= thetaNow;
	  OutputManager.TreeDoubleV["engyMC"]	= clusterClassMap[armNow][particleId]->EngyMC;
	  OutputManager.TreeDoubleV["thetaMC"]	= clusterClassMap[armNow][particleId]->ThetaMC;
	  OutputManager.TreeDoubleV["phi"]	= phiNow;
	  OutputManager.TreeDoubleV["rzStart"]	= rzstartNow;
	  OutputManager.TreeDoubleV["Xglob"]    = xglob;
	  OutputManager.TreeDoubleV["Yglob"]    = yglob;
	  OutputManager.TreeDoubleV["Zglob"]    = zglob;
	  //--
	  OutputManager.FillRootTree("LumiRecoParticleTree");
	}
      }

      /* fill in a flag entry (all values = -1) in the tree in the case that no arm has any clusters
      if(clusterInFlag == 0) {

          OutputManager.TreeIntV["nEvt"]	= NumEvt;
	  OutputManager.TreeIntV["outFlag"]	=  1;
	  OutputManager.TreeDoubleV["engy"]	= -1;
	  OutputManager.TreeDoubleV["theta"]	= -1;
	  OutputManager.TreeDoubleV["phi"]	= -1;
	  OutputManager.TreeDoubleV["rzStart"]	= -1;
	  OutputManager.TreeDoubleV["Xglob"]    = -1;
	  OutputManager.TreeDoubleV["Yglob"]    = -1;
	  OutputManager.TreeDoubleV["Zglob"]    = -1;

	OutputManager.TreeIntV["sign"]	= 1;
	//OutputManager.TreeMap["bhabhaSelectionTree"] -> Fill();
	OutputManager.FillRootTree(treeName);

	OutputManager.TreeIntV["sign"]	= -1;
	//	OutputManager.TreeMap["bhabhaSelectionTree"] -> Fill();
  	OutputManager.FillRootTree(treeName);
    }
      */
      // generated MC particles ( within fiducial volume of LumiCal

    int mcparticleInFlag = 0;
try
  {
   
    const double beta = tan( GlobalMethods.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle]/2.);
    /*
    const double gamma = 1./sqrt( 1. - beta*beta );
    const double betagamma = beta*gamma;
    */
    //(BP)  as in Mokka PrimaryGeneratorAction
    const double betagamma = beta;
    const double gamma = sqrt( 1. + beta*beta );

    const double thmin = GlobalMethods.GlobalParamD[GlobalMethodsClass::ThetaMin];
    const double thmax = GlobalMethods.GlobalParamD[GlobalMethodsClass::ThetaMax];
    const double LcalZstart = GlobalMethods.GlobalParamD[GlobalMethodsClass::ZStart]; 

    LCCollection * particles = evt->getCollection( "MCParticle" );

    if ( particles ){
      numMCParticles = particles->getNumberOfElements();
      for ( int jparticle=0; jparticle<numMCParticles; jparticle++ ){
	MCParticle * particle = dynamic_cast<MCParticle*>( particles->getElementAt(jparticle) );
	// only primary particles wanted
       	if( !particle->isCreatedInSimulation() ){
	  int pdg = particle->getPDG();
	  // skip neutrinos
	  if( abs(pdg) == 12 || abs(pdg) == 14 || abs(pdg) == 16 ) continue; 
	  double* p = (double*)particle->getMomentum();
	  double engy = particle->getEnergy();
       	  double pxlocal = -betagamma*engy + gamma*p[0];
	  double pTheta = atan( sqrt( sqr(pxlocal) + sqr(p[1]) )/fabs(p[2]));
	  // check if within Lcal acceptance
          if( fabs(pTheta) > thmin  && fabs(pTheta) < thmax ){
	    double *pos = (double*)particle->getVertex();
	    double *endPoint = (double*)particle->getEndpoint();
	    // did it reach at least LCal face ?
	    if( fabs( endPoint[2] ) >= LcalZstart  ){
	      // energy above min 
	      if( engy >= GlobalMethods.GlobalParamD[GlobalMethodsClass::MinClusterEngyGeV] ){
		mcparticleInFlag++;
		double phi = atan2( p[1], pxlocal);
		phi = ( phi >=0 )? phi : phi+2.*M_PI;
		int sign = int(p[2]/fabs( p[2]));
		OutputManager.TreeIntV["nEvt"]	= NumEvt;
		OutputManager.TreeIntV["sign"]	= sign;
		OutputManager.TreeIntV["pdg"]	= pdg;
		OutputManager.TreeDoubleV["engy"]	= engy;
		OutputManager.TreeDoubleV["theta"]	= pTheta;
		OutputManager.TreeDoubleV["phi"]	= phi; 
		// position at the face of LCal
		OutputManager.TreeDoubleV["vtxX"]	= pos[0]+LcalZstart*tan(p[0]/fabs(p[2]));
		OutputManager.TreeDoubleV["vtxY"]	= pos[1]+LcalZstart*tan(p[1]/fabs(p[2]));
		OutputManager.TreeDoubleV["vtxZ"]	= double(sign)*LcalZstart;
		OutputManager.TreeDoubleV["endX"]    = endPoint[0];
		OutputManager.TreeDoubleV["endY"]    = endPoint[1];
		OutputManager.TreeDoubleV["endZ"]    = endPoint[2];
		// fill tree
		OutputManager.FillRootTree("LumiMCParticleTree" );
	      } // energy 
	    } // if endPoint
	  } // if pTheta
       	} // if !isCreatedinSimulation 
      }// for jparticle
    }//if particles
  }// try
 catch ( DataNotAvailableException &e ){
   streamlog_out(WARNING)<< "MCParticle data not available for event "<< NumEvt << std::endl;
 }
   hisName = "NumClustersIn_numMCParticlesIn";
   OutputManager.HisMap2D[hisName] -> Fill ( mcparticleInFlag, clusterInFlag );

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
	 cleanUp
	 -------------------------------------------------------------------------- */
      for(int armNow = -1; armNow < 2; armNow += 2) {
	clusterClassMapIterator = clusterClassMap[armNow].begin();
	numMCParticles          = clusterClassMap[armNow].size();
	for (int MCParticleNow = 0; MCParticleNow < numMCParticles; MCParticleNow++, clusterClassMapIterator++) {
	  particleId = (int)(*clusterClassMapIterator).first;

	  delete	clusterClassMap[armNow][particleId];
	}
      }


      /* --------------------------------------------------------------------------
	 write to the root tree
	 -------------------------------------------------------------------------- */
      OutputManager.WriteToRootTree("" , NumEvt);

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
						std::map < int , std::map < int , ClusterClass * > > & clusterClassMap,
						EVENT::LCEvent * evt ) {


    std::map < int , ClusterClass * > :: iterator		clusterClassMapIterator;

    std::map < int , std::vector<int> >  ::const_iterator	clusterIdToCellIdIterator;
    std::vector < myMCP > mcParticlesVec;

    LCCollection * particles = evt->getCollection( "MCParticle" );

    if ( particles ){
      int numMCParticles = particles->getNumberOfElements();
      for ( int jparticle=0; jparticle<numMCParticles; jparticle++ ){
	MCParticle * particle = dynamic_cast<MCParticle*>( particles->getElementAt(jparticle) );
	// only primaries
	if( particle->isCreatedInSimulation() ) continue;
	double* mom = (double *)particle->getMomentum();
	double theta = atan( sqrt( sqr(mom[0]) + sqr(mom[1]))/fabs( mom[2]) ); 
	myMCP p = { particle, theta};
	mcParticlesVec.push_back( p );

      }
    }

    for(int armNow = -1; armNow < 2; armNow += 2) {
      clusterIdToCellIdIterator = clusterIdToCellId.at(armNow).begin();
      int numClusters               = clusterIdToCellId.at(armNow).size();
      for(int superClusterNow = 0; superClusterNow < numClusters; superClusterNow++, clusterIdToCellIdIterator++) {
	int clusterId = (int)(*clusterIdToCellIdIterator).first;

	clusterClassMap[armNow][clusterId] = new ClusterClass(clusterId );
#if _CLUSTER_RESET_STATS_DEBUG == 1
      streamlog_out(DEBUG) << "Initial Reset Stats for LumiCal arm = " << armNow << "  ...... " << std::endl;
#endif
	clusterClassMap[armNow][clusterId] -> SetStatsMC();
	clusterClassMap[armNow][clusterId] -> SignMC = armNow;

	int numElementsInCluster = clusterIdToCellId.at(armNow).at(clusterId).size();
	for(int cellNow = 0; cellNow < numElementsInCluster; cellNow++) {

	  int cellId  = clusterIdToCellId.at(armNow).at(clusterId).at(cellNow);
	  double engyHit = cellIdToCellEngy.at(armNow).at(clusterId).at(cellNow);

	  clusterClassMap[armNow][clusterId] -> FillHit(cellId , engyHit);
	}
      }
    }



#if _CLUSTER_RESET_STATS_DEBUG == 1
    streamlog_out(DEBUG) << "Transfering information into ClusterClass objects ......  "
			 << std::endl;
#endif

    /* --------------------------------------------------------------------------
       calculate the energy and position of each cluster
       -------------------------------------------------------------------------- */
    for(int armNow = -1; armNow < 2; armNow += 2) {

      clusterClassMapIterator = clusterClassMap[armNow].begin();
      int numClusters             = clusterClassMap[armNow].size();
      for (int MCParticleNow = 0; MCParticleNow < numClusters; MCParticleNow++, clusterClassMapIterator++) {
	int clusterId = (int)(*clusterClassMapIterator).first;
	int	resetStatsFlag = clusterClassMap[armNow][clusterId] -> ResetStats();
	_sortAgainstThisTheta = clusterClassMap[armNow][clusterId] -> Theta;
	sort( mcParticlesVec.begin(), mcParticlesVec.end(), ThetaCompAsc );
	clusterClassMap[armNow][clusterId]->SetStatsMC( mcParticlesVec[0].mcp );

  if(resetStatsFlag == 1) {
    continue; 
  } 

#if _CLUSTER_RESET_STATS_DEBUG == 1
	if(clusterClassMap[armNow][clusterId] -> SignMC != armNow) continue;

	streamlog_out( MESSAGE4 ) << "\tParticle Out ("
		  << clusterClassMap[armNow][clusterId] -> OutsideReason << "):   " << clusterId << std::endl
		  << "\t\t pdg, parentId , NumMCDaughters = "
		  << "\t" << clusterClassMap[armNow][clusterId] -> Pdg
		  << "\t" << clusterClassMap[armNow][clusterId] -> ParentId
		  << "\t" << clusterClassMap[armNow][clusterId] -> NumMCDaughters << std::endl
		  << "\t\t vertex -> endPoint X, Y, Z = "
		  << "\t" << clusterClassMap[armNow][clusterId] -> VtxX
		  << "\t" << clusterClassMap[armNow][clusterId] -> VtxY
		  << "\t" << clusterClassMap[armNow][clusterId] -> VtxZ << "   ->   "
		  << "\t" << clusterClassMap[armNow][clusterId] -> EndPointX
		  << "\t" << clusterClassMap[armNow][clusterId] -> EndPointY
		  << "\t" << clusterClassMap[armNow][clusterId] -> EndPointZ << std::endl
		  << "\t\t engy, theta, phi (mc):\t "
		  << clusterClassMap[armNow][clusterId] -> EngyMC  << "\t"
		  << clusterClassMap[armNow][clusterId] -> ThetaMC << "\t"
		  << clusterClassMap[armNow][clusterId] -> PhiMC
		  << std::endl << std::endl;
#endif
      }
    }

    return;

  }
