
#include "MarlinLumiCalClusterer.h"

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

/* -----------------------------------------------------------------------------------------*/
// (BP) defs for sorting

double _sortAgainstThisTheta;
double _sortAgainstThisX;
double _sortAgainstThisY;

struct myMCP{
  MCParticle * mcp;
  double theta;
  double x;
  double y;
};

bool ThetaCompAsc( myMCP a, myMCP b ){
  return fabs( _sortAgainstThisTheta - a.theta) < fabs( _sortAgainstThisTheta -b.theta) ;
}

bool PositionCompAsc( myMCP a, myMCP b ){
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

      std::string hisName, optionName;
      int numMCParticles, particleId;
      double	engyNow, thetaNow, phiNow, rzstartNow;

      std::map < int , std::map < int , ClusterClass * > >	clusterClassMap;
      std::map < int , ClusterClass * > :: iterator	clusterClassMapIterator;

      int	numClusters, clusterId;
      double	weightNow;

      csbx[-1] = cos( M_PI - GlobalMethods.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle]/2.);
      csbx[ 1] = cos( GlobalMethods.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle]/2.);
      snbx[-1] = sin( M_PI - GlobalMethods.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle]/2.);
      snbx[ 1] = sin( GlobalMethods.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle]/2.);
      /* --------------------------------------------------------------------------
	 create clusters using: LumiCalClustererClass
	 -------------------------------------------------------------------------- */

      if ( !LumiCalClusterer.processEvent(evt) ) return;


       
// instantiate a clusterClass object for each mcParticle which was created
// infront of LumiCal and was destroyed after lumical.
// (BP) this can optional
     CreateClusters( LumiCalClusterer._superClusterIdToCellId,
		      LumiCalClusterer._superClusterIdToCellEngy,
		      //		      LumiCalClusterer._superClusterIdClusterInfo,
		      clusterClassMap,
		      evt );


      /* --------------------------------------------------------------------------
	 flag the highest-energy cluster in each arm
	 -------------------------------------------------------------------------- */
      
      LCCollectionVec* LCalClusterCol = new LCCollectionVec(LCIO::CLUSTER);
      LCCollectionVec* LCalRPCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
      std::map < int , std::vector<int> > ::const_iterator  clusterIdToCellIdIterator;
      for(int armNow = -1; armNow < 2; armNow += 2) {
	clusterIdToCellIdIterator = LumiCalClusterer._superClusterIdToCellId[armNow].begin();
	numClusters               = LumiCalClusterer._superClusterIdToCellId[armNow].size();

	//	numClusters             = clusterClassMap[armNow].size();
	//	clusterClassMapIterator = clusterClassMap[armNow].begin();
	int highestEnergyClusterId = -1;
	double theHighestEnergy = 0.;
	//	for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterClassMapIterator++) {
	for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterIdToCellIdIterator++) {
	  clusterId = (int)(*clusterIdToCellIdIterator).first;

	  //	  ClusterClass const* thisCluster = clusterClassMap[armNow][clusterId];
	  LCCluster const& thisClusterInfo = LumiCalClusterer._superClusterIdClusterInfo[armNow][clusterId];

	  ClusterImpl* cluster = new ClusterImpl;
	  const double thetaCluster = thisClusterInfo.getTheta();
	  const double energyCluster = GlobalMethods.SignalGevConversion(GlobalMethodsClass::Signal_to_GeV , thisClusterInfo.getE());
	  const double phiCluster = thisClusterInfo.getPhi();
	  cluster->setEnergy( energyCluster );
          const float xloc =  float(thisClusterInfo.getX());
          const float yloc =  float(thisClusterInfo.getY());
          const float zloc =  fabs(float(thisClusterInfo.getZ()));
	  //(BP) local -> global rot.
	  float xglob = csbx[ armNow ]*xloc + snbx[ armNow ]*zloc;
          float zglob = snbx[ armNow ]*xloc + csbx[ armNow ]*zloc;
	  const float clusterPosition[3] = { xglob, yloc, zglob };

#pragma message ("FIXME: Link Calohits to the Cluster")
	  //Get the cellID from the container Hit of the cluster and find the
	  //matching LumiCal SimCalorimeterHit in the event collection, probably
	  //have to loop over all of them until it is found

#if _BP_DEBUG
	  //(BP)  cluster position from ClusterInfo are inconsistent with Thet, Phi values !
	  // i.e. atan( sqrt( xloc**2 + yloc**2 )/zloc != Theta
	  // and are dfferent then those calculated as energy log weighted average of hits positions belonging to this cluster
	  // What is this cluster position ?!
	  streamlog_out(MESSAGE)<< "\n ----------------------------------------------------------------------------------\n";
	  streamlog_out(MESSAGE)<< " Cluster "<<clusterId<<"\t Weight'n method: "<< thisClusterInfo.getMethod() <<"\n";
	  streamlog_out(MESSAGE)<< " Position from clusterInfo : X = "<< xloc <<"\t"<<thisCluster->clusterPosition[0]<<"\n"
				<< "                             Y = "<< yloc <<"\t"<<thisCluster->clusterPosition[1]<<"\n"
				<< "                             Z = "<< zloc <<"\t"<<thisCluster->clusterPosition[2]<<"\n";

	  double eneInfo = GlobalMethods.SignalGevConversion(GlobalMethodsClass::Signal_to_GeV , thisClusterInfo.getE());
	  double thetaInfo = atan( sqrt( sqr(xloc) + sqr(yloc) )/fabs( zloc ));
	  double phiInfo = atan2 ( yloc, xloc ); 
	  streamlog_out(MESSAGE) << " energy = "<< eneInfo                    <<"\t"<< energyCluster <<"\n"
				 << " theta  = "<< thisClusterInfo.getTheta() <<" ("<< thetaInfo <<" )" <<"\t"<< thetaCluster  <<"\n"
				 << " phi    = "<< thisClusterInfo.getPhi()   <<" ("<< phiInfo   <<" )" <<"\t"<< phiCluster    <<"\n";
#endif	  

	  cluster->setPosition( clusterPosition );
	  // (BP) take care about components sign
          double px = energyCluster * sin ( thetaCluster ) * cos ( phiCluster ) * float(armNow);
          double py = energyCluster * sin ( thetaCluster ) * sin ( phiCluster );
	  double pz = energyCluster * cos ( thetaCluster ) * float(armNow);
	  // (BP) do boost 
	  double pxlab = _betagamma * energyCluster + _gamma*px; 
	  float momentumCluster[3] = { float( pxlab ), float( py ), float( pz )  };

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
	  //	  thetaNow = atan( sqrt( sqr(xglob) + sqr(yglob) )/fabs(zglob) );
	  int mFlag = clusterClassMap[armNow][particleId]->Pdg;
	  int nHits = clusterClassMap[armNow][particleId]->NumHits;
	  int highE = clusterClassMap[armNow][particleId]->HighestEnergyFlag;
	  //-------------
	  OutputManager.TreeIntV["nEvt"]	= NumEvt;
	  OutputManager.TreeIntV["outFlag"]	= clusterClassMap[armNow][particleId]->OutsideFlag;
	  OutputManager.TreeIntV["mFlag"]	= mFlag;
	  OutputManager.TreeIntV["highE"]	= highE;
	  OutputManager.TreeIntV["sign"]	= armNow;
	  OutputManager.TreeIntV["nHits"]	= nHits;
	  OutputManager.TreeDoubleV["engy"]	= engyNow;
	  OutputManager.TreeDoubleV["theta"]	= thetaNow;
	  OutputManager.TreeDoubleV["engyMC"]	= clusterClassMap[armNow][particleId]->EngyMC;
	  OutputManager.TreeDoubleV["thetaMC"]	= clusterClassMap[armNow][particleId]->ThetaMC;
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
   //
   // generated MC particles ( within fiducial volume of LumiCal

    int mcparticleInFlag = 0;

try
  {
   

    const double thmin = GlobalMethods.GlobalParamD[GlobalMethodsClass::ThetaMin];
    const double thmax = GlobalMethods.GlobalParamD[GlobalMethodsClass::ThetaMax];
    const double LcalZstart = GlobalMethods.GlobalParamD[GlobalMethodsClass::ZStart]; 

    LCCollection * particles = evt->getCollection( "MCParticle" );

    if ( particles ){
      numMCParticles = particles->getNumberOfElements();
      for ( int jparticle=0; jparticle<numMCParticles; jparticle++ ){
	MCParticle * particle = dynamic_cast<MCParticle*>( particles->getElementAt(jparticle) );
	// only primary particles wanted
       	if( particle->isCreatedInSimulation() ) continue;
	  int pdg = particle->getPDG();
	  // skip neutrinos
	  if( abs(pdg) == 12 || abs(pdg) == 14 || abs(pdg) == 16 ) continue; 
	  double* p = (double*)particle->getMomentum();
	  double engy = particle->getEnergy();
       	  double pxlocal = -_betagamma*engy + _gamma*p[0];
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
		// position at the face of LCal (global coordinates)
		OutputManager.TreeDoubleV["begX"]	= pos[0]+double(sign)*LcalZstart*tan(p[0]/(p[2]));
		OutputManager.TreeDoubleV["begY"]	= pos[1]+double(sign)*LcalZstart*tan(p[1]/(p[2]));
		OutputManager.TreeDoubleV["begZ"]	= double(sign)*LcalZstart;
		OutputManager.TreeDoubleV["vtxX"]	= pos[0];
		OutputManager.TreeDoubleV["vtxY"]	= pos[1];
		OutputManager.TreeDoubleV["vtxZ"]	= pos[2];
		OutputManager.TreeDoubleV["endX"]    = endPoint[0];
		OutputManager.TreeDoubleV["endY"]    = endPoint[1];
		OutputManager.TreeDoubleV["endZ"]    = endPoint[2];
		// fill tree
		OutputManager.FillRootTree("LumiMCParticleTree" );
	      } // energy 
	    } // if endPoint
	  } // if pTheta
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

    const double LcalZstart = GlobalMethods.GlobalParamD[GlobalMethodsClass::ZStart]; 

#if _CLUSTER_RESET_STATS_DEBUG == 1
	streamlog_out(MESSAGE4) << "------------------------------------------------------------------------------" << std::endl;
	streamlog_out(MESSAGE4) << "CreateClusters: event = " <<NumEvt  << std::endl;
#endif
    if ( particles ){
      int numMCParticles = particles->getNumberOfElements();
      for ( int jparticle=0; jparticle<numMCParticles; jparticle++ ){
	MCParticle * particle = dynamic_cast<MCParticle*>( particles->getElementAt(jparticle) );
	// only primaries
	if( particle->isCreatedInSimulation() ) continue;
	int pdg = particle->getPDG();
	// skip neutrinos
	if( abs(pdg) == 12 || abs(pdg) == 14 || abs(pdg) == 16 ) continue;
	double *endPoint = (double*)particle->getEndpoint();
	// did it reach at least LCal face ?
	if( fabs( endPoint[2] ) >= LcalZstart  ){
	double engy = particle->getEnergy();
	if( engy < GlobalMethods.GlobalParamD[GlobalMethodsClass::MinClusterEngyGeV] ) continue;
	double* mom = (double *)particle->getMomentum();
	double sign = double( int( mom[2]/fabs(mom[2]))); 
	double* pos = (double *)particle->getVertex();
	// particle position at LCal face ( global )
	double  x0  = pos[0] + sign*LcalZstart*tan(mom[0]/mom[2]);
	double  y0  = pos[1] + sign*LcalZstart*tan(mom[1]/mom[2]);
	// unbooost momentum as clusters theta are
       	double pxlocal = _betagamma*engy -_gamma*mom[0];
	double theta = atan( sqrt( sqr(pxlocal) + sqr(mom[1]))/fabs( mom[2]) ); 
	myMCP p = { particle, theta, x0, y0};
	mcParticlesVec.push_back( p );
	}
      }
    }

    for(int armNow = -1; armNow < 2; armNow += 2) {
      clusterIdToCellIdIterator = clusterIdToCellId.at(armNow).begin();
      int numClusters               = clusterIdToCellId.at(armNow).size();
      for(int superClusterNow = 0; superClusterNow < numClusters; superClusterNow++, clusterIdToCellIdIterator++) {
	int clusterId = (int)(*clusterIdToCellIdIterator).first;
	//	LCCluster const& thisClusterInfo = LumiCalClusterer._superClusterIdClusterInfo[armNow][clusterId];

	// create a new cluster and put it on map
	clusterClassMap[armNow][clusterId] = new ClusterClass(clusterId );

      //	clusterClassMap[armNow][clusterId] -> SetStatsMC(); // done in creator
	clusterClassMap[armNow][clusterId] -> SignMC = armNow;

	int numElementsInCluster = clusterIdToCellId.at(armNow).at(clusterId).size();
	for(int cellNow = 0; cellNow < numElementsInCluster; cellNow++) {

	  int cellId  = clusterIdToCellId.at(armNow).at(clusterId).at(cellNow);
	  double engyHit = cellIdToCellEngy.at(armNow).at(clusterId).at(cellNow);

	  clusterClassMap[armNow][clusterId] -> FillHit(cellId , engyHit);
	  // calculate energy, position for the cluster
	}
#if _CLUSTER_RESET_STATS_DEBUG == 1
	streamlog_out(MESSAGE4) << "Initial Reset Stats for LumiCal arm = " << armNow <<"\t cluster "<< clusterId<< "  ...... " << std::endl;
#endif
	  clusterClassMap[armNow][clusterId] -> ResetStats();
	  /*
	  clusterClassMap[armNow][clusterId] -> clusterPosition[0] = thisClusterInfo.getX();
	  clusterClassMap[armNow][clusterId] -> clusterPosition[1] = thisClusterInfo.getY();
	  clusterClassMap[armNow][clusterId] -> clusterPosition[2] = thisClusterInfo.getZ();
	  clusterClassMap[armNow][clusterId] -> Theta     = thisClusterInfo.getTheta();
	  clusterClassMap[armNow][clusterId] -> Phi       = thisClusterInfo.getPhi();
	  clusterClassMap[armNow][clusterId] -> RZStart   = LcalZstart * tan( clusterClassMap[armNow][clusterId] -> Theta );
	  */
      }
    }



#if _CLUSTER_RESET_STATS_DEBUG == 1
    streamlog_out(MESSAGE4) << "Transfering MC information into ClusterClass objects ......  "
			 << std::endl;
#endif

    for(int armNow = -1; armNow < 2; armNow += 2) {

      clusterClassMapIterator = clusterClassMap[armNow].begin();
      int numClusters             = clusterClassMap[armNow].size();
      for (int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterClassMapIterator++) {
	int clusterId = (int)(*clusterClassMapIterator).first;
	//	int	resetStatsFlag = clusterClassMap[armNow][clusterId] -> ResetStats();
	// try to match MC true particle with cluster
	_sortAgainstThisTheta = clusterClassMap[armNow][clusterId] -> Theta;

	double phi = clusterClassMap[armNow][clusterId] -> Phi;
	double RZStart = clusterClassMap[armNow][clusterId] -> RZStart;
	double xglob = RZStart*cos(phi)*csbx[armNow] + LcalZstart*snbx[armNow];
	double yglob = RZStart*sin(phi);
	_sortAgainstThisX = xglob;
	_sortAgainstThisY = yglob;

	if( mcParticlesVec.size() ){
	  //	  sort( mcParticlesVec.begin(), mcParticlesVec.end(), ThetaCompAsc );
	  sort( mcParticlesVec.begin(), mcParticlesVec.end(), PositionCompAsc );
	  double det = fabs( _sortAgainstThisTheta - mcParticlesVec[0].theta );
	  // 10 mrad tolerance
	  if(  det < 0.005 ) clusterClassMap[armNow][clusterId]->SetStatsMC( mcParticlesVec[0].mcp );
	}

	if( clusterClassMap[armNow][clusterId] -> Pdg == 0)    continue;  

#if _CLUSTER_RESET_STATS_DEBUG == 1
	if(clusterClassMap[armNow][clusterId] -> SignMC != armNow) continue;

	streamlog_out( MESSAGE4 ) << "\tParticle Out ("
		  << clusterClassMap[armNow][clusterId] -> OutsideReason << "):   " << clusterId << std::endl
		  << "\t\t side, pdg, parentId , NumMCDaughters = "
		  << "\t" << clusterClassMap[armNow][clusterId] -> SignMC
		  << "\t" << clusterClassMap[armNow][clusterId] -> Pdg
		  << "\t" << clusterClassMap[armNow][clusterId] -> ParentId
		  << "\t" << clusterClassMap[armNow][clusterId] -> NumMCDaughters << std::endl
		  << "\t\t vertex   X, Y, Z = "
		  << "\t" << clusterClassMap[armNow][clusterId] -> VtxX
		  << "\t" << clusterClassMap[armNow][clusterId] -> VtxY
		  << "\t" << clusterClassMap[armNow][clusterId] -> VtxZ << std::endl
		  << "\t\t endPoint X, Y, Z = "
		  << "\t" << clusterClassMap[armNow][clusterId] -> EndPointX
		  << "\t" << clusterClassMap[armNow][clusterId] -> EndPointY
		  << "\t" << clusterClassMap[armNow][clusterId] -> EndPointZ << std::endl
		  << "\t\t engy, theta, phi (mc): \t "
		  << clusterClassMap[armNow][clusterId] -> EngyMC  << "\t"
		  << clusterClassMap[armNow][clusterId] -> ThetaMC << "\t"
		  << clusterClassMap[armNow][clusterId] -> PhiMC << std::endl
		  << "\t\t engy, theta, phi (reco):\t "
		  << GlobalMethods.SignalGevConversion(GlobalMethodsClass::Signal_to_GeV ,clusterClassMap[armNow][clusterId] -> Engy) << "\t"
		  << clusterClassMap[armNow][clusterId] -> Theta << "\t"
		  << clusterClassMap[armNow][clusterId] -> Phi
		  << std::endl << std::endl;
#endif
      }
    }

    return;

  }

