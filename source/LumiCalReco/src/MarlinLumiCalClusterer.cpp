
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
#include <iomanip>


/* >> */ 

/* -----------------------------------------------------------------------------------------*/
// (BP) defs for sorting

double _sortAgainstThisTheta;
double _sortAgainstThisX;
double _sortAgainstThisY;

struct myMCP{
  MCParticle * mcp;
  double engy;
  double theta;
  double x;
  double y;
};
// variable _sortAgainstThisTheta must be set before sorting 
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
      int         numMCParticles, particleId;
      double	  engyNow, thetaNow, phiNow, rzstartNow;

      std::map < int , std::map < int , ClusterClass * > >	clusterClassMap;
      std::map < int , ClusterClass * > :: iterator	clusterClassMapIterator;

      int	numClusters, clusterId;

      double ThetaMid = (GlobalMethods.GlobalParamD[GlobalMethodsClass::ThetaMin] + GlobalMethods.GlobalParamD[GlobalMethodsClass::ThetaMax])/2.;
      double ThetaTol = (GlobalMethods.GlobalParamD[GlobalMethodsClass::ThetaMax] - GlobalMethods.GlobalParamD[GlobalMethodsClass::ThetaMin])/2.;


      csbx[-1] = cos( M_PI - GlobalMethods.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle]/2.);
      csbx[ 1] = cos( GlobalMethods.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle]/2.);
      snbx[-1] = sin( M_PI - GlobalMethods.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle]/2.);
      snbx[ 1] = sin( GlobalMethods.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle]/2.);
      /* --------------------------------------------------------------------------
	 create clusters using: LumiCalClustererClass
	 -------------------------------------------------------------------------- */

      if ( !LumiCalClusterer.processEvent(evt) ) return;




      
      LCCollectionVec* LCalClusterCol = new LCCollectionVec(LCIO::CLUSTER);
      LCCollectionVec* LCalRPCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
      std::map < int , std::vector<int> > ::const_iterator  clusterIdToCellIdIterator;

#if _CREATE_CLUSTERS_DEBUG == 1
      std::cout << " Transfering reco results to LCalClusterCollection....."<<std::endl;
#endif 
  
      for(int armNow = -1; armNow < 2; armNow += 2) {
	clusterIdToCellIdIterator = LumiCalClusterer._superClusterIdToCellId[armNow].begin();
	numClusters               = LumiCalClusterer._superClusterIdToCellId[armNow].size();

	//	numClusters             = clusterClassMap[armNow].size();
	//	clusterClassMapIterator = clusterClassMap[armNow].begin();
	//	for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterClassMapIterator++) {
#if _CREATE_CLUSTERS_DEBUG == 1
       	std::cout<<" Arm  "<< std::setw(4)<< armNow << "\t Number of clusters: "<< numClusters <<std::endl;
#endif
	for(int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterIdToCellIdIterator++) {
	  clusterId = (int)(*clusterIdToCellIdIterator).first;

	  //	  ClusterClass const* thisCluster = clusterClassMap[armNow][clusterId];
	  LCCluster const& thisClusterInfo = LumiCalClusterer._superClusterIdClusterInfo[armNow][clusterId];

	  const double clusterEnergy = GlobalMethods.SignalGevConversion(GlobalMethodsClass::Signal_to_GeV , thisClusterInfo.getE());
	  if( clusterEnergy < _minClusterEngy ) continue;

	  const double clusterTheta = thisClusterInfo.getTheta();
	  if( fabs ( clusterTheta - ThetaMid ) >  ThetaTol ) continue;
 
	  const double clusterPhi = thisClusterInfo.getPhi();
	  const float  xloc =  float(thisClusterInfo.getX());
          const float  yloc =  float(thisClusterInfo.getY());
          const float  zloc =  fabs(float(thisClusterInfo.getZ())); // take abs , since we kept wrong local minus z in getHts
#if _CREATE_CLUSTERS_DEBUG == 1
	  std::cout<<std::setw(8) << std::fixed << std::setprecision(3)
		   << " Cluster "<< clusterId << ": Position X,Y,Z [mm] ( "<< xloc <<", "<< yloc <<", "<< zloc 
		   <<")\t E( "<< clusterEnergy <<" GeV)" 
		   <<std::endl;
#endif
	  //(BP) local -> global rot.
	  float xglob = csbx[ armNow ]*xloc + snbx[ armNow ]*zloc;
          float zglob = snbx[ armNow ]*xloc + csbx[ armNow ]*zloc;
	  const float clusterPosition[3] = { xglob, yloc, zglob };
	  
	  ClusterImpl* cluster = new ClusterImpl;
	  cluster->setEnergy( clusterEnergy );
 	  cluster->setPosition( clusterPosition );
#pragma message ("FIXME: Link Calohits to the Cluster")
	  //Get the cellID from the container Hit of the cluster and find the
	  //matching LumiCal SimCalorimeterHit in the event collection, probably
	  //have to loop over all of them until it is found

	  ReconstructedParticleImpl* particle = new ReconstructedParticleImpl;
	  const float mass = 0.0;
	  const float charge = 1e+19;
	  // (BP) take care about components sign
          double px = clusterEnergy * sin ( clusterTheta ) * cos ( clusterPhi ) * float(armNow);
          double py = clusterEnergy * sin ( clusterTheta ) * sin ( clusterPhi );
	  double pz = clusterEnergy * cos ( clusterTheta ) * float(armNow);
	  // (BP) do boost 
	  double pxlab = _betagamma * clusterEnergy + _gamma*px; 
	  float clusterMomentum[3] = { float( pxlab ), float( py ), float( pz )  };

	  particle->setMass( mass ) ;
	  particle->setCharge( charge ) ;
	  particle->setMomentum ( clusterMomentum) ;
	  particle->setEnergy ( clusterEnergy ) ;
	  particle->addCluster( cluster ) ;

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

// instantiate a clusterClass object for each mcParticle which was created
// infront of LumiCal and was destroyed after lumical.
     
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
	  OutputManager.TreeDoubleV["distTheta"]= clusterClassMap[armNow][particleId]->DiffTheta;
	  OutputManager.TreeDoubleV["distXY"]   = clusterClassMap[armNow][particleId]->DiffPosXY;
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
    const double rmin = LcalZstart*tan(thmin);
    const double rmax = LcalZstart*tan(thmax);
    const double r0   = LcalZstart*tan(_BeamCrossingAngle);
    LCCollection * particles = evt->getCollection( "MCParticle" );

    if ( particles ){
      numMCParticles = particles->getNumberOfElements();
      for ( int jparticle=0; jparticle<numMCParticles; jparticle++ ){
	MCParticle *particle = dynamic_cast<MCParticle*>( particles->getElementAt(jparticle) );
	// only primary particles wanted
       	if( particle->isCreatedInSimulation() ) continue;                           //<---- is primary ? 
	  int pdg = particle->getPDG();
	  // skip neutrinos
	  if( abs(pdg) == 12 || abs(pdg) == 14 || abs(pdg) == 16 ) continue;                //<---- detectable ?
	  // energy above min 
	  double engy = particle->getEnergy();
	  if( engy < GlobalMethods.GlobalParamD[GlobalMethodsClass::MinClusterEngyGeV] ) continue; 
	  double* p = (double*)particle->getMomentum();
	  int sign = int ( p[2]/fabs( p[2] ));
	  // check if within Lcal acceptance
	  double *endPoint = (double*)particle->getEndpoint();
	  // did it reach at least LCal face ?
	  if( fabs( endPoint[2] ) < LcalZstart && fabs( endPoint[2] ) > 0. ) continue;       //<---- endPoint set and particle did not make to Lcal 
	  double *vx = (double*)particle->getVertex();
	  double begX = vx[0]+double(sign)*LcalZstart*tan(p[0]/(p[2]));
	  double begY = vx[1]+double(sign)*LcalZstart*tan(p[1]/(p[2]));
	  double rt = sqrt( sqr(begX-r0) + sqr(begY) );
	  if( rt < rmin || rt > rmax ) continue;                                           //<---- within geo  acceptance ?
       	  double pxlocal = -_betagamma*engy + _gamma*p[0];
	  double pTheta = atan( sqrt( sqr(pxlocal) + sqr(p[1]) )/fabs(p[2]));
          if( fabs(pTheta) < thmin  || fabs(pTheta) > thmax ) continue;         //<--- within theta range ? 
		mcparticleInFlag++;
		double phi = atan2( p[1], pxlocal);

		OutputManager.TreeIntV["nEvt"]	= NumEvt;
		OutputManager.TreeIntV["sign"]	= sign;
		OutputManager.TreeIntV["pdg"]	= pdg;
		OutputManager.TreeDoubleV["engy"]	= engy;
		OutputManager.TreeDoubleV["theta"]	= pTheta;
		OutputManager.TreeDoubleV["phi"]	= phi; 
		// position at the face of LCal (global coordinates)
		OutputManager.TreeDoubleV["begX"]	= begX;
		OutputManager.TreeDoubleV["begY"]	= begY;
		OutputManager.TreeDoubleV["begZ"]	= double(sign)*LcalZstart;
		OutputManager.TreeDoubleV["vtxX"]	= vx[0];
		OutputManager.TreeDoubleV["vtxY"]	= vx[1];
		OutputManager.TreeDoubleV["vtxZ"]	= vx[2];
		OutputManager.TreeDoubleV["endX"]    = endPoint[0];
		OutputManager.TreeDoubleV["endY"]    = endPoint[1];
		OutputManager.TreeDoubleV["endZ"]    = endPoint[2];
		// fill tree
		OutputManager.FillRootTree("LumiMCParticleTree" );
	  
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
	 write to the root tree
	 -------------------------------------------------------------------------- */
      OutputManager.WriteToRootTree("" , NumEvt);

      /* --------------------------------------------------------------------------
	 clean ClusterClassMap
	 -------------------------------------------------------------------------- */
      for(int armNow = -1; armNow < 2; armNow += 2) {
	clusterClassMapIterator = clusterClassMap[armNow].begin();
	numClusters          = clusterClassMap[armNow].size();
	for (int numC = 0; numC < numClusters; numC++, clusterClassMapIterator++) {
	  particleId = (int)(*clusterClassMapIterator).first;

	  delete clusterClassMap[armNow][particleId];
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
						std::map < int , std::map < int , ClusterClass * > > & clusterClassMap,
						EVENT::LCEvent * evt ) {


    std::map < int , ClusterClass * > :: iterator		clusterClassMapIterator;
    std::map < int , std::vector<int> >  ::const_iterator	clusterIdToCellIdIterator;

    std::vector < myMCP > mcParticlesVecPos;
    std::vector < myMCP > mcParticlesVecNeg;


    const double LcalZstart = GlobalMethods.GlobalParamD[GlobalMethodsClass::ZStart];
    const double Rmin = GlobalMethods.GlobalParamD[GlobalMethodsClass::RMin]; 
    const double Rmax = GlobalMethods.GlobalParamD[GlobalMethodsClass::RMax]; 
    const double thmin = GlobalMethods.GlobalParamD[GlobalMethodsClass::ThetaMin];
    const double thmax = GlobalMethods.GlobalParamD[GlobalMethodsClass::ThetaMax];

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
	MCParticle * particle = dynamic_cast<MCParticle*>( particles->getElementAt(jparticle) );
	// only primary particles wanted
	if( particle->isCreatedInSimulation() ) continue;                           //<--- is primary ? 
	  int pdg = particle->getPDG();
	  // skip neutrinos
	  if( abs(pdg) == 12 || abs(pdg) == 14 || abs(pdg) == 16 ) continue;                        //<--- detectable ?
	  // energy above min 
	  double engy = particle->getEnergy();
	  if( engy < GlobalMethods.GlobalParamD[GlobalMethodsClass::MinClusterEngyGeV] ) continue;  //<--- Energy >  Emin ?
	  double* pp = (double*)particle->getMomentum();
	  int sign = int ( pp[2]/fabs( pp[2] ));
          // check if within Lcal acceptance
	  double *endPoint = (double*)particle->getEndpoint();
         	// did it reach at least LCal face
       	  if( fabs( endPoint[2] ) < LcalZstart && fabs( endPoint[2] )>0. ) continue;                 //<--- reached LCal ?
       	  double pxloc = -_betagamma*engy + _gamma*pp[0];
	  double *vx = (double*)particle->getVertex();
        	// particle position at LCal face ( local )
	  double begX = vx[0]+ LcalZstart*tan(pxloc/(pp[2]));
	  double begY = vx[1]+double(sign)*LcalZstart*tan(pp[1] /(pp[2]));
	  double rt = sqrt( sqr(begX) + sqr(begY) );
       	  if( rt < Rmin || rt > Rmax ) continue;                                                    //<--- within geo  acceptance ?
        	//  polar angle ( local )
	  double theta = atan( sqrt( sqr(pxloc) + sqr(pp[1]) )/fabs(pp[2]));
          if( fabs(theta) < thmin  || fabs(theta) > thmax ) continue;                  //<--- within theta range ? 

 
	  myMCP p = { particle, engy, theta, begX, begY};
	  //	  std::cout << " pdg, engy, theta "<< pdg << ", "<< engy <<", "<< theta << std::endl;
	if ( sign > 0 ) mcParticlesVecPos.push_back( p );
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
	clusterIdToCellIdIterator = clusterIdToCellId.at(armNow).begin();
	int numClusters               = clusterIdToCellId.at(armNow).size();
	double EngyMax =  0.;
	int  EngyMaxID = -1 ;
	for(int superClusterNow = 0; superClusterNow < numClusters; superClusterNow++, clusterIdToCellIdIterator++){
	  int clusterId = (int)(*clusterIdToCellIdIterator).first;
	// create a new cluster and put it on map
	  clusterClassMap[armNow][clusterId] = new ClusterClass(clusterId );
	  clusterClassMap[armNow][clusterId] -> SignMC = armNow;

	  int numElementsInCluster = clusterIdToCellId.at(armNow).at(clusterId).size();
	  double engySum = 0.;
	  for(int cellNow = 0; cellNow < numElementsInCluster; cellNow++) {
	    int cellId     = clusterIdToCellId.at(armNow).at(clusterId).at(cellNow);
	    double engyHit = cellIdToCellEngy.at(armNow).at(clusterId).at(cellNow);
	    clusterClassMap[armNow][clusterId] -> FillHit(cellId , engyHit);
	    engySum += engyHit;
	  }
	  if( engySum > EngyMax ) {
	    EngyMax   = engySum;
	    EngyMaxID = clusterId;
	  }
	
	  clusterClassMap[armNow][clusterId] -> ResetStats(); // calculate energy, position for the cluster
	
#if _CREATE_CLUSTERS_DEBUG == 1
       	LCCluster const& thisClusterInfo = LumiCalClusterer._superClusterIdClusterInfo[armNow][clusterId];
	streamlog_out(MESSAGE) << "arm =   " << armNow <<"\t cluster "<< clusterId<< "  ...... " << std::endl
			       << std::setw(20) << "X, Y, Z:            " << std::endl
			       << std::setw(20) << "ClusterClass        "
			       << std::fixed << std::setprecision(3)
			       << std::setw(13) << clusterClassMap[armNow][clusterId]-> clusterPosition[0]
			       << std::setw(13) << clusterClassMap[armNow][clusterId]-> clusterPosition[1]
			       << std::setw(13) << clusterClassMap[armNow][clusterId]-> clusterPosition[2] << std::endl
			       << std::setw(20) << "CalusterInfo        " 
			       << std::setw(13) << thisClusterInfo.getX()
			       << std::setw(13) << thisClusterInfo.getY()
			       << std::setw(13) << thisClusterInfo.getZ() << std::endl
			       << std::setw(20) << "Energy, Theta, Phi: " << std::endl
			       << std::setw(20) << "ClusterClass        "
			       << std::setw(13) << clusterClassMap[armNow][clusterId]-> Engy
			       << std::setw(13) << clusterClassMap[armNow][clusterId]-> Theta
			       << std::setw(13) << clusterClassMap[armNow][clusterId]-> Phi  << std::endl
			       << std::setw(20) << "ClusterInfo         " 
			       << std::setw(13) << thisClusterInfo.getEnergy()
			       << std::setw(13) << thisClusterInfo.getTheta()
			       << std::setw(13) << thisClusterInfo.getPhi()  << std::endl
			       << std::endl;
#endif
	}
	// set flag for highest energy cluster found in this event
	if ( EngyMax > 0.0 ) clusterClassMap[armNow][EngyMaxID]->HighestEnergyFlag = 1;
      }



#if _CREATE_CLUSTERS_DEBUG == 1
      if( mcParticlesVecNeg.size() || mcParticlesVecPos.size() ) {
	streamlog_out(MESSAGE4) << "Transfering MC information into ClusterClass objects ......  "  << std::endl;
	streamlog_out(MESSAGE4) << "MC particles arm z-: "<< mcParticlesVecNeg.size() << "\t arm z+: "<< mcParticlesVecPos.size() << std::endl;
      }else{
	streamlog_out(MESSAGE4) << "Transfering MC information into ClusterClass objects: No primary MCParticles enetring LumiCal found!"  << std::endl;
      }
#endif
    
      if( mcParticlesVecNeg.size() || mcParticlesVecPos.size() ) {
	for(int armNow = -1; armNow < 2; armNow += 2) {

	  std::vector < myMCP > mcParticlesVec;
	  if ( armNow == -1 ) mcParticlesVec = mcParticlesVecNeg ;
	  else                mcParticlesVec = mcParticlesVecPos ;

	  clusterClassMapIterator = clusterClassMap[armNow].begin();
	  int numClusters             = clusterClassMap[armNow].size();
	  for (int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterClassMapIterator++) {
	    int clusterId = (int)(*clusterClassMapIterator).first;


	    double eneCL = GlobalMethods.SignalGevConversion(GlobalMethodsClass::Signal_to_GeV ,clusterClassMap[armNow][clusterId] -> Engy);
	    double phiCL = clusterClassMap[armNow][clusterId] -> Phi;
	    double RZStart = clusterClassMap[armNow][clusterId] -> RZStart;
	    double xs = RZStart*cos(phiCL);
	    double ys = RZStart*sin(phiCL);
	    _sortAgainstThisX = xs;
	    _sortAgainstThisY = ys;
	    _sortAgainstThisTheta = clusterClassMap[armNow][clusterId] -> Theta;

	    if( mcParticlesVec.size() ){
	// try to match MC true particle with cluster by comparing positions at Lumical entry
	  //	  sort( mcParticlesVec.begin(), mcParticlesVec.end(), ThetaCompAsc );
	      sort( mcParticlesVec.begin(), mcParticlesVec.end(), PositionCompAsc );

	      double dTheta  = fabs( _sortAgainstThisTheta - mcParticlesVec[0].theta );
	      double dPos0    = sqrt( ( sqr(xs - mcParticlesVec[0].x) + sqr( ys - mcParticlesVec[0].y)));
	      double dEne0    = fabs( mcParticlesVec[0].engy - eneCL);
	      clusterClassMap[armNow][clusterId] -> DiffTheta = dTheta;
	      if( mcParticlesVec.size() > 1 ) {
		double dPos1    = sqrt( ( sqr(xs - mcParticlesVec[1].x) + sqr( ys - mcParticlesVec[1].y)));
                double dEne1    = fabs( mcParticlesVec[1].engy - eneCL);
		if( dPos1 <  GlobalMethods.GlobalParamD[GlobalMethodsClass::MinSeparationDist] && dEne1 < dEne0 ) {
		  clusterClassMap[armNow][clusterId]->SetStatsMC( mcParticlesVec[1].mcp );
		  clusterClassMap[armNow][clusterId] -> DiffPosXY = dPos1;
		  mcParticlesVec.erase( mcParticlesVec.begin()+1 ); 
		}else{
		  clusterClassMap[armNow][clusterId]->SetStatsMC( mcParticlesVec[0].mcp );
		  clusterClassMap[armNow][clusterId] -> DiffPosXY = dPos0;
		  mcParticlesVec.erase( mcParticlesVec.begin() ); 
		}
	      }else{
		clusterClassMap[armNow][clusterId] -> DiffPosXY = dPos0;
		if( dPos0 <  GlobalMethods.GlobalParamD[GlobalMethodsClass::MinSeparationDist] ){  
		  clusterClassMap[armNow][clusterId]->SetStatsMC( mcParticlesVec[0].mcp );
		  mcParticlesVec.erase( mcParticlesVec.begin() ); 
		}
	      }
	    }

#if _CREATE_CLUSTERS_DEBUG == 1
	    if( clusterClassMap[armNow][clusterId] -> Pdg != 0) {   
	      double Ereco = GlobalMethods.SignalGevConversion(GlobalMethodsClass::Signal_to_GeV ,clusterClassMap[armNow][clusterId] -> Engy); 
	streamlog_out( MESSAGE4 ) << "\tParticle Out ("
		  << clusterClassMap[armNow][clusterId] -> OutsideReason << "):   " << clusterId << std::endl
		  << "\t\t side(arm), pdg, parentId , NumMCDaughters = "
		  << "\t" << clusterClassMap[armNow][clusterId] -> SignMC<<"("<<armNow<<")"
		  << "\t" << clusterClassMap[armNow][clusterId] -> Pdg
		  << "\t" << clusterClassMap[armNow][clusterId] -> ParentId
		  << "\t" << clusterClassMap[armNow][clusterId] -> NumMCDaughters << std::endl
		  << "\t\t MCPos    X, Y, Z: "
		  << "\t" << std::setw(13) << clusterClassMap[armNow][clusterId] -> mcpPosition[0]
		  << "\t" << std::setw(13) << clusterClassMap[armNow][clusterId] -> mcpPosition[1]
		  << "\t" << std::setw(13) << clusterClassMap[armNow][clusterId] -> mcpPosition[2]<< std::endl
		  << "\t\t Cluster  X, Y, Z: "
		  << "\t" << std::setw(13) << clusterClassMap[armNow][clusterId] -> clusterPosition[0]
		  << "\t" << std::setw(13) << clusterClassMap[armNow][clusterId] -> clusterPosition[1]
		  << "\t" << std::setw(13) << clusterClassMap[armNow][clusterId] -> clusterPosition[2] << std::endl
		  << "\t\t engy, theta, phi (mc): "
		  <<  std::setw(13)<< clusterClassMap[armNow][clusterId] -> EngyMC  
		  << "\t"<< std::setw(13)<<clusterClassMap[armNow][clusterId] -> ThetaMC 
		  << "\t"<< std::setw(13)<<clusterClassMap[armNow][clusterId] -> PhiMC << std::endl
		  << "\t\t engy, theta, phi (rec):"
		  << std::setw(13)<<Ereco
		  << "\t"<< std::setw(13)<<clusterClassMap[armNow][clusterId] -> Theta 
		  << "\t"<< std::setw(13)<<clusterClassMap[armNow][clusterId] -> Phi
		  << std::endl << std::endl;
	    }else{
	      streamlog_out( MESSAGE4 ) << "Arm: "<< armNow << "\t Cluster: " << clusterId << " does not match MC particle !\n";
	      clusterClassMap[armNow][clusterId]->PrintInfo();
	    }
#endif
	  }
	}
      }
    }
    return;

  }

