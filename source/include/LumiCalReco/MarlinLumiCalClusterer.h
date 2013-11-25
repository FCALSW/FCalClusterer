// check if this file has already been called:
#ifndef MarlinLumiCalClusterer_
#define MarlinLumiCalClusterer_

class LumiCalClustererClass;

// Marlin classes:
#include "marlin/Processor.h"
#include <marlin/Global.h>
#include "lcio.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimCalorimeterHit.h>
#include "IMPL/CalorimeterHitImpl.h"
#include "IMPL/SimCalorimeterHitImpl.h"
#include "IMPL/ClusterImpl.h"

// gear
#include <gear/GEAR.h>
#include <gear/CalorimeterParameters.h>
#include <gear/LayerLayout.h>


// cout color definitios
#define coutDefault		"\033[0m"
#define coutRed			"\033[31m"
#define coutGreen		"\033[32m"
#define coutBlue		"\033[34m"
#define coutPurple		"\033[35m"
#define coutUnderLine		"\033[4;30m"
#define coutWhiteOnBlack	"\33[40;37;1m"


// my standard includes:
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <assert.h>

// ROOT libs
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TColor.h>

// namespaces:
using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace TMath ;




namespace MarlinLumiCalClusterer {


#define _VERY_SMALL_NUMBER 1e-9
#define _CLUSTER_MERGE_DEBUG 1
#define _GLOBAL_COUNTERS_UPDATE_DEBUG 1
#define _BHABHA_SELECTION_CUTS_DEBUG 1
#define _CLUSTER_RESET_STATS_DEBUG 1

#include <GlobalMethodsClass.h>
#include <ClusterClass.h>
#include <SortingClass.h>
#include <OutputManagerClass.h>

  //#include <LumiCalClusterer/LumiCalClusterer.cc>



class MarlinLumiCalClusterer : public Processor {

public:

	// returns the processor instance
	virtual Processor * newProcessor() {return new MarlinLumiCalClusterer ;}
	
	// Constructor
	MarlinLumiCalClusterer() ;
	
	// initialization routine - Called at the begining of the job.
	virtual void init() ;
	
	// pre-run actions - Called for every run
	virtual void processRunHeader( LCRunHeader * run  ) ;
	
	// main actions in each event -Called for every event - the working horse. 
	virtual void processEvent( LCEvent * evt ) ;
	
	// final action after last event analysis is over.
	virtual void end() ;

//protected:

	// Processor Parameters
	string	LumiName ;
	int	SkipNEvents, MaxRecordNumber;

	// global counters
	int	NumRun, NumEvt, EvtNumber;
	string	OutDirName;

	int	NumEventsTree;

	OutputManagerClass	* OutputManager;
	GlobalMethodsClass	* GlobalMethods;
	LumiCalClustererClass	* LumiCalClusterer;

	void	TryMarlinLumiCalClusterer(LCEvent * evt);

	void	CreateClusters( map < int , map < int , vector<int> > > 	clusterIdToCellId,
				map < int , map < int , vector<double> > > 	cellIdToCellEngy,
				map < int , map < int , ClusterClass * > >	* clusterClassMapP );


};




}  // namespace MarlinLumiCalClusterer

#endif
