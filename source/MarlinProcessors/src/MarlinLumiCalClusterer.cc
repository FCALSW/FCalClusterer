#include "Global.hh"
#include "MarlinLumiCalClusterer.h"

#include <EVENT/LCEvent.h>

#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/CalorimeterParameters.h>
#include <gear/LayerLayout.h>

#include <string>
#include <map>

///////////////////////////////////////////////////////////////////////////////
// GENERAL NOTES:
// ============================================================================
///////////////////////////////////////////////////////////////////////////////



/* ============================================================================
   Register all Parameters that are activated in the steering file:
   ========================================================================= */
// create an instance of the MarlinLumiCalClusterer class
MarlinLumiCalClusterer aMarlinLumiCalClusterer ;
// constructor for MarlinLumiCalClusterer - MarlinLumiCalClusterer is an instance of Processor
MarlinLumiCalClusterer::MarlinLumiCalClusterer() : Processor("MarlinLumiCalClusterer"),
						   LumiInColName(""),
						   LumiClusterColName(""),
						   LumiRecoParticleColName(""),
						   _rMoliere(16.),
						   _minClusterEngy(2.),
						   _minHitEnergy(5.e-6),
						   _logWeigthConstant(6.),
						   _ElementsPercentInShowerPeakLayer(0.03),
						   _MiddleEnergyHitBoundFrac(0.01),
						   _WeightingMethod(-1),
						   _ClusterMinNumHits(15),
						   _NumOfNearNeighbor(6),
						   SkipNEvents(0),
						   MaxRecordNumber(0),
						   NumRun(0),
						   NumEvt(0),
						   EvtNumber(0),
						   OutDirName(""),
						   NumEventsTree(0),
						   OutputManager(),
						   GlobalMethods(),
						   LumiCalClusterer(LumiInColName)
						   
{
  _description = "whatever..." ;

  //---------------------------------------------------------
  // registerProcessorParameter(	 "LumiCal_Collection" ,
  //					 "some description" ,
  //					 LumiName ,
  //					 std::string("LumiCalCollection") ) ;

  registerInputCollection(LCIO::SIMCALORIMETERHIT,
			  "LumiCal_Collection" ,
			  "Collection Containing the Hits in the LumiCal" ,
			  LumiInColName ,
			  std::string("LumiCalCollection") ) ;

  registerOutputCollection(LCIO::CLUSTER,
			   "LumiCal_Clusters" ,
			   "Collection of Cluster found in the LumiCal" ,
			   LumiClusterColName ,
			   std::string("LumiCalClusters") ) ;

  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE,
			   "LumiCal_RecoParticles" ,
			   "Collection of Reconstructed Particles found in the LumiCal" ,
			   LumiRecoParticleColName ,
			   std::string("LumiCalRecoParticles") ) ;


  //---------------------------------------------------------
  registerProcessorParameter(	"SkipNEvents" ,
				" Number of events to skip at the begining of the LCIO file" ,
				SkipNEvents,
				0 );
  //---------------------------------------------------------
  registerProcessorParameter(	"MaxRecordNumber" ,
				" Number of event to work with" ,
				MaxRecordNumber,
				10 );
  //---------------------------------------------------------
  registerProcessorParameter(	"OutDirName" ,
				"Name of output directory" ,
				OutDirName,
				std::string("rootOut") );
  //---------------------------------------------------------
  // Processor Clustering Parameters
  //---------------------------------------------------------
  registerProcessorParameter(  "LogWeigthConstant",
                               " Sets minimum for logarithmic energy weights",
                               _logWeigthConstant,
			       6. );
  registerProcessorParameter(  "MoliereRadius",
                               " Moliere radius, controls clusters separation distance [mm]",
                               _rMoliere,
                               16. );
  registerProcessorParameter(  "MinClusterEngy",
                               " Sets minimum energy deposit for cluster to be accepted [GeV]",
                               _minClusterEngy,
                               2. );
  registerProcessorParameter(  "WeightingMethod",
                               " Choose cluster hit position weights LogMthod=-1 or simple EnergyMethod=1 ",
			       _WeightingMethod,
                               -1 );
  registerProcessorParameter(  "MinHitEnergy",
                               " Hit energy cut [Gev] ",
			       _minHitEnergy,
                               5.e-6 );
  registerProcessorParameter(  "ClusterMinNumHits",
                               " Minimal number of hits in cluster",
			       _ClusterMinNumHits,
                               15 );
  registerProcessorParameter(  "ElementsPercentInShowerPeakLayer",
                               " BP: Not sure what it is",
			       _ElementsPercentInShowerPeakLayer,
                               0.03 );
  registerProcessorParameter(  "MiddleEnergyHitBoundFrac",
                               " BP: see explanation in LumiCalClusterer.cpp",
			       _MiddleEnergyHitBoundFrac,
                               0.01 );
   registerProcessorParameter(  "NumOfNearNeighbor",
                               "Number of neighbor hits to consider ",
			       _NumOfNearNeighbor,
                                6 );
}



/* ============================================================================
   initial action before first event analysis starts:
   Called at the begining of the job before anything is read.
   ========================================================================= */
void MarlinLumiCalClusterer::init(){

  	printParameters();

  // global vars
  NumRun = 0 ;
  NumEvt = SkipNEvents;


  //  GlobalMethods = new GlobalMethodsClass();
  
  GlobalMethods.SetConstants();

  //LumiCalClusterer = new LumiCalClustererClass(LumiInColName);
  LumiCalClusterer.setLumiCollectionName(LumiInColName);
  LumiCalClusterer.init(GlobalMethods.GlobalParamI, GlobalMethods.GlobalParamD);

  NumEventsTree = 500;
  //OutputManager = new OutputManagerClass();
  OutputManager.Initialize(SkipNEvents , NumEventsTree, OutDirName.c_str());

  /* --------------------------------------------------------------------------
     Print out Processor Parameters
     -------------------------------------------------------------------------- */
  streamlog_out(MESSAGE) << std::endl << "Global parameters for "<< name() << " Processor:" << std::endl;
  GlobalMethods.PrintAllParameters();

  streamlog_out(MESSAGE) << std::endl;

}


/* ============================================================================
   pre run action:
   Called for every run, e.g.
   ========================================================================= */
void MarlinLumiCalClusterer::processRunHeader( LCRunHeader * /* run */ ) {
  NumRun++ ;
}

/* ============================================================================
   main actions in each event:
   Called for every event - the working horse.
   ========================================================================= */
void MarlinLumiCalClusterer::processEvent( EVENT::LCEvent * evt ) {

  // increment / initialize global variables
  NumEvt++;
  EvtNumber = evt->getEventNumber();

  streamlog_out( DEBUG ) << std::endl
	    << "Run MarlinLumiCalClusterer::processEvent - event counter: NumEvt = " << NumEvt
	    << " ( event index " << EvtNumber << " )"
	    << std::endl ;


  OutputManager.NumEventsTree = 500;	
  TryMarlinLumiCalClusterer(evt);

}


/* ============================================================================
   final action after last event analysis is over:
   Called after data processing for clean up in the inverse order of the
   init() method so that resources allocated in the first processor also will
   be available for all following processors.
   ========================================================================= */
void MarlinLumiCalClusterer::end(){

  streamlog_out( DEBUG ) << "Run MarlinLumiCalClusterer::end() "
			 << "Went through " << NumEvt << " events from " << NumRun << " file(s)" 
			 << std::endl;

  // write to the root tree
  OutputManager.WriteToRootTree("forceWrite" , NumEvt);


  // write out the counter map
  int numCounters = OutputManager.Counter.size();

  if(numCounters > 0)
    std::cout << std::endl << "Global counters:"  << std::endl;

  OutputManager.CounterIterator = OutputManager.Counter.begin();
  for(int hisNow = 0; hisNow < numCounters; hisNow++ , OutputManager.CounterIterator++) {
    std::string counterName = (std::string)(*OutputManager.CounterIterator).first;
    std::cout << "\t" << OutputManager.Counter[counterName] << "  \t <->  " << counterName << std::endl;
  }

  OutputManager.CleanUp();

}
