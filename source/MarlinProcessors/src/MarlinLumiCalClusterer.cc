#include "MarlinLumiCalClusterer.h"

#include <EVENT/LCEvent.h>
#include <EVENT/LCIO.h>

#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>

#include <iostream>
#include <map>
#include <string>
#include <utility>

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
						   _BeamCrossingAngle(0.0),
						   _zLayerPhiOffset(3.75),
						   _rMoliere(16.),
						   _minClusterEngy(2.),
						   _minHitEnergy(5.e-6),
						   _logWeigthConstant(6.),
						   _ElementsPercentInShowerPeakLayer(0.03),
						   _MiddleEnergyHitBoundFrac(0.01),
						   _EnergyCalibConst(0.0105),
						   _WeightingMethod("LogMethod"),
						   _ClusterMinNumHits(15),
						   _NumOfNearNeighbor(6),
						   SkipNEvents(0),
						   MaxRecordNumber(0),
						   NumRun(0),
						   NumEvt(0),
						   EvtNumber(0),
						   OutDirName(""),
						   OutRootFileName(""),
						   NumEventsTree(0),
						   MemoryResidentTree(0),
						   OutputManager(),
						   gmc(),
                                                   LumiCalClusterer(LumiInColName)
						   
{
  _description = "Reconstruction of clusters in the LumiCal detector" ;

  //---------------------------------------------------------
  // registerProcessorParameter(	 "LumiCal_Collection" ,
  //					 "some description" ,
  //					 LumiName ,
  //					 std::string("LumiCalCollection") ) ;

  registerInputCollection(LCIO::SIMCALORIMETERHIT,
			  "LumiCal_Collection" ,
			  "Collection Containing Hits in the LumiCal" ,
			  LumiInColName ,
			  std::string("LumiCalCollection") ) ;

  registerOutputCollection(LCIO::CALORIMETERHIT, "LumiCal_Hits", "Collection of CalorimtersHits from the LumiCal",
                           LumiOutColName, LumiOutColName);

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
   registerProcessorParameter(	"OutRootFileName" ,
				"Name of output ROOT file ( without suffix) NO DEFAULT. If no name provided root file will not be generated" ,
				OutRootFileName,
				std::string("") );
  registerProcessorParameter(   "NumEventsTree",
				"Number of events in memory resident ROOT tree.",
				NumEventsTree,
				500 );
  registerProcessorParameter(   "MemoryResidentTree",
				" Place for ROOT tree memory(1) or disk(0)",
				MemoryResidentTree,
				0 );
  //---------------------------------------------------------
  // Processor Clustering Parameters
  //---------------------------------------------------------
  registerProcessorParameter(  "LogWeigthConstant",
                               " Sets minimum for logarithmic energy weights",
                               _logWeigthConstant,
			       6. );
  registerProcessorParameter(  "ZLayerPhiOffset",
                               " Relative offset of LCal z-layers [deg] default is half of the phi sector size",
                               _zLayerPhiOffset,
                               3.75 );
  registerProcessorParameter(  "EnergyCalibConst",
                               " Calibration const E_dep = EnergyCalibConst*E_primary ( default for LCal ILD) ",
                               _EnergyCalibConst,
                               0.0105 );
  registerProcessorParameter(  "MoliereRadius",
                               " Moliere radius, controls clusters separation distance [mm]",
                               _rMoliere,
                               16. );
  registerProcessorParameter(  "MinClusterEngy",
                               " Sets minimum energy deposit for cluster to be accepted [GeV]",
                               _minClusterEngy,
                               2. );
  registerProcessorParameter(  "WeightingMethod",
                               " Hit positions weighting method (LogMthod, EnergyMethod ) ",
			       _WeightingMethod,
                               std::string("LogMethod") );
  registerProcessorParameter(  "MinHitEnergy",
                               " Hit energy cut [Mev] ",
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
  registerProcessorParameter(  "CutOnFiducuialVolume",
                               "Whether to cut clusters outside of the fiducial volume or not",
                               _cutOnFiducialVolume,
                               false );
}



/* ============================================================================
   initial action before first event analysis starts:
   Called at the begining of the job before anything is read.
   ========================================================================= */
void MarlinLumiCalClusterer::init(){


  // global vars
  NumRun = 0 ;
  NumEvt = SkipNEvents;
  //  LumiInColName = "LumiCalSimHits";

  gmc.SetConstants( this );

  _BeamCrossingAngle = gmc.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle]/2.;

  printParameters();
  /* --------------------------------------------------------------------------
     Print out Processor Parameters
     -------------------------------------------------------------------------- */
  streamlog_out(MESSAGE) << std::endl << "Global parameters for Processor:"<< name() <<"\t"<< type() << std::endl;
  gmc.PrintAllParameters();
  streamlog_out(MESSAGE) << std::endl;

  //LumiCalClusterer = new LumiCalClustererClass(LumiInColName);
  LumiCalClusterer.setLumiCollectionName(LumiInColName);
  LumiCalClusterer.setLumiOutCollectionName(LumiOutColName);
  LumiCalClusterer.init( gmc );
  LumiCalClusterer.setCutOnFiducialVolume(_cutOnFiducialVolume);

  //OutputManager = new OutputManagerClass();
  OutputManager.Initialize(MemoryResidentTree, SkipNEvents , NumEventsTree, OutDirName, OutRootFileName);



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
  NumEvt++;                             // number of processed events 
  EvtNumber = evt->getEventNumber();    // event number 

  streamlog_out( DEBUG ) << std::endl
	    << "Run MarlinLumiCalClusterer::processEvent - event counter: NumEvt = " << NumEvt
	    << " ( event index " << EvtNumber << " )"
	    << std::endl ;


  //  OutputManager.NumEventsTree = 500;	
  TryMarlinLumiCalClusterer(evt);

}

void MarlinLumiCalClusterer::check( EVENT::LCEvent * /* evt */ ) {
  /*----  NOP for now ---- */
}

/* ============================================================================
   final action after last event analysis is over:
   Called after data processing for clean up in the inverse order of the
   init() method so that resources allocated in the first processor also will
   be available for all following processors.
   ========================================================================= */
void MarlinLumiCalClusterer::end(){

  streamlog_out( MESSAGE ) << "Run MarlinLumiCalClusterer::end() "
			 << "Went through " << NumEvt << " events from " << NumRun << " file(s)" 
			 << std::endl;



  // write out the counter map
  int numCounters = OutputManager.Counter.size();

  if(numCounters > 0)
    std::cout << std::endl << "Global counters:"  << std::endl;

  OutputManager.CounterIterator = OutputManager.Counter.begin();
  for(int hisNow = 0; hisNow < numCounters; hisNow++ , OutputManager.CounterIterator++) {
    std::string counterName = (std::string)(*OutputManager.CounterIterator).first;
    std::cout << "\t" << OutputManager.Counter[counterName] << "  \t <->  " << counterName << std::endl;
  }

  // write to the root tree
  OutputManager.WriteToRootTree("forceWrite" , NumEvt);

  // (BP) It is done by destructor!
  //OutputManager.CleanUp();       

}

