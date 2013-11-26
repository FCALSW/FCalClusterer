#define _LC_DEBUG

// include my header file + whatever LCIO classes i use below

#include "MarlinLumiCalClusterer.h"
#include "MarlinLumiCalClusterer_auxiliary.h"

#include <EVENT/LCEvent.h>


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
						   LumiOutColName(""),
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
			   LumiOutColName ,
			   std::string("LumiCalClusters") ) ;

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

}



/* ============================================================================
   initial action before first event analysis starts:
   Called at the begining of the job before anything is read.
   ========================================================================= */
void MarlinLumiCalClusterer::init(){

  //	printParameters();

  // global vars
  NumRun = 0 ;
  NumEvt = SkipNEvents;

  //	GlobalMethods = new GlobalMethodsClass();
  GlobalMethods.SetConstants();

  //LumiCalClusterer = new LumiCalClustererClass(LumiInColName);
  LumiCalClusterer.init(GlobalMethods.GlobalParamI, GlobalMethods.GlobalParamD);

  NumEventsTree = 500;
  //OutputManager = new OutputManagerClass();
  OutputManager.Initialize(SkipNEvents , NumEventsTree, OutDirName.c_str());

  /* --------------------------------------------------------------------------
     Print out Processor Parameters
     -------------------------------------------------------------------------- */
  std::cout << std::endl << "Global parameters for MarlinLumiCalClusterer Processor:" << std::endl;
  GlobalMethods.PrintAllParameters();

  std::cout << std::endl;

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

#ifdef _LC_DEBUG
  std::cout << std::endl
	    << "Run MarlinLumiCalClusterer::processEvent - event counter: NumEvt = " << NumEvt
	    << " ( event index " << EvtNumber << " )"
	    << std::endl ;
#endif


  OutputManager.NumEventsTree = 500;	TryMarlinLumiCalClusterer(evt);

}


/* ============================================================================
   final action after last event analysis is over:
   Called after data processing for clean up in the inverse order of the
   init() method so that resources allocated in the first processor also will
   be available for all following processors.
   ========================================================================= */
void MarlinLumiCalClusterer::end(){

#ifdef _LC_DEBUG
  std::cout << std::endl << "Run MarlinLumiCalClusterer::end() "
	    << std::endl << std::endl
	    << "Went through " << NumEvt << " events from " << NumRun << " file(s)" << std::endl << std::endl;
#endif

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

}
