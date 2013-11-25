#define _LC_DEBUG


// include my header file + whatever LCIO classes i use below
#include "MarlinLumiCalClusterer.h"

#include "LumiCalClusterer.h"

#include "MarlinLumiCalClusterer_auxiliary.h"



namespace MarlinLumiCalClusterer {
	
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
MarlinLumiCalClusterer::MarlinLumiCalClusterer() : Processor("MarlinLumiCalClusterer") {
	_description = "whatever..." ;

//---------------------------------------------------------
	registerProcessorParameter(	 "LumiCal_Collection" ,
                                   	 "some description" ,
                                   	 LumiName ,
                                   	 string("LumiCalCollection") ) ;
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
					string("rootOut") );
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

	GlobalMethods = new GlobalMethodsClass();
	GlobalMethods -> SetConstants();

	LumiCalClusterer = new LumiCalClustererClass(LumiName);
	LumiCalClusterer -> init(GlobalMethods-> GlobalParamI, GlobalMethods->GlobalParamD);
	
	NumEventsTree = 500;
	OutputManager = new OutputManagerClass();
	OutputManager -> Initialize(SkipNEvents , NumEventsTree, TString(OutDirName.c_str()));
	
	/* --------------------------------------------------------------------------
	   Print out Processor Parameters
	-------------------------------------------------------------------------- */
	map < TString , int > ::iterator	GlobalParamIIterator;
	map < TString , double > ::iterator	GlobalParamDIterator;
	map < TString , TString > ::iterator	GlobalParamSIterator;
	int numParam;

	cout	<< endl << coutUnderLine << coutBlue << "Global parameters for MarlinLumiCalClusterer Processor:" << coutDefault << endl;

	GlobalParamIIterator = GlobalMethods->GlobalParamI.begin();
	numParam             = GlobalMethods->GlobalParamI.size();
	for(int paramNow = 0; paramNow < numParam; paramNow++, GlobalParamIIterator++){
		TString paramName = (TString)(*GlobalParamIIterator).first;
		int	paramVal  = (int)(*GlobalParamIIterator).second;

		cout	<< " - (int)     " << paramName << "  =  " << paramVal << endl;
	}

	GlobalParamDIterator = GlobalMethods->GlobalParamD.begin();
	numParam             = GlobalMethods->GlobalParamD.size();
	for(int paramNow = 0; paramNow < numParam; paramNow++, GlobalParamDIterator++){
		TString paramName = (TString)(*GlobalParamDIterator).first;
		double	paramVal  = (double)(*GlobalParamDIterator).second;

		cout	<< " - (double)  " << paramName << "  =  " << paramVal << endl;
	}

	GlobalParamSIterator = GlobalMethods->GlobalParamS.begin();
	numParam             = GlobalMethods->GlobalParamS.size();
	for(int paramNow = 0; paramNow < numParam; paramNow++, GlobalParamSIterator++){
		TString paramName = (TString)(*GlobalParamSIterator).first;
		TString	paramVal  = (TString)(*GlobalParamSIterator).second;

		cout	<< " - (TString) " << paramName << "  =  " << paramVal << endl;
	}


	cout	<< endl;

}


/* ============================================================================
   pre run action:
   Called for every run, e.g. 
========================================================================= */
void MarlinLumiCalClusterer::processRunHeader( LCRunHeader * run ) {
	NumRun++ ;
}

/* ============================================================================
   main actions in each event:
   Called for every event - the working horse. 
========================================================================= */
void MarlinLumiCalClusterer::processEvent( LCEvent * evt ) {

	// increment / initialize global variables
	NumEvt++;
	EvtNumber = evt->getEventNumber();

	#ifdef _LC_DEBUG
	cout	<< endl << coutWhiteOnBlack
		<< "Run MarlinLumiCalClusterer::processEvent - event counter: NumEvt = " << NumEvt << " ( event index " << EvtNumber << " )"
		<< coutDefault << endl ;
	#endif


	OutputManager -> NumEventsTree = 500;	TryMarlinLumiCalClusterer(evt);

}


/* ============================================================================
   final action after last event analysis is over:
   Called after data processing for clean up in the inverse order of the
   init() method so that resources allocated in the first processor also will
   be available for all following processors. 
========================================================================= */
void MarlinLumiCalClusterer::end(){

	#ifdef _LC_DEBUG
	cout	<< endl << coutWhiteOnBlack << "Run MarlinLumiCalClusterer::end() "
		<< coutDefault << endl << endl
		<< "Went through " << NumEvt << " events from " << NumRun << " file(s)" << endl << endl;	
	#endif

	// write to the root tree
	OutputManager->WriteToRootTree("forceWrite" , NumEvt);

	
	// write out the counter map
	int numCounters = OutputManager->Counter.size();
	
	if(numCounters > 0)
		cout 	<< endl << coutRed << "Global counters:"  << coutDefault << endl;

	OutputManager->CounterIterator = OutputManager->Counter.begin();
	for(int hisNow = 0; hisNow < numCounters; hisNow++ , OutputManager->CounterIterator++) {
		TString counterName = (TString)(*OutputManager->CounterIterator).first;
			cout << "\t" << OutputManager->Counter[counterName] << "  \t <->  " << counterName << endl;
	}

	delete	GlobalMethods;
	delete	LumiCalClusterer;
	delete	OutputManager;
}


} // namespace MarlinLumiCalClusterer
