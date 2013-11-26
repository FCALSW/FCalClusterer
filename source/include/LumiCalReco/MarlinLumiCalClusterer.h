// check if this file has already been called:
#ifndef MarlinLumiCalClusterer_
#define MarlinLumiCalClusterer_

#include "Global.hh"

// Marlin classes:
#include <marlin/Processor.h>
#include <marlin/Global.h>
#include <lcio.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/SimCalorimeterHitImpl.h>
#include <IMPL/ClusterImpl.h>

namespace EVENT{
  class LCEvent;
}

// gear
#include <gear/GEAR.h>
#include <gear/CalorimeterParameters.h>
#include <gear/LayerLayout.h>


// // cout color definitios
// #define coutDefault          "\033[0m"
// #define coutRed                      "\033[31m"
// #define coutGreen            "\033[32m"
// #define coutBlue             "\033[34m"
// #define coutPurple           "\033[35m"
// #define coutUnderLine                "\033[4;30m"
// #define coutWhiteOnBlack     "\33[40;37;1m"


// my standard includes:
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

// ROOT libs
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TColor.h>

#include "GlobalMethodsClass.h"
#include "ClusterClass.h"
#include "SortingClass.h"
#include "OutputManagerClass.h"

#include "LumiCalClusterer.h"



  class MarlinLumiCalClusterer : public marlin::Processor {

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
    virtual void processEvent( EVENT::LCEvent * evt ) ;

    // final action after last event analysis is over.
    virtual void end() ;

    //protected:

    // Processor Parameters
    std::string	LumiInColName, LumiOutColName ;
    int	SkipNEvents, MaxRecordNumber;

    // global counters
    int	NumRun, NumEvt, EvtNumber;
    std::string	OutDirName;

    int	NumEventsTree;

    OutputManagerClass	OutputManager;
    GlobalMethodsClass	GlobalMethods;
    LumiCalClustererClass	LumiCalClusterer;

    void TryMarlinLumiCalClusterer(EVENT::LCEvent * evt);

    void CreateClusters( std::map < int , std::map < int , std::vector<int> > > const& clusterIdToCellId,
			 std::map < int , std::map < int , std::vector<double> > > const& cellIdToCellEngy,
			 std::map < int , std::map < int , ClusterClass * > > & clusterClassMapP );


  };

#endif
