#ifndef MarlinLumiCalClusterer_
#define MarlinLumiCalClusterer_
#include <marlin/Global.h>
#include "Global.hh"

#include "GlobalMethodsClass.h"
#include "LumiCalClusterer.h"
#include "OutputManagerClass.h"

class ClusterClass;

// Marlin classes:
#include <marlin/Processor.h>

namespace EVENT{
  class LCEvent;
}


// // cout color definitios
// #define coutDefault          "\033[0m"
// #define coutRed                      "\033[31m"
// #define coutGreen            "\033[32m"
// #define coutBlue             "\033[34m"
// #define coutPurple           "\033[35m"
// #define coutUnderLine                "\033[4;30m"
// #define coutWhiteOnBlack     "\33[40;37;1m"




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

    // Processor Parameters
    std::string	LumiInColName, LumiClusterColName, LumiRecoParticleColName ;
    double _rMoliere,_minClusterEngy, _minHitEnergy, _logWeigthConstant;
    double _ElementsPercentInShowerPeakLayer, _MiddleEnergyHitBoundFrac;
    int    _WeightingMethod, _ClusterMinNumHits, _NumOfNearNeighbor;
         

    // global counters
    int	SkipNEvents, MaxRecordNumber;
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
