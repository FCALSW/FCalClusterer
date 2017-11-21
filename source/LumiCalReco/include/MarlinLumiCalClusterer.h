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

namespace IMPL {
  class ClusterImpl;
  class ReconstructedParticleImpl;
}

typedef std::map < int , ClusterClass * > MapIntPClusterClass;
typedef std::map < int , std::vector<int> >  MapIntVInt;

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

    // actions after each event - place to produce control plots/printouts ....
    virtual void check( EVENT::LCEvent * evt ) ;

    // final action after last event analysis is over.
    virtual void end() ;

    // Processor Parameters
    std::string	LumiInColName;
    std::string LumiOutColName = "LumiCalHits";
    std::string LumiClusterColName;
    std::string LumiRecoParticleColName ;
    double _BeamCrossingAngle;
    double _zLayerPhiOffset;
    double _rMoliere,_minClusterEngy, _minHitEnergy, _logWeigthConstant;
    double _ElementsPercentInShowerPeakLayer, _MiddleEnergyHitBoundFrac;
    double _EnergyCalibConst;
    std::string _WeightingMethod;
    int     _ClusterMinNumHits, _NumOfNearNeighbor;
         

    // global counters
    int	SkipNEvents, MaxRecordNumber;
    int	NumRun, NumEvt, EvtNumber;
    std::string	OutDirName;
    std::string	OutRootFileName;
    int	NumEventsTree;
    int MemoryResidentTree;

    OutputManagerClass	OutputManager;
    GlobalMethodsClass	gmc;
    LumiCalClustererClass	LumiCalClusterer;
    bool _cutOnFiducialVolume=false;

    void TryMarlinLumiCalClusterer(EVENT::LCEvent * evt);

    void CreateClusters( std::map < int , std::map < int , std::vector<int> > > const& clusterIdToCellId,
			 std::map < int , std::map < int , std::vector<double> > > const& cellIdToCellEngy,
			 //			 std::map < int , std::map < int , LCCluster > > const& superClusterIdClusterInfo,
			 std::map < int , MapIntPClusterClass > & clusterClassMapP,
			 EVENT::LCEvent * evt);
    std::tuple<ClusterImpl*, ReconstructedParticleImpl*> getLCIOObjects(LCCluster const& clusterInfo) const;
    void writeRootInfo(EVENT::LCEvent* evt);

    inline double sqr( double a){ return a*a;};
    inline float sqr( float a){ return a*a;};
    inline int sqr( int a){ return a*a;};

    void storeMCParticleInfo( LCEvent *evt, int clusterInFlag  );

  };

#endif

