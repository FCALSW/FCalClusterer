#ifndef ReadBeamCal_h
#define ReadBeamCal_h 1

#include <BCPadEnergies.hh>

#include <string>
#include <vector>

#include <marlin/Processor.h>
#include <lcio.h>


class TH1D;
class TH2D;
class TH3D;
class TFile;
class TRandom3;
class TTree;


class ReadBC_tower : public marlin::Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new ReadBC_tower ; }
  
  
  ReadBC_tower() ;
  ~ReadBC_tower() {};

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

  /** Input collection name.
   */
  std::string m_colNameBCal, m_nameOutputFile, m_nameFinalOutputFile, m_nameInputFile;
  //std::string m_pdftitle;
  int m_nRun, m_nEvt ;
  double m_probFactor, m_threshold;

  TRandom3 *m_random3;

  BCPadEnergies* m_padEnergiesLeft;
  BCPadEnergies* m_padEnergiesRight;
  TTree *tree;

  BeamCalGeo* m_bcg;
  BCPCuts* m_bcpCuts;

private://to shut the warnings up
  ReadBC_tower(const ReadBC_tower&);
  ReadBC_tower& operator=(const ReadBC_tower&);

} ;


#endif



