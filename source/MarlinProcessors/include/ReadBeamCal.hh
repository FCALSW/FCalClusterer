#ifndef ReadBeamCal_h
#define ReadBeamCal_h 1

#include <string>

#include <marlin/Processor.h>

class BCPadEnergies;
class BeamCalGeo;

class TRandom3;

namespace EVENT {
  class LCEvent;
  class LCRunHeader;
}

class ReadBeamCal : public marlin::Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new ReadBeamCal ; }
  
  
  ReadBeamCal() ;

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
  int m_nRun ;
  int m_nEvt ;
  double m_probFactor;

  TRandom3 *m_random3;

  BCPadEnergies* m_padEnergiesLeft;
  BCPadEnergies* m_padEnergiesRight;

  BeamCalGeo* m_bcg;
  bool m_usingDD4HEP;

private://to shut the warnings up
  ReadBeamCal(const ReadBeamCal&);
  ReadBeamCal& operator=(const ReadBeamCal&);

} ;


#endif



