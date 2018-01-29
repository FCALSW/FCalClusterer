#ifndef DrawBeamCalFromDD4hep_h
#define DrawBeamCalFromDD4hep_h 1

#include <marlin/Processor.h>

#include <DD4hep/DetElement.h>
#include <DD4hep/Segmentations.h>

#include <string>
#include <map>

class TFile;
class TTree;

namespace EVENT {
  class LCEvent;
  class LCRunHeader;
}

typedef std::map< unsigned int, double> MapIdVal;

class DrawBeamCalFromDD4hep : public marlin::Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new DrawBeamCalFromDD4hep ; }
  
  
  DrawBeamCalFromDD4hep() ;
  DrawBeamCalFromDD4hep(const DrawBeamCalFromDD4hep&) = delete ;
  DrawBeamCalFromDD4hep & operator = (const DrawBeamCalFromDD4hep&) = delete ;

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
  

  void drawSegmentation();
  void drawPolarGridRPhi2();
  void drawCartesianGridXY();

 protected:

  /** Input collection name.
   */
  std::string m_colNameBCal, m_nameOutputFile, m_nameFinalOutputFile, m_nameInputFile;
  //std::string m_pdftitle;
  int m_nRun ;
  int m_nEvt ;
  bool m_drawDensities;

  dd4hep::DetElement m_BeamCal{};
  dd4hep::Segmentation m_seg{};
  TFile * m_file;
  TTree * m_tree;
  double m_x, m_y, m_z, m_energy;
private://to shut the warnings up

  MapIdVal hitEnergies{};

} ;


#endif



