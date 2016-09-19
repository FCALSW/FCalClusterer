#ifndef DrawBeamCalFromDD4hep_h
#define DrawBeamCalFromDD4hep_h 1

#include <BCPadEnergies.hh>

#include <marlin/Processor.h>
#include <lcio.h>

#include <DD4hep/LCDD.h>
#include <DD4hep/Detector.h>

#include <string>
#include <vector>
#include <map>


class TH1D;
class TH2D;
class TH3D;
class TFile;
class TRandom3;
class TTree;

typedef std::map< unsigned int, double> MapIdVal;

class DrawBeamCalFromDD4hep : public marlin::Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new DrawBeamCalFromDD4hep ; }
  
  
  DrawBeamCalFromDD4hep() ;

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

  BeamCalGeo* m_bcg;
  DD4hep::Geometry::DetElement m_BeamCal;
  DD4hep::Geometry::Segmentation m_seg;
  TFile * m_file;
  TTree * m_tree;
  double m_x, m_y, m_z, m_energy;
private://to shut the warnings up
  DrawBeamCalFromDD4hep(const DrawBeamCalFromDD4hep&);
  DrawBeamCalFromDD4hep& operator=(const DrawBeamCalFromDD4hep&);

  MapIdVal hitEnergies;

} ;


#endif



