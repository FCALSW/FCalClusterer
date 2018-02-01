#include "ReadBeamCal.hh"
#include "ProcessorUtilities.hh"

#include "BCPadEnergies.hh"
#include "BCUtilities.hh"
#include "BeamCal.hh"
#include "BeamCalGeo.hh"

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCIO.h>
#include <EVENT/SimCalorimeterHit.h>
#include <UTIL/CellIDDecoder.h>

#include <Exceptions.h>

// ----- include for verbosity dependend logging ---------
#include <streamlog/baselevels.h>
#include <streamlog/loglevels.h>
#include <streamlog/logstream.h>
#include <streamlog/streamlog.h>

#include <marlin/ProcessorEventSeeder.h>
#include <marlin/Global.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TH2.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TTree.h>

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace lcio ;
using namespace marlin ;

ReadBeamCal aReadBeamCal ;

// #pragma GCC diagnostic push
// #pragma GCC diagnostic ignored "-Wshadow"

ReadBeamCal::ReadBeamCal()
    : Processor("ReadBeamCal"),
      m_colNameBCal(""),
      m_nameOutputFile(""),
      m_nameFinalOutputFile(""),
      m_nameInputFile(""),
      m_nRun(0),
      m_nEvt(0),
      m_probFactor(0.0),
      m_random3(nullptr),
      m_padEnergiesLeft(nullptr),
      m_padEnergiesRight(nullptr),
      m_bcg(nullptr),
      m_usingDD4HEP(false) {
  // modify processor description
  _description = "ReadBeamCal reads the simulation for the pairs and creates two std::vector<double> in a tree, which can then be used later on for Overlay, calculation of fluctiuations, etc." ;


  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::SIMCALORIMETERHIT,
			   "BeamCalCollectionName" ,
			   "Name of BeamCal Collection"  ,
			   m_colNameBCal ,
			   std::string("BeamCalCollection") ) ;

  registerProcessorParameter ("OutputFileBackground",
			      "Root OutputFile ",
			      m_nameOutputFile,
			      std::string("BeamCal.root") ) ;

  registerProcessorParameter ("ProbabilityFactor",
			      "Probability [0,100%] for particle to be added to event. Allows the scaling of"\
			      "the background to smaller background rate",
			      m_probFactor,
			      double(100.0) ) ;

}
//#pragma GCC diagnostic pop

void ReadBeamCal::init() {

  Global::EVENTSEEDER->registerProcessor(this);
  m_random3 = new TRandom3();

  // usually a good idea to
  printParameters() ;

  m_nEvt = 0;
  m_bcg = ProcessorUtilities::getBeamCalGeo(m_usingDD4HEP);
  m_padEnergiesLeft = new BCPadEnergies(m_bcg);
  m_padEnergiesRight = new BCPadEnergies(m_bcg);
 
}//init

void ReadBeamCal::processRunHeader( LCRunHeader* ) {
  //  streamlog_out (DEBUG) << "Runnumber "<< _nRun << std::endl;
  //   if(_nRun % 4 == 0) {
}

void ReadBeamCal::processEvent( LCEvent * evt ) {
  m_random3->SetSeed(Global::EVENTSEEDER->getSeed(this));
  
  if ( m_random3->Uniform() > m_probFactor/100. ) return;
 
  LCCollection  *colBCal;

  try {
    colBCal = evt->getCollection( m_colNameBCal ) ;
  } catch (Exception &e) {
    colBCal = 0;
  }
  m_nEvt ++ ;
  if( not colBCal ) return;

  CellIDDecoder<SimCalorimeterHit> mydecoder(colBCal);
  int nHits = colBCal->getNumberOfElements();
  for(int i=0; i < nHits; i++) {
    SimCalorimeterHit *bcalhit = static_cast<SimCalorimeterHit*>(colBCal->getElementAt(i));
    int side, layer, cylinder, sector;
    BCUtil::DecodeCellID(mydecoder, bcalhit, side, layer, cylinder, sector, m_usingDD4HEP);
    const float energy = bcalhit->getEnergy();
    try{
      (side == 0 ) ?
	m_padEnergiesLeft->addEnergy(layer, cylinder, sector, energy):
	m_padEnergiesRight->addEnergy(layer, cylinder, sector, energy);
    } catch (std::out_of_range &e) {
      streamlog_out( DEBUG4 ) 
	<< "Caught exception  " <<  e.what()
	<< std::setw(10) << layer
	<< std::setw(10) << cylinder
	<< std::setw(10) << sector
	<< std::setw(20) << bcalhit->getPosition()[0]
	<< std::setw(20) << bcalhit->getPosition()[1]
	<< std::setw(20) << bcalhit->getPosition()[2]
	<< std::setw(20) << energy
	<< std::endl;
    }
  }//for all entries in the collection
    
  return;
}//processEvent



void ReadBeamCal::check( LCEvent * ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ReadBeamCal::end(){

  streamlog_out ( MESSAGE ) << __PRETTY_FUNCTION__ << " " << name()
			    << " processed " << m_nEvt << " events."
			    << std::endl ;
  //Do the average for every bin, and calculate the maximal difference to the mean, which means, we have to loop twice.

  TTree *tree = new TTree("bcTree","bcTree");
  // tree->Branch("h3BC_left",h3BeamCalDeposits_left);
  // tree->Branch("h3BC_right",h3BeamCalDeposits_right);
  tree->Branch("vec_right",m_padEnergiesRight->getEnergies());
  tree->Branch("vec_left",m_padEnergiesLeft->getEnergies());
  tree->Fill();

  TFile* rootfile = TFile::Open(m_nameOutputFile.c_str(), "RECREATE");

  // for(unsigned int k = 0; k < bgruns.size(); k++) {
  //   bgruns[k]->Write();
  //   //     bgpos[k]->Write();
  // }

  if( streamlog::out.write< DEBUG >() ) {

    BeamCal bc( *m_bcg );

    gStyle->SetOptStat(0);
    bc.SetLogz(1);
    bc.SetAxisMax(10);
    bc.SetAxisMin(1e-3);
    bc.SetBeamCalHisto(m_padEnergiesRight);
    TCanvas c1("c1","c1");
    c1.SetRightMargin(0.16);
    TH2F frame("frame","BeamCal",10,-160,160,10,-160,160);
    bc.BeamCalDraw(&c1, &frame);
    bc.SetBeamCalHisto(m_padEnergiesLeft);
    bc.BeamCalDraw(&c1, &frame);

  }//Only draw all the things in DEBUG Mode


  tree->Write();
  // background.Write();
  // fluctuation.Write();
  rootfile->Write();
  rootfile->Close();
  delete rootfile;
  rootfile = nullptr;
  delete m_random3;
  delete m_bcg;

}
