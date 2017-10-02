#include "DrawBeamCalFromDD4hep.hh"

#include <BeamCal.hh>
#include <BeamCalGeoGear.hh>
#include <BeamCalGeoCached.hh>
#include <BCUtilities.hh>

#ifdef FCAL_USE_CustomRoot
#include <CustomRoot.hh>
#endif


#include "DD4hep/DD4hepUnits.h"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include <UTIL/CellIDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>

#include <marlin/ProcessorEventSeeder.h>
#include <marlin/Global.h>

#include <TTree.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TPaletteAxis.h>
#include <TCrown.h>

#include <iostream>
#include <iomanip>


using namespace lcio ;
using namespace marlin ;

DrawBeamCalFromDD4hep aDrawBeamCalFromDD4hep;

// #pragma GCC diagnostic push
// #pragma GCC diagnostic ignored "-Wshadow"


int getColor(double energy) {
	Double_t minenergy=1e-3, maxenergy=1e2;
	Int_t color1, color2, color;
	double wlmin = TMath::Log10(minenergy);
	double wlmax = TMath::Log10(maxenergy);
	Int_t ncolors = gStyle->GetNumberOfColors();
	Int_t ndivz   = 99;
	Double_t scale = ndivz/(wlmax - wlmin);
	double zc = TMath::Log10(energy);
	if (zc < wlmin) zc = wlmin;
	if (zc > wlmax) zc = wlmax;
	color1 = Int_t(0.01+(zc-wlmin)*scale);
	if(color1 > 98) color1 = 98;
	color2 = Int_t((color1+0.99)*Double_t(ncolors)/Double_t(ndivz));
	color = gStyle->GetColorPalette(color2);

	return color;

}

DrawBeamCalFromDD4hep::DrawBeamCalFromDD4hep() : Processor("DrawBeamCalFromDD4hep"),
			     m_colNameBCal(""),
			     m_nameOutputFile(""),
			     m_nameFinalOutputFile(""),
			     m_nameInputFile(""),
			     m_nRun(0),
			     m_nEvt(0),
			     m_bcg(NULL),
						 m_file(NULL),
						 m_tree(NULL),
						 m_x(0.0),
						 m_y(0.0),
						 m_z(0.0),
						 m_energy(0.0)

						 
{

  // modify processor description
  _description = "DrawBeamCalFromDD4hep draws segmentation from dd4hep with energy deposits in the lcio files" ;


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

}
//#pragma GCC diagnostic pop

void DrawBeamCalFromDD4hep::init() {


  // usually a good idea to
  printParameters() ;

  m_nEvt = 0;

  m_file = TFile::Open(m_nameOutputFile.c_str(), "RECREATE");

  m_tree = new TTree("hitBC","hitBC");
  m_tree->Branch("x",&m_x);
  m_tree->Branch("y",&m_y);
  m_tree->Branch("z",&m_z);
  m_tree->Branch("e",&m_energy);


  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();
  
  m_BeamCal = theDetector.detector("BeamCal");

  //  dd4hep::PlacedVolume place = m_BeamCal.placement();
  dd4hep::Position local(4.0, 0.0, 0.0);
  dd4hep::Position global(0.0, 0.0, 0.0);
  m_BeamCal.child("BeamCal01").nominal().localToWorld(local, global);
  std::cout << global  << std::endl;
  std::cout << local  << std::endl;
  dd4hep::DetElement::Children children = m_BeamCal.children();
  std::cout << children.size()  << std::endl;

  for (dd4hep::DetElement::Children::iterator it = children.begin(); it != children.end(); ++it) {
    std::cout << it->first  << std::endl;
  }
 
  m_BeamCal.child("BeamCal01").nominal().worldTransformation().Print();

  dd4hep::Detector::HandleMap readouts = theDetector.readouts();
  for (dd4hep::Detector::HandleMap::iterator it = readouts.begin(); it != readouts.end(); ++it) {
    std::cout << it->first  << std::endl;
  }

  dd4hep::Readout myReadout = theDetector.readout("BeamCalHits");
  m_seg = myReadout.segmentation();
  std::cout << m_seg.type()  << std::endl;
  std::cout <<"FieldDef: " << m_seg.segmentation()->fieldDescription()  << std::endl;

  dd4hep::DDSegmentation::Parameters pars = m_seg.parameters();
  for (dd4hep::DDSegmentation::Parameters::iterator it = pars.begin(); it != pars.end(); ++it) {
    std::cout << (*it)->toString()
      //	      << " " << (*it)->value()  
	      << std::endl;
  }

  dd4hep::DDSegmentation::CellID cid = m_seg.cellID(local, global, 0);
  
  std::cout << "CellID:" <<  cid  << std::endl;
  std::cout <<  "Some position: " << m_seg.position(cid)  << std::endl;

}//init

void DrawBeamCalFromDD4hep::processRunHeader( LCRunHeader* ) {
  //  streamlog_out (DEBUG) << "Runnumber "<< _nRun << std::endl;
  //   if(_nRun % 4 == 0) {
}

void DrawBeamCalFromDD4hep::processEvent( LCEvent * evt ) {
 
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
    int // side,
      layer, cylinder, sector;
    //BCUtil::DecodeCellID(mydecoder, bcalhit, side, layer, cylinder, sector);

    // std::cout << "PosX: "<< bcalhit->getPosition()[0]  << std::endl;
    // std::cout << "PosY: "<< bcalhit->getPosition()[1]  << std::endl;
    // std::cout << "PosZ: "<< bcalhit->getPosition()[2]  << std::endl;

    dd4hep::Position global(double(bcalhit->getPosition()[0])*dd4hep::mm,
				      double(bcalhit->getPosition()[1])*dd4hep::mm,
				      double(bcalhit->getPosition()[2])*dd4hep::mm );
    //    dd4hep::Position local(0.0,0.0,0.0);


    // (global.z() > 0 ) ?
    //   m_BeamCal.child("BeamCal01").worldToLocal(global, local):
    //   m_BeamCal.child("BeamCal02").worldToLocal(global, local);

    // std::cout << "Local Cell ID: " << global << " -- > "  << local << std::endl;
    // std::cout << "  cid " << m_seg.cellID(local, global, 0)
    // 	      << std::endl;



    dd4hep::DDSegmentation::CellID cID = (((unsigned long long)(bcalhit->getCellID0()) ) +
					  ((unsigned long long)(bcalhit->getCellID1()) << 32));
    // cID = m_seg.cellID(local, global, 0);

    dd4hep::DDSegmentation::Vector3D pos = m_seg.position(cID);

    // std::cout << std::setw(25) << cID << ":" 
    // 	      << std::setw(15) << pos.x()/dd4hep::mm
    // 	      << std::setw(15) << pos.y()/dd4hep::mm
    // 	      << std::setw(15) << pos.z()/dd4hep::mm
    // 	      << std::endl;

    const float energy = bcalhit->getEnergy();

    m_x = pos.x();
    m_y = pos.y();
    m_z = bcalhit->getPosition()[2] >0 ? 1 : -1;
    m_energy = energy;
    m_tree->Fill();

    hitEnergies[bcalhit->getCellID1()] += energy;

    try{
      // (side == 0 ) ?
	// m_padEnergiesLeft->addEnergy(layer, cylinder, sector, energy):
	// m_padEnergiesRight->addEnergy(layer, cylinder, sector, energy);
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



void DrawBeamCalFromDD4hep::check( LCEvent * ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void DrawBeamCalFromDD4hep::end(){

  drawSegmentation();


  std::cout << "Tree entries: " << m_tree->GetEntries()  << std::endl;
  //  m_tree->Write();
  m_file->Write();
  m_file->Close();
  delete m_file;
  m_file = NULL;
  streamlog_out ( MESSAGE ) << __PRETTY_FUNCTION__ << " " << name()
			    << " processed " << m_nEvt << " events."
			    << std::endl ;
  //Do the average for every bin, and calculate the maximal difference to the mean, which means, we have to loop twice.


}


void DrawBeamCalFromDD4hep::drawSegmentation() {


#ifdef FCAL_USE_CustomRoot
  CustomRoot::SetCDRStyle2();
#endif


  if (m_seg.type() == "PolarGridRPhi2" ){
    drawPolarGridRPhi2();
  } else if (m_seg.type() == "CartesianGridXY") {
    drawCartesianGridXY();
  }

  return;
}
  
void DrawBeamCalFromDD4hep::drawCartesianGridXY() {

  TCanvas c1("xy", "xy", 800, 800);

  typedef dd4hep::DDSegmentation::TypedSegmentationParameter< double > ParDou;

  ParDou* parGridX = static_cast<ParDou*>(m_seg.segmentation()->parameter("grid_size_x"));
  ParDou* parGridY = static_cast<ParDou*>(m_seg.segmentation()->parameter("grid_size_y"));

  double gridX = parGridX->typedValue();
  double gridY = parGridY->typedValue();

  //aim for 20 by 20 cm size, but need multiple of bins
  const double size = 20;
  const int nBinsX = int(2*size/gridX+0.5);
  const double xRange = nBinsX*gridX*0.5+gridX*0.5;

  const int nBinsY = int(2*size/gridY+0.5);
  const double yRange = nBinsY*gridY*0.5+gridY*0.5;


  std::cout << "Bins In X: " << nBinsX << "   range in X: " << xRange  << std::endl;

  TH2D g("g", "g;x [cm]; y [cm]", nBinsX+1, -xRange, xRange, nBinsY+1, -yRange, yRange); 
  for (MapIdVal::iterator it = hitEnergies.begin(); it != hitEnergies.end(); ++it) {
    unsigned int cellid = it->first;
    double energy = it->second;
    auto const& decoder = *m_seg.segmentation()->decoder();
    auto llCellID = decoder.toLong(0, cellid);
    int xBin = decoder.get(llCellID, "x");
    int yBin = decoder.get(llCellID, "y");
    g.Fill( xBin*gridX, yBin*gridY, energy );
  }

  c1.SetLogz();
  g.SetAxisRange(1e-5, 1e2, "Z");
  g.Draw("col");

  c1.SaveAs("Segmentation.eps");

}


void DrawBeamCalFromDD4hep::drawPolarGridRPhi2() {

  TCanvas c1("rphi2", "rphi2", 800, 800);
  TH2D dummy("dummy", "dummy", 20, 0, 20, 200, 0, 200); //16 bits to get colours from pal->GetValueColor
  TH2D g("g", "g;x [cm]; y [cm]", 200, -20, 20, 200, -20, 20); //16 bits to get colours from pal->GetValueColor
  typedef dd4hep::DDSegmentation::TypedSegmentationParameter< std::vector<double> > ParVec;
  typedef dd4hep::DDSegmentation::TypedSegmentationParameter< double > ParDou;
  ParVec* rPar = static_cast<ParVec*>(m_seg.segmentation()->parameter("grid_r_values"));
  ParVec* pPar = static_cast<ParVec*>(m_seg.segmentation()->parameter("grid_phi_values"));

  ParDou* oPPar = static_cast<ParDou*>(m_seg.segmentation()->parameter("offset_phi"));


  std::vector<double> rValues = rPar->typedValue();
  std::vector<double> pValues = pPar->typedValue();
  double offsetPhi = oPPar->typedValue();

  for (MapIdVal::iterator it = hitEnergies.begin(); it != hitEnergies.end(); ++it) {
    unsigned int cellid = it->first;
    double energy = it->second;
    auto const& decoder = *m_seg.segmentation()->decoder();
    auto llCellID = decoder.toLong(0, cellid);
    int rBin = decoder.get(llCellID, "r");
    int pBin = decoder.get(llCellID, "phi");

    dummy.Fill( rBin, pBin, energy);
  }

  dummy.Draw("colz");
  gPad->SetLogz();
  dummy.SetAxisRange(1e-5, 1e2, "Z");
  gPad->Update();
  //  TPaletteAxis *pal = (TPaletteAxis*)dummy.GetListOfFunctions()->FindObject("palette");

  c1.SaveAs("dummy.eps");

  TCanvas c2("Segmentation", "Segmentation", 800, 800);
  c2.cd();
  g.Draw("AXIS");

  for (MapIdVal::iterator it = hitEnergies.begin(); it != hitEnergies.end(); ++it) {
    unsigned int cellid = it->first;
    double energy = it->second;
    auto const& decoder = *m_seg.segmentation()->decoder();
    auto llCellID = decoder.toLong(0, cellid);
    int rBin = decoder.get(llCellID, "r");
    int pBin = decoder.get(llCellID, "phi");

    const double offset = offsetPhi;
    const double rI = rValues[rBin];
    const double rO = rValues[rBin+1];
    const double pI = pValues[rBin]*pBin    +offset;
    const double pO = pValues[rBin]*(pBin+1)+offset;
  
    if (rBin == 3 && pBin > 47) continue;

    TCrown *tc = new TCrown(0,0, rI, rO, pI*TMath::RadToDeg(), pO*TMath::RadToDeg());
    std::cout 
      << std::setw(6) << rBin
      << std::setw(6) << pBin
      << std::setw(14) << rI
      << std::setw(14) << rO
      << std::setw(14) << pI*TMath::RadToDeg()
      << std::setw(14) << pO*TMath::RadToDeg()
      << std::setw(14) << energy
      << std::endl;
    tc->SetLineColor( getColor(energy) );
    tc->SetLineWidth(0);
    //    tc->SetFillColor( pal->GetValueColor(energy) );
    tc->SetFillColor( getColor(energy) );

    tc->Draw();
  }

  c2.SaveAs("Segmentation.eps");

}
