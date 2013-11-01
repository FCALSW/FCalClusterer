#include "BeamCalReco/BeamCalGeo.hh"
#include "BeamCalReco/BeamCalGeoGear.hh"
#include "BeamCalReco/BeamCalGeoCached.hh"
#include "BeamCalReco/BeamCal.hh"
#include "BeamCalReco/BCPadEnergies.hh"

//GEAR
#include <gearxml/GearXML.h>
#include <gear/GearMgr.h>

//ROOT
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TStyle.h>

#include <string>
#include <iostream>

void ReadRootFile(std::string const& fileName, std::vector<BCPadEnergies>& newPads);

void ReadRootFile(std::string const& fileName, std::vector<BCPadEnergies>& newPads) {

  TFile* file = TFile::Open(fileName.c_str());
  if ( not file ) {
    std::cout << "File not found: "<< fileName  << std::endl;
    return;
  }
  TTree* tree;
  file->GetObject("bcTree", tree);
  if ( not tree ) {
    std::cout << "Tree not found in file " << fileName  << std::endl;
    return;
  }

  std::vector<double> *depositsLeft(NULL), *depositsRight(NULL);
  tree->SetBranchAddress("vec_left" , &depositsLeft);
  tree->SetBranchAddress("vec_right", &depositsRight);
  tree->GetEntry(0);
  newPads[0].addEnergies(*depositsLeft);
  newPads[1].addEnergies(*depositsRight);

  file->Close();
  delete file;
  return;
}

int main (int argn, char **argc) {
  
  if ( argn < 3 ) {
    std::cout << "Not enough parameters"  << std::endl;
    std::cout << "DrawBeamCals <GearFile1> <BCPRootFile> "  << std::endl;
    std::exit(1);
  } 

  //Create the gearmanager from the gearfile
  std::string gearFile ( argc[1] );
  std::string rootFile ( argc[2]);
  gear::GearXML gearXML( gearFile ) ;
  gear::GearMgr* gearMgr = gearXML.createGearMgr() ;
  BeamCalGeo* geoCache = new BeamCalGeoCached (gearMgr);
  //  BeamCalGeo* geoGear = new BeamCalGeoGear (gearMgr);
  std::vector<BCPadEnergies> signalBeamCals(2, geoCache);
  ReadRootFile(rootFile, signalBeamCals);
  //Draw The BeamCal Energies::
  gStyle->SetOptStat(0);
  TCanvas c("c","c", 1000, 800);
  c.SetRightMargin(0.3);
  TH2F frame (" frame", "frame", 20, -200, 200, 20, -200, 200);
  BeamCal bc (gearMgr);
  //bc.SetLogz(true);
  bc.SetBeamCalHisto( &(signalBeamCals[0]) );
  bc.SetAxisMin(1e-6);
  bc.SetAxisMax(10.4);
  bc.SetLogz(true);
  bc.BeamCalDraw( &c, &frame );
  c.SaveAs("BeamCalCLIC.eps");

  return 0;
}
