#include "BeamCalGeo.hh"
#include "BeamCalGeoGear.hh"
#include "BeamCalGeoCached.hh"
#include "BeamCal.hh"
#include "BCPadEnergies.hh"
#include "BCRootUtilities.hh"

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


int main (int argn, char **argc) {

  if ( argn < 2 ) {
    std::cout << "Not enough parameters"  << std::endl;
    std::cout << "DrawBeamCals <GearFile1> [PRootFile] "  << std::endl;
    std::exit(1);
  } 
  const bool haveEnergy(argn > 2);

  //Create the gearmanager from the gearfile
  std::string gearFile ( argc[1] );
  gear::GearXML gearXML( gearFile ) ;
  gear::GearMgr* gearMgr = gearXML.createGearMgr() ;
  BeamCalGeo* geoCache = new BeamCalGeoCached (gearMgr);
  //  BeamCalGeo* geoGear = new BeamCalGeoGear (gearMgr);
  std::vector<BCPadEnergies> signalBeamCals(2, geoCache);
  if (haveEnergy) {
    std::string rootFile ( argc[2]);
    int wrong(0);
    try {
      BCUtil::ReadRootFile(rootFile, signalBeamCals);
    } catch (std::invalid_argument &e) {
      ++wrong;
    }
    try {
      BCUtil::ReadBecasFile(rootFile, signalBeamCals);
    } catch (std::invalid_argument &e) {
      ++wrong;
    }
    if (wrong==2) {
      std::cerr << "This file has neither BecAs nor BCPadEnergies as a tree" << std::endl;
      return 1;
    }
  }
  //Draw The BeamCal Energies::
  gStyle->SetOptStat(0);
  TCanvas c("c","c", 800, 800);
  //  c.SetRightMargin(0.3);
  TH2F frame ("frame", "frame", 20, -200, 200, 20, -200, 200);
  gStyle->SetOptTitle(0);
  BeamCal bc (gearMgr);
  //bc.SetLogz(true);
  bc.SetBeamCalHisto( &(signalBeamCals[0]) );
  bc.SetAxisMin(1e-6);
  bc.SetAxisMax(10.4);
  bc.SetLogz(true);
  bc.BeamCalDraw( &c, &frame );
  c.SaveAs("BeamCalLayout.eps");

  return 0;
}
