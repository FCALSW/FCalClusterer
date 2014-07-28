#include "BeamCalReco/BeamCal.hh"
#include "BeamCalReco/BCPadEnergies.hh"
#include "BeamCalReco/BeamCalGeoCached.hh"
#include "BeamCalReco/BeamCalGeoGear.hh"
#include "BeamCalReco/BeamCalCluster.hh"
#include "BeamCalReco/BCPCuts.hh"
#include "BeamCalReco/BCRootUtilities.hh"

//GEAR
#include <gearxml/GearXML.h>
#include <gear/GearMgr.h>

//ROOT
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TStyle.h>

#include <vector>
#include <iostream>
#include <iomanip>

int main (int argn, char **argc) {
 
  if ( argn < 4 ) {
    std::cout << "Not enough parameters"  << std::endl;
    std::cout << "ReconstructBeCaS GearFile SignalFile backgroundSigmaFile"  << std::endl;
    std::exit(1);
  } 

  std::string gearFile(argc[1]);
  std::string signalFile(argc[2]);
  std::string backgroundSigmaFile(argc[3]);

  BCPCuts cuts;
  cuts.setStartLayer(7).setSigmaCut(1.0).setMinimumTowerSize(5);

  ////////////////////////////////////////////////////////////////////////////////
  //Create the gearmanager from the gearfile
  gear::GearXML gearXML( gearFile ) ;
  gear::GearMgr* gearMgr = gearXML.createGearMgr() ;
  BeamCalGeo* geo = new BeamCalGeoCached (gearMgr);

  ////////////////////////////////////////////////////////////////////////////////
  // Read Becas files into BCPadEnergiesx
  std::cout << "Read Files"  << std::endl;
  std::vector<BCPadEnergies> signalBeamCals(2, geo), backgroundBeamCals(2, geo);

  BCUtil::ReadBecasFile(signalFile, signalBeamCals);

  BCUtil::ReadBecasFile(backgroundSigmaFile, backgroundBeamCals);
  //subtract sigma from signal
  signalBeamCals[0].subtractEnergies(backgroundBeamCals[0]);
  signalBeamCals[1].subtractEnergies(backgroundBeamCals[1]);

  ////////////////////////////////////////////////////////////////////////////////
  // Reconstruct based on the sigma criterium
  std::cout << "Reconstructing"  << std::endl;

  bool detailedPrintout;
  std::vector<BeamCalCluster> clusters0 = signalBeamCals[0].lookForNeighbouringClustersOverSigma(backgroundBeamCals[0], cuts, detailedPrintout=false);
  std::vector<BeamCalCluster> clusters1 = signalBeamCals[1].lookForNeighbouringClustersOverSigma(backgroundBeamCals[1], cuts, detailedPrintout=false);

  ///Output the reconstructed objects 
  for (std::vector<BeamCalCluster>::iterator it = clusters0.begin(); it != clusters0.end(); ++it) {
    std::cout << *it  << std::endl;
  }

  for (std::vector<BeamCalCluster>::iterator it = clusters1.begin(); it != clusters1.end(); ++it) {
    std::cout << *it  << std::endl;
  }


  //signalBeamCals[0].addEnergy(1234, 0.001);

  //Draw The BeamCal Energies::
  gStyle->SetOptStat(0);
  TCanvas c("c","c", 1000, 800);
  c.SetRightMargin(0.3);
  TH2F frame (" frame", "frame", 20, -200, 200, 20, -200, 200);
  BeamCal bc (gearMgr);
  //bc.SetLogz(true);

  bc.SetBeamCalHisto( &(signalBeamCals[1]) );
  bc.SetAxisMin(1e-6);
  bc.SetAxisMax(0.4);
  bc.SetLogz(true);
  bc.BeamCalDraw( &c, &frame );
  c.SaveAs("BeamCal.eps");

  // for (int layer = 1; layer <= geo->getBCLayers(); ++layer) {
  //   bc.BeamCalDraw( &c, &frame, layer);
  //   c.SaveAs(Form("BeamCal_L%i.png", layer));
  // }

  return 0;
}

