#include "BeamCalReco/BeamCal.hh"
#include "BeamCalReco/BCPadEnergies.hh"
#include "BeamCalReco/BeamCalGeoCached.hh"
#include "BeamCalReco/BeamCalGeoGear.hh"
#include "BeamCalReco/BeamCalCluster.hh"
#include "BeamCalReco/BCPCuts.hh"

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

void ReadBecasFile(std::string filename, BeamCalGeo* geo, std::vector<BCPadEnergies>& newPads );
void ReadBecasFromTree(TFile* becasFile, BeamCalGeo* geo, std::vector<BCPadEnergies>& newPads );

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

  ReadBecasFile(signalFile, geo, signalBeamCals);

  ReadBecasFile(backgroundSigmaFile, geo, backgroundBeamCals);
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



void ReadBecasFile(std::string filename, BeamCalGeo* geo, std::vector<BCPadEnergies>& newPads) {

  TFile *becasFile = TFile::Open(TString(filename));
  if( not becasFile ) {
    std::cerr << "File not found! " << filename  << std::endl;
    exit(1);
  }

  //check if this histogram exists, if it does we read the energy from histograms
  TH2F *testhisto(NULL);
  becasFile->GetObject("EnergyDepositionRPhiLr0", testhisto);
  ReadBecasFromTree( becasFile, geo, newPads );

  becasFile->Close();

  return;

}


void ReadBecasFromTree(TFile* becasFile, BeamCalGeo* geo, std::vector<BCPadEnergies>& newPads ) {

  TTree *becasTree(NULL);
  becasFile->GetObject("tSegment",becasTree);

  int position[3];
  double energy;
  becasTree->SetBranchAddress("sPos",position);
  becasTree->SetBranchAddress("sEdep",&energy);

  for (int i = 0 ; i < becasTree->GetEntries(); ++i)  {
    becasTree->GetEntry(i);
    int layer(position[0]), ring(position[1]), phiPad(position[2]);
    const int globalPadIndex = geo->getPadIndex(abs(layer), ring, phiPad);
    BCPadEnergies::BeamCalSide_t side = ( layer < 0 ) ? BCPadEnergies::kLeft : BCPadEnergies::kRight ;
    newPads[side].addEnergy( globalPadIndex, energy);

  }//for all entries

  

  return;

}// ReadBecasFromTree
