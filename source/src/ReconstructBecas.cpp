#include "BeamCal.hh"
#include "BCPadEnergies.hh"
#include "BeamCalGeoCached.hh"
#include "BeamCalCluster.hh"
#include "BCPCuts.hh"
#include "BCRootUtilities.hh"

//GEAR
#include <GEAR.h>
#include <gearxml/GearXML.h>

//ROOT
#include <TCanvas.h>
#include <TH2.h>
#include <TStyle.h>

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

class BeamCalGeo;

namespace gear {
  class GearMgr;
}

int reconstructBecas (int argn, char **argc) {
 
  if ( argn < 4 ) {
    throw std::invalid_argument("Not enough parameters\nReconstructBeCaS GearFile SignalFile backgroundSigmaFile");
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

  const std::vector<BeamCalCluster>& clusters0 = signalBeamCals[0].lookForNeighbouringClustersOverSigma(backgroundBeamCals[0], cuts, /*detaildPrintout=*/false);
  const std::vector<BeamCalCluster>& clusters1 = signalBeamCals[1].lookForNeighbouringClustersOverSigma(backgroundBeamCals[1], cuts, /*detaildPrintout=*/false);

  ///Output the reconstructed objects 
  std::cout << "cluster0"  << std::endl;
  for (std::vector<BeamCalCluster>::const_iterator it = clusters0.begin(); it != clusters0.end(); ++it) {
    std::cout << *it  << std::endl;
  }

  std::cout << "cluster1"  << std::endl;
  for (std::vector<BeamCalCluster>::const_iterator it = clusters1.begin(); it != clusters1.end(); ++it) {
    std::cout << *it  << std::endl;
  }


  //signalBeamCals[0].addEnergy(1234, 0.001);

  //Draw The BeamCal Energies::
  gStyle->SetOptStat(0);
  TCanvas c("c","c", 1000, 800);
  c.SetRightMargin(0.3);
  TH2F frame (" frame", "frame", 20, -200, 200, 20, -200, 200);
  BeamCal bc (*geo);
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

int main (int argn, char **argc) {

  try {
    return reconstructBecas(argn, argc);
  } catch (std::out_of_range &e) {
    std::cerr << "Geometry does not agree with energy in the trees:" << e.what()
	      << std::endl;
    return 1;
  } catch (std::invalid_argument &e) {
    std::cerr << e.what() << std::endl;
    return 1;
  } catch (gear::ParseException &e) {
    std::cerr << e.what();
    return 1;
  }

  return 0;
}
