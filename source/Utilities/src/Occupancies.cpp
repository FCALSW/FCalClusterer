/**
 *  Occupancies calculates the occupancy in the LumiCal or BeamCal given a list
 *  of background files produced by the ReadBeamCal processor
 * Arguments: [BeamCal|LumiCal] <Threshold[GeV]> <CompactFile> <Background1.root> [<Background2.root ...>]
 * The larger the number of background files the better the averaging of course
 */


#include <BeamCalGeoDD.hh>
#include <BCPadEnergies.hh>
#include <BeamCal.hh>
#include <RootUtils.hh>

#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TPaletteAxis.h>

#include <DD4hep/Detector.h>

#include <string>
#include <iostream>
#include <iomanip>
#include <vector>

using VD=std::vector<double>;

void calculateOccupancy(std::string& detectorName, std::string& compactFile, double threshold,
                        std::vector<std::string> const& backgroundFiles);
void readBackgroundFile(std::string const& backgroundFile, VD& countLeft, VD& countRight, double threshold, int& nBX, double& totalEnergyBX);
void drawOccupancy(std::string const& detectorName, BeamCalGeo const& bcg,
                   std::string const& name, BCPadEnergies const& bcp);


int main (int argc, char **args) {

  if(argc < 4) {
    std::cout << "Not enough arguments "
              << "Occupancies [BeamCal|LumiCal] <Threshold[GeV]> compactFile Background1.root [Background2.root ...]"
              << std::endl;
    return 1;
  }

  std::string detectorName = std::string(args[1]);
  std::string compactFile =  std::string(args[3]);
  double threshold = std::atof(args[2]);

  std::vector<std::string> backgroundFiles{};
  for (int i = 4; i < argc ;++i) {
    backgroundFiles.emplace_back(args[i]);
  }


  RootUtils::SetStyle();

  try{
    calculateOccupancy(detectorName, compactFile, threshold, backgroundFiles);
  } catch(std::exception &e) {
    std::cout << "Exception " << e.what()  << std::endl;
    return 1;
  }

}


void calculateOccupancy(std::string& detectorName, std::string& compactFile, double threshold,
                        std::vector<std::string> const& backgroundFiles) {

  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();
  theDetector.fromCompact(compactFile);

  BeamCalGeoDD bcg(theDetector, detectorName, detectorName+"Collection");

  int nBX(0);
  int nPads = bcg.getPadsPerBeamCal();
  VD countLeft(nPads, 0), countRight(nPads, 0);
  double totalEnergy(0.0);

  std::cout << "Have " << backgroundFiles.size() << " Files"  << std::endl;
  for (auto const& backgroundFile: backgroundFiles ) {
    double totalEnergyBX = 0.0;
    readBackgroundFile(backgroundFile, countLeft, countRight, threshold, nBX, totalEnergyBX);
    totalEnergy += totalEnergyBX;
  }

  std::cout << "Have " << nBX << " bunch crossings"  << std::endl;
  std::cout << "Have " << countLeft.size() << " pads"  << std::endl;

  //Calculate average occupancy given threshold, loop over all pads and all
  //bunch crossings and count values above threshold why store all the values
  //and not directly count things above threshold?

  totalEnergy /= (double(nBX)*2.0); //counting both sides
  for (size_t nPad = 0; nPad < countLeft.size();++nPad) {
    countLeft[nPad] /= double(nBX);
    countRight[nPad] /= double(nBX);
  }

  BCPadEnergies left(bcg), right(bcg);

  left.setEnergies(countLeft);
  right.setEnergies(countRight);

  drawOccupancy(detectorName, bcg, "Left", left);
  drawOccupancy(detectorName, bcg, "Right", right);

  std::cout << "Total Energy deposit per BX in _one_ Detector: " << totalEnergy << " GeV" << std::endl;

}

void drawOccupancy(std::string const& detectorName, BeamCalGeo const& bcg,
                   std::string const& name, BCPadEnergies const& bcp) {
  TCanvas canv("cB", "cB", 800, 700);
  canv.SetRightMargin(0.21);
  canv.SetLeftMargin(0.15);
  TH2D occupancyHisto("hOcc","Occupancy;Layer;Radial Pad;Pad Occupancy [1/BX]", bcg.getBCLayers(), 0.5, bcg.getBCLayers()+0.5,
                      bcg.getBCRings(), -0.5, bcg.getBCRings());

  for (int layer = 0; layer < bcg.getBCLayers(); ++layer) {
    for (int ring = 0; ring < bcg.getBCRings() ; ++ring) {
      double averageOccupancy(0);
      for (int pad = 0; pad <  bcg.getPadsInRing(ring); ++pad) {
        averageOccupancy += bcp.getEnergy(layer, ring, pad);
      }
      averageOccupancy /= double(bcg.getPadsInRing(ring));
      occupancyHisto.Fill(layer+1, ring, averageOccupancy);
    }//for each ring
  }//for each layer
  occupancyHisto.Draw("colz");
  canv.SetLogz();
  canv.Update();
  RootUtils::MovePaletteHorizontally(&occupancyHisto, 0.0);
  occupancyHisto.GetZaxis()->SetLabelOffset(-0.01);
  gPad->Update();
  canv.Modified();
  canv.Update();
  canv.SaveAs(Form("%s_%sOccupancies.eps", detectorName.c_str(), name.c_str()));
  canv.SaveAs(Form("%s_%sOccupancies.C", detectorName.c_str(), name.c_str()));

}


void readBackgroundFile(std::string const& backgroundFile, VD& countLeft, VD& countRight, double threshold, int& nBX, double& totalEnergyBX){
  TTree* tree;
  VD *depLeft=NULL;
  VD *depRight=NULL;

  TFile* file = TFile::Open(backgroundFile.c_str());
  if (not file) {
    std::cerr << "File not found: "<< backgroundFile  << std::endl;
    throw std::runtime_error("Failed to find file");
  }

  file->GetObject("bcTree", tree);
  if (not tree) {
    std::cerr << "Tree not found in file " << backgroundFile  << std::endl;
    file->Close();
    delete file;
    throw std::runtime_error("Failed to find tree in file");
  }

  tree->SetBranchAddress("vec_left" , &depLeft);
  tree->SetBranchAddress("vec_right", &depRight);
  for (int i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    nBX++;
    for (size_t nPad = 0; nPad < countLeft.size();++nPad) {
      totalEnergyBX += (*depLeft)[nPad];
      totalEnergyBX += (*depRight)[nPad];
      if((*depLeft)[nPad] > threshold) {
        countLeft[nPad] += 1;
      }
      if((*depRight)[nPad] > threshold) {
        countRight[nPad] += 1;
      }
    }
  }
  file->Close();
}
