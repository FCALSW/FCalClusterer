#include "OutputManagerClass.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>

#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>

#include <algorithm>
#include <cassert>
#include <map>
#include <string>



OutputManagerClass::OutputManagerClass():
  HisMap1D(),
  HisMap1DIterator(),
  HisMap2D(),
  HisMap2DIterator(),
  TreeMap(),
  TreeMapIterator(),
  TreeIntV(),
  TreeDoubleV(),
  OutputRootFileName(""), OutDirName(""),
  OutputRootFile(NULL),
  Counter(),
  CounterIterator(),
  SkipNEvents(0), WriteRootTrees(0), NumEventsTree(0)
{
}

OutputManagerClass::~OutputManagerClass() {

  CleanUp();
}


void OutputManagerClass::Initialize(int skipNEventsNow, int numEventsTreeNow, std::string outDirNameNow){


  OutDirName = outDirNameNow;
  //	OutDirName = "rootOut";

  SkipNEvents   = skipNEventsNow;
  NumEventsTree = numEventsTreeNow;

  //	Counter["High engy cluster In positive side"] = 0;
  //	Counter["High engy cluster In negative side"] = 0;
  //	Counter["Bhabhas after selection cuts"]   = 0;



  TH1F	* his1;
  TH2F	* his2;
  TTree	* tree;
  Int_t	markerCol = kRed-4, lineCol = kAzure+3, fillCol = kBlue-8;
  std::string hisName;  int numBins1;  double hisRange1[2];  int numBins2;  double hisRange2[2];

  WriteRootTrees = int(SkipNEvents/double(NumEventsTree));


  // check that the output directory exists
  std::string  outDirName = "cd "; outDirName += OutDirName;
  assert( !system(outDirName.c_str()) );


  hisRange1[0] = 0;	hisRange1[1] = 150;	numBins1 = 3000;

  hisName = "totEnergyOut";
  his1 = new TH1F(hisName.c_str(),hisName.c_str(),numBins1,hisRange1[0],hisRange1[1]); HisMap1D[hisName] = his1;

  his1->SetLineColor(lineCol);
  his1->SetFillColor(fillCol);
  his1->GetXaxis()->SetTitle("E [GeV]");
  his1->GetXaxis()->SetTitleSize(0.05);
  his1->GetXaxis()->SetTitleOffset(1.);
  his1->GetYaxis()->SetTitle("Entries");
  his1->GetYaxis()->SetTitleSize(0.05);
  his1->GetYaxis()->SetTitleOffset(1.);

  HisMap1D[hisName] = his1;


  hisRange1[0] = 0;	hisRange1[1] = 300;	numBins1 = 300;

  hisName = "totEnergyIn";
  his1 = new TH1F(hisName.c_str(),hisName.c_str(),numBins1,hisRange1[0],hisRange1[1]);

  his1->SetLineColor(lineCol);
  his1->SetFillColor(fillCol);
  his1->GetXaxis()->SetTitle("E [GeV]");
  his1->GetXaxis()->SetTitleSize(0.05);
  his1->GetXaxis()->SetTitleOffset(1.);
  his1->GetYaxis()->SetTitle("Entries");
  his1->GetYaxis()->SetTitleSize(0.05);
  his1->GetYaxis()->SetTitleOffset(1.);

  HisMap1D[hisName] = his1;

  hisName = "higestEngyParticle_Engy";
  his1 = new TH1F(hisName.c_str(),hisName.c_str(),numBins1,hisRange1[0],hisRange1[1]); HisMap1D[hisName] = his1;

  his1->SetLineColor(lineCol);
  his1->SetFillColor(fillCol);
  his1->GetXaxis()->SetTitle("E [GeV]");
  his1->GetXaxis()->SetTitleSize(0.05);
  his1->GetXaxis()->SetTitleOffset(1.);
  his1->GetYaxis()->SetTitle("Entries");
  his1->GetYaxis()->SetTitleSize(0.05);
  his1->GetYaxis()->SetTitleOffset(1.);

  HisMap1D[hisName] = his1;

  hisRange1[0] = 0;	hisRange1[1] = .1;	numBins1 = 200;
  hisName = "higestEngyParticle_Theta";
  his1 = new TH1F(hisName.c_str(),hisName.c_str(),numBins1,hisRange1[0],hisRange1[1]); HisMap1D[hisName] = his1;

  his1->SetLineColor(lineCol);
  his1->SetFillColor(fillCol);
  his1->GetXaxis()->SetTitle("#theta [rad]");
  his1->GetXaxis()->SetTitleSize(0.05);
  his1->GetXaxis()->SetTitleOffset(1.);
  his1->GetYaxis()->SetTitle("Entries");
  his1->GetYaxis()->SetTitleSize(0.05);
  his1->GetYaxis()->SetTitleOffset(1.);

  HisMap1D[hisName] = his1;


  hisRange1[0] = 0;	hisRange1[1] = 300;	numBins1 = 300;
  hisRange2[0] = 0;	hisRange2[1] = .1;	numBins2 = 300;
  hisName = "thetaEnergyOut_DepositedEngy";
  his2 = new TH2F(hisName.c_str(),hisName.c_str(),numBins1,hisRange1[0],hisRange1[1],numBins2,hisRange2[0],hisRange2[1]);

  his2->SetMarkerColor(markerCol);
  his2->SetMarkerStyle(5);
  his2->SetMarkerSize(0.7);
  his2->GetXaxis()->SetTitle("E [GeV]");
  his2->GetXaxis()->SetTitleSize(0.05);
  his2->GetXaxis()->SetTitleOffset(1.);
  his2->GetYaxis()->SetTitle("#theta [rad]");
  his2->GetYaxis()->SetTitleSize(0.05);
  his2->GetYaxis()->SetTitleOffset(1.);

  HisMap2D[hisName] = his2;


  hisName = "bhabhaSelectionTree";
  tree = new TTree(hisName.c_str(),hisName.c_str());
  tree -> Branch("Event", &TreeIntV["nEvt"], "nEvt/I");
  tree -> Branch("Sign", &TreeIntV["sign"], "sign/I");
  tree -> Branch("Energy", &TreeDoubleV["engy"], "engy/D");
  tree -> Branch("Theta", &TreeDoubleV["theta"], "theta/D");
  tree -> Branch("Phi", &TreeDoubleV["phi"], "phi/D");

  TreeMap["bhabhaSelectionTree"] = tree;



  hisName = "mcParticleTree";
  tree = new TTree(hisName.c_str(),hisName.c_str());
  tree -> Branch("Event", &TreeIntV["nEvt"], "nEvt/I");
  tree -> Branch("Sign", &TreeIntV["sign"], "sign/I");
  tree -> Branch("PDG", &TreeIntV["pdg"], "pdg/I");
  tree -> Branch("Energy", &TreeDoubleV["engy"], "engy/D");
  tree -> Branch("Theta", &TreeDoubleV["theta"], "theta/D");
  tree -> Branch("VtxX", &TreeDoubleV["vtxX"], "vtxX/D");
  tree -> Branch("VtxY", &TreeDoubleV["vtxY"], "vtxY/D");
  tree -> Branch("VtxZ", &TreeDoubleV["vtxZ"], "vtxZ/D");
  tree -> Branch("EndX", &TreeDoubleV["endX"], "endX/D");
  tree -> Branch("EndY", &TreeDoubleV["endY"], "endY/D");
  tree -> Branch("EndZ", &TreeDoubleV["endZ"], "endZ/D");

  TreeMap["mcParticleTree"] = tree;


  return;
}

void OutputManagerClass::CleanUp(){

  // 1D histogram std::map
  HisMap1DIterator  = HisMap1D.begin();
  for(; HisMap1DIterator != HisMap1D.end(); ++HisMap1DIterator) {
    delete HisMap1DIterator->second;
  }

  // 2D histogram std::map
  HisMap2DIterator  = HisMap2D.begin();
  for(;HisMap2DIterator != HisMap2D.end(); ++HisMap2DIterator) {
    delete HisMap2DIterator->second;
  }

  // tree std::map
  TreeMapIterator   = TreeMap.begin();
  for(; TreeMapIterator != TreeMap.end(); ++TreeMapIterator) {
    delete TreeMapIterator->second;
  }

  TreeMap.clear();
  HisMap2D.clear();
  HisMap1D.clear();

  return;
}


void OutputManagerClass::WriteToRootTree(std::string optName, int nEvtNow){

  if( int(nEvtNow/double(NumEventsTree)) ==  WriteRootTrees || optName == "forceWrite") {
    OutputRootFileName = "./";  OutputRootFileName += OutDirName;  OutputRootFileName += "/output_";
    OutputRootFileName += WriteRootTrees;  OutputRootFileName += ".root";
    OutputRootFile = new TFile(OutputRootFileName.c_str(),"RECREATE");

    // 1D histogram std::map
    HisMap1DIterator  = HisMap1D.begin();
  for(; HisMap1DIterator != HisMap1D.end(); ++HisMap1DIterator) {
      HisMap1DIterator->second -> Write();
      HisMap1DIterator->second -> Reset();
    }

    // 2D histogram map
    HisMap2DIterator  = HisMap2D.begin();
    for(;HisMap2DIterator != HisMap2D.end();++HisMap2DIterator) {
      HisMap2DIterator->second -> Write();
      HisMap2DIterator->second -> Reset();
    }

    // tree map
    TreeMapIterator   = TreeMap.begin();
    for(; TreeMapIterator != TreeMap.end(); ++TreeMapIterator) {
      TreeMapIterator->second -> Write();
      TreeMapIterator->second -> Reset();
    }

    OutputRootFile -> Close();
    delete	OutputRootFile;
    WriteRootTrees++;

    streamlog_out(DEBUG) << std::endl << "Went through  " << nEvtNow << "  events..." << std::endl;
  }

  return;
}
