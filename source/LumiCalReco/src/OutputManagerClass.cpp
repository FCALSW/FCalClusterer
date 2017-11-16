
#include "OutputManagerClass.h"
#include "Global.hh"
#include <TROOT.h>
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
#include <sstream>



OutputManagerClass::OutputManagerClass():
  HisMap1D(),
  HisMap1DIterator(),
  HisMap2D(),
  HisMap2DIterator(),
  TreeMap(),
  TreeMapIterator(),
  TreeIntV(),
  TreeDoubleV(),
  OutputRootFileName("LcalOut"), OutDirName("rootOut"),
  OutputRootFile(NULL),
  Counter(),
  CounterIterator(),
  SkipNEvents(0), WriteRootTrees(1), NumEventsTree(0),
  MemoryResidentTree(0)
{
}

OutputManagerClass::~OutputManagerClass() {

  CleanUp();
}


void OutputManagerClass::Initialize(int treeLocOptNow,
                                    int skipNEventsNow,
                                    int numEventsTreeNow,
                                    std::string outDirNameNow,
                                    std::string outFileNameNow){


  OutDirName = outDirNameNow;
  OutputRootFileName = outFileNameNow;

  SkipNEvents   = skipNEventsNow;
  NumEventsTree = numEventsTreeNow;

  //	Counter["High engy cluster In positive side"] = 0;
  //	Counter["High engy cluster In negative side"] = 0;
  //	Counter["Bhabhas after selection cuts"]   = 0;
  OutputRootFile = NULL;
  if ( outFileNameNow == "" ) return ; // nothing to do

  TH1F	* his1;
  TH2F	* his2;
  TTree	* tree;
  Int_t	markerCol = kRed-4, lineCol = kAzure+3, fillCol = kBlue-8;
  std::string hisName;  int numBins1;  double hisRange1[2];  int numBins2;  double hisRange2[2];
  MemoryResidentTree  = treeLocOptNow;
  WriteRootTrees = int(SkipNEvents/double(NumEventsTree))+1;


  /* check that the output directory exists
  std::string  outDirName = "cd "; outDirName += OutDirName;
  assert( !system(outDirName.c_str()) );
  */
  // it is probably more convenient to check and create if needed (BP) 
  std::string command = "mkdir -p ";
  command += OutDirName;
  if( system( command.c_str() ) ) {
    exit( EXIT_FAILURE );
  }
  if( MemoryResidentTree == 0 ) {
  // Open Root file before booking histos/Trees to have diskresident dir
    std::string outFileName = OutDirName;
    outFileName += "/";
    outFileName += OutputRootFileName;
    outFileName += ".root";
    OutputRootFile = new TFile( outFileName.c_str(),"RECREATE"); 
  }else{
    // make sure we are at top root dir (BP)
    gROOT->cd(); 
  }
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

  hisRange1[0] = -0.5;	hisRange1[1] = 19.5;	numBins1 = 20;
  hisRange2[0] = -0.5;	hisRange2[1] = 19.5;	numBins2 = 20;
  hisName = "NumClustersIn_numMCParticlesIn";
  his2 = new TH2F(hisName.c_str(),hisName.c_str(),numBins1,hisRange1[0],hisRange1[1],numBins2,hisRange2[0],hisRange2[1]);

  his2->SetMarkerColor(markerCol);
  his2->SetMarkerStyle(5);
  his2->SetMarkerSize(0.7);
  his2->GetXaxis()->SetTitle("Number of MC Particles");
  his2->GetXaxis()->SetTitleSize(0.05);
  his2->GetXaxis()->SetTitleOffset(1.);
  his2->GetYaxis()->SetTitle("Number of Reconstructed Particles");
  his2->GetYaxis()->SetTitleSize(0.05);
  his2->GetYaxis()->SetTitleOffset(1.);

  HisMap2D[hisName] = his2;


  hisName = "LumiRecoParticleTree";
  TreeIntV["nEvt"] = 0;
  TreeIntV["sign"] = 0;
  TreeIntV["pdg"] = 0;
  TreeIntV["outFlag"] = 0;
  TreeIntV["mFlag"] = 0;
  TreeIntV["nHits"] = 0;
  TreeIntV["highE"] = 0;

  TreeDoubleV["distTheta"]=0.;
  TreeDoubleV["distXY"]=0.;
  TreeDoubleV["engy"]=0.;
  TreeDoubleV["theta"]=0.;
  TreeDoubleV["phi"]=0.;
  TreeDoubleV["rzStart"] = 0.;
  TreeDoubleV["Xglob"] = 0.;
  TreeDoubleV["Yglob"] = 0.;
  TreeDoubleV["Zglob"] = 0.;

  tree = new TTree(hisName.c_str(),hisName.c_str());
  tree -> Branch("Event", &TreeIntV["nEvt"], "nEvt/I");
  tree -> Branch("OutFlag", &TreeIntV["outFlag"], "outFlag/I");
  tree -> Branch("MatchFlag", &TreeIntV["mFlag"], "mFlag/I");
  tree -> Branch("distTheta", &TreeDoubleV["distTheta"], "distTheta/D");
  tree -> Branch("distXY", &TreeDoubleV["distXY"], "distXY/D");
  tree -> Branch("HighEngy", &TreeIntV["highE"], "highE/I");
  tree -> Branch("Sign", &TreeIntV["sign"], "sign/I");
  tree -> Branch("NHits", &TreeIntV["nHits"], "nHits/I");
  tree -> Branch("Energy", &TreeDoubleV["engy"], "engy/D");
  tree -> Branch("EngyMC", &TreeDoubleV["engyMC"], "engyMC/D");
  tree -> Branch("Theta", &TreeDoubleV["theta"], "theta/D");
  tree -> Branch("ThetaMC", &TreeDoubleV["thetaMC"], "thetaMC/D");
  tree -> Branch("Phi", &TreeDoubleV["phi"], "phi/D");
  tree -> Branch("RZstart", &TreeDoubleV["rzStart"], "rzStart/D");
  tree -> Branch("Xstart", &TreeDoubleV["Xglob"], "Xglob/D");
  tree -> Branch("Ystart", &TreeDoubleV["Yglob"], "Yglob/D");
  tree -> Branch("Zstart", &TreeDoubleV["Zglob"], "Zglob/D");

  TreeMap[ hisName] = tree;
  TreeDoubleV["vtxX"] = 0.;
  TreeDoubleV["vtxY"] = 0.;
  TreeDoubleV["vtxZ"] = 0.;
  TreeDoubleV["begX"] = 0.;
  TreeDoubleV["begY"] = 0.;
  TreeDoubleV["begZ"] = 0.;
  TreeDoubleV["endX"] = 0.;
  TreeDoubleV["endY"] = 0.;
  TreeDoubleV["endZ"] = 0.;
 
  hisName = "LumiMCParticleTree";
  tree = new TTree(hisName.c_str(),hisName.c_str());
  tree -> Branch("Event", &TreeIntV["nEvt"], "nEvt/I");
  tree -> Branch("Sign", &TreeIntV["sign"], "sign/I");
  tree -> Branch("PDG", &TreeIntV["pdg"], "pdg/I");
  tree -> Branch("Energy", &TreeDoubleV["engy"], "engy/D");
  tree -> Branch("Theta", &TreeDoubleV["theta"], "theta/D");
  tree -> Branch("Phi", &TreeDoubleV["phi"], "phi/D");
  tree -> Branch("VX", &TreeDoubleV["vtxX"], "vtxX/D");
  tree -> Branch("VY", &TreeDoubleV["vtxY"], "vtxY/D");
  tree -> Branch("VZ", &TreeDoubleV["vtxZ"], "vtxZ/D");
  tree -> Branch("BegX", &TreeDoubleV["begX"], "begX/D");
  tree -> Branch("BegY", &TreeDoubleV["begY"], "begY/D");
  tree -> Branch("BegZ", &TreeDoubleV["begZ"], "begZ/D");
  tree -> Branch("EndX", &TreeDoubleV["endX"], "endX/D");
  tree -> Branch("EndY", &TreeDoubleV["endY"], "endY/D");
  tree -> Branch("EndZ", &TreeDoubleV["endZ"], "endZ/D");

  TreeMap[hisName] = tree;
  

  return;
}

void OutputManagerClass::CleanUp(){

  // 1D histogram std::map
  streamlog_out(DEBUG)<<"OutputManagerClass::CleanUp: cleaning 1D map .........\n";
  HisMap1DIterator  = HisMap1D.begin();
  for(; HisMap1DIterator != HisMap1D.end(); ++HisMap1DIterator) {
    if ( HisMap1DIterator->second ) delete HisMap1DIterator->second;
  }

  // 2D histogram std::map
  streamlog_out(DEBUG)<<"OutputManagerClass::CleanUp: cleaning 2D map .........\n";
  HisMap2DIterator  = HisMap2D.begin();
  for(;HisMap2DIterator != HisMap2D.end(); ++HisMap2DIterator) {
    delete HisMap2DIterator->second;
  }

  // tree std::map
  streamlog_out(DEBUG)<<"OutputManagerClass::CleanUp: cleaning Tree map .......\n";
  TreeMapIterator   = TreeMap.begin();
  for(; TreeMapIterator != TreeMap.end(); ++TreeMapIterator) {
    delete TreeMapIterator->second;
  }

  TreeMap.clear();
  HisMap2D.clear();
  HisMap1D.clear();

  return;
}

void OutputManagerClass::FillRootTree( const std::string &treeName){
  TreeMap[ treeName ]->Fill();
}

void OutputManagerClass::WriteToRootTree(std::string optName, int nEvtNow){

  if ( !OutputRootFile ) return;

    if ( MemoryResidentTree == 1 ){ 
    //memory resident tree 
    if( int(nEvtNow/double(NumEventsTree)) ==  WriteRootTrees || optName == "forceWrite") {
      // OutputRootFileName = "./";  one migth want to write in some other place (BP)
    std::string outFileName = OutDirName;
    outFileName += "/";
    outFileName += OutputRootFileName;
    outFileName += "_";
    //    OutputRootFileName += WriteRootTrees;  <---  bad idea (BP)
     std::stringstream fnum;
     fnum << WriteRootTrees-1;
     outFileName += fnum.str();
     outFileName += ".root";
     OutputRootFile = new TFile( outFileName.c_str(),"RECREATE");
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
      TTree* pTreeNow = TreeMapIterator->second;
      if( pTreeNow ) {
	pTreeNow -> Write();
        pTreeNow -> Reset();
      }
    }

    OutputRootFile -> Close();
    delete	OutputRootFile;
    streamlog_out(MESSAGE) <<"\n"<< "OutputManagerClass::WriteToRootTree: Went through  " << nEvtNow << "  events..." << "\n";
    streamlog_out(MESSAGE) <<"OutputManagerClass::WriteToRootTree: Root file  " << outFileName << " written" << "\n";
    WriteRootTrees++;
  }
  }else{
      // disk resident hits/trees
    if( optName == "forceWrite" ){
      OutputRootFile -> Write();
      OutputRootFile -> Close();
      streamlog_out(MESSAGE) <<"\n"<< "OutputManagerClass::WriteToRootTree: Went through  " << nEvtNow << "  events..." << "\n";
      streamlog_out(MESSAGE) <<"OutputManagerClass::WriteToRootTree: Root file  " << OutputRootFile->GetName() << " written" << "\n";
      delete	OutputRootFile;
    }
  }

  return;
}
