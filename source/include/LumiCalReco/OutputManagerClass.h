
#include <TH1F.h>


#include <map>
#include <string>



class OutputManagerClass {

public:
  OutputManagerClass();
  ~OutputManagerClass();


  std::map < std::string , TH1 * >		HisMap1D;
  std::map < std::string , TH1 * > :: iterator	HisMap1DIterator;

  std::map < std::string , TH2 * >		HisMap2D;
  std::map < std::string , TH2 * > :: iterator	HisMap2DIterator;


  std::map < std::string , TTree * >		TreeMap;
  std::map < std::string , TTree * > :: iterator	TreeMapIterator;

  std::map < std::string , int > TreeIntV;
  std::map < std::string , double > TreeDoubleV;


  std::string	OutputRootFileName, OutDirName;
  TFile *OutputRootFile;

  std::map < std::string , int >		Counter;
  std::map < std::string , int >::iterator	CounterIterator;

  int	SkipNEvents, WriteRootTrees, NumEventsTree;

  void	WriteToRootTree(std::string optName, int nEvtNow);
  void	Initialize(int skipNEventsNow , int numEventsTreeNow, std::string outDirNameNow);

private:
  // private method that may only be called ONCE in the destructor
  void	CleanUp();

  OutputManagerClass& operator=(OutputManagerClass const& rhs);
  OutputManagerClass(OutputManagerClass const&);
};

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

  int numHis;

  // 1D histogram std::map
  numHis            = HisMap1D.size();
  HisMap1DIterator  = HisMap1D.begin();
  for(int hisNow = 0; hisNow < numHis; hisNow++ , HisMap1DIterator++) {
    std::string hisName = (std::string)(*HisMap1DIterator).first;

    delete HisMap1D[hisName];
  }

  // 2D histogram std::map
  numHis            = HisMap2D.size();
  HisMap2DIterator  = HisMap2D.begin();
  for(int hisNow = 0; hisNow < numHis; hisNow++ , HisMap2DIterator++) {
    std::string hisName = (std::string)(*HisMap2DIterator).first;

    delete HisMap2D[hisName];
  }

  // tree std::map
  numHis            = TreeMap.size();
  TreeMapIterator   = TreeMap.begin();
  for(int hisNow = 0; hisNow < numHis; hisNow++ , ++TreeMapIterator) {
    std::string hisName = (std::string)(*TreeMapIterator).first;

    delete TreeMap[hisName];
  }


  return;
}


void OutputManagerClass::WriteToRootTree(std::string optName, int nEvtNow){

  int numHis;

  if( int(nEvtNow/double(NumEventsTree)) ==  WriteRootTrees || optName == "forceWrite") {
    OutputRootFileName = "./";  OutputRootFileName += OutDirName;  OutputRootFileName += "/output_";
    OutputRootFileName += WriteRootTrees;  OutputRootFileName += ".root";
    OutputRootFile = new TFile(OutputRootFileName.c_str(),"RECREATE");

    // 1D histogram std::map
    numHis            = HisMap1D.size();
    HisMap1DIterator  = HisMap1D.begin();
    for(int hisNow = 0; hisNow < numHis; hisNow++ , HisMap1DIterator++) {
      std::string hisName = (std::string)(*HisMap1DIterator).first;
      HisMap1D[hisName] -> Write();
      HisMap1D[hisName] -> Reset();
    }

    // 2D histogram map
    numHis            = HisMap2D.size();
    HisMap2DIterator  = HisMap2D.begin();
    for(int hisNow = 0; hisNow < numHis; hisNow++ , HisMap2DIterator++) {
      std::string hisName = (std::string)(*HisMap2DIterator).first;
      HisMap2D[hisName] -> Write();
      HisMap2D[hisName] -> Reset();
    }

    // tree map
    numHis            = TreeMap.size();
    TreeMapIterator   = TreeMap.begin();
    for(int hisNow = 0; hisNow < numHis; hisNow++ , ++TreeMapIterator) {
      std::string hisName = (std::string)(*TreeMapIterator).first;
      TreeMap[hisName] -> Write();
      TreeMap[hisName] -> Reset();
    }

    OutputRootFile -> Close();
    delete	OutputRootFile;
    WriteRootTrees++;

    streamlog_out(DEBUG) << std::endl << "Went through  " << nEvtNow << "  events..." << std::endl;
  }

  return;
}
