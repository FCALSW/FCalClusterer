#ifndef OUTPUTMANAGERCLASS_H
#define OUTPUTMANAGERCLASS_H 1
// NOTE: Memory resident trees are not properly witten to disk
// this need to be fixed if such option is wanted
#define LCAL_MEMORY_RESIDENT_TREE 1

class TFile;
class TH1;
class TH2;
class TTree;

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

  int SkipNEvents, WriteRootTrees, NumEventsTree;
  int MemoryResidentTree;

  void  FillRootTree( const std::string & treeName );
  void	WriteToRootTree(std::string optName, int nEvtNow);
  void	Initialize(int treeLocNow, int skipNEventsNow , int numEventsTreeNow, std::string outDirNameNow, std::string outFileName);

  void	CleanUp();

private:
  // private method that may only be called ONCE in the destructor

  OutputManagerClass& operator=(OutputManagerClass const& rhs);
  OutputManagerClass(OutputManagerClass const&);
};

#endif
