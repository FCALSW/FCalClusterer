#ifndef OUTPUTMANAGERCLASS_H
#define OUTPUTMANAGERCLASS_H 1

class TFile;
class TH1;
class TH2;
class TTree;

#include <map>
#include <string>
#include <vector>

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

  void	CleanUp();

private:
  // private method that may only be called ONCE in the destructor

  OutputManagerClass& operator=(OutputManagerClass const& rhs);
  OutputManagerClass(OutputManagerClass const&);
};

#endif
