#ifndef BCRootUtilities_hh
#define BCRootUtilities_hh 1

#include "BCPadEnergies.hh"

class TFile;
class BeamCalGeo;

#include <TError.h>

#include <string>
#include <vector>


namespace BCUtil{

  void ReadRootFile(std::string const& fileName, std::vector<BCPadEnergies>& newPads);
  void ReadBecasFile(std::string const& fileName, std::vector<BCPadEnergies>& newPads,
		     std::string treeName = "tSegment",
		     std::string energyField="sEdep",
		     bool isFromMokka = false
		     );


  class IgnoreRootError {
  public:
    IgnoreRootError(bool enable=true): orgErrorLevel(gErrorIgnoreLevel) {
      if( enable ) { gErrorIgnoreLevel=kError+1; }
    }
    ~IgnoreRootError(){
      gErrorIgnoreLevel=orgErrorLevel;
    }
  private:
    int orgErrorLevel;
  };

}//end namespace

#endif // BCRootUtilities_hh
