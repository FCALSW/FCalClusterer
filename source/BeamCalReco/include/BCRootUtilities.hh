#ifndef BCRootUtilities_hh
#define BCRootUtilities_hh 1

#include "BCPadEnergies.hh"

class TFile;
class BeamCalGeo;

#include <string>
#include <vector>


namespace BCUtil{

  void ReadRootFile(std::string const& fileName, std::vector<BCPadEnergies>& newPads);
  void ReadBecasFile(std::string const& fileName, std::vector<BCPadEnergies>& newPads,
		     std::string treeName = "tSegment",
		     std::string energyField="sEdep",
		     bool isFromMokka = false
		     );

}//end namespace

#endif // BCRootUtilities_hh
