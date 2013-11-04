#ifndef BCRootUtilities_hh
#define BCRootUtilities_hh 1

#include "BCPadEnergies.hh"

class TFile;
class BeamCalGeo;

#include <string>
#include <vector>


namespace BCUtil{

  void ReadRootFile(std::string const& fileName, std::vector<BCPadEnergies>& newPads);
  void ReadBecasFile(std::string const& fileName, std::vector<BCPadEnergies>& newPads );

}//end namespace

#endif // BCRootUtilities_hh
