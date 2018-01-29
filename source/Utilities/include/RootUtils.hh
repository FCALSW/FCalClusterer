#ifndef RootUtils_hh
#define RootUtils_hh 1

#include <RtypesCore.h>

#include <vector>

class TEfficiency;
class TGraph;
class TGraphErrors;
class TH1;


namespace RootUtils {

  class Colors {
  private:
    std::vector<Color_t> m_colors;
  public:
    Colors();
    Color_t GetColor();
    void SetColors();
  };//Class Colors

  void SetAllColors( TH1* object, Color_t color );
  void SetAllColors( TEfficiency* object, Color_t color );
  void SetAllColors( TGraph* object, Color_t color );
  void SetAllColors( TGraphErrors* object, Color_t color );


}//namespace

#endif // RootUtils_hh
