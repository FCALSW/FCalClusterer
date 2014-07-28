#include "BCUtilities.hh"


void BCUtil::DecodeCellID(lcio::CellIDDecoder<lcio::SimCalorimeterHit> &mydecoder, const lcio::SimCalorimeterHit* hit, int& side, int& layer, int& cylinder, int& sector){

  try{
    side     = mydecoder( hit )[ "S-1" ] ;
    cylinder = mydecoder( hit )["I"] ;
    sector   = mydecoder( hit )["J"] ;
    layer    = mydecoder( hit )["K"] ;
  } catch (Exception &e) {
    std::cout << e.what()  << std::endl;
  }

  return;
}//DecodeCellID
