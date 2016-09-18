#include "BCUtilities.hh"


void BCUtil::DecodeCellID(lcio::CellIDDecoder<lcio::SimCalorimeterHit> &mydecoder, const lcio::SimCalorimeterHit* hit, int& side, int& layer, int& cylinder, int& sector, bool usingDD4HEP){

  if( usingDD4HEP ){
    try{
      side     = mydecoder( hit )[ "barrel" ] - 1; // 1 and 2 originally, change to 0 and 1
      cylinder = mydecoder( hit )["r"] ;
      sector   = mydecoder( hit )["phi"] ;
      layer    = mydecoder( hit )["layer"] - 1 ; //starting at 1
    } catch (Exception &e) {
      std::cout << "Exception in BCUtil with DD4hep:" << e.what()  << std::endl;
    }
  } else {
    try{
      side     = mydecoder( hit )[ "S-1" ] ;
      cylinder = mydecoder( hit )["I"] ;
      sector   = mydecoder( hit )["J"] ;
      layer    = mydecoder( hit )["K"] ;
    } catch (Exception &e) {
      std::cout << "Exception in BCUtil without DD4hep:" << e.what()  << std::endl;
    }

  }

  return;
}//DecodeCellID
