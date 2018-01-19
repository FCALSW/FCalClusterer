#ifndef BCutilities_hh
#define BCutilities_hh 1

#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/LayerLayout.h>
#include <gear/CalorimeterParameters.h>

#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <UTIL/CellIDDecoder.h>

namespace BCUtil{

  template<class t> inline t vectorprod(const t *a,const  t *b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  }

  template<class t> inline t vectorabs(const t *a) {
    return sqrt(vectorprod(a, a));
  }
  
  template<class t> inline t vectorangle(const t *a, const t *b) {
    return acos(vectorprod(a, b)/(vectorabs(a)*vectorabs(b)));
  }

  template<class t> inline void vectoradd(t *a, t *b) {
    a[0] += b[0];
    a[1] += b[1];
    a[2] += b[2];
    a[3] += b[3];
  }

  template<class t> inline void RotateToBeamCal(const t* vector, t* rotated, double angle) {
    const double CosAngle = cos(double(angle)/1000.0);
    const double SinAngle = sin(double(angle)/1000.0);
    rotated[0] = CosAngle*vector[0]-SinAngle*vector[2];
    rotated[1] = vector[1];
    rotated[2] = SinAngle*vector[0]+CosAngle*vector[2];
  }

  template<class t> inline t AngleToBeamCal(const t* vec, double angle) {
    const double CosAngle = cos(double(angle)/1000);
    const double SinAngle = sin(double(angle)/1000);
    if(vec[2] > 0) {
      const t zP[3] = {SinAngle, 0, CosAngle};
      return vectorangle(zP, vec);
    } else {
      //-Cos, because the beamcal is on the other side, but still rotated to plus x
      const t zP[3] = {SinAngle, 0, -CosAngle};
      return vectorangle(zP, vec);
    }  
  }

  template<class t> void RotateFromBeamCal(const t* vector, t* rotated, double angle) {
    const double CosAngle = cos(double(angle)/1000.0);
    const double SinAngle = sin(double(angle)/1000.0);
    rotated[0] = CosAngle*vector[0]+SinAngle*vector[2];
    rotated[1] = vector[1];
    rotated[2] = -SinAngle*vector[0]+CosAngle*vector[2];
  }


  template<class t> void RotateToLabFrame(const t* vector, t* rotated, double angle) {
    if(vector[2] > 0) {
      BCUtil::RotateFromBeamCal(vector, rotated, angle);
    } else {
      BCUtil::RotateToBeamCal(vector, rotated, angle);
    }
  }

  template<class t> void RotateToBCFrame(const t* vector, t* rotated, double angle) {
    if(vector[2] < 0) {
      BCUtil::RotateFromBeamCal(vector, rotated, angle);
    } else {
      BCUtil::RotateToBeamCal(vector, rotated, angle);
    }
  }

  template<class T, class U>
  void DecodeCellID(T &mydecoder, const U* hit, int& side, int& layer, int& cylinder, int& sector, bool usingDD4HEP=false){
    if (usingDD4HEP) {
      try{
        side     = mydecoder( hit )[ "barrel" ] - 1; // 1 and 2 originally, change to 0 and 1
        cylinder = mydecoder( hit )["r"] ;
        sector   = mydecoder( hit )["phi"] ;
        layer    = mydecoder( hit )["layer"]; // starting at 0
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


  // inline const gear::CalorimeterParameters& BCPs() { return marlin::Global::GEAR->getBeamCalParameters(); }

  // //Wrappters around Gear Interface:
  // inline double                getBCInnerRadius(){ return BCPs().getExtent()[0]; }
  // inline double                getBCOuterRadius(){ return BCPs().getExtent()[1]; }
  // inline int                   getBCLayers()     { return BCPs().getLayerLayout().getNLayers(); } 
  // inline int                   getBCRings()      { return BCPs().getDoubleVals("phi_segmentation").size(); }
  // inline std::vector<double>   getSegmentation() { return BCPs().getDoubleVals("phi_segmentation"); }
  // inline std::vector<int>      getNSegments()    { return BCPs().getIntVals("nPhi_segmentation"); }
  // inline double                getCutout()       { return BCPs().getDoubleVal("dead_area_outer_r")+0.1; }
  // inline double getBCZDistanceToIP() {    static const double globalBeamCalDist = 3350; return globalBeamCalDist; }


  inline double BXperYear () { return 134784000000.; } //100 * 24 * 60 * 60 * 50 * 312; // 


  //Functions from BCPadEnergies.hh
  //They should be made to only rely on the GearInterface wrappers from above...
  //Or we write a class which either uses Gear or DD4hep accesses and point to that...


  inline bool areCloseTogether(double theta1, double phi1, double theta2, double phi2) {
  const double DegToRad = M_PI/180.0;
  return (
	  ( fabs( theta1 - theta2 ) < 5.0 )
	  and (( fabs ( sin(phi1*DegToRad) - sin(phi2*DegToRad) ) < 0.35 ))
	  and (( fabs ( cos(phi1*DegToRad) - cos(phi2*DegToRad) ) < 0.35 )));
}

}//Namespace

#endif // BCutilities_hh
