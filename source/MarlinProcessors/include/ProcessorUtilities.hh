#ifndef ProcessorUtilities_HH
#define ProcessorUtilities_HH 1

#include <BeamCalGeoCached.hh>
#include <BeamCalGeoGear.hh>

#include <marlin/Global.h>

//DD4hep (optional for now)
#ifdef FCAL_WITH_DD4HEP
#include <DD4hep/LCDD.h>
namespace DD4hep{
  namespace Geometry {
    struct DetElement;
    struct LCDD;
  }
}
#include <BeamCalGeoDD.hh>
#endif


namespace ProcessorUtilities {

  ///Creates the BeamCalGeometry either from DD4hep if compiled with DD4hep and
  ///the geometry is available or from GearFile in all other cases
  inline BeamCalGeo* getBeamCalGeo(bool& usingDD4HEP) {

#ifdef FCAL_WITH_DD4HEP
    DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();
    try {
      const DD4hep::Geometry::DetElement& beamcal = lcdd.detector("BeamCal");
      if (beamcal.isValid()){
	streamlog_out(DEBUG) << "Creating DD4hep Based geometry" << std::endl;
	usingDD4HEP = true;
	return new BeamCalGeoDD(lcdd);
      }
    } catch( std::runtime_error &e ) {
      streamlog_out(ERROR) << " Failed to created BeamCalGeometry from DD4hep: "
			   << e.what()
			   << std::endl;
      streamlog_out(ERROR) << " Falling back to using GEAR as geometry source."
			   << std::endl;
    } catch (...) {
      streamlog_out(ERROR) << " Falling back to using GEAR as geometry source."
			   << std::endl;
    }
#endif

    streamlog_out(DEBUG) << "Creating GEAR based geometry" << std::endl;
    return new BeamCalGeoCached(marlin::Global::GEAR);
  }

} //end namespace

#endif
