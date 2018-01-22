#include "BeamCalBkg.hh"
#include "BeamCalBkgAverage.hh"
#include "BeamCalBkgEmpty.hh"
#include "BeamCalBkgGauss.hh"
#include "BeamCalBkgParam.hh"
#include "BeamCalBkgPregen.hh"
#include "BeamCalGeo.hh"

#include <streamlog/streamlog.h>

#include <iostream>
#include <string>

BeamCalBkg* BeamCalBkg::Factory(std::string const& backgroundMethod, BeamCalGeo const* BCG) {
  if (std::string("Pregenerated") == backgroundMethod) {
    return new BeamCalBkgPregen(backgroundMethod, BCG);
  } else if (std::string("Gaussian") == backgroundMethod) {
    return new BeamCalBkgGauss(backgroundMethod, BCG);
  } else if (std::string("Parametrised") == backgroundMethod) {
    return new BeamCalBkgParam(backgroundMethod, BCG);
  } else if (std::string("Averaged") == backgroundMethod) {
    return new BeamCalBkgAverage(backgroundMethod, BCG);
  } else if (std::string("Empty") == backgroundMethod) {
    return new BeamCalBkgEmpty(backgroundMethod, BCG);
  } else {
    streamlog_out(ERROR) << "== Error From BeamCalBackground == "
                         << "Unknown BeamCal method of background esitmation \"" << backgroundMethod << "\"" << std::endl;
    throw std::runtime_error("Unknown BeamCal background method");
  }
}
