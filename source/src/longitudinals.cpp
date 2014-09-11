/****************************************************************************
 *
 * longitudinals.cpp
 *
 * Create longitudinal projections of multiple BeamCal histograms
 * and store them in a file.
 *
 *  Strahinja Lukic, June 2014
 *
 ****************************************************************************/

#include "BeamCalGeo.hh"
#include "BeamCalGeoGear.hh"
#include "BeamCalGeoCached.hh"
#include "BeamCal.hh"
#include "BCPadEnergies.hh"
#include "BCRootUtilities.hh"
#include "ProfileExtractor.hh"

//GEAR
#include <gearxml/GearXML.h>
#include <gear/GearMgr.h>

//ROOT
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TStyle.h>

#include <string>
#include <iostream>
#include <sstream>



int main (int argn, char **argc) {

  if ( argn < 3 ) {
    std::cout << "Not enough parameters"  << std::endl;
    std::cout << "longitudinals <GearFile1> <PRootFile> [<profileType> <parameter> ] "  << std::endl;
    std::exit(1);
  } 

  ProfileExtractor* profiler;
  if (argn < 4 ) profiler = new ProfileAll;
  else if (TString(argc[3]).Contains("top", TString::kIgnoreCase))
  {
	  if(argn < 5)
	  {
		std::cout << "Not enough parameters for the ""top"" profile taker. "  << std::endl;
		std::cout << "longitudinals <GearFile1> <PRootFile> top <threshold> "  << std::endl;
		std::exit(1);
	  }
	  if(atof(argc[4]) <= 0. )
	  {
		std::cout << "Please set the threshold above zero." << std::endl;
		std::cout << "longitudinals <GearFile1> <PRootFile> top <threshold> "  << std::endl;
		std::exit(1);
	  }
	  profiler = new ProfileTopAndNNNTowers(atof(argc[4]));
  }
  else if (TString(argc[3]).Contains("cm", TString::kIgnoreCase))
  {
	  if(argn < 5)
	  {
		std::cout << "Not enough parameters for the ""cm"" profile taker. "  << std::endl;
		std::cout << "longitudinals <GearFile1> <PRootFile> cm <radius> "  << std::endl;
		std::exit(1);
	  }
	  profiler = new ProfileInRadiusFromCM(atof(argc[4]));
  }
  else profiler = new ProfileAll;

//  const double ares = 0.03; //Constant relating the width of the deposit distribution in a cell to the square root of the deposition
  const double e0 = 1.e-3; // Effective "unit energy" to express the cell deposit distribution as Poissonian
  const double threshold = 1.e-3;

  //Create the gearmanager from the gearfile
  std::string gearFile ( argc[1] );
  gear::GearXML gearXML( gearFile ) ;
  gear::GearMgr* gearMgr = gearXML.createGearMgr() ;
  BeamCalGeo* geoCache = new BeamCalGeoCached (gearMgr);
  //  BeamCalGeo* geoGear = new BeamCalGeoGear (gearMgr);
  std::vector<BCPadEnergies> *signalBeamCals = NULL;

  std::string rootFile ( argc[2]);
  int wrong(0);
  try {
    	signalBeamCals = BCUtil::ReadMultiRootFile(rootFile, geoCache);
    } catch (std::invalid_argument &e) {
      ++wrong;
  }
  if (!signalBeamCals) {
      std::cerr << "This file has no BCPadEnergies as a tree" << std::endl;
      return 1;
  }

  TString outName("longitudinals."); outName += rootFile;
  TFile histograms(outName.Data(), "RECREATE");

  int nEntries = signalBeamCals->size();

  TH1D *BeamCalProfile = profiler->extractProfile(&(signalBeamCals->at(0)));
  TH1D BCaverage(*BeamCalProfile); BCaverage.SetName("Layers_average");
  BCaverage.Reset();
  int nAdded=0;
  if(BeamCalProfile->Integral() > threshold)
	  { BeamCalProfile->Write("Layers_0"); BCaverage.Add(BeamCalProfile); nAdded++;}

// Loop to create the average profile

  for( int iEvt=1; iEvt<nEntries; iEvt++)
  {
	  if(BeamCalProfile) delete BeamCalProfile;
	  BeamCalProfile = profiler->extractProfile(&(signalBeamCals->at(iEvt)));
	  if(BeamCalProfile->Integral() > threshold)
	  {
		  BeamCalProfile->Write(Form("Layers_%i", iEvt));
		  BCaverage.Add(BeamCalProfile);
		  nAdded++;
	  }
  }

  BCaverage.Scale(1./double(nAdded));
  TH1D BCrms(BCaverage); BCrms.SetName("Layers_RMSdev");
  BCrms.Reset();
  TH1D BCtemp(BCrms); BCtemp.SetName("Layers_temp");
  TH1D BCsq(BCrms); BCsq.SetName("Layers_square");

// Loop to create the RMS of the individual bins
  for( int iEvt=0; iEvt<nEntries; iEvt++)
  {
	  if(BeamCalProfile) delete BeamCalProfile;
	  BeamCalProfile = profiler->extractProfile(&(signalBeamCals->at(iEvt)));
	  if(BeamCalProfile->Integral() > threshold)
	  {
		  BCtemp.Add(BeamCalProfile, &BCaverage, 1., -1.);
		  BCsq.Multiply(&BCtemp,&BCtemp);
		  BCrms.Add(&BCsq);
	  }
  }
  BCrms.Scale(1./double(nAdded));

// Set the errors of the average and the individual profiles.
  //Calculate the bin factor
  TH1D haRes(BCrms);
  for( int ibin=1; ibin<=BCaverage.GetNbinsX(); ibin++)
  {
	  if(BCaverage.GetBinContent(ibin) > 0)
		  haRes.SetBinContent(ibin, sqrt(BCrms.GetBinContent(ibin) / BCaverage.GetBinContent(ibin)));
	  BCaverage.SetBinError(ibin, sqrt(BCrms.GetBinContent(ibin) / double(nAdded)));
  }

  BCaverage.Write("Layers_average");
  haRes.Write("ares");

  for( int iEvt=0; iEvt<nEntries; iEvt++)
  {
    std::stringstream hname; hname << "Layers_" << iEvt;
    histograms.GetObject(hname.str().c_str(), BeamCalProfile);
    if(BeamCalProfile)
    {
      BeamCalProfile->Scale(1./e0);
      BeamCalProfile->Sumw2();
      BeamCalProfile->Scale(e0);
	  BeamCalProfile->Write(hname.str().c_str(), TObject::kOverwrite);
    }
  }

//  histograms.Write();


  return 0;
}
