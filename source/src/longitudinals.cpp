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

class profileTaker
{
public:
	profileTaker();
	~profileAll();
	TH1D* takeProfile(BCPadEnergies *bcpads) const = 0;
};

class profileAll : profileTaker
{
public:
	profileAll();
	~profileAll();

	TH1D* takeProfile(BCPadEnergies *bcpads) const {return bcpads->longitudinalProfile();}
};

class profileTopAndNNNTowers : profileTaker
{
public:
	profileTopAndNNNTowers(double threshold = 0.001) {m_threshold = threshold;}
	~profileTopAndNNNTowers();
	TH1D* takeProfile(BCPadEnergies *bcpads) const
	{
	  BCPadEnergies::TowerIndexList *towerlist = bcpads->getTopAndNNNeighbourTowers(m_threshold);
	  TH1D *profile = bcpads->longitudinalProfile(towerlist);
	  delete towerlist;
	  return profile;
	}
private:
	double m_threshold;
};

class profileInRadiusFromCM : profileTaker
{
public:
	profileInRadiusFromCM(double radius = 1.) {m_radius = radius;};
	~profileInRadiusFromCM();
	TH1D* takeProfile(BCPadEnergies *bcpads) const
	{
	  double z, rho, phi;
	  bcpads->getGlobalCM(z, rho, phi);
	  BCPadEnergies::TowerIndexList *towerlist = bcpads->getTowersWithinRadiusFromPoint(rho, phi, m_radius);
	  TH1D *profile = bcpads->longitudinalProfile(towerlist);
	  delete towerlist;
	  return profile;
	}
private:
	double m_radius;
};



int main (int argn, char **argc) {

  if ( argn < 3 ) {
    std::cout << "Not enough parameters"  << std::endl;
    std::cout << "longitudinals <GearFile1> <PRootFile> [profileType] "  << std::endl;
    std::exit(1);
  } 

  profileTaker* profiler;
  if (argn < 4 ) profiler = new profileAll;
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
	  profiler = new profileTopAndNNNTowers(atof(argc[4]));
  }
  else if (TString(argc[3]).Contains("cm", TString::kIgnoreCase))
  {
	  if(argn < 5)
	  {
		std::cout << "Not enough parameters for the ""cm"" profile taker. "  << std::endl;
		std::cout << "longitudinals <GearFile1> <PRootFile> cm <radius> "  << std::endl;
		std::exit(1);
	  }
	  profiler = new profileInRadiusFromCM(atof(argc[4]));
  }
  else profiler = new profileAll;

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

  TH1D *BeamCalProfile = takeProfile(&(signalBeamCals->at(0)));
  TH1D BCaverage(*BeamCalProfile);
  BCaverage.Reset();
  int nAdded=0;
  if(BeamCalProfile->Integral() > threshold)
	  { BeamCalProfile->Write("Layers_0"); BCaverage.Add(BeamCalProfile); nAdded++;}

// Loop to create the average profile

  for( int iEvt=1; iEvt<nEntries; iEvt++)
  {
	  delete BeamCalProfile; BeamCalProfile = takeProfile(&(signalBeamCals->at(iEvt)));
	  if(BeamCalProfile->Integral() > threshold)
	  {
		  BeamCalProfile->Write(Form("Layers_%i", iEvt));
		  BCaverage.Add(BeamCalProfile);
		  nAdded++;
	  }
  }

  BCaverage.Scale(1./double(nAdded));
  TH1D BCrms(BCaverage);
  BCrms.Reset();
  TH1D BCtemp(BCrms);
  TH1D BCsq(BCrms);

// Loop to create the RMS of the individual bins
  for( int iEvt=0; iEvt<nEntries; iEvt++)
  {
	  delete BeamCalProfile; BeamCalProfile = takeProfile(&(signalBeamCals->at(iEvt)));
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
