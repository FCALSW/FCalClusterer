#include "BeamCalClusterReco.hh"

#include "BCPCuts.hh"
#include "BCPadEnergies.hh"
#include "BCRecoObject.hh"
#include "BeamCal.hh"
#include "BeamCalCluster.hh"
#include "BCUtilities.hh"
#include "BeamCalGeoCached.hh"

//LCIO
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/SimCalorimeterHit.h>
#include <UTIL/CellIDDecoder.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/ClusterImpl.h>

// ----- include for verbosity dependent logging ---------
#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>

#include <marlin/ProcessorEventSeeder.h>
#include <marlin/Global.h>

//ROOT
#include <TCanvas.h>
#include <TChain.h>
#include <TCrown.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TRandom3.h>

//STDLIB
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>

using namespace lcio ;
using namespace marlin ;

BeamCalClusterReco aBeamCalClusterReco ;

//Just for formatting /////////////////////////////////////////////////////////////////
#define LONGSTRING "                                                                  "
///////////////////////////////////////////////////////////////////////////////////////

class LackingFilesException :public std::runtime_error {
public:
LackingFilesException(std::string error) : std::runtime_error(error) { }
};

class WrongParameterException :public std::runtime_error {
public:
WrongParameterException(std::string error) : std::runtime_error(error) { }
};


BeamCalClusterReco::BeamCalClusterReco() : Processor("BeamCalClusterReco"),
                                           m_colNameMC(""),
                                           m_colNameBCal(""),
                                           m_files(),
                                           m_nEvt(0),
                                           m_specialEvent(-1),
                                           m_nBXtoOverlay(0),
                                           m_eventSide(-1),
                                           m_minimumTowerSize(0),
                                           m_startLookingInLayer(0),
                                           m_usePadCuts(true),
                                           m_createEfficienyFile(false),
                                           m_sigmaCut(1.0),
                                           m_startingRings(),
                                           m_requiredRemainingEnergy(),
                                           m_requiredClusterEnergy(),
                                           m_BeamCalDepositsLeft(NULL),
                                           m_BeamCalDepositsRight(NULL),
                                           m_BeamCalAverageLeft(NULL),
                                           m_BeamCalAverageRight(NULL),
                                           m_BeamCalErrorsLeft(NULL),
                                           m_BeamCalErrorsRight(NULL),
                                           m_random3(NULL),
                                           m_backgroundBX(0),
                                           m_BCG(NULL),
                                           m_bcpCuts(NULL),
                                           m_thetaEfficieny(NULL),
                                           m_phiEfficiency(NULL),
                                           m_twoDEfficiency(NULL),
                                           m_phiFake(NULL),
                                           m_thetaFake(NULL),
                                           m_originalParticles(0),
                                           m_BCalClusterColName(""),
                                           m_BCalRPColName(""),
                                           m_EfficiencyFileName("")
{

// modify processor description
_description = "BeamCalClusterReco takes a list of beamcal background files from the ReadBeamCal"\
  "processor, adds NumberOfBX of them together and puts the signal hits from the"\
  "lcio input file on top of that, and then clustering is attempted." ;


// register steering parameters: name, description, class-variable, default value

registerInputCollection( LCIO::MCPARTICLE,
			   "MCParticle Collection Name, only needed and used to estimate efficiencies" ,
			   "Name of the MCParticle Collection"  ,
			   m_colNameMC ,
			   std::string("MCParticle") ) ;


registerInputCollection( LCIO::SIMCALORIMETERHIT,
			   "BeamCalCollectionName" ,
			   "Name of BeamCal Collection"  ,
			   m_colNameBCal ,
			   std::string("BeamCalCollection") ) ;

std::vector<std::string> defaultFile;
defaultFile.push_back("BeamCal.root");

registerProcessorParameter ("InputFileBackgrounds",
			      "Root Inputfile(s)",
			      m_files,
			      defaultFile ) ;

registerProcessorParameter ("NumberOfBX",
			      "Number of Bunch Crossings of Background",
			      m_nBXtoOverlay,
			      int(1) ) ;


std::vector<float> startingRing, padCut, clusterCut;
startingRing.push_back(0.0);  padCut.push_back(0.5);  clusterCut.push_back(3.0);
startingRing.push_back(1.0);  padCut.push_back(0.3);  clusterCut.push_back(2.0);
startingRing.push_back(2.0);  padCut.push_back(0.2);  clusterCut.push_back(1.0);

registerProcessorParameter ("StartingRing",
			      "Rings from which onwards the outside Thresholds are used",
			      m_startingRings,
			      startingRing ) ;

registerProcessorParameter ("ETPad",
			      "Energy in a Pad, after subtraction of background required to consider it for signal",
			      m_requiredRemainingEnergy,
			      padCut ) ;

registerProcessorParameter ("ETCluster",
			      "Energy in a Cluster to consider it an electron",
			      m_requiredClusterEnergy,
			      clusterCut) ;

registerProcessorParameter ("MinimumTowerSize",
			      "Minimum number of pads in a single tower to be considered for signal",
			      m_minimumTowerSize,
			      int(4) ) ;

registerProcessorParameter ("StartLookingInLayer",
			      "Layer (inclusive) from which on we start looking for signal clusters",
			      m_startLookingInLayer,
			      int(10) ) ;


registerProcessorParameter ("UseConstPadCuts",
			      "Use the cuts for the pads specified in ETPad, if false, the variance in each pad is used times SigmaPad Factor",
			      m_usePadCuts,
			      true ) ;

registerProcessorParameter ("SigmaCut",
			      "If not using ConstPadCuts, each pad SigmaCut*variance is considered for clusters",
			      m_sigmaCut,
			      double(1.0) ) ;


// registerProcessorParameter ("PDFFile",
//			      "Title of Output PDF only if verbosity is DEBUG!",
//			      m_pdftitle,
//			      std::string("BeamCalClusterReco.pdf") ) ;

registerProcessorParameter ("PrintThisEvent",
			    "Number of Event that should be printed to PDF File",
			    m_specialEvent,
			    int(-1) ) ;

// registerProcessorParameter ("LimitAreaOfBeamCal",
//			      "Skips events with electrons at the edge of the BeamCal or near the Keyhole cutout.",
//			      m_LimitedAreaOfBeamCal,
//			      bool(true) ) ;

registerProcessorParameter ("CreateEfficiencyFile",
			    "Flag to create the TEfficiency for fast tagging library",
			    m_createEfficienyFile,
			    false ) ;

registerProcessorParameter ("EfficiencyFilename",
			    "The name of the rootFile which will contain the TEfficiency objects",
			    m_EfficiencyFileName,
			    std::string("TaggingEfficiency.root") ) ;


 registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "RecoParticleCollectionname" ,
			   "Name of the Reconstructed Particle collection"  ,
			   m_BCalRPColName,
			   std::string("BCalRecoParticle") ) ;

 registerOutputCollection( LCIO::CLUSTER,
			   "RecoClusterCollectionname" ,
			   "Name of the Reconstructed Cluster collection"  ,
			   m_BCalClusterColName ,
			   std::string("BCalClusters") ) ;





}

void BeamCalClusterReco::init() {

  Global::EVENTSEEDER->registerProcessor(this);
  m_random3 = new TRandom3();

  // usually a good idea to
  if( streamlog::out.write< MESSAGE3 >() ) {
    printParameters() ;
  }

  m_nEvt = 0;

#pragma message "FIXME: Allow for possibility to use precomputed sigma file"
#pragma message "FIXME: Need energy calibration"

  if ( (m_startingRings.size() != m_requiredClusterEnergy.size() ) ||
       (m_requiredClusterEnergy.size() != m_requiredRemainingEnergy.size())){
    throw WrongParameterException("== Error From BeamCalClusterReco == The number of starting rings and required" \
				  "cluster energy or pad energy are not the same!");
  }
  if( m_startingRings[0] != 0 ) {
    throw WrongParameterException("== Error From BeamCalClusterReco == startingRings must always start with 0");
  }


  //Open the Files given as the list into a TChain...
  m_backgroundBX = new TChain("bcTree");

  //mix up the files, because the random numbers are ordered to avoid repeating
  std::random_shuffle(m_files.begin(), m_files.end());

  for (std::vector<std::string>::iterator file = m_files.begin(); file != m_files.end(); ++file) {
    streamlog_out(DEBUG1) << *file << std::endl;
    m_backgroundBX->Add(TString(*file));
  }

  //Ready the energy deposit vectors for the tree
  m_BeamCalDepositsLeft=NULL;
  m_BeamCalDepositsRight=NULL;

  m_backgroundBX->SetBranchAddress("vec_left" , &m_BeamCalDepositsLeft);
  m_backgroundBX->SetBranchAddress("vec_right", &m_BeamCalDepositsRight);

  streamlog_out(DEBUG2) << "We have " << m_backgroundBX->GetEntries() << " background BXs" << std::endl;

  //Create an Average BeamCal, needed to setup the BCPadEnergies
  m_BCG = new BeamCalGeoCached(marlin::Global::GEAR);

  m_BeamCalAverageLeft  =  new BCPadEnergies(m_BCG);
  m_BeamCalAverageRight =  new BCPadEnergies(m_BCG);


  //Fill BCPCuts object with cuts from the processor parameters
  m_bcpCuts = new BCPCuts(m_startingRings,
			  m_requiredRemainingEnergy, m_requiredClusterEnergy,
			  m_minimumTowerSize,
			  m_startLookingInLayer,
			  m_usePadCuts,
			  m_sigmaCut);


  // Create the average BCPadEnergies objects used to subtract average backgrounds
  std::vector<BCPadEnergies*> listOfBunchCrossingsRight, listOfBunchCrossingsLeft;

  for (int i = 0; i < 10; ++i) {
    listOfBunchCrossingsLeft.push_back( new BCPadEnergies(m_BCG) );
    listOfBunchCrossingsRight.push_back( new BCPadEnergies(m_BCG) );
  }

  std::set<int> randomNumbers;
  const unsigned int nBackgroundBX = m_backgroundBX->GetEntries();

  //Check that we have
  if( int(nBackgroundBX) < m_nBXtoOverlay*10 ) {
    throw LackingFilesException("There are not Enough BeamCal Background files to calculate a proper average!");
  }

  //we use a set so no duplication occurs
  const int numberForAverage = 10;
  while( int(randomNumbers.size()) < m_nBXtoOverlay*numberForAverage ){//do it ten times as often
    randomNumbers.insert( int(m_random3->Uniform(0, nBackgroundBX)) );
  }

  int counter = 0;


  for (std::set<int>::iterator it = randomNumbers.begin(); it != randomNumbers.end();++it) {
    streamlog_out(DEBUG1) << std::setw(5) << *it << std::flush;
    m_backgroundBX->GetEntry(*it);
    m_BeamCalAverageLeft-> addEnergies(*m_BeamCalDepositsLeft);
    m_BeamCalAverageRight->addEnergies(*m_BeamCalDepositsRight);
    listOfBunchCrossingsLeft.at(counter/m_nBXtoOverlay)-> addEnergies(*m_BeamCalDepositsLeft);
    listOfBunchCrossingsRight.at(counter/m_nBXtoOverlay)-> addEnergies(*m_BeamCalDepositsRight);
    ++counter;
  }

  //Now divide by ten, to get the average distributions...
  m_BeamCalAverageLeft ->scaleEnergies(1./double(numberForAverage));
  m_BeamCalAverageRight->scaleEnergies(1./double(numberForAverage));
  Double_t totalEnergyMean = m_BeamCalAverageLeft->getTotalEnergy();
  Double_t varEn(0.0);
  for (int l = 0; l < numberForAverage;++l) {
    varEn += (listOfBunchCrossingsRight[l]->getTotalEnergy() - totalEnergyMean) * (listOfBunchCrossingsRight[l]->getTotalEnergy() - totalEnergyMean);
  }//histograms
  varEn /= double(numberForAverage);
  varEn = sqrt(varEn);
  streamlog_out(MESSAGE4) << "Total Energy " << totalEnergyMean << " +- "
			  << varEn << " GeV/" << m_nBXtoOverlay <<"BX" << std::endl;

  // m_h3BeamCalAverageLeft = bc.getBeamCalHistogram("AverageLeft");
  // m_h3BeamCalAverageRight = bc.getBeamCalHistogram("AverageRight");

  //And now calculate the error for every bin....
  m_BeamCalErrorsLeft  = getBeamCalErrors(m_BeamCalAverageLeft,  listOfBunchCrossingsLeft, numberForAverage);
  m_BeamCalErrorsRight = getBeamCalErrors(m_BeamCalAverageRight, listOfBunchCrossingsRight, numberForAverage);

  //Add one sigma to the averages -- > just do it once here
  m_BeamCalAverageLeft ->addEnergies( m_BeamCalErrorsLeft );
  m_BeamCalAverageRight->addEnergies( m_BeamCalErrorsRight);


  streamlog_out(DEBUG1) << std::endl;

  for (int i = 0; i < numberForAverage;++i) {
    delete listOfBunchCrossingsLeft[i];
    delete listOfBunchCrossingsRight[i];
  }

  //Create Efficiency Objects if required
  if(m_createEfficienyFile) {
    const double //angles in mrad
      minAngle(0.9*m_BCG->getBCInnerRadius()/m_BCG->getBCZDistanceToIP()*1000), 
      maxAngle(1.1*m_BCG->getBCOuterRadius()/m_BCG->getBCZDistanceToIP()*1000); 
    const int bins = 50;
    m_thetaEfficieny = new TEfficiency("thetaEff","Efficiency vs. #Theta", bins, minAngle, maxAngle);
    m_phiEfficiency  = new TEfficiency("phiEff","Efficiency vs. #Phi", 72, 0, 360);
    m_twoDEfficiency = new TEfficiency("TwoDEff","Efficiency vs. #Theta and #Phi", bins, minAngle, maxAngle, 72, 0, 360);
    m_phiFake  = new TEfficiency("phiFake","Fake Rate vs. #Phi", 72, 0, 360);
    m_thetaFake = new TEfficiency("thetaFake","Fake Rate vs. #Theta", bins, minAngle, maxAngle);
  }//Creating Efficiency objects



}//init

void BeamCalClusterReco::processRunHeader( LCRunHeader*) {
  //  streamlog_out (DEBUG) << "Runnumber "<< _nRun << std::endl;
  //   if(_nRun % 4 == 0) {
}

void BeamCalClusterReco::processEvent( LCEvent * evt ) {

  m_random3->SetSeed(m_nEvt+Global::EVENTSEEDER->getSeed(this));

  LCCollection  *colBCal;

  try {
    colBCal = evt->getCollection( m_colNameBCal ) ;
  } catch (Exception &e) {
    colBCal = NULL;
  }

  BCPadEnergies padEnergiesLeft(m_BCG, BCPadEnergies::kLeft);
  BCPadEnergies padEnergiesRight(m_BCG, BCPadEnergies::kRight);

  ////////////////////////////////////////////////////////
  // Prepare the randomly chosen Background BeamCals... //
  ////////////////////////////////////////////////////////
  std::set<int> randomNumbers;
  unsigned int nBackgroundBX = m_backgroundBX->GetEntries();
  while( int(randomNumbers.size()) < m_nBXtoOverlay ){
    randomNumbers.insert( int(m_random3->Uniform(0, nBackgroundBX)) );
  }

  ////////////////////////
  // Sum them all up... //
  ////////////////////////
  for (std::set<int>::iterator it = randomNumbers.begin(); it != randomNumbers.end();++it) {
    m_backgroundBX->GetEntry(*it);
    padEnergiesRight.addEnergies(*m_BeamCalDepositsRight);
    padEnergiesLeft.addEnergies(*m_BeamCalDepositsLeft);
  }

  streamlog_out(MESSAGE1) << "*************** Event " << std::setw(6) << m_nEvt << " ***************" << std::endl;

  // Some event classification variables
  double depositedEnergy(0);
  double maxDeposit(0.0);
  int maxLayer(0);

  // add the energy in the event to the background/average energy
  if(colBCal) {
    CellIDDecoder<SimCalorimeterHit> mydecoder(colBCal);
    int nHits = colBCal->getNumberOfElements();
    for(int i=0; i < nHits; i++) {
      SimCalorimeterHit *bcalhit = static_cast<SimCalorimeterHit*>(colBCal->getElementAt(i));
      int side, layer, ring, sector;
      BCUtil::DecodeCellID(mydecoder, bcalhit, side, layer, ring, sector);
      const float energy = bcalhit->getEnergy();
      depositedEnergy += energy;

      if(maxDeposit < energy) {
	maxDeposit = energy;
	maxLayer = layer;
      }

      try{
	if(side == BCPadEnergies::kLeft) {
	  padEnergiesLeft.addEnergy(layer, ring, sector, energy);
	} else if(side == BCPadEnergies::kRight) {
	  padEnergiesRight.addEnergy(layer, ring, sector, energy);
	}
      } catch (std::out_of_range &e){
	streamlog_out(DEBUG4) << "Filling from signal: " << e.what()
			      << std::setw(10) << layer
			      << std::setw(10) << ring
			      << std::setw(10) << sector
			      << std::endl;
      }

    }//for all entries in the collection
  }//if there were hits from the signal

  // Run the clustering
  std::vector<BCRecoObject*> LeftSide  (FindClusters(padEnergiesLeft,  m_BeamCalAverageLeft ,  m_BeamCalErrorsLeft,  "Sig 6 L"));
  std::vector<BCRecoObject*> RightSide (FindClusters(padEnergiesRight, m_BeamCalAverageRight, m_BeamCalErrorsRight, "Sig 6 R") );

  //merge the two list of clusters so that we can run in one loop
  LeftSide.insert( LeftSide.end(), RightSide.begin(), RightSide.end() );



  if( (streamlog::out.write< DEBUG3 >() && m_nEvt == m_specialEvent ) ) {
    printBeamCalEventDisplay(padEnergiesLeft, padEnergiesRight, maxLayer, maxDeposit, depositedEnergy, LeftSide);
  }//DEBUG

  if(m_createEfficienyFile) {
    findOriginalMCParticles(evt);
    fillEfficiencyObjects(LeftSide);
  }

  /////////////////////////////////////////////////////////
  // Add the found objects to the RecoParticleCollection //
  /////////////////////////////////////////////////////////

  LCCollectionVec* BCalClusterCol = new LCCollectionVec(LCIO::CLUSTER);
  LCCollectionVec* BCalRPCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

  for (std::vector<BCRecoObject*>::iterator it = LeftSide.begin(); it != LeftSide.end(); ++it) {

    // Create Reconstructed Particles and Clusters from the BCRecoObjects" )
    const float energyCluster((*it)->getEnergy());
    const float thetaCluster((*it)->getThetaRad());
    const float phiCluster((*it)->getPhi()*TMath::DegToRad());

    const float mass = 0.0;
    const float charge = 1e+19;
    const float mmBeamCalDistance( m_BCG->getBCZDistanceToIP() );

    float position[3] = { float(mmBeamCalDistance * sin ( thetaCluster ) * cos ( phiCluster )),
			  float(mmBeamCalDistance * sin ( thetaCluster ) * sin ( phiCluster )),
			  float(mmBeamCalDistance * cos ( thetaCluster )) };

    float momentumCluster[3] = { float(energyCluster * sin ( thetaCluster ) * cos ( phiCluster )),
				 float(energyCluster * sin ( thetaCluster ) * sin ( phiCluster )),
				 float(energyCluster * cos ( thetaCluster )) };
#pragma message "FIXME: Rotate position and momentum of cluster to lab system"

    ClusterImpl* cluster = new ClusterImpl;
    cluster->setEnergy( energyCluster );
    cluster->setPosition( position );

    ReconstructedParticleImpl* particle = new ReconstructedParticleImpl;
    particle->setMass( mass ) ;
    particle->setCharge( charge ) ;
    particle->setMomentum ( momentumCluster ) ;
    particle->setEnergy ( energyCluster ) ;
    particle->addCluster( cluster ) ;

    BCalClusterCol->addElement(cluster);
    BCalRPCol->addElement(particle);

    //CleanUp
    delete *it;
  }//for all found clusters

  ///////////////////////////////////////
  // Done with Running Reco Clustering //
  ///////////////////////////////////////
  
  if( BCalClusterCol->getNumberOfElements() != 0 ) {
    evt->addCollection(BCalClusterCol, m_BCalClusterColName);
    evt->addCollection(BCalRPCol, m_BCalRPColName);
  } else {
    delete BCalClusterCol;
    delete BCalRPCol;
  }

  m_nEvt++ ;


}//processEvent


void BeamCalClusterReco::fillEfficiencyObjects(const std::vector<BCRecoObject*>& RecoedObjects) {

  if( not m_createEfficienyFile ) return;

  //Try to match all clusters and particles, only then can we fill our efficiencies,
  //this should allow one in principle to estimate efficiencies when there are
  //multiple particles or cluster, but should be checked!!

  //See if we have a cluster matching one of the MCParticles
  for (std::vector<BCRecoObject*>::const_iterator it = RecoedObjects.begin(); it != RecoedObjects.end(); ++it) {
    BCRecoObject* bco = *it;

    bool hasRightCluster = false;
    const double theta(bco->getThetaMrad());
    const double phi(bco->getPhi());

    for (std::vector<OriginalMC>::iterator mcIt = m_originalParticles.begin(); mcIt != m_originalParticles.end(); ++mcIt) {
      if (BCUtil::areCloseTogether(theta, phi, (*mcIt).m_theta, (*mcIt).m_phi )) {
	hasRightCluster = true;
	(*mcIt).m_wasFound = true;
	break;
      }
    }

    bco->setHasRightCluster(hasRightCluster);
    streamlog_out(MESSAGE2) << "Have we found a cluster matching a particle? "
			    << std::boolalpha << hasRightCluster 
			    << std::endl;
  }


  //Here we fill the efficiency for reconstructing MCParticles
  for (std::vector<OriginalMC>::iterator mcIt = m_originalParticles.begin(); mcIt != m_originalParticles.end(); ++mcIt) {
    OriginalMC const& omc = (*mcIt);
    m_thetaEfficieny->Fill( omc.m_wasFound, omc.m_theta);
    m_phiEfficiency->Fill ( omc.m_wasFound, omc.m_phi);
    m_twoDEfficiency->Fill( omc.m_wasFound, omc.m_theta, omc.m_phi);
    streamlog_out(MESSAGE2) << "Particle was found? " << std::boolalpha << omc.m_wasFound 
			    << std::setw(13) << omc.m_theta
			    << std::setw(13) << omc.m_phi
			    << std::setw(13) << omc.m_energy
			    << std::endl;
  }

  //Here we fill the fake rate, for reconstructed clusters that do not have an MCParticle
  bool foundFake(false);
  for (std::vector<BCRecoObject*>::const_iterator it = RecoedObjects.begin(); it != RecoedObjects.end(); ++it) {
    BCRecoObject* bco = *it;
    const bool hasRightCluster(bco->hasRightCluster());
    const double theta(bco->getThetaMrad());
    const double phi(bco->getPhi());
    if (not hasRightCluster) {
      m_thetaFake->Fill(true, theta);
      m_phiFake->Fill(true, phi);
      foundFake = true;
    }
  }

  if (not foundFake) {
    const int nbins   = m_thetaFake->GetTotalHistogram()->GetNbinsX();
    const double low  = m_thetaFake->GetTotalHistogram()->GetBinLowEdge(1);
    const double high = m_thetaFake->GetTotalHistogram()->GetBinLowEdge(nbins+1);
    const double step = (high-low)/double(nbins);
    for (int i = 1; i <= nbins ;++i) {
      m_thetaFake->Fill( false, i * step + step/2.0 + low );
    }
  }

}//fillEfficiencyObjects




void BeamCalClusterReco::check( LCEvent * ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void BeamCalClusterReco::end(){

  streamlog_out ( MESSAGE4 ) << __PRETTY_FUNCTION__ << " " << name()
			     << " processed " << m_nEvt << " events."
			     << std::endl ;


  if(m_createEfficienyFile) {
    TFile *effFile = TFile::Open(m_EfficiencyFileName.c_str(),"RECREATE");
    m_thetaEfficieny->Write();
    m_phiEfficiency->Write();
    m_twoDEfficiency->Write();
    m_phiFake->Write();
    m_thetaFake->Write();
    effFile->Close();
  }

  delete m_BeamCalAverageLeft;
  delete m_BeamCalAverageRight;
  delete m_BeamCalErrorsLeft;
  delete m_BeamCalErrorsRight;

  delete m_random3;
  delete m_BCG;
  delete m_bcpCuts;

}

/// Calculate the one sigma errors for the given selected background bunch crossings
BCPadEnergies* BeamCalClusterReco::getBeamCalErrors(const BCPadEnergies *averages, const std::vector<BCPadEnergies*> singles, int numberForAverage ) {

  BCPadEnergies * BCPErrors = new BCPadEnergies(m_BCG);
  for (int i = 0; i < m_BCG->getPadsPerBeamCal()  ;++i) {
    Double_t mean(averages->getEnergy(i));
    Double_t variance(0);
    Double_t nHistos(numberForAverage);
    for (int l = 0; l < nHistos;++l) {
      Double_t energy( singles[l]->getEnergy(i) );
      variance += ( energy - mean ) * ( energy - mean );
    }//histograms
    variance /= nHistos;
    variance = sqrt(variance);
    BCPErrors->setEnergy(i, variance);

  }//for all pads

  return BCPErrors;

}//getBeamCalErrors



std::vector<BCRecoObject*> BeamCalClusterReco::FindClusters(const BCPadEnergies& signalPads,
							    const BCPadEnergies& backgroundPads,
							    const BCPadEnergies& backgroundSigma,
							    const TString& title) {

  std::vector<BCRecoObject*> recoVec;

  //////////////////////////////////////////
  // This calls the clustering function!
  //////////////////////////////////////////
  const std::vector<BeamCalCluster> &bccs =
    signalPads.lookForNeighbouringClustersOverWithVetoAndCheck(backgroundPads, backgroundSigma, *m_bcpCuts);
  const bool isRealParticle = false; //always false here, decide later

  for (std::vector<BeamCalCluster>::const_iterator it = bccs.begin(); it != bccs.end(); ++it) {

    streamlog_out(MESSAGE2) << title;
    if(signalPads.getSide() == BCPadEnergies::kRight) streamlog_out(MESSAGE2) << LONGSTRING;
    streamlog_out(MESSAGE2) << " " << (*it);
    if(signalPads.getSide() == BCPadEnergies::kLeft) streamlog_out(MESSAGE2) << LONGSTRING;

    //Apply cuts on the reconstructed clusters, then calculate angles
    if ( ( it->getNPads() > 2 ) && m_bcpCuts->isClusterAboveThreshold( (*it) ) ) {

      double theta(it->getTheta());
      double phi  (it->getPhi());

      streamlog_out(MESSAGE2) << " found something "
			      << std::setw(10) << theta
			      << std::setw(10) << phi
	;//ending the streamlog!

      recoVec.push_back( new BCRecoObject(isRealParticle, true, theta, phi, it->getEnergy(), it->getNPads(), signalPads.getSide() ) );

    }//if we have enough pads and energy in the clusters

    //Finish the output line
    streamlog_out(MESSAGE2) << std::endl;

  }//clusterloop

  return recoVec;

}//tryReco6


void BeamCalClusterReco::printBeamCalEventDisplay(BCPadEnergies& padEnergiesLeft, BCPadEnergies& padEnergiesRight,
						  int maxLayer, double maxDeposit, double depositedEnergy,
						  const std::vector<BCRecoObject*> & RecoedObjects) const {

  BCPadEnergies *padEnergies, *padErrors, *padAverages;
  if( m_eventSide == BCPadEnergies::kLeft ) {
    padEnergies = &padEnergiesLeft;
    padAverages = m_BeamCalAverageLeft;
    padErrors = m_BeamCalErrorsLeft;
  } else {
    padEnergies = &padEnergiesRight;
    padAverages = m_BeamCalAverageRight;
    padErrors = m_BeamCalErrorsRight;
  }

  ///////////////////////////////////////
  // Deal with the Canvas and its pads //
  ///////////////////////////////////////
  TCanvas canv("canv","canv", 3200, 1600);
  TPad *pads[9];
  Double_t padSizeX = canv.GetWindowWidth()/4;
  Double_t padSizeY = (canv.GetWindowHeight()-200)/2;
  for (int i = 0; i < 8 ;++i) {
    Double_t lowX, lowY, highX, highY;

    lowX  = padSizeX *	double(int(i%4))    / 3200.;
    highX = padSizeX *	double(int(i%4)+1)  / 3200.;

    lowY  = (padSizeY *	double(int(i/4))   + 200.) / 1600.;
    highY = (padSizeY *	double(int(i/4)+1) + 200.) / 1600.;

    pads[i] = new TPad(Form("Pad%i",i+1), Form("Pad%i",i+1), lowX, lowY, highX, highY,
		       kWhite, short(0), short(0));
    pads[i]->SetNumber( ( i+1 < 5) ? i+1 + 4 : i+1 - 4 );//to start at the top... DONT TOUCH THIS EVER!
    pads[i]->Draw();
  }
  pads[8] = new TPad("TextPad","Text", 0, 0., 1., 1./9., kWhite, 0, 0);
  pads[8]->SetNumber(9);
  pads[8]->Draw();

  //Draw one before and two after maxLayer
  if(maxDeposit < 0.1) maxLayer = 10;
  const int startLayer = (maxLayer - 1 <= 37) ? maxLayer - 1 : 37 ;
  ////////////////////////
  // Deal with the data //
  ////////////////////////

  BeamCal bc(marlin::Global::GEAR);
  bc.SetLogz(1);
  bc.SetBeamCalHisto(padEnergies,"tempLeft");

  Double_t ymax = 5;


  for (int layer = startLayer; layer < startLayer + 4; ++layer) {

    TH2F frame("frame",Form("BeamCal Layer %i", layer), 160, -160, 160, 160, -160, 160);

    const int pad1 = ( layer - startLayer ) + 1;
    const int pad2 = ( layer - startLayer ) + 5;

    //------------------------------
    // Pad 1
    //------------------------------
    canv.cd(pad1);
    canv.cd(pad1)->SetRightMargin(0.18);
    canv.cd(pad1)->SetLeftMargin (0.18);
    bc.BeamCalDraw((TPad*)canv.GetPad(pad1), &frame, layer);
    gPad->Update();
    DrawElectronMarkers( RecoedObjects );

    {
      TCanvas tCanv("temp1",Form("BeamCal Layer %i", layer));
      tCanv.cd();
      bc.BeamCalDraw((TPad*)tCanv.GetPad(0), &frame, layer);
      tCanv.SaveAs(Form("SpecialEvent_Pad%i.eps", pad1));
    }

    //------------------------------
    // Pad 2
    //------------------------------
    canv.cd(pad2);
    bc.SetBeamCalHisto(padAverages, padErrors);
    bc.DrawPhiDistributions((TPad*)canv.GetPad(pad2), layer, "dotted,errors");

    bc.SetBeamCalHisto(padEnergies,"padLeft");
    bc.DrawPhiDistributions((TPad*)canv.GetPad(pad2), layer, "same");
    static_cast<TH1*>(((TPad*)canv.GetPad(pad2))->GetListOfPrimitives()->At(0))->SetMaximum(ymax);
    DrawLineMarkers( RecoedObjects );

    {
      TCanvas tCanv("temp1",Form("PhiDistribution Layer %i", layer));
      tCanv.cd();
      bc.SetBeamCalHisto(padAverages, padErrors);
      bc.DrawPhiDistributions(&tCanv, layer, "dotted,errors");
      bc.SetBeamCalHisto(padEnergies,"padLeft");
      bc.DrawPhiDistributions(&tCanv, layer, "same");
      static_cast<TH1*>(tCanv.GetListOfPrimitives()->At(0))->SetMaximum(ymax);
      //DrawLineMarkers( RecoedObjects );
      tCanv.SaveAs(Form("SpecialEvent_Pad%i.eps", pad2));
    }

  }//run over layers

  //Write some information
  canv.cd(9);
  TPaveText text(0.0, 0.0, 1.0, 1.0);
  text.SetFillColor(kWhite);
  text.AddText(Form("ImpactAngle #theta: %2.1f, Energy Deposit: %2.1f GeV", -1.0, depositedEnergy));
  text.AddText(Form("Max Deposit: %2.1f GeV in Layer %i", maxDeposit, maxLayer));
  text.Draw();

  canv.SaveAs(Form("Event%i.eps", m_nEvt));

}


void BeamCalClusterReco::DrawElectronMarkers ( const std::vector<BCRecoObject*> & RecoedObjects ) const {

  const double BeamCalDist = m_BCG->getBCZDistanceToIP();

  for( std::vector<BCRecoObject*>::const_iterator it = RecoedObjects.begin();
       it != RecoedObjects.end(); ++it) {

    double radius = BeamCalDist*tan((*it)->getThetaRad() );
    double circX = radius*cos((*it)->getPhi()*TMath::DegToRad());
    double circY = radius*sin((*it)->getPhi()*TMath::DegToRad());
    // electron.SetNextPoint(radius*cos(m_impactAnglePhi*TMath::DegToRad()),
    //			    radius*sin(m_impactAnglePhi*TMath::DegToRad()));

    TCrown electron(circX, circY, 25, 30);
    electron.SetLineColor(kRed);
    electron.SetFillColor(kRed);
    // electron.SetFillStyle(4000);
    // electron.SetMarkerStyle(kOpenCircle);
    // electron.SetMarkerSize(5);
    // electron.SetMarkerColor(kRed);
    //      Double_t ymin = 0, ymax = 35;
    electron.Draw();

  }

  return;
}

void BeamCalClusterReco::DrawLineMarkers( const std::vector<BCRecoObject*> & RecoedObjects ) const {

  Double_t ymin = 0, ymax = 5;

  for( std::vector<BCRecoObject*>::const_iterator it = RecoedObjects.begin();
       it != RecoedObjects.end(); ++it) {
    TLine line((*it)->getPhi(),ymin,(*it)->getPhi(),ymax);
    line.SetLineStyle(kDashed);
    line.SetLineColor(kRed);
    line.SetLineWidth(0);
    line.Draw();//on gPad whatever active
  }

  return;
}


void BeamCalClusterReco::findOriginalMCParticles(LCEvent *evt) {
  m_originalParticles.clear();
  const double maxAngle(1.1*m_BCG->getBCOuterRadius()/m_BCG->getBCZDistanceToIP()*1000); 
  try {
    LCCollection* colMC = evt->getCollection ( m_colNameMC );
    for (int i =0; i < colMC->getNumberOfElements();++i) {
      MCParticle *tempmc = dynamic_cast<MCParticle*>(colMC->getElementAt(i));

      //only use particles generated by generator
      //this does not work, if we used the particle gun, in which case everything was created in simulation
      // in this case we take only the first
      if (not ( tempmc->getGeneratorStatus() == 1 || (tempmc->getGeneratorStatus() == 0 && i == 0 ) ) ) continue;

      const double *momentum = tempmc->getMomentum();
      double momentum2[3];
    
      //Check if it inside the BeamCalAcceptance, and has sufficientmomentum(1GeV);
      const double absMom = sqrt(momentum[0]*momentum[0]+momentum[1]*momentum[1]+momentum[2]*momentum[2]);
      if ( absMom < 10 ) continue;
      //Rotate to appropriate beamCal System
#pragma message "Fix the rotation thing"
      if(momentum[2] > 0) {
	BCUtil::RotateToBeamCal<10>(momentum, momentum2);
	m_eventSide = 0;
      } else {
	BCUtil::RotateFromBeamCal<10>(momentum, momentum2);
	m_eventSide = 1;
      }
      //    double radius = sqrt(momentum2[0]*momentum2[0]+momentum2[1]*momentum2[1]);
      double impactTheta = BCUtil::AngleToBeamCal<10>(momentum)*1000;//mrad
      if ( impactTheta > maxAngle ) continue; //only particles that are inside the BeamCal acceptance

      double impactPhi   = TMath::ATan2(momentum2[1], momentum2[0]) * TMath::RadToDeg();
      if(impactPhi < 0) impactPhi += 360;
      
      OriginalMC omc(impactTheta, impactPhi, absMom);
      m_originalParticles.push_back(omc);
      streamlog_out(MESSAGE2) << "Found MCParticle: (Theta, Phi, Eng) "
			      << std::setw(13) << impactTheta
			      << std::setw(13) << impactPhi
			      << std::setw(13) << absMom
			      << std::endl;
    }
  } catch (Exception &e) {
  }

  return;
}//FindOriginalMCParticle
