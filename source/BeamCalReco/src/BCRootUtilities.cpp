#include "BCRootUtilities.hh"
#include "BeamCalGeo.hh"


#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <cmath>
#include <stdexcept>
#include <cstdlib>


void BCUtil::ReadRootFile(std::string const& fileName, std::vector<BCPadEnergies>& newPads) {

  TFile* file = TFile::Open(fileName.c_str());
  if ( not file ) {
    std::cerr << "File not found: "<< fileName  << std::endl;
    return;
  }
  TTree* tree;
  file->GetObject("bcTree", tree);
  if ( not tree ) {
    std::cerr << "Tree not found in file " << fileName  << std::endl;
    file->Close();
    delete file;
    throw std::invalid_argument("The file does not contain a BCPadenergies Rootfile");
  }

  std::vector<double> *depositsLeft(NULL), *depositsRight(NULL);
  tree->SetBranchAddress("vec_left" , &depositsLeft);
  tree->SetBranchAddress("vec_right", &depositsRight);
  tree->GetEntry(0);
  newPads[0].addEnergies(*depositsLeft);
  newPads[1].addEnergies(*depositsRight);

  file->Close();
  delete file;
  return;
}



std::vector<BCPadEnergies>* BCUtil::ReadMultiRootFile(std::string const& fileName, BeamCalGeo* geoCache) {

  TFile* file = TFile::Open(fileName.c_str());
  if ( not file ) {
    std::cerr << "File not found: "<< fileName  << std::endl;
    return NULL;
  }
  TTree* tree;
  file->GetObject("bcTree", tree);
  if ( not tree ) {
    std::cerr << "Tree not found in file " << fileName  << std::endl;
    file->Close();
    delete file;
    throw std::invalid_argument("The file does not contain a BCPadenergies Rootfile");
  }

  std::vector<double> *depositsLeft(NULL), *depositsRight(NULL);
  tree->SetBranchAddress("vec_left" , &depositsLeft);
  tree->SetBranchAddress("vec_right", &depositsRight);
  int nEntries = tree->GetEntries();

  std::vector<BCPadEnergies> *newPads = new std::vector<BCPadEnergies>(2*nEntries, geoCache);

  for( int iEvt=0; iEvt<nEntries; iEvt++)
  {
	  tree->GetEntry(iEvt);

	  std::cout << "Processing entry #" << iEvt << "\r";

	  newPads->at(iEvt*2).setEnergies(*depositsLeft);
	  newPads->at(iEvt*2+1).setEnergies(*depositsRight);
  }

  file->Close();
  delete file;
  return newPads;
}



void BCUtil::ReadBecasFile(std::string const& fileName, std::vector<BCPadEnergies>& newPads) {

  TFile *becasFile = TFile::Open(TString(fileName));
  if( not becasFile ) {
    std::cerr << "File not found! " << fileName  << std::endl;
    return;
  }

  TTree *becasTree(NULL);
  becasFile->GetObject("tSegment",becasTree);
  if ( not becasTree ) {
    std::cerr << "Tree not found in file " << fileName  << std::endl;
    becasFile->Close();
    delete becasFile;
    throw std::invalid_argument("The file does not contain a BeCas file");
  }

  int position[3];
  double energy;
  becasTree->SetBranchAddress("sPos",position);
  becasTree->SetBranchAddress("sEdep",&energy);

  const BeamCalGeo* geo = &(newPads[0].m_BCG);

  for (int i = 0 ; i < becasTree->GetEntries(); ++i)  {
    becasTree->GetEntry(i);
    int layer(position[0]), ring(position[1]), phiPad(position[2]);
    const int globalPadIndex = geo->getPadIndex( std::abs(layer), ring, phiPad);
    BCPadEnergies::BeamCalSide_t side = ( layer < 0 ) ? BCPadEnergies::kLeft : BCPadEnergies::kRight ;
    newPads[side].addEnergy( globalPadIndex, energy);

  }//for all entries

  becasFile->Close();
  return;

}
