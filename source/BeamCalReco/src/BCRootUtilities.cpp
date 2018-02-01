#include "BCRootUtilities.hh"
#include "BCPadEnergies.hh"
#include "BeamCalGeo.hh"

#include <TFile.h>
#include <TString.h>
#include <TTree.h>

#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>


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



void BCUtil::ReadBecasFile(std::string const& fileName, std::vector<BCPadEnergies>& newPads,
			   std::string treeName,
			   std::string energyField,
			   bool hasPairMonitor
			   ) {

  TFile *becasFile = TFile::Open(TString(fileName));
  if( not becasFile ) {
    std::cerr << "File not found! " << fileName  << std::endl;
    throw std::invalid_argument("File does not exist!");
  }

  TTree *becasTree(NULL);
  becasFile->GetObject(treeName.c_str(), becasTree);
  if ( not becasTree ) {
    std::cerr << "Tree not found in file " << fileName  << std::endl;
    becasFile->Close();
    delete becasFile;
    throw std::invalid_argument("The file does not contain a BeCas tree");
  }

  int position[3];
  double energy;
  becasTree->SetBranchAddress("sPos",position);
  becasTree->SetBranchAddress(energyField.c_str(), &energy);

  const BeamCalGeo* geo = &(newPads[0].m_BCG);

  for (int i = 0 ; i < becasTree->GetEntries(); ++i)  {
    becasTree->GetEntry(i);
    int layer(position[0]), ring(position[1]), phiPad(position[2]);
    if(hasPairMonitor and std::abs(layer) > geo->getBCLayers() ){
      continue;
    }

  try {
    const int globalPadIndex = geo->getPadIndex( std::abs(layer), ring, phiPad);
    BCPadEnergies::BeamCalSide_t side = ( layer < 0 ) ? BCPadEnergies::kLeft : BCPadEnergies::kRight ;
    newPads[side].addEnergy( globalPadIndex, energy);
  } catch (std::exception &e)  {
    std::cout << "Exception " << e.what()  << std::endl;
    std::cout << "MaxValues (Ring, Layer, Pads)"
	      << std::setw(10) << geo->getBCLayers()
	      << std::setw(10) << geo->getBCRings()
	      << std::setw(10) << geo->getNSegments()[ring]
	      << std::endl;
    std::cout << "Found (Ring, Layer, Pads)    "
	      << std::setw(10) << layer
	      << std::setw(10) << ring
	      << std::setw(10) << phiPad
	      << std::endl;
    throw;
  }

  }//for all entries

  becasFile->Close();
  return;

}
