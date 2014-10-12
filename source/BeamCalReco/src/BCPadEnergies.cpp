#include "BCPadEnergies.hh"
#include "BeamCalCluster.hh"
#include "BCPCuts.hh"
#include "BeamCalGeo.hh"

#include <TH2D.h>
#include <TROOT.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

// ----- include for verbosity dependend logging ---------
#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>

//typedef std::map<int, int> PadIndexList;
// typedef BCPadEnergies::PadIndexList;
// typedef BCPadEnergies::TowerIndexList;



bool value_comparer( const BCPadEnergies::TowerIndexList::value_type &i1, const BCPadEnergies::TowerIndexList::value_type &i2)
{
  return i1.second<i2.second;
}

BCPadEnergies::BCPadEnergies(const BeamCalGeo &bcg, BeamCalSide_t side):
  m_PadEnergies(bcg.getPadsPerBeamCal()),
  m_side(side),
  m_BCG(bcg)
{

}

BCPadEnergies::BCPadEnergies(const BeamCalGeo *bcg, BeamCalSide_t side):
  m_PadEnergies(bcg->getPadsPerBeamCal()),
  m_side(side),
  m_BCG(*bcg)
{

}

BCPadEnergies::BCPadEnergies(const BCPadEnergies &bcp):
  m_PadEnergies(bcp.m_PadEnergies),
  m_side(bcp.m_side),
  m_BCG(bcp.m_BCG)
{
}

BCPadEnergies::BCPadEnergies(const BCPadEnergies *bcp):
  m_PadEnergies(bcp->m_PadEnergies),
  m_side(bcp->m_side),
  m_BCG(bcp->m_BCG)
{
}

std::vector<double>* BCPadEnergies::getEnergies() { return &m_PadEnergies; }

double BCPadEnergies::getEnergy(int layer, int ring, int pad) const {
  return m_PadEnergies[ m_BCG.getPadIndex(layer, ring, pad) ];
}

double BCPadEnergies::getEnergy(int padIndex) const {
  return m_PadEnergies[ padIndex ];
}


void BCPadEnergies::setEnergy(int layer, int ring, int pad, double energy){
  m_PadEnergies[ m_BCG.getPadIndex(layer, ring, pad) ] = energy;
}

void BCPadEnergies::setEnergy(int padIndex, double energy){
  m_PadEnergies[ padIndex ] = energy;
}

void BCPadEnergies::resetEnergies(){
  for (int i = 0; i < m_BCG.getPadsPerBeamCal();++i) {
    m_PadEnergies[i] = 0 ;
  }
}

void BCPadEnergies::scaleEnergies(double factor){
  for (int i = 0; i < m_BCG.getPadsPerBeamCal();++i) {
    m_PadEnergies[i] *= factor;
  }
}//Scale

double BCPadEnergies::getTotalEnergy() const{
  double sum(0.0);
  for (int i = 0; i < m_BCG.getPadsPerBeamCal();++i) {
    sum += m_PadEnergies[i];
  }
  return sum;
}



void BCPadEnergies::setEnergies(const std::vector<double> &energies){
  if( (int)energies.size() != m_BCG.getPadsPerBeamCal() ) {
    std::stringstream errorMessage;
    errorMessage << "Energies vector has wrong size! " << energies.size() << " vs. "  <<  m_BCG.getPadsPerBeamCal();
    throw std::out_of_range( errorMessage.str()  );
  }
  for (int i = 0; i < m_BCG.getPadsPerBeamCal();++i) {
    m_PadEnergies[i] = energies[i] ;
  }
  //  m_PadEnergies = energies;
}

void BCPadEnergies::setEnergies(const BCPadEnergies &bcp){
  if( bcp.m_BCG.getPadsPerBeamCal() != m_BCG.getPadsPerBeamCal()) throw std::out_of_range("BCPadEnergies has wrong size!");
  for (int i = 0; i < m_BCG.getPadsPerBeamCal();++i) {
    m_PadEnergies[i] = bcp.m_PadEnergies[i] ;
  }
  //  m_PadEnergies = energies;
}


void BCPadEnergies::addEnergies(const std::vector<double> &energies){
  if( (int)energies.size() != m_BCG.getPadsPerBeamCal() ) {
    std::stringstream errorMessage;
    errorMessage << "Energies vector has wrong size! " << energies.size() << " vs. "  <<  m_BCG.getPadsPerBeamCal();
    throw std::out_of_range( errorMessage.str()  );
  }
  for (int i = 0; i < m_BCG.getPadsPerBeamCal();++i) {
    m_PadEnergies[i] += energies[i] ;
  }
}

void BCPadEnergies::addEnergies(const BCPadEnergies &bcp){
  if( bcp.m_BCG.getPadsPerBeamCal() != m_BCG.getPadsPerBeamCal() ) throw std::out_of_range("BCPadEnergies has wrong size!");
  for (int i = 0; i < m_BCG.getPadsPerBeamCal();++i) {
    m_PadEnergies[i] += bcp.m_PadEnergies[i] ;
  }
}



void BCPadEnergies::subtractEnergies(const std::vector<double> &energies){
  if( (int)energies.size() != m_BCG.getPadsPerBeamCal() ) {
    std::stringstream errorMessage;
    errorMessage << "Energies vector has wrong size! " << energies.size() << " vs. "  <<  m_BCG.getPadsPerBeamCal();
    throw std::out_of_range( errorMessage.str()  );
  }
  for (int i = 0; i < m_BCG.getPadsPerBeamCal();++i) {
    m_PadEnergies[i] -= energies[i] ;
  }
}


void BCPadEnergies::subtractEnergies(const BCPadEnergies &bcp){
  if( bcp.m_BCG.getPadsPerBeamCal() != m_BCG.getPadsPerBeamCal()) throw std::out_of_range("BCPadEnergies has wrong size!");
  for (int i = 0; i < m_BCG.getPadsPerBeamCal();++i) {
    m_PadEnergies[i] -= bcp.m_PadEnergies[i] ;
  }
}



/**
 * Checks if a certain percentage of the pads in the first ring of layer 10 are
 * significantly above the average energy deposits of the background, and in
 * this case, we reduce the background by one sigma. If there are 25 or more
 * pads below -0.9 sigma, we increase in addEnergiesWithCheck
 */
void BCPadEnergies::subtractEnergiesWithCheck(const BCPadEnergies &bcp, const BCPadEnergies &sigma){
  int tooMuchAbove = 0, tooMuchBelow = 0;
  if( bcp.m_BCG.getPadsPerBeamCal() != m_BCG.getPadsPerBeamCal()) throw std::out_of_range("BCPadEnergies has wrong size!");

  for (int i = 0; i < m_BCG.getPadsPerBeamCal();++i) {

    if( &bcp != &sigma ) {
      m_PadEnergies[i] -= bcp.m_PadEnergies[i] ;
    } else {
      m_PadEnergies[i] -= 0.10 * bcp.m_PadEnergies[i] ;
    }

    if( (m_BCG.getLayer(i) == 10) &&
	(m_BCG.getRing(i) == 0) ) {
      if ( m_PadEnergies[i] > 0.9 * sigma.m_PadEnergies[i])  {
	tooMuchAbove++;
      } else if( m_PadEnergies[i] < -0.9 * sigma.m_PadEnergies[i])  {
	tooMuchBelow++;
      }
    }
  }
  //if it is above the limit, do it again
  if( tooMuchAbove >= 5 ) {
    subtractEnergiesWithCheck(sigma, sigma);
  } else if (tooMuchBelow >= 25) {
    addEnergiesWithCheck(sigma, sigma);
  }

}//subtractEnergiesWithCheck

/**
 * Increase by 10% of sigma the energy deposits
 */
void BCPadEnergies::addEnergiesWithCheck(const BCPadEnergies &bcp, const BCPadEnergies &sigma){
  int tooMuchBelow = 0;
  if( bcp.m_BCG.getPadsPerBeamCal() != m_BCG.getPadsPerBeamCal()) throw std::out_of_range("BCPadEnergies has wrong size!");

  for (int i = 0; i < m_BCG.getPadsPerBeamCal();++i) {

    m_PadEnergies[i] += 0.10 * bcp.m_PadEnergies[i] ;

    if( (m_BCG.getLayer(i) == 10) &&
	(m_BCG.getRing(i) == 0) &&
	( m_PadEnergies[i] < -0.9 * sigma.m_PadEnergies[i]) ) {
      tooMuchBelow++;
    }
  }
  //if it is above the limit, do it again
  if( tooMuchBelow >= 25 ) {
    addEnergiesWithCheck(sigma, sigma);
  }
}//addEnergiesWithCheck


void BCPadEnergies::addEnergy(int layer, int ring, int pad, double energy){
  m_PadEnergies[ m_BCG.getPadIndex(layer, ring, pad) ] += energy;
}

void BCPadEnergies::addEnergy(int padIndex, double energy){
  m_PadEnergies[ padIndex ] += energy;
}


//Take padEnergies and subtract the averages and errors from it, then print all
//remaining clusters with energy above 0.5 GeV behind layer10.  Why not just
//subtract averages, and the user can supply averages and errors like he wants?
//The function makes a copy of the pads so different approaches can be tried
BeamCalCluster BCPadEnergies::lookForClustersOver(const BCPadEnergies &background, const BCPCuts &cuts) const {
  BCPadEnergies testPads(*this);

  //remove the average
  testPads.subtractEnergies(background);

  const PadIndexList& myPadIndices = getPadsAboveThresholds(testPads, cuts);

  BeamCalCluster BCCluster = getClusterFromAcceptedPads( testPads,  myPadIndices, cuts );

  return BCCluster;

} // lookForClusters

//The function makes a copy of the pads so different approaches can be tried
//For every position at which we have a cluster we count the number of pads, then we take those where we have the most
BeamCalCluster BCPadEnergies::lookForAlignedClustersOver(const BCPadEnergies &background, const BCPCuts &cuts) const {
  BCPadEnergies testPads(*this);

  //remove the average
  testPads.subtractEnergies(background);


  PadIndexList myPadIndices = getPadsAboveThresholds(testPads, cuts);
  TowerIndexList myTowerIndices = getTowersFromPads( this->m_BCG, myPadIndices );

  //Get Pointer to max element...
  TowerIndexList::iterator largestTower = std::max_element(myTowerIndices.begin(), myTowerIndices.end(), value_comparer);

  //if we dont have at least 4 of them aligned we leave this place
  if( largestTower->second < cuts.getMinimumTowerSize() )  return BeamCalCluster();

  //Now throw away all but the ones in the largestTower
  for (TowerIndexList::iterator it = myTowerIndices.begin(); it != myTowerIndices.end(); ++it) {
    if( it->first == largestTower->first) continue;
    removeTowerFromPads( this->m_BCG, myPadIndices, it->first );
  }

  BeamCalCluster BCCluster = getClusterFromAcceptedPads(testPads, myPadIndices, cuts);

  //set the padIndexIn the layer to let BeamCal calculate the angles from it
  BCCluster.setPadIndexInLayer(largestTower->first);

  return BCCluster;

} // lookForAlignedClustersOver


//The function makes a copy of the pads so different approaches can be tried
BeamCalCluster BCPadEnergies::lookForNeighbouringClustersOver(const BCPadEnergies &background, const BCPCuts &cuts) const {
  BCPadEnergies testPads(*this);

  //remove the average
  testPads.subtractEnergies(background);

  //here cuts are applied on the pads
  PadIndexList myPadIndices = getPadsAboveThresholds(testPads, cuts);

  TowerIndexList myTowerIndices = getTowersFromPads( this->m_BCG, myPadIndices );

  // only if we have just 1 tower we can just use the other function...
  if( myTowerIndices.size() == 1 )  return lookForAlignedClustersOver(background, cuts);

  //Compare the largest deposit to the others and check if they are neighbours
  TowerIndexList::iterator largestTower = std::max_element(myTowerIndices.begin(), myTowerIndices.end(),
							   value_comparer);

  for (TowerIndexList::iterator it = myTowerIndices.begin(); it != myTowerIndices.end(); ++it) {
    if( it->first == largestTower->first) continue;

    const bool isNeighbour = m_BCG.arePadsNeighbours(largestTower->first, it->first);
    if( not isNeighbour ) {
      //remove padIndices from myPadIndices
      removeTowerFromPads ( this->m_BCG, myPadIndices, it->first);
    }

  }// find neighbouring towers/pads

  BeamCalCluster BCCluster = getClusterFromAcceptedPads(testPads, myPadIndices, cuts);
  //#pragma "FIXME:Use average over neighbouring towers"
  //Should probably use an average over all the towers
  BCCluster.setPadIndexInLayer(largestTower->first);
  return BCCluster;
} // lookForNeighbouringClustersOver


BeamCalCluster BCPadEnergies::lookForNeighbouringClustersOverWithVeto(const BCPadEnergies &background, const BCPCuts &cuts) const {
  BCPadEnergies testPads(*this);
  testPads.subtractEnergies(background);

  //here cuts are applied on the pads
  PadIndexList myPadIndices = getPadsAboveThresholds(testPads, cuts);
  TowerIndexList myTowerIndices = getTowersFromPads( this->m_BCG, myPadIndices );

  // only if we have just 1 tower we can just use the other function...
  if( myTowerIndices.size() == 1 )  return lookForAlignedClustersOver(background, cuts);

  //Compare the largest deposit to the others and check if they are neighbours
  TowerIndexList::iterator largestTower = std::max_element(myTowerIndices.begin(), myTowerIndices.end(),
							   value_comparer);

  if(largestTower->second < cuts.getMinimumTowerSize()) return BeamCalCluster();

  for (TowerIndexList::iterator it = myTowerIndices.begin(); it != myTowerIndices.end(); ++it) {
    if( it->first == largestTower->first) continue;

    const bool isNeighbour = m_BCG.arePadsNeighbours(largestTower->first, it->first);

    if( not isNeighbour ) {
      //remove padIndices from myPadIndices
      removeTowerFromPads ( this->m_BCG, myPadIndices, it->first);
    } else if( it->second < cuts.getMinimumTowerSize() ){
      removeTowerFromPads ( this->m_BCG, myPadIndices, it->first);
    }
  }// find neighbouring towers/pads

  BeamCalCluster BCCluster = getClusterFromAcceptedPads(testPads, myPadIndices, cuts);
  BCCluster.setPadIndexInLayer(largestTower->first);

  return BCCluster;
} // lookForNeighbouringClustersOverWithVeto



BCPadEnergies::BeamCalClusterList BCPadEnergies::listOfNeighbouringClustersOverWithVeto(const BCPadEnergies &background,
												 const BCPadEnergies &backgroundSigma,
												 const BCPCuts &cuts) const {
  BCPadEnergies::BeamCalClusterList BeamCalClusters;

  //We make a copy, because we might want to apply different clustering on the same pads
  BCPadEnergies testPads(*this);
  testPads.subtractEnergies(background);
  //here cuts are applied on the pads
  PadIndexList myPadIndices = ( cuts.useConstPadCuts() ) ?
    getPadsAboveThresholds(testPads, cuts) :
    testPads.getPadsAboveSigma(backgroundSigma, cuts);

  ClusterNextToNearestNeighbourTowers(testPads, myPadIndices, cuts, BeamCalClusters);

  return BeamCalClusters;
} // listOfNeighbouringClustersOverWithVeto



BCPadEnergies::BeamCalClusterList BCPadEnergies::lookForNeighbouringClustersOverWithVetoAndCheck(const BCPadEnergies &background,
												 const BCPadEnergies &backgroundSigma,
												 const BCPCuts &cuts) const {
  BCPadEnergies::BeamCalClusterList BeamCalClusters;

  //We make a copy, because we might want to apply different clustering on the same pads
  BCPadEnergies testPads(*this);
  testPads.subtractEnergiesWithCheck(background, backgroundSigma);
  //here cuts are applied on the pads
  PadIndexList myPadIndices = ( cuts.useConstPadCuts() ) ?
    getPadsAboveThresholds(testPads, cuts) :
    testPads.getPadsAboveSigma(backgroundSigma, cuts);

  ClusterNextToNearestNeighbourTowers(testPads, myPadIndices, cuts, BeamCalClusters);

  return BeamCalClusters;
} // lookForNeighbouringClustersOverWithVetoAndCheck



BCPadEnergies::BeamCalClusterList BCPadEnergies::lookForNeighbouringClustersOverSigma( const BCPadEnergies &backgroundSigma,
										       const BCPCuts &cuts,
										       bool detailedPrintout) const {
  BCPadEnergies::BeamCalClusterList BeamCalClusters;

  //here cuts are applied on the pads
  PadIndexList myPadIndices = ( cuts.useConstPadCuts() ) ?
    getPadsAboveThresholds(cuts) :
    getPadsAboveSigma(backgroundSigma, cuts);

  ClusterNextToNearestNeighbourTowers(*this, myPadIndices, cuts, BeamCalClusters, detailedPrintout);

  return BeamCalClusters;
} // lookForNeighbouringClustersOverWithVetoAndCheck


// Cluster towers as far as they satisfy the size cut + one neighborhood level outside
void ClusterNextToNearestNeighbourTowers( BCPadEnergies const& testPads, BCPadEnergies::PadIndexList& myPadIndices, const BCPCuts &cuts, BCPadEnergies::BeamCalClusterList& BeamCalClusters, bool DetailedPrintout)
{
  BeamCalClusters.clear();
  // Use only towers that pass cuts
  BCPadEnergies::TowerIndexList allTowersInBeamCal = BCPadEnergies::getContiguousTowersFromPads( testPads.m_BCG, myPadIndices);

  while( not allTowersInBeamCal.empty() )
  {
    //
    //Compare the largest deposit to the others and check if they are neighbours
    //
    BCPadEnergies::TowerIndexList::iterator largestTower =
    		std::max_element( allTowersInBeamCal.begin(), allTowersInBeamCal.end(), value_comparer);
    // If largest tower too small, we will not find anything else.
    if (largestTower->second < cuts.getMinimumTowerSize()) return;

    if (DetailedPrintout) {
      std::cout << "Largest Tower PadID " << largestTower->first << " : " << std::setw(3) << largestTower->second
		<< testPads.streamPad(largestTower->first)
		<< std::endl;
    }

    //We keep all the padIDs making up this tower in here, and we store the neighbors we have to check for additional neighbors
    BCPadEnergies::TowerIndexList towersInThisCluster, checkNextNeighborsList;
    towersInThisCluster.insert( *largestTower );

    do
    {
      if( not checkNextNeighborsList.empty() ) {
    	  largestTower = checkNextNeighborsList.begin();
      }

      //compare largest tower with all other towers
      for (BCPadEnergies::TowerIndexList::iterator it = allTowersInBeamCal.begin(); it != allTowersInBeamCal.end(); ++it)
      {
        // If the tower is our largest tower do nothing
		if( it->first == largestTower->first ) continue;
		// if the tower is already in the list we do nothing
		if ( towersInThisCluster.find(it->first) != towersInThisCluster.end() ) continue;

		const bool isNeighbour = testPads.m_BCG.arePadsNeighbours(largestTower->first, it->first);
		if ( isNeighbour ) {
		  if ( DetailedPrintout ) {
			std::cout << "Found a neighbor " << std::setw(6) << it->first << " : " << std::setw(3) << it->second
				  << testPads.streamPad(it->first)
				  << std::endl;
		  }//debug output

		  towersInThisCluster.insert( *it );
		  // Next level of neighbours is sought only if this tower is big enough
		  if (it->second >= cuts.getMinimumTowerSize()) checkNextNeighborsList.insert( *it );

		} //if they are neighbours
      }// find neighbouring towers/pads

      //remove the tower we just used from the checkNextNeighborsList
      checkNextNeighborsList.erase( largestTower->first );

    } while ( not checkNextNeighborsList.empty() );

    // Remove the towers from the towerIndexList for the next cluster round
    // and perform the cut on the cluster size.
    // The cluster size is calculated only from towers that pass the single size cut themselves
    int clusterSize=0;
    for (BCPadEnergies::TowerIndexList::iterator it = towersInThisCluster.begin(); it != towersInThisCluster.end(); ++it)
    {
      allTowersInBeamCal.erase(it->first);
      if(it->second >= cuts.getMinimumTowerSize()) clusterSize += it->second;
    }
    // Factor two makes sense for showers centered between neighboring pads
    if (clusterSize < 2*cuts.getMinimumTowerSize()) continue;


    //Create Cluster from the selected pads
    BCPadEnergies::PadIndexList padsForThisCluster(BCPadEnergies::getPadsFromTowers( testPads.m_BCG, myPadIndices, towersInThisCluster ));
    BeamCalClusters.push_back( testPads.getClusterFromAcceptedPads( testPads, padsForThisCluster, cuts) );
    BeamCalClusters.back().setPadIndexInLayer(max_element(towersInThisCluster.begin(), towersInThisCluster.end(), value_comparer)->first);

  }//while there are towers

}//ClusterNextToNearestNeighbourTowers



// Cluster contiguous towers in unlimited levels of neighborhood.
// Use only towers that pass the size cut
void ClusterContiguousTowers( BCPadEnergies const& testPads, BCPadEnergies::PadIndexList& myPadIndices, const BCPCuts &cuts, BCPadEnergies::BeamCalClusterList& BeamCalClusters, bool DetailedPrintout)
{
  BeamCalClusters.clear();
  // Use only towers that pass cuts
  BCPadEnergies::TowerIndexList allTowersInBeamCal = BCPadEnergies::getContiguousTowersFromPads( testPads.m_BCG, myPadIndices, cuts );

  while( not allTowersInBeamCal.empty() )
  {
    //
    //Compare the largest deposit to the others and check if they are neighbours
    //
    BCPadEnergies::TowerIndexList::iterator largestTower =
    		std::max_element( allTowersInBeamCal.begin(), allTowersInBeamCal.end(), value_comparer);

    if (DetailedPrintout) {
      std::cout << "Largest Tower PadID " << largestTower->first << " : " << std::setw(3) << largestTower->second
		<< testPads.streamPad(largestTower->first)
		<< std::endl;
    }

    //We keep all the padIDs making up this tower in here, and we store the neighbors we have to check for additional neighbors
    BCPadEnergies::TowerIndexList towersInThisCluster, checkNextNeighborsList;
    towersInThisCluster.insert( *largestTower );

    do
    {
      if( not checkNextNeighborsList.empty() ) {
    	  largestTower = checkNextNeighborsList.begin();
      }

      //compare largest tower with all other towers
      for (BCPadEnergies::TowerIndexList::iterator it = allTowersInBeamCal.begin(); it != allTowersInBeamCal.end(); ++it)
      {
        // If the tower is our largest tower do nothing
		if( it->first == largestTower->first ) continue;
		// if the tower is already in the list we do nothing
		if ( towersInThisCluster.find(it->first) != towersInThisCluster.end() ) continue;

		const bool isNeighbour = testPads.m_BCG.arePadsNeighbours(largestTower->first, it->first);
		if ( isNeighbour ) {
		  if ( DetailedPrintout ) {
			std::cout << "Found a neighbor " << std::setw(6) << it->first << " : " << std::setw(3) << it->second
				  << testPads.streamPad(it->first)
				  << std::endl;
		  }//debug output

		  towersInThisCluster.insert( *it );
		  checkNextNeighborsList.insert( *it );

		} //if they are neighbours
      }// find neighbouring towers/pads

      //remove the tower we just used from the checkNextNeighborsList
      checkNextNeighborsList.erase( largestTower->first );

    } while ( not checkNextNeighborsList.empty() );

    //Remove the towers from the towerIndexList for the next cluster round
    for (BCPadEnergies::TowerIndexList::iterator it = towersInThisCluster.begin(); it != towersInThisCluster.end(); ++it) {
      allTowersInBeamCal.erase(it->first);
    }

    //Create Cluster from the selected pads
    BCPadEnergies::PadIndexList padsForThisCluster(BCPadEnergies::getPadsFromTowers( testPads.m_BCG, myPadIndices, towersInThisCluster ));
    BeamCalClusters.push_back( testPads.getClusterFromAcceptedPads( testPads, padsForThisCluster, cuts) );
    BeamCalClusters.back().setPadIndexInLayer(max_element(towersInThisCluster.begin(), towersInThisCluster.end(), value_comparer)->first);

  }//while there are towers

}//ClusterNextToNearestNeighbourTowers




BCPadEnergies::BCPadEnergies::PadIndexList BCPadEnergies::getPadsAboveThreshold(double threshold) const{
  PadIndexList myPadIndices;
  for (int k = 0; k < m_BCG.getPadsPerBeamCal();++k) {
    double padEnergy(getEnergy(k));
	if(padEnergy > threshold) myPadIndices.push_back(k);
  }//all pads
  return myPadIndices;
}


BCPadEnergies::BCPadEnergies::PadIndexList BCPadEnergies::getPadsAboveThresholds(const BCPadEnergies& testPads, const BCPCuts& cuts) const{
  PadIndexList myPadIndices;
  for (int k = 0; k < testPads.m_BCG.getPadsPerBeamCal();++k) {
    if( testPads.m_BCG.getLayer(k) >= cuts.getStartingLayer() ) {
      double padEnergy(testPads.getEnergy(k));
      int padRing = testPads.m_BCG.getRing(k);
      if( cuts.isPadAboveThreshold(padRing, padEnergy) ) {
	myPadIndices.push_back(k);
      }
    }
  }//all pads
  return myPadIndices;
}


BCPadEnergies::BCPadEnergies::PadIndexList BCPadEnergies::getPadsAboveThresholds(const BCPCuts& cuts) const{
  PadIndexList myPadIndices;
  for (int k = 0; k < m_BCG.getPadsPerBeamCal();++k) {
    if( m_BCG.getLayer(k) >= cuts.getStartingLayer() ) {
      double padEnergy(getEnergy(k));
      int padRing = m_BCG.getRing(k);
      if( cuts.isPadAboveThreshold(padRing, padEnergy) ) {
	myPadIndices.push_back(k);
      }
    }
  }//all pads
  return myPadIndices;
}


BCPadEnergies::PadIndexList BCPadEnergies::getPadsAboveSigma(const BCPadEnergies& sigma,
							     const BCPCuts& cuts) const {
  PadIndexList myPadIndices;
  int nAbove=0;
  for (int k = 0; k < this->m_BCG.getPadsPerBeamCal();++k) {
    if( this->m_BCG.getLayer(k) >=  cuts.getStartingLayer() ) {
      double padEnergy(this->getEnergy(k));
      double padSigma(sigma.getEnergy(k));
      if( padEnergy > cuts.getPadSigmaCut() * padSigma ) {
	myPadIndices.push_back(k); nAbove++;
      }
    }
  }//all pads
  streamlog_out(DEBUG) << double(nAbove)/m_BCG.getPadsPerBeamCal()*100. << " % pads are above " << cuts.getPadSigmaCut() << " sigma.\n";
  return myPadIndices;
}


/// Gets a list of indices which make up the cluster myPadIndices
/// Sums up all the energy of this cluster, calculates the average position of the cluster
/// operates on *this
BeamCalCluster BCPadEnergies::getClusterFromPads(const PadIndexList& myPadIndices) const {
  BeamCalCluster BCCluster;
  double phi(0.0), ringAverage(0.0), totalEnergy(0.0);
  double thetaAverage(0.0);

  //Averaging an azimuthal angle is done via sine and cosine
  double yStore(0.0), xStore(0.0);

  //now take all the indices and add them to a cluster
  for (PadIndexList::const_iterator it = myPadIndices.begin(); it != myPadIndices.end(); ++it) {
    //Threshold was applied to get the padIndices
    const double energy(getEnergy(*it));
    const int ring = m_BCG.getRing(*it);
    const double thisPhi = m_BCG.getPadPhi(*it)* M_PI / 180.0; //Degrees to Radian
    BCCluster.addPad(*it, energy);
    //    phi+= thisPhi*energy;
    ringAverage += double( ring ) * energy;
    totalEnergy += energy;
    // Quasi-Cartesian coordinates in units of rad
    yStore += energy * sin( thisPhi ) * m_BCG.getThetaFromRing( ring ) ;
    xStore += energy * cos( thisPhi ) * m_BCG.getThetaFromRing( ring ) ;
//    thetaAverage += m_BCG.getThetaFromRing( ring ) * energy;

  }
  if(totalEnergy > 0.0) {
    phi =  atan2(yStore,  xStore) * 180.0 / M_PI ;
    if(phi < 0) phi += 360;
    //    phi /= totalEnergy;
    ringAverage /= totalEnergy;
    thetaAverage = sqrt(xStore*xStore + yStore*yStore) / totalEnergy;

    BCCluster.setPhi(phi);
    BCCluster.setRing(ringAverage);
    //    BCCluster.setTheta( m_BCG.getThetaFromRing( ringAverage ) * 1000 );
    BCCluster.setTheta( thetaAverage * 1000 );
  }


  return BCCluster;
}//getClusterFromPads



/// Gets a list of padEnergies testPads, and a list of indices which make up the cluster myPadIndices
/// Sums up all the energy of this cluster, calculates the average position of the cluster
/// Could be static except for m_BCG, should be m_BCG function
BeamCalCluster BCPadEnergies::getClusterFromAcceptedPads(const BCPadEnergies& testPads, const PadIndexList& myPadIndices, const BCPCuts& ) const {
  BeamCalCluster BCCluster;
  double phi(0.0), ringAverage(0.0), totalEnergy(0.0);
  double thetaAverage(0.0);

  //Averaging an azimuthal angle is done via sine and cosine
  double yStore(0.0), xStore(0.0);

  //now take all the indices and add them to a cluster
  for (PadIndexList::const_iterator it = myPadIndices.begin(); it != myPadIndices.end(); ++it) {
    //Threshold was applied to get the padIndices
    const double energy(testPads.getEnergy(*it));
    const int ring = m_BCG.getRing(*it);
    const double thisPhi = m_BCG.getPadPhi(*it)* M_PI / 180.0; //Degrees to Radian
    BCCluster.addPad(*it, energy);
    ringAverage += double( ring ) * energy;
    totalEnergy += energy;
    // Quasi-Cartesian coordinates in units of rad
    yStore += energy * sin( thisPhi ) * m_BCG.getThetaFromRing( ring+1 ) ;
    xStore += energy * cos( thisPhi ) * m_BCG.getThetaFromRing( ring+1 ) ;
  }
  if(totalEnergy > 0.0) {
    phi =  atan2(yStore,  xStore) * 180.0 / M_PI ;
    if(phi < 0) phi += 360;
    //    phi /= totalEnergy;
    ringAverage /= totalEnergy;
    thetaAverage = sqrt(xStore*xStore + yStore*yStore) / totalEnergy;

    BCCluster.setPhi(phi);
    BCCluster.setRing(ringAverage);
    BCCluster.setTheta( thetaAverage * 1000 );
  }

  return BCCluster;
}//getClusterFromAcceptedPads




BCPadEnergies::TowerIndexList BCPadEnergies::getTowersFromPads( BeamCalGeo const& geo, const PadIndexList& myPadIndices) {
  TowerIndexList myTowerIndices;
  for (PadIndexList::const_iterator it = myPadIndices.begin(); it != myPadIndices.end(); ++it) {
    myTowerIndices[ *it % geo.getPadsPerLayer() ] += 1;
  }
  return myTowerIndices;
}


BCPadEnergies::TowerIndexList BCPadEnergies::getTowersFromPads( BeamCalGeo const& geo, const PadIndexList& myPadIndices, const BCPCuts& cuts) {
  TowerIndexList myTowerIndices, finalTowerIndices;
  for (PadIndexList::const_iterator it = myPadIndices.begin(); it != myPadIndices.end(); ++it) {
    myTowerIndices[ *it % geo.getPadsPerLayer() ] += 1;
  }
  for (TowerIndexList::const_iterator it = myTowerIndices.begin(); it != myTowerIndices.end(); ++it)
  {
	  if(it->second >= cuts.getMinimumTowerSize()) finalTowerIndices.insert(*it);
  }
  return finalTowerIndices;
}



BCPadEnergies::TowerIndexList BCPadEnergies::getContiguousTowersFromPads( BeamCalGeo const& geo, const PadIndexList& myPadIndices)
{
  // Declare and initialize tower maps
  std::map<int, std::vector<bool> > towerMaps;
  for(int iTower=0; iTower<geo.getPadsPerLayer(); iTower++)
  {
	  std::vector<bool> *newTowerMap = new std::vector<bool>;
	  for(int iLayer=0; iLayer<geo.getBCLayers(); iLayer++)
		  newTowerMap->push_back(false);
	  towerMaps.insert(std::pair<int, std::vector<bool> >(iTower, *newTowerMap));
  }

  // Fill tower maps
  for (PadIndexList::const_iterator it = myPadIndices.begin(); it != myPadIndices.end(); ++it)
  {
	  int iTower = (*it) % geo.getPadsPerLayer();
	  int iLayer = (*it) / geo.getPadsPerLayer();
	  towerMaps[iTower].at(iLayer) = true;
  }

  // Find maximum uninterrupted series of hits in each tower
  // and fill the list
  TowerIndexList contiguousTowers;
  for(int iTower=0; iTower<geo.getPadsPerLayer(); iTower++)
  {
	  int maxSize=0, lastSize=0;
	  for(int iLayer=0; iLayer < geo.getBCLayers(); iLayer++)
	  {
		  if(towerMaps[iTower].at(iLayer)) lastSize++;
		  else
		  {
			  if(maxSize<lastSize) maxSize=lastSize;
			  lastSize=0;
		  }
	  }
	  contiguousTowers[iTower]=maxSize;
  }

  return contiguousTowers;
}

BCPadEnergies::TowerIndexList BCPadEnergies::getContiguousTowersFromPads( BeamCalGeo const& geo, const PadIndexList& myPadIndices, const BCPCuts& cuts) {

  // Declare and initialize tower maps
  std::map<int, std::vector<bool> > towerMaps;
  for(int iTower=0; iTower<geo.getPadsPerLayer(); iTower++)
  {
	  std::vector<bool> *newTowerMap = new std::vector<bool>;
	  for(int iLayer=0; iLayer<geo.getBCLayers(); iLayer++)
		  newTowerMap->push_back(false);
	  towerMaps.insert(std::pair<int, std::vector<bool> >(iTower, *newTowerMap));
  }

  // Fill tower maps
  for (PadIndexList::const_iterator it = myPadIndices.begin(); it != myPadIndices.end(); ++it)
  {
	  int iTower = (*it) % geo.getPadsPerLayer();
	  int iLayer = (*it) / geo.getPadsPerLayer();
	  towerMaps[iTower].at(iLayer) = true;
  }

  // Find maximum uninterrupted series of hits in each tower
  // and fill the list
  TowerIndexList contiguousTowers, finalTowers;
  for(int iTower=0; iTower<geo.getPadsPerLayer(); iTower++)
  {
	  int maxSize=0, lastSize=0;
	  for(int iLayer=0; iLayer < geo.getBCLayers(); iLayer++)
	  {
		  if(towerMaps[iTower].at(iLayer)) lastSize++;
		  else
		  {
			  if(maxSize<lastSize) maxSize=lastSize;
			  lastSize=0;
		  }
	  }
	  contiguousTowers[iTower]=maxSize;
  }

  // Keep only towers that pass the cuts
  for (TowerIndexList::const_iterator it = contiguousTowers.begin(); it != contiguousTowers.end(); ++it)
  {
	  if(it->second >= cuts.getMinimumTowerSize()) finalTowers.insert(*it);
  }
  return finalTowers;
}


void BCPadEnergies::removeTowerFromPads( BeamCalGeo const& geo, BCPadEnergies::PadIndexList& myPadIndices, int towerNumber) {
  BCPadEnergies::PadIndexList::iterator iter = myPadIndices.begin();
  while (  iter != myPadIndices.end() )  {
    if ( ( *iter % geo.getPadsPerLayer() ) == towerNumber ) {
      iter = myPadIndices.erase( iter );
    } else {
      ++iter;
    }
  }
  return;
}


void BCPadEnergies::removeTowersFromPads( BeamCalGeo const& geo,
					  BCPadEnergies::PadIndexList& myPadIndices,
					  BCPadEnergies::TowerIndexList towerNumbers) {

  BCPadEnergies::TowerIndexList::iterator towerIt = towerNumbers.begin();
  while ( towerIt != towerNumbers.end() ) {
    const int towerNumber = towerIt->first;
    BCPadEnergies::removeTowerFromPads( geo, myPadIndices, towerNumber);
    ++towerIt;
  }//all towerIndices
  return;
}




std::string BCPadEnergies::streamPad(int padID) const {
  std::stringstream out;
  int  layer, ring, pad;
  this->m_BCG.getLayerRingPad(padID, layer, ring, pad);
  out << "  Ring:" << std::setw(3) << ring
      << "  Pad:" << std::setw(4) << pad
      << "  Phi:" << std::setw(10) << this->m_BCG.getPadPhi(ring, pad);
  return out.str();
}


BCPadEnergies::PadIndexList BCPadEnergies::getPadsFromTowers ( BeamCalGeo const& geo,
							       BCPadEnergies::PadIndexList const&  allPads,
							       BCPadEnergies::TowerIndexList const& towers) {

  BCPadEnergies::PadIndexList padsInTowers;
  for (BCPadEnergies::PadIndexList::const_iterator it = allPads.begin(); it != allPads.end(); ++it) {
    if ( towers.find( (*it) % geo.getPadsPerLayer() ) != towers.end() ) {
      padsInTowers.push_back( *it );
    }
  }//for all pads

  return padsInTowers;

}



BCPadEnergies::TowerIndexList* BCPadEnergies::getTopAndNeighbourTowers(double threshold) const{

  TowerIndexList *outputlist = new TowerIndexList;

  //here cuts are applied on the pads
  PadIndexList myPadIndices = getPadsAboveThreshold(threshold);

  TowerIndexList myTowerIndices = getTowersFromPads( this->m_BCG, myPadIndices );

  //Compare the largest deposit to the others and check if they are neighbours
  TowerIndexList::iterator largestTower = std::max_element(myTowerIndices.begin(), myTowerIndices.end(),
							   value_comparer);
  outputlist->insert(std::pair<int, int>(largestTower->first, largestTower->second));

  for (TowerIndexList::iterator it=myTowerIndices.begin(); it!=myTowerIndices.end(); it++) {
	if (it==largestTower) continue;
	const bool isNeighbour = this->m_BCG.arePadsNeighbours(largestTower->first, it->first);
	if (isNeighbour) outputlist->insert(std::pair<int, int>(it->first, it->second));
  }

  return outputlist;
}



BCPadEnergies::TowerIndexList* BCPadEnergies::getTopAndNNNeighbourTowers(double threshold) const{

  TowerIndexList *outputlist = new TowerIndexList;

  TowerIndexList firstNeighbours = *(this->getTopAndNeighbourTowers(threshold));
  for(TowerIndexList::iterator it=firstNeighbours.begin(); it!=firstNeighbours.end(); it++){
	  outputlist->insert(std::pair<int, int>(it->first, it->second));
  }

  //here cuts are applied on the pads
  PadIndexList myPadIndices = getPadsAboveThreshold(threshold);

  TowerIndexList myTowerIndices = getTowersFromPads( this->m_BCG, myPadIndices );

  for (TowerIndexList::iterator it=myTowerIndices.begin(); it!=myTowerIndices.end(); it++) {
	bool isNeighbour = false;
	for(TowerIndexList::iterator jt=firstNeighbours.begin(); jt!=firstNeighbours.end(); jt++){
		if(m_BCG.arePadsNeighbours(jt->first, it->first))
			{ isNeighbour=true; break;}
	}
	// The insert method checks for repeated indices itself.
	if (isNeighbour) outputlist->insert(std::pair<int, int>(it->first, it->second));
  }

  return outputlist;
}


int BCPadEnergies::maxDepositTower() const{

  double maxDepo = 0.;
  int iMaxTower = 0;
  for(int iTower=1; iTower<=m_BCG.getPadsPerLayer(); iTower++)
  {
	  double depo=0.;
	  for(int iLayer=0; iLayer<m_BCG.getPadsPerLayer(); iLayer++)
		  depo += this->getEnergy(iTower+iLayer*m_BCG.getPadsPerLayer());

	  if(depo>maxDepo) { maxDepo=depo; iMaxTower = iTower;}
  }
  return iMaxTower;
}


BCPadEnergies::TowerIndexList* BCPadEnergies::getMaxTowerAndWithinRadius(double radius) const
{
  TowerIndexList *outputlist = new TowerIndexList;
//  TowerIndexList::iterator maxTower

  int iMaxTower = maxDepositTower();
  outputlist->insert(std::pair<int,int>(iMaxTower, m_BCG.getBCLayers()));

  for(int iTower=1; iTower<=m_BCG.getPadsPerLayer(); iTower++)
  {
	  if (m_BCG.getTransversalDistancePads(iMaxTower, iTower) < radius)
		  outputlist->insert(std::pair<int,int>(iTower, m_BCG.getBCLayers()));
  }

  return outputlist;
}


void BCPadEnergies::getGlobalCM(double &z, double &rho, double &phi)
{
	  z = phi = rho = 0.0;
	  double totalEnergy(0.0);

	  //Averaging is done via cartesian coordinates
	  double xStore(0.0), yStore(0.0), zStore(0.0);
	  double extents[6];

	  //now take all the indices and add them to a cluster
	  for (int iPad=1; iPad<m_BCG.getPadsPerBeamCal(); iPad++) {
		int layer, ring, pad;
		m_BCG.getLayerRingPad(iPad, layer, ring, pad);
	    m_BCG.getPadExtents(ring, pad, extents);
	    const double energy(getEnergy(iPad));
	    xStore += energy * extents[4]*cos( extents[5] * M_PI/180 );
	    yStore += energy * extents[4]*sin( extents[5] * M_PI/180 );
	    zStore += energy * double(layer);
	    totalEnergy += energy;
	  }

	  if(totalEnergy > 0.0) {
		double x = xStore/totalEnergy;
		double y = yStore/totalEnergy;
		z = zStore/totalEnergy;
	    phi =  atan2(y, x) * 180.0 / M_PI ;
	    rho = sqrt(x*x + y*y);
	  }

	  return;
}



BCPadEnergies::TowerIndexList* BCPadEnergies::getTowersWithinRadiusFromPoint(double rho, double phi, double radius) const
{
  TowerIndexList *outputlist = new TowerIndexList;
//  TowerIndexList::iterator maxTower

  for(int iTower=1; iTower<=m_BCG.getPadsPerLayer(); iTower++)
  {
	  if (m_BCG.getTransversalDistancePadToPoint(iTower, rho, phi) < radius)
		  outputlist->insert(std::pair<int,int>(iTower, m_BCG.getBCLayers()));
  }

  return outputlist;
}


// Return the pointer to a histogram containing the longitudinal profile of energy deposits in BeamCal
TH1D* BCPadEnergies::longitudinalProfile() const
{
  TH1D *profile = new TH1D("BClongProfile", "BeamCal longitudinal profile; layer; energy (a.u.)", m_BCG.getBCLayers(), 0.5, double(m_BCG.getBCLayers())+0.5);

  for (int i = 0; i < m_BCG.getPadsPerBeamCal();++i) {
	  profile->Fill(double(m_BCG.getLayer(i)), this->getEnergy(i));
  }

  return profile;
}

// Return the pointer to a histogram containing the longitudinal profile of energy deposits in BeamCal
// Restricted to pads in the padlist
TH1D* BCPadEnergies::longitudinalProfile(PadIndexList* padlist) const
{
  TH1D *profile = new TH1D("BClongProfile", "BeamCal longitudinal profile; layer; energy (a.u.)", m_BCG.getBCLayers(), 0.5, double(m_BCG.getBCLayers())+0.5);

  for (PadIndexList::iterator it=padlist->begin(); it!=padlist->end(); it++) {
	  profile->Fill(double(m_BCG.getLayer(*it)), this->getEnergy(*it));
  }

  return profile;
}


// Return the pointer to a histogram containing the longitudinal profile of energy deposits in BeamCal
// Restricted to pads in the towerlist
TH1D* BCPadEnergies::longitudinalProfile(TowerIndexList* towerlist) const
{
  TH1D *profile = dynamic_cast<TH1D*> (gROOT->FindObject("BClongProfile"));
  if(profile) delete profile;
  profile = new TH1D("BClongProfile", "BeamCal longitudinal profile; layer; energy (a.u.)", m_BCG.getBCLayers(), 0.5, double(m_BCG.getBCLayers())+0.5);

  for (TowerIndexList::iterator it=towerlist->begin(); it!=towerlist->end(); it++) {
	  for(int iLayer=0; iLayer < m_BCG.getBCLayers(); iLayer++){
		  // iLayer is actually the layer number minus 1
		  int iPad = it->first + iLayer*m_BCG.getPadsPerLayer();
		  profile->Fill(double(iLayer+1), this->getEnergy(iPad));
	  }
  }

  return profile;
}


// Return the pointer to a histogram containing the radial profile of energy deposits in BeamCal
// Restricted to pads in the towerlist
TH1D* BCPadEnergies::radialProfile(TowerIndexList* towerlist) const
{
  TH1D *profile = dynamic_cast<TH1D*> (gROOT->FindObject("BCradProfile"));
  if(profile) delete profile;
  profile = new TH1D("BCradProfile", "BeamCal radial profile; cylinder; energy (a.u.)", m_BCG.getBCRings(), 0.5, double(m_BCG.getBCRings())+0.5);

  for (TowerIndexList::iterator it=towerlist->begin(); it!=towerlist->end(); it++) {
	  for(int iLayer=0; iLayer < m_BCG.getBCLayers(); iLayer++){
		  // iLayer is actually the layer number minus 1
		  int iPad = it->first + iLayer*m_BCG.getPadsPerLayer();
		  int layer, ring, pad;
		  m_BCG.getLayerRingPad(iPad, layer, ring, pad);
		  profile->Fill(double(ring)+.5, this->getEnergy(iPad));
	  }
  }

  return profile;
}


// Return the pointer to a histogram containing the azimuthal profile of energy deposits in BeamCal
// Restricted to pads in the given ring
TH1D* BCPadEnergies::azimuthalProfile(int ring) const
{
  TH1D *profile = dynamic_cast<TH1D*> (gROOT->FindObject("BCphiProfile"));
  if(profile) delete profile;
  profile = new TH1D("BCphiProfile", "BeamCal azimuthal profile; pad; energy (a.u.)", m_BCG.getPadsInRing(ring), 0.5, double(m_BCG.getPadsInRing(ring))+0.5);

  for (int iPad = 1; iPad <= m_BCG.getPadsInRing(ring); iPad++) {
	  for(int iLayer=1; iLayer <= m_BCG.getBCLayers(); iLayer++){
		  // iLayer is actually the layer number minus 1
		  int iGlobalPad = m_BCG.getPadIndex(iLayer, ring, iPad);
		  profile->Fill(double(iPad)+.5, this->getEnergy(iGlobalPad));
	  }
  }

  return profile;
}



