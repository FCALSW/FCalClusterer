#include "BCPadEnergies.hh"
#include "BeamCalCluster.hh"
#include "BCPCuts.hh"
#include "BeamCalGeo.hh"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>


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


int BCPadEnergies::getTowerEnergies(int padIndex, std::vector<double> & te) const
{
  te.clear();
  if (padIndex < 0 || padIndex > m_BCG.getPadsPerLayer() )
  {
    std::cout << "BCPadEnergies::getTowerEnergies bad padIndex requested:" << padIndex
      << std::endl;
  }

  int pc(padIndex);
  while ( pc < m_BCG.getPadsPerBeamCal() ){
    te.push_back(m_PadEnergies.at(pc));
    pc+= m_BCG.getPadsPerLayer();
  }

  return te.size();
}

double BCPadEnergies::getTowerEnergy(int padIndex, int startLayer) const
{
  if (padIndex < 0 || padIndex > m_BCG.getPadsPerLayer() )
  {
    std::cout << "BCPadEnergies::getTowerEnergies bad padIndex requested:" << padIndex
      << std::endl;
  }

  int pc(padIndex + startLayer*m_BCG.getPadsPerLayer());
  double te(0.);
  while ( pc < m_BCG.getPadsPerBeamCal() ){
    te+= m_PadEnergies.at(pc);
    pc+= m_BCG.getPadsPerLayer();
  }

  return te;
}

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
      m_PadEnergies[i] -= 0.10 * sigma.m_PadEnergies[i] ;
    }

    if( (m_BCG.getLayer(i) == 10) &&
	(m_BCG.getRing(i) == 0) ) {
      if ( m_PadEnergies[i] > 0.9 * sigma.m_PadEnergies[i] && sigma.m_PadEnergies[i] > 1e-9 )  {
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

int BCPadEnergies::addEnergy(int layer, int ring, int pad, double energy) {
  const int padID = m_BCG.getPadIndex(layer, ring, pad);
  m_PadEnergies[padID] += energy;
  return padID;
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



BCPadEnergies::BeamCalClusterList BCPadEnergies::lookForNeighbouringClustersOverWithVetoAndCheck(const BCPadEnergies &background,
												 const BCPadEnergies &backgroundSigma,
												 const BCPCuts &cuts) const {
  BCPadEnergies::BeamCalClusterList BeamCalClusters;

  //We make a copy, because we might want to apply different clustering on the same pads
  BCPadEnergies testPads(*this);
//   testPads.subtractEnergies(background);
  testPads.subtractEnergiesWithCheck(background, backgroundSigma);
  //here cuts are applied on the pads
  PadIndexList myPadIndices = ( cuts.useConstPadCuts() ) ?
    getPadsAboveThresholds(testPads, cuts) :
    testPads.getPadsAboveSigma(backgroundSigma, cuts);

  testPads.clusterNextToNearestNeighbourTowers(myPadIndices, cuts, BeamCalClusters);

  return BeamCalClusters;
} // lookForNeighbouringClustersOverWithVetoAndCheck



BCPadEnergies::BeamCalClusterList BCPadEnergies::lookForNeighbouringClustersOverSigma( const BCPadEnergies &backgroundSigma,
										       const BCPCuts &cuts,
										       bool detailedPrintout) const {
  BCPadEnergies::BeamCalClusterList BeamCalClusters;

  //here cuts are applied on the pads
  PadIndexList myPadIndices = getPadsAboveSigma(backgroundSigma, cuts);

  this->clusterNextToNearestNeighbourTowers(myPadIndices, cuts, BeamCalClusters, detailedPrintout);

  return BeamCalClusters;
} // lookForNeighbouringClustersOverWithVetoAndCheck


void BCPadEnergies::clusterNextToNearestNeighbourTowers( const BCPadEnergies::PadIndexList &myPadIndices,
					  const BCPCuts &cuts,
					  BCPadEnergies::BeamCalClusterList &BeamCalClusters,
					  bool DetailedPrintout) const {

  BCPadEnergies::TowerIndexList allTowersInBeamCal = BCPadEnergies::getTowersFromPads( this->m_BCG, myPadIndices );
  while( not allTowersInBeamCal.empty() ) {

    //
    //Compare the largest deposit to the others and check if they are neighbours
    //

    BCPadEnergies::TowerIndexList::iterator largestTower = std::max_element( allTowersInBeamCal.begin(),
									     allTowersInBeamCal.end(),
									     value_comparer);

    if (DetailedPrintout) {
      std::cout << "Largest Tower PadID " << largestTower->first << " : " << std::setw(3) << largestTower->second
		<< this->streamPad(largestTower->first)
		<< std::endl;
    }

    //If the largest element in smaller than the minimal size we wont find anything else
    if( largestTower->second < cuts.getMinimumTowerSize() ) { return; }

    //We keep all the padIDs making up this tower in here, and we store the neighbors we have to check for additional neighbors
    BCPadEnergies::TowerIndexList towersInThisCluster, checkNextNeighborsList;
    towersInThisCluster.insert( *largestTower );

    do {

      if( not checkNextNeighborsList.empty() ) {
	largestTower = checkNextNeighborsList.begin();
      }

      //compare largest tower with all other towers
      for (BCPadEnergies::TowerIndexList::iterator it = allTowersInBeamCal.begin(); it != allTowersInBeamCal.end(); ++it) {
	if( it->first == largestTower->first ) continue;

	//if the tower is already in the list we do nothing
	if ( towersInThisCluster.find(it->first) != towersInThisCluster.end() ) continue;

	const bool isNeighbour = this->m_BCG.arePadsNeighbours(largestTower->first, it->first);
	if ( isNeighbour ) {
	  if ( DetailedPrintout ) {
	    std::cout << "Found a neighbor " << std::setw(6) << it->first << " : " << std::setw(3) << it->second
		      << this->streamPad(it->first)
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
    BCPadEnergies::PadIndexList padsForThisCluster(BCPadEnergies::getPadsFromTowers( this->m_BCG, myPadIndices, towersInThisCluster ));
    BeamCalClusters.push_back( this->getClusterFromAcceptedPads( *this, padsForThisCluster, cuts) );
    BeamCalClusters.back().setPadIndexInLayer(max_element(towersInThisCluster.begin(), towersInThisCluster.end(), value_comparer)->first);

  }//while there are towers

}//clusterNextToNearestNeighbourTowers

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


BCPadEnergies::PadIndexList BCPadEnergies::getPadsAboveSigma(const BCPadEnergies& sigma,
							     const BCPCuts& cuts) const {
  PadIndexList myPadIndices;
  for (int k = 0; k < this->m_BCG.getPadsPerBeamCal();++k) {
    if( this->m_BCG.getLayer(k) >=  cuts.getStartingLayer() ) {
      const double padEnergy(this->getEnergy(k));
      const double padSigma(sigma.getEnergy(k));
      const double cutValue = std::max(cuts.getPadSigmaCut() * padSigma, double(cuts.getMinPadEnergy()));
      if( padEnergy > cutValue ) {
	myPadIndices.push_back(k);
      }
    }
  }//all pads
  return myPadIndices;
}


/// Gets a list of padEnergies testPads, and a list of indices which make up the cluster myPadIndices
/// Sums up all the energy of this cluster, calculates the average position of the cluster
/// Could be static except for m_BCG, should be m_BCG function
BeamCalCluster BCPadEnergies::getClusterFromAcceptedPads(const BCPadEnergies& testPads, const PadIndexList& myPadIndices,
                                                         const BCPCuts& cuts) const {
  BeamCalCluster BCCluster;

  //Averaging an azimuthal angle is done via sine and cosine
  double sinStore(0.0), cosStore(0.0);
  double phi(0.0), ringAverage(0.0), thetaAverage(0.0), zAverage(0.0), totalEnergy(0.0), totalWeight(0.0);

  //get total cluster energy needed for log-weighting
  for (PadIndexList::const_iterator it = myPadIndices.begin(); it != myPadIndices.end(); ++it) {
    totalEnergy += testPads.getEnergy(*it);
  }

  const double logWeight       = cuts.getLogWeighting();
  const double firstRingRadius = m_BCG.getPadMiddleR(0, 0);

  //now take all the indices and add them to a cluster
  for (auto padID : myPadIndices) {
    //Threshold was applied to get the padIndices
    const double energy  = testPads.getEnergy(padID);
    const int    ring    = m_BCG.getRing(padID);
    const int    layer   = m_BCG.getLayer(padID);
    const double thisPhi = m_BCG.getPadPhi(padID) * M_PI / 180.0;  //Degrees to Radian
    BCCluster.addPad(padID, energy);
    const double weight =
        (logWeight < 0) ? energy : std::max(0.0,
                                            (logWeight + std::log(testPads.getEnergy(padID) / totalEnergy)) *
                                                firstRingRadius / m_BCG.getPadMiddleR(ring, 0));
    if (not(weight > 0.0)) {
      continue;
    }
    totalWeight += weight;
    ringAverage += double(ring) * weight;
    sinStore += weight * sin(thisPhi);
    cosStore += weight * cos(thisPhi);
    thetaAverage += m_BCG.getThetaFromRing(layer, ring) * weight;
    zAverage += m_BCG.getLayerZDistanceToIP(layer) * weight;
  }
  if (totalWeight > 0.0) {
    phi = atan2(sinStore / totalWeight, cosStore / totalWeight) * 180.0 / M_PI;
    if(phi < 0) phi += 360;
    ringAverage /= totalWeight;
    thetaAverage /= totalWeight;
    zAverage /= totalWeight;
    BCCluster.setPhi(phi);
    BCCluster.setRing(ringAverage);
    BCCluster.setTheta( thetaAverage * 1000 );
    BCCluster.setZ(zAverage);
  }

  // correct the reconstructed Phi for the "Right" side
  // beamcal is rotated and phi goes the other way in global coordinates
  if( m_side == kRight ) {
    phi = 360 - phi;
    while(phi < 0)   phi += 360;
    while(phi > 360) phi -= 360;
    BCCluster.setPhi(phi);
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
