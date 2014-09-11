/*
 * profileExtractor.hh
 *
 * Classes derived from virtul class ProfileExtractor perform extraction of the
 * longitudinal energy profile from a BCPadEnergies object.
 * Different derived classes extract the profile in different ways, from the whole
 * BeamCal, or restricted to towers or pads specified in TowerIndexList or PadIndexList
 * objects.
 * The purpose is to allow one to choose the desired extraction method at runtime.
 *
 *  Created on: Sep 11, 2014
 *      Author: strahinja
 */

#ifndef PROFILEEXTRACTOR_HH_
#define PROFILEEXTRACTOR_HH_


class ProfileExtractor
{
public:
	ProfileExtractor() {};
	virtual ~ProfileExtractor() {};
	virtual TH1D* extractProfile(BCPadEnergies *bcpads) const = 0;
};

class ProfileAll : public ProfileExtractor
{
public:
	ProfileAll() {};
	~ProfileAll() {};

	TH1D* extractProfile(BCPadEnergies *bcpads) const {return bcpads->longitudinalProfile();}
};

class ProfileTopAndNNNTowers : public ProfileExtractor
{
public:
	ProfileTopAndNNNTowers(double threshold = 0.001) : m_threshold(threshold) {};
	~ProfileTopAndNNNTowers() {};
	TH1D* extractProfile(BCPadEnergies *bcpads) const
	{
	  BCPadEnergies::TowerIndexList *towerlist = bcpads->getTopAndNNNeighbourTowers(m_threshold);
	  TH1D *profile = bcpads->longitudinalProfile(towerlist);
	  delete towerlist;
	  return profile;
	}
protected:
	double m_threshold;
};


class ProfileInRadiusFromCM : public ProfileExtractor
{
public:
	ProfileInRadiusFromCM(double radius = 1.) : m_radius(radius){};
	~ProfileInRadiusFromCM() {};
	TH1D* extractProfile(BCPadEnergies *bcpads) const
	{
	  double z, rho, phi;
	  bcpads->getGlobalCM(z, rho, phi);
	  BCPadEnergies::TowerIndexList *towerlist = bcpads->getTowersWithinRadiusFromPoint(rho, phi, m_radius);
	  TH1D *profile = bcpads->longitudinalProfile(towerlist);
	  delete towerlist;
	  return profile;
	}
protected:
	double m_radius;
};



#endif /* PROFILEEXTRACTOR_HH_ */
