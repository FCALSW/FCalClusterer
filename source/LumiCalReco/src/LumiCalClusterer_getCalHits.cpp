// Local
#include "Global.hh"
#include "GlobalMethodsClass.h"
#include "LumiCalClusterer.h"
#include "LumiCalHit.hh"

#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>
//LCIO
#include <EVENT/CalorimeterHit.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCIO.h>
#include <EVENT/LCObject.h>
#include <EVENT/LCParameters.h>
#include <EVENT/SimCalorimeterHit.h>
#include <Exceptions.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/BitField64.h>
#include <UTIL/CellIDDecoder.h>
// Stdlib
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

// IWYU pragma: no_include <bits/shared_ptr.h>

/* --------------------------------------------------------------------------
   Loop over al hits in the LCCollection and write the hits into vectors
   of CalorimeterHitImpl. Hits are split in two vectors, one for each arm
   of LumiCal.
   -------------------------------------------------------------------------- */
int LumiCalClustererClass::getCalHits(	EVENT::LCEvent * evt,
					MapIntMapIntVCalHit & calHits) {
  EVENT::LCCollection* col = nullptr;
  try {
    col = evt->getCollection(_lumiName.c_str());
  } catch (EVENT::DataNotAvailableException& e) {
    streamlog_out(WARNING) << "Event does not have the '" << _lumiName << "' collection" << std::endl;
    return 0;
  }

  streamlog_out(DEBUG6) << std::endl << "Getting hit information .... event: " << evt->getEventNumber() << std::endl;


    const int nHitsCol = col->getNumberOfElements();
    if ( nHitsCol < _clusterMinNumHits ) return 0;

    // figure out if we have a CalorimeterHit or SimCalorimeterHit collection,
    // and create CalorimeterHit collection if necessary
    if (dynamic_cast<EVENT::SimCalorimeterHit*>(col->getElementAt(0)) != nullptr) {
      col = createCaloHitCollection(col);
      evt->addCollection(col, _lumiOutName.c_str());
    }

    if (not _mydecoder) {
      _mydecoder =
        std::unique_ptr<UTIL::CellIDDecoder<EVENT::CalorimeterHit>>(new UTIL::CellIDDecoder<EVENT::CalorimeterHit>(col));
    }

    for (int i=0; i<nHitsCol; ++i) {
      
      int arm(0), layer(0);
      int rCell(0), phiCell(0);

      auto* calHitIn = static_cast<EVENT::CalorimeterHit*>(col->getElementAt(i));

      const double engyHit = (double)calHitIn -> getEnergy();

      if(engyHit < _hitMinEnergy)	continue;

      // ???????? DECIDE/FIX - the global coordinates are going to change in new Mokka
      // versions, so the z must be extracted from the cellId instead ... ????????
      /// APS: Can be taken from the position of the calohit, when it is stored

      // get parameters from the input IMPL::CalorimeterHitImpl
     

      //using Mokka simulated files
      if( not _useDD4hep ) {
	arm     = (*_mydecoder)( calHitIn )["S-1"]; // from 0
	rCell   = (*_mydecoder)( calHitIn )["I"]; // from 0
	phiCell = (*_mydecoder)( calHitIn )["J"]; // from 0
	layer   = (*_mydecoder)( calHitIn )["K"]; // counts from 1
	// detector layer  - count layers from zero and not from one
	layer -= 1 ;
	// determine the side (arm) of the hit -> (+,-)1
	arm = (( arm == 0 ) ? -1 : 1);

	if( arm < 0 ) {
	  //for rotation around the Y-axis, so that the phiCell increases counter-clockwise for z<0
	  phiCell = _cellPhiMax - phiCell;
	  phiCell += _cellPhiMax/2;
	  if(phiCell >= _cellPhiMax) phiCell -= _cellPhiMax;
	}

      } else {
	arm =(*_mydecoder)( calHitIn )["barrel"]; // from 1 and 2
	if( arm == 2 ) arm = -1;

	phiCell = (*_mydecoder)( calHitIn )["phi"]; // goes from -phiMax/2 to +phiMax/2-1
	if( arm < 0 ) {
	  //for rotation around the Y-axis, so that the phiCell increases counter-clockwise for z<0
          phiCell *= -1;
          //This is not needed any more because we no currently do not calculate pad positions from IDs any longer
          // now it is just important that we have the correct phiID range
	  //LumiCal is (or not, if we fix it) rotated by pi around Z for negative side
	  //phiCell += int(_gmc._backwardRotationPhi/(2.0*M_PI)*_cellPhiMax+0.5);
	}
        //limit to range 0 to _cellPhiMax-1
	if(phiCell >= _cellPhiMax) phiCell -= _cellPhiMax;
	if(phiCell < 0 ) phiCell += _cellPhiMax;

	rCell = (*_mydecoder)( calHitIn )["r"];
	layer = (*_mydecoder)( calHitIn )["layer"]; // counts from 0
      }

      //Calculate internal cellID
      int cellId = GlobalMethodsClass::CellIdZPR(layer, phiCell, rCell, arm);

      // skip this hit if the following conditions are met
      if(layer >= _maxLayerToAnalyse || layer < 0 )	continue;

      const float* Pos = calHitIn->getPosition();
      double       locPos[3] = {0.0, 0.0, 0.0};
      _gmc.rotateToLumiCal(Pos, locPos);

      streamlog_message(DEBUG2,
                        std::stringstream p;
                        p << std::scientific << std::setprecision(3)
                        << "\t Arm, CellId, Pos(x,y,z), hit energy [MeV]: "
                        << std::setw(5) << arm
                        << std::setw(13)
                        << cellId << "\t ("
                        << std::setw(13) << locPos[0] << ", "
                        << std::setw(13) << locPos[1] << ", "
                        << std::setw(13) << locPos[2] << "), "
                        << 1000.*engyHit
                        << "Layer/Phi/R/Arm"
                        << std::setw(5) << layer
                        << std::setw(5) << phiCell
                        << std::setw(5) << rCell
                        << std::setw(5) << arm
                        <<std::endl; ,
                        p.str(););

        // create a new LumiCalHit
        auto calHitNew = std::make_shared<LumiCalHit>();

        // write the parameters to the new LumiCalHit
        calHitNew->setCellID0(cellId);
        calHitNew->setEnergy(engyHit);
        calHitNew->setPosition(locPos);
        calHitNew->addHit(calHitIn);

        // add the LumiCalHit to a vector according to the detector
        // arm, and sum the total collected energy at either arm
        calHits[arm][layer].push_back(std::move(calHitNew));
        _numHitsInArm[arm]++;
        _totEngyArm[arm] += engyHit;
    }//for all simHits

    streamlog_out(DEBUG4) << "Energy deposit: "<< _totEngyArm[-1] << "\t" << _totEngyArm[1] <<"\n"
                          << "Number of hits: "<< _numHitsInArm[-1] << "\t" << _numHitsInArm[1] << std::endl;

    if (((_numHitsInArm[-1] < _clusterMinNumHits) || (_totEngyArm[-1] < _minClusterEngyGeV)) &&
        ((_numHitsInArm[1] < _clusterMinNumHits) || (_totEngyArm[1] < _minClusterEngyGeV))) {
      return 0;
    }else{ return 1; }
}

EVENT::LCCollection* LumiCalClustererClass::createCaloHitCollection(EVENT::LCCollection* simCaloHitCollection) const {
  streamlog_out(DEBUG7) << "Creating the CalorimeterHit collection with dummy digitization" << std::endl;

  const double calibrationFactor = _gmc.getCalibrationFactor();

  auto caloHitCollection = new IMPL::LCCollectionVec(EVENT::LCIO::CALORIMETERHIT);
  caloHitCollection->parameters().setValue(EVENT::LCIO::CellIDEncoding,
                                           simCaloHitCollection->getParameters().getStringVal(EVENT::LCIO::CellIDEncoding));
  caloHitCollection->setFlag(simCaloHitCollection->getFlag());
  for (int i = 0; i < simCaloHitCollection->getNumberOfElements(); ++i) {
    auto simHit = static_cast<EVENT::SimCalorimeterHit*>(simCaloHitCollection->getElementAt(i));
    auto calHit = new IMPL::CalorimeterHitImpl();

    calHit->setCellID0(simHit->getCellID0());
    calHit->setCellID1(simHit->getCellID1());
    calHit->setEnergy(simHit->getEnergy() * calibrationFactor);
    calHit->setPosition(simHit->getPosition());

    caloHitCollection->addElement(calHit);
  }

  return caloHitCollection;
}
