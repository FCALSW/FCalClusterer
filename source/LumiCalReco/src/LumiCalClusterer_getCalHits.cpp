#ifndef LumiCalClusterer_getCalHits_h
#define LumiCalClusterer_getCalHits_h 1
// Local
#include "LumiCalClusterer.h"
//LCIO
#include <Exceptions.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/SimCalorimeterHitImpl.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <UTIL/CellIDDecoder.h>
// Stdlib
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>

/* --------------------------------------------------------------------------
   Loop over al hits in the LCCollection and write the hits into vectors
   of CalorimeterHitImpl. Hits are split in two vectors, one for each arm
   of LumiCal.
   -------------------------------------------------------------------------- */
int LumiCalClustererClass::getCalHits(	EVENT::LCEvent * evt,
					MapIntMapIntVCalHit & calHits) {

  EVENT::LCCollection * col = NULL;
  try {
    col = evt->getCollection(_lumiName.c_str());
  } // try
  // if an exception has been thrown (no *col for this event) then do....
  catch( EVENT::DataNotAvailableException &e){
#if _GENERAL_CLUSTERER_DEBUG == 1
    streamlog_out( ERROR ) << "Event has a SimCalorimeterHitImpl exception"<< std::endl;
#endif
    return 0;
  }

#if _GENERAL_CLUSTERER_DEBUG == 1
    streamlog_out( DEBUG ) << std::endl  << "Getting hit information ...."
	      << std::endl << std::endl;
#endif

    if (not _mydecoder) {
      _mydecoder = new CellIDDecoder<SimCalorimeterHit> (col);
    }
    const int nHitsCol = col->getNumberOfElements();
    for (int i=0; i<nHitsCol; ++i) {
      
      int arm(0), layer, rCell, phiCell;
      //int side(0);

      // get the hit from the LCCollection with index i
      IMPL::SimCalorimeterHitImpl * calHitIn = static_cast<IMPL::SimCalorimeterHitImpl*> (col->getElementAt(i));

      // ???????? DECIDE/FIX - the global coordinates are going to change in new Mokka
      // versions, so the z must be extracted from the cellId instead ... ????????
      /// APS: Can be taken from the position of the calohit, when it is stored

      // get parameters from the input IMPL::CalorimeterHitImpl
      const double engyHit = (double)calHitIn -> getEnergy();
      //      const int cellid = (int)calHitIn->getCellID0();
      layer   = (*_mydecoder)( calHitIn )["K"] ;    // counts from 1
      phiCell = (*_mydecoder)( calHitIn )["J"] ;    //        from 0
      rCell   = (*_mydecoder)( calHitIn )["I"] ;    //        from 0
      //      side    = (*_mydecoder)( calHitIn )["S-1"] ;  //        from 0

      ///APS Calculate the cellID with the function from here
      const int cellId = GlobalMethodsClass::CellIdZPR(layer, phiCell, rCell, arm);
      // detector layer (ring) - count layers from zero and not from one
      layer -= 1 ;


      // skip this hit if the following conditions are met
      if(layer >= _maxLayerToAnalyse)	continue;
      if(engyHit < _hitMinEnergy)	continue;


      const double rHit = (rCell+0.5) * _rCellLength + _rMin;
      // BP: add relative layer offset
      //const double phiHit = (phiCell+0.5) * _phiCellLength;
      const double phiHit = (phiCell+0.5) * _phiCellLength + double( layer%2 )*_zLayerPhiOffset;

      const float zHit = (double)calHitIn -> getPosition()[2];
      // determine the side (arm) of the hit -> (+,-)1
      arm = (( zHit < 0 ) ? -1 : 1);

      // compute x.y hit coordinates ( local )
      const float xHit = rHit * cos(phiHit);
      const float yHit = rHit * sin(phiHit);
      // write x,y,z to an array
      float hitPosV[3] = {xHit, yHit, zHit};

      // create a new IMPL::CalorimeterHitImpl
      IMPL::CalorimeterHitImpl *calHitNew = new IMPL::CalorimeterHitImpl();

      // write the parameters to the new IMPL::CalorimeterHitImpl
      calHitNew -> setCellID0(cellId);
      calHitNew -> setEnergy(engyHit);
      calHitNew -> setPosition(hitPosV);

      // add the IMPL::CalorimeterHitImpl to a vector according to the detector
      // arm, and sum the total collected energy at either arm
      calHits[arm /* arm is variable */ ][layer].push_back( calHitNew );
      _numHitsInArm[arm]++;
      _totEngyArm[arm] += engyHit;
    }//for all simHits

    if(    (( _numHitsInArm[-1] < _clusterMinNumHits) || (_totEngyArm[-1] < _minClusterEngySignal))
	   && (( _numHitsInArm[ 1] < _clusterMinNumHits) || (_totEngyArm[ 1] < _minClusterEngySignal)) ){
      calHits.clear();
      return 0;
    }else{ return 1; }
}

#endif // LumiCalClusterer_getCalHits_h
