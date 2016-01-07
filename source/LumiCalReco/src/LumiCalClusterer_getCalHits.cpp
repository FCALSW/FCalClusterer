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
  // (BP) LumiCal rotation ( local -> global )
  const double CSBX[2] = { cos( M_PI - _beamCrossingAngle ),
                           cos( _beamCrossingAngle )        };
  const double SNBX[2] = { sin( M_PI - _beamCrossingAngle ),
                           sin( _beamCrossingAngle  )       };
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
      //      const int cellid = (int)calHitIn->getCellID0();
      layer   = (*_mydecoder)( calHitIn )["K"] ;    // counts from 1
      phiCell = (*_mydecoder)( calHitIn )["J"] ;    //        from 0
      rCell   = (*_mydecoder)( calHitIn )["I"] ;    //        from 0
      arm     = (*_mydecoder)( calHitIn )["S-1"] ;  //        from 0

      const double engyHit = (double)calHitIn -> getEnergy();

      ///APS Calculate the cellID with the function from here
      const int cellId = GlobalMethodsClass::CellIdZPR(layer, phiCell, rCell, arm);
      // detector layer (ring) - count layers from zero and not from one
      layer -= 1 ;

       // skip this hit if the following conditions are met
      if(layer >= _maxLayerToAnalyse || layer < 0 )	continue;
      if(engyHit < _hitMinEnergy)	continue;

      
      /*(BP) it is not safe - in case non-zero crossing angle
            - phi sectors numbering order changes on -ve side
            - in some models there is layers relative phi offset  
             also in local system zHit is +ve always
      const double rHit = (rCell+0.5) * _rCellLength + _rMin;
      // BP: add relative layer offset
      //const double phiHit = (phiCell+0.5) * _phiCellLength;
      const double phiHit = (phiCell+0.5) * _phiCellLength + double( layer%2 )*_zLayerPhiOffset;

      // compute x.y hit coordinates ( local )
      const float xHit = rHit * cos(phiHit);
      const float yHit = rHit * sin(phiHit);
     
      // write x,y,z to an array
      float hitPosV[3] = {xHit, yHit, zHit};
      */
      const float xHit = (double)calHitIn -> getPosition()[0];
      const float yHit = (double)calHitIn -> getPosition()[1];
      const float zHit = (double)calHitIn -> getPosition()[2];



      // write x,y,z to an array
      // rotate to local Lcal
 	float hitPosV[3];
	hitPosV[0] =  xHit*CSBX[arm] - zHit*SNBX[arm];
	hitPosV[1] =  yHit;
	hitPosV[2] =  xHit*SNBX[arm] + zHit*CSBX[arm];

     // determine the side (arm) of the hit -> (+,-)1
      arm = (( arm == 0 ) ? -1 : 1);
      // (BP) keep wrong sign of zHit local as it was  
      // in case it  matters
      //     hitPosV[2] *= arm;
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
