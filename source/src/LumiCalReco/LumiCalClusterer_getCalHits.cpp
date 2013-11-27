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
void LumiCalClustererClass::getCalHits(	EVENT::LCEvent * evt,
					std::map < int, std::map < int, std::vector <IMPL::CalorimeterHitImpl*> > >& calHits) {


  try {
    EVENT::LCCollection	* col = evt->getCollection(_lumiName.c_str());

    int	nHitsCol;

    std::map < int , std::vector <IMPL::CalorimeterHitImpl*> > :: iterator		calHitsIterator;

#if _GENERAL_CLUSTERER_DEBUG == 1
    streamlog_out( DEBUG ) << std::endl  << "Getting hit information ...."
	      << std::endl << std::endl;
#endif

    nHitsCol = col->getNumberOfElements() ;
    if (not _mydecoder) {
      _mydecoder = new CellIDDecoder<SimCalorimeterHit> (col);
    }
    for (int i=0; i<nHitsCol; ++i) {

      // local-loop variables
      IMPL::SimCalorimeterHitImpl	* calHitIn;      // input CalorimeterHitImpl
      IMPL::CalorimeterHitImpl	* calHitNew;      // new   CalorimeterHitImpl

      int	cellId, arm, layer, rCell, phiCell;
      double	engyHit, xHit, yHit, zHit, rHit, phiHit;
      float	hitPosV[3];

      // get the hit from the LCCollection with index i
      calHitIn = static_cast<IMPL::SimCalorimeterHitImpl*> (col->getElementAt(i));


      // ???????? DECIDE/FIX - the global coordinates are going to change in new Mokka
      // versions, so the z must be extracted from the cellId instead ... ????????
      /// APS: Can be taken from the position of the calohit, when it is stored

      // get parameters from the input IMPL::CalorimeterHitImpl
      engyHit = (double)calHitIn -> getEnergy();

      layer   = (*_mydecoder)( calHitIn )["K"] ;
      phiCell = (*_mydecoder)( calHitIn )["J"] ;
      rCell   = (*_mydecoder)( calHitIn )["I"] ;

      ///APS Calculate the cellID with the function from here
      cellId = GlobalMethodsClass::CellIdZPR(layer, phiCell, rCell);
      // detector layer (ring) - count layers from zero and not from one
      layer -= 1 ;


      // skip this event if the following conditions are not met
      if(layer >= _maxLayerToAnalyse)	continue;
      if(engyHit < _hitMinEnergy)	continue;


      rHit	= (rCell+0.5)   * _rCellLength + _rMin;
      phiHit	= (phiCell+0.5) * _phiCellLength;

      zHit      = (double)calHitIn -> getPosition()[2];

      // compute x.y hit coordinates
      xHit = rHit * cos(phiHit);
      yHit = rHit * sin(phiHit);
      // write x,y,z to an array
      hitPosV[0] = xHit; hitPosV[1] = yHit; hitPosV[2] = zHit;

      // determine the side (arm) of the hit -> (+,-)1
      arm = (( zHit < 0 ) ? -1 : 1);

      // create a new IMPL::CalorimeterHitImpl
      calHitNew = new IMPL::CalorimeterHitImpl();

      // write the parameters to the new IMPL::CalorimeterHitImpl
      calHitNew -> setCellID0(cellId);
      calHitNew -> setEnergy(engyHit);
      calHitNew -> setPosition(hitPosV);

      // add the IMPL::CalorimeterHitImpl to a vector according to the detector
      // arm, and sum the total collected energy at either arm
      calHits[arm /* arm is variable */ ][layer].push_back( calHitNew );
      _totEngyArm[arm] += engyHit;
    }

  } // try


  // if an exception has been thrown (no *col for this event) then do....
  catch( EVENT::DataNotAvailableException &e){
#if _GENERAL_CLUSTERER_DEBUG == 1
    streamlog_out( ERROR ) << "Event has a SimCalorimeterHitImpl exception"<< std::endl;
#endif
  }



  return;
}

#endif // LumiCalClusterer_getCalHits_h
