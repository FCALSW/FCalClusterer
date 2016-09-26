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
  } catch ( EVENT::DataNotAvailableException &e ) {
    streamlog_out( ERROR ) << "Event has a SimCalorimeterHitImpl exception"<< std::endl;
    return 0;
  }

#if _GENERAL_CLUSTERER_DEBUG == 1
  streamlog_out( MESSAGE4 ) << std::endl  << "Getting hit information .... event: "<< evt->getEventNumber() << std::endl;
#endif

    if (not _mydecoder) {
      _mydecoder = new CellIDDecoder<SimCalorimeterHit> (col);
    }
    const int nHitsCol = col->getNumberOfElements();
    if ( nHitsCol < _clusterMinNumHits ) return 0;

    for (int i=0; i<nHitsCol; ++i) {
      
      int arm(0), layer(0), rCell(0), phiCell(0);

      // get the hit from the LCCollection with index i
      IMPL::SimCalorimeterHitImpl * calHitIn = static_cast<IMPL::SimCalorimeterHitImpl*> (col->getElementAt(i));

      const double engyHit = (double)calHitIn -> getEnergy();
     if(engyHit < _hitMinEnergy)	continue;

      // get cell parameters from the input IMPL::CalorimeterHitImpl

#ifdef FCAL_WITH_DD4HEP
     try{     
       arm     = (*_mydecoder)( calHitIn )["barrel"] - 1 ;   //        from 0
       layer   = (*_mydecoder)( calHitIn )["layer"] ;        // counts from 1
       rCell   = (*_mydecoder)( calHitIn )["r"] ;            //        from 0
       phiCell = (*_mydecoder)( calHitIn )["phi"] ;          //        from 0
     } catch ( Exception &e) {
       streamlog_out( ERROR )<<" Exception: " << e.what() << " in LumiCal_getCalHits with DD4HEP while encoding cellID: " << std::endl; 
     } 
#else
     try{
       arm     = (*_mydecoder)( calHitIn )["S-1"] ;          //        from 0
       layer   = (*_mydecoder)( calHitIn )["K"] ;            // counts from 1
       rCell   = (*_mydecoder)( calHitIn )["I"] ;            //        from 0
       phiCell = (*_mydecoder)( calHitIn )["J"] ;            //        from 0
     } catch ( Exception &e ) {
       streamlog_out( ERROR )<<" Exception: " << e.what() << " in LumiCal_getCalHits without DD4HEP while encoding cellID " << std::endl; 
     }
#endif

     //          


      ///APS Calculate the cellID with the function from here
       const int cellId = GlobalMethodsClass::CellIdZPR(layer, phiCell, rCell, arm);
       /*  
         const int cellId_test = calHitIn->getCellID0(); 
         if( cellId_test != cellId )   std::cout<<" mydecoder says S,I,J,K "<< arm <<", "<< rCell <<", "<< phiCell <<", "<< layer <<std::endl;
       */
      // detector layer  - count layers from zero and not from one
      layer -= 1 ;
     // determine the side (arm) of the hit -> (+,-)1
      arm = (( arm == 0 ) ? -1 : 1);

       // skip this hit if the following conditions are met
      if(layer >= _maxLayerToAnalyse || layer < 0 )	continue;

      
      /*(BP) it is not safe - in case non-zero crossing angle
            - phi sectors numbering order changes on -ve side
            - in some models there is layers relative phi offset  
            - also in local system zHit is +ve always
     
      const double rHit = (rCell+0.5) * _rCellLength + _rMin;
      //const double phiHit = (phiCell+0.5) * _phiCellLength;
      const double phiHit =  phiCell * _phiCellLength + double( (layer+1)%2 )*_zLayerPhiOffset;

      // compute x.y.z hit coordinates ( local )
      const float xHit = rHit * cos(phiHit);
      const float yHit = rHit * sin(phiHit);
      const float zHit = _zFirstLayer + float( layer ) * _zLayerThickness;

      // write x,y,z to an array
      float hitPosV[3] = {xHit, yHit, zHit};
      */
         
      // rotate to local Lcal reference system 
      const float* Pos = calHitIn->getPosition();     
      float xl =  Pos[0]*RotMat[arm]["cos"] - Pos[2]*RotMat[arm]["sin"];
      float yl =  Pos[1];
      float zl =  Pos[0]*RotMat[arm]["sin"] + Pos[2]*RotMat[arm]["cos"];

      // (BP) keep wrong sign of zHit local as it was  
      // in case it  matters ( still some guys might want to determine arm using sign of z )
       zl *= arm;
       float locPos[3] = { xl, yl, zl };

#if _GENERAL_CLUSTERER_DEBUG == 1
  	std::cout << "\t CellId, Pos(x,y,z), hit energy [MeV]: "
		  << std::setw(8)
		  << cellId << "\t ("
		  << locPos[0] << ", " 
		  << locPos[1] << ", " 
		  << locPos[2] << "), " 
		  << std::scientific << std::setprecision(3)
		  << 1000.*engyHit
		  <<std::endl;
#endif    
     // 
      // create a new IMPL::CalorimeterHitImpl
      IMPL::CalorimeterHitImpl *calHitNew = new IMPL::CalorimeterHitImpl();

      // write the parameters to the new IMPL::CalorimeterHitImpl
      calHitNew -> setCellID0(cellId);
      calHitNew -> setEnergy(engyHit);
      calHitNew -> setPosition(locPos);

      // add the IMPL::CalorimeterHitImpl to a vector according to the detector
      // arm, and sum the total collected energy at either arm
      calHits[arm][layer].push_back( calHitNew );
      _numHitsInArm[arm]++;
      _totEngyArm[arm] += engyHit;
    }//for all simHits

#if _GENERAL_CLUSTERER_DEBUG == 1
    streamlog_out( MESSAGE4 ) << std::endl  << "Energy deposit: "<< _totEngyArm[-1] << "\t" << _totEngyArm[1] <<"\n"
			   << "Number of hits: "<< _numHitsInArm[-1] << "\t" << _numHitsInArm[1] << "\n\n";
#endif

    if(    (( _numHitsInArm[-1] < _clusterMinNumHits) || (_totEngyArm[-1] < _minClusterEngySignal))
	   && (( _numHitsInArm[ 1] < _clusterMinNumHits) || (_totEngyArm[ 1] < _minClusterEngySignal)) ){
      calHits.clear();
      return 0;
    }else{ return 1; }
}

#endif // LumiCalClusterer_getCalHits_h
