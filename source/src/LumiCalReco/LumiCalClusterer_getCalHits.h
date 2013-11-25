
/* --------------------------------------------------------------------------
   Loop over al hits in the LCCollection and write the hits into vectors
   of CalorimeterHitImpl. Hits are split in two vectors, one for each arm
   of LumiCal.
   -------------------------------------------------------------------------- */
void LumiCalClustererClass::getCalHits(	LCEvent								* evt,
					map < int , map < int , vector <CalorimeterHitImpl*> > >	* calHitsP ) {


  try{
    LCCollection	* col = evt->getCollection(_lumiName.c_str());

    int	nHitsCol;

    map < int , map < int , vector <CalorimeterHitImpl*> > >	calHits = * calHitsP;
    map < int , vector <CalorimeterHitImpl*> > :: iterator		calHitsIterator;

#if _GENERAL_CLUSTERER_DEBUG == 1
    cout	<< endl << coutBlue << "Getting hit information ...."
		<< coutDefault << endl << endl;
#endif

    if(col) nHitsCol = col->getNumberOfElements() ;
    if(col) for (int i=0; i<nHitsCol; ++i) {

	// local-loop variables
	SimCalorimeterHitImpl	* calHitIn;      // input CalorimeterHitImpl
	CalorimeterHitImpl	* calHitNew;      // new   CalorimeterHitImpl

	int	cellId, arm, layer, rCell, phiCell;
	double	engyHit, xHit, yHit, zHit, rHit, phiHit;
	float	hitPosV[3];

	// get the hit from the LCCollection with index i
	calHitIn = dynamic_cast<SimCalorimeterHitImpl*> (col->getElementAt(i));

	// create a new CalorimeterHitImpl
	calHitNew = new CalorimeterHitImpl();

	// ???????? DECIDE/FIX - the global coordinets are going to change in new Mokka
	// versions, so the z must be extracted from the cellId instead ... ????????

	// get parameters from the input CalorimeterHitImpl
	cellId  = (int)calHitIn    -> getCellID0();
	engyHit = (double)calHitIn -> getEnergy();

	layer   = ( cellId >> 0  ) & 0xff ;
	phiCell = ( cellId >> 8  ) & 0xff ;
	rCell	= ( cellId >> 16 ) & 0xff ;

	// detector layer (ring) - count laters from zero and not from one
	layer-- ;

	rHit	= (rCell+0.5)   * _rCellLength + _rMin;
	phiHit	= (phiCell+0.5) * _phiCellLength;

	zHit    = (double)calHitIn -> getPosition()[2];

	// compute x.y hit coordinets
	xHit = rHit * Cos(phiHit);
	yHit = rHit * Sin(phiHit);
	// write x,y,z to an array
	hitPosV[0] = xHit; hitPosV[1] = yHit; hitPosV[2] = zHit;

	// determine the side (arm) of the hit -> (+,-)1
	arm = int(zHit/Abs(zHit));

	// skip this event if the following conditions are not met
	if(layer >= _maxLayerToAnalyse)	continue;
	if(engyHit < _hitMinEnergy)	continue;

	// write the parameters to the new CalorimeterHitImpl
	calHitNew -> setCellID0(cellId);
	calHitNew -> setEnergy(engyHit);
	calHitNew -> setPosition(hitPosV);

	// add the CalorimeterHitImpl to a vector according to the detector
	// arm, and sum the total collected energy at either arm
	if(arm == -1)	calHits[arm][layer].push_back( calHitNew );
	else		calHits[arm][layer].push_back( calHitNew );

	_totEngyArm[arm] += engyHit;
      }


    * calHitsP = calHits;

  } // try


  // if an exception has been thrown (no *col for this event) than do....
  catch( DataNotAvailableException &e){
#if _GENERAL_CLUSTERER_DEBUG == 1
    cout << "Event has a SimCalorimeterHitImpl exception"<< endl;
#endif
  }



  return;
}
