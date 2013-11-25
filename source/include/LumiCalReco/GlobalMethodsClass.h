class GlobalMethodsClass {

public:

  GlobalMethodsClass();
  ~GlobalMethodsClass();

  void	SetConstants();

  map < TString , int >		GlobalParamI;
  map < TString , double >	GlobalParamD;
  map < TString , TString >	GlobalParamS;

  double	SignalGevConversion( TString optName , double valNow );

  double	Distance2DPolar( map <TString , double> pos1 , map <TString , double> pos2 );

  void	ThetaPhiCell(int cellId , map <TString , double> * thetaPhiCellP);

  int	CellIdZPR(int cellZ, int cellPhi, int cellR);
  int	CellIdZPR(int cellId, string ZPR);

};


GlobalMethodsClass :: GlobalMethodsClass(){
}

GlobalMethodsClass :: ~GlobalMethodsClass(){
}



/* --------------------------------------------------------------------------
   (1):	return a cellId for a given Z (layer), R (cylinder) and Phi (sector)
   (2):	return Z (layer), R (cylinder) and Phi (sector) for a given cellId
   -------------------------------------------------------------------------- */
int GlobalMethodsClass::CellIdZPR(int cellZ, int cellPhi, int cellR) {

  int cellId = 0;

  cellId |= ( cellZ   << 0  ) ;
  cellId |= ( cellPhi << 8  ) ;
  cellId |= ( cellR   << 16 ) ;

  return cellId;
}

int GlobalMethodsClass::CellIdZPR(int cellId, string ZPR) {

  int cellZ, cellPhi, cellR;
  int result = 0;

  // compute Z,Phi,R coordinets acording to the cellId
  cellZ   = (cellId >> 0 ) & 0xff ;
  cellPhi = (cellId >> 8 ) & 0xff ;
  cellR   = (cellId >> 16 ) & 0xff ;

  if(ZPR == "Z") result = cellZ;
  if(ZPR == "R") result = cellR;
  if(ZPR == "P") result = cellPhi;

  return result;
}



void GlobalMethodsClass::SetConstants() {

  // fiducial volume of LumiCal - minimal and maximal polar angles [rad]
  GlobalParamD["ThetaMin"] = 41e-3;
  GlobalParamD["ThetaMax"] = 69e-3;

  // geometry parameters of LumiCal

  // starting position [mm]
  GlobalParamD["ZStart"]   = 2270.;
  // inner and outer radii [mm]
  GlobalParamD["RMin"] = 80;
  GlobalParamD["RMax"] = 190;
  // cell division numbers
  GlobalParamI["NumCellsR"] = 64;
  GlobalParamI["NumCellsPhi"] = 48;
  GlobalParamI["NumCellsZ"] = 30;

  GlobalParamD["RCellLengt"] = (GlobalParamD["RMax"] - GlobalParamD["RMin"]) / double(GlobalParamI["NumCellsR"]);
  // layer thickness
  GlobalParamD["ZLayerThickness"] = 4.5;

  // logarithmic constant for position reconstruction
  GlobalParamD["LogWeightConstant"] = 6.;

  // Moliere radius of LumiCal [mm]
  GlobalParamD["MoliereRadius"] = 14;

  // minimal separation distance between any pair of clusters [mm]
  GlobalParamD["MinSeparationDist"] = GlobalParamD["MoliereRadius"];

  // minimal energy of a single cluster
  GlobalParamD["MinClusterEngy"] = 20;  // value in GeV
  GlobalParamD["MinClusterEngy"] = SignalGevConversion("GeV_to_Signal" , GlobalParamD["MinClusterEngy"]);  // conversion to "detector Signal"


}


/* --------------------------------------------------------------------------
   merge pairs of particles with separation/energy low cuts
   -------------------------------------------------------------------------- */
double GlobalMethodsClass::Distance2DPolar( map <TString , double> pos1 , map <TString , double> pos2 ) {

  return Sqrt(	pos1["r"]*pos1["r"] + pos2["r"]*pos2["r"]
		- 2 * pos1["r"]*pos2["r"] * Cos(pos1["phi"] - pos2["phi"]) );

}


/* --------------------------------------------------------------------------
   ccccccccc
   -------------------------------------------------------------------------- */
double GlobalMethodsClass::SignalGevConversion( TString optName , double valNow ){

  double	returnVal = -1;

  if(optName == "GeV_to_Signal")
    returnVal = valNow * .0105;
  //		returnVal = valNow * .0105 + .0013;

  if(optName == "Signal_to_GeV")
    returnVal = valNow / .0105;
  //		returnVal = (valNow-.0013) / .0105;

  return returnVal;
}



void GlobalMethodsClass::ThetaPhiCell(int cellId , map <TString , double> * thetaPhiCellP) {

  map <TString , double>  thetaPhiCellV = * thetaPhiCellP;

  int	cellIdR, cellIdZ, cellIdPhi;
  double	rCell, zCell, thetaCell, phiCell;

  // compute Z,Phi,R coordinets acording to the cellId
  cellIdZ   = (cellId >> 0 ) & 0xff ;
  cellIdPhi = (cellId >> 8 ) & 0xff ;
  cellIdR   = (cellId >> 16 ) & 0xff ;

  // theta
  rCell      = GlobalParamD["RMin"] + (cellIdR + .5) * GlobalParamD["RCellLengt"];
  zCell      = Abs(GlobalParamD["ZStart"]) + GlobalParamD["ZLayerThickness"] * (cellIdZ - 1);
  thetaCell  = ATan(rCell / zCell);

  // phi
  phiCell   = TwoPi() * (double(cellIdPhi) + .5) / double(GlobalParamI["NumCellsPhi"]);


  // fill output container
  thetaPhiCellV["theta"] = thetaCell;
  thetaPhiCellV["phi"]   = phiCell;
  thetaPhiCellV["r"]     = rCell;

  * thetaPhiCellP = thetaPhiCellV;

  return;
}
