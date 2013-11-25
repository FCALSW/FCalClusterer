/* --------------------------------------------------------------------------
   class ....
   -------------------------------------------------------------------------- */
class ClusterClass : public GlobalMethodsClass {

public:
  ClusterClass(int idNow);
  ~ClusterClass();

  void	FillHit(int cellNow, double engyNow);
  int	ResetStats();
  void	SetStatsMC(MCParticle * mcParticle);
  void	SetStatsMC();

  int	Id, Pdg, SignMC, ParentId, NumMCDaughters;
  int	OutsideFlag, MergedFlag, HighestEnergyFlag, ModifiedFlag;
  double	Engy, Theta, Phi, RZStart;
  double	VtxX, VtxY, VtxZ, EndPointX, EndPointY, EndPointZ;
  double	EngyMC, ThetaMC, PhiMC;

  vector <int> MergedV;

  TString	OutsideReason;

  map < int , double >	Hit;

};


/* --------------------------------------------------------------------------
   class ....
   -------------------------------------------------------------------------- */
ClusterClass::ClusterClass(int idNow) {


  OutsideFlag		= 0;
  MergedFlag		= 0;
  HighestEnergyFlag	= 0;
  ModifiedFlag		= 0;

  OutsideReason = "";

  Engy = 0.;

  Id = ParentId = idNow;

  Hit.clear();
  MergedV.clear();

  // inherited method from BeamPipe class, to set all global constants
  SetConstants();
}

ClusterClass::~ClusterClass() {

}

void ClusterClass::SetStatsMC(MCParticle * mcParticle) {


  double	MCmomX, MCmomY, MCmomZ, MCr;
  MCParticle	* mcParticleParent;

  Pdg = (int)mcParticle->getPDG();
  Id  = (int)mcParticle -> id();

  VtxX    = (double)mcParticle->getVertex()[0];
  VtxY    = (double)mcParticle->getVertex()[1];
  VtxZ    = (double)mcParticle->getVertex()[2];

  EndPointX    = (double)mcParticle->getEndpoint()[0];
  EndPointY    = (double)mcParticle->getEndpoint()[1];
  EndPointZ    = (double)mcParticle->getEndpoint()[2];

  MCmomX    = (double)mcParticle->getMomentum()[0];
  MCmomY    = (double)mcParticle->getMomentum()[1];
  MCmomZ    = (double)mcParticle->getMomentum()[2];

  SignMC = int(MCmomZ / Abs(MCmomZ));
  MCr     = Sqrt(Power(MCmomX,2) + Power(MCmomY,2));
  ThetaMC = ATan(MCr / Abs(MCmomZ));

  // correction for the polar angle for cases where the particle does not
  // come from the IP (according to a projection on the face of LumiCal)
  if( Abs(VtxZ) > _VERY_SMALL_NUMBER ) {
    //				cout <<"   ---   "<< Id << "\t" << ThetaMC << "\t -> \t" ;
    MCr = Tan(ThetaMC) * (GlobalParamD["ZStart"] - Abs(VtxZ));
    ThetaMC = ATan( MCr / GlobalParamD["ZStart"] );
    //				cout	<< ThetaMC << "\t r,z:  " << MCr <<"\t" << Abs(vtxZ) <<  "\t"
    //					<<  mcParticle->getDaughters().size() << endl ;
    //				mcParticleParent = mcParticle -> getDaughters()[0];
    //				cout	<< "\t\t\t"<< mcParticleParent->getVertex()[2] <<  "\t"<< mcParticleParent->id() <<  "\t"<< endl ;
  }

  PhiMC = ATan(Abs(MCmomY / MCmomX));
  if(MCmomZ > 0) {
    if(MCmomX>0 && MCmomY>0) PhiMC  =		PhiMC;
    if(MCmomX<0 && MCmomY>0) PhiMC  = Pi() -	PhiMC;
    if(MCmomX<0 && MCmomY<0) PhiMC  = Pi() +	PhiMC;
    if(MCmomX>0 && MCmomY<0) PhiMC  = TwoPi() -		PhiMC;
  } else {
    if(MCmomX>0 && MCmomY>0) PhiMC  = Pi() -	PhiMC;
    if(MCmomX<0 && MCmomY>0) PhiMC  =		PhiMC;
    if(MCmomX<0 && MCmomY<0) PhiMC  = TwoPi() -		PhiMC;
    if(MCmomX>0 && MCmomY<0) PhiMC  = Pi() +	PhiMC;
  }

  NumMCDaughters = (int)mcParticle -> getDaughters().size();
  EngyMC        = (double)mcParticle->getEnergy();


  // find the original parent of the particle
  mcParticleParent = mcParticle;
  ParentId = Id;
  while(1) {
    if( !(mcParticleParent -> isCreatedInSimulation()) )
      break;

    mcParticleParent = mcParticleParent -> getParents()[0];
    ParentId = (int)mcParticleParent -> id();
  }


  return;
}



void ClusterClass::SetStatsMC() {

  ParentId = -1;
  Pdg = Id = SignMC = NumMCDaughters = 0;
  VtxX = VtxY = VtxZ = 0.;
  EndPointX = EndPointY = EndPointZ = 0.;
  ThetaMC = PhiMC = EngyMC = 0.;


  return;
}



void ClusterClass::FillHit(int cellNow, double engyNow) {

  Hit[cellNow] += engyNow;
  Engy         += engyNow;

  return;
}


int ClusterClass::ResetStats() {

  int	cellId, numHits;
  double	weightNow, engyNow;
  double	thetaSum, phiSum, weightSum, engySum;

  map < int , double > :: iterator	hitIterator;
  map < TString , double >		thetaPhiCellV;

  // initialize modification flag
  ModifiedFlag = 0;

  // initialization of summation variables
  engySum = weightSum = thetaSum = phiSum = 0.;

  // re-sum energy from all cells (just in case)
  numHits     = Hit.size();
  hitIterator = Hit.begin();
  for(int cellNow = 0; cellNow < numHits; cellNow++, hitIterator++) {
    cellId = (int)(*hitIterator).first;

    engySum += Hit[cellId];
  }
  Engy = engySum;

  // if the particle has no deposits in LumiCal
  if(Engy < _VERY_SMALL_NUMBER) {
    OutsideFlag = 1;
    OutsideReason = "No energy deposits at all";

    return 0;
  }

  numHits     = Hit.size();
  hitIterator = Hit.begin();
  for(int cellNow = 0; cellNow < numHits; cellNow++, hitIterator++) {
    cellId = (int)(*hitIterator).first;

    ThetaPhiCell(cellId , &thetaPhiCellV);

    engyNow = Hit[cellId];
    weightNow = GlobalParamD["LogWeightConstant"] + Log(engyNow/engySum);
    if(weightNow < 0) continue;

    thetaSum  += thetaPhiCellV["theta"] * weightNow;
    phiSum    += thetaPhiCellV["phi"]   * weightNow;

    weightSum += weightNow;
  }

  if(weightSum < 1e-9) {
    OutsideFlag = 1;
    OutsideReason = "No energy deposits above the minimal threshold";

    return 0;
  }


  Theta    = thetaSum / weightSum;
  RZStart  = ATan(Abs(Theta)) * GlobalParamD["ZStart"];
  Phi      = phiSum   / weightSum;

#if _CLUSTER_RESET_STATS_DEBUG == 1
  cout	<< coutBlue << "\tReset stats:   " << Id << endl
	<< "\t\tsign, engyHits, engyMC:  " << SignMC << "\t"
	<< SignalGevConversion("Signal_to_GeV" , Engy) << "\t" << EngyMC << endl
	<< "\t\ttheta mc,rec , phi mc,rec:   " << ThetaMC << "\t" << Theta << "\t" << PhiMC << "\t" << Phi << endl
	<< "\t\tvertex X,Y,Z:                " << VtxX << "\t" << VtxY << "\t" << VtxZ << endl
	<< coutDefault << endl;
#endif

  return 1;
}
