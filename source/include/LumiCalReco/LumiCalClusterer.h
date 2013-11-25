class LumiCalClustererClass {

	
public:

	// Constructor
	LumiCalClustererClass( string lumiNameNow ) ;
	
	// initialization routine - Called at the begining of the job.
	void init(	map < TString , int >		GlobalParamI,
			map < TString , double >	GlobalParamD );


	// main actions in each event -Called for every event - the working horse. 
	void processEvent( EVENT::LCEvent * evt ) ;

	map < int , map < int , vector<int> > > 	_superClusterIdToCellId;
	map < int , map < int , vector<double> > > 	_superClusterIdToCellEngy;
	

protected:

	// Processor Parameters
	string 	_lumiName;
	int    	_clusterMinNumHits;
	double 	_hitMinEnergy;

	
	// global variables
	int	_numEventsPerTree, _resetRootTrees;
	int	_maxLayerToAnalyse;
	double	_zFirstLayer, _zLayerThickness, _rMin, _rMax, _rCellLength, _phiCellLength;
	double	_elementsPercentInShowerPeakLayer;
	double	_logWeightConst;
	int	_nNearNeighbor ;
	int	_cellRMax, _cellPhiMax ;
	double	_middleEnergyHitBoundFrac;
	string	_methodCM;
	double	_moliereRadius;
	double	_thetaContainmentBouds[2];
	double	_minSeparationDistance, _minClusterEngyGeV;

	map < int , double >	_totEngyArm;
	vector < int >		_armsToCluster;


	// methods:
	void	getCalHits(	EVENT::LCEvent 							* evt,
				map < int , map < int , vector <CalorimeterHitImpl*> > >	* calHitsP );


	int	buildClusters(	map < int , vector <CalorimeterHitImpl*> >	&calHits,
				map < int , CalorimeterHitImpl* > 		* calHitsCellIdGlobalP,
				map < int , vector<int> > 			* superClusterIdToCellIdP,
				map < int , vector<double> >  			* superClusterIdToCellEngyP,
				map < int , vector<double> > 			* superClusterCMP );
		
	int	initialClusterBuild( map < int , CalorimeterHitImpl* > 	calHitsCellId,
				     map < int , int > 			* cellIdToClusterIdP,
				     map < int , vector<int> > 		* clusterIdToCellIdP,
				     map < int , vector<double> > 	* clusterCMP,
				     vector < int > 			controlVar );

	int	initialLowEngyClusterBuild( map < int , CalorimeterHitImpl* > 	calHitsSmallEngyCellId,
					    map < int , CalorimeterHitImpl* > 	* calHitsCellIdP,
					    map < int , int > 			* cellIdToClusterIdP,
					    map < int , vector<int> > 		* clusterIdToCellIdP,
					    map < int , vector<double> > 	* clusterCMP );


	int	virtualCMClusterBuild( map < int , CalorimeterHitImpl* > 	calHitsCellId,
				       map < int , int > 			* cellIdToClusterIdP,
				       map < int , vector<int> > 		* clusterIdToCellIdP,
				       map < int , vector<double> > 		* clusterCMP,
				       map < int , vector<double> > 		virtualClusterCM );
	
	int	virtualCMPeakLayersFix(	map < int , CalorimeterHitImpl* > 	calHitsCellId,
					map < int , int > 			* cellIdToClusterIdP,
					map < int , vector<int> > 		* clusterIdToCellIdP,
					map < int , vector<double> > 		* clusterCMP,
					map < int , vector<double> > 		virtualClusterCM );

	int	buildSuperClusters ( map <int , CalorimeterHitImpl* > 			* calHitsCellIdGlobalP,
				     vector < map < int , CalorimeterHitImpl* > > 	calHitsCellId,
				     vector < map < int , vector<int> > > 		clusterIdToCellId,
				     vector < map < int , vector<double> > > 		clusterCM,
				     vector < map < int , vector<double> > > 		virtualClusterCM, 
				     map < int , int > 					* cellIdToSuperClusterIdP,
				     map < int , vector<int> > 				* superClusterIdToCellIdP,
				     map < int , vector<double> > 			* superClusterCMP );

	int	engyInMoliereCorrections ( map <int , CalorimeterHitImpl* > 		calHitsCellIdGlobal,
					   map < int , vector <CalorimeterHitImpl*> > 	calHits,
					   vector < map < int , CalorimeterHitImpl* > > calHitsCellIdLayer,
					   vector < map < int , vector<int> > > 	* clusterIdToCellIdP,
					   vector < map < int , vector<double> > > 	* clusterCMP,
					   vector < map < int , int > > 		* cellIdToClusterIdP,
					   map < int , int > 				* cellIdToSuperClusterIdP,
					   map < int , vector<int> > 			* superClusterIdToCellIdP,
					   map < int , vector<double> > 		* superClusterCMP,
					   double 					middleEnergyHitBound,
					   int 						detectorArm );



	void	energyCorrections (	map < int , vector<int> > 		* superClusterIdToCellIdP,
					map < int , vector<double> > 		* superClusterIdToCellEngyP,
					map < int , vector<double> >		* superClusterCMP,
					map < int , CalorimeterHitImpl* >	calHitsCellIdGlobal ) ;


	void	clusterMerger (	      map < int , vector<double> > 	* clusterIdToCellEngyP,
				      map < int , vector<int> > 	* clusterIdToCellIdP,
				      map < int , vector<double> >	* clusterCMP,
				      map < int , CalorimeterHitImpl* >	calHitsCellIdGlobal ) ;


	void	fiducialVolumeCuts (	map < int , vector<int> > 		* superClusterIdToCellIdP,
					map < int , vector<double> > 		* superClusterIdToCellEngyP,
					map < int , vector<double> >		* superClusterCMP ) ;


	void	getThetaPhiZCluster( map < int , CalorimeterHitImpl* >	calHitsCellId,
				     vector <int> 			clusterIdToCellId,
				     double 				totEngy,
				     double 				* output );
		
	int	getNeighborId( int	cellId,
			       int	neighborIndex );
	
	int	idToZPR( int 	cellId,
			 string	ZPR );
	
	int	ZPRToId( int	cellZ,
			 int	cellPhi,
			 int	cellR );
	
	double	posWeight( CalorimeterHitImpl	* calHit ,
			   string 		method );

	double	posWeightTureCluster( CalorimeterHitImpl	* calHit,
				      double 			cellEngy,
				      string 			method );

	double	posWeight( CalorimeterHitImpl	* calHit,
			   double 		totEngy,
			   string 		method );

	double	posWeight( CalorimeterHitImpl	* calHit,
			   double 		totEngy,
			   string 		method,
			   double 		logWeightConstNow );

	double	distance2D( double * pos1,
			    double * pos2 );

	double	distance2DPolar( double * pos1,
				 double * pos2 );

	double	thetaPhiCell( int	cellId,
			      string	output );
	
	void	getEngyPosCMValues( vector <int>			cellIdV,
				    map < int , CalorimeterHitImpl* > 	calHitsCellId,
				    double 				* engyPosCM,
				    string 				method );

	void	calculateEngyPosCM( vector <int> 			cellIdV,
				    map < int , CalorimeterHitImpl* > 	calHitsCellId,
				    map < int , vector<double> > 	* clusterCMp,
				    int 				clusterId,
				    string 				method );

	void	calculateEngyPosCM_EngyV(      vector <int> 				cellIdV,
					       vector <double> 				cellEngyV,
					       map < int , CalorimeterHitImpl* >	calHitsCellId,
					       map < int , vector<double> > 		* clusterCMp,
					       int 					clusterId,
					       string 					method );

	void	updateEngyPosCM( CalorimeterHitImpl	* calHit,
				 vector<double> 	* clusterCMp );

	int	checkClusterMergeCM( int 				clusterId1,
				     int 				clusterId2,
				     map < int , vector<int> > 		clusterIdToCellId, 
				     map < int , CalorimeterHitImpl* > 	calHitsCellId,
				     double 				distanceAroundCM,
				     double 				percentOfEngyAroungCM,
				     string 				method );

	double	getDistanceAroundCMWithEnergyPercent( vector < double >			clusterCM,
						      vector < int > 			clusterIdToCellId,
						      map < int , CalorimeterHitImpl* >	calHitsCellId,
						      double 				engyPercentage );

	double	getMoliereRadius( map < int , CalorimeterHitImpl* >	calHitsCellId,
				  vector <int>				clusterIdToCellId, 
				  vector <double>			clusterCM );

	double	getEngyInMoliereFraction( map < int , CalorimeterHitImpl* >	calHitsCellId, 
					  vector <int> 				clusterIdToCellId, 
					  vector <double> 			clusterCM,
					  double 				moliereFraction );

	double	getEngyInMoliereFraction( map < int , CalorimeterHitImpl* >	calHitsCellId, 
					  vector < int >			clusterIdToCellId, 
					  vector < double >			clusterCM,
					  double				moliereFraction,
					  map < int , int >			* flagP );

	double	engySignalGeV( double 	engy,
			       TString 	transformMethod );

	

};

