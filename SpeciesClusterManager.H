#ifndef _SPECIES_CLUSTER_MNGR_
#define _SPECIES_CLUSTER_MNGR_
/*This class reads in a set of clusters for each of the species. Then it maps each gene
 * to its orthogroup. If there are multiple representations of a gene in an orthogroup, then
 * we consider only one of them. 
 * */
#include <map>
#include <string>
#include "CommonTypes.H"
using namespace std;
class Expert;
class GeneExpManager;
class MappedOrthogroupReader;
class GammaManager;
class GeneNameMapper;
typedef map<int,Expert*> CLUSTERSET;
typedef map<int,Expert*>::iterator CLUSTERSET_ITER;

class SpeciesClusterManager
{
	public:
	SpeciesClusterManager();
	~SpeciesClusterManager();
	int setOrthogroupReader(MappedOrthogroupReader*);
	int setGammaManager(GammaManager*);
	//SK: function to set the test set GammaManager object when running in CVMode
	int setGammaManager_Test(GammaManager*);
	int setSrcSpecies(const char*);
	int readSpeciesData(const char*);
	int readSpeciesData(const char*,const char*);
	int initExperts(bool);
	int getExtantClusterAssignment(map<int,map<string,int>*>&);
	int getExtantClusterAssignment(map<int,map<string,int>*>&,map<string,int>&);
	map<string,CLUSTERSET*>& getExtantSpeciesClusters();
	//Do an EM: Could do one step or could do until convergence
	int estimateExpertParameters(const char*);
	double getScore();
	int showClusters(const char*);
	int showClusters_Extant(const char*);
	int showClusters_Ancestral(const char*);
	int showMeans(const char* outputDir);
	int showMeans(const char* outputDir,int);
	int dumpAllInferredClusterAssignments(const char* outputDir);
	int dumpAllInferredClusterAssignments(const char* outputDir,int);
	int dumpAllInferredClusters_Scerwise(const char* outputDir,vector<string>& speciesList);
	int dumpAllInferredClusters_ScerwiseGrouped(const char* outputDir,vector<string>& speciesList);
	int dumpAllInferredClusters_LCA(const char* outputDir,vector<string>& speciesList, string& lcaName, const char* treefile);	// repop
	int dumpAllInferredClusterGammas(const char* outputDir,vector<string>& speciesList);
	//int dumpAllInferredClusters_Scerwise_Dup(const char* outputDir,vector<string>& speciesList);
	map<string,int>* getGenesForSpecies(string&);
	int setRandom(bool);
	int setRandSeed(int);
	int setMaxClusterCnt(int);
	int generateData(const char*,string&,vector<string>&);
        int setConstCov(double val); // FIXED COVARIANCE; Sets constant cov value (DC Add)
	//SK: set the number of iteractions to allow in the learning function if the algorithm doesn't converge
	int setNMaxIterations(int n);
	//SK: set the conversion threshold for the learning proceedure, the algorithm has convererged if the likelihood score changes by less than this amount
	int setConversionThreshold(double t);
	//SK: function to read in the lists of the test and training OGIDS sets for cross validation
	int setTrainingAndTestSet(map <int,int>& foldids,int f);
	//SK: function to set the boolean CVMode class variable, which is initialised to false in the constructor, and this function turns on the CVMode. 
	int setCVMode(bool);
	//SK: function for getting the score of the test data in cross validation mode.
	double getScore_Test();
	//SK: function for writting our score information
	int writeScores(const char*);
	//SK: function to set the secondStage boolean variable, see below.
	int setSecondStageOption(bool);
	int setPredictionInputMode(bool);
	int setMergedOGIDSet(map <int,int>& idset);
	int executePredictionMode(const char* outputDir, vector <string>, const char* treeFName);	// repop
	private:
	//SK: boolean class variable to change the mode of the objects functioning if cross validation is being done instead
	bool CVMode;
	//SK: map listing the OGIDS that belong to the training set. The first int in the map is the OGID for each entry
	map <int,int> trainingSet;
	//SK: map listing the OGIDS belonging to the test set, first int is the OGID value for each entry
	map <int,int> testSet;
	//SK: map object listing the orthogroups from which genes were read in with cluster assignments in readClusters.
	//SK: vector of likelihood scores for trainging data
	vector <double> scores;
	//SK: vector of likelihood scores for trainging data
        vector <double> testScores;
	map <int,int> workingOrthoGroups;
	//SK: the readClusters function was slightly updated for the cross validation step.
	int readClusters(string&, const char*);
	int maximizationStep();
	int expectationStep(bool estimate=true);
	int expectationStep_Species(string& specName, CLUSTERSET* expertSet);
	//SK: function to mimic the expectation step, but on the test OGIDS set for the cross validation proceedure
	int expectationStep_Test();
	//SK: function that is called for each species in the above function, and only for the cross validation step of scoring the test OGIDS (data)
        int expectationStep_Test_Species(string& specName, CLUSTERSET* expertSet);
	int estimateMeanCov(Expert* e, string& specName, int clusterID);
	int displaySpeciesClusters(const char* outFName,string& specName);
	int displaySpeciesClusterAssignments(const char* outFName,string& specName);
	int displaySpeciesClusterAssignments_NonScerNames(const char* outFName,string& specName);
	int assignGenesToExperts();
	int assignGenesToExperts_FromMap();
	map<string,CLUSTERSET*> speciesExpertSet;
	map<string,map<string,int>*> speciesClusterSet_Genewise;
	//Key is the ogid and the subsequent key is the name of the species.
	map<int,map<string,int>*> mappedClusterAssignment;
	map<int,map<string,int>*> allClusterAssignmentsGrouped;
	//First key is the species name. Second level key is the ogid
	map<int,map<string,map<int,double>*>*> mappedClusterGamma;
	map<string,map<int,int>*> speciesClusterSet_OGwise;
	map<string,int> testGenes;
	map<string,GeneExpManager*> speciesExprSet;
	MappedOrthogroupReader* mor;
	GammaManager* gammaMgr;
	//SK: additional gamma manager oject for the test data set. To be used only in cross validation mode
	GammaManager* gammaMgrTst;//SK: additional gamma manager oject for the test data set. TO be used only in cross validation mode
	//SK: set the number of clusters in the analysis
	int maxClusterCnt;
	bool randFlag;
	gsl_rng* r;
	GeneNameMapper gnm;
	char srcSpecies[256];
	int rseed;
	map<int,int> initRandAssign;
	//SK: maximum number of iterations to allow before ending the learning step regardless if the likelihood score has converged to the desired threshold level

        // FIXED COVARIANCE START //
        double constCov; // Fixed value for covariance (DC)
        bool fixCov; // If true, uses fixed value (DC)
        // FIXED COVARIANCE END //

	int nMaxIterations;
	//SK: threshold for the change in the likelihood score to meet the convergence condition.
	//SK: In order to converge the change in the score between iterations in the learning algorithm should be less than this value.
	double convThresh;
	//SK: a boolean to control if the second maximization phase should be done in SpeciesClusterManager::dumpAllInferredClusterAssignments() 
	bool secondStage;
	//SK: a boolean to control if the cluster means should be initialized from the source species. 
	bool predictionInputMode;
	map <int,int> mergedOGIDS;
	map <int,int> predictionOGIDS;
	map <int,int> initialOGIDS;
	int initPredictionOGIDS();
	int initPredictionGammas();
	int updateWorkingOGIDs();
	map <string,string> readSpeciesTree(const char*);	// repop
	map <string,string> children;	// repop
};
#endif
