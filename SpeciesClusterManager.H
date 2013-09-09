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
	int setSrcSpecies(const char*);
	int readSpeciesData(const char*);
	int initExperts();
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
	int dumpAllInferredClusters_LCA(const char* outputDir,vector<string>& speciesList, string& lcaName);
	int dumpAllInferredClusterGammas(const char* outputDir,vector<string>& speciesList);
	//int dumpAllInferredClusters_Scerwise_Dup(const char* outputDir,vector<string>& speciesList);
	map<string,int>* getGenesForSpecies(string&);
	int setRandom(bool);
	int setRandSeed(int);
	int setMaxClusterCnt(int);
	int generateData(const char*,string&,vector<string>&);
	private:
	
	int readClusters(string&, const char*);
	int maximizationStep();
	int expectationStep();
	int expectationStep_Species(string& specName, CLUSTERSET* expertSet);
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
	//First key is the species name. Second level key is the ogid
	map<int,map<string,map<int,double>*>*> mappedClusterGamma;
	map<string,map<int,int>*> speciesClusterSet_OGwise;
	map<string,int> testGenes;
	map<string,GeneExpManager*> speciesExprSet;
	MappedOrthogroupReader* mor;
	GammaManager* gammaMgr;
	int maxClusterCnt;
	bool randFlag;
	gsl_rng* r;
	GeneNameMapper gnm;
	char srcSpecies[256];
	int rseed;
	map<int,int> initRandAssign;
};

#endif
