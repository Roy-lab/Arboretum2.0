#ifndef _FRAMEWOaRK_
#define _FRAMEWORK_
#include "gsl/gsl_randist.h"
/*The way we will do this is that we will assume that we have some mapping of the genes across
 * species. That is cluster k in species 1 will correspond to cluster k in another species.*/
class SpeciesDistManager;
class SpeciesClusterManager;
class GeneTreeManager;
class MappedOrthogroupReader;
class Framework
{
	public:
	Framework();
	~Framework();
	int readSpeciesData(const char*,const char*);
	int readSpeciesTree(int,const char*);
	int readOrthology(const char* specOrder, const char* orthomapfile);
	int setSrcSpecies(const char*);
	int startClustering(const char*);
	int generateData(const char*);
	int redisplay(const char*);
	int setPdiagonalLeaf(double);
	int setPdiagonalNonLeaf(double);
	int setClusterTransProb(double);
	int setClusterTransProb(const char*);
	private:
	int learnAncestralClusters();
	int initClusterTransitionProb();
	int initClusterTransitionProb(SpeciesDistManager::Species* anode);
	int initTransitionProb(Matrix* m,double);
	//This is all part of the expectation step
	int inferAncestralClusters(map<int,map<string,int>*>&);
	int inferExtantClusters(map<int,map<string,int>*>&);
	int estimateClusterTransProb(map<int,map<string,int>*>& clusterAssignments);
	int estimateClusterTransProb(SpeciesDistManager::Species* parent, SpeciesDistManager::Species* child, map<int,map<string,int>*>& clusterAssignments);
	int estimateTransitionMatrix(string& parentname,string& childname, SpeciesDistManager::Species* child,map<int,map<string,int>*>& clusterAssignments);
	//Maximization step. Not sure if we should estimate the cluster transition probabilites here or not.
	int maximizationStep();

	bool checkConvergence(map<int,map<string,int>*>& clusterAssignments,map<int,map<string,int>*>& oldclusterAssignments);
	int saveClusterAssignments(map<int,map<string,int>*>& clusterAssignments,map<int,map<string,int>*> &oldclusterAssignments);

	int sampleAncestralCluster(gsl_rng* r,SpeciesDistManager::Species* root);
	int sampleChildCluster(gsl_rng* r,SpeciesDistManager::Species*,int,map<string,int>& clusterAssign);
	int sampleChildCluster(gsl_rng* r, int parentClusterId,Matrix* params,vector<int>*);
	int sortIndices(SpeciesDistManager::Species*);
	int sortIndices(Matrix* values,int rowID, vector<int>*);
	int dumpInferredAssignments(const char* outputDir, const char* suffix);

	//This is very much like the SpeciesDist. Only it works with general potentials rather than pairwise potentials.
	SpeciesDistManager sdMgr;
	//This is the wrapper over the species-specific clusterings
	SpeciesClusterManager scMgr;
	GammaManager gammaMgr;
	//For now assume in each species we store the ML estimate of the cluster ID
	map<string,map<string,int>*> clusterAssignments;
	MappedOrthogroupReader mor;
	int maxClusterCnt;
	int learnIter;
	char outputDir[1024];
	gsl_rng* randnum;
	char srcSpecies[256];
	map<string,double> clusterTransitionProb;
	double p_diagonal_leaf;
	double p_diagonal_nonleaf;
};
#endif
