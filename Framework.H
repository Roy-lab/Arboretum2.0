#ifndef _FRAMEWOaRK_
#define _FRAMEWORK_

#include <vector>
#include "Graph.H"
#include "Node.H"
#include "GeneManager.H"
#include "ProteinManager.H"
#include "InteractionManager.H"
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "MotifManager.H"
#include "BFGSWrapperData.H"
#include "BFGSWrapper.H"
#include "TestData.H"
#include "MotifRegressor.H"



using namespace std;

#define CONT 1
#define DISC 2

#include "gsl/gsl_randist.h"
/*The way we will do this is that we will assume that we have some mapping of the genes across
 * species. That is cluster k in species 1 will correspond to cluster k in another species.*/

class Error; //SK: for learnMoE component
class VariableManager;//SK: for learnMoE component
class EvidenceManager;//SK: for learnMoE component
class MotifManager;//SK: for learnMoE component
class MotifRegressor;//SK: for learnMoE compoment
class TestData; //SK: for learnMoE component
class InteractionManager; //SK: added for the fgconverter component.
class ProteinManager; //SK: added for the fgconverter component
class Graph; //SK: added for the fgconverter component.
class GeneManager; //SK: added for the fgconverter component. 
class SpeciesDistManager;
class SpeciesClusterManager;
class GeneTreeManager;
class MappedOrthogroupReader;
class Framework
{
	public:
		Framework();
		~Framework();
		int readSpeciesData(const char*,const char*,const char*);
		int readSpeciesTree(int,const char*);
		int readOrthology(const char* specOrder, const char* orthomapfile);
		int setSrcSpecies(const char*);
		int startClustering(const char*);
		int startCrossValidationCheck(const char*, int cvF,const char*,const char*,const char*, const char*);//function for runnign the cross validation study
		int generateData(const char*);
		int fgconverter(const char*,const char*,int,int);
		int learnMoE(const char*,const char*,int);
		int genSpeciesClusters(const char*,const char*);
		int genSpeciesClustersNonSrc(const char*,const char*);
		int redisplay(const char*);
		int setPdiagonalLeaf(double);
		int setPdiagonalNonLeaf(double);
		int setClusterTransProb(double);
		int setClusterTransProb(const char*);
		int setNMaxIterations(int);
		int setConversionThreshold(double);
		//SK: function to set the option to run the second optimization stage in SpeciesClusterManager::dumpAllInferredClusterAssignments(const char* outputDir) or not.
		int setSecondStageOption(bool);
		int mergeData();
		int readExpressionData(const char* aFName,const char* bFName);
		int readExpressionDataNonSrc(const char* aFName,const char* bFName);
		int setPreClustering(bool);
		int countOGIDSNonSrc(const char* aFName);
		int setInputMode(bool);
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
		//SK: class variable for the number of maximum iterations to allow
		int nMaxIterations;
		//SK: class variable for the convergence threhsold
		//SK: the likelihood score must not change by more than this value between iterations for the learning algorithm to converge
		double convThresh;
		//SK: gsl_rng object to be used for random number generation to partition the OGIDS sets into folds for the cross validation proceedure
		gsl_rng* r_og;
		//SK: A boolean variable to run the second stage of analysis in SpeciesClusterManager::dumpAllInferredClusterAssignments(const char* outputDir) or not, initialized to true.  
		bool secondStage;
		bool preClustering;
		//SK: data merging and write out components
		char outFName[256];
		char pnlGraphFName[256];
		int vType;
		Graph phyNw;
		vector<double> mult_dist;
                InteractionManager pdMgr;
                InteractionManager ppMgr;
                ProteinManager protMgr;
                GeneManager geneMgr;
                int expCntPerGene;
                vector<int> expIds;
		int createGraph(map<int,int>&, map<string,int>&,map<string,int>&);
		int writeToFile(vector<Node*>&);
		int writeToFile_MeanStd(vector<Node*>& topOrder);
		//SK: learnMoE components
		VariableManager varManager;
                EvidenceManager evManager;
                MotifManager motifManager;
                TestData testData;
                double epsThreshold;
                char actionStr[256];
                double penalty;
                MotifRegressor::InitType initTypeL;
                int cvfolds;
                char gmmClusterFName[1024];
                char untransformedData[1024];
		//SK: genClusterAssign components
		map <string,int> speciesList;
		map<string,int> geneClusterAssignment;
		map<int,int> excludedOGList;
		map<int,int> ogidSet;
		bool expandInputMode;
};
#endif
