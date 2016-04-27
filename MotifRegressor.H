#ifndef _MOTIF_REGRESSOR_
#define _MOTIF_REGRESSOR_
class VariableManager;
class EvidenceManager;
class MotifManager;
class ExpertL;
class Matrix;
class BFGSWrapper;
class MotifRegressor
{
	public:
		MotifRegressor();
		~MotifRegressor();
		int setVariableManager(VariableManager*);
		int setEvidenceManager(EvidenceManager*);
		int setMotifManager(MotifManager*);
		int setOutputDir(const char*);
		int setExpertCnt(int);
		int learnMoE();
		int learnMoE_CrossValidation(int neededfoldCnt);
		int showMoE();
		int showMoEParameters(const char*);
		int showGenatomyModule();
		int showClusterAssignment();
		int dispTFsPerCluster();
		int dispTFsPerCluster(int);
		int generateData(const char*,int);
		typedef enum
		{
			KMEANS,
			RAND,
			GMM
		}InitType;
		int setInitType(InitType);
		int setInitClusterFile(const char*);
		int predictTestData(map<string,map<string,int>*>&);
		int predictTestExpression();
		double getTestDataLikelihood(double&);
		int setUntransformedData(const char*);
	private:
		int initExperts();
		int initExpertsKMeans();
		int initExpertsGMM();
		int initVarsWithPredictors();

		int initPredictorsWithVars();
		int split(int);
		int cleanUp();
		int expectationStep();
		int maximizationStep();
		int normalizeWeights(map<int,double>&);
		int populateRandIntegers(vector<int>& randIndex,int evidCnt);
		int estimateMeanCov(vector<int>& randIndex,int startInd,int endInd,ExpertL* e);
		int estimateMeanCov(int expertid,ExpertL* e);
		int estimateGlobalCov(ExpertL* e);
		int estimateMean(vector<int>& randIndex,int startInd,int endInd,ExpertL* e);
		int setInitFeatureWeights(ExpertL* e);
		int populateBFGSData();
		int populateBFGSGammas();
		int updateFeatureWeights();
		double getScore();
		int matrixifyData();
		int assignGenesToExperts();
		int assignGenesToExperts_Exclusive();
		int assignGenesToExperts_Exclusive_ML();
		int getFeatureEnrichment(ExpertL* e,INTDBLMAP& featEnrich);
		int showGammas(const char*);
		int readClusterMembership();
		int showUntransformedFName();
		
		EvidenceManager* evMgr;
		VariableManager* varMgr;
		MotifManager* motifMgr;
		int expertCnt;
		map<int,ExpertL*> expertSet;
		BFGSWrapper bfgs;
		
		map<int,int> varsWithPredictors_Train;
		map<int,int> varsWithPredictors_Test;
		map<string,int> varsWithPredictors_Train_NameID;
		map<int,double> predictionCorr;
		map<int,double> predictionProb;
		map<int,int> varsWithPredictors;
		map<int,int> predsWithVars;
		map<int,map<int,double>*> gammas;
		map<int,Matrix*> exprProfileSet;
		map<int,int> featIDMatidMap;
		map<int,int> featMatidIDMap;
		gsl_rng* r;
		map<int,int> matIdvIdMap;
		char outputDir[1024];
		InitType initType;
		int foldCnt;
		map<int,double> holdoutLL;
		map<int,double> trainLL;
		map<int,double> holdout_unpenLL;
		char clusterFName[1024];
		char untransformedFName[1024];
		map<int,map<string,int>*> gmmClusters;
};
#endif
