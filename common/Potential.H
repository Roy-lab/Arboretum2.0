#ifndef _POTENTIAL_
#define _POTENTIAL_
//The goal of this class is to learn a decision tree given a target node and its neighbourhood.
//Learning is done in two phases. In the first phase we perform consider one neighbour variable
//at a time and then split the node which improves the score the most.

#include <queue>
#include <list>
using namespace std;
#define PI 3.1472
#include "gsl/gsl_randist.h"
#include "CommonTypes.H"
#include "RegressionTree.H"
class Evidence;
class EvidenceManager;
class Variable;
class Potential
{
	public:
	Potential();
	~Potential();
	typedef enum
	{
		FACTOR,
		MARKOV_BNKT
	}VariableRole;
	int setAssocVariable(Variable*,VariableRole);
	VSET& getAssocVariables();
	int potZeroInit();
	int calculateConditionalEntropy();
	double getConditionalEntropy();
	double getJointEntropy();
	int setEvidenceManager(EvidenceManager*);
	int setMinLeafSize(int);
	int populateMe(double);
	int clearMe();
	double getJointPotValueFor(INTINTMAP&);
	double getJointPotValueForConf(string& varConf);
	int generateSample(INTINTMAP& jointConf, int vId,gsl_rng* r);
	double predictSample(map<int,Evidence*>* jointConf, int vId);
	double getMSE();
	int dumpPotential(ostream&);
	int makeValidJPD();
	RegressionTree* getRegressionTree();
	int showTree();
	int prune();
	int getAssocVariables_PostPruning(INTINTMAP&);

	private:
	int getPartitions(int vId, double& splitValue, double marginalEntropy,INTINTMAP& dataSet, double& infoGain);
	int getPartitions_Cached(int vId, double& splitValue, double marginalEntropy,INTINTMAP& dataSet, double& infoGain);
	int clearCache();
	int estimateMarginal(RegressionTree*);
	int sortAttrVals(DBLVECT&, INTVECT&);
	int getSubEntropy(INTVECT&,int,int,double&,double&,double&);
	int getSubEntropy(double,double,double&);
	double computeEntropyIfDeleted(Rule* r,int delMe,int&);
	int reducePruneSet(Rule* r);
	int computeCodingLength(RegressionTree*);
	int removeDuplicateCoverage();
	RegressionTree* dtree;
	EvidenceManager* evMgr;
	map<int,int> done;
	queue<RegressionTree*> currentLeafNodes;
	map<int,PARTITION*> allPartitions;
	map<int,INTDBLMAP*> allMeans;
	map<int,INTDBLMAP*> allVariances;
	map<int,INTDBLMAP*> allMarginalEntropy;
	VSET varSet;
	INTINTMAP factorVariables;
	INTINTMAP markovBlnktVariables;
	double conditionalEntropy;
	double jointEntropy;
	vector<Rule*> prunedRuleSet;
	double lambda;
	INTINTMAP pruneDataSet;
	INTDBLMAP minValMap;
	INTDBLMAP maxValMap;
	map<int,INTINTMAP*> dataRuleCoverage;
	int minLeafSize;
};
#endif
