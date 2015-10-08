#ifndef _BIONETWORK
#define _BIONETWORK
class BioNode;
class ProteinManager;
class GeneManager;
class InteractionManager;

class BioNetwork
{
	public:
		BioNetwork();
		~BioNetwork();
		int setProteinManager(ProteinManager*);
		int setGeneManager(GeneManager*);
		int setPPInteractionManager(InteractionManager*);
		int setPDInteractionManager(InteractionManager*);
		int createNetwork();
		int showNetworkStat();
		int getReachability();
		int getGraphInfo(const char*,int,map<int,int>&, map<string,int>&, map<string,int>&);
		int getDiffExpGenes(vector<int>&,int);
		int getAllNodeIDs(map<int,int>&);
	private:
		ProteinManager* pMgr;
		GeneManager* gMgr;
		InteractionManager* ppMgr;
		InteractionManager* pdMgr;
		map<string,BioNode*> bioNodeSet;
};
#endif
