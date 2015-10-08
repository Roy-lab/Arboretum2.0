/*The goal of this class is to interface the factor graph classes with the kmeans or other
 * clustering programs. This class first checks if there any existing clusters generated
 * via non-negative matrix factorizations or previous run of kmeans. If not it gets
 * data from the evidence manager and then uses kmeans to generate the clusters*/

#include <vector>
using namespace std;

class VariableManager;
class EvidenceManager;
class Kmeans;

class ClusterManager
{
	public:
	ClusterManager();
	~ClusterManager();
	int setVariableManager(VariableManager* aPtr);
	int setEvidenceManager(EvidenceManager* aPtr);
	int setOutputDir(const char*);
	int setClusterCnt(int);
	int setMaxIter(int);
	vector<VSET*>& getClusters();
	vector<INTINTMAP*>& getClusters(map<int,INTDBLMAP*>& kmeansData);
	vector<VSET*>& getRandomClusters();
	private:
	int readClusters();
	int generateClusters();
	int generateClusters_Randomrestarts();
	int generateClusters_Randomrestarts(map<int,INTDBLMAP*>&);
	int maxIter;
	Kmeans kmeans;
	int clusterCnt;
	char outputDir[1024];
	VariableManager* vMgr;
	EvidenceManager* evidMgr;
	vector<VSET*> variableSubsets;
	vector<INTINTMAP*> variableSubsets_IDs;
	vector<VSET*> randomVariableSubsets;
};



