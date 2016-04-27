#ifndef _GENETREEMGR_
#define _GENETREEMGR_
#include <map>
using namespace std;
class GeneTree;
class SpeciesDistManager;
class MappedOrthogroup;
class GeneTreeManager
{
	public:
	GeneTreeManager();
	~GeneTreeManager();
	int setSpeciesTree(SpeciesDistManager::Species*);
	int setSpeciesDistManager(SpeciesDistManager*);
	int setGeneTreeDir(const char*);
	GeneTree* getGeneTree(int);	
	GeneTree* getGeneTree(MappedOrthogroup*);	
	GeneTree* generateTree(MappedOrthogroup*);
	GeneTree* generateTreeFromFile(MappedOrthogroup*);
	inline int getTotalAdditionalParents()
	{
		return totalAdditionalParents;
	}
	int showTotalAdditionalParentDist();
	int countLossForTree(int);
	private:
	int showTree(GeneTree*);
	int pruneTree(GeneTree*);
	int insertTree(GeneTree* node);
	int setNodeType(GeneTree* node);
	int identifyAncestorWithDuplicates(SpeciesDistManager::Species* speciesnode,MappedOrthogroup* mor,SpeciesDistManager::Species** duplicateSpecies);
	const char* getSpeciesName(GeneTree* node);
	int getSpeciesName(GeneTree* node,map<string,int>&);
	GeneTree* createTree(SpeciesDistManager::Species* duplicateSpecies,SpeciesDistManager::Species* node, MappedOrthogroup* mor, int copyno);
	GeneTree* readTreeFromFile(const char* treePath);
	int addKpolGene(GeneTree*,string&);
	GeneTree* isKpolAnc(GeneTree* gtree,string& parentName);
	int countLossForTree(GeneTree*);
	map<int,GeneTree*> genetreeSet;
	SpeciesDistManager::Species* speciestreeRoot;
	SpeciesDistManager* sdMgr;
	//SK: variable for gene tree direrctory, to be used if readung in from files instead of the MappedOrthogroupReader object
	char dirName[1024];
	int totalAdditionalParents;
	map<string,int> addedParents;
};
#endif
