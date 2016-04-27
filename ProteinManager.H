#ifndef _PROTEIN_MANAGER
#define _PROTEIN_MANAGER
#include <vector>
#include <map>
#include <string>
using namespace std;
class Protein;
class GeneManager;

typedef vector<Protein*> PROTSET;
typedef map<string,Protein*> PROTMAP;
typedef map<string,Protein*>::iterator PROTMAP_ITER;

class ProteinManager
{
	public:
		ProteinManager();
		~ProteinManager();
		int setGeneMgr(GeneManager*);
		int readProteinData(const char*);
		//The argument is a single line of the protein expression file
		int addProteinNode(const char*);
		Protein* getProteinNode(int id);
		Protein* getProteinWithName(const char*);
		int getTotalNumberOfProteins();
		int discretizeProteinExpr();
		int rankOrderProteinExpr();
		int scaleProteinExpr();
		int standardizeProteinExpr();
	private:
		PROTSET proteinSet;
		PROTMAP proteinMap;
		int proteinIds;
		GeneManager* geneMgr;
};
#endif
