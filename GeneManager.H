#ifndef _GENE_MANAGER
#define _GENE_MANAGER
#include <vector>
#include <map>
#include <string>
using namespace std;
class Gene;

typedef vector<Gene*> GENESET;
typedef map<string,Gene*> GENEMAP;
typedef map<string,Gene*>::iterator GENEMAP_ITER;

class GeneManager
{
	public:
		GeneManager();
		~GeneManager();
		int readGeneData(const char*, int logtrans=0);
		int readGeneData(const char*, int datastart, int logtrans=0);
		//The parameter is a single line of the protein expression file
		int addGeneNode(const char*,int logtrans=0);
		Gene* getGeneNode(int id);
		Gene* getGeneWithName(const char*);
		int getTotalNumberOfGenes();
		int discretizeGeneExpr();
		int rankOrderGeneExpr();
		int scaleGeneExpr();
		int standardizeGeneExpr();
		int zeroMeanGeneExpr();

		int assocProtWithCodingGene(const char*, int);
	private:
		GENESET geneSet;
		GENEMAP geneMap;
		int geneIds;
		int totalMeasurements;
		double maxExp;
		double minExp;
};
#endif
