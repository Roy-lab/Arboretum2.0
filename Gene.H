#ifndef _GENE
#define _GENE
#include <vector>

using namespace std;

#include  "BioNode.H"

class Gene : public BioNode
{
	public:
		Gene();
		~Gene();
		int setCodedProteinName(const char* protName);
		const char* getCodedProteinName();
	private:
		char assocProtName[GENENAME_LEN];
};
#endif
