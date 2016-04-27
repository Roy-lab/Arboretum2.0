#ifndef _GENETREE_
#define _GENETREE_
#include <string>
using namespace std;
class GeneTree
{
	public:
	GeneTree();
	~GeneTree();
	GeneTree* leftchild;
	GeneTree* rightchild;
	GeneTree* parent;
	string species;
	string name;
	int nodeType;
	int nodeID;
	bool exclude;
};
#endif

