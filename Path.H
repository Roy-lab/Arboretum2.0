//This stores a directed path of BioNodes. Some of the edges in this path are undirected.
//In principle I could store names of nodes, but this would mean I am storing strings which is
//expensive. Instead storing pointers to BioNodes is better
#ifndef _PATH
#define _PATH

#include <vector>
#include <fstream>

using namespace std;

class BioNode;

class Path
{
	public:
		Path();
		~Path();
		int setStartNode(BioNode*);
		int addNode(BioNode*);
		int getPathLength();
		BioNode* getNodeAt(int);
		Path* copy();	
		int showPath(ostream&);
		int hasBeenCopied();
		bool isCopied();
	private:
		vector<BioNode*> nodeSet;
		bool copyFlag;
};

#endif
