#ifndef _GRAPH
#define _GRAPH
#include <map>
#include <string>
#include <vector>
#include "Node.H"
using namespace std;
typedef map<string,Node*> VERTEX_SET;
typedef map<string,Node*>::iterator VERTEX_SET_ITER;

class Node;
class Graph
{
	public:
		Graph();
		~Graph();
		int addNode(const char*,Node::NodeType,int,vector<int>&);
		int addEdge(const char*,const char*);
		int doTopologicalSort();
		vector<Node*>& getTopologicalSort();
		int genPNLInputFormat(const char*);
	private:
		int dfsVisit(Node*,int);
		VERTEX_SET nodeSet;
		vector<Node*> topSorted_R;
		vector<Node*> topSorted;
		map<string,int> visitFlag;
};
#endif
