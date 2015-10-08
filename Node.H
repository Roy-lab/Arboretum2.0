#ifndef _NODE_
#define _NODE_

#include <vector>
using namespace std;
#include "CommonDefs.H"


class Node
{
	public:
		Node();
		~Node();
		int setName(const char*);
		const char* getName();
		
		//This applies primarily to the gene and protein nodes
		int setIdInData(int id);
		int getIdInData();
		//Set the id according to topsort
		int setTopSortId(int id);
		int getTopSortId();
		typedef enum
		{
			GENE,
			PROTEIN,
			STAT_PG,
			STAT_PP,
			DYN_PG,
			DYN_PP
		}NodeType;
		
		int setNodeType(NodeType);
		NodeType getNodeType();

		int setValues(vector<int>&);
		vector<int>& getValues();
		//This is the string version of the value vector above
		const char* getValueStr();

		int addChild(Node*);
		int addParent(Node*);
		vector<Node*>& getChildren();
		vector<Node*>& getParents();
		int getNeighbourCnt();
	private:
		char nodeName[GENENAME_LEN];
		int nodeId;
		int topsortId;
		NodeType ntype;
		vector<Node*> children;
		vector<Node*> parents;
		vector<int> values;
		char valueStr[VALSTR_LEN];
};
#endif
