#include <iostream>
#include <fstream>
#include <string.h>

#include "Node.H"
#include "Graph.H"

Graph::Graph()
{
}

Graph::~Graph()
{
}

int 
Graph::addNode(const char* nName,Node::NodeType ntype,int idInData, vector<int>& values)
{
	Node* aNode=new Node;
	aNode->setName(nName);
	aNode->setIdInData(idInData);
	aNode->setNodeType(ntype);
	aNode->setValues(values);
	string nodeKey(nName);
	nodeSet[nodeKey]=aNode;
	return 0;
}

int
Graph::addEdge(const char* node1, const char* node2)
{
	string key1(node1);
	string key2(node2);
	Node* v1=nodeSet[key1];
	Node* v2=nodeSet[key2];
	v1->addChild(v2);
	v2->addParent(v1);
	return 0;
}

int 
Graph::doTopologicalSort()
{
	//Do a DFS to obtain top-sort as described in CLR
	//First make a vector of ints to make sure that we have visited all the vertices
	//0 -> white
	//1 -> gray
	//2 -> black
	for(VERTEX_SET_ITER vsIter=nodeSet.begin();vsIter!=nodeSet.end();vsIter++)
	{
		visitFlag[vsIter->first]=0;
	}
	int visitCnt=0;
	int currTime=0;
	//Outer loop is required just to make sure I have visited all the vertices
	for(VERTEX_SET_ITER vsIter=nodeSet.begin();vsIter!=nodeSet.end();vsIter++)
	{
		if(visitFlag[vsIter->first]==0)
		{
			//Do the DFS-visit, i.e. set its color to gray and push it on to the queue
			Node* aNode=vsIter->second;
			currTime=dfsVisit(aNode,currTime);
		}
	}
	//Show the top sorted vector
	for(int i=topSorted_R.size()-1;i>=0;i--)
	{
		topSorted.push_back(topSorted_R[i]);
	}
	/*for(int i=0;i<topSorted.size();i++)
	{
		cout << topSorted[i]->getName()<< endl;
	}*/
	return 0;
}

vector<Node*>&
Graph::getTopologicalSort()
{
	return topSorted;
}

int
Graph::dfsVisit(Node* aNode, int currTime)
{
	string myStr(aNode->getName());
	visitFlag[myStr]=1;
	int endTime=currTime;
	vector<Node*> children=aNode->getChildren();
	for(int i=0;i<children.size();i++)
	{
		endTime++;
		string childStr(children[i]->getName());
		if(visitFlag[childStr]==0)
		{
			dfsVisit(children[i],currTime);
		}
	}
	topSorted_R.push_back(aNode);
	visitFlag[myStr]=2;
	return endTime;
}
/*The information required by PNL is the 
 * type of a node i.e whether it is discrete or continuous.
 * Number of neighbours of each node
 * For each neighbour, whether it a parent or a child
 * We are going to generate a file that has the information
 * in the order in which the PNL toolkit expects it to be.
 * */

int
Graph::genPNLInputFormat(const char* aFName)
{
	ofstream oFile(aFName);
	//First assign ids for each node in the topsort
	//While doing this also get the continuous nodes and the discrete nodes
	vector<int> contNodeID;
	vector<int> discreteNodeID;
	for(int i=0;i<topSorted.size();i++)
	{
		Node* aNode=topSorted[i];
		aNode->setTopSortId(i);
		if( (strstr(aNode->getName(),"pp")!=NULL) || 
		    (strstr(aNode->getName(),"pg")!=NULL) ||
		    (strstr(aNode->getName(),"gg")!=NULL) ||
		    (strstr(aNode->getName(),"gp")!=NULL))
		{
			discreteNodeID.push_back(i);
		}
		else
		{
			contNodeID.push_back(i);
		}
	}
	//Dump the number of nodes
	oFile <<"NodeCnt\t"<< topSorted.size() << endl;
	//Dump the continuous and discrete node ids, ideally we want only the discrete nodes or the
	//continuous nodes but we will use both
	oFile<<"ContinuousNodes";
	for(int i=0;i<contNodeID.size();i++)
	{
		oFile <<"\t"<< contNodeID[i];
	}
	oFile << endl;
	
	oFile << "DiscreteNodes";
	for(int i=0;i<discreteNodeID.size();i++)
	{
		oFile <<"\t" << discreteNodeID[i];
	}
	oFile << endl;
	
	//Now for each node output the number of neighbours
	/*oFile <<"Neighbours";
	for(int i=0;i<topSorted.size();i++)
	{
		Node* aNode=topSorted[i];
		oFile <<"\t" << aNode->getNeighbourCnt();
	}
	oFile << endl;*/

	//Now specify the actual neighbours first the parents and then the children
	for(int i=0;i<topSorted.size();i++)
	{
		Node* aNode=topSorted[i];
		oFile <<"NodeName=" << aNode->getName() <<"\tNodeID="<<aNode->getTopSortId() <<"\tParents=";
		vector<Node*>& parents=aNode->getParents();
		for(int j=0;j<parents.size();j++)
		{
			oFile<<parents[j]->getTopSortId();
			if(j!=(parents.size()-1))
			{
				oFile << ",";
			}
		}
		oFile << "\tChildren=";
		vector<Node*>& children=aNode->getChildren();
		for(int j=0;j<children.size();j++)
		{
			oFile << children[j]->getTopSortId();
			if(j!=(children.size()-1))
			{
				oFile << ",";
			}
		}
		oFile << "\tValues="<< aNode->getValueStr();
		oFile << endl;
	}
	oFile.close();
	return 0;
}
