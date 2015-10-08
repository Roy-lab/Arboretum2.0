#include "BioNode.H"
#include "Path.H"


Path::Path()
{
	copyFlag=false;
}

Path::~Path()
{
	nodeSet.clear();
}

int 
Path::setStartNode(BioNode* aNode)
{
	nodeSet.push_back(aNode);
	return 0;
}

int 
Path::addNode(BioNode* aNode)
{
	nodeSet.push_back(aNode);
	return 0;
}
	
int 
Path::getPathLength()
{
	return nodeSet.size();
}


BioNode*
Path::getNodeAt(int nodeId)
{
	return nodeSet[nodeId];
}

Path*
Path::copy()
{
	Path* newPath=new Path;
	for(int i=0;i<nodeSet.size();i++)
	{
		newPath->addNode(nodeSet[i]);
	}
	return newPath;
}

int 
Path::showPath(ostream& oFile)
{
	for(int i=0;i<nodeSet.size();i++)
	{
		if(i>0)
		{
			oFile <<" ";
		}
		oFile << nodeSet[i]->getName();
	}
	oFile << endl;
	return 0;
}


int 
Path::hasBeenCopied()
{
	copyFlag=true;
	return 0;
}

bool 
Path::isCopied()
{
	return copyFlag;
}
