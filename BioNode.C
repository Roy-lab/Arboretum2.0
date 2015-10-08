#include <iostream>
#include <queue>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "Path.H"
#include "BioNode.H"

BioNode::BioNode()
{
	expMin=10000;
}

BioNode::~BioNode()
{
}

int 
BioNode::setName(const char* aName)
{
	strcpy(name,aName);
	return 0;
}

int 
BioNode::setID(int aId)
{
	nId=aId;
	return 0;
}

const char* 
BioNode::getName()
{
	return name;
}

int 
BioNode::getID()
{
	return nId;
}

int 
BioNode::setType(BioNode::BioNodeType bnType)
{
	nodeType=bnType;
	return 0;
}

BioNode::BioNodeType 
BioNode::getType()
{
	return nodeType;
}

int 
BioNode::addExpLevel(double aVal)
{
	if((expMin>aVal) && (aVal>0))
	{
		expMin=aVal;
	}

	expLevels.push_back(aVal);
	return 0;
}

double 
BioNode::getExpLevelAt(int expId)
{
	if((expId<0) || (expId>=expLevels.size()))
	{
		cout << "Bad ID for experiment id  "<< endl;
		return -1;
	}
	double aVal=expLevels[expId];
	if(aVal<=0)
	{
	//	cout <<"Updating " << aVal << " to " << expMin << " for " << name << endl;
	//	aVal=expMin;
	}
	return aVal;
}

double
BioNode::getRankedExpLevelAt(int expId)
{
	if((expId<0) || (expId>=rankedExpLevels.size()))
	{
		cout << "Bad ID for experiment id  "<< endl;
		return -1;
	}
	return rankedExpLevels[expId];
}

//Right now we will discretize into the three levels by computing mean
//and stdev of the expression levels. We'll use 2sd difference from mean
//to get the levels
int 
BioNode::discretizeLevels()
{
	double mExp=getMean();
	double sdExp=getStdev();
	for(int i=0;i<expLevels.size();i++)
	{
		double diff=expLevels[i]-mExp;
		if(diff<0)
		{
			if(fabs(diff)>=(1*sdExp))
			{
				discreteExpLevels.push_back(BioNode::LOW);
			}
			else
			{
				discreteExpLevels.push_back(BioNode::MEDIUM);
			}
		}
		else if(diff>0)
		{
			if(fabs(diff)>=(1*sdExp))
			{
				discreteExpLevels.push_back(BioNode::HIGH);
			}
			else
			{
				discreteExpLevels.push_back(BioNode::MEDIUM);
			}
		}
		else
		{
			discreteExpLevels.push_back(BioNode::MEDIUM);
		}
	}
	return 0;
}

//Use rankorder or quantile method. This will allow a more or less
//equal presence of the different values
//For now we will hard code the number of levels
int
BioNode::rankOrder()
{
	int level=3;
	INTVECT ranks;
	sort(ranks);
	int maxRank=ranks.size();
	double binSize=((double)maxRank)/3.0;
	for(int g=0;g<ranks.size();g++)
	{
		int modRank=(int) floor(ranks[g]/binSize);
		discreteExpLevels.push_back(modRank);
	}
	return 0;
}

int
BioNode::rankOrder_Store()
{
	INTVECT ranks;
	sort(ranks);
	double maxRank=(double)ranks.size();
	for(int g=0;g<ranks.size();g++)
	{
		double r=(ranks[g]+0.5)/maxRank;
		rankedExpLevels.push_back(r);
	}
	return 0;
}


int
BioNode::discretizeEqualBins()
{
	double max=-1;
	double min=-1;
	getMaxMin(max,min);
	double binSize=(max-min)/3;
	for(int g=0;g<expLevels.size();g++)
	{
		int dVal=0;
		if(binSize>0)
		{
			dVal=(int) floor(expLevels[g]/binSize);
		}
		if(dVal>2)
		{
			dVal=2;
		}
		discreteExpLevels.push_back(dVal);
	}
	
	return 0;
}

int
BioNode::discretizeEqualBins_Global(double min, double max)
{
	double binSize=(max-min)/3;
	for(int g=0;g<expLevels.size();g++)
	{
		if(expLevels[g]>200)
		{
			cout <<"found biggie" << endl;
		}
		int dVal=0;
		if(binSize>0)
		{
			//dVal=(int) floor((expLevels[g]-min)/binSize);
			dVal =(int) floor(3*((log(expLevels[g])-log(min))/(log(max)-log(min))));
		}
		if(dVal>2)
		{
			dVal=2;
		}
		discreteExpLevels.push_back(dVal);
	}
	
	return 0;
}

int
BioNode::scaleExpression()
{
	double max=-1;
	double min=-1;
	getMaxMin(max,min);
	for(int g=0;g<expLevels.size();g++)
	{
		expLevels[g]=expLevels[g]/max;
	}
	return 0;
}

int
BioNode::standardizeData()
{
	calculateMean();	
	calculateStdev(mean);	
	double mean=getMean();
	double std=getStdev();
	for(int g=0;g<expLevels.size();g++)
	{
		expLevels[g]=exp((log(expLevels[g]+1)-mean)/std);
		if(isnan(expLevels[g]))
		{
			cout <<"Stop! Found nan for " << name << endl;
		}
	}
	return 0;
}


int
BioNode::zeroMeanData()
{
	calculateMean();	
	double mean=getMean();
	for(int g=0;g<expLevels.size();g++)
	{
		expLevels[g]=exp((log(expLevels[g]+1)-mean));
		if(isnan(expLevels[g]))
		{
			cout <<"Stop! Found nan for " << name << endl;
		}
	}
	return 0;
}

int
BioNode::getDiscreteExpLevelAt(int expId)
{
	return discreteExpLevels[expId];
}

int
BioNode::calculateMean()
{
	mean=0;
	for(int i=0;i<expLevels.size();i++)
	{
		mean=mean+log(expLevels[i]+1);
	}
	mean=mean/(double)expLevels.size();
	return 0;
}

int
BioNode::calculateStdev(double mean)
{
	stdev=0;
	for(int i=0;i<expLevels.size();i++)
	{
		double diff=log(expLevels[i]+1)-mean;
		stdev=stdev+(diff*diff);
	}
	stdev=sqrt(stdev/(double) (expLevels.size()-1));
	
	return 0;
}

double
BioNode::getMean()
{
	return mean;
}

double
BioNode::getStdev()
{
	return stdev;
}

double
BioNode::getMaxMin(double& max, double& min)
{
	for(int i=0;i<expLevels.size();i++)
	{
		if(max==-1 || expLevels[i]>max)
		{
			max=expLevels[i];
		}
		if(min==-1 || expLevels[i]<min)
		{
			min=expLevels[i];
		}
	}
	return 0;
}

int
BioNode::addNeighbour(BioNode* bNode, BioNode::NeighbourType ntype)
{
	string key(bNode->getName());
	if(ntype==BioNode::IN_NEIGHBOUR)
	{
		if(inNeighbours.find(key)==inNeighbours.end())
		{
			inNeighbours[key]=bNode;
		}
		else
		{
			//cout <<"Already added in-neighbour " << key.c_str() << " to " << name << endl;
		}
	}
	else if (ntype==BioNode::OUT_NEIGHBOUR)
	{
		if(outNeighbours.find(key)==outNeighbours.end())
		{
			outNeighbours[key]=bNode;
		}
		else
		{
			//cout <<"Already added out-neighbour " << key.c_str()<< " to " << name << endl;
		}
		//out neighbours are reachable from this node using a single edge
		reachableNeighbours[key]=1;
	}
		
	return 0;
}

BNODE_MAP& 
BioNode::getOutNeighbours()
{
	return outNeighbours;
}

BNODE_MAP& 
BioNode::getInNeighbours()
{
	return inNeighbours;
}

int
BioNode::findReachableNodes()
{
	//lazy bum me is going to STL queue of vertices
	queue<BioNode*> reachableQ;
	//Initialize reachableQ with immediate neighbours
	for(BNODE_MAP_ITER bnIter=outNeighbours.begin();bnIter!=outNeighbours.end();bnIter++)
	{
		reachableQ.push(bnIter->second);
		Path* newPath=new Path;
		newPath->setStartNode(this);
		newPath->addNode(bnIter->second);
		pathSet[bnIter->first]=newPath;
	}

	while(!reachableQ.empty())
	{
		BioNode* v=reachableQ.front();
		//get the neighbours of aNode
		BNODE_MAP& vNeighbours=v->getOutNeighbours();
		//This node, u, is reachable to all the neighbours of v by 1 + dist(u,v)
		string vkey(v->getName());
		if(reachableNeighbours.find(vkey)==reachableNeighbours.end())
		{
			cout <<"Expected a node in reachable list, did not find " << endl;
			exit(0);
		}
		int dist_uv=reachableNeighbours[vkey];
		map<string,Path*> tempPathSet;
		for(BNODE_MAP_ITER bnIter=vNeighbours.begin();bnIter!=vNeighbours.end();bnIter++)
		{
			//Check for cycles
			if(reachableNeighbours.find(bnIter->first)==reachableNeighbours.end())
			{
				reachableNeighbours[bnIter->first]=dist_uv+1;
				if(pathSet.find(vkey)==pathSet.end())
				{
					cout <<"Did not find existing path  for " << vkey.c_str() << endl;
					return -1;
				}
					
				Path* existingPath=pathSet[vkey];
				existingPath->hasBeenCopied();
				Path* newPath=existingPath->copy();
				newPath->addNode(bnIter->second);
				tempPathSet[bnIter->first]=newPath;
				reachableQ.push(bnIter->second);
			}
		}
		reachableQ.pop();
		//update pathSet
		for(map<string,Path*>::iterator pIter=pathSet.begin();pIter!=pathSet.end();pIter++)
		{
			Path* aPath=pIter->second;
			if(aPath->isCopied())
			{
				delete aPath;
				pathSet.erase(pIter);
			}
		}
		for(map<string,Path*>::iterator pIter=tempPathSet.begin();pIter!=tempPathSet.end();pIter++)
		{
			pathSet[pIter->first]=pIter->second;	
		}
		tempPathSet.clear();
	}
	buildReachableMap();
	buildNodePathAssoc();
	return 0;
}

int
BioNode::getReachableNodeCnt()
{
	return reachableNeighbours.size();
}

map<string,Path*>&
BioNode::getPaths()
{
	return pathSet;
}


//Store the path length and the set of nodes that are reachable from this node at this length
int 
BioNode::buildReachableMap()
{
	for(STRINTMAP_ITER aIter=reachableNeighbours.begin();aIter!=reachableNeighbours.end();aIter++)
	{
		int length=aIter->second;
		if(neighbourSets.find(length)==neighbourSets.end())
		{
			STRVECT* sVect=new STRVECT;
			neighbourSets[length]=sVect;
			sVect->push_back(aIter->first);
		}
		else
		{
			STRVECT* sVect=neighbourSets[length];
			sVect->push_back(aIter->first);
		}
	}
	return 0;
}

//Store the path in which a node and this node have their corresponding shortest path
int 
BioNode::buildNodePathAssoc()
{
	for(map<string,Path*>::iterator pIter=pathSet.begin();pIter!=pathSet.end();pIter++)
	{
		Path* apath=pIter->second;
		int pathLen=apath->getPathLength();
		for(int i=0;i<pathLen;i++)
		{
			BioNode* anode=apath->getNodeAt(i);
			string aKey(anode->getName());
			nodePathAssoc[aKey]=pIter->first;
		}
	}
	return 0;
}

int
BioNode::dumpPaths()
{
	cout <<"Total number of paths" << pathSet.size() << endl;
	cout <<"Path length\tNumber of nodes " << endl;
	for(map<int,STRVECT*>::iterator aIter=neighbourSets.begin();aIter!=neighbourSets.end();aIter++)
	{
		cout << aIter->first <<"\t" << aIter->second->size() <<endl;
	}
	for(map<string,Path*>::iterator pIter=pathSet.begin();pIter!=pathSet.end();pIter++)
	{
	//	pIter->second->showPath(cout);
	}
	return 0;
}

int
BioNode::showInNeighbours()
{
	for(BNODE_MAP_ITER bIter=inNeighbours.begin();bIter!=inNeighbours.end();bIter++)
	{
		cout << " " << bIter->first.c_str();
	}
	cout << endl;
	return 0;
}

int
BioNode::showPathTo(const char* nodeName)
{
	string nodeKey(nodeName);
	if(nodePathAssoc.find(nodeKey)!=nodePathAssoc.end())
	{
		string& pathKey=nodePathAssoc[nodeKey];
		Path* apath=pathSet[pathKey];
		apath->showPath(cout);
	}
	return 0;
}


int 
BioNode::sort(INTVECT& rank)
{
	for(int i=0;i<expLevels.size();i++)
	{
		rank.push_back(i);
	}
	for(int i=0;i<expLevels.size();i++)
	{
		for(int j=i+1;j<expLevels.size();j++)
		{
			if(expLevels[i]>expLevels[j])
			{
				double temp=expLevels[i];
				expLevels[i]=expLevels[j];
				expLevels[j]=temp;
				int indTemp=rank[i];
				rank[i]=rank[j];
				rank[j]=indTemp;
			}
		}
	}
	return 0;
}
