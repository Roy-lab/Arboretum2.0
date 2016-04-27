#include <iostream>
#include <string.h>
using namespace std;

#include "Path.H"
#include "BioNode.H"
#include "Gene.H"
#include "Protein.H"
#include "ProteinManager.H"
#include "GeneManager.H"
#include "Interaction.H"
#include "InteractionManager.H"
#include "BioNetwork.H"

BioNetwork::BioNetwork()
{
}

BioNetwork::~BioNetwork()
{
}

int
BioNetwork::setProteinManager(ProteinManager* mPtr)
{
	pMgr=mPtr;
	return 0;
}

int
BioNetwork::setGeneManager(GeneManager* mPtr)
{
	gMgr=mPtr;
	return 0;
}

int
BioNetwork::setPPInteractionManager(InteractionManager* imPtr)
{
	ppMgr=imPtr;
	return 0;
}

int
BioNetwork::setPDInteractionManager(InteractionManager* imPtr)
{
	pdMgr=imPtr;
	return 0;
}

int
BioNetwork::createNetwork()
{
	//Get all the nodes from the protein manager and gene manager
	int geneNodeCnt=gMgr->getTotalNumberOfGenes();
	int proteinNodeCnt=pMgr->getTotalNumberOfProteins();
	for(int i=0;i<geneNodeCnt;i++)
	{
		//cout << "Read gene at " << i << endl;
		Gene* gene=gMgr->getGeneNode(i);
		string key(gene->getName());
		BioNode* genebNode=(BioNode*) gene;
		bioNodeSet[key]=genebNode;
		if(proteinNodeCnt==0)
		{
			continue;
		}
		const char* codedPName=gene->getCodedProteinName();
		BioNode* proteinbNode=(BioNode*)pMgr->getProteinWithName(codedPName);

		genebNode->addNeighbour(proteinbNode,BioNode::OUT_NEIGHBOUR);
		proteinbNode->addNeighbour(genebNode,BioNode::IN_NEIGHBOUR);
	}
	for(int j=0;j<proteinNodeCnt;j++)
	{
		Protein* protein=pMgr->getProteinNode(j);
		string key(protein->getName());
		bioNodeSet[key]=(BioNode*) protein;
	}
	//Now we will look the protein dna interactions and add the in and out neighbours of the gene and protein nodes
	//respectively
	int pdIntrCnt=pdMgr->getNumberOfIntr();
	for(int i=0;i<pdIntrCnt;i++)
	{
		Interaction* aIntr=pdMgr->getInteractionAt(i);
		Protein* p1=pMgr->getProteinWithName(aIntr->getFirstMember());
		Gene* g2=gMgr->getGeneWithName(aIntr->getSecondMember());
		p1->addNeighbour((BioNode*)g2,BioNode::OUT_NEIGHBOUR);
		g2->addNeighbour((BioNode*)p1,BioNode::IN_NEIGHBOUR);
	}
	int ppIntrCnt=ppMgr->getNumberOfIntr();
	for(int i=0;i<ppIntrCnt;i++)
	{
		Interaction* aIntr=ppMgr->getInteractionAt(i);
		string key1(aIntr->getFirstMember());
		string key2(aIntr->getSecondMember());
		if(bioNodeSet.find(key1)==bioNodeSet.end())
		{
			cout <<"No node named " << key1.c_str() << endl;
			return -1;
		}
		if(bioNodeSet.find(key2)==bioNodeSet.end())
		{
			cout <<"No node named " << key2.c_str() << endl;
			return -1;
		}
		BioNode* b1=bioNodeSet[key1];
		BioNode* b2=bioNodeSet[key2];
		b1->addNeighbour(b2,BioNode::OUT_NEIGHBOUR);
		b1->addNeighbour(b2,BioNode::IN_NEIGHBOUR);
		b2->addNeighbour(b1,BioNode::OUT_NEIGHBOUR);
		b2->addNeighbour(b1,BioNode::IN_NEIGHBOUR);
	}

	return 0;
}


int 
BioNetwork::showNetworkStat()
{
	cout <<"Number of bionodes "<< bioNodeSet.size() << endl;
	for(map<string,BioNode*>::iterator aIter=bioNodeSet.begin();aIter!=bioNodeSet.end();aIter++)
	{
		//cout<< aIter->first.c_str() << "\t" << aIter->second->getPaths().size() << endl;
		if(aIter->second->getType()==BioNode::GENE)
		{
			cout << aIter->first.c_str() << "\t";
			aIter->second->showInNeighbours();
		}
	}
	return 0;
}


int
BioNetwork::getReachability()
{
	//For reachability, we will call findReachable nodes on every BioNode.
	//We will store the paths (shortest) for each node alongwith its length in every bionode
	for(map<string,BioNode*>::iterator aIter=bioNodeSet.begin();aIter!=bioNodeSet.end();aIter++)
	{
		BioNode* bnode=aIter->second;
		bnode->findReachableNodes();
	}
	return 0;
}

int
BioNetwork::getGraphInfo(const char* hubGene,int pathLenLimit, map<int,int>& geneIDs, 
			 map<string,int>& pdIntr, map<string,int>& ppIntr)
{
	string key(hubGene);
	if(bioNodeSet.find(key)==bioNodeSet.end())
	{
		cout <<"Hub gene " << hubGene << " not found " << endl;
		return -1;
	}
	BioNode* bnode=bioNodeSet[key];
	//bnode->dumpPaths();
	geneIDs[bnode->getID()]=0;
	
	map<string, Path*>& pathSet=bnode->getPaths();
	map<string, int> shownPaths;
	for(map<string,Path*>::iterator pIter=pathSet.begin();pIter!=pathSet.end();pIter++)
	{
		//We want all nodes on this path which are of length equal or less than pathLenLimit
		Path* apath=pIter->second;
		//Just for display purposes
		string pathStr;
		BioNode* prevnode=bnode;
		pathStr.append(prevnode->getName());
		int pid=1;
		char intrID[256];
		while((pid<=pathLenLimit) && (pid<apath->getPathLength()))
		{
			BioNode* nextnode=apath->getNodeAt(pid);
			geneIDs[nextnode->getID()]=0;
			pathStr.append("-");
			pathStr.append(nextnode->getName());

			//Add the interaction
			if((prevnode->getType()==BioNode::PROTEIN) &&(nextnode->getType()==BioNode::PROTEIN))
			{
				sprintf(intrID,"%d-%d",prevnode->getID(),nextnode->getID());
				ppIntr[intrID]=0;
				
			}
			else if((prevnode->getType()==BioNode::PROTEIN) && (nextnode->getType()==BioNode::GENE))
			{
				sprintf(intrID,"%d-%d",prevnode->getID(),nextnode->getID());
				pdIntr[intrID]=0;
			}
			else if((prevnode->getType()==BioNode::GENE) && (nextnode->getType()==BioNode::PROTEIN))
			{
				//cout <<"Ignoring gene-protein interaction " << endl;	
			}
			else
			{
				cout <<"Wrong type of interaction between " << prevnode->getName() 
					<< " and " << nextnode->getName() << endl;
				return -1;
			}
			prevnode=nextnode;
			pid++;
		}
		if(shownPaths.find(pathStr)==shownPaths.end())
		{
			cout << pathStr.c_str() << endl;
			shownPaths[pathStr]=0;
		}
	}

	return 0;
}


int 
BioNetwork::getDiffExpGenes(vector<int>& expIds, int pertGeneId)
{
	cout << "The differentially expressed genes" << endl;
	BioNode* pertNode=NULL;
	for(int i=0;i<expIds.size();i++)
	{
		vector<string> diffExpGenes;
		for(map<string,BioNode*>::iterator aIter=bioNodeSet.begin();aIter!=bioNodeSet.end();aIter++)
		{
			BioNode* anode=aIter->second;
			if(anode->getID()==pertGeneId)
			{
				pertNode=anode;
			}
			if(anode->getType()==BioNode::GENE)
			{
				int level=anode->getDiscreteExpLevelAt(expIds[i]);
				if(level!=BioNode::MEDIUM)
				{
					diffExpGenes.push_back(aIter->first);
				}
			}
		}
		cout <<"Exp ID " << expIds[i];
		for(int i=0;i<diffExpGenes.size();i++)
		{
			cout << " " << diffExpGenes[i].c_str();
		}
		cout << endl;
		cout <<"Paths to " << pertGeneId << endl;
		for(int i=0;i<diffExpGenes.size();i++)
		{
			if(strcmp(pertNode->getName(),diffExpGenes[i].c_str())!=0)
			{
				pertNode->showPathTo(diffExpGenes[i].c_str());
			}
		}
	}
	return 0;
}

int 
BioNetwork::getAllNodeIDs(map<int,int>& geneIDs)
{
	for(map<string,BioNode*>::iterator aIter=bioNodeSet.begin();aIter!=bioNodeSet.end();aIter++)
	{
		BioNode* anode=aIter->second;
		if(anode->getType()==BioNode::GENE)
		{
			geneIDs[anode->getID()]=0;
		}
	}
	return 0;
}
