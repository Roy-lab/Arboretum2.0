#include <fstream>
#include <iostream>
#include <string.h>

#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"
#include "EvidenceManager.H"
#include "Randomizer.H"
#include "Kmeans.H"
#include "ClusterManager.H"

ClusterManager::ClusterManager()
{
	maxIter=50;
}

ClusterManager::~ClusterManager()
{
}

int 
ClusterManager::setVariableManager(VariableManager* aPtr)
{
	vMgr=aPtr;
	return 0;
}

int 
ClusterManager::setEvidenceManager(EvidenceManager* aPtr)
{
	evidMgr=aPtr;
	return 0;
}

int 
ClusterManager::setOutputDir(const char* aFName)
{
	strcpy(outputDir,aFName);
	return 0;
}

int
ClusterManager::setClusterCnt(int aCnt)
{
	clusterCnt=aCnt;
	return 0;
}

int 
ClusterManager::setMaxIter(int aIter)
{
	maxIter=aIter;
	return 0;
}

vector<VSET*>& 
ClusterManager::getClusters()
{
	readClusters();
	if(variableSubsets.size()==0)
	{
		//generateClusters();
		generateClusters_Randomrestarts();
	}
	return variableSubsets;
}

vector<INTINTMAP*>&
ClusterManager::getClusters(map<int,INTDBLMAP*>& kmeansData)
{
	generateClusters_Randomrestarts(kmeansData);
	return variableSubsets_IDs;
}

vector<VSET*>&
ClusterManager::getRandomClusters()
{
	readClusters();
	if(variableSubsets.size()==0)
	{
		cout <<"Did not find clusters!! " << endl;
	}
	VSET& variableSet=vMgr->getVariableSet();
	Randomizer rand;
	rand.initialize(0,variableSet.size()-1);
	char aFName[1024];
	sprintf(aFName,"%s/km_randclusters_%d.txt",outputDir,clusterCnt);
	ofstream oFile(aFName);
	map<int,int> usedVariables;
	for(int i=0;i<variableSubsets.size();i++)
	{
		VSET* trueCluster=variableSubsets[i];
		VSET* rndCluster=new VSET;
		randomVariableSubsets.push_back(rndCluster);
		oFile <<"ClusterId " << i << endl;
		if(i==(variableSubsets.size()-1))
		{
			for(VSET_ITER vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
			{
				if(usedVariables.find(vIter->first)==usedVariables.end())
				{
					(*rndCluster)[vIter->first]=vIter->second;
				}
			}
		}
		else
		{
			int csize=trueCluster->size();
			while(rndCluster->size()<trueCluster->size())
			{
				int rndIndex=rand.getRandomNumber();
				while(usedVariables.find(rndIndex)!=usedVariables.end())
				{
					rndIndex=rand.getRandomNumber();
				}
				Variable* aVar=variableSet[rndIndex];
				(*rndCluster)[rndIndex]=aVar;
				oFile <<" " << rndIndex;
				usedVariables[rndIndex]=0;
			}
			oFile << endl;
		}
	}
	oFile.close();
	return randomVariableSubsets;
}

int 
ClusterManager::readClusters()
{
	char aFName[1024];
	sprintf(aFName,"%s/km_clusters_%d.txt",outputDir,clusterCnt);
	ifstream inFile(aFName);
	string strbuffer;
	char* buffer=NULL;
	unsigned int currBuffSize=0;
	VSET* varCluster=NULL;
	VSET& variableSet=vMgr->getVariableSet();
	while(inFile.good())
	{
		getline(inFile,strbuffer);
		if(strbuffer.length()<=0)
		{
			continue;
		}
		if(currBuffSize<(strbuffer.length()+1))
		{
			if(buffer!=NULL)
			{
				delete[] buffer;
			}
			buffer=new char[strbuffer.length()+1];
			currBuffSize=strbuffer.length()+1;
		}
		strcpy(buffer,strbuffer.c_str());
		if(strstr(buffer,"ClusterId")!=NULL)
		{
			varCluster=new VSET;
			variableSubsets.push_back(varCluster);
		}
		else
		{
			char* tok=strtok(buffer," ");
			while(tok!=NULL)
			{
				int vId=atoi(tok);
				(*varCluster)[vId]=variableSet[vId];
				tok=strtok(NULL," ");
			}
		}
	}
	inFile.close();
	return 0;
}

int 
ClusterManager::generateClusters()
{
	map<int,INTDBLMAP*> kmeansData;
	//evidMgr->generateMatrixFormat(kmeansData);
	kmeans.setConvergenceThreshold(0.01);
	kmeans.setClusterCnt(clusterCnt);
	kmeans.setDistanceType(Kmeans::MI);
	kmeans.setMaxIter(maxIter);
	kmeans.cluster(&kmeansData);
	char aFName[1024];
	sprintf(aFName,"%s/km_clusters_%d.txt",outputDir,clusterCnt);
	kmeans.showClusters(aFName);
	CLUSTERS& clusters=kmeans.getClusters();
	VSET& variableSet=vMgr->getVariableSet();
	for(map<int,MEMBER*>::iterator cIter=clusters.begin();cIter!=clusters.end();cIter++)
	{
		MEMBER* members=cIter->second;
		VSET* varCluster=new VSET;
		variableSubsets.push_back(varCluster);
		for(unsigned int i=0;i<members->size();i++)
		{
			int vId=(*members)[i];
			(*varCluster)[vId]=variableSet[vId];
		}
	}

	return 0;
}

int 
ClusterManager::generateClusters_Randomrestarts()
{
	map<int,INTDBLMAP*> kmeansData;
	//evidMgr->generateMatrixFormat(kmeansData);
	vector<Kmeans*> kmset;
	double minDist=-1;
	double maxDist=-1;
	int bestkm_min=0;
	int bestkm_max=0;
	int restartCnt=10;
	char rstFName[1024];
	sprintf(rstFName,"%s/rand_restart.txt",outputDir);
	ofstream oFile(rstFName);
	for(int i=0;i<restartCnt;i++)
	{
		Kmeans* k=new Kmeans;
		k->setConvergenceThreshold(0.01);
		k->setClusterCnt(clusterCnt);
		k->setDistanceType(Kmeans::MI);
		k->setMaxIter(maxIter);
		k->cluster(&kmeansData);
		double dist=k->getWithinClusterDist();
		if(minDist==-1)
		{
			minDist=dist;
			maxDist=dist;
			bestkm_min=kmset.size();
			bestkm_max=kmset.size();
		}
		else 
		{
			if(minDist>dist)
			{
				minDist=dist;
				bestkm_min=kmset.size();
			}
			if(maxDist<dist)
			{
				maxDist=dist;
				bestkm_max=kmset.size();
			}
		}
		kmset.push_back(k);
		oFile <<i <<"\t" << dist << endl;
		sprintf(rstFName,"%s/restarts/clusters_%d_%d.txt",outputDir,clusterCnt,i);
		k->showClusters(rstFName);
	}
	Kmeans* bestkmeans=kmset[bestkm_max];
	cout << "Best kmeans found at " << bestkm_max <<" with dist " << maxDist << endl;
	char aFName[1024];
	sprintf(aFName,"%s/km_clusters_%d.txt",outputDir,clusterCnt);
	bestkmeans->showClusters(aFName);
	CLUSTERS& clusters=bestkmeans->getClusters();
	VSET& variableSet=vMgr->getVariableSet();
	for(map<int,MEMBER*>::iterator cIter=clusters.begin();cIter!=clusters.end();cIter++)
	{
		MEMBER* members=cIter->second;
		VSET* varCluster=new VSET;
		variableSubsets.push_back(varCluster);
		for(unsigned int i=0;i<members->size();i++)
		{
			int vId=(*members)[i];
			(*varCluster)[vId]=variableSet[vId];
		}
	}
	for(int i=0;i<kmset.size();i++)
	{
		Kmeans* k=kmset[i];
		k->clear();
		delete k;
	}
	kmset.clear();
	return 0;
}

int 
ClusterManager::generateClusters_Randomrestarts(map<int,INTDBLMAP*>& kmeansData)
{
	vector<Kmeans*> kmset;
	double minDist=-1;
	double maxDist=-1;
	int bestkm_min=0;
	int bestkm_max=0;
	int restartCnt=10;
	char rstFName[1024];
	sprintf(rstFName,"%s/rand_restart.txt",outputDir);
	ofstream oFile(rstFName);
	for(int i=0;i<restartCnt;i++)
	{
		Kmeans* k=new Kmeans;
		k->setConvergenceThreshold(0.01);
		k->setClusterCnt(clusterCnt);
		k->setDistanceType(Kmeans::MI);
		k->setMaxIter(maxIter);
		k->cluster(&kmeansData);
		double dist=k->getWithinClusterDist();
		if(minDist==-1)
		{
			minDist=dist;
			maxDist=dist;
			bestkm_min=kmset.size();
			bestkm_max=kmset.size();
		}
		else 
		{
			if(minDist>dist)
			{
				minDist=dist;
				bestkm_min=kmset.size();
			}
			if(maxDist<dist)
			{
				maxDist=dist;
				bestkm_max=kmset.size();
			}
		}
		kmset.push_back(k);
		oFile <<i <<"\t" << dist << endl;
		sprintf(rstFName,"%s/restarts/clusters_%d_%d.txt",outputDir,clusterCnt,i);
		k->showClusters(rstFName);
	}
	Kmeans* bestkmeans=kmset[bestkm_max];
	cout << "Best kmeans found at " << bestkm_max <<" with dist " << maxDist << endl;
	char aFName[1024];
	sprintf(aFName,"%s/km_clusters_%d.txt",outputDir,clusterCnt);
	bestkmeans->showClusters(aFName);
	CLUSTERS& clusters=bestkmeans->getClusters();
	for(map<int,MEMBER*>::iterator cIter=clusters.begin();cIter!=clusters.end();cIter++)
	{
		MEMBER* members=cIter->second;
		INTINTMAP* varCluster=new INTINTMAP;
		variableSubsets_IDs.push_back(varCluster);
		for(unsigned int i=0;i<members->size();i++)
		{
			int vId=(*members)[i];
			(*varCluster)[vId]=0;
		}
	}
	for(int i=0;i<kmset.size();i++)
	{
		Kmeans* k=kmset[i];
		k->clear();
		delete k;
	}
	kmset.clear();
	return 0;
}
