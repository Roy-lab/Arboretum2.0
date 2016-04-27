/*
Arboretum: An algorithm to cluster functional genomesomics data from multiple species
    Copyright (C) 2013 Sushmita Roy sushroy@gmail.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <iostream>
#include <fstream>
#include <iostream>
#include <string.h>
#include "Matrix.H"
#include "SpeciesDistManager.H"

SpeciesDistManager::SpeciesDistManager()
{
	root=NULL;
	maxClusterCnt=0;
}

SpeciesDistManager::~SpeciesDistManager()
{
}

int
SpeciesDistManager::setMaxClusters(int k)
{
	maxClusterCnt=k;
	return 0;
}

int 
SpeciesDistManager::readSpeciesTree(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		if(strchr(buffer,'#')!=NULL)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string childSpeciesName;
		string parentSpeciesName;
		string childtype;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				childSpeciesName.append(tok);	
			}
			else if(tokCnt==1)
			{
				childtype.append(tok);
			}
			else if(tokCnt==2)
			{
				parentSpeciesName.append(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		SpeciesDistManager::Species* childSpecies=NULL;
		SpeciesDistManager::Species* parentSpecies=NULL;
		if(speciesSet.find(childSpeciesName)==speciesSet.end())
		{
			childSpecies=new SpeciesDistManager::Species;
			childSpecies->name.append(childSpeciesName.c_str());
			childSpecies->parent=NULL;
			childSpecies->leftchild=NULL;
			childSpecies->rightchild=NULL;
			speciesSet[childSpeciesName]=childSpecies;
			childSpecies->conditional=new Matrix(maxClusterCnt,maxClusterCnt);
		}
		else
		{
			childSpecies=speciesSet[childSpeciesName];
			if(childSpecies->conditional==NULL)
			{
				childSpecies->conditional=new Matrix(maxClusterCnt,maxClusterCnt);
			}
		}
		if(speciesSet.find(parentSpeciesName)==speciesSet.end())
		{
			parentSpecies=new SpeciesDistManager::Species;
			parentSpecies->name.append(parentSpeciesName.c_str());
			parentSpecies->parent=NULL;
			parentSpecies->leftchild=NULL;
			parentSpecies->rightchild=NULL;
			speciesSet[parentSpeciesName]=parentSpecies;
		}
		else
		{
			parentSpecies=speciesSet[parentSpeciesName];
		}
		childSpecies->parent=parentSpecies;
		if(strcmp(childtype.c_str(),"left")==0)
		{
			parentSpecies->leftchild=childSpecies;
		}	
		else if(strcmp(childtype.c_str(),"right")==0)
		{
			parentSpecies->rightchild=childSpecies;
		}
		else
		{
			cout <<"Unknown edge type " << childtype.c_str() << endl;
		}
		if(parentSpecies->parent==NULL)
		{
			cout << "Root set to " << parentSpecies->name << endl;
			root=parentSpecies;
		}
	}
	root->conditional=new Matrix(1,maxClusterCnt);

	inFile.close();
	return 0;
}


int
SpeciesDistManager::assignLevel()
{
	levelFromRoot[root->name]=0;
	assignLevel(root->leftchild,1);
	assignLevel(root->rightchild,1);
	return 0;
}

int
SpeciesDistManager::assignLevel(SpeciesDistManager::Species* node,int level)
{
	levelFromRoot[node->name]=level;
	if(node->leftchild!=NULL)
	{
		assignLevel(node->leftchild,level+1);
	}
	if(node->rightchild!=NULL)
	{
		assignLevel(node->rightchild,level+1);
	}
	return 0;
}


int
SpeciesDistManager::getLevelFromRoot(const char* nodeName)
{
	string nKey(nodeName);
	if(levelFromRoot.find(nKey)==levelFromRoot.end())
	{	
		cout <<"No node with name " << nKey << endl;
		return -1;
	}
	return levelFromRoot[nKey];
}


SpeciesDistManager::Species*
SpeciesDistManager::getRoot()
{
	return root;
}

int 
SpeciesDistManager::getSpeciesListPrefix(vector<string>& specList)
{
	getSpeciesListPrefix(root,specList);
	return 0;
}


map<string,SpeciesDistManager::Species*>& 
SpeciesDistManager::getAllSpecies()
{	
	return speciesSet;
}

SpeciesDistManager::Species*
SpeciesDistManager::getSpecies(string& specKey)
{
	if(speciesSet.find(specKey)==speciesSet.end())
	{
		return NULL;
	}
	return speciesSet[specKey];
}

int
SpeciesDistManager::getSpeciesListPrefix(SpeciesDistManager::Species* node, vector<string>& specList)
{
	if(node->leftchild!=NULL)
	{
		getSpeciesListPrefix(node->leftchild,specList);
	}
	if(node->rightchild!=NULL)
	{	
		getSpeciesListPrefix(node->rightchild,specList);
	}
	specList.push_back(node->name);
	return 0;
}

int 
SpeciesDistManager::resetTransitionProbability()
{
	//root->conditional->setAllValues(0.001);
	root->conditional->setAllValues(0);
	resetTransitionProbability(root->leftchild);
	resetTransitionProbability(root->rightchild);
	return 0;
}

int
SpeciesDistManager::resetTransitionProbability(SpeciesDistManager::Species* species)
{
	species->conditional->setAllValues(0.001);
	if(species->leftchild!=NULL)
	{	
		resetTransitionProbability(species->leftchild);
	}
	if(species->rightchild!=NULL)
	{
		resetTransitionProbability(species->rightchild);
	}
	return 0;
}

int 
SpeciesDistManager::normalizeTransitionMatrix()
{
	normalizeTransitionMatrix(root);
	return 0;
}

int
SpeciesDistManager::normalizeTransitionMatrix(SpeciesDistManager::Species* node)
{
	Matrix* conditional=node->conditional;
	double globalsum=0;
	cout <<"Transition for " << node->name << " before norm" << endl;
	node->conditional->showMatrix();
	for(int r=0;r<conditional->getRowCnt();r++)
	{
		double sum=0;
		for(int c=0;c<conditional->getColCnt();c++)
		{
			double aval=conditional->getValue(r,c);
			sum=sum+aval;
		}
		globalsum=globalsum+sum;
		for(int c=0;c<conditional->getColCnt();c++)
		{
			double prob=conditional->getValue(r,c);
			prob=prob/sum;
			conditional->setValue(prob,r,c);
		}
	}
	cout <<"Transition for " << node->name << endl;
	node->conditional->showMatrix();
	if(node->leftchild!=NULL)
	{
		normalizeTransitionMatrix(node->leftchild);
	}
	if(node->rightchild!=NULL)
	{
		normalizeTransitionMatrix(node->rightchild);
	}
	return 0;
}

int
SpeciesDistManager::initTransitionMatrix_ML()
{
	//Start at the root and move downwards
	root->conditional_ml=new Matrix(1,maxClusterCnt);
	root->conditional_ml->setAllValues(0.001);
	if(root->leftchild!=NULL)
	{
		initTransitionMatrix_ML(root->leftchild);
	}
	if(root->rightchild!=NULL)
	{
		initTransitionMatrix_ML(root->rightchild);
	}
	return 0;
}

int
SpeciesDistManager::initTransitionMatrix_ML(SpeciesDistManager::Species* node)
{
	node->conditional_ml=new Matrix(maxClusterCnt,maxClusterCnt);
	node->conditional_ml->setAllValues(0.001);
	if(node->leftchild!=NULL)
	{
		initTransitionMatrix_ML(node->leftchild);
	}
	if(node->rightchild!=NULL)
	{	
		initTransitionMatrix_ML(node->rightchild);
	}
	return 0;
}

int 
SpeciesDistManager::normalizeTransitionMatrix_ML()
{
	normalizeTransitionMatrix_ML(root);
	return 0;
}

int
SpeciesDistManager::normalizeTransitionMatrix_ML(SpeciesDistManager::Species* node)
{
	Matrix* conditional=node->conditional_ml;
	double globalsum=0;
	cout <<"Transition for " << node->name << " before norm" << endl;
	conditional->showMatrix();
	for(int r=0;r<conditional->getRowCnt();r++)
	{
		double sum=0;
		for(int c=0;c<conditional->getColCnt();c++)
		{
			double aval=conditional->getValue(r,c);
			sum=sum+aval;
		}
		globalsum=globalsum+sum;
		for(int c=0;c<conditional->getColCnt();c++)
		{
			double prob=conditional->getValue(r,c);
			prob=prob/sum;
			conditional->setValue(prob,r,c);
		}
	}
	cout <<"Transition for " << node->name << endl;
	conditional->showMatrix();
	if(node->leftchild!=NULL)
	{
		normalizeTransitionMatrix_ML(node->leftchild);
	}
	if(node->rightchild!=NULL)
	{
		normalizeTransitionMatrix_ML(node->rightchild);
	}
	return 0;
}


Matrix*
SpeciesDistManager::getTransitionMatrix_ML(string& spName)
{
	if(speciesSet.find(spName)==speciesSet.end())
	{
		cout <<"No species with name " << spName.c_str() << endl;
		exit(0);
	}
	Species* species=speciesSet[spName];
	return species->conditional_ml;
}




Matrix*
SpeciesDistManager::getConditional(string& spName)
{
	if(speciesSet.find(spName)==speciesSet.end())
	{
		cout <<"No species with name " << spName.c_str() << endl;
		exit(0);
	}
	Species* species=speciesSet[spName];
	return species->conditional;
}

//The probability that an edge is maintained in a child given that the edge is present in the ancestor
double 
SpeciesDistManager::getConditionalProb(string& spName, int parentCluster,int childCluster)
{
	if(speciesSet.find(spName)==speciesSet.end())
	{
		cout <<"No species with name " << spName.c_str() << endl;
		exit(0);
	}
	Species* species=speciesSet[spName];
	double prob=species->conditional->getValue(parentCluster,childCluster);
	return prob;
}


double
SpeciesDistManager::getEdgeStatusProb(map<string,int>& edgeStatus)
{
	/*for(map<string,int>::iterator eIter=edgeStatus.begin();eIter!=edgeStatus.end();eIter++)
	{
		cout <<eIter->first <<"\t" << eIter->second << endl;
	}*/
	double edgeprior=0;
	//Start with the root
	for(int i=0;i<maxClusterCnt;i++)
	{
		double leftScore=getSubTree(i,root->leftchild,edgeStatus);
		double rightScore=getSubTree(i,root->rightchild,edgeStatus);
		double score=0.5*leftScore*rightScore;
		edgeprior=edgeprior+score;
	}
	return edgeprior;
}


int 
SpeciesDistManager::showInferredConditionals(const char* outputDir)
{
	showConditionals(outputDir,root);
	return 0;
}

int
SpeciesDistManager::showConditionals(const char* outputDir,SpeciesDistManager::Species* species)
{
	char output[1024];
	sprintf(output,"%s/%s",outputDir,species->name.c_str());
	ofstream oFile(output);	
	species->conditional->showMatrix(oFile);
	if(species->leftchild!=NULL)
	{
		showConditionals(outputDir,species->leftchild);
		showConditionals(outputDir,species->rightchild);
	}
	return 0;
}

int 
SpeciesDistManager::showInferredConditionals_ML(const char* outputDir)
{
	showConditionals_ML(outputDir,root);
	return 0;
}

int
SpeciesDistManager::showConditionals_ML(const char* outputDir,SpeciesDistManager::Species* species)
{
	char output[1024];
	sprintf(output,"%s/ml_%s",outputDir,species->name.c_str());
	ofstream oFile(output);	
	species->conditional_ml->showMatrix(oFile);
	if(species->leftchild!=NULL)
	{
		showConditionals_ML(outputDir,species->leftchild);
		showConditionals_ML(outputDir,species->rightchild);
	}
	return 0;
}


double
SpeciesDistManager::getSubTree(int parentCluster, Species* child, map<string,int>& edgeStatus)
{
	double score=0;
	if(child->leftchild==NULL && child->rightchild==NULL)
	{
		//This is a leaf node
		int status=edgeStatus[child->name];
		score=getConditionalProb(child->name,parentCluster,status);
	}
	else
	{
		//This is an interior node
		for(int i=0;i<maxClusterCnt;i++)
		{
			double leftScore=getSubTree(i,child->leftchild,edgeStatus);
			double rightScore=getSubTree(i,child->rightchild,edgeStatus);
			double prob=getConditionalProb(child->name,parentCluster,edgeStatus[child->name]);
			score=score+(prob*leftScore*rightScore);
		}
	}
	return score;
}

double
SpeciesDistManager::getAncestralClustering(map<string,int>& extantClusterAssign,map<string,int>& ancestralClusterAssign)
{
	//Start with the root
	int maxParentCluster=-1;
	double maxScore=0;
	map<int,double> assignProb;
	map<string,int> tempAssign;
	for(int i=0;i<maxClusterCnt;i++)
	{
		double leftScore=maxSubTree(i,root->leftchild,extantClusterAssign,tempAssign);
		double rightScore=maxSubTree(i,root->rightchild,extantClusterAssign,tempAssign);
		double prior=root->conditional->getValue(i,0);
		double score=prior*leftScore*rightScore;
		assignProb[i]=score;
		if(score>maxScore)	
		{
			maxScore=score;
			maxParentCluster=i;
			for(map<string,int>::iterator aIter=tempAssign.begin();aIter!=tempAssign.end();aIter++)
			{
				ancestralClusterAssign[aIter->first]=aIter->second;
			}
		}
		/*cout << score <<" " << maxScore;
		for(map<string,int>::iterator aIter=ancestralClusterAssign.begin();aIter!=ancestralClusterAssign.end();aIter++)
		{
			cout <<" " << aIter->first<<"=" << aIter->second;
		}
		cout << endl;*/
	}
	tempAssign.clear();
	for(map<int,double>::iterator aIter=assignProb.begin();aIter!=assignProb.end();aIter++)
	{
		if(aIter->second>=maxScore && aIter->first!=maxParentCluster)
		{
			cout<< root->name<<" " <<aIter->first <<"=" << aIter->second << endl;
		}
	}
	ancestralClusterAssign[root->name]=maxParentCluster;
	return maxScore;
}

double 
SpeciesDistManager::getExtantClustering(map<string,int>& ancestralClusterAssign, map<string,int>& extantClusterAssign)
{
	int maxParentCluster=-1;
	assignExtantClustering(ancestralClusterAssign,root,extantClusterAssign);
	return 0;
}

int
SpeciesDistManager::assignExtantClustering(map<string,int>& ancestralClusterAssign,Species* node, map<string,int>& extantClusterAssign)
{
	if(node->leftchild==NULL)
	{
		Species* parent=node->parent;
		int clusterid=ancestralClusterAssign[parent->name];
		int maxchildclusterid=getMaxClusterAssignForChild(clusterid,node);
		extantClusterAssign[node->name]=maxchildclusterid;
	}
	else
	{
		assignExtantClustering(ancestralClusterAssign,node->leftchild,extantClusterAssign);
		assignExtantClustering(ancestralClusterAssign,node->rightchild,extantClusterAssign);
	}
	return 0;
}

int
SpeciesDistManager::getMaxClusterAssignForChild(int parentClustID,Species* child)
{
	double maxchildprob=0;
	int maxchildcluster=-1;
	for(int i=0;i<maxClusterCnt;i++)
	{
		double prob=getConditionalProb(child->name,parentClustID,i);
		if(prob>maxchildprob)
		{
			maxchildprob=prob;
			maxchildcluster=i;
		}
	}
	return maxchildcluster;
}

/*
int
SpeciesDistManager::enumerateAllAncestors(map<string,int>& edgeStatus)
{
	bool edgeInParent[]={false,true};
	map<string,int> jointAssign;
	for(map<string,int>::iterator eIter=edgeStatus.begin();eIter!=edgeStatus.end();eIter++)
	{
		jointAssign[eIter->first]=eIter->second;
	}
	vector<map<string,int>*> assignSet;
	map<string,int> path;
	for(int i=0;i<2;i++)
	{
		path[root->name]=edgeInParent[i];
		enumerateChild(root->leftchild,path,edgeStatus,assignSet);			
		enumerateChild(root->rightchild,path,edgeStatus,assignSet);			
		for(int i=0;i<assignSet.size();i++)
		{
			map<string,int>* assignment=assignSet[i];
			if(assignment->size()<(edgeStatus.size() - 1))
			{
				continue;
			}
			for(map<string,int>::iterator cIter=edgeStatus.begin();cIter!=edgeStatus.end();cIter++)
			{
				(*assignment)[cIter->first]=cIter->second;
			}
			double ass=scoreAssignment(*assignment);
			cout << ass;
			for(map<string,int>::iterator aIter=assignment->begin();aIter!=assignment->end();aIter++)
			{
				cout << " " << aIter->first<<"=" << aIter->second;
			}
			cout << endl;
			assignment->clear();
			delete assignment;
		}
		assignSet.clear();
	}
	return 0;
}

int
SpeciesDistManager::enumerateChild(Species* child,map<string,int>& path, map<string,int>& leafAssign, vector<map<string,int>*>& assignmentSet)
{
	bool edgeInChild[]={false,true};
	//For every parent configuration get the children configuration
	if(child->leftchild==NULL)
	{
		//cout << child->name << "=" << leafAssign[child->name];
		//map<string,int>* assignment=new map<string,int>;
		//assignmentSet.push_back(assignment);
	//	path[child->name]=leafAssign[child->name];
	}
	else
	{
		for(int i=0;i<2;i++)
		{
			//cout <<child->name <<"=" << edgeInChild[i];
			map<string,int>* assignment=new map<string,int>;
			for(map<string,int>::iterator aIter=path.begin();aIter!=path.end();aIter++)
			{
				(*assignment)[aIter->first]=aIter->second;
			}
			(*assignment)[child->name]=edgeInChild[i];
			assignmentSet.push_back(assignment);
			enumerateChild(child->leftchild,path,leafAssign,assignmentSet);
			enumerateChild(child->rightchild,path,leafAssign,assignmentSet);
		}
	}
	return 0;
}*/

double
SpeciesDistManager::scoreAssignment(map<string,int>& jointAssign)
{
	double assignpdf=1;
	for(map<string,int>::iterator aIter=jointAssign.begin();aIter!=jointAssign.end();aIter++)
	{
		Species* node=speciesSet[aIter->first];
		int clusterid=aIter->second;
		if(node->parent==NULL)
		{
			assignpdf=assignpdf*root->conditional->getValue(clusterid,0);;
		}
		else
		{
			Species* parent=node->parent;
			if(jointAssign.find(parent->name)==jointAssign.end())
			{
				cout <<"No assignment for " << parent->name << endl;
			}
			int parentAssign=jointAssign[parent->name];
			double score=getConditionalProb((string&)aIter->first,parentAssign,aIter->second);
			assignpdf=assignpdf*score;
		}
	}
	return assignpdf;
}


double
SpeciesDistManager::maxSubTree(int parentCluster, Species* child, map<string,int>& extantCluster,map<string,int>& ancestralCluster)
{
	double score=0;
	if(child->leftchild==NULL && child->rightchild==NULL)
	{
		//This is a leaf node
		int clusterassign=extantCluster[child->name];
		score=getConditionalProb(child->name,parentCluster,clusterassign);
	}
	else
	{
		//This is an interior node
		int maxParentCluster=-1;
		double maxScore=0;
		map<int,double> assignProb;
		for(int i=0;i<maxClusterCnt;i++)
		{
			double leftScore=maxSubTree(i,child->leftchild,extantCluster,ancestralCluster);
			double rightScore=maxSubTree(i,child->rightchild,extantCluster,ancestralCluster);
			double prob=getConditionalProb(child->name,parentCluster,i);
			double currScore=prob*leftScore*rightScore;
			if(currScore>maxScore)
			{
				maxScore=currScore;
				maxParentCluster=i;
			}
			assignProb[i]=currScore;
		}
		for(map<int,double>::iterator aIter=assignProb.begin();aIter!=assignProb.end();aIter++)
		{
			/*if((aIter->second>=maxScore) && (aIter->first!=maxParentCluster))
			{
				cout<< "Found clash score: " << aIter->second << " currmax=" << maxParentCluster << " parent:"  <<child->parent->name << "=" << parentCluster 
				<< " current: " << child->name<<"="<< aIter->first;
				if(child->leftchild->leftchild==NULL)
				{
					cout << "left: "<< child->leftchild->name<<"="<< extantCluster[child->leftchild->name];
					cout << " right: "<< child->rightchild->name<<"="<< extantCluster[child->rightchild->name] << endl;
				}
				else
				{
					cout << "left: "<< child->leftchild->name<<"="<< ancestralCluster[child->leftchild->name];
					cout << " right: "<< child->rightchild->name<<"="<< ancestralCluster[child->rightchild->name] << endl;
				}
			}*/
		}
		score=maxScore;
		if(strcmp(child->name.c_str(),"Anc3")==0)
		{
	//		cout <<"Update: " << child->parent->name<<"="<< parentCluster<<" "<< child->name<<"=" << maxParentCluster <<" Score=" <<maxScore << endl;
		}
		ancestralCluster[child->name]=maxParentCluster;
		assignProb.clear();
	}
	return score;
}

/*
int 
SpeciesDistManager::estimateNormalizationConstants(map<string,map<int,double>*>& dataPtProb, map<string,double>& normConstants)
{
	double nodeNormConstant=0;
	Matrix* conditional=root->conditional;
	int clusterCnt=conditional->getRowCnt();
	for(int cId=0;cId<clusterCnt;cId++)
	{
		double leftNorm=estimateNormalizationConstant_Node(dataPtProb,normConstants,root->leftchild,cId);
		double rightNorm=estimateNormalizationConstant_Node(dataPtProb,normConstants,root->rightchild,cId);
		double prior=condition->getValueAt(cId,0);
		double score=prior*leftNorm*rightNorm;
		nodeNormConstant=nodeNormConstant+score;
	}
	normConstants[root->name]=nodeNormConstant;
	return nodeNormConstant;
}

int 
SpeciesDistManager::estimateNormalizationConstant_Node(map<string,map<int,double>*>& dataPtProb, map<string,double>& normConstants,SpeciesDistManager::Species* node,int parentClustID)
{
	double nodeNormConstant=0;
	Matrix* params=node->conditional;
	int clusterCnt=params->getClusterCnt();
	if(node->leftchild==NULL)
	{
		if(dataPtProb.find(node->name)==dataPtProb.end())
		{
			return nodeNormCosntant;
		}
		map<int,double>* obsProb=dataPtProb[node->name];
		for(int i=0;i<clusterCnt;i++)
		{
			double prior=params->getValueAt(parentClustID,i);
			double dprob=0;
			if(obsProb->find(i)==obsProb->ebd())
			{
				dprob=(*obsProb)[i];
			}
			nodeNormConstant=nodeNormConstant+(dprob*prior);
		}
	}
	else
	{
		for(int i=0;i<clusterCnt;i++)
		{
			double prior=params->getValueAt(parentClustID,i);
			double leftChildValue=estimateNormalizationConstant_Node(dataPtProb,nodeNormConstant,node->leftchild,i);
			double rightChildValue=estimateNormalizationConstant_Node(dataPtProb,nodeNormConstant,node->leftchild,i);
			nodeNormConstant=nodeNormConstant+(prior*leftChildValue*rightChildValue);
		}
	}
	return nodeNormConstant;
}*/
