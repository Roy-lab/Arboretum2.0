
#include <math.h>
#include <string.h>
#include <gsl/gsl_randist.h>
#include "SpeciesDistManager.H"
#include "GeneTree.H"
#include "GeneTreeManager.H"
#include "Matrix.H"
#include "Gamma.H"
#include "CommonTypes.H"
#include "MappedOrthogroup.H"
#include "MappedOrthogroupReader.H"
#include "GammaManager.H"


GammaManager::GammaManager()
{
}

GammaManager::~GammaManager()
{
}

int
GammaManager::showTotalUpdatedParents()
{
	gtMgr.showTotalAdditionalParentDist();
	return 0;
}


int 
GammaManager::setMaxClusterCnt(int k)
{
	maxClusterCnt=k;
	return 0;
}

int 
GammaManager::setOrthogroupReader(MappedOrthogroupReader* aPtr)
{
	mgrpReader=aPtr;
	return 0;
}

int 
GammaManager::setSpeciesDistManager(SpeciesDistManager* aPtr)
{
	spdistMgr=aPtr;
	gtMgr.setSpeciesTree(spdistMgr->getRoot());
	gtMgr.setSpeciesDistManager(spdistMgr);
	return 0;
}


int 
GammaManager::initGamma(int ogid, string& geneName, string& specName,int clustID)
{
	if(ogid==1002|| ogid==13692)
	{
		cout <<"Found OGID "<< ogid << endl;
	}
	if(strcmp(geneName.c_str(),"CAGL0B03899g")==0)
	{
		cout <<"Stop here " << endl;
	}
	Gamma* gamma=NULL;
	MappedOrthogroup* mor=mgrpReader->getMappedOrthogroup(geneName.c_str(),specName.c_str());
	if(gammaSet.find(ogid)==gammaSet.end())
	{
		gamma=new Gamma;
		gamma->setMaxClusterCnt(maxClusterCnt);
		gammaSet[ogid]=gamma;	
		GeneTree* gtree=gtMgr.getGeneTree(mor);
		gamma->initUsingGeneTree(gtree);
		if(gamma->getDupAncestor().length()>0)
		{
			ogDupAncMap[ogid]=gamma->getDupAncestor();
		}
		if(ogid==4597)
		{
			cout <<"Now showing " << ogid << endl;
			gamma->showTree();
		}
	}
	else
	{
		gamma=gammaSet[ogid];
	}
	//gamma->addGene(geneName,specName,maxClusterCnt);
	gamma->setGeneClusterID(geneName,specName,clustID);
	return 0;
}

int 
GammaManager::estimateLeafGamma(int ogid,map<int,double>& prob,string& geneName, string& specName)
{
	if(ogid==4595 && strcmp(specName.c_str(),"Scer")==0)
	{
	//	cout << "Stop here" << endl;	
	}
	
	Gamma* gamma_og=gammaSet[ogid];
	Gamma::Node* m=gamma_og->getGeneNode(geneName,specName);
	if(m->normTerm==NULL)
	{
		m->normTerm=new Matrix(1,maxClusterCnt);
	}
	Matrix* conditional=spdistMgr->getConditional(specName);
	for(int r=0;r<m->gamma->getRowCnt();r++)
	{
		double anormTerm=0;
		for(map<int,double>::iterator pIter=prob.begin();pIter!=prob.end();pIter++)
		{
			double pval=pIter->second;
			double prior=conditional->getValue(r,pIter->first);
			double post=pval*prior;
			if(isnan(post) || isinf(post))
			{
				cout <<"Posterior is nan/inf for  "<< geneName <<" species " << specName << "cluster id "  << pIter->first
				<< " pval=" << pval << " prior=" << prior << endl;
			}
			anormTerm=anormTerm+post;
			m->gamma->setValue(post,r,pIter->first);
		}
		for(int c=0;c<m->gamma->getColCnt();c++)
		{
			double aval=m->gamma->getValue(r,c);
			aval=aval/anormTerm;
			m->gamma->setValue(aval,r,c);
		}
		m->normTerm->setValue(anormTerm,0,r);
	}
	/*if(ogid==1585|| ogid==2178|| ogid==2894)
	{
		cout << "Leaf gamma: " <<specName <<" " << m->name << endl;
		m->gamma->showMatrix();
	}*/
	return 0;
}

//At this point all the leaf nodes have been initialized with whatever the data told us
int 
GammaManager::estimateNonLeafPosterior()
{
	for(map<int,Gamma*>::iterator gIter=gammaSet.begin();gIter!=gammaSet.end();gIter++)
	{
		Gamma* gamma_og=gIter->second;
		showGammas=false;
		if(gIter->first==20057|| gIter->first==13692)
		{
			//showGammas=true;
			cout <<"Stop here" << endl;
		}
		estimateNonLeafPosterior(gamma_og->getRoot());
	}
	return 0;
}

int
GammaManager::estimateNonLeafPosterior(Gamma::Node* node)
{	
	if(node==NULL)
	{
		return 0;
	}
	if(node->leftchild==NULL && node->rightchild==NULL)
	{
		if(showGammas)
		{
			cout << "Gamma for " << node->name << endl;
			node->gamma->showMatrix();
		}
		return 0;
	}
	else
	{
		estimateNonLeafPosterior(node->leftchild);
		estimateNonLeafPosterior(node->rightchild);
		Matrix* conditional=spdistMgr->getConditional(node->species);
		int rowcnt=conditional->getRowCnt();
		if(node->normTerm==NULL)
		{
			node->normTerm=new Matrix(1,rowcnt);
		}
		for(int r=0;r<node->gamma->getRowCnt();r++)
		{
			double anormTerm=0;
			for(int c=0;c<node->gamma->getColCnt();c++)
			{
				double leftval=1;
				if(node->leftchild!=NULL && node->leftchild->normTerm!=NULL)
				{
					leftval=node->leftchild->normTerm->getValue(0,c);	
				}
				double rightval=1;	
				if(node->rightchild!=NULL && node->rightchild->normTerm!=NULL)
				{
					rightval=node->rightchild->normTerm->getValue(0,c);
				}
				double unnormVal=leftval*rightval;
				if(node->nodeType==1)
				{
					double prior=conditional->getValue(r,c);
					unnormVal=unnormVal*prior;
					if(unnormVal<1e-180)
					{
						unnormVal=1e-180;
					}
					if(isnan(unnormVal) || isinf(unnormVal))
					{
						cout <<"Posterior is nan/inf for  "<< node->name <<" species " << node->species<< " cluster id "  << c
						<< " pval=" << unnormVal << " prior=" << prior << endl;
					}
				}
				node->gamma->setValue(unnormVal,r,c);
				anormTerm=anormTerm+unnormVal;
			}
			if(anormTerm==0)
			{
				cout <<"normTerm is 0 for  "<< node->name <<" species " << node->species<< " cluster id "  << r << endl;
			}
			node->normTerm->setValue(anormTerm,0,r);
			for(int c=0;c<node->gamma->getColCnt();c++)
			{
				double post=node->gamma->getValue(r,c);
				node->gamma->setValue(post/anormTerm,r,c);
			}
		}
		if(showGammas)
		{
			cout << "Gamma for " << node->name << endl;
			node->gamma->showMatrix();
		}
	}
	return 0;
}


map<int,double>*
GammaManager::getLeafLikelihood_store(int ogid,string& geneName)
{
	Gamma* gamma_og=gammaSet[ogid];
	map<int,double>* llstore=gamma_og->getDataLL(geneName);
	return llstore;
}


double
GammaManager::getAllNodeScore()
{
	double totalScore=0;
	for(map<int,Gamma*>::iterator gIter=gammaSet.begin();gIter!=gammaSet.end();gIter++)
	{
		double score=0;
		score=gIter->second->root->normTerm->getValue(0,0);
		totalScore=totalScore+log(score);
		score=getNodeScore(gIter->second->root);
	//	totalScore=totalScore+score;
	}
	return totalScore;
}

double
GammaManager::getNodeScore(Gamma::Node* anode)
{
	double nodeContrib=0;
	if(anode==NULL)
	{
		return 0;
	}
	if(anode->nodeType==1)
	{
		if(anode->leftchild==NULL && anode->rightchild==NULL)
		{
		//	cout << "Reached a leaf node " << endl;
		}
		Matrix* gamma=anode->gamma;
		Matrix* conditional=spdistMgr->getConditional(anode->species);
		for(int r=0;r<gamma->getRowCnt();r++)
		{
			for(int c=0;c<gamma->getColCnt();c++)
			{
				double gval=gamma->getValue(r,c);
				double prior=conditional->getValue(r,c);
				double pval=0; 
				if(anode->hasData)
				{	
					map<int,double>* dataLL=anode->dataLL;
					double dpll=(*dataLL)[c];
					pval=gval*(log(dpll)+log(prior));
				}
				else
				{
					pval=gval*log(prior);
				}
				nodeContrib=nodeContrib+pval;
			}
		}
	}
	double leftScore=getNodeScore(anode->leftchild);
	double rightScore=getNodeScore(anode->rightchild);
	double score=nodeContrib+leftScore+rightScore;
	return score;
}


int
GammaManager::getAllClusterAssignments(map<int,map<string,int>*>& allClusterAssignments, bool revisit)
{
	//First clean up
	for(map<int,map<string,int>*>::iterator aIter=allClusterAssignments.begin();aIter!=allClusterAssignments.end();aIter++)
	{
		aIter->second->clear();
		delete aIter->second;
	}
	for(map<string,Matrix*>::iterator aIter=clusterFlipCnt.begin();aIter!=clusterFlipCnt.end();aIter++)
	{
		delete aIter->second;
	}
	allClusterAssignments.clear();
	clusterFlipCnt.clear();
	//Init the ML estimates
	spdistMgr->initTransitionMatrix_ML();
	for(map<int,Gamma*>::iterator gIter=gammaSet.begin();gIter!=gammaSet.end();gIter++)
	{
		int ogid=gIter->first;
		if(gIter->first==20057|| gIter->first==25639)
		{
			cout <<"Stop here" << endl;
		}
		Gamma* gamma=gIter->second;
		if(gIter->first==1325)
		{
			gamma->showTree();
		}
		map<string,int>* maximalAssignment=new map<string,int>;	
		int dupID=1;
		getMaximalAssignment(gamma->root,maximalAssignment,dupID);
		map<string,int> parentCAssign;
		map<string,int> extantCAssign;
		for(map<string,int>::iterator aIter=maximalAssignment->begin();aIter!=maximalAssignment->end();aIter++)
		{
			if(strstr(aIter->first.c_str(),"Anc")==NULL)
			{
				continue;
			}
			char ancName[256];
			strcpy(ancName,aIter->first.c_str());
			char* pos=strchr(ancName,':');
			*pos='\0';
			string key(pos+1);
			parentCAssign[key]=aIter->second;
			extantCAssign[key]=aIter->second;
		}
		spdistMgr->getExtantClustering(parentCAssign,extantCAssign);
		updateTransitionMatrix_ML(gamma->root,-1);
		bool dup=hasDuplicate(gamma->root);
		if(!dup && revisit)
		{
			revisitLeafAssignments(gamma->root, extantCAssign,maximalAssignment);
		}
		allClusterAssignments[gIter->first]=maximalAssignment;
	}
	spdistMgr->normalizeTransitionMatrix_ML();
	return 0;
}


int
GammaManager::getAllClusterAssignments_Grouped(map<int,map<string,int>*>& allClusterAssignments)
{
	//First clean up
	for(map<int,map<string,int>*>::iterator aIter=allClusterAssignments.begin();aIter!=allClusterAssignments.end();aIter++)
	{
		aIter->second->clear();
		delete aIter->second;
	}
	for(map<string,Matrix*>::iterator aIter=clusterFlipCnt.begin();aIter!=clusterFlipCnt.end();aIter++)
	{
		delete aIter->second;
	}
	allClusterAssignments.clear();
	for(map<int,Gamma*>::iterator gIter=gammaSet.begin();gIter!=gammaSet.end();gIter++)
	{
		int ogid=gIter->first;
		if(gIter->first==20057|| gIter->first==25639)
		{
			cout <<"Stop here" << endl;
		}
		Gamma* gamma=gIter->second;
		if(gIter->first==1325)
		{
			gamma->showTree();
		}
		map<string,int>* maximalAssignment=new map<string,int>;	
		int dupID=1;
		getMaximalAssignment_Grouped(gamma->root,maximalAssignment,dupID);
		allClusterAssignments[gIter->first]=maximalAssignment;
	}
	return 0;
}

int
GammaManager::getAllClusterAssignments_Conditional(map<int,map<string,int>*>& allClusterAssignments)
{
	for(map<int,map<string,int>*>::iterator aIter=allClusterAssignments.begin();aIter!=allClusterAssignments.end();aIter++)
	{
		aIter->second->clear();
		delete aIter->second;
	}
	allClusterAssignments.clear();
	//Init the ML estimates
	spdistMgr->initTransitionMatrix_ML();
	for(map<int,Gamma*>::iterator gIter=gammaSet.begin();gIter!=gammaSet.end();gIter++)
	{
		int ogid=gIter->first;
		if(gIter->first==20057|| gIter->first==25639)
		{
			cout <<"Stop here" << endl;
		}
		Gamma* gamma=gIter->second;
		if(gIter->first==1325)
		{
			gamma->showTree();
		}
		map<string,int>* maximalAssignment=new map<string,int>;	
		int dupID=1;
		getMaximalConditionalAssignment(gamma->root,maximalAssignment,dupID,-1);
		map<string,int> parentCAssign;
		map<string,int> extantCAssign;
		for(map<string,int>::iterator aIter=maximalAssignment->begin();aIter!=maximalAssignment->end();aIter++)
		{
			if(strstr(aIter->first.c_str(),"Anc")==NULL)
			{
				continue;
			}
			char ancName[256];
			strcpy(ancName,aIter->first.c_str());
			char* pos=strchr(ancName,':');
			*pos='\0';
			string key(pos+1);
			parentCAssign[key]=aIter->second;
			extantCAssign[key]=aIter->second;
		}
		spdistMgr->getExtantClustering(parentCAssign,extantCAssign);
		updateTransitionMatrix_ML(gamma->root,-1);
		bool dup=hasDuplicate(gamma->root);
		if(!dup)
		{
			revisitLeafAssignments(gamma->root, extantCAssign,maximalAssignment);
		}
		allClusterAssignments[gIter->first]=maximalAssignment;
	}
	spdistMgr->normalizeTransitionMatrix_ML();
	return 0;
}

int
GammaManager::getAllClusterGammas(map<int,map<string,map<int,double>*>*>& allClusterGammas)
{
	for(map<int,Gamma*>::iterator gIter=gammaSet.begin();gIter!=gammaSet.end();gIter++)
	{
		Gamma* gamma=gIter->second;
		map<string,map<int,double>*>* maximalAssignment=new map<string,map<int,double>*>;	
		getClusterProb(gamma->root,maximalAssignment);
		allClusterGammas[gIter->first]=maximalAssignment;
	}
	return 0;
}

int
GammaManager::getMaximalAssignment(Gamma::Node* anode,map<string,int>* assignment)
{
	if(anode==NULL)
	{
		return 0;
	}
	if(strstr(anode->name.c_str(),"orf19.5809")!=NULL)
	{
		cout << "Stop here " << endl;
	}
	if(anode->nodeType==1)
	{
		if(anode->normTerm!=NULL)
		{
			double maxScore=0;
			int maxClusterID=-1;
			Matrix* margp=anode->normTerm;
			if(margp->getColCnt()==1)
			{
				margp=anode->gamma;
			}
			for(int c=0;c<margp->getColCnt();c++)
			{
				double val=margp->getValue(0,c);
				if(val>maxScore)
				{
					maxScore=val;
					maxClusterID=c;
				}
			}
			string key;
			if(strstr(anode->name.c_str(),"Anc")!=NULL)
			{
				char temp[256];
				sprintf(temp,"%s_1:%s",anode->species.c_str(),anode->species.c_str());
				key.append(temp);
				if(assignment->find(key)==assignment->end())
				{
					(*assignment)[key]=maxClusterID;
				}
				else
				{
					key.clear();
					sprintf(temp,"%s_2:%s",anode->species.c_str(),anode->species.c_str());
					key.append(temp);
					if(assignment->find(key)==assignment->end())
					{
						(*assignment)[key]=maxClusterID;
					}
					else
					{
						cout <<"Found more instances than expected for "<< anode->name << endl;	
					}
				}	
			}
			else
			{
				key.append(anode->name);
				key.append(":");
				key.append(anode->species);
				(*assignment)[key]=maxClusterID;
			}
		}
	}
	getMaximalAssignment(anode->leftchild,assignment);
	getMaximalAssignment(anode->rightchild,assignment);
	return 0;
}

int
GammaManager::updateTransitionMatrix_ML(Gamma::Node* node, int parentClusterID)
{
	if(node==NULL)
	{
		return 0;	
	}
	if(node->nodeType==1)
	{
		//Now update the ML transition matrix based on the assignments
		Matrix* conditional_ml=spdistMgr->getTransitionMatrix_ML(node->species);
		if(node->normTerm!=NULL)
		{
			double maxScore=0;
			int maxClusterID=-1;
			Matrix* margp=node->normTerm;
			if(margp->getColCnt()==1)
			{
				margp=node->gamma;
			}
			for(int c=0;c<margp->getColCnt();c++)
			{
				double val=margp->getValue(0,c);
				if(val>maxScore)
				{
					maxScore=val;
					maxClusterID=c;
				}
			}
			if(parentClusterID==-1)
			{
				double currval=conditional_ml->getValue(0,maxClusterID);
				conditional_ml->setValue(currval+1,0,maxClusterID);
			}
			else
			{
				double currval=conditional_ml->getValue(parentClusterID,maxClusterID);
				conditional_ml->setValue(currval+1,parentClusterID,maxClusterID);
			}
			updateTransitionMatrix_ML(node->leftchild,maxClusterID);
			updateTransitionMatrix_ML(node->rightchild,maxClusterID);
		}
	}
	else 
	{
		updateTransitionMatrix_ML(node->leftchild,parentClusterID);
		updateTransitionMatrix_ML(node->rightchild,parentClusterID);
	}
	return 0;
}

const char* 
GammaManager::getDupAncestor(int ogid)
{
	const char* dupAnc=NULL;
	if(ogDupAncMap.find(ogid)!=ogDupAncMap.end())
	{
		dupAnc=ogDupAncMap[ogid].c_str();
	}
	return dupAnc;
}

int
GammaManager::getMaximalAssignment(Gamma::Node* anode,map<string,int>* assignment,int dupID)
{
	if(anode==NULL)
	{
		return 0;
	}
	if(strcmp(anode->name.c_str(),"Kwal6.520")==0)
	{	
		cout << "Stop here " << endl;
	}
	int localDupID=dupID;
	if(anode->nodeType==1)
	{
		localDupID=dupID;
		int maxClusterID=-1;
		if(anode->normTerm!=NULL)
		{
			double maxScore=0;
			Matrix* margp=anode->normTerm;
			if(margp->getColCnt()==1)
			{
				margp=anode->gamma;
			}
			for(int c=0;c<margp->getColCnt();c++)
			{
				double val=margp->getValue(0,c);
				if(val>maxScore)
				{
					maxScore=val;
					maxClusterID=c;
				}
			}
		}
		string key;
		if(strstr(anode->name.c_str(),"Anc")!=NULL)
		{
			char temp[256];
			sprintf(temp,"%s_%d:%s",anode->species.c_str(),dupID,anode->species.c_str());
			key.append(temp);
			if(assignment->find(key)==assignment->end())
			{
				(*assignment)[key]=maxClusterID;
			}
			else
			{
				cout <<"Found more instances than expected for "<< anode->name << endl;	
			}	
		}
		else
		{
			key.append(anode->name);
			key.append(":");
			key.append(anode->species);
			(*assignment)[key]=maxClusterID;
		}
		getMaximalAssignment(anode->leftchild,assignment,localDupID);
		getMaximalAssignment(anode->rightchild,assignment,localDupID);
	}
	else 
	{
		getMaximalAssignment(anode->leftchild,assignment,localDupID);
		getMaximalAssignment(anode->rightchild,assignment,localDupID+1);
	}

	return 0;
}


//In addition to the ancestral species, we also add an ID to the extant species and put them all
//in a set. 
int
GammaManager::getMaximalAssignment_Grouped(Gamma::Node* anode,map<string,int>* assignment,int dupID)
{
	if(anode==NULL)
	{
		return 0;
	}
	if(strcmp(anode->name.c_str(),"Kwal6.520")==0)
	{	
		cout << "Stop here " << endl;
	}
	int localDupID=dupID;
	if(anode->nodeType==1)
	{
		localDupID=dupID;
		int maxClusterID=-1;
		if(anode->normTerm!=NULL)
		{
			double maxScore=0;
			Matrix* margp=anode->normTerm;
			if(margp->getColCnt()==1)
			{
				margp=anode->gamma;
			}
			for(int c=0;c<margp->getColCnt();c++)
			{
				double val=margp->getValue(0,c);
				if(val>maxScore)
				{
					maxScore=val;
					maxClusterID=c;
				}
			}
		}
		string key;
		if(strstr(anode->name.c_str(),"Anc")!=NULL)
		{
			char temp[256];
			sprintf(temp,"%d:%s:%s",localDupID,anode->species.c_str(),anode->species.c_str());
			key.append(temp);
			if(assignment->find(key)==assignment->end())
			{
				(*assignment)[key]=maxClusterID;
			}
			else
			{
				cout <<"Found more instances than expected for "<< anode->name << endl;	
			}	
		}
		else
		{
			char temp[256];
			sprintf(temp,"%d:%s:%s",localDupID,anode->species.c_str(),anode->name.c_str());
			key.append(temp);
			(*assignment)[key]=maxClusterID;
		}
		getMaximalAssignment_Grouped(anode->leftchild,assignment,localDupID);
		getMaximalAssignment_Grouped(anode->rightchild,assignment,localDupID);
	}
	else 
	{
		getMaximalAssignment_Grouped(anode->leftchild,assignment,localDupID);
		getMaximalAssignment_Grouped(anode->rightchild,assignment,localDupID+1);
	}
	return 0;
}


int
GammaManager::getMaximalConditionalAssignment(Gamma::Node* anode,map<string,int>* assignment,int dupID, int parentClusterID)
{
	if(anode==NULL)
	{
		return 0;
	}
	int localDupID=dupID;
	if(anode->nodeType==1)
	{
		localDupID=dupID;
		double maxScore=0;
		int maxClusterID=-1;
		if(anode->normTerm!=NULL)
		{
			Matrix* condp=anode->gamma;
			int localParentClusterID=-1;
			if(condp->getRowCnt()==1)
			{
				localParentClusterID=0;
			}
			else
			{
				localParentClusterID=parentClusterID;
			}
			Matrix* conditional_ml=spdistMgr->getTransitionMatrix_ML(anode->species);
			for(int c=0;c<condp->getColCnt();c++)
			{
				double val=condp->getValue(localParentClusterID,c);
				if(val>maxScore)
				{
					maxScore=val;
					maxClusterID=c;
				}
			}
			if(parentClusterID==-1)
			{
				double currval=conditional_ml->getValue(0,maxClusterID);
				conditional_ml->setValue(currval+1,0,maxClusterID);
			}
			else
			{
				double currval=conditional_ml->getValue(parentClusterID,maxClusterID);
				conditional_ml->setValue(currval+1,parentClusterID,maxClusterID);
			}
			string key;
			if(strstr(anode->name.c_str(),"Anc")!=NULL)
			{
				char temp[256];
				sprintf(temp,"%s_%d:%s",anode->species.c_str(),dupID,anode->species.c_str());
				key.append(temp);
				if(assignment->find(key)==assignment->end())
				{
					(*assignment)[key]=maxClusterID;
				}
				else
				{
					cout <<"Found more instances than expected for "<< anode->name << endl;	
				}	
			}
			else
			{
				key.append(anode->name);
				key.append(":");
				key.append(anode->species);
				(*assignment)[key]=maxClusterID;
			}
		}
		getMaximalConditionalAssignment(anode->leftchild,assignment,localDupID,maxClusterID);
		getMaximalConditionalAssignment(anode->rightchild,assignment,localDupID,maxClusterID);
	}
	else 
	{
		getMaximalConditionalAssignment(anode->leftchild,assignment,localDupID,parentClusterID);
		getMaximalConditionalAssignment(anode->rightchild,assignment,localDupID+1,parentClusterID);
	}

	return 0;
}

int
GammaManager::revisitLeafAssignments(Gamma::Node* anode, map<string,int>& clusterAssignment_prior, map<string,int>* clusterAssignment_posterior)
{
	if(strstr(anode->species.c_str(),"Anc")!=NULL)
	{
		if(anode->leftchild!=NULL)
		{
			revisitLeafAssignments(anode->leftchild,clusterAssignment_prior,clusterAssignment_posterior);
		}
		if(anode->rightchild!=NULL)
		{
			revisitLeafAssignments(anode->rightchild,clusterAssignment_prior,clusterAssignment_posterior);
		}
		return 0;
	}
	if(anode->nodeType==2 || anode->normTerm==NULL)
	{
		return 0;
	}
	//Get ancestry
	map<int,int> ancestralClusterCnt;
	Gamma::Node* ancestor=anode->parent;
	while(ancestor!=NULL)
	{
		if(clusterAssignment_prior.find(ancestor->species)==clusterAssignment_prior.end())
		{
			cout <<"No species by name " << ancestor->species << endl;
			continue;
		}
		int ancCID=clusterAssignment_prior[ancestor->species];
		if(ancestralClusterCnt.find(ancCID)==ancestralClusterCnt.end())
		{
			ancestralClusterCnt[ancCID]=1;
		}
		else
		{
			ancestralClusterCnt[ancCID]=ancestralClusterCnt[ancCID]+1;
		}
		ancestor=ancestor->parent;
	}
	//If my cluster assignment does not agree with any of the ancestors, pick something from the prior or something from the ancestor. the majority vote of ancestor
	//should agree with the prior
	string key;
	key.append(anode->name);
	key.append(":");
	key.append(anode->species);
	if(strcmp(anode->name.c_str(),"Calb")==0)
	{
		cout <<"Stop here " << endl;
	}
	if(clusterAssignment_posterior->find(key)==clusterAssignment_posterior->end())
	{
		return 0;
	}
	int c_postClusterAssign=(*clusterAssignment_posterior)[key];
	if(ancestralClusterCnt.find(c_postClusterAssign)==ancestralClusterCnt.end())
	{
		//Changing the cluster assignment of this gene
		int c_MajorityAncestor=-1;
		int maxCnt=0;
		for(map<int,int>::iterator mIter=ancestralClusterCnt.begin();mIter!=ancestralClusterCnt.end();mIter++)
		{
			if(mIter->second>maxCnt)
			{
				maxCnt=mIter->second;
				c_MajorityAncestor=mIter->first;
			}
		}
		int c_Prior=clusterAssignment_prior[anode->species];
		//if(c_Prior==c_MajorityAncestor || ancestralClusterCnt.size()==1)
		if(ancestralClusterCnt.size()==1)
		//if(c_Prior==c_MajorityAncestor) 
		{
			(*clusterAssignment_posterior)[key]=c_MajorityAncestor;
			//If I am changing things, I should also update the gamma of this node
			Matrix* gvals=anode->gamma;
			for(int r=0;r<gvals->getRowCnt();r++)
			{
				double src=gvals->getValue(r,c_postClusterAssign);
				double dest=gvals->getValue(r,c_MajorityAncestor);
				gvals->setValue(src,r,c_MajorityAncestor);
				gvals->setValue(dest,r,c_postClusterAssign);
			}
			//Also keep a tab of which clusters may switch per species
			Matrix* cnt=NULL;
			if(clusterFlipCnt.find(anode->species)==clusterFlipCnt.end())
			{
				cnt=new Matrix(gvals->getRowCnt(),gvals->getColCnt());
				cnt->setAllValues(0);
				clusterFlipCnt[anode->species]=cnt;
			}
			else	
			{
				cnt=clusterFlipCnt[anode->species];
			}
			double acnt=cnt->getValue(c_postClusterAssign,c_MajorityAncestor);
			cnt->setValue(acnt+1,c_postClusterAssign,c_MajorityAncestor);
		}
		/*if(c_Prior==c_MajorityAncestor)
		{
			(*clusterAssignment_posterior)[key]=c_Prior;
		}
		else if(ancestralClusterCnt.size()==1)
		{
			(*clusterAssignment_posterior)[key]=c_MajorityAncestor;
		}*/
	}
	ancestralClusterCnt.clear();
	return 0;
}


int
GammaManager::getClusterProb(Gamma::Node* anode,map<string,map<int,double>*>* prob)
{
	if(anode==NULL)
	{
		return 0;
	}
	if(anode->nodeType==1)
	{
		if(anode->normTerm!=NULL)
		{
			map<int,double>* pmap=new map<int,double>;
			int maxClusterID=-1;
			Matrix* margp=anode->normTerm;
			if(margp->getColCnt()==1)
			{
				margp=anode->gamma;
			}
			for(int c=0;c<margp->getColCnt();c++)
			{
				double val=margp->getValue(0,c);
				(*pmap)[c]=val;
			}
			string key(anode->name);
			key.append(":");
			key.append(anode->species);
			(*prob)[key]=pmap;
		}
	}
	getClusterProb(anode->leftchild,prob);
	getClusterProb(anode->rightchild,prob);
	return 0;
}

int
GammaManager::estimateTransitionProbability()
{
	spdistMgr->resetTransitionProbability();
	for(map<int,Gamma*>::iterator gIter=gammaSet.begin();gIter!=gammaSet.end();gIter++)
	{
		if(gIter->first==13692)
		{
			cout <<"Stop here" << endl;
		}
		//cout <<"Estimating transition for OG : " << gIter->first << endl;
		Gamma* gamma=gIter->second;
		Gamma::Node* root=gamma->root;
		updateTransitionMatrix(root);
	}
	spdistMgr->normalizeTransitionMatrix();
	return 0;
}

/*int
GammaManager::updateTransitionMatrix(Gamma::Node* node)
{
	if(node->normTerm==NULL)
	{
		return 0;
	}
	Matrix* condmatrix=spdistMgr->getConditional(node->species);
	Matrix* gval=node->gamma;
	
	//use the norm term to store the marginals of this node, which will be used in the child node
	Matrix* parentmarg=NULL;
	if(node->parent!=NULL)
	{
		parentmarg=node->parent->normTerm;
		if(parentmarg==NULL)
		{
			cout << "Parent normterm null for " << node->name << "\t" << node->parent->species << endl;
		}
		//cout <<"Estimating normterm for " << node->name << endl;
		if(parentmarg->getColCnt()==1)
		{
			if(node->parent->parent==NULL)
			{
				//cout <<"Check is correct" << endl;
			}
			parentmarg=node->parent->gamma;
		}
		for(int c=0;c<gval->getColCnt();c++)
		{
			double anormterm=0;
			for(int r=0;r<gval->getRowCnt();r++)
			{	
				double g=gval->getValue(r,c);
				double margp=parentmarg->getValue(0,r);
				double jointp=g*margp;
				if(isnan(jointp) || isinf(jointp))
				{
					cout <<"Posterior is nan/inf for  "<< node->name <<" species " << node->species<< " cluster id "  << c
					<< " gamma=" << g<< " marg=" << margp<< endl;
				}
				gval->setValue(jointp,r,c);
				anormterm=anormterm+jointp;
			}
			node->normTerm->setValue(anormterm,0,c);
		}
		//node->normTerm->showMatrix();
		//node->gamma->showMatrix();
	}
	for(int c=0;c<gval->getColCnt();c++)
	{
		for(int r=0;r<gval->getRowCnt();r++)
		{	
			double currval=condmatrix->getValue(r,c);
			double g=gval->getValue(r,c);
			condmatrix->setValue(currval+g,r,c);
			//change the gamma back to the conditional
			if(parentmarg!=NULL)
			{
				double margp=parentmarg->getValue(0,r);
				//double condp=g/margp;
				//gval->setValue(condp,r,c);
			}
		}
	}
	if(node->leftchild!=NULL)
	{
		updateTransitionMatrix(node->leftchild);
	}
	if(node->rightchild!=NULL)
	{
		updateTransitionMatrix(node->rightchild);
	}
	return 0;
}*/

int
GammaManager::updateTransitionMatrix(Gamma::Node* node)
{
	if(node->normTerm==NULL)
	{
		return 0;
	}
	if(node->nodeType==2)
	{
		//Then the marginal is simply the marginal of its parent
		if(node->parent==NULL)
		{
			return 0;
		}
		Matrix* amarg=node->parent->normTerm;
		if(amarg->getColCnt()==1)
		{
			amarg=node->parent->gamma;
		}
		//If I am a node of type 2, the just inherit the marginal from the parent
		for(int r=0;r<amarg->getColCnt();r++)
		{
			node->normTerm->setValue(amarg->getValue(0,r),0,r);
		}
		if(node->leftchild!=NULL)
		{
			updateTransitionMatrix(node->leftchild);
		}
		if(node->rightchild!=NULL)
		{
			updateTransitionMatrix(node->rightchild);
		}
		return 0;
	}
	Matrix* condmatrix=spdistMgr->getConditional(node->species);
	Matrix* gval=node->gamma;
	Matrix* parentmarg=NULL;
	//use the norm term to store the marginals of this node, which will be used in the child node
	if(node->parent!=NULL)
	{
		parentmarg=node->parent->normTerm;
		if(parentmarg->getColCnt()==1)
		{
			parentmarg=node->parent->gamma;
		}
		for(int c=0;c<gval->getColCnt();c++)
		{
			double anormterm=0;
			for(int r=0;r<gval->getRowCnt();r++)
			{	
				double g=gval->getValue(r,c);
				double margp=parentmarg->getValue(0,r);
				double jointp=g*margp;
				if(isnan(jointp) || isinf(jointp))
				{
					cout <<"Posterior is nan/inf for  "<< node->name <<" species " << node->species<< " cluster id "  << c
					<< " gamma=" << g<< " marg=" << margp<< endl;
				}
				gval->setValue(jointp,r,c);
				anormterm=anormterm+jointp;
			}
			node->normTerm->setValue(anormterm,0,c);
		}
		//node->normTerm->showMatrix();
		//node->gamma->showMatrix();
	}
	for(int c=0;c<gval->getColCnt();c++)
	{
		for(int r=0;r<gval->getRowCnt();r++)
		{	
			double currval=condmatrix->getValue(r,c);
			double g=gval->getValue(r,c);
			condmatrix->setValue(currval+g,r,c);
		}
	}
	if(node->leftchild!=NULL)
	{
		updateTransitionMatrix(node->leftchild);
	}
	if(node->rightchild!=NULL)
	{
		updateTransitionMatrix(node->rightchild);
	}
	return 0;
}


int 
GammaManager::reestimateTransitionProbability()
{
	spdistMgr->resetTransitionProbability();
	for(map<int,Gamma*>::iterator gIter=gammaSet.begin();gIter!=gammaSet.end();gIter++)
	{
		if(gIter->first==3494)
		{
			cout <<"Stop here" << endl;
		}
		Gamma* gamma=gIter->second;
		Gamma::Node* root=gamma->root;
		updateTransitionMatrix_GammaFix(root);
	}
	spdistMgr->normalizeTransitionMatrix();
	
	return 0;
}


int
GammaManager::updateTransitionMatrix_GammaFix(Gamma::Node* node)
{
	if(node->normTerm==NULL)
	{
		return 0;
	}
	if(node->nodeType==2)
	{
		//Then the marginal is simply the marginal of its parent
		if(node->parent==NULL)
		{
			return 0;
		}
		Matrix* amarg=node->parent->normTerm;
		if(amarg->getColCnt()==1)
		{
			amarg=node->parent->gamma;
		}
		//If I am a node of type 2, the just inherit the marginal from the parent
		for(int r=0;r<amarg->getColCnt();r++)
		{
			node->normTerm->setValue(amarg->getValue(0,r),0,r);
		}
		if(node->leftchild!=NULL)
		{
			updateTransitionMatrix_GammaFix(node->leftchild);
		}
		if(node->rightchild!=NULL)
		{
			updateTransitionMatrix_GammaFix(node->rightchild);
		}
		return 0;
	}
	Matrix* condmatrix=spdistMgr->getConditional(node->species);
	Matrix* gval=node->gamma;
	
	//use the norm term to store the marginals of this node, which will be used in the child node
	Matrix* parentmarg=NULL;
	if(node->parent!=NULL)
	{
		parentmarg=node->parent->normTerm;
		//cout <<"Estimating normterm for " << node->name << endl;
		if(parentmarg->getColCnt()==1)
		{
			if(node->parent->parent==NULL)
			{
				//cout <<"Check is correct" << endl;
			}
			parentmarg=node->parent->gamma;
		}
		for(int c=0;c<gval->getColCnt();c++)
		{
			double anormterm=0;
			for(int r=0;r<gval->getRowCnt();r++)
			{	
				double jointp=gval->getValue(r,c);
				if(isnan(jointp) || isinf(jointp))
				{
					cout <<"joint is nan/inf for  "<< node->name <<" species " << node->species<< " cluster id "  << c
					<< " gamma=" << jointp<< endl;
				}
				gval->setValue(jointp,r,c);
				anormterm=anormterm+jointp;
			}
			node->normTerm->setValue(anormterm,0,c);
		}
	}
	for(int c=0;c<gval->getColCnt();c++)
	{
		for(int r=0;r<gval->getRowCnt();r++)
		{	
			double currval=condmatrix->getValue(r,c);
			double g=gval->getValue(r,c);
			condmatrix->setValue(currval+g,r,c);
		}
	}
	if(node->leftchild!=NULL)
	{
		updateTransitionMatrix_GammaFix(node->leftchild);
	}
	if(node->rightchild!=NULL)
	{
		updateTransitionMatrix_GammaFix(node->rightchild);
	}
	return 0;
}



int 
GammaManager::showClusterFlipCnts()
{
	if(clusterFlipCnt.size()==0)
	{
		return 0;
	}
	cout << "Cluster Flips per species" << endl;
	for(map<string,Matrix*>::iterator aIter=clusterFlipCnt.begin();aIter!=clusterFlipCnt.end();aIter++)
	{
		cout << aIter->first << endl;
		Matrix* m=aIter->second;
		m->showMatrix();
	}
	return 0;
}



Matrix* 
GammaManager::getGamma(int ogid,string& geneName, string& specName)
{
	Gamma* gamma_og=gammaSet[ogid];
	Matrix* g=gamma_og->getGamma(geneName,specName);
	return g;
}

Matrix* 
GammaManager::getNormTerm(int ogid,string& geneName, string& specName)
{
	Gamma* gamma_og=gammaSet[ogid];
	Matrix* g=gamma_og->getNormTerm(geneName,specName);
	return g;
}



double
GammaManager::getPrior(string& speciesname,int parentClusterID,int nodeClusterID)
{
	Matrix* conditional=spdistMgr->getConditional(speciesname);
	double condprob=conditional->getValue(parentClusterID,nodeClusterID);
	return condprob;
}

bool
GammaManager::hasDuplicate(Gamma::Node* anode)
{
	if(anode->nodeType==2)
	{
		return true;
	}
	bool leftStatus=false;
	if(anode->leftchild!=NULL)
	{
		leftStatus=hasDuplicate(anode->leftchild);
	}
	bool rightStatus=false;
	if(anode->rightchild!=NULL)
	{
		rightStatus=hasDuplicate(anode->rightchild);
	}
	if(leftStatus==true || rightStatus==true)
	{
		return true;
	}
	return false;
}

int
GammaManager::sampleAllClusterAssignments(map<int,map<string,int>*>& allClusterAssignments)
{
	r=gsl_rng_alloc(gsl_rng_default);
	int rseed=getpid();
	gsl_rng_set(r,rseed);
	for(map<int,Gamma*>::iterator gIter=gammaSet.begin();gIter!=gammaSet.end();gIter++)
	{
		int ogid=gIter->first;
		Gamma* gamma=gIter->second;
		map<string,int>* assignment=new map<string,int>;	
		int dupID=1;
		drawAssignment(gamma->root,assignment,dupID,-1);
		allClusterAssignments[gIter->first]=assignment;
	}
	return 0;
}

int
GammaManager::drawAssignment(Gamma::Node* anode, map<string,int>* assignment,int dupID,int parentID)
{
	if(anode==NULL)
	{
		return 0;
	}
	int localDupID=dupID;
	if(anode->nodeType==1)
	{
		localDupID=dupID;
		int clusterID=-1;
		SpeciesDistManager::Species* species=spdistMgr->getSpecies(anode->species);
		if(parentID==-1)
		{
			clusterID=sampleAncestralCluster(r,species);
		}
		else
		{
			clusterID=sampleChildCluster(r,species,parentID);
		}
		string key;
		if(strstr(anode->name.c_str(),"Anc")!=NULL)
		{
			char temp[256];
			sprintf(temp,"%d:%s:%s",localDupID,anode->species.c_str(),anode->species.c_str());
			key.append(temp);
			if(assignment->find(key)==assignment->end())
			{
				(*assignment)[key]=clusterID;
			}
			else
			{
				cout <<"Found more instances than expected for "<< anode->name << endl;	
			}	
		}
		else
		{
			char temp[256];
			sprintf(temp,"%d:%s:%s",localDupID,anode->species.c_str(),anode->name.c_str());
			key.append(temp);
			(*assignment)[key]=clusterID;
		}
		drawAssignment(anode->leftchild,assignment,localDupID,clusterID);
		drawAssignment(anode->rightchild,assignment,localDupID,clusterID);
	}
	else
	{
		drawAssignment(anode->leftchild,assignment,localDupID,parentID);
		drawAssignment(anode->rightchild,assignment,localDupID+1,parentID);
	}
	return 0;
}


int 
GammaManager::sampleAncestralCluster(gsl_rng* r,SpeciesDistManager::Species* root)
{
	int clustId=-1;
	double pval=gsl_ran_flat(r,0,1);
	Matrix* params=root->getParams();
	vector<int>* sortedClustIDs=root->getSortedClusterIDs(0);
	if(sortedClustIDs==NULL)
	{
		sortedClustIDs=new vector<int>;
		for(int i=0;i<params->getColCnt();i++)
		{
			sortedClustIDs->push_back(i);
		}
		for(int i=0;i<params->getColCnt();i++)
		{
			for(int j=i+1;j<params->getColCnt();j++)
			{
				double v1=params->getValue(0,(*sortedClustIDs)[i]);
				double v2=params->getValue(0,(*sortedClustIDs)[j]);
				if(v1<v2)
				{
					int oldval=(*sortedClustIDs)[i];
					(*sortedClustIDs)[i]=(*sortedClustIDs)[j];
					(*sortedClustIDs)[j]=oldval;
				}
			}
		}
		root->setSortedClusterIDs(0,sortedClustIDs);
	}
	double cdf=params->getValue(0,(*sortedClustIDs)[0]);
	clustId=0;
	while(pval>cdf)
	{
		clustId++;
		cdf=cdf+params->getValue(0,(*sortedClustIDs)[clustId]);
	}
	int actualClustID=(*sortedClustIDs)[clustId];
	return actualClustID;
}

int 
GammaManager::sampleChildCluster(gsl_rng* r,SpeciesDistManager::Species* child,int parentID)
{
	vector<int>* sortedClustIDs=child->getSortedClusterIDs(parentID);
	if(sortedClustIDs==NULL)
	{
		sortedClustIDs=new vector<int>;
		sortIndices(child->getParams(),parentID,sortedClustIDs);
		child->setSortedClusterIDs(parentID,sortedClustIDs);
	}
	int childID=sampleChildCluster(r,parentID,child->getParams(),sortedClustIDs);
	if(child->leftchild!=NULL)
	{
		sampleChildCluster(r,child->leftchild,childID);
		sampleChildCluster(r,child->rightchild,childID);
	}
	return childID;
}

int 
GammaManager::sampleChildCluster(gsl_rng* r, int parentClusterId,Matrix* params,vector<int>* sortedClustIDs)
{
	int clustId=-1;
	double pval=gsl_ran_flat(r,0,1);
	double cdf=params->getValue(parentClusterId,(*sortedClustIDs)[0]);
	clustId=0;
	while(pval>cdf)
	{
		clustId++;
		cdf=cdf+params->getValue(parentClusterId,(*sortedClustIDs)[clustId]);
	}
	int actualClustID=(*sortedClustIDs)[clustId];
	return actualClustID;
}


int
GammaManager::sortIndices(Matrix* params,int parentID,vector<int>* sortedClustIDs)
{
	int colCnt=params->getColCnt();
	for(int i=0;i<colCnt;i++)
	{
		sortedClustIDs->push_back(i);
	}
	for(int i=0;i<params->getColCnt();i++)
	{
		for(int j=i+1;j<params->getColCnt();j++)
		{
			double v1=params->getValue(parentID,(*sortedClustIDs)[i]);
			double v2=params->getValue(parentID,(*sortedClustIDs)[j]);
			if(v1<v2)
			{
				int oldval=(*sortedClustIDs)[i];
				(*sortedClustIDs)[i]=(*sortedClustIDs)[j];
				(*sortedClustIDs)[j]=oldval;
			}
		}
	}
	return 0;
}
