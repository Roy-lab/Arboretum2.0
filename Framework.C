#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string.h>
#include "GeneMap.H"
#include "MappedOrthogroup.H"
#include "MappedOrthogroupReader.H"
#include "Matrix.H"
#include "SpeciesDistManager.H"
#include "GeneTree.H"
#include "GeneTreeManager.H"
#include "Expert.H"
#include "Gamma.H"
#include "GammaManager.H"
#include "GeneNameMapper.H"
#include "SpeciesClusterManager.H"
#include "Framework.H"


Framework::Framework()
{
}

Framework::~Framework()
{
}

int 
Framework::readSpeciesData(const char* aFName, const char* rand)
{
	if(strcmp(rand,"none")==0)
	{
		scMgr.setRandom(false);
	}
	else if(strcmp(rand,"yes")==0)
	{
		scMgr.setRandom(true);
	}
	else if(isdigit(rand[0]))
	{
		scMgr.setRandom(true);
		scMgr.setRandSeed(atoi(rand));
	}
	scMgr.setMaxClusterCnt(maxClusterCnt);
	scMgr.readSpeciesData(aFName);
	randnum=gsl_rng_alloc(gsl_rng_default);
	return 0;
}

int 
Framework::readSpeciesTree(int clusterCnt, const char* aFName)
{
	maxClusterCnt=clusterCnt;
	sdMgr.setMaxClusters(clusterCnt);
	sdMgr.readSpeciesTree(aFName);
	sdMgr.assignLevel();
	gammaMgr.setSpeciesDistManager(&sdMgr);
	return 0;
}

int 
Framework::readOrthology(const char* specOrder, const char* orthomapfile)
{
	mor.readSpeciesMapping(specOrder);
	mor.readFile(orthomapfile);
	scMgr.setOrthogroupReader(&mor);
	gammaMgr.setOrthogroupReader(&mor);
	gammaMgr.setMaxClusterCnt(maxClusterCnt);
	scMgr.setGammaManager(&gammaMgr);
	return 0;
}

int
Framework::setSrcSpecies(const char* specName)
{
	strcpy(srcSpecies,specName);
	scMgr.setSrcSpecies(specName);
	return 0;
}

int
Framework::startClustering(const char* aDir)
{
	strcpy(outputDir,aDir);
	scMgr.initExperts();
	cout <<"Total updated parent nodes "<< gammaMgr.getTotalUpdatedParentCnt() << endl; 
	gammaMgr.showTotalUpdatedParents();
	initClusterTransitionProb();
	scMgr.estimateExpertParameters(outputDir);
	double newScore=scMgr.getScore();
	scMgr.dumpAllInferredClusterAssignments(outputDir);
	double newScore_PP=scMgr.getScore();
	cout <<"Score before PP " << newScore << "\t" << " Score after PP " << newScore_PP << endl;
	scMgr.showClusters_Extant(outputDir);
	scMgr.showClusters_Ancestral(outputDir);
	scMgr.showMeans(outputDir);
	//This is only for visualization purposes
	vector<string> speciesList;
	sdMgr.getSpeciesListPrefix(speciesList);
	for(int i=0;i<speciesList.size();i++)
	{	
		cout << speciesList[i] << endl;
	}
	scMgr.dumpAllInferredClusters_ScerwiseGrouped(outputDir,speciesList);
	scMgr.dumpAllInferredClusters_LCA(outputDir,speciesList,sdMgr.getRoot()->name);
	scMgr.dumpAllInferredClusterGammas(outputDir,speciesList);
	sdMgr.showInferredConditionals(outputDir);
	sdMgr.showInferredConditionals_ML(outputDir);
	return 0;
}


int
Framework::generateData(const char* outputDir)
{
	scMgr.initExperts();
	//Lets do one round of learning shall we.
	initClusterTransitionProb();
	scMgr.estimateExpertParameters(outputDir);
	sdMgr.showInferredConditionals(outputDir);
	scMgr.showMeans(outputDir);
	char dirName[1024];
	sprintf(dirName,"mkdir -p %s/samples",outputDir);
	system(dirName);
	sprintf(dirName,"%s/samples",outputDir);
	vector<string> speciesList;
	sdMgr.getSpeciesListPrefix(speciesList);

	scMgr.generateData(dirName,sdMgr.getRoot()->name,speciesList);
	return 0;
	char fName[1024];
	//Assume we want to generate the same number of genes as in the original data
	SpeciesDistManager::Species* root=sdMgr.getRoot();
	map<string,CLUSTERSET*>& extantSpeciesSet=scMgr.getExtantSpeciesClusters();
	map<string,ofstream*> filePtrSet;
	map<string,ofstream*> filePtrClusteredSet;
	map<string,ofstream*> clusterFilePtrSet;
	map<string,map<int,int>*> speciesClusterDist;
	map<string,map<int,map<string,int>*>*> speciesClusterMembers;
	for(map<string,CLUSTERSET*>::iterator sIter=extantSpeciesSet.begin();sIter!=extantSpeciesSet.end();sIter++)
	{
		sprintf(fName,"%s/%s_samples.txt",outputDir,sIter->first.c_str());
		ofstream* oFile=new ofstream(fName);
		filePtrSet[sIter->first]=oFile;
		sprintf(dirName,"mkdir -p %s/%s",outputDir,sIter->first.c_str());
		system(dirName);
		sprintf(fName,"%s/%s/clusterassign.txt",outputDir,sIter->first.c_str());
		ofstream* cFile=new ofstream(fName);
		clusterFilePtrSet[sIter->first]=cFile;
	}
	map<int,MappedOrthogroup*>& ogSet=mor.getMappedOrthogroups();
	vector<double> sampleValues;
	char clusterAssignmentFName[1024];
	sprintf(clusterAssignmentFName,"%s/clusterassign_multspecies.txt",outputDir);
	ofstream caFile(clusterAssignmentFName);
	string scer(srcSpecies);
	map<string,int>* scerGenes=scMgr.getGenesForSpecies(scer);
	int shown=0;
	for(map<string,int>::iterator gIter=scerGenes->begin();gIter!=scerGenes->end();gIter++)
	{
		int cId=sampleAncestralCluster(randnum,root);
		map<string,int>* clusterAssign=new map<string,int>;
		clusterAssignments[gIter->first]=clusterAssign;
		(*clusterAssign)[root->name]=cId;
		sampleChildCluster(randnum,root->leftchild,cId,*clusterAssign);
		sampleChildCluster(randnum,root->rightchild,cId,*clusterAssign);
		caFile << gIter->first;
		for(map<string,int>::iterator cIter=clusterAssign->begin();cIter!=clusterAssign->end();cIter++)
		{
			cout<<" " <<cIter->first <<"=" << cIter->second;
			caFile<<"\t" <<cIter->first <<"=" << cIter->second;
			map<int,int>* clusterCnt=NULL;
			if(speciesClusterDist.find(cIter->first)==speciesClusterDist.end())
			{
				clusterCnt=new map<int,int>;
				speciesClusterDist[cIter->first]=clusterCnt;
			}
			else
			{
				clusterCnt=speciesClusterDist[cIter->first];
			}
			if(clusterCnt->find(cIter->second)==clusterCnt->end())
			{
				(*clusterCnt)[cIter->second]=1;
			}
			else
			{
				(*clusterCnt)[cIter->second]=(*clusterCnt)[cIter->second]+1;
			}
		}
		caFile<< endl;
		cout << endl;
		int ogid=mor.getMappedOrthogroupID(gIter->first.c_str(),srcSpecies);
		MappedOrthogroup* mgrp=ogSet[ogid];
		//Then for all extant species draw the expression vector
		for(map<string,CLUSTERSET*>::iterator sIter=extantSpeciesSet.begin();sIter!=extantSpeciesSet.end();sIter++)
		{
			CLUSTERSET* clusterSet=sIter->second;
			int specClustId=(*clusterAssign)[sIter->first];
			Expert* e=(*clusterSet)[specClustId];
			e->generateSample(randnum,sampleValues);
			GeneMap* gMap=mgrp->getSpeciesHits(sIter->first.c_str());
			ofstream* oFile=filePtrSet[sIter->first];
			ofstream* cFile=clusterFilePtrSet[sIter->first];
			const string& geneName=gMap->getGeneSet().begin()->first;
			(*cFile) <<geneName <<"\t" << specClustId <<endl;
			(*oFile) << geneName;
			for(int j=0;j<sampleValues.size();j++)
			{
				(*oFile) <<"\t" << sampleValues[j];
			}
			(*oFile) << endl;
			sampleValues.clear();
		}
		shown++;
		//clusterAssign.clear();
	}
	for(map<string,ofstream*>::iterator fIter=filePtrSet.begin();fIter!=filePtrSet.end();fIter++)
	{
		fIter->second->close();
		ofstream* cFile=clusterFilePtrSet[fIter->first];
		cFile->close();
	}
	caFile.close();
	scMgr.showClusters_Extant(outputDir);
	
	//Cluster size dist
	for(map<string,map<int,int>*>::iterator sIter=speciesClusterDist.begin();sIter!=speciesClusterDist.end();sIter++)
	{
		cout <<sIter->first;
		map<int,int>* sizeDist=sIter->second;
		for(map<int,int>::iterator aIter=sizeDist->begin();aIter!=sizeDist->end();aIter++)
		{
			cout <<" " << aIter->first<<":"<< aIter->second;
		}
		cout << endl;
	}

	return 0;
}


int
Framework::redisplay(const char* outputDir)
{
	scMgr.initExperts();
	scMgr.showClusters(outputDir);
	//Only redisplay the data
	return 0;
}


int 
Framework::setPdiagonalLeaf(double aval)
{
	p_diagonal_leaf=aval;
	return 0;
}

int 
Framework::setPdiagonalNonLeaf(double aval)
{
	p_diagonal_nonleaf=aval;
	return 0;
}

//Need to initialize a conditional distribution for every branch which species
//the probability of transitioning from cluster k to cluster j
//TODO: initialize using the ribosomal clusters or some known cluster membership across species
//TODO: Use the rate and the branch length information
int 
Framework::initClusterTransitionProb()
{
	SpeciesDistManager::Species* root=sdMgr.getRoot();
	Matrix* conditional=root->getParams();
	int colcnt=conditional->getColCnt();
	for(int i=0;i<colcnt;i++)
	{
		double aval=1/((double)colcnt);
		conditional->setValue(aval,0,i);
	}
	initClusterTransitionProb(root->leftchild);
	initClusterTransitionProb(root->rightchild);	
	return 0;
}

int
Framework::initClusterTransitionProb(SpeciesDistManager::Species* anode)
{
	cout <<"Transitions for " << anode->name << endl;
	if(anode->leftchild==NULL)
	{
		if(clusterTransitionProb.find(anode->name)==clusterTransitionProb.end())
		{
			cout <<"No cluster transition prob for  " << anode->name << endl;
			exit(0);
		}
		double pval=clusterTransitionProb[anode->name];
		initTransitionProb(anode->getParams(),pval);
	}
	else
	{
		if(clusterTransitionProb.find(anode->name)==clusterTransitionProb.end())
		{
			cout <<"No cluster transition prob for  " << anode->name << endl;
			exit(0);
		}
		double pval=clusterTransitionProb[anode->name];
		initTransitionProb(anode->getParams(),pval);
	}
	if(anode->leftchild!=NULL)
	{
		cout <<"Transitions for " << anode->leftchild->name << endl;
		initClusterTransitionProb(anode->leftchild);
		cout <<"Transitions for " << anode->rightchild->name << endl;
		initClusterTransitionProb(anode->rightchild);
	}
	return 0;
}

//The matrix is supposed to be a transition matrix of going from cluster k to cluster l
int
Framework::initTransitionProb(Matrix* m,double initval)
{
	int rowcnt=m->getRowCnt();
	int colcnt=m->getColCnt();
	for (int i=0;i<rowcnt;i++)
	{
		double s=0;
		for(int j=0;j<colcnt;j++)
		{
			double aval=0;
			if(i==j)
			{
				aval=initval;
			}
			else
			{
				aval=(1-initval)/((double) (rowcnt-1));
			}
			double err=gsl_ran_flat(randnum,0,0.01);
			m->setValue(aval+err,i,j);
			s=s+err+aval;
		}
		for(int j=0;j<colcnt;j++)
		{
			double aval=m->getValue(i,j);	
			aval=aval/s;
			m->setValue(aval,i,j);
		}
	}
	m->showMatrix(1e-5);
	return 0;
}


int 
Framework::inferAncestralClusters(map<int,map<string,int>*>& clusterAssignments)
{
	cout <<"Inferring ancestral clusters" << endl;
	map<string,map<int,int>*> speciesClusterSizeDist;
	int disp=0;
	for(map<int,map<string,int>*>::iterator cIter=clusterAssignments.begin();cIter!=clusterAssignments.end();cIter++)
	{
		map<string,int>* extantClustering=cIter->second;
		map<string,int> ancestralClustering;
		sdMgr.getAncestralClustering(*extantClustering,ancestralClustering);
		for(map<string,int>::iterator aIter=ancestralClustering.begin();aIter!=ancestralClustering.end();aIter++)
		{
			(*extantClustering)[aIter->first]=aIter->second;
		}
		for(map<string,int>::iterator aIter=extantClustering->begin();aIter!=extantClustering->end();aIter++)
		{
			map<int,int>* clusterCnts=NULL;
			if(speciesClusterSizeDist.find(aIter->first)==speciesClusterSizeDist.end())
			{
				clusterCnts=new map<int,int>;
				speciesClusterSizeDist[aIter->first]=clusterCnts;
			}
			else
			{
				clusterCnts=speciesClusterSizeDist[aIter->first];
			}
			if(clusterCnts->find(aIter->second)==clusterCnts->end())
			{
				(*clusterCnts)[aIter->second]=1;
			}
			else
			{
				(*clusterCnts)[aIter->second]=(*clusterCnts)[aIter->second]+1;
			}
		}
		
		/*if(cIter==clusterAssignments.begin())
		{
			cout<<"Gene";
			for(map<string,int>::iterator aIter=extantClustering->begin();aIter!=extantClustering->end();aIter++)
			{
				cout <<" " << aIter->first;
			}
			cout << endl;
		}*/
		if(disp<10)
		{
			cout<<cIter->first;
			for(map<string,int>::iterator aIter=extantClustering->begin();aIter!=extantClustering->end();aIter++)
			{	
				cout <<" " << aIter->first<<"="<< aIter->second;
			} 
			cout << endl;
		}
		disp++;
	
	}
	for(map<string,map<int,int>*>::iterator sIter=speciesClusterSizeDist.begin();sIter!=speciesClusterSizeDist.end();sIter++)
	{
		cout <<sIter->first;
		map<int,int>* csize=sIter->second;
		for(map<int,int>::iterator cIter=csize->begin();cIter!=csize->end();cIter++)
		{
			cout <<" " << cIter->first <<":"<< cIter->second;
		}
		cout << endl;
	}
	return 0;
}

int
Framework::inferExtantClusters(map<int,map<string,int>*>& clusterAssignments)
{
	cout <<"Inferring extant clusters" << endl;
	map<string,map<int,int>*> speciesClusterSizeDist;
	int disp=0;
	for(map<int,map<string,int>*>::iterator cIter=clusterAssignments.begin();cIter!=clusterAssignments.end();cIter++)
	{
		map<string,int>* ancestralClustering=cIter->second;
		map<string,int> extantClustering;
		sdMgr.getExtantClustering(*ancestralClustering,extantClustering);
		for(map<string,int>::iterator aIter=extantClustering.begin();aIter!=extantClustering.end();aIter++)
		{
			(*ancestralClustering)[aIter->first]=aIter->second;
		}
		for(map<string,int>::iterator aIter=ancestralClustering->begin();aIter!=ancestralClustering->end();aIter++)
		{
			map<int,int>* clusterCnts=NULL;
			if(speciesClusterSizeDist.find(aIter->first)==speciesClusterSizeDist.end())
			{
				clusterCnts=new map<int,int>;
				speciesClusterSizeDist[aIter->first]=clusterCnts;
			}
			else
			{
				clusterCnts=speciesClusterSizeDist[aIter->first];
			}
			if(clusterCnts->find(aIter->second)==clusterCnts->end())
			{
				(*clusterCnts)[aIter->second]=1;
			}
			else
			{
				(*clusterCnts)[aIter->second]=(*clusterCnts)[aIter->second]+1;
			}
		}
		/*if(cIter==clusterAssignments.begin())
		{
			cout<<"Gene";
			for(map<string,int>::iterator aIter=ancestralClustering->begin();aIter!=ancestralClustering->end();aIter++)
			{
				cout <<" " << aIter->first;
			}
			cout << endl;
		}*/
		if(disp<10)
		{
			cout<<cIter->first;
			for(map<string,int>::iterator aIter=ancestralClustering->begin();aIter!=ancestralClustering->end();aIter++)
			{	
				cout <<" " << aIter->first<<"="<< aIter->second;
			} 
			cout << endl;
		}
		disp++;
	}
	for(map<string,map<int,int>*>::iterator sIter=speciesClusterSizeDist.begin();sIter!=speciesClusterSizeDist.end();sIter++)
	{
		cout <<sIter->first;
		map<int,int>* csize=sIter->second;
		for(map<int,int>::iterator cIter=csize->begin();cIter!=csize->end();cIter++)
		{
			cout <<" " << cIter->first <<":"<< cIter->second;
		}
		cout << endl;
	}
	return 0;
}

int 
Framework::estimateClusterTransProb(map<int,map<string,int>*>& clusterAssignments)
{
	SpeciesDistManager::Species* root=sdMgr.getRoot();
	//Estimate prior of the root
	Matrix* rootparam=root->conditional;
	rootparam->setAllValues(1e-5);
	for(map<int,map<string,int>*>::iterator aIter=clusterAssignments.begin();aIter!=clusterAssignments.end();aIter++)
	{
		int cid=(*aIter->second)[root->name];
		double aval=rootparam->getValue(cid,0);
		aval=aval+1;
		rootparam->setValue(aval,cid,0);
	}
	for(int i=0;i<rootparam->getRowCnt();i++)
	{
		double aval=rootparam->getValue(i,0)/clusterAssignments.size();
		rootparam->setValue(aval,i,0);
	}
	rootparam->showMatrix(1e-5);
	estimateClusterTransProb(root,root->leftchild,clusterAssignments);
	estimateClusterTransProb(root,root->rightchild,clusterAssignments);
	return 0;
}

int
Framework::estimateClusterTransProb(SpeciesDistManager::Species* parent, SpeciesDistManager::Species* child, map<int,map<string,int>*>& clusterAssignments)
{
	estimateTransitionMatrix(parent->name,child->name,child,clusterAssignments);
	if(child->leftchild!=NULL)
	{
		estimateClusterTransProb(child,child->leftchild,clusterAssignments);
		estimateClusterTransProb(child,child->rightchild,clusterAssignments);
	}
	return 0;
}

int
Framework::estimateTransitionMatrix(string& parentname,string& childname, SpeciesDistManager::Species* child, map<int,map<string,int>*>& clusterAssignments)
{
	Matrix* param=child->getParams();
	param->setAllValues(1);
	for(map<int,map<string,int>*>::iterator cIter=clusterAssignments.begin();cIter!=clusterAssignments.end();cIter++)
	{
		map<string,int>* assignments=cIter->second;
		int ancid=(*assignments)[parentname];
		int childid=(*assignments)[childname];
		double currval=param->getValue(ancid,childid);
		currval=currval+1;
		param->setValue(currval,ancid,childid);
	}
	cout <<"New Params for " << childname << " before norm " << endl;
	param->showMatrix(1e-5);
	//Row is for ancestral cluster id. Cols in a row must add to 1
	for(int i=0;i<param->getRowCnt();i++)
	{
		double sum=0;
		for(int j=0;j<param->getColCnt();j++)
		{
			sum=sum+param->getValue(i,j);
		}
		for(int j=0;j<param->getColCnt();j++)
		{
			double prob=param->getValue(i,j);
			prob=prob/sum;
			param->setValue(prob,i,j);
		}
	}
	cout <<"After normalzation" << childname << endl;
	param->showMatrix(1e-5);
	return 0;
}


bool
Framework::checkConvergence(map<int,map<string,int>*>& clusterAssignments,map<int,map<string,int>*>& oldclusterAssignments)
{
	bool convergence=false;
	int changes=0;
	for(map<int,map<string,int>*>::iterator cIter=clusterAssignments.begin();cIter!=clusterAssignments.end();cIter++)
	{
		map<string,int>* newassign=cIter->second;
		map<string,int>* oldassign=oldclusterAssignments[cIter->first];
		
		for(map<string,int>::iterator nIter=newassign->begin();nIter!=newassign->end();nIter++)
		{
			if(oldassign->find(nIter->first)==oldassign->end())
			{
				cout <<"Fatal error " << endl;
				exit(0);
			}
			int currassignval=(*oldassign)[nIter->first];
			if(currassignval!=nIter->second)
			{
				changes++;
			}
		}
	}
	if(changes==0)
	{
		convergence=true;
	}
	return convergence;

}

int
Framework::saveClusterAssignments(map<int,map<string,int>*>& clusterAssignments,map<int,map<string,int>*> &oldclusterAssignments)
{
	for(map<int,map<string,int>*>::iterator aIter=oldclusterAssignments.begin();aIter!=oldclusterAssignments.end();aIter++)
	{
		//delete aIter->second;
		//aIter->second->clear();
	}
	//oldclusterAssignments.clear();
	for(map<int,map<string,int>*>::iterator bIter=clusterAssignments.begin();bIter!=clusterAssignments.end();bIter++)
	{
		map<string,int>* newassign=NULL;
		if(oldclusterAssignments.find(bIter->first)==oldclusterAssignments.end())
		{
			newassign=new map<string,int>;
			oldclusterAssignments[bIter->first]=bIter->second;
		}
		else
		{
			newassign=oldclusterAssignments[bIter->first];
		}
		map<string,int>* oldassign=bIter->second;
		for(map<string,int>::iterator dIter=oldassign->begin();dIter!=oldassign->end();dIter++)
		{
			//cout <<"updating " << bIter->first <<" "<< dIter->first << " " << dIter->second << endl;
			(*newassign)[dIter->first]=dIter->second;
		}
		oldclusterAssignments[bIter->first]=newassign;
	}
	//clusterAssignments.clear();
	return 0;

}


int 
Framework::sampleAncestralCluster(gsl_rng* r,SpeciesDistManager::Species* root)
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
Framework::sampleChildCluster(gsl_rng* r,SpeciesDistManager::Species* child,int,map<string,int>& clusterAssign)
{
	int parentID=-1;
	if(clusterAssign.find(child->parent->name)==clusterAssign.end())
	{
		cout <<"No parent cluster id for " << child->parent->name << endl;
		exit(0);
	}
	if((strcmp(child->name.c_str(),"Scer")==0) || (strcmp(child->name.c_str(),"Sbay")==0))
	{
	//	cout <<"At child " << child->name << endl;
	}
	parentID=clusterAssign[child->parent->name];
	vector<int>* sortedClustIDs=child->getSortedClusterIDs(parentID);
	if(sortedClustIDs==NULL)
	{
		sortedClustIDs=new vector<int>;
		sortIndices(child->getParams(),parentID,sortedClustIDs);
		child->setSortedClusterIDs(parentID,sortedClustIDs);
	}
	int childID=sampleChildCluster(r,parentID,child->getParams(),sortedClustIDs);
	clusterAssign[child->name]=childID;
	if(child->leftchild!=NULL)
	{
		sampleChildCluster(r,child->leftchild,childID,clusterAssign);
		sampleChildCluster(r,child->rightchild,childID,clusterAssign);
	}
	return 0;
}

int 
Framework::sampleChildCluster(gsl_rng* r, int parentClusterId,Matrix* params,vector<int>* sortedClustIDs)
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
Framework::sortIndices(Matrix* params,int parentID,vector<int>* sortedClustIDs)
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


int 
Framework::dumpInferredAssignments(const char* outputDir,const char* suffix )
{
	char aFName[1024];
	sprintf(aFName,"%s/clusterassign_multspecies_%s_%d.txt",outputDir,suffix,learnIter);
	ofstream oFile(aFName);
	for(map<string,map<string,int>*>::iterator oIter=clusterAssignments.begin();oIter!=clusterAssignments.end();oIter++)
	{
		oFile <<oIter->first;
		map<string,int>* assignments=oIter->second;
		for(map<string,int>::iterator aIter=assignments->begin();aIter!=assignments->end();aIter++)
		{
			oFile <<"\t" << aIter->first<<"=" << aIter->second;
		}
		oFile << endl;
	}
	oFile.close();
	return 0;
}


int 
Framework::setClusterTransProb(double ctransProb)
{
	vector<string> speciesList;
	sdMgr.getSpeciesListPrefix(speciesList);
	for(int s=0;s<speciesList.size();s++)
	{
		clusterTransitionProb[speciesList[s]]=ctransProb;	
	}
	return 0;
}

int 
Framework::setClusterTransProb(const char* aFName)
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
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string speciesName;
		double pval;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				speciesName.append(tok);
			}
			else if(tokCnt==1)
			{
				pval=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		clusterTransitionProb[speciesName]=pval;
	}
	inFile.close();
	return 0;
}

int
main(int argc, const char** argv)
{
	if(argc!=12)
	{
		cout <<"Usage: infAncCluster specorder orthogroup maxk speciestree clusterassignments rand[rseed|none] outputDir mode[learn|generate|visualize] srcSpecies inittype[uniform|branchlength] p_diagonal_nonleaf" << endl;
		cout <<"Inittype is for specifying how the cluster transition probabilities will be initialized. If inittype is uniform then set to everything to the same and if branchlength then set everything from the file" << endl;
		return 0;
	}
	Framework fw;
	fw.readSpeciesTree(atoi(argv[3]),argv[4]);
	fw.readOrthology(argv[1],argv[2]);
	fw.readSpeciesData(argv[5],argv[6]);
	fw.setSrcSpecies(argv[9]);
	if(strcmp(argv[10],"uniform")==0)
	{
		fw.setClusterTransProb(atof(argv[11]));
	}
	else 
	{
		fw.setClusterTransProb(argv[11]);
	}
	if(strcmp(argv[8],"learn")==0)
	{
		fw.startClustering(argv[7]);
	}
	else if(strcmp(argv[8],"generate")==0)
	{	
		fw.generateData(argv[7]);
	}
	else 
	{
		fw.redisplay(argv[7]);
	}
	return 0;
}
