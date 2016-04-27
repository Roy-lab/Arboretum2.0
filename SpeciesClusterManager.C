
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
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include "GeneMap.H"
#include "MappedOrthogroup.H"
#include "MappedOrthogroupReader.H"
#include "Matrix.H"
#include "Expert.H"
#include "GeneExpManager.H"
#include "SpeciesDistManager.H"
#include "GeneTree.H"
#include "GeneTreeManager.H"
#include "Gamma.H"
#include "GammaManager.H"
#include "GeneNameMapper.H"
#include "SpeciesClusterManager.H"

SpeciesClusterManager::SpeciesClusterManager()
{
	randFlag=false;
	CVMode=false;
	secondStage=true;
	r=NULL;
	gnm.readGeneNames();
	rseed=-1;
}

SpeciesClusterManager::~SpeciesClusterManager()
{
}

int
SpeciesClusterManager::setOrthogroupReader(MappedOrthogroupReader* aPtr)
{
	mor=aPtr;
	return 0;
}

int 
SpeciesClusterManager::setGammaManager(GammaManager* aPtr)
{
	gammaMgr=aPtr;
	return 0;
}

//SK: new function to set the GammaManagerObject for test data when working in cross validation mode
int 
SpeciesClusterManager::setGammaManager_Test(GammaManager* aPtr)
{
        gammaMgrTst=aPtr;
        return 0;
}

int
SpeciesClusterManager::setSrcSpecies(const char* aName)
{
	strcpy(srcSpecies,aName);
	return 0;
}

int
SpeciesClusterManager::setMaxClusterCnt(int k)
{
	maxClusterCnt=k;
	return 0;
}

int
SpeciesClusterManager::setRandom(bool flag)
{
	randFlag=flag;
	return 0;
}

int
SpeciesClusterManager::setRandSeed(int aseed)
{
	rseed=aseed;
	return 0;
}

int 
SpeciesClusterManager::readSpeciesData(const char* clusterFName)
{
	ifstream inFile(clusterFName);
	char buffer[1024];
	r=gsl_rng_alloc(gsl_rng_default);
	if(rseed==-1)
	{
		rseed=getpid();
	}
	gsl_rng_set(r,rseed);
	cout << rseed << endl;
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
		string clusterFName;
		string expressionFName;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				speciesName.append(tok);
			}
			else if(tokCnt==1)
			{
				clusterFName.append(tok);
			}
			else if(tokCnt==2)
			{
				expressionFName.append(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		GeneExpManager* exprManager=new GeneExpManager;
		exprManager->readExpression(expressionFName.c_str());
		speciesExprSet[speciesName]=exprManager;
		readClusters(speciesName,clusterFName.c_str());
	}
	inFile.close();
	cout <<"Read clusterings for " << speciesClusterSet_Genewise.size() << " species" << endl;
	return 0;
}

int
SpeciesClusterManager::readSpeciesData(const char* clusterFName,const char* Dir)
{
        ifstream inFile(clusterFName);
        char buffer[1024];
        r=gsl_rng_alloc(gsl_rng_default);
        if(rseed==-1)
        {
                rseed=getpid();
        }
        gsl_rng_set(r,rseed);
        cout << rseed << endl;
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
                string clusterFName;
                string expressionFName;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				speciesName.append(tok);
			}
			else if(tokCnt==1)
			{
				expressionFName.append(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		clusterFName.append(Dir);
		clusterFName.append("/");
		clusterFName.append(speciesName.c_str());
		clusterFName.append("_initial_clusterassign.txt");
		cout << clusterFName << endl;
                GeneExpManager* exprManager=new GeneExpManager;
                exprManager->readExpression(expressionFName.c_str());
                speciesExprSet[speciesName]=exprManager;
                readClusters(speciesName,clusterFName.c_str());
        }
        inFile.close();
        cout <<"Read clusterings for " << speciesClusterSet_Genewise.size() << " species" << endl;
        return 0;
}

int
SpeciesClusterManager::initExperts()
{
	for(map<string,CLUSTERSET*>::iterator aIter=speciesExpertSet.begin();aIter!=speciesExpertSet.end();aIter++)
        {
                if(strcmp(aIter->first.c_str(),srcSpecies)!=0)
                {
                        //continue;
                }
                CLUSTERSET* speciesExpert=aIter->second;
                map<int,int> deleteme;
                for(CLUSTERSET_ITER cIter=speciesExpert->begin();cIter!=speciesExpert->end();cIter++)
                {
                        Expert* e=cIter->second;
			int sample=estimateMeanCov(e,(string&)aIter->first,cIter->first);
                        if(sample<10)
                        {
				if(strcmp(aIter->first.c_str(),srcSpecies)==0)
                		{
					cout << srcSpecies << " cluster " << cIter->first << " has " << sample << " genes." << endl;
					cout << "Exiting program in SpeciesClusterManager::initExperts()" << endl;
                			exit(EXIT_FAILURE);
				}
                                deleteme[cIter->first]=0;
                        }
                }
                for(map<int,int>::iterator cIter=deleteme.begin();cIter!=deleteme.end();cIter++)
                {
                        cout << "Deleting " << cIter->first << endl;
                        CLUSTERSET_ITER dIter=speciesExpert->find(cIter->first);
                        delete dIter->second;
                        speciesExpert->erase(dIter);
                }
        }
        string baseSpecies;
        //SK: set map of Exper objects for base species
        baseSpecies.append(srcSpecies);
        CLUSTERSET* baseSet=speciesExpertSet.find(baseSpecies)->second;
	//SK: check for missing clusters (meaning an Epert object that isn't initialized but should be) in all species other than the base species. 
        for(map<string,CLUSTERSET*>::iterator aIter=speciesExpertSet.begin();aIter!=speciesExpertSet.end();aIter++)
        {
                //SK: ignore if looking at the base species and continue
                if(aIter->first==baseSpecies)
                {
                        cout << "Continuing: source species" << endl;
                        continue;
                }
                //SK: if not a base species check if the number of Expert object is not the number of clusters, meaning one cluster id is missing from the input clusterassignments for that species
                CLUSTERSET* speciesExpert=aIter->second;
                if(speciesExpert->size()<maxClusterCnt)
                {
                        for(int c=0;c<maxClusterCnt;c++)
                        {
                                //SK: check for each cluster id that may be missing an initilization
                                cout << "Cluster " << c << " for " << aIter->first << endl;
                                if(speciesExpert->find(c)==speciesExpert->end())
                                {
                                        //SK:get expert from baske species and copy to a new object for the current species
                                        Expert* i = baseSet->find(c)->second;
                                        Expert* j = new Expert();
                                        //SK: initialize means for the new Expert object from the genes of the base species for this cluster.
                                        estimateMeanCov(j,(string&)baseSpecies,c);
                                        speciesExpert->insert(std::pair<int,Expert*>(c,j));

                                }
                        }
                }
        }
	return 0;
}

/*
int 
SpeciesClusterManager::getExtantClusterAssignment(map<int,map<string,int>*>& clusterAssignments)
{
	for(map<int,map<string,INTDBLMAP*>*>::iterator gIter=gammas.begin();gIter!=gammas.end();gIter++)
	{
		map<string,INTDBLMAP*>* gValSet=gIter->second;
		map<string,int>* clustermembers=NULL;
		if(clusterAssignments.find(gIter->first)==clusterAssignments.end())
		{
			clustermembers=new map<string,int>;
			clusterAssignments[gIter->first]=clustermembers;
		}
		else
		{
			clustermembers=clusterAssignments[gIter->first];
		}
		for(map<string,INTDBLMAP*>::iterator sIter=gValSet->begin();sIter!=gValSet->end();sIter++)
		{
			//For each species get the best cluster id for this og
			INTDBLMAP* gamma_i=sIter->second;
			double maxPval=-1;
			int maxexpertID=-1;
			for(INTDBLMAP_ITER eIter=gamma_i->begin();eIter!=gamma_i->end();eIter++)
			{
				if(eIter->second>maxPval)
				{
					maxPval=eIter->second;
					maxexpertID=eIter->first;
				}
			}
			(*clustermembers)[sIter->first]=maxexpertID;
		}
	}
	return 0;
}*/


map<string,CLUSTERSET*>& 
SpeciesClusterManager::getExtantSpeciesClusters()
{
	return speciesExpertSet;
}

int
SpeciesClusterManager::estimateExpertParameters(const char* outputDir)
{
	double currScore=0;
	bool convergence=false;
	int iter=0;
	//SK: use nMaxIterations to control this while loop
	//SK: this is the maximum number of iterations to run without convergence
	while((!convergence)&&(iter<nMaxIterations))
	//while((!convergence)&&(iter<2))
	{
		expectationStep();
		maximizationStep();
		//dumpAllInferredClusterAssignments(outputDir,iter);
		double newScore=0;
		newScore=getScore();
		double diff=fabs(newScore-currScore);
		//SK: check if the change in the likelihood score is less than the threshold value
		if((iter>0) && (diff<convThresh))
		{
			//SK: if it is, note that the algorith has converged
			convergence=true;
		}
		cout <<"Iter: " << iter << " score: " << newScore << " diff " << diff << endl;
		currScore=newScore;
		iter++;
	}
	if(convergence)
	{
		cout <<"Convergence at iteration:" << iter << endl;
	}
	return 0;
}


int
SpeciesClusterManager::readClusters(string& specName, const char* fName)
{
	ifstream inFile(fName);
	CLUSTERSET* cset=new CLUSTERSET;
	map<string,int>* geneset=new map<string,int>;
	speciesExpertSet[specName]=cset;
	speciesClusterSet_Genewise[specName]=geneset;
	char buffer[1024];
	GeneExpManager* geMgr=speciesExprSet[specName];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		int clustid=0;
		string genename;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				genename.append(tok);
			}	
			else if(tokCnt==1)
			{
				clustid=atoi(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		int ogid=mor->getMappedOrthogroupID(genename.c_str(),specName.c_str());
		if(ogid==-1)
		{
			cout <<"OGID -1 for " << genename << ":" << specName << endl;
			continue;
		}
		if(geMgr->getExp(genename)==NULL)
		{
			continue;
		}
		Expert* expert=NULL;
		if(randFlag)
		{
			if(initRandAssign.find(ogid)==initRandAssign.end())
			{
				//don't use this clustid but some other one
				double pval=gsl_ran_flat(r,0,1);
				double step=1.0/(double)maxClusterCnt;
				int newclustid=(int)(floor(pval/step));
				if(newclustid>=maxClusterCnt)
				{
					newclustid=maxClusterCnt-1;
				}
				clustid=newclustid;
				//initRandAssign[ogid]=clustid;
			}
			else{
				clustid=initRandAssign[ogid];
			}
		}
		//SK: enter this ogid into the list of working orthogroups
		workingOrthoGroups[ogid]=0;
		//SK: if CVMode is turned on and the orthogroup is in the test set...
		if(CVMode==true && testSet.find(ogid)!=testSet.end())
		{
			//SK: initialize this OGID group in the test GammaManager object
                        gammaMgrTst->initGamma(ogid,genename,specName,clustid);
			//SK: continue and apply remaining initialization only to the training set OGIDS
                       	continue;
                }
		if(cset->find(clustid)==cset->end())
		{
			expert=new Expert;
			(*cset)[clustid]=expert;
		}
		else{
			expert=(*cset)[clustid];
		}
		(*geneset)[genename]=clustid;
		expert->assignGeneToExpert(genename.c_str());
		gammaMgr->initGamma(ogid,genename,specName,clustid);
	}
	cout <<"Read " << geneset->size() << " genes in " << specName << endl;
	inFile.close();
	return 0;
}

int 
SpeciesClusterManager::maximizationStep()
{
	for(map<string,CLUSTERSET*>::iterator aIter=speciesExpertSet.begin();aIter!=speciesExpertSet.end();aIter++)
	{
		CLUSTERSET* speciesExpert=aIter->second;
		for(CLUSTERSET_ITER cIter=speciesExpert->begin();cIter!=speciesExpert->end();cIter++)
		{
			Expert* e=cIter->second;
			estimateMeanCov(e,(string&)aIter->first,cIter->first);
		}
	}
	return 0;
}


int
SpeciesClusterManager::expectationStep()
{
	//First get the gammas for each leaf node
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		CLUSTERSET* speciesClusters=sIter->second;
		expectationStep_Species((string&)sIter->first,speciesClusters);
	}
	//The use the speciesdist manager's conditional distributions to infer the rest of the gammas.
	gammaMgr->estimateNonLeafPosterior();
	gammaMgr->estimateTransitionProbability();
	return 0;
}

double 
SpeciesClusterManager::getScore()
{
	double totalComplexity=0;
	//SK: count number of degrees of freedom for extant species.
	int ndfExtant=0;
	for(map<string,CLUSTERSET*>::iterator cIter=speciesExpertSet.begin();cIter!=speciesExpertSet.end();cIter++)
	{
		//SK: count number of data measurements per gene in this analysis 
		int NData=0;
		GeneExpManager* speciesMgr=speciesExprSet[cIter->first];
		CLUSTERSET* expertSet=cIter->second;
		map<string,int>* speciesGenes=speciesClusterSet_Genewise[cIter->first];
		int paramCnt=(expertSet->size()*2*(speciesGenes->size()));
		double modelComplexity=(paramCnt/2)*log(speciesGenes->size());
		totalComplexity=totalComplexity+modelComplexity;
		for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
		{
			eIter->second->resetClip();
		}
		for(map<string,int>::iterator vIter=speciesGenes->begin();vIter!=speciesGenes->end();vIter++)
		{
			int ogid=mor->getMappedOrthogroupID(vIter->first.c_str(),cIter->first.c_str());
			//SK: check if the SpeciesClusterManager object is in CVMode and if the represented OGID is not in the training set or in the working set
			if(CVMode==true && (trainingSet.find(ogid)==trainingSet.end()|| workingOrthoGroups.find(ogid)==workingOrthoGroups.end()))
			{
				//SK: if the OGID is not in the training or working set, continue
				continue;
			}
			vector<double>* exprProf=speciesMgr->getExp(vIter->first);
			//SK: set number of data points
			if(NData==0)
			{
				NData=exprProf->size();
			}
			map<int,double>* mixOutProbs=gammaMgr->getLeafLikelihood_store(ogid,(string&)vIter->first);
			for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
			{
				Expert* e=eIter->second;
				double pdf=e->getOutputPDF(exprProf);
				if(isnan(pdf))
				{
					cout <<"PDF is nan for " << cIter->first << " " << vIter->first << " for expert " << eIter->first << endl;
				}	
				(*mixOutProbs)[eIter->first]=pdf;
			}
		}
		for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
		{
			cout <<cIter->first<<":" << eIter->first <<" " << eIter->second->getClip()<< " minpdf " << eIter->second->getMinPDF() << endl;;
		}
		//SK: Add to the total number of degrees of freedom for this species
                ndfExtant+=(NData*maxClusterCnt)+(maxClusterCnt*(maxClusterCnt-1));
	}
	//get number of extant species
	int NSpecies=speciesExpertSet.size();
        //SK: Add number of degrees of freedom for correction, for most species
        double ndfAnc=(NSpecies-1)*(maxClusterCnt*(maxClusterCnt-1));
        //SK: and consider for LCA root species with k-1 parameters.
        ndfAnc += maxClusterCnt-1;
        //SK: Count the number of OGIDS to consider in this instance
        int NOrthos=0;
        //SK: loop over test OGIDS
	//SK: if not in CVMOde
	if(CVMode==false)
	{
		//SK: simply get the number of working ortho groups. 
		NOrthos=workingOrthoGroups.size();
	}
	//SK: if in CV mode
        if(CVMode==true)
	{
		//SK: loop over training set OGIDS
		for(map<int,int>::iterator itr=trainingSet.begin();itr!=trainingSet.end();itr++)
        	{
                	//check if it is one of the working orthos.
                	if(workingOrthoGroups.find(itr->first)!=workingOrthoGroups.end())
                	{
                        	NOrthos+=1;
                	}
        	}
	}
        //SK: calculate full correction for the penalized likelihood score.
        cout << ndfExtant << "\t" << ndfAnc << "\t" << NOrthos << "\t" << log(NOrthos) << endl;
	double penaltyBIC = (double) ((double )ndfExtant+ndfAnc)*((double) 1/2)*log(NOrthos);
	double penaltyAIC = (double) ((double )ndfExtant+ndfAnc)*((double) 1/2)*2;
	cout << "PenaltyBIC: " << penaltyBIC << endl;
	double netLL_unpen=gammaMgr->getAllNodeScore();
	double netLL=netLL_unpen+totalComplexity;
	//SK: set corrected score
        double netLL_correctedBIC=netLL_unpen-penaltyBIC;
	double netLL_correctedAIC=netLL_unpen-penaltyAIC;
        //print out all three score values
        cout <<"Unpenalized score= " << netLL_unpen << " Penalized score=" << netLL << " Corrected BIC score=" << netLL_correctedBIC << "Corrected AIC score=" << netLL_correctedAIC << endl;
        //SK: add score values to score vector
	scores.clear();
        scores.push_back(netLL_unpen);
        scores.push_back(netLL);
        scores.push_back(netLL_correctedBIC);
	scores.push_back(netLL_correctedAIC);
	return netLL_unpen;
}

int
SpeciesClusterManager::estimateMeanCov(Expert* e, string& specName, int clusterID)
{
	GeneExpManager* exprMgr=speciesExprSet[specName];
	map<string,int>* speciesGenes=speciesClusterSet_Genewise[specName];
	vector<double>* expr=exprMgr->getExp(speciesGenes->begin()->first);
	//cout << speciesGenes->begin()->first << endl;
	CLUSTERSET* expertSet=speciesExpertSet[specName];
	int dim=expr->size();
	//cout << dim << endl;
	Matrix* mean=e->getMean();
	if(mean==NULL)
	{
		mean= new Matrix(1,dim);
	}
	double totaldp=0;
	for(int d=0;d<dim;d++)
	{
		double sum=0;
		double meanval=0;
		for(map<string,int>::iterator dIter=speciesGenes->begin();dIter!=speciesGenes->end();dIter++)
		{
			//Assuming that each evidence is actually a joing assignment to RVs. What we
			//want is a mean and variance for each experiment
			vector<double>* geneExpr=exprMgr->getExp(dIter->first);
			if(geneExpr==NULL)
			{
				cout <<"No gene by name " << dIter->first << " in species " << specName << endl;
			}
			int ogid=mor->getMappedOrthogroupID(dIter->first.c_str(),specName.c_str());
			//SK: check if operating in the CV mode and if the OGID is not in the training set
			if(CVMode==true && (trainingSet.find(ogid)==trainingSet.end() || workingOrthoGroups.find(ogid)==workingOrthoGroups.end()))
			{
				//SK: continue if the orthogorup is not in the training and working set
				continue;
			}
			//This gamma matrix is a kXk matrix which stores the joint probability of being in cluster i
			//given its ancestor was in some other species. 
			Matrix* gamma_i_k_s=gammaMgr->getGamma(ogid,(string&)dIter->first,specName);
			//Need to sum over all possible ways in which we could end up in clusterID
			for(int r=0;r<gamma_i_k_s->getRowCnt();r++)
			{
				double g_i=gamma_i_k_s->getValue(r,clusterID);
				if(g_i>1)
				{
					cout <<"Weird gamma found for " << dIter->first << " at row " << r << endl;
				}
				sum=sum+g_i;
				if(sum>10000)
				{
					cout <<"Weird sum " << sum << " found at " << dIter->first << " at row " << r << endl;
				}

				double eval=(*geneExpr)[d];
				meanval=meanval+(g_i*eval);
			}
		}
		if(sum==0)
		{
			cout <<"No members in cluster " << clusterID << endl;
		}
		meanval=meanval/sum;
		mean->setValue(meanval,0,d);
		totaldp=sum;
	}
	Matrix* covariance=e->getCovariance();
	if(covariance==NULL)
	{
		covariance=new Matrix(dim,dim);
	}
	covariance->setAllValues(0);
	for(int i=0;i<dim;i++)
	{
		double m1=mean->getValue(0,i);
		//for(int j=i;j<dim;j++)
		for(int j=i;j<i+1;j++)
		{
			double m2=mean->getValue(0,j);
			double cov=0;
			if(i==j)
			{
				cov=0.00001;
			}
			double sum=0;
			for(map<string,int>::iterator dIter=speciesGenes->begin();dIter!=speciesGenes->end();dIter++)
			{
				vector<double>* geneExpr=exprMgr->getExp(dIter->first);
				int ogid=mor->getMappedOrthogroupID(dIter->first.c_str(),specName.c_str());
				//SK: check if CV mode is active and if the represented in OGIDS is not in the training set
				if(CVMode==true && (trainingSet.find(ogid)==trainingSet.end() || workingOrthoGroups.find(ogid)==workingOrthoGroups.end()))
				{
					//SK: continue if the orthogroups is not in the training and working set
					continue;
				}
				Matrix* gamma_i_k_s=gammaMgr->getGamma(ogid,(string&)dIter->first,specName);
				//Need to sum over all possible ways in which we could end up in clusterID
				for(int r=0;r<gamma_i_k_s->getRowCnt();r++)
				{
					double g_i=gamma_i_k_s->getValue(r,clusterID);
					sum=sum+g_i;
					double diff=g_i*((*geneExpr)[i]-m1)*((*geneExpr)[j]-m2);
					cov=cov+diff;
				}
			}
			cov=cov/sum;
			covariance->setValue(cov,i,j);
			covariance->setValue(cov,j,i);
		}
	}
	cout <<"Mean estimated from "<< totaldp << " in " << specName <<":" << clusterID<< endl;
//	mean->showMatrix();
	if(e->getMean()==NULL)
	{
		e->setMean(mean);
	}
	cout <<"Covariance estimated from " << totaldp << " in " << specName << ":" << clusterID << endl;
//	covariance->showMatrix();
	if(e->getCovariance()==NULL)
	{
		e->setCovariance(covariance);
	}
	else
	{
		e->updateCovariance();
	}
	//cout <<"Estimated params for " << specName << ":"<< clusterID << " prior="<< prior << endl;
	return totaldp;
}

int
SpeciesClusterManager::expectationStep_Species(string& specName, CLUSTERSET* expertSet)
{
	GeneExpManager* exprMgr=speciesExprSet[specName];
	map<string,int>* speciesGenes=speciesClusterSet_Genewise[specName];
	for(map<string,int>::iterator vIter=speciesGenes->begin();vIter!=speciesGenes->end();vIter++)
	{
		int ogid=mor->getMappedOrthogroupID(vIter->first.c_str(),specName.c_str());
		//SK: check if operating in CV mode, and if the represented OGID is in the training set and a working orthogroups
		if(CVMode==true && (trainingSet.find(ogid)==trainingSet.end() || workingOrthoGroups.find(ogid)==workingOrthoGroups.end()))
		{
			//SK: continue if the OGID is not in the training and working set
                        continue;
                }
		vector<double>* exprProf=exprMgr->getExp(vIter->first);
		map<int,double> mixOutProbs;
		double sum=0;
		if(strcmp(vIter->first.c_str(),"orf19.5809")==0 || strcmp(vIter->first.c_str(),"Q0080")==0)
		{
			cout << "Found " << vIter->first << endl;
		}
		for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
		{
			Expert* e=eIter->second;
			double pdf=e->getOutputPDF(exprProf);
			if(isnan(pdf))
			{
				cout <<"PDF is nan for " << specName << " " << vIter->first << " for expert " << eIter->first << endl;
			}
			mixOutProbs[eIter->first]=pdf;
			sum=sum+pdf;
		}
		gammaMgr->estimateLeafAlpha(ogid,mixOutProbs,(string&)vIter->first,specName);
		mixOutProbs.clear();
	}
	return 0;
}
//SK: function for applying the expectation step to the test OGIDS data, as part of the CV mode, much like the expectationStep function, and is called in the getScore_Test function
int
SpeciesClusterManager::expectationStep_Test()
{
        //First get the gammas for each leaf node
        for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
                CLUSTERSET* speciesClusters=sIter->second;
		//SK: call the expectationStep_Test_Species function for the individual species.
                expectationStep_Test_Species((string&)sIter->first,speciesClusters);
        }
        //The use the speciesdist manager's conditional distributions to infer the rest of the gammas.
	//SK: calculate the non-leaf posterior probabilities for the test data 
        gammaMgrTst->estimateNonLeafPosterior();
	//SK: calculate the transition Probability for the test data
        gammaMgrTst->estimateTransitionProbability();
        return 0;
}

//SK: function called in expectationStep_Test as part of getScore_Test
int
SpeciesClusterManager::expectationStep_Test_Species(string& specName, CLUSTERSET* expertSet)
{
	//SK: call the expression data for the given species
        GeneExpManager* exprMgr=speciesExprSet[specName];
        /*map<string,int>* speciesGenes=speciesClusterSet_Genewise[specName];
        for(map<string,int>::iterator vIter=speciesGenes->begin();vIter!=speciesGenes->end();vIter++)
	{
		int ogid=mor->getMappedOrthogroupID(vIter->first.c_str(),specName.c_str());//Get the associated OGID for the given gene
                if(CVMode==true && testSet.find(ogid)==testSet.end())//SK: check if the gene is in in the test OGIDS
                {
                        continue;
                }
                vector<double>* exprProf=exprMgr->getExp(vIter->first);
                map<int,double> mixOutProbs;
                double sum=0;
                if(strcmp(vIter->first.c_str(),"orf19.5809")==0 || strcmp(vIter->first.c_str(),"Q0080")==0)
		{
                        cout << "Found " << vIter->first << endl;
                }
                for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
		{
                        Expert* e=eIter->second;
                        double pdf=e->getOutputPDF(exprProf);
                        if(isnan(pdf))
			{
                        	cout <<"PDF is nan for " << specName << " " << vIter->first << " for expert " << eIter->first << endl;
                        }
                        mixOutProbs[eIter->first]=pdf;
                        sum=sum+pdf;//SK: add pdf variable to sum
                }
                gammaMgrTst->estimateLeafGamma(ogid,mixOutProbs,(string&)vIter->first,specName);
                mixOutProbs.clear();
        }*/
	//SK: call mapped orthogroup set
	map<int,MappedOrthogroup*>& OMap= mor->getMappedOrthogroups();
	//SK: loop over all test set OGIDS
	for(map<int,int>::iterator itr=testSet.begin();itr!=testSet.end();itr++)
	{
		//SK: init variable for the OGID for clarity in the lines below
		int ogid=itr->first;
		//make sure the orthogroup is in the set of groups being worked. 
		if(workingOrthoGroups.find(ogid)==workingOrthoGroups.end())
		{
			continue;
		}
		//SK: define grp variable for the associated MappedOrthogroup object
		MappedOrthogroup* grp = OMap.find(ogid)->second;
		//SK: call the gene map for this species from the MappedOrthogroup onject
		GeneMap* geneMap=grp->getSpeciesHits(specName.c_str());
		//SK: check of the GeneMap is null and ortherwise continue
                if(geneMap==NULL)
                {
			continue;
			//cout << "Warning: GeneMap Null in expectationStep_Test_Species" << endl;
                }
		//SK: otherwise loop over the species genes in this orthogroup
                else
                {
			//SK: get gene map, based on lines near ~847 below in this file
                	map<string,map<string,STRINTMAP*>*>& geneset=geneMap->getGeneSet();
			//SK: loop over gene set
                        for(map<string,map<string,STRINTMAP*>*>::iterator sIter=geneset.begin();sIter!=geneset.end();sIter++)
                        {
				//SK: print OGID and gene name for sanity check
                        	//cout << ogid << "\t" << sIter->first << endl;
				//SK: obtain the expression profile for this gene
				vector<double>* exprProf=exprMgr->getExp(sIter->first);
				//SK: ckeck if the expression profile is null and otherwise continue
                		if(exprProf==NULL)
				{
					continue;
				}
				//SK: define map for the mix out probabiity values
				map<int,double> mixOutProbs;
				//SK: define variable for sum of "pdf" probabilities
                		double sum=0;
				//SK: loop over the expert set for the clusters in this species
				for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
				{
					//SK: set e variable to an Expert object for a given cluster
					Expert* e=eIter->second;
					//SK: get the ouput PDF probability for this Expert (cluster) based on the expression profile of this gene
					double pdf=e->getOutputPDF(exprProf);
					//SK: check if this pdf probability is a nan value
					if(isnan(pdf))
					{
						cout <<"PDF is nan for " << specName << " " << sIter->first << " for expert " << eIter->first << endl;
					}
					//SK: set mixOutProbs variable for this cluster 
					mixOutProbs[eIter->first]=pdf;
					//SK: add pdf variable to sum
					sum=sum+pdf;
				}
				//SK: call the gamma manager object to estimate the leafe node gamma matrix with the map of the mix out probabilities 
				gammaMgrTst->estimateLeafGamma(ogid,mixOutProbs,(string&)sIter->first,specName);
				//SK: clear the mixOutProbs map object
				mixOutProbs.clear();	
			}
                }

	}
        return 0;
}

//SK: this is the penultimate function to call for caluculating the score for the test data
double
SpeciesClusterManager::getScore_Test()
{
	//SK: check if CVMode is set and if there are test OGIDS.
	if(!CVMode || testSet.size()==0)
	{
		return 0;
	}
	//SK: check if CVMode has been initialized
	if(CVMode==false)
	{
		//SK: given warning statement 
		cout << "Cannot calculate test data score because the crossvalidation analysis is not initialized; SpeciesClusterManager::exiting getScore_Test()" << endl;
		//exit function.
		return 0;
	}
	//SK: Initialized gamma manager for test data in readClusters in startup.
	//SK: apply expectation step.
	expectationStep_Test();
	//SK: now calculate score for test data
 	//SK: define the varible for calculating complexity of the model on the test data set
	double totalComplexity=0;
	//SK: Count the number of data degrees of freedom in the extand species.
	int ndfExtant=0;
	//SK count number of orthogroups included
	
	//SK: the following line commented out below are for calculating the score from the original getScore function, which this function is emulating.
        /*for(map<string,CLUSTERSET*>::iterator cIter=speciesExpertSet.begin();cIter!=speciesExpertSet.end();cIter++)
	{
                GeneExpManager* speciesMgr=speciesExprSet[cIter->first];
                CLUSTERSET* expertSet=cIter->second;
                //map<string,int>* speciesGenes=speciesClusterSet_Genewise[cIter->first];
                int paramCnt=(expertSet->size()*2*(speciesGenes->size()));
                double modelComplexity=(paramCnt/2)*log(speciesGenes->size());
                totalComplexity=totalComplexity+modelComplexity;
                for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
		{
                        eIter->second->resetClip();
                }
                for(map<string,int>::iterator vIter=speciesGenes->begin();vIter!=speciesGenes->end();vIter++)
		{
                        int ogid=mor->getMappedOrthogroupID(vIter->first.c_str(),cIter->first.c_str());
			if(CVMode==true && testSet.find(ogid)==testSet.end())
			{
				continue;
			}
                        vector<double>* exprProf=speciesMgr->getExp(vIter->first);
                        map<int,double>* mixOutProbs=gammaMgrTst->getLeafLikelihood_store(ogid,(string&)vIter->first);
                        for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
			{
                                Expert* e=eIter->second;
                                double pdf=e->getOutputPDF(exprProf);
                                if(isnan(pdf))
				{
                                        cout <<"PDF is nan for " << cIter->first << " " << vIter->first << " for expert " << eIter->first << endl;
                                }
                                (*mixOutProbs)[eIter->first]=pdf;
                        }
                }
                for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
		{
                        cout <<cIter->first<<":" << eIter->first <<" " << eIter->second->getClip()<< " minpdf " << eIter->second->getMinPDF() << endl;
                }
        }*/
	//SK: access the mapp of the ortho group objects
	map<int,MappedOrthogroup*>& OMap= mor->getMappedOrthogroups();
	//SK: first loop over the species represented in the analysis via the speciesExpertSet map object
	for(map<string,CLUSTERSET*>::iterator cIter=speciesExpertSet.begin();cIter!=speciesExpertSet.end();cIter++)
        {
		//SK: define variable for number of genes in the species
		int NSpeciesGenes=0;
		//SK: count the number of measurements for genes in this species.
		int NData=0;
		//SK: call species expression data map
		GeneExpManager* speciesMgr=speciesExprSet[cIter->first];
		//SK: call expert set for species
                CLUSTERSET* expertSet=cIter->second;		
		//SK: reset expert sets as in original example
		for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
                {
                        eIter->second->resetClip();
                }
		//SK: loop over the test OGIDS to loop over the species  
        	for(map<int,int>::iterator itr=testSet.begin();itr!=testSet.end();itr++)
        	{
			//SK: set ogid varible for clarity below
                	int ogid=itr->first;
			//SK: check that this orthogroup is in the set of working orthogroups that were initialized in the GammaManager objects in readClusters
			if(workingOrthoGroups.find(ogid)==workingOrthoGroups.end()){
				//SK: skip this orthogroup because it isn't associated with genes that were initialized with cluster assignments in readClusters
				continue;
			}
			//SK: set the grp varible to the right Mapped Orthogroup object from the set in OMap
                	MappedOrthogroup* grp = OMap.find(ogid)->second;
			//SK: get the geneMap for this orthogroup
                	GeneMap* geneMap=grp->getSpeciesHits(cIter->first.c_str());
			//SK: check that the GeneMap object isn't null
                	if(geneMap==NULL)
                	{
                        	continue;
                	}
			//SK: otherwise continue to loop over the genes for this species
                	else
                	{
				//SK: get the geneset for the relevant species 
                        	map<string,map<string,STRINTMAP*>*>& geneset=geneMap->getGeneSet();
				//loop over the gene set
                        	for(map<string,map<string,STRINTMAP*>*>::iterator sIter=geneset.begin();sIter!=geneset.end();sIter++)
                        	{
					//SK: print ogid and gene name for sanity check
                                	//cout << ogid << sIter->first << endl;
					//SK: obtain the expression profile for this gene
                                	vector<double>* exprProf=speciesMgr->getExp(sIter->first);
                                	//SK: ckeck if the expression profile is null and continue if it is
                                	if(exprProf==NULL)
                                	{
                                        	continue;
                                	}
					//SK: increment the number of genes being counted in this species by one
                                        NSpeciesGenes+=1;
					//SK: set number of data entries for each gene
					if(NData==0)
					{
						NData=exprProf->size();
					}
					//SK: get mixOutProbs for this gene
					map<int,double>* mixOutProbs=gammaMgrTst->getLeafLikelihood_store(ogid,(string&)sIter->first);
					//SK: loop over the Expert objects, as modeled above
                        		for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
                        		{
						//SK: take Expert object and define the variable for the pdf value
                                		Expert* e=eIter->second;
                                		double pdf=e->getOutputPDF(exprProf);
						//SK: check if that probability variable is set to  "nan"
                                		if(isnan(pdf))
                                		{
							//SK: print out warning about "nan" pdf value
                                        		cout <<"PDF is nan for " << cIter->first << " " << cIter->first << " for expert " << eIter->first << endl;
                                		}
						//SK: add the pdf value to the mixOutProbs map
                                		(*mixOutProbs)[eIter->first]=pdf;
					//SK: end loop over Experts
                        		}
				//SK: end loop over the gene set	
				}
				//SK: loo over the expert set one more time
                		for(map<int,Expert*>::iterator eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
                		{
					//SK: print out the pdf values
                        		cout <<cIter->first<<":" << eIter->first <<" " << eIter->second->getClip()<< " minpdf " << eIter->second->getMinPDF() << endl;
				//SK: end loop over Expert objects
                		}
			//SK: complete else statement if the GeneMap is not NUll, end loop over genes or this orthogroup id in orther words
			}
		//SK: end the loop over OGIDs for orthogroups in the test set
		}
		//SK: calculate total complexity of the model for this species, first count the number of parameters
		int paramCnt=(expertSet->size()*2*NSpeciesGenes); //(speciesGenes->size()));
		//SK: calculate the model complexity score
                double modelComplexity=(paramCnt/2)*log(NSpeciesGenes); //speciesGenes->size());
		//SK: add it to the total complexity score
                totalComplexity=totalComplexity+modelComplexity;
		//SK: Add to the total number of degrees of freedom for this species
		ndfExtant+=(NData*maxClusterCnt)+(maxClusterCnt*(maxClusterCnt-1));
	}
	//SK: get total score for test data from the GammaManager for the test set
        double netLL_unpen=gammaMgrTst->getAllNodeScore();
	//SK: Add model score and complexity score
        double netLL=netLL_unpen+totalComplexity;
	//SK: make further correction
	//SK: get number of extant species. 
	int NSpecies=speciesExpertSet.size();
	//SK: Add number of degrees of freedom for correction, for most species
	double ndfAnc=(NSpecies-1)*(maxClusterCnt*(maxClusterCnt-1));
	//SK: and consider for LCA root species with k-1 parameters.
	ndfAnc += maxClusterCnt-1;
	//SK: Count the number of OGIDS to consider in this instance
	int NOrthos=0;
	//SK: loop over test OGIDS
	for(map<int,int>::iterator itr=testSet.begin();itr!=testSet.end();itr++)
	{
		//check if it is one of the working orthos.
		if(workingOrthoGroups.find(itr->first)!=workingOrthoGroups.end())
		{
			NOrthos+=1;
		}
	}
	//SK: calculate full correction for the penalized likelihood score.
	cout << ndfExtant << "\t" << ndfAnc << "\t" << NOrthos << "\t" << log(NOrthos) << endl;
	int ndf=ndfExtant+ndfAnc;
	double hold=(double) ndf / 2;
	double penaltyBIC = (double) hold*log(NOrthos);
	double penaltyAIC = (double) hold*2;
	cout << "Penalty BIC: " << penaltyBIC << endl;
	//SK: set corrected score
	double netLL_correctedBIC=netLL_unpen-penaltyBIC;
	double netLL_correctedAIC=netLL_unpen-penaltyAIC;
	//print out values
        cout <<"Unpenalized score= " << netLL_unpen << " Penalized score=" << netLL << " Corrected BIC score=" << netLL_correctedBIC << "Corrected AIC score=" << netLL_correctedAIC << endl;
	//SK: add score values to score vector
	testScores.push_back(netLL_unpen);
	testScores.push_back(netLL);
	testScores.push_back(netLL_correctedBIC);
	testScores.push_back(netLL_correctedAIC);
	//SK: complete the function and return the score value
        return netLL_unpen;
}
int
SpeciesClusterManager::showClusters_Extant(const char* outputDir)
{
	char aFName[1024];
	sprintf(aFName,"%s/rseed.txt",outputDir);
	ofstream oFile(aFName);
	oFile << rseed << endl;
	oFile.close();
	assignGenesToExperts_FromMap();
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		char outFName[1024];
		sprintf(outFName,"%s/%s_exprtab.txt",outputDir,sIter->first.c_str());
		displaySpeciesClusters(outFName,(string&)sIter->first);
		sprintf(outFName,"%s/%s_clusterassign.txt",outputDir,sIter->first.c_str());
		displaySpeciesClusterAssignments(outFName,(string&)sIter->first);
		sprintf(outFName,"%s/%s_speciesspecnames_clusterassign.txt",outputDir,sIter->first.c_str());
		displaySpeciesClusterAssignments_NonScerNames(outFName,(string&)sIter->first);
	}
	return 0;
}

int
SpeciesClusterManager::showMeans(const char* outputDir)
{
	char outFName[1024];
	sprintf(outFName,"%s/clustermeans.txt",outputDir);
	ofstream oFile(outFName);
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		CLUSTERSET* clustering=sIter->second;
		for(CLUSTERSET_ITER cIter=clustering->begin();cIter!=clustering->end();cIter++)
		{
			oFile <<sIter->first<<"_"<<cIter->first;
			Expert* expert=cIter->second;
			Matrix* m=expert->getMean();
			for(int c=0;c<m->getColCnt();c++)
			{
				oFile <<"\t" << m->getValue(0,c);
			}
			oFile << endl;
		}
	}
	oFile.close();
	return 0;
}


int
SpeciesClusterManager::showMeans(const char* outputDir, int iter)
{
	char outFName[1024];
	sprintf(outFName,"%s/clustermeans_%d.txt",outputDir,iter);
	ofstream oFile(outFName);
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		CLUSTERSET* clustering=sIter->second;
		for(CLUSTERSET_ITER cIter=clustering->begin();cIter!=clustering->end();cIter++)
		{
			oFile <<sIter->first<<"_"<<cIter->first;
			Expert* expert=cIter->second;
			Matrix* m=expert->getMean();
			for(int c=0;c<m->getColCnt();c++)
			{
				oFile <<"\t" << m->getValue(0,c);
			}
			oFile << endl;
		}
	}
	oFile.close();
	return 0;
}

int
SpeciesClusterManager::showClusters(const char* outputDir)
{
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		char outFName[1024];
		sprintf(outFName,"%s/%s_exprtab.txt",outputDir,sIter->first.c_str());
		displaySpeciesClusters(outFName,(string&)sIter->first);
		sprintf(outFName,"%s/%s_clusterassign.txt",outputDir,sIter->first.c_str());
		displaySpeciesClusterAssignments(outFName,(string&)sIter->first);
	}
	return 0;
}




int
SpeciesClusterManager::showClusters_Ancestral(const char* outputDir)
{
	map<string,ofstream*> filePtrs;
	char geneName[32];
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	/*for(map<int,map<string,int>*>::iterator ogIter=mappedClusterAssignment.begin();ogIter!=mappedClusterAssignment.end();ogIter++)
	{
		map<string,int>* clusterAssign=ogIter->second;
		MappedOrthogroup* mgrp=allOgs[ogIter->first];
		for(map<string,int>::iterator cIter=clusterAssign->begin();cIter!=clusterAssign->end();cIter++)
		{
			strcpy(geneName,cIter->first.c_str());
			char* pos=strchr(geneName,':');
			if(pos==NULL)
			{
				cout <<"Bad format " << endl;	
				exit(0);
			}
			*pos='\0';
			string specName(pos+1);
			if(speciesExpertSet.find(specName)!=speciesExpertSet.end())
			{	
				continue;
			}
			ofstream* oFile=NULL;
			if(filePtrs.find(specName)==filePtrs.end())
			{	
				char aFName[1024];
				sprintf(aFName,"%s/%s_clusterassign.txt",outputDir,specName.c_str());
				oFile=new ofstream(aFName);
				filePtrs[specName]=oFile;
			}
			else
			{	
				oFile=filePtrs[specName];
			}
			GeneMap* geneMap=mgrp->getSpeciesHits(srcSpecies);
			if(geneMap==NULL)
			{
				(*oFile)<< "OG"<<ogIter->first<<":"<< geneName <<"\t" << cIter->second << endl;
			}
			else
			{
				//Get the scer ortholog of this gene
				map<string,map<string,STRINTMAP*>*>& geneset=geneMap->getGeneSet();
				//Fear here is that the same scer gene is repeated within and between clusters
				for(map<string,map<string,STRINTMAP*>*>::iterator sIter=geneset.begin();sIter!=geneset.end();sIter++)
				{
					(*oFile) <<sIter->first <<"\t" << cIter->second << endl;
				}
			}
		}
	}*/

	for(map<int,map<string,int>*>::iterator ogIter=allClusterAssignmentsGrouped.begin();ogIter!=allClusterAssignmentsGrouped.end();ogIter++)
	{
		int ogid=ogIter->first;
		
		map<string,int>* clusterAssign=ogIter->second;
		int groupID=0;
		map<int,map<string,int>*> geneAssignSet;
		map<int,map<string,string>*> geneNamesSet;
		char geneName[1024];
		for(map<string,int>::iterator cIter=clusterAssign->begin();cIter!=clusterAssign->end();cIter++)
		{
			char tempBuff[1024];
			strcpy(tempBuff,cIter->first.c_str());
			char* tok=strtok(tempBuff,":");
			int tokCnt=0;
			string speciesName;
			string aName;
			int dup=-1;
			while(tok!=NULL)
			{
				if(tokCnt==0)
				{
					dup=atoi(tok);
				}
				else if(tokCnt==1)
				{
					speciesName.append(tok);
				}	
				else
				{
					aName.append(tok);
				}
				tok=strtok(NULL,":");
				tokCnt++;
			}
			map<string,int>* geneSet=NULL;
			map<string,string>* geneName=NULL;
			if(geneAssignSet.find(dup)==geneAssignSet.end())
			{
				geneSet=new map<string,int>;
				geneName=new map<string,string>;
				geneAssignSet[dup]=geneSet;
				geneNamesSet[dup]=geneName;	
			}
			else
			{
				geneSet=geneAssignSet[dup];
				geneName=geneNamesSet[dup];
			}
			(*geneSet)[speciesName]=cIter->second;
			(*geneName)[speciesName]=aName;
		}
		map<string,int> localKey;
		map<string,int> localKeyDup;
		bool toShowOG=false;
		for(map<int,map<string,int>*>::iterator hIter=geneAssignSet.begin();hIter!=geneAssignSet.end();hIter++)
		{
			map<string,string>* gname=geneNamesSet[hIter->first];
			map<string,int>* gset=hIter->second;
			int ancClusterAssign=-3;
			//Now use the scername to show this og
			string dispName;
			if(gset->find(srcSpecies)==gset->end())
			{
				char tempName[1024];	
				sprintf(tempName,"OG%d_%d",ogid,hIter->first);
				dispName.append(tempName);
			}	
			else
			{
				string geneName((*gname)[srcSpecies]);
				//dispName.append(gnm.getCommonName((*gname)[srcSpecies].c_str()));
				dispName.append(geneName);
			}

			for(map<string,int>::iterator sIter=gset->begin();sIter!=gset->end();sIter++)
			{
				int geneclusterID=-2;
				string specName(sIter->first);
				geneclusterID=sIter->second;
				if(geneclusterID<0)
				{
					continue;
				}
				ofstream* oFile=NULL;
				if(filePtrs.find(specName)==filePtrs.end())
				{	
					char aFName[1024];
					sprintf(aFName,"%s/%s_clusterassign.txt",outputDir,specName.c_str());
					oFile=new ofstream(aFName);
					filePtrs[specName]=oFile;
				}
				else
				{	
					oFile=filePtrs[specName];
				}
				(*oFile) << dispName << "\t"<< geneclusterID << endl;;
			}
		}
		
	}
	for(map<string,ofstream*>::iterator fIter=filePtrs.begin();fIter!=filePtrs.end();fIter++)
	{
		fIter->second->close();
	}
	return 0;
}




/*
int
SpeciesClusterManager::dumpAllInferredClusters_Scerwise_Dup(const char* outputDir,vector<string>& speciesList)
{	
	char outputFName[1024];
	sprintf(outputFName,"%s/allspecies_duplicates_clusterassign_brk.txt",outputDir);
	string specKey(srcSpecies);
	map<string,int>* speciesGenes=speciesClusterSet_Genewise[specKey];
	ofstream oFile(outputFName);
	CLUSTERSET* expertSet=speciesExpertSet[specKey];
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	map<string,int> shownGenes;
	for(CLUSTERSET_ITER cIter=expertSet->begin();cIter!=expertSet->end();cIter++)
	{
		Expert* e=cIter->second;
		map<string,int>& geneSet=e->getGeneSet();
		for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			if(shownGenes.find(gIter->first)!=shownGenes.end())
			{
				continue;
			}
			MappedOrthogroup* mgrp=mor->getMappedOrthogroup(gIter->first.c_str(),srcSpecies);
			if(mgrp->getCnt()<2)
			{
				continue;
			}
			map<string,int>* clusterAssign=mappedClusterAssignment[ogid];
			GeneMap* dupspecies_GeneMap=mgrp->getSpeciesHits(srcName);
			if(dupspecies_GeneMap==NULL)
			{
				cout <<"Warning! No " << srcSpecies << " gene in OGid "<<  << endl;
				continue;
			}
			bool dupInScer=false;
			//get the number of species that have 2 genes
			map<string,int> speciesWithDupgenes;
			for(map<string,int>::iterator aIter=speciesList.begin();aIter!=speciesList.end();aIter++)
			{
				GeneMap* specRep=mgrp->getSpeciesHits(aIter->first.c_str());
				if(specRep==NULL)
				{
					continue;
				}	
				if(specRep->getGeneSet().size()>=2)
				{
					speciesWithDupGenes[aIter->first]=0;
				}
			}
			if(speciesWithDupGenes.find(srcSpecies)!=speciesWithDupGenes.end())
			{
				dupInScer=true;
			}
			string dupSpecies(srcName);
			if(!dupInScer)
			{
				dupSpecies.clear();
				dupSpecies.append(speciesWithDupgenes.begin()->first.c_str());
				dupspecies_GeneMap=mgrp->getSpeciesHits(dupSpecies.c_str());
			}
			for(map<string,int>::iterator cIter=clusterAssign->begin();cIter!=clusterAssign->end();cIter++)
			{
				strcpy(geneName,cIter->first.c_str());
				char* pos=strchr(geneName,':');
				if(pos==NULL)
				{
					cout <<"Bad format " << endl;	
					exit(0);
				}
				*pos='\0';
				string specName(pos+1);
				geneAssign[geneName]=cIter->second;
			}
			map<string,map<string,STRINTMAP*>*>& dupGenes=dupSpecies_GeneMap->getGeneSet();
			for(map<string,map<string,STRINTMAP*>*>::iterator hIter=dupGenes.begin();hIter!=dupGenes.end();hIter++)
			{
				oFile << gnm.getCommonName(hIter->first.c_str());
				shownGenes[hIter->first]=0;
				map<string,int> specAssign;
				int clusterID=-1;
				if(geneAssign.find(hIter->first)!=geneAssign.end())
				{
					clusterID=geneAssign[hIter->first];
				}
				specAssign[srcSpecies]=clusterID;
				map<string,STRINTMAP*>* hitsInOthers=hIter->second;
				for(int s=0;s<speciesList.size();s++)
				{
					if(strcmp(speciesLis[s].c_str(),"Scer")==0)
					{
						continue;
					}
					//This means the gene is altogether missing in the species
					int clusterID=-2;
					if(hitsInOthers->find(speciesList[s])==hitsInOthers->end())
					{
						clusterID=-2;
					}
					else
					{
						map<string,int>* genesInOthers=(*hitsInOthers)[speciesList[s]];
						if(genesInOthers->size()>1)
						{
							cout <<"Warning! Multiple hits to " << hIter->first << " from " << speciesList[s] << endl;
						}
						string& orthogene=genesInOthers->begin();
						if(geneAssign.find(orthogene)!=geneAssignAssign.end())
						{	
							clusterID=geneAssign[orthogene];
						}	
						else
						{
							clusterID=-1;
						}
						shownGenes[orthogene]=1;
					}
					specAssign[speciesList[s]]=clusterID;
				}
				for(map<string,int>::iterator sIter=specAssign.begin();sIter!=specAssign.end();sIter++)
				{
					oFile <<"\t" << sIter->second;
				}
				oFile << endl;	
				specAssign.clear();
			}
			
			geneAssign.clear();
		}
		oFile <<"Dummy" << cIter->first;
		for(int s=0;s<speciesList.size();s++)
		{
			oFile <<"\t" <<"-100";
		}
		oFile <<endl;
	}
	oFile.close();
	return 0;
}*/


int
SpeciesClusterManager::displaySpeciesClusterAssignments(const char* outFName,string& specName)
{
	map<string,int>* speciesGenes=speciesClusterSet_Genewise[specName];
	ofstream oFile(outFName);
	CLUSTERSET* expertSet=speciesExpertSet[specName];
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	for(CLUSTERSET_ITER cIter=expertSet->begin();cIter!=expertSet->end();cIter++)
	{
		Expert* e=cIter->second;
		map<string,int>& geneSet=e->getGeneSet();
		for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			MappedOrthogroup* mgrp=mor->getMappedOrthogroup(gIter->first.c_str(),specName.c_str());
			const string& genename=gIter->first;
			if(strcmp(specName.c_str(),srcSpecies)==0)
			{
				oFile<< genename <<"\t" << cIter->first << endl;
			}
			else
			{
				//Get the scer ortholog of this gene
				map<string,int>* scerOrtholog=mor->getOrtholog(specName.c_str(),genename.c_str(),srcSpecies);
				if(scerOrtholog==NULL)
				{
					oFile <<"OG"<< mgrp->getID() << ":"<< genename<<"\t"<< cIter->first << endl;
				}
				else
				{
					for(map<string,int>::iterator sIter=scerOrtholog->begin();sIter!=scerOrtholog->end();sIter++)
					{
						oFile <<sIter->first<<"\t"<< cIter->first << endl;
					}
				}
			}
		}
	}
	oFile.close();
	return 0;
}

//Ugh if only I knew how to name my functions better! The main differenc here and the displaySpeciesClusterAssignments
//is that here we display the original species-specific name of the species
int
SpeciesClusterManager::displaySpeciesClusterAssignments_NonScerNames(const char* outFName,string& specName)
{
	map<string,int>* speciesGenes=speciesClusterSet_Genewise[specName];
	ofstream oFile(outFName);
	CLUSTERSET* expertSet=speciesExpertSet[specName];
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	for(CLUSTERSET_ITER cIter=expertSet->begin();cIter!=expertSet->end();cIter++)
	{
		Expert* e=cIter->second;
		map<string,int>& geneSet=e->getGeneSet();
		for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			oFile<< gIter->first<<"\t" << cIter->first << endl;
		}
	}
	oFile.close();
	return 0;
}


int
SpeciesClusterManager::assignGenesToExperts()
{
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		CLUSTERSET* expertSet=sIter->second;
		for(CLUSTERSET_ITER eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
		{
			eIter->second->resetAssignedGenes();
		}
	}
	for(map<string,GeneExpManager*>::iterator sIter=speciesExprSet.begin();sIter!=speciesExprSet.end();sIter++)
	{
		map<string,int>* speciesGenes=speciesClusterSet_Genewise[sIter->first];
		CLUSTERSET* expertSet=speciesExpertSet[sIter->first];	
		for(map<string,int>::iterator gIter=speciesGenes->begin();gIter!=speciesGenes->end();gIter++)
		{
			int ogid=mor->getMappedOrthogroupID(gIter->first.c_str(),sIter->first.c_str());
			//This gamma matrix is a kXk matrix which stores the joint probability of being in cluster i
			//given its ancestor was in some other species. 
			Matrix* normterm=gammaMgr->getNormTerm(ogid,(string&)gIter->first,(string&)sIter->first);
			int maxClusterID=-1;	
			double maxProb=0;
			for(int c=0;c<normterm->getColCnt();c++)
			{
				double prob=normterm->getValue(0,c);
				if(prob> maxProb)
				{
					maxProb=prob;
					maxClusterID=c;
				}
			}
			if(maxClusterID==-1)
			{
				cout <<"Gene " << gIter->first << " cannot be assigned to any cluster of " << sIter->first.c_str() << endl;
				continue;
			}
			Expert* e=(*expertSet)[maxClusterID];
			e->assignGeneToExpert(gIter->first.c_str());
		}	
	}
	return 0;
}


int
SpeciesClusterManager::assignGenesToExperts_FromMap()
{
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		CLUSTERSET* expertSet=sIter->second;
		for(CLUSTERSET_ITER eIter=expertSet->begin();eIter!=expertSet->end();eIter++)
		{
			eIter->second->resetAssignedGenes();
		}
	}
	char geneName[1024];
	for(map<int,map<string,int>*>::iterator ogIter=mappedClusterAssignment.begin();ogIter!=mappedClusterAssignment.end();ogIter++)
	{
		map<string,int>* geneAssignments=ogIter->second;
		int ogid=ogIter->first;
		for(map<string,int>::iterator gIter=geneAssignments->begin();gIter!=geneAssignments->end();gIter++)
		{
			strcpy(geneName,gIter->first.c_str());
			char* pos=strchr(geneName,':');
			if(pos==NULL)
			{
				cout <<"Bad format" << endl;
				exit(0);
			}
			if(ogid==4597)
			{
				cout << gIter->first << "\t" << gIter->second << endl;
			}
			*pos='\0';
			if(strcmp(geneName,"orf19.5809")==0)
			{	
				cout << "Stop here " << geneName << " " << gIter->first << " " << gIter->second<< endl;
			}
			string specName(pos+1);
			if(speciesExpertSet.find(specName)==speciesExpertSet.end())
			{
				continue;
			}
			CLUSTERSET* expertSet=speciesExpertSet[specName];
			if(expertSet->find(gIter->second)==expertSet->end())
			{
				cout << "No cluster " << gIter->second<< " for " << gIter->first  << endl;
				continue;
			}
			Expert* e=(*expertSet)[gIter->second];
			e->assignGeneToExpert(geneName);
		}	
	}
	return 0;
}



int
SpeciesClusterManager::displaySpeciesClusters(const char* outFName,string& specName)
{
	map<string,int>* speciesGenes=speciesClusterSet_Genewise[specName];
	GeneExpManager* expMgr=speciesExprSet[specName];
	ofstream oFile(outFName);
	CLUSTERSET* expertSet=speciesExpertSet[specName];
	vector<double>* expr=expMgr->getExp(speciesGenes->begin()->first);
	int dim=expr->size();
	oFile <<"Gene";
	for(int i=0;i<dim;i++)
	{
		oFile <<"\tExp"<<i;
	}
	oFile << endl;
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	for(int clusterID=0;clusterID<maxClusterCnt;clusterID++)
	{
		if(expertSet->find(clusterID)==expertSet->end())
		{
			oFile <<"Dummy" << clusterID;
			for(int i=0;i<dim;i++)
			{
				oFile <<"\t" <<"-100";
			}
			oFile <<endl;
			continue;
		}
		Expert* e=(*expertSet)[clusterID];
		map<string,int>& geneSet=e->getGeneSet();
		for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			if(strcmp(gIter->first.c_str(),"YDR398W")==0)
			{
				cout <<"Stop here found "<< gIter->first << endl;
			}
			oFile<< gIter->first;
			vector<double>* expr=expMgr->getExp(gIter->first);
			if(expr==NULL)
			{
				cout <<"No expression for " << gIter->first << endl;
				continue;
			}
			for(int i=0;i<expr->size();i++)
			{
				oFile << "\t"<<(*expr)[i];
			}
			oFile << endl;
		}
		oFile <<"Dummy" << clusterID;
		for(int i=0;i<dim;i++)
		{
			oFile <<"\t" <<"-100";
		}
		oFile <<endl;
	}
	oFile.close();
	return 0;
}


int
SpeciesClusterManager::dumpAllInferredClusters_Scerwise(const char* outputDir,vector<string>& speciesList)
{	
	char outputFName[1024];
	sprintf(outputFName,"%s/allspecies_clusterassign_brk.txt",outputDir);
	ofstream oFile(outputFName);
	sprintf(outputFName,"%s/allspecies_duplicates_clusterassign_brk.txt",outputDir);
	ofstream dFile(outputFName);
	string specKey(srcSpecies);
	map<string,int>* speciesGenes=speciesClusterSet_Genewise[specKey];
	CLUSTERSET* expertSet=speciesExpertSet[specKey];
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	//for(CLUSTERSET_ITER cIter=expertSet->begin();cIter!=expertSet->end();cIter++)
	map<string,int> shownGenes;
	for(int clusterID=0;clusterID<maxClusterCnt;clusterID++)
	{
		if(expertSet->find(clusterID)==expertSet->end())
		{
			oFile <<"Dummy" << clusterID;
			dFile <<"Dummy" << clusterID;
			for(int s=0;s<speciesList.size();s++)
			{
				oFile <<"\t" <<"-100";
				dFile <<"\t" <<"-100";
			}
			oFile <<endl;
			dFile <<endl;
			continue;
		}
		Expert* e=(*expertSet)[clusterID];
		map<string,int>& geneSet=e->getGeneSet();
		for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			if(shownGenes.find(gIter->first)!=shownGenes.end())
			{	
				continue;
			}
			MappedOrthogroup* mgrp=mor->getMappedOrthogroup(gIter->first.c_str(),srcSpecies);
			int ogid=mgrp->getID();
			if(ogid==3865)
			{
				cout <<"OG " << ogid << endl;
			}
			map<string,int>* clusterAssign=mappedClusterAssignment[ogid];
			map<string,int> geneAssign;
			char geneName[1024];
			for(map<string,int>::iterator cIter=clusterAssign->begin();cIter!=clusterAssign->end();cIter++)
			{
				strcpy(geneName,cIter->first.c_str());
				char* pos=strchr(geneName,':');
				if(pos==NULL)
				{
					cout <<"Bad format " << endl;	
					exit(0);
				}
				*pos='\0';
				string specName(pos+1);
				string geneNameKey(geneName);
				if(strstr(specName.c_str(),"Anc")==NULL)
				{
					geneAssign[geneNameKey]=cIter->second;
				}
				else	
				{
					geneAssign[geneNameKey]=cIter->second;
				}
			}
			const char* dupAnc=gammaMgr->getDupAncestor(ogid);
			bool wgd=false;
			if(dupAnc!=NULL && strcmp(dupAnc,"Anc5")==0)
			{
				wgd=true;
			}
			map<int,map<string,string>*>& geneSets=mgrp->getGeneSets();
			map<string,int> localKey;
			for(map<int,map<string,string>*>::iterator hIter=geneSets.begin();hIter!=geneSets.end();hIter++)
			{
				map<string,string>* gset=hIter->second;
				if(gset->find(srcSpecies)==gset->end())
				{
					oFile << "OG" << ogid <<"_" << hIter->first;
					if(wgd)
					{
						oFile <<"*";
					}
					if(geneSets.size()>1)
					{
						dFile << "OG" << ogid <<"_" << hIter->first;
						if(wgd)
						{
							dFile <<"*";
						}
					}
				}	
				else
				{
					string dispName(gnm.getCommonName((*gset)[srcSpecies].c_str()));
					if(localKey.find(dispName)!=localKey.end())
					{
						dispName.append("_2");
					}
					localKey[dispName]=0;
					oFile << dispName;
					if(wgd)
					{
						oFile <<"*";
					}
					if(geneSets.size()>1)
					{
						dFile << dispName;
						if(wgd)
						{
							dFile <<"*";
						}
					}
				}
				for(int s=0;s<speciesList.size();s++)
				{
					int geneclusterID=-2;
					if(strstr(speciesList[s].c_str(),"Anc")!=NULL)
					{
						char temp[1024];
						sprintf(temp,"%s_%d",speciesList[s].c_str(),hIter->first+1);
						string tempKey(temp);
						//if(geneAssign.find(speciesList[s])!=geneAssign.end())
						if(geneAssign.find(tempKey)!=geneAssign.end())
						{
							geneclusterID=geneAssign[tempKey];
						}
								
					}
					else
					{
						if(gset->find(speciesList[s])==gset->end())
						{
							geneclusterID=-2;
						}
						else 
						{	
							string& geneName=(*gset)[speciesList[s]];
							if(geneAssign.find(geneName)!=geneAssign.end())
							{
								geneclusterID=geneAssign[geneName];
							}
							shownGenes[geneName]=0;
						}
					}
					oFile<< "\t"<< geneclusterID;
					if(geneSets.size()>1)
					{
						dFile<< "\t"<< geneclusterID;
					}
				}
				oFile << endl;	
				if(geneSets.size()>1)
				{
					dFile<< endl;
				}
			}
			localKey.clear();
		}
		oFile <<"Dummy" << clusterID;
		dFile <<"Dummy" << clusterID;
		for(int s=0;s<speciesList.size();s++)
		{
			oFile <<"\t" <<"-100";
			dFile <<"\t" <<"-100";
		}
		oFile <<endl;
		dFile <<endl;
	}
	oFile.close();
	dFile.close();
	return 0;
}


int
SpeciesClusterManager::dumpAllInferredClusters_ScerwiseGrouped(const char* outputDir,vector<string>& speciesList)
{	
	char outputFName[1024];
	sprintf(outputFName,"%s/allspecies_clusterassign_brk.txt",outputDir);
	ofstream oFile(outputFName);
	sprintf(outputFName,"%s/allspecies_duplicates_clusterassign_brk.txt",outputDir);
	ofstream dFile(outputFName);
	string specKey(srcSpecies);
	map<string,int>* speciesGenes=speciesClusterSet_Genewise[specKey];
	CLUSTERSET* expertSet=speciesExpertSet[specKey];
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	map<string,int> shownGenes;
	gammaMgr->getAllClusterAssignments_Grouped(allClusterAssignmentsGrouped);
	sprintf(outputFName,"%s/allspecies_clusterassign_tree.txt",outputDir);
	gammaMgr->showAssignments(outputFName);
	map<int,string> ogDupSpeciesMap;
	for(int clusterID=0;clusterID<maxClusterCnt;clusterID++)
	{
		if(expertSet->find(clusterID)==expertSet->end())
		{
			oFile <<"Dummy" << clusterID;
			dFile <<"Dummy" << clusterID;
			for(int s=0;s<speciesList.size();s++)
			{
				oFile <<"\t" <<"-100";
				dFile <<"\t" <<"-100";
			}
			oFile <<endl;
			dFile <<endl;
			continue;
		}
		Expert* e=(*expertSet)[clusterID];
		map<string,int>& geneSet=e->getGeneSet();
		for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			if(shownGenes.find(gIter->first)!=shownGenes.end())
			{	
				continue;
			}
			MappedOrthogroup* mgrp=mor->getMappedOrthogroup(gIter->first.c_str(),srcSpecies);
			int ogid=mgrp->getID();
			map<string,int>* clusterAssign=allClusterAssignmentsGrouped[ogid];
			map<string,int>* clusterAssignOrig=mappedClusterAssignment[ogid];
			if(ogid==4597)
			{
				cout <<"OG " << ogid << endl;
				for(map<string,int>::iterator cIter=clusterAssign->begin();cIter!=clusterAssign->end();cIter++)
				{
					cout << cIter->first << " " << cIter->second << endl;
				}
				cout <<"Original assignment" << endl;
				for(map<string,int>::iterator cIter=clusterAssignOrig->begin();cIter!=clusterAssignOrig->end();cIter++)
				{
					cout << cIter->first << " " << cIter->second << endl;
				}
			}

			int groupID=0;
			map<int,map<string,int>*> geneAssignSet;
			map<int,map<string,string>*> geneNamesSet;
			char geneName[1024];
			for(map<string,int>::iterator cIter=clusterAssign->begin();cIter!=clusterAssign->end();cIter++)
			{
				char tempBuff[1024];
				strcpy(tempBuff,cIter->first.c_str());
				char* tok=strtok(tempBuff,":");
				int tokCnt=0;
				string speciesName;
				string aName;
				int dup=-1;
				while(tok!=NULL)
				{
					if(tokCnt==0)
					{
						dup=atoi(tok);
					}
					else if(tokCnt==1)
					{
						speciesName.append(tok);
					}	
					else
					{
						aName.append(tok);
					}
					tok=strtok(NULL,":");
					tokCnt++;
				}
				map<string,int>* geneSet=NULL;
				map<string,string>* geneName=NULL;
				if(geneAssignSet.find(dup)==geneAssignSet.end())
				{
					geneSet=new map<string,int>;
					geneName=new map<string,string>;
					geneAssignSet[dup]=geneSet;
					geneNamesSet[dup]=geneName;	
				}
				else
				{
					geneSet=geneAssignSet[dup];
					geneName=geneNamesSet[dup];
				}
				(*geneSet)[speciesName]=cIter->second;
				(*geneName)[speciesName]=aName;
				char buffer[1024];
				string key;
				if(strstr(cIter->first.c_str(),"Anc")!=NULL)
				{
					sprintf(buffer,"%s_%d:%s",speciesName.c_str(),dup,speciesName.c_str());
				}	
				else
				{
					sprintf(buffer,"%s:%s",aName.c_str(),speciesName.c_str());
				}
				key.append(buffer);
				int otherassign=-1;
				if(clusterAssignOrig->find(key)==clusterAssignOrig->end())
				{
					cout << " No cluster for "<< ogid << " " << key.c_str() << endl;
				}
				else
				{
					otherassign=(*clusterAssignOrig)[key];
					if(otherassign!=cIter->second)
					{
						cout <<"Cluster asssignment mismatch for  " << ogid<<  " " << key.c_str() << " OLD " << otherassign << " " << cIter->second<< endl;
					}
				}
			}
			const char* dupAnc=gammaMgr->getDupAncestor(ogid);
			bool wgd=false;
			if(dupAnc!=NULL)
			{
				string dupAncStr(dupAnc);
				ogDupSpeciesMap[ogid]=dupAncStr;	
			}
			if(dupAnc!=NULL && strcmp(dupAnc,"Anc5")==0)
			{
				wgd=true;
			}
			map<int,map<string,string>*>& geneSets=mgrp->getGeneSets();
			if(geneSets.size()!=geneAssignSet.size())
			{
				cout <<"Missed a duplicate for  " << ogid << endl;
			}
			map<string,int> localKey;
			map<string,int> localKeyDup;
			for(map<int,map<string,int>*>::iterator hIter=geneAssignSet.begin();hIter!=geneAssignSet.end();hIter++)
			{
				map<string,int>* gset=hIter->second;
				map<string,string>* gname=geneNamesSet[hIter->first];
				int srcClusterAssign=-3;
				if(gset->find(srcSpecies)!=gset->end())
				{
					srcClusterAssign=(*gset)[srcSpecies];
				}
				
				if(gset->find(srcSpecies)==gset->end())
				{
					if(srcClusterAssign==clusterID)
					{
						oFile << "OG" << ogid <<"_" << hIter->first;
						if(wgd)
						{
							oFile <<"*";
						}
					}
					//if(geneSets.size()>1)
					if(geneAssignSet.size()>1)
					{
						dFile << "OG" << ogid <<"_" << hIter->first;
						if(wgd)
						{
							dFile <<"*";
						}
					}
				}	
				else
				{
					if(srcClusterAssign==clusterID)
					{
						string& geneName=(*gname)[srcSpecies];
						shownGenes[geneName]=0;
						string dispName(gnm.getCommonName((*gname)[srcSpecies].c_str()));
						if(localKey.find(dispName)!=localKey.end())
						{
							dispName.append("_2");
						}
						oFile << dispName;
						if(wgd)
						{
							oFile <<"*";
						}
						localKey[dispName]=0;
					}
					//if(geneSets.size()>1)
					if(geneAssignSet.size()>1)
					{
						string dispName(gnm.getCommonName((*gname)[srcSpecies].c_str()));
						if(localKeyDup.find(dispName)!=localKeyDup.end())
						{
							dispName.append("_2");
						}
						dFile << dispName;
						if(wgd)
						{
							dFile <<"*";
						}
						localKeyDup[dispName]=0;
					}
				}

				for(int s=0;s<speciesList.size();s++)
				{
					int geneclusterID=-2;
					if(gset->find(speciesList[s])==gset->end())
					{
						geneclusterID=-2;
					}
					else 
					{	
						string& geneName=(*gname)[speciesList[s]];
						geneclusterID=(*gset)[speciesList[s]];
				//		shownGenes[geneName]=0;
					}
					if(srcClusterAssign==clusterID)
					{
						oFile<< "\t"<< geneclusterID;
					}
					//if(geneSets.size()>1)
					if(geneAssignSet.size()>1)
					{
						dFile<< "\t"<< geneclusterID;
					}
				}
				if(srcClusterAssign==clusterID)
				{
					oFile << endl;	
				}
				//if(geneSets.size()>1)
				if(geneAssignSet.size()>1)
				{
					dFile<< endl;
				}
			}
		}
		oFile <<"Dummy" << clusterID;
		dFile <<"Dummy" << clusterID;
		for(int s=0;s<speciesList.size();s++)
		{
			oFile <<"\t" <<"-100";
			dFile <<"\t" <<"-100";
		}
		oFile <<endl;
		dFile <<endl;
	}
	oFile.close();
	dFile.close();
	sprintf(outputFName,"%s/dupinfo_inferred.txt",outputDir);
	ofstream dupFile(outputFName);
	for(map<int,string>::iterator ogIter=ogDupSpeciesMap.begin();ogIter!=ogDupSpeciesMap.end();ogIter++)
	{
		dupFile << ogIter->first <<"\t" << ogIter->second << endl;
	}
	dupFile.close();
	return 0;
}


int
SpeciesClusterManager::dumpAllInferredClusters_LCA(const char* outputDir,vector<string>& speciesList, string& lcaName)
{	
	char outputFName[1024];
	sprintf(outputFName,"%s/allspecies_clusterassign_lca_brk.txt",outputDir);
	ofstream oFile(outputFName);
	sprintf(outputFName,"%s/allspecies_duplicates_clusterassign__lcabrk.txt",outputDir);
	ofstream dFile(outputFName);
	sprintf(outputFName,"%s/genemembers_perog.txt",outputDir);
	ofstream geneMembersFile(outputFName);
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	map<string,int> shownGenes;
	map<int,string> ogDupSpeciesMap;
	int badOG_DupLCA=0;
	int badOG_NoLCA=0;
	map<int,int> shownOGs;
	for(int clusterID=0;clusterID<maxClusterCnt;clusterID++)
	{
		//Iterate over the allClusterAssignments but show only those groups for which the Anc14 cluster assingment matches
		//the clusterID
		for(map<int,map<string,int>*>::iterator ogIter=allClusterAssignmentsGrouped.begin();ogIter!=allClusterAssignmentsGrouped.end();ogIter++)
		{
			int ogid=ogIter->first;
			
			map<string,int>* clusterAssign=ogIter->second;
			map<string,int>* clusterAssignOrig=mappedClusterAssignment[ogid];
			int groupID=0;
			map<int,map<string,int>*> geneAssignSet;
			map<int,map<string,string>*> geneNamesSet;
			char geneName[1024];
			for(map<string,int>::iterator cIter=clusterAssign->begin();cIter!=clusterAssign->end();cIter++)
			{
				char tempBuff[1024];
				strcpy(tempBuff,cIter->first.c_str());
				char* tok=strtok(tempBuff,":");
				int tokCnt=0;
				string speciesName;
				string aName;
				int dup=-1;
				while(tok!=NULL)
				{
					if(tokCnt==0)
					{
						dup=atoi(tok);
					}
					else if(tokCnt==1)
					{
						speciesName.append(tok);
					}	
					else
					{
						aName.append(tok);
					}
					tok=strtok(NULL,":");
					tokCnt++;
				}
				map<string,int>* geneSet=NULL;
				map<string,string>* geneName=NULL;
				if(geneAssignSet.find(dup)==geneAssignSet.end())
				{
					geneSet=new map<string,int>;
					geneName=new map<string,string>;
					geneAssignSet[dup]=geneSet;
					geneNamesSet[dup]=geneName;	
				}
				else
				{
					geneSet=geneAssignSet[dup];
					geneName=geneNamesSet[dup];
				}
				(*geneSet)[speciesName]=cIter->second;
				(*geneName)[speciesName]=aName;
				char buffer[1024];
				string key;
				if(strstr(cIter->first.c_str(),"Anc")!=NULL)
				{
					sprintf(buffer,"%s_%d:%s",speciesName.c_str(),dup,speciesName.c_str());
				}	
				else
				{
					sprintf(buffer,"%s:%s",aName.c_str(),speciesName.c_str());
				}
				key.append(buffer);
				int otherassign=-1;
				if(clusterAssignOrig->find(key)==clusterAssignOrig->end())
				{
					cout << " No cluster for "<< ogid << " " << key.c_str() << endl;
				}
				else
				{
					otherassign=(*clusterAssignOrig)[key];
					if(otherassign!=cIter->second)
					{
						cout <<"Cluster asssignment mismatch for  " << ogid<<  " " << key.c_str() << " OLD " << otherassign << " " << cIter->second<< endl;
					}
				}
			}
			const char* dupAnc=gammaMgr->getDupAncestor(ogid);
			bool wgd=false;
			if(dupAnc!=NULL)
			{
				string dupAncStr(dupAnc);
				ogDupSpeciesMap[ogid]=dupAncStr;	
			}
			if(dupAnc!=NULL && strcmp(dupAnc,"Anc5")==0)
			{
				wgd=true;
			}
			map<string,int> localKey;
			map<string,int> localKeyDup;
			//Also check if Anc14 has more than 1 copy
			//Show the orthogroup if and only if one of the copies has the cluster ID as current cluster ID.
			bool toShowOG=false;
			for(map<int,map<string,int>*>::iterator hIter=geneAssignSet.begin();hIter!=geneAssignSet.end();hIter++)
			{
				map<string,int>* gset=hIter->second;
				int ancClusterAssign=-3;
				if(gset->find(lcaName)!=gset->end())
				{
					ancClusterAssign=(*gset)[lcaName];
					if(ancClusterAssign==clusterID)
					{
						toShowOG=true;
					}
				}
				
			}
			if(!toShowOG)
			{
				continue;
			}
			shownOGs[ogid]=0;
			int lcacnt=0;
			for(map<int,map<string,int>*>::iterator hIter=geneAssignSet.begin();hIter!=geneAssignSet.end();hIter++)
			{
				map<string,string>* gname=geneNamesSet[hIter->first];
				map<string,int>* gset=hIter->second;
				int ancClusterAssign=-3;
				if(gset->find(lcaName)!=gset->end())
				{
					ancClusterAssign=(*gset)[lcaName];
					lcacnt++;
				}
				//Now use the scername to show this og
				string dispName;
				if(gset->find(srcSpecies)==gset->end())
				{
					char tempName[1024];	
					sprintf(tempName,"OG%d_%d",ogid,hIter->first);
					dispName.append(tempName);
					if(wgd)
					{
						dispName.append("*");
					}
				}	
				else
				{
					string geneName((*gname)[srcSpecies]);
					shownGenes[geneName]=0;
					dispName.append(gnm.getCommonName((*gname)[srcSpecies].c_str()));
					if(localKey.find(dispName)!=localKey.end())
					{
						dispName.append("_2");
					}
					localKey[dispName]=0;
					if(wgd)
					{
						dispName.append("*");
					}
				}
				//if(ancClusterAssign==clusterID)
				//{
					oFile <<dispName;
					if(geneAssignSet.size()>1)
					{
						dFile << dispName;
					}
					//Show the gene members
				//}
				string members;

				for(int s=0;s<speciesList.size();s++)
				{
					int geneclusterID=-2;
					if(gset->find(speciesList[s])==gset->end())
					{
						geneclusterID=-2;
					}
					else 
					{	
						geneclusterID=(*gset)[speciesList[s]];
					}
					if(members.length()>0)
					{
						members.append(";");
					}	
					if(gname->find(speciesList[s])==gname->end())
					{
						members.append("NULL");
					}
					else
					{
						members.append((*gname)[speciesList[s]]);
					}		
					//if(ancClusterAssign==clusterID)
				//	{
						oFile<< "\t"<< geneclusterID;
				//	}
					if(geneAssignSet.size()>1)
					{
						dFile<< "\t"<< geneclusterID;
					}
				}
				geneMembersFile <<"OG"<<ogid <<"_" << hIter->first <<"\t" << members << endl;
				//if(ancClusterAssign==clusterID)
				//{
					oFile << endl;	
				//}
				if(geneAssignSet.size()>1)
				{
					dFile<< endl;
				}
			}
			if(lcacnt>1)
			{
				badOG_DupLCA++;
				cout <<"LCA has multiple copies in "<< ogid << endl;
			}
		}
		oFile <<"Dummy" << clusterID;
		dFile <<"Dummy" << clusterID;
		for(int s=0;s<speciesList.size();s++)
		{
			oFile <<"\t" <<"-100";
			dFile <<"\t" <<"-100";
		}
		oFile <<endl;
		dFile <<endl;
	}
	oFile.close();
	dFile.close();
	geneMembersFile.close();
	sprintf(outputFName,"%s/dupinfo_inferred.txt",outputDir);
	ofstream dupFile(outputFName);
	for(map<int,string>::iterator ogIter=ogDupSpeciesMap.begin();ogIter!=ogDupSpeciesMap.end();ogIter++)
	{
		dupFile << ogIter->first <<"\t" << ogIter->second << endl;
	}
	dupFile.close();
	cout <<" BadOG: Dup_LCA_OGs " << badOG_DupLCA  << endl;
	cout <<"Shown " << shownOGs.size () << " out of " << allClusterAssignmentsGrouped.size() <<  " total OGs"  << endl;
	return 0;
}




/*
int
SpeciesClusterManager::dumpAllInferredClusters_Scerwise_Dup(const char* outputDir,vector<string>& speciesList)
{	
	char outputFName[1024];
	sprintf(outputFName,"%s/allspecies_duplicates_clusterassign_brk.txt",outputDir);
	string specKey(srcSpecies);
	map<string,int>* speciesGenes=speciesClusterSet_Genewise[specKey];
	ofstream oFile(outputFName);
	CLUSTERSET* expertSet=speciesExpertSet[specKey];
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	map<string,int> shownGenes;
	for(CLUSTERSET_ITER cIter=expertSet->begin();cIter!=expertSet->end();cIter++)
	{
		Expert* e=cIter->second;
		map<string,int>& geneSet=e->getGeneSet();
		for(map<string,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			if(shownGenes.find(gIter->first)!=shownGenes.end())
			{
				continue;
			}
			MappedOrthogroup* mgrp=mor->getMappedOrthogroup(gIter->first.c_str(),srcSpecies);
			if(mgrp->getCnt()<2)
			{
				continue;
			}
			map<string,int>* clusterAssign=mappedClusterAssignment[ogid];
			GeneMap* dupspecies_GeneMap=mgrp->getSpeciesHits(srcName);
			if(dupspecies_GeneMap==NULL)
			{
				cout <<"Warning! No " << srcSpecies << " gene in OGid "<<  << endl;
				continue;
			}
			bool dupInScer=false;
			//get the number of species that have 2 genes
			map<string,int> speciesWithDupgenes;
			for(map<string,int>::iterator aIter=speciesList.begin();aIter!=speciesList.end();aIter++)
			{
				GeneMap* specRep=mgrp->getSpeciesHits(aIter->first.c_str());
				if(specRep==NULL)
				{
					continue;
				}	
				if(specRep->getGeneSet().size()>=2)
				{
					speciesWithDupGenes[aIter->first]=0;
				}
			}
			if(speciesWithDupGenes.find(srcSpecies)!=speciesWithDupGenes.end())
			{
				dupInScer=true;
			}
			string dupSpecies(srcName);
			if(!dupInScer)
			{
				dupSpecies.clear();
				dupSpecies.append(speciesWithDupgenes.begin()->first.c_str());
				dupspecies_GeneMap=mgrp->getSpeciesHits(dupSpecies.c_str());
			}
			for(map<string,int>::iterator cIter=clusterAssign->begin();cIter!=clusterAssign->end();cIter++)
			{
				strcpy(geneName,cIter->first.c_str());
				char* pos=strchr(geneName,':');
				if(pos==NULL)
				{
					cout <<"Bad format " << endl;	
					exit(0);
				}
				*pos='\0';
				string specName(pos+1);
				geneAssign[geneName]=cIter->second;
			}
			map<string,map<string,STRINTMAP*>*>& dupGenes=dupSpecies_GeneMap->getGeneSet();
			for(map<string,map<string,STRINTMAP*>*>::iterator hIter=dupGenes.begin();hIter!=dupGenes.end();hIter++)
			{
				oFile << gnm.getCommonName(hIter->first.c_str());
				shownGenes[hIter->first]=0;
				map<string,int> specAssign;
				int clusterID=-1;
				if(geneAssign.find(hIter->first)!=geneAssign.end())
				{
					clusterID=geneAssign[hIter->first];
				}
				specAssign[srcSpecies]=clusterID;
				map<string,STRINTMAP*>* hitsInOthers=hIter->second;
				for(int s=0;s<speciesList.size();s++)
				{
					if(strcmp(speciesLis[s].c_str(),"Scer")==0)
					{
						continue;
					}
					//This means the gene is altogether missing in the species
					int clusterID=-2;
					if(hitsInOthers->find(speciesList[s])==hitsInOthers->end())
					{
						clusterID=-2;
					}
					else
					{
						map<string,int>* genesInOthers=(*hitsInOthers)[speciesList[s]];
						if(genesInOthers->size()>1)
						{
							cout <<"Warning! Multiple hits to " << hIter->first << " from " << speciesList[s] << endl;
						}
						string& orthogene=genesInOthers->begin();
						if(geneAssign.find(orthogene)!=geneAssignAssign.end())
						{	
							clusterID=geneAssign[orthogene];
						}	
						else
						{
							clusterID=-1;
						}
						shownGenes[orthogene]=1;
					}
					specAssign[speciesList[s]]=clusterID;
				}
				for(map<string,int>::iterator sIter=specAssign.begin();sIter!=specAssign.end();sIter++)
				{
					oFile <<"\t" << sIter->second;
				}
				oFile << endl;	
				specAssign.clear();
			}
			
			geneAssign.clear();
		}
		oFile <<"Dummy" << cIter->first;
		for(int s=0;s<speciesList.size();s++)
		{
			oFile <<"\t" <<"-100";
		}
		oFile <<endl;
	}
	oFile.close();
	return 0;
}*/




int
SpeciesClusterManager::dumpAllInferredClusterAssignments(const char* outputDir)
{
	char outputFName[1024];
	sprintf(outputFName,"%s/clusterassign_multspecies.txt",outputDir);
	ofstream oFile(outputFName);
	//gammaMgr->getAllClusterAssignments_Conditional(mappedClusterAssignment);
	gammaMgr->getAllClusterAssignments(mappedClusterAssignment,secondStage);
	gammaMgr->showClusterFlipCnts();
	gammaMgr->reestimateTransitionProbability();
	char outputDirCmd[1024];
	sprintf(outputDirCmd,"mkdir %s/prior_pp/",outputDir);
	system(outputDirCmd);
	char outputSubDir[1024];
	sprintf(outputSubDir,"%s/prior_pp",outputDir);
	showClusters_Extant(outputSubDir);
	int iter=0;
	double currScore=0;
	bool convergence=false;
	//while(iter<50 && !convergence)
	while(secondStage && iter<nMaxIterations && !convergence)
	{
		showMeans(outputDir,iter);
		if(iter>0)
		{
			expectationStep();
		}
		maximizationStep();
		double newScore=getScore();
		double diff=fabs(newScore-currScore);
		if((iter>0) && (diff<convThresh))
		{
			convergence=true;
		}
		cout <<"PPIter: " << iter << " score: " << newScore << " diff " << diff << endl;
		currScore=newScore;
		iter++;
	}
	gammaMgr->getAllClusterAssignments(mappedClusterAssignment,false);
	//gammaMgr->showClusterFlipCnts();
	//gammaMgr->reestimateTransitionProbability();
	//gammaMgr->getAllClusterAssignments_Conditional(mappedClusterAssignment);
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	for(map<int,map<string,int>*>::iterator ogIter=mappedClusterAssignment.begin();ogIter!=mappedClusterAssignment.end();ogIter++)
	{
		//oFile<< ogIter->first;
		MappedOrthogroup* mgrp=allOgs[ogIter->first];
		GeneMap* geneMap=mgrp->getSpeciesHits(srcSpecies);
		if(geneMap==NULL)
		{
			oFile <<" --";
		}
		else
		{
			map<string,map<string,STRINTMAP*>*>& geneSet=geneMap->getGeneSet();
			for(map<string,map<string,STRINTMAP*>*>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
			{
				if(gIter==geneSet.begin())
				{
					oFile << gIter->first;
				}
				else
				{
					oFile <<","<< gIter->first;
				}
			}
		}
		map<string,int>* clusterAss=ogIter->second;
		for(map<string,int>::iterator cIter=clusterAss->begin();cIter!=clusterAss->end();cIter++)
		{
			const char* pos=strchr(cIter->first.c_str(),':');
			if(pos==NULL)
			{
				cout <<"Bad format " << endl;	
				exit(0);
			}
			string specName(pos+1);
			oFile << "\t"<<specName<<"=" << cIter->second;
		}
		oFile << endl;
	}
	oFile.close();
	return 0;
}


int
SpeciesClusterManager::dumpAllInferredClusterAssignments(const char* outputDir,int clusterID)
{
	char outputFName[1024];
	sprintf(outputFName,"%s/clusterassign_multspecies_%d.txt",outputDir,clusterID);
	ofstream oFile(outputFName);
	gammaMgr->getAllClusterAssignments(mappedClusterAssignment,true);
	//gammaMgr->getAllClusterAssignments_Conditional(mappedClusterAssignment);
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	for(map<int,map<string,int>*>::iterator ogIter=mappedClusterAssignment.begin();ogIter!=mappedClusterAssignment.end();ogIter++)
	{
		//oFile<< ogIter->first;
		MappedOrthogroup* mgrp=allOgs[ogIter->first];
		GeneMap* geneMap=mgrp->getSpeciesHits(srcSpecies);
		if(geneMap==NULL)
		{
			oFile <<" --";
		}
		else
		{
			map<string,map<string,STRINTMAP*>*>& geneSet=geneMap->getGeneSet();
			for(map<string,map<string,STRINTMAP*>*>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
			{
				if(gIter==geneSet.begin())
				{
					oFile << gIter->first;
				}
				else
				{
					oFile <<","<< gIter->first;
				}
			}
		}
		map<string,int>* clusterAss=ogIter->second;
		for(map<string,int>::iterator cIter=clusterAss->begin();cIter!=clusterAss->end();cIter++)
		{
			const char* pos=strchr(cIter->first.c_str(),':');
			if(pos==NULL)
			{
				cout <<"Bad format " << endl;	
				exit(0);
			}
			string specName(pos+1);
			oFile << "\t"<<specName<<"=" << cIter->second;
		}
		oFile << endl;
	}
	oFile.close();
	return 0;
}

//Show all the cluster assignments per ancestral species too
int
SpeciesClusterManager::dumpAllInferredClusterGammas(const char* outputDir,vector<string>& speciesList)
{
	char outputFName[1024];
	gammaMgr->getAllClusterGammas(mappedClusterGamma);
	sprintf(outputFName,"%s/clustergamma_multspecies.txt",outputDir);
	ofstream oFile(outputFName);
	map<int,MappedOrthogroup*>& allOgs=mor->getMappedOrthogroups();
	for(map<int,map<string,map<int,double>*>*>::iterator ogIter=mappedClusterGamma.begin();ogIter!=mappedClusterGamma.end();ogIter++)
	{
		//oFile<< ogIter->first;
		MappedOrthogroup* mgrp=allOgs[ogIter->first];
		int ogid=mgrp->getID();
		map<string,map<int,double>*> geneAssign;
		map<string,map<int,double>*>* clusterAssign=mappedClusterGamma[ogid];
		char geneName[1024];
		for(map<string,map<int,double>*>::iterator cIter=clusterAssign->begin();cIter!=clusterAssign->end();cIter++)
		{
			strcpy(geneName,cIter->first.c_str());
			char* pos=strchr(geneName,':');
			if(pos==NULL)
			{
				cout <<"Bad format " << endl;	
				exit(0);
			}
			*pos='\0';
			string specName(pos+1);
			string geneNameKey(geneName);
			if(strstr(specName.c_str(),"Anc")==NULL)
			{
				geneAssign[geneNameKey]=cIter->second;
			}
			else	
			{
				geneAssign[specName]=cIter->second;
			}
		}
		
		map<int,map<string,string>*>& geneSets=mgrp->getGeneSets();
		for(map<int,map<string,string>*>::iterator hIter=geneSets.begin();hIter!=geneSets.end();hIter++)
		{
			map<string,string>* gset=hIter->second;
			if(gset->find(srcSpecies)==gset->end())
			{
				continue;
			}
			oFile << gnm.getCommonName((*gset)[srcSpecies].c_str());
			for(int s=0;s<speciesList.size();s++)
			{
				map<int,double>* geneclusterID=NULL;
				if(strstr(speciesList[s].c_str(),"Anc")!=NULL)
				{
					if(geneAssign.find(speciesList[s])!=geneAssign.end())
					{
						geneclusterID=geneAssign[speciesList[s]];
					}
								
				}
				else
				{
					if(gset->find(speciesList[s])==gset->end())
					{
						geneclusterID=NULL;
					}
					else 
					{	
						string& geneName=(*gset)[speciesList[s]];
						if(geneAssign.find(geneName)!=geneAssign.end())
						{
							geneclusterID=geneAssign[geneName];
						}
					}
				}
				if(geneclusterID!=NULL)
				{
					for(map<int,double>::iterator aIter=geneclusterID->begin();aIter!=geneclusterID->end();aIter++)
					{
						oFile <<"\t"<< aIter->second;
					}
				}
				else
				{
					for(int k=0;k<maxClusterCnt;k++)
					{
						oFile <<"\t-1";
					}
				}
			}
			oFile << endl;
		}
	}
	return 0;
}


map<string,int>* 
SpeciesClusterManager::getGenesForSpecies(string& specName)
{
	if(speciesClusterSet_Genewise.find(specName)==speciesClusterSet_Genewise.end())
	{	
		return NULL;
	}
	return speciesClusterSet_Genewise[specName];
}

int
SpeciesClusterManager::generateData(const char* outputDir, string& lcaName,vector<string>& speciesList)
{
	map<int,map<string,int>*> clusterAssignments;
	gammaMgr->sampleAllClusterAssignments(clusterAssignments);
	char fName[1024];
	char dirName[1024];
	map<string,ofstream*> filePtrSet;
	map<string,ofstream*> filePtrClusteredSet;
	map<string,ofstream*> clusterFilePtrSet;
	for(map<string,CLUSTERSET*>::iterator sIter=speciesExpertSet.begin();sIter!=speciesExpertSet.end();sIter++)
	{
		sprintf(fName,"%s/%s_samples.txt",outputDir,sIter->first.c_str());
		ofstream* oFile=new ofstream(fName);
		filePtrSet[sIter->first]=oFile;
		sprintf(fName,"%s/%s_clusterassign.txt",outputDir,sIter->first.c_str());
		ofstream* cFile=new ofstream(fName);
		clusterFilePtrSet[sIter->first]=cFile;
	}
	vector<double> sampleValues;
	char clusterAssignmentFName[1024];
	char allassignfname[1024];
	sprintf(allassignfname,"%s/allspecies_clusterassign.txt",outputDir);
	ofstream allAssignFile(allassignfname);
	map<string,int> shownGenes;
	for(map<int,map<string,int>*>::iterator gIter=clusterAssignments.begin();gIter!=clusterAssignments.end();gIter++)
	{
		
		int ogid=gIter->first;
		map<string,int>* assign=gIter->second;
		map<int,map<string,int>*> geneAssignSet;
		map<int,map<string,string>*> geneNamesSet;
		for(map<string,int>::iterator cIter=assign->begin();cIter!=assign->end();cIter++)
		{
			char tempBuff[1024];
			strcpy(tempBuff,cIter->first.c_str());
			char* tok=strtok(tempBuff,":");
			int tokCnt=0;
			string speciesName;
			string aName;
			int dup=-1;
			while(tok!=NULL)
			{
				if(tokCnt==0)
				{
					dup=atoi(tok);
				}
				else if(tokCnt==1)
				{
					speciesName.append(tok);
				}	
				else
				{
					aName.append(tok);
				}
				tok=strtok(NULL,":");
				tokCnt++;
			}
			map<string,int>* geneSet=NULL;
			map<string,string>* geneName=NULL;
			if(geneAssignSet.find(dup)==geneAssignSet.end())
			{
				geneSet=new map<string,int>;
				geneName=new map<string,string>;
				geneAssignSet[dup]=geneSet;
				geneNamesSet[dup]=geneName;	
			}
			else
			{
				geneSet=geneAssignSet[dup];
				geneName=geneNamesSet[dup];
			}
			(*geneSet)[speciesName]=cIter->second;
			(*geneName)[speciesName]=aName;
		}
		for(map<int,map<string,int>*>::iterator gIter=geneAssignSet.begin();gIter!=geneAssignSet.end();gIter++)
		{
			map<string,int>* assignment=gIter->second;
			map<string,string>* geneName=geneNamesSet[gIter->first];
			for(map<string,int>::iterator sIter=assignment->begin();sIter!=assignment->end();sIter++)
			{
				if(clusterFilePtrSet.find(sIter->first)==clusterFilePtrSet.end())
				{
					sprintf(fName,"%s/%s_clusterassign.txt",outputDir,sIter->first.c_str());
					ofstream* cFile=new ofstream(fName);
					clusterFilePtrSet[sIter->first]=cFile;
				}
				ofstream* cFile=clusterFilePtrSet[sIter->first];
				string& aName=(*geneName)[sIter->first];
				if(strstr(sIter->first.c_str(),"Anc")!=NULL)
				{
					aName=(*geneName)[srcSpecies];
				}
				(*cFile) <<aName <<"\t" << sIter->second<<endl;
				if(speciesExpertSet.find(sIter->first)==speciesExpertSet.end())
				{
					continue;	
				}
				CLUSTERSET* clusterSet=speciesExpertSet[sIter->first];
				Expert* e=(*clusterSet)[sIter->second];
				e->generateSample(r,sampleValues);
				ofstream* oFile=filePtrSet[sIter->first];
				(*oFile) << aName;
				for(int j=0;j<sampleValues.size();j++)
				{
					(*oFile) <<"\t" << sampleValues[j];
				}
				(*oFile) << endl;
				sampleValues.clear();
			}
		}
		int lcacnt=0;
		map<string,int> localKey;
		const char* dupAnc=gammaMgr->getDupAncestor(ogid);
		bool wgd=false;
		if(dupAnc!=NULL && strcmp(dupAnc,"Anc5")==0)
		{
			wgd=true;
		}
		for(map<int,map<string,int>*>::iterator hIter=geneAssignSet.begin();hIter!=geneAssignSet.end();hIter++)
		{
			map<string,string>* gname=geneNamesSet[hIter->first];
			map<string,int>* gset=hIter->second;
			int ancClusterAssign=-3;
			if(gset->find(lcaName)!=gset->end())
			{
				ancClusterAssign=(*gset)[lcaName];
				lcacnt++;
			}
			//Now use the scername to show this og
			string dispName;
			if(gset->find(srcSpecies)==gset->end())
			{
				char tempName[1024];	
				sprintf(tempName,"OG%d_%d",ogid,hIter->first);
				dispName.append(tempName);
				if(wgd)
				{
					dispName.append("*");
				}
			}	
			else
			{
				string geneName((*gname)[srcSpecies]);
				shownGenes[geneName]=0;
				dispName.append(gnm.getCommonName((*gname)[srcSpecies].c_str()));
				if(localKey.find(dispName)!=localKey.end())
				{
					dispName.append("_2");
				}
				localKey[dispName]=0;
				if(wgd)
				{
					dispName.append("*");
				}
			}
			allAssignFile <<dispName;
			string members;
			for(int s=0;s<speciesList.size();s++)
			{
				int geneclusterID=-2;
				if(gset->find(speciesList[s])==gset->end())
				{
					geneclusterID=-2;
				}
				else 
				{	
					geneclusterID=(*gset)[speciesList[s]];
				}
				if(members.length()>0)
				{
					members.append(";");
				}	
				if(gname->find(speciesList[s])==gname->end())
				{
					members.append("NULL");
				}
				else
				{
					members.append((*gname)[speciesList[s]]);
				}		
				allAssignFile<< "\t"<< geneclusterID;
			}
			allAssignFile << endl;	
		}
	}
	for(map<string,ofstream*>::iterator fIter=clusterFilePtrSet.begin();fIter!=clusterFilePtrSet.end();fIter++)
	{
		fIter->second->close();
		if(filePtrSet.find(fIter->first)==filePtrSet.end())
		{
			continue;
		}
		ofstream* file=filePtrSet[fIter->first];
		file->close();
	}
	
	return 0;
}

//set the maximum number of iterations to be done if the expectation/maximization step for the learning function doesn't meet the convergence condition.
int
SpeciesClusterManager::setNMaxIterations(int n)//set the maximum number of iteratnction doesn't meet the convergence condition.
{
        //SK: set the class object equal to the input variable n
        nMaxIterations=n;
        //SK:exit the function
        return 0;
}

//SK: function to set the convergence threshold condition
int
SpeciesClusterManager::setConversionThreshold(double t)
{
        //SK: set the class variable to the argument variable t
        convThresh=t;
        return 0;
}

//SK: function to set training and test ogids
int
SpeciesClusterManager::setTrainingAndTestSet(map <int,int>& foldids,int f)
{
        //SK: loop over elements in foldids
        for(map<int,int>::iterator itr=foldids.begin();itr!=foldids.end();itr++)
        {
                //SK: set orthogroup variable
                int ogid=itr->first;
                //SK: set fold id variable
                int fid=itr->second;
                //SK: check if the id is for the identified fold
                if(fid==f)
                {
                        //SK: enter ogid into test set if it is
                        testSet[ogid]=0;
                }
                else
                {
                        //SK: otherwise place in training set
                        trainingSet[ogid]=0;
                }
        }
        return 0;
}

//SK: function to set the boolean CVMode class variable
//SK: this varible is initialised to false in the constructor object.
int
SpeciesClusterManager::setCVMode(bool in)
{
        //SK: set the variable
        CVMode=in;
        //SK: exit the function
        return 0;
}

//SK: function for writing out test score information
int
SpeciesClusterManager::writeScores(const char* outputDir)
{
        char outputFName[1024];
        sprintf(outputFName,"%s/likelihood.txt",outputDir);
        ofstream oFile(outputFName);
        oFile << scores[0] << "\t" << scores[1] << "\t" << scores[2] << "\t" << scores[3] << endl;
        oFile.close();
        if(CVMode==true && testScores.size()>0)
        {
                char outputFNameT[1024];
                sprintf(outputFNameT,"%s/likelihoodTest.txt",outputDir);
                ofstream oFileT(outputFNameT);
                oFileT << testScores[0] << "\t" << testScores[1] << "\t" << testScores[2] << "\t" << testScores[3] << endl;
                oFileT.close();
        }
        return 0;
}

//SK: function to set the second stage option
int
SpeciesClusterManager::setSecondStageOption(bool in)
{
        secondStage=in;
        return 0;
}
