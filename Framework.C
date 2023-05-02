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
#include <vector>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstring>
#include <getopt.h>
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

//SK: for fgconverter component
#include "Gene.H"
#include "GeneManager.H"
#include "Protein.H"
#include "ProteinManager.H"
#include "Interaction.H"
#include "InteractionManager.H"
#include "Node.H"
#include "Graph.H"
#include "BioNetwork.H"

//SK: for learnMoE component
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "MotifManager.H"
#include "BFGSWrapperData.H"
#include "BFGSWrapper.H"
#include "TestData.H"
#include "MotifRegressor.H"

bool GENFROMFILE=false;//SK: define two global variables for the gene tree mode. First a boolean vaiable if the gene trees should be generated from the gene tree files, or the mapped ortho group reader object.
char* GENETREEDIR; //SK: Second, the name of the directory that the .tre files are in.

Framework::Framework()
{
	convThresh=0.5;
	secondStage=true;
	vType=CONT;
	epsThreshold=-1;
        untransformedData[0]='\0';
	predictionInputMode=false;
	sourceInit=false;
}

Framework::~Framework()
{
}

int
Framework::setSourceInitOption(bool in)
{
	sourceInit=in;
	return 0;
}

int 
Framework::readSpeciesData(const char* aFName, const char* rand,const char* Dir)
{
	if(predictionInputMode==true)
        {
                scMgr.setPredictionInputMode(true);
	}
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
	scMgr.setNMaxIterations(nMaxIterations); 
        scMgr.setConversionThreshold(convThresh);
	scMgr.setMaxClusterCnt(maxClusterCnt);
	if(!preClustering)
	{
		scMgr.readSpeciesData(aFName);
	}
	else
	{
		scMgr.readSpeciesData(aFName,Dir);
	}
	randnum=gsl_rng_alloc(gsl_rng_default);
	return 0;
}

int 
Framework::readSpeciesTree(int clusterCnt, const char* aFName)
{
	maxClusterCnt=clusterCnt;
	sdMgr.setMaxClusters(clusterCnt);
	sdMgr.readSpeciesTree(aFName);
	vector<string> speciesList;
	sdMgr.getSpeciesListPrefix(speciesList);
	for(int i=0;i<speciesList.size();i++)
	{	
		cout << speciesList[i] << endl;
	}
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
//Framework::startClustering(const char* aDir)
Framework::startClustering(const char* aDir, const char* treeFName)	// repop
{
	//SK: set maximum number of iterations to allow
        scMgr.setNMaxIterations(nMaxIterations);
        //SK: set convesion threshold
        scMgr.setConversionThreshold(convThresh);
	//SK: set the second stage parameter
        scMgr.setSecondStageOption(secondStage);
	strcpy(outputDir,aDir);
	scMgr.initExperts(sourceInit);
	cout <<"Total updated parent nodes "<< gammaMgr.getTotalUpdatedParentCnt() << endl;
        gammaMgr.showTotalUpdatedParents();
	if(predictionInputMode==false)
	{
		initClusterTransitionProb();
		scMgr.setMergedOGIDSet(mergedOgidSet);
		vector<string> speciesList;
		sdMgr.getSpeciesListPrefix(speciesList);
		//SK: initialize cluster assignments for genes in orthogroups that can be added from the non-source species merged data in prediction mode
		//scMgr.setPredictionInputMode(true);
		//scMgr.executePredictionMode(outputDir,speciesList,treeFName);	// repop
		//scMgr.setPredictionInputMode(false);
		//SK:at this point the cluster assingments should be initialized for all genes in orthogrous that can be used.
		//SK:these clusters assingments will be in the prediction subdirectory of the output directory.  
		scMgr.estimateExpertParameters(outputDir);
		double newScore=scMgr.getScore();
		//SK: set the second stage parameter
		scMgr.setSecondStageOption(secondStage);
		scMgr.dumpAllInferredClusterAssignments(outputDir);
		double newScore_PP=scMgr.getScore();
		cout <<"Score before PP " << newScore << "\t" << " Score after PP " << newScore_PP << endl;
		scMgr.showClusters_Extant(outputDir);
		scMgr.showClusters_Ancestral(outputDir);
		scMgr.showMeans(outputDir);
		//SK: write out scores
		scMgr.writeScores(outputDir);
		//This is only for visualization purposes
		for(int i=0;i<speciesList.size();i++)
		{	
			cout << speciesList[i] << endl;
		}
		scMgr.dumpAllInferredClusters_ScerwiseGrouped(outputDir,speciesList);
		//scMgr.dumpAllInferredClusters_LCA(outputDir,speciesList,sdMgr.getRoot()->name);
		scMgr.dumpAllInferredClusters_LCA(outputDir,speciesList,sdMgr.getRoot()->name,treeFName);	// repop
		scMgr.showClusters_Ancestral(outputDir);
		scMgr.dumpAllInferredClusterGammas(outputDir,speciesList);
		sdMgr.showInferredConditionals(outputDir);
		sdMgr.showInferredConditionals_ML(outputDir);
	}
	else if(predictionInputMode==true)
        {       
                cout << "Initializing prediction mode" << endl;
                scMgr.setPredictionInputMode(true);
                scMgr.setMergedOGIDSet(mergedOgidSet);
		setTransitionMatrices();
		vector<string> speciesList;
                sdMgr.getSpeciesListPrefix(speciesList);
                for(int i=0;i<speciesList.size();i++)
                {
                        cout << speciesList[i] << endl;
                }
		scMgr.executePredictionMode(outputDir,speciesList,treeFName);	// repop
        }
	return 0;
}

int
Framework::readTransitionMatrices(const char* cFile)
{
	//SK: first read in fourth column of configuration file for prediction mode. 
	ifstream inFile(cFile);
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
                string transitionFName;
                while(tok!=NULL)
                {
                        if(tokCnt==0)
                        {
                                speciesName.append(tok);
                        }
                        else if(tokCnt==3)
                        {
                                transitionFName.append(tok);
                        }
                        tok=strtok(NULL,"\t");
                        tokCnt++;
                }
		//SK: for species speciesName, read in the transitionMatrix
		Matrix* toAdd = new Matrix(maxClusterCnt,maxClusterCnt);
		cout << speciesName << "\t" << transitionFName << endl;
		ifstream inFile2(transitionFName.c_str());
		int r=0;
		char buffer2[1024];
		while(inFile2.good())
		{
			inFile2.getline(buffer2,1023);
			if(strlen(buffer2)<=0)
			{
				continue;
			}
			char* tok2=strtok(buffer2," ");
			int c=0;
			while(tok2!=NULL)
			{
				double v=atof(tok2);
				toAdd->setValue(v,r,c);
				tok2=strtok(NULL," ");
				c++;
			}
			r++;
		}
		cout << "Input transition matrix for " << speciesName << endl;
		toAdd->showMatrix();
		transitionMatrixMap[speciesName]=toAdd;
	}	
	return 0;
}

int
Framework::setTransitionMatrices()
{
        SpeciesDistManager::Species* root=sdMgr.getRoot();
	Matrix* forRoot = new Matrix(1,maxClusterCnt);
	for(int i=0;i<maxClusterCnt;i++)
	{
		double v=transitionMatrixMap.find(root->name)->second->getValue(0,i);
		forRoot->setValue(v,0,i);
	}
	root->setParams(forRoot);
       	setTransitionMatrices(root->leftchild);
	setTransitionMatrices(root->rightchild);
        return 0;
}

int
Framework::setTransitionMatrices(SpeciesDistManager::Species* anode)
{
	if(transitionMatrixMap.find(anode->name)!=transitionMatrixMap.end())
	{
		anode->setParams(transitionMatrixMap.find(anode->name)->second);
	}
        if(anode->leftchild!=NULL)
        {
        	setTransitionMatrices(anode->leftchild);
        }
        if(anode->rightchild!=NULL)
        {      
                setTransitionMatrices(anode->rightchild);
        }
        return 0;
}

//SK: function for cross validation proceedure
int 
//Framework::startCrossValidationCheck(const char* aDir,int cvF,const char* rand,const char* confFName,const char* specOrder, const char* orthomapfile)
Framework::startCrossValidationCheck(const char* aDir,int cvF,const char* rand,const char* confFName,const char* specOrder, const char* orthomapfile, const char* treeFName)	// repop
{
	cout << "Starting cross validation check" << endl;//SK: State that the code is starting cross validation
	r_og=gsl_rng_alloc(gsl_rng_default);//SK: define the random number generator
        int rseed=getpid();//SK: set set vairable
        gsl_rng_set(r_og,rseed);//SK: initialize the randome number generator
        cout << "Random seed for partitioning orthogroups is "<< rseed << endl;//SK: print out seed setting
	//SK: partition the OGID groups
	map <int,int> FoldIDs;//SK: define map of of fold id's for the orthogroups
	//SK: get the map of mapperorthogroup objects from the class mor object
	map<int,MappedOrthogroup*>& ogSet=mor.getMappedOrthogroups();
	//SK: loop over each orthogroups
	for(map<int,MappedOrthogroup*>::iterator oIter=ogSet.begin();oIter!=ogSet.end();oIter++)
	{
		//SK: set ogid number
		int og=oIter->first;
		//SK: get random number from 0 to 1 from the random number generator
		double pval=gsl_ran_flat(r_og,0,1);
		//SK: the inverse of the number of cross validation fold we're dividing the data into.
                double step=1.0/(double)cvF;
		//SK: obtain an integer number from 0 to cvF-1
                int foldid=(int)(floor(pval/step));
		//SK: set the fold id the FildIDs map object.
		FoldIDs[og]=foldid;
		//cout << foldid << endl;
	}
	cout << "Partitioned orthogroups" << endl;
	//now apply the clustering in the cross validation context
	for(int f=0;f<cvF;f++)
	{
		if(cvF==1)
		{
			f=1;
		}
		//SK: set up output directory for training data results
		//SK: define the varible for the  output directory to be used for writing out the results from analyzing the training data. 
		char outDirF[1024]; 
		//SK: set the output directory name
		sprintf(outDirF,"%s/fold%i",aDir,f);
		//SK: set the varible for the command to make the output directory 
		char outputDirCmdF[1024];
		//SK: fill the string for the command to make the directory
                sprintf(outputDirCmdF,"mkdir -p  %s",outDirF);
		//SK: call the system to make the output directory.
                system(outputDirCmdF);
		//SK: print out the output directory name 
		cout << outDirF << " made for fold " << f << " training" << endl;
		
		//SK: set up application of Arboretum to training data, here the set up the SpeciesClusterManager object and the GammaManager object(s) is the main task. 
		
		//SK: define the species cluster manager object to do the cross validation                 
		SpeciesClusterManager scMgrTrn;
		//SK: print  out the fold id for which the cross validation is being done
		cout << "Generated fold " << f << " training OGIDS." << endl;
		//SK: set the mor object in scMgrTrn
        	scMgrTrn.setOrthogroupReader(&mor);
		//SK: Define the Gamma manager object for the training OGIDS set
		GammaManager gammaMgrTrn;
		//SK: set mor in gammaMgrTrn
        	gammaMgrTrn.setOrthogroupReader(&mor);
		//SK: set the number of clusters in gammaMgrTrn
        	gammaMgrTrn.setMaxClusterCnt(maxClusterCnt);
		//SK: set the species distance manager  object
        	gammaMgrTrn.setSpeciesDistManager(&sdMgr);
		//SK: initialize a spearate GammaManager object for the test set data in the SpeciesClusterManager object
                GammaManager gammaMgrTst;
                gammaMgrTst.setOrthogroupReader(&mor);
                gammaMgrTst.setMaxClusterCnt(maxClusterCnt);
                gammaMgrTst.setSpeciesDistManager(&sdMgr);
		//SK: the following lines emulate those in Framework::readSpeciesData above, but for scMgrTrn		

		//SK: first set of randomization variable information
		if(strcmp(rand,"none")==0)
                {
                        //SK set randomization option to false
                        scMgrTrn.setRandom(false);
                }
                //SK: check if randomization option is set to true
                else if(strcmp(rand,"yes")==0)
                {
                        //SK: set randomization option
                        scMgrTrn.setRandom(true);
                }
                //SK: check if randomization option is set to true
                else if(isdigit(rand[0]))
                {
                        //SK: set randomization to true in scMgrTrn
                        scMgrTrn.setRandom(true);
                        //SK: set randomization seed number
                        scMgrTrn.setRandSeed(atoi(rand));
                }
		//SK: set the CVMode object to true
                scMgrTrn.setCVMode(true);
                //SK: set training and test set OGIDS
                scMgrTrn.setTrainingAndTestSet(FoldIDs,f);
                //SK: set gamma manager for the spcMgrTrn
                scMgrTrn.setGammaManager(&gammaMgrTrn);
                //SK: set the GammaManager for the test OGIDS, since we're in CV mode
                scMgrTrn.setGammaManager_Test(&gammaMgrTst);
		//SK: set the CVMode object to true
                scMgrTrn.setCVMode(true);
                //SK: set secondStage option.
                scMgrTrn.setSecondStageOption(secondStage);
                //SK: set training and test set OGIDS
                scMgrTrn.setTrainingAndTestSet(FoldIDs,f);
                //SK: set gamma manager for the spcMgrTrn
                scMgrTrn.setGammaManager(&gammaMgrTrn);
                //SK: set the GammaManager for the test OGIDS, since we're in CV mode
                scMgrTrn.setGammaManager_Test(&gammaMgrTst);
		//SK: set the source species in scMgrTrn
                scMgrTrn.setSrcSpecies(srcSpecies);
                //SK: set maximum number of iterations to allow
                scMgrTrn.setNMaxIterations(nMaxIterations);
                //SK: set convesion threshold
                scMgrTrn.setConversionThreshold(convThresh);
                //SK: set the number of clusters in scMgrTrn
                scMgrTrn.setMaxClusterCnt(maxClusterCnt);
                //SK: read in the species data to scMgeTrn
                if(!preClustering)
		{
			scMgrTrn.readSpeciesData(confFName);
		}
		else
		{
			scMgrTrn.readSpeciesData(confFName,aDir);
		}	
		//SK: set randnum variable, necessary for initClusterTransitionProb() to complete
                randnum=gsl_rng_alloc(gsl_rng_default);

		//SK: the following follows what is done in Framework::startClustering(const char* aDir)

		//SK: initialize expert parameters
                scMgrTrn.initExperts(sourceInit);
		//SK: needs to be done from Framework::startClustering(const char* aDir)
		gammaMgrTrn.showTotalUpdatedParents();
		//SK: initialize the SpeciesDistManager
                initClusterTransitionProb();
		//apply Aboretum to trainging data
        	scMgrTrn.estimateExpertParameters(outDirF);
		//SK: get score for training data analysis
        	double newScore=scMgrTrn.getScore();
		//SK: writeout cluster assignments from training data analysis
        	scMgrTrn.dumpAllInferredClusterAssignments(outDirF);
		//SK: get corrected score for the training data
        	double newScore_PP=scMgrTrn.getScore();//SK: get corrected score for the training data
		//write out scores for the training data analysis
        	//cout <<"Score before PP " << newScore << "\t" << " Score after PP " << newScore_PP << endl;
		//SK: write out training data results as exemplified in the Framework::startClustering function above
		//SK: show and write out means from model
        	//scMgrTrn.showMeans(outDirF);
        	//SK: This is only for purposes of writting out the results
        	//vector<string> speciesList;//SK: define vector species list
        	//sdMgr.getSpeciesListPrefix(speciesList);//SK: fill species names into the speciesList vector
        	//for(int i=0;i<speciesList.size();i++)//SK: loop over all species names
		//{
                //	cout << speciesList[i] << endl;//SK: print out species names
        	//}
		//SK: write out the clusterassignment matrix files
        	//scMgrTrn.dumpAllInferredClusters_ScerwiseGrouped(outDirF,speciesList);
		//SK: write out cluster assignment matrix files
        	//scMgrTrn.dumpAllInferredClusters_LCA(outDirF,speciesList,sdMgr.getRoot()->name);
		//SK: write out gamma matrices for each species
        	//scMgrTrn.dumpAllInferredClusterGammas(outDirF,speciesList);
		//SK: write out conditionals from the species dist manager
        	//sdMgr.showInferredConditionals(outDirF);
        	//sdMgr.showInferredConditionals_ML(outDirF);
		//SK:  cluster assignments for the extant species
                //scMgrTrn.showClusters_Extant(outDirF);
                //SK: cluster assigments for ancestral species
                //scMgrTrn.showClusters_Ancestral(outDirF);
		cout <<"Score before PP " << newScore << "\t" << " Score after PP " << newScore_PP << endl;
		scMgrTrn.showClusters_Extant(outDirF);
		scMgrTrn.showClusters_Ancestral(outDirF);
		scMgrTrn.showMeans(outDirF);
		//SK: write out scores
		scMgrTrn.writeScores(outDirF);
		//This is only for visualization purposes
		vector<string> speciesList;
		sdMgr.getSpeciesListPrefix(speciesList);
		for(int i=0;i<speciesList.size();i++)
		{
			cout << speciesList[i] << endl;
		}
		scMgrTrn.dumpAllInferredClusters_ScerwiseGrouped(outDirF,speciesList);
		//scMgrTrn.dumpAllInferredClusters_LCA(outDirF,speciesList,sdMgr.getRoot()->name);
		scMgrTrn.dumpAllInferredClusters_LCA(outDirF,speciesList,sdMgr.getRoot()->name,treeFName);	// repop
		scMgrTrn.showClusters_Ancestral(outDirF);
		scMgrTrn.dumpAllInferredClusterGammas(outDirF,speciesList);
		sdMgr.showInferredConditionals(outDirF);
		sdMgr.showInferredConditionals_ML(outDirF);
		//SK: finally apply the get score function for the test data in scMgrTrn, which is unique to this cross validation part of the code. 
		double testDataScore = scMgrTrn.getScore_Test();
		//SK: write out scores.
		scMgrTrn.writeScores(outDirF);
	}
}

int
Framework::generateData(const char* outputDir)
{
	scMgr.initExperts(sourceInit);
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
	scMgr.initExperts(sourceInit);
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
Framework::setNMaxIterations(int n)
{
	nMaxIterations=n;
	return 0;
}

int 
Framework::setConversionThreshold(double t)
{
	convThresh=t;
	return 0;
}

//SK: function to set the option for running the second optimization step in SpeciesClusterManager::dumpAllInferredClusterAssignments(const char* outputDir)
int
Framework::setSecondStageOption(bool in)
{
	secondStage=in;
	return 0;
}

void print_usage() 
{
	//SK: if the number of arguments is incorrect so now print out the following copyright and usage information.
	//SK: print copyright information
	cout <<"Arboretum: Copyright (C) 2013 Sushmita Roy" << endl;
	cout <<"This program comes with ABSOLUTELY NO WARRANTY; for details type ./arboretum. "
		<<"This is free software, and you are welcome to redistribute it"
		<<" under certain conditions;" <<endl << endl;
	//SK: print example usage infromation
	cout <<"Usage: ./arboretum -s specorder -e orthogroup -k maxk -t speciestree -c clusterassignments -r rand -o outputDir -m learn -b baseSpecies -i inittype -p p_diagonal_nonleaf" << endl;
	cout << "Required options:" << endl;
	//SK: print information on the required variables that are needed as input for running a command. 
	cout << "-s\t\tFile listing species in the orthology." << endl;
	cout << "-e\t\tFile listing the orthology relationships of genes across species." << endl;
	cout << "-k\t\tNumber of clusters in each species." <<  endl;
	cout << "-t\t\tSpecies tree file of species which are represented by the input data in the analysis" << endl;
	cout << "-c\t\tThe input cluster assignment and expression data for each species in the analysis." << endl;
	cout << "This file of the form: species <tab> cluster_assign_file <tab> expression_data_file <endl>" << endl;
	cout << "-r\t\tOption to randomize input cluster assignments rseed|none|yes." << endl;
	cout << "-o\t\tOutput directory." << endl;
	cout << "-m\t\tDefines the mode in which the algorithm is to be used; learn|generate|visualize|crossvalidation|prediction are the options." << endl;
	cout << "When in prediction mode the format of the configuration file (pointed to with the -c option) should be : species <tab> cluster_assign_file <tab> expression_data_file <tab> transition_matrix_file <end>" << endl;
	cout << "-b\t\tA well annotated species to which gene names of other species will be mapped in the *_clusterassign.txt output." << endl;
	cout << "-i\t\tInitialization method for transition probabilities for cluster membership across species, uniform|branchlength." << endl;
	cout << "-p\t\tInitial transition probability values, either a single value or a file defining values every specie tree branch." << endl;
	//SK: print comment on the initialization type and probability value type
	cout <<"The init--type setting is for specifying how the diagonal cluster membership transition probabilities across species will be initialized. If inittype is uniform then the -p argument should be the default initial value for diagonal transition probabilities on all branches of the species tree. If the branchlength option is used, then non-uniform transition probabilities will be taken from a file set by the -p option." << endl;
	//SK: now also print out the information about the non-required variables. 
	cout << "\n\n\nThe following are non-required options." << endl;
	cout << "-g\t\tDirectory containing gene trees, which can be used to directly define the GeneTreeManager gene mappings." << endl;
	cout << "-d,--conv-thresh\t\tThis is the change in the likelihood score between iterations that defines the convergence of the algorithm." << endl;
	cout << "-n\t\tThe maximum number of iterations if convergence condition is not met." << endl;
	cout << "-v\t\tThis option turns on a cross validation test of the clustering; and defines the number of partitions to use." << endl;
	cout << "-w\t\tA true|false option to overwrite the output directory if it already exists, the default is false." << endl;
	cout << "-l\t\tA true|false option to run the last EM optimization step in SpeciesClusterManager::dumpAllInferredClusterAssignments(). Default is True." << endl;
	cout << "-f\t\tA true|false option to run the merged clustering of the data to generate input clusterassignments. Default is False." << endl;
	cout << " When this option is used the input configuration (-c) file should be in the form: species <tab> expression_data_file <endl>" << endl;
	cout << "-u\t\tA true|false option for updating the cluster means from the source species." << endl;
	cout << "-h,--help\t\tReports this usage information." << endl;
	cout << "A note on prediction mode: this is most commonly used when you have an arboretum clustering result to which you have applied a reordering, and hence need to reinfer ancestral assignments for the reordered clusters." << endl;
	return;
}

int
Framework::readExpressionData(const char* aFName,const char* bFName)
{
	ofstream oFile(bFName);
	map < string , map <string , vector <string> > > MAPS;
	vector <string> specOrder;
	map <int,int> ogidSet;
        map <string,int> geneCnts;
	map <string,int> dataCnts;
	int uniformGenes;
	int uniformGenesWithExpr;
	int genesWithExpr=0;
        ifstream inFile(aFName);
	int maxGeneCnt=2;
        char buffer[1024];
	speciesList.clear();
	cout << "Reading " << aFName << endl;
	//SK: write out header of output merged file.
	//SK: read in data and proceed 
        while(inFile.good())
        {
                inFile.getline(buffer,1023);
                if(strlen(buffer)<=0)
                {
                        continue;
                }
                char* tok=strtok(buffer,"\t");
                int tokCnt=0;
		if(strncmp(buffer,"Anc",3)==0)
                {
                        continue;
                }
                string specName;
                string fileName;
                while(tok!=NULL)
                {
                        if(tokCnt==0)
                        {
                                specName.append(tok);
                        }
                        else if(tokCnt==1)
                        {
                                fileName.append(tok);
                        }
                        tok=strtok(NULL,"\t");
                        tokCnt++;
                }
                cout << specName << endl;
                cout << fileName << endl;
		speciesList[specName]=0;
		specOrder.push_back(specName);
		map <string , vector <string> > gexp;
		ifstream inFileS(fileName.c_str());
        	char bufferS[1024];
        	while(inFileS.good())
        	{
                	inFileS.getline(bufferS,1023);
                	if(strlen(bufferS)<=0)
                	{
                        	continue;
                	}
                	char* tokS=strtok(bufferS,"\t");
                	int tokCntS=0;
                	string geneName;
                	vector <string> data;
                	while(tokS!=NULL)
                	{
                        	if(tokCntS==0)
                        	{
                                	geneName.append(tokS);
                        	}
				else if(tokCntS>0)
				{
					string tmp;
					tmp.append(tokS);
					data.push_back(tmp);
				}
				tokS=strtok(NULL,"\t");
                        	tokCntS++;
			}
			gexp[geneName]=data;
			dataCnts[specName]=data.size();
		}
                MAPS[specName]=gexp;
		//cout << specName << "\t" << gexp.size() << endl;
		//cout << gexp.begin()->second.size() << endl;
		inFileS.close();
        }
        inFile.close();
	//SK: write out header of output merged file.
        oFile << "Name";
        for(int s=0;s<specOrder.size();s++)
        {
                int exp=dataCnts.find(specOrder[s])->second;
                for(int i=0;i<exp;i++)
                {
                        oFile << "\t" << specOrder[s] << "_Exp" << i;
                }
        }
        oFile << endl;
	//now aggregate the data across epcies as in mergeExpr_Uniform
	string baseSpecies;
	cout << MAPS.size() << endl;
	baseSpecies.append(srcSpecies);
	for(map < string , vector <string> >::iterator gIter=MAPS.find(baseSpecies)->second.begin();gIter!=MAPS.find(baseSpecies)->second.end();gIter++)
	{
		//cout << gIter->first << "\t" << baseSpecies << endl;
		int ogid=mor.getMappedOrthogroupID(gIter->first.c_str(),baseSpecies.c_str());
                MappedOrthogroup* agrp=mor.getMappedOrthogroup(gIter->first.c_str(),baseSpecies.c_str());
		if(agrp==NULL)
                {
			cout << "No orthogroup" << endl;
                        continue;
                }
                if(agrp->getOrthoMembers().size()<maxGeneCnt)
                {
			cout << "Orthogroup " << ogid << " only has " << agrp->getOrthoMembers().size() << " genes, continuing." << endl; 
                        continue;
                }
                /*if((filteredOGSet.size()>0) && (filteredOGSet.find(ogid)==filteredOGSet.end()))
                {
                        continue;
                }*/
                GeneMap* selfhit=agrp->getSpeciesHits(baseSpecies.c_str());
		bool isUniform=true;
                bool hasExpression=true;
                int speciesWithGeneExpr=0;
                for(map < string , map <string , vector <string> > >::iterator sIter=MAPS.begin();sIter!=MAPS.end();sIter++)
                {
                        //cout << sIter->first << endl;
                        if(strcmp(sIter->first.c_str(),baseSpecies.c_str())==0)
                        {
                                continue;
                        }
                        GeneMap* specieshit=agrp->getSpeciesHits(sIter->first.c_str());
                        if(specieshit==NULL)
                        {
                                isUniform=false;
                                continue;
                        }
                        if(specieshit->getGeneSet().size()>maxGeneCnt)
                        {
                                isUniform=false;
                        }
                }
                if(isUniform)
                {
                        //cout << "Is uniform" << endl;
                        uniformGenes++;
                }
                for(map <string, map <string , vector <string> > >::iterator sIter=MAPS.begin();sIter!=MAPS.end();sIter++)
                {
                        if(strcmp(sIter->first.c_str(),baseSpecies.c_str())==0)
                        {
                               continue;
                        }
                        //cout << "Species: " << sIter->first << endl;
                        GeneMap* specieshit=agrp->getSpeciesHits(sIter->first.c_str());
                        if(specieshit==NULL)
                        {
                                hasExpression=false;
                        }
                        else
                        {
                                map<string,int>* orthogenes=mor.getOrtholog(baseSpecies.c_str(),gIter->first.c_str(),sIter->first.c_str());
                                if(orthogenes==NULL)
                                {
					cout << "orthogenes NULL" << endl;
                                        hasExpression=false;
                                }
				else
                                {
                                        for(map<string,int>::iterator oIter=orthogenes->begin();oIter!=orthogenes->end();oIter++)
                                        {
						//cout << oIter->first << endl;
						if(MAPS.find(sIter->first)->second.find(oIter->first)!=MAPS.find(sIter->first)->second.end())
						{
                                                	vector <string> exprVect=MAPS.find(sIter->first)->second.find(oIter->first)->second;
							if(exprVect.size()>=1)
							{
                                                        	//cout << oIter->first << " has expression data." << endl;
                                                        	hasExpression=true;
                                                        	speciesWithGeneExpr++;
							}
                                                }
                                        }
                                }
                        }
		}
                if(speciesWithGeneExpr==0)
                {
                        //cout << "No expression data found\n";
                        continue;
                }
                if(isUniform)
                {
                        //cout << "Uniform Gene with expression" << endl;
			uniformGenesWithExpr++;
		}
                genesWithExpr++;
                vector<string> scerExpr=gIter->second;
                if(geneCnts.find(baseSpecies)==geneCnts.end())
                {
                        geneCnts[baseSpecies]=1;
                }
                else
                {
                        geneCnts[baseSpecies]=geneCnts[baseSpecies]+1;
                }
		vector <string> entries;
		for(int i=0;i<scerExpr.size();i++)
                {
                	entries.push_back(scerExpr[i]);
                }
		int numEmpty=0;
                for(int s=0;s<specOrder.size();s++)
                {
                        if(strcmp(specOrder[s].c_str(),baseSpecies.c_str())==0)
                        {
                               continue;
                        }
                        if(MAPS.find(specOrder[s])==MAPS.end())
                        {
                                continue;
                        }
                        map < string , vector < string > > specExprMap=MAPS.find(specOrder[s])->second;
                        int scnt=dataCnts.find(specOrder[s])->second;
                        GeneMap* specieshit=agrp->getSpeciesHits(specOrder[s].c_str());
                        vector <string> speciesExpr;
                        map<string,int>* orthogenes=mor.getOrtholog(baseSpecies.c_str(),gIter->first.c_str(),specOrder[s].c_str());
                        if(orthogenes!=NULL)
                        {
                                for(map<string,int>::iterator oIter=orthogenes->begin();oIter!=orthogenes->end();oIter++)
                                {
					if(MAPS.find(specOrder[s])->second.find(oIter->first)!=MAPS.find(specOrder[s])->second.end())
					{
						speciesExpr=MAPS.find(specOrder[s])->second.find(oIter->first)->second;
					}
                                }
                        }
                        if(speciesExpr.size()>=1)
                        {
				//SK: add values to vector if data is found
                                for(int i=0;i<speciesExpr.size();i++)
                                {
                                	entries.push_back(speciesExpr[i]);
                                }
                                if(geneCnts.find(specOrder[s])==geneCnts.end())
				{
				        geneCnts[specOrder[s]]=1;
                                }
                                if(geneCnts.find(specOrder[s])==geneCnts.end())
                                {
                                        geneCnts[specOrder[s]]=1;
                                }
                                else
                                {
                                        geneCnts[specOrder[s]]=geneCnts[specOrder[s]]+1;
                                }
                        }
                        else
                        {
				//SK: otherwise insert empty values
				string in;
                                in.append("<nodata>");
                                for(int i=0;i<dataCnts.find(specOrder[s])->second;i++)
                                {
					entries.push_back(in);
					numEmpty++;		
                                }
                        }
                }
		//SK: define here the fraction of missing values allowed. 
		if(numEmpty/entries.size()>2/specOrder.size())
		{
			//SK: continue to the next base species gene if that fraction is higher than the allowed threshold
			continue;
		}
		double total=0;
		string in;
                in.append("<nodata>");
		for(int i=0;i<entries.size();i++)
		{
			if(entries[i]!=in)
			{
				total += atof(entries[i].c_str());
			}
		}
		int dat=entries.size()-numEmpty;
		double ave=(double)total/(double)dat;
		//SK: now fill a double vector object with the average of present values for missing values, selecting  
		vector <double> values;
		for(int i=0;i<entries.size();i++)
                {
		
                        if(entries[i]!="<nodata>")
                        {
                                values.push_back(atof(entries[i].c_str()));
                        }
			else
			{
				values.push_back(ave);
			}
                }
                //SK: write out first line of the matrix file.
		oFile << gIter->first;
		for(int i=0;i<values.size();i++)
                {
			oFile << "\t" << values[i];
		}
		oFile << endl;
		expCntPerGene=values.size();
	}
        return 0;
}

int
Framework::readExpressionDataNonSrc(const char* aFName,const char* bFName)
{
        ofstream oFile(bFName);
        map < string , map <string , vector <string> > > MAPS; //object containing species and gene names and data in string format
	//SK: order of species
        vector <string> specOrder;
	//SK: set of OGIDS relevant to analysis
        map <int,int> ogidSet;
	//SK: set of numbers of genes relevant to the analysis for orthogroups
        map <string,int> geneCnts;
	//SK: set of numbers of measurements per gene for each species expression data set 
        map <string,int> dataCnts;
	//SK: interger to count if a given gene/orthogroup is uniform
        int uniformGenes;
        //SK interger to count the number of genes 
	int uniformGenesWithExpr; //integer to count if that orthogorups also has data for sufficient number of genes
	//SK: integer to count the number of genes with expression
        int genesWithExpr=0;
	//SK: input expression data configuration file name 
        ifstream inFile(aFName); 
        int maxGeneCnt=2;
        char buffer[1024];
        speciesList.clear();
        cout << "Reading " << aFName << endl;
	//SK: read in configuration file name 
        while(inFile.good())
        {
                inFile.getline(buffer,1023);
                if(strlen(buffer)<=0)
                {
                        continue;
                }
		if(strncmp(buffer,"Anc",3)==0)
                {
                        continue;
                }
                char* tok=strtok(buffer,"\t");
                int tokCnt=0;
                string specName;
                string fileName;
		//SK: for each line of the configuration file read in the species name and the data file name.
                while(tok!=NULL)
                {
                        if(tokCnt==0)
                        {
                                specName.append(tok);
                        }
                        else if(tokCnt==1 && !predictionInputMode)
                        {
                                fileName.append(tok);
                        }
			else if(tokCnt==2 && predictionInputMode)
                        {
                                fileName.append(tok);
                        }
                        tok=strtok(NULL,"\t");
                        tokCnt++;
                }
		//print out that species name and data file name, and read in that data file to teh MAPS object.
                cout << specName << endl;
                cout << fileName << endl;
		//SK: add species name to list of species of interest
                speciesList[specName]=0;
                specOrder.push_back(specName);
		//SK : define map object for reading in the expression data information 
                map <string , vector <string> > gexp;
                ifstream inFileS(fileName.c_str());
                char bufferS[1024];
		//SK: read in expression data file for this species and save the expression data in string format for each gene in this species data sets
                while(inFileS.good())
                {
                        inFileS.getline(bufferS,1023);
                        if(strlen(bufferS)<=0)
                        {
                                continue;
                        }
                        char* tokS=strtok(bufferS,"\t");
                        int tokCntS=0;
			//SK: define variables for gene name and expression data for this entry in the species data matrix
			string geneName;
                        vector <string> data;
                        while(tokS!=NULL)
                        {
                                if(tokCntS==0)
                                {
					//SK: set gene name 
                                        geneName.append(tokS); //read in gene name 
                                }
                                else if(tokCntS>0)
                                {
					//SK: fill vector object for data
                                        string tmp;
                                        tmp.append(tokS);
                                        data.push_back(tmp); //read in data as string.
                                }
                                tokS=strtok(NULL,"\t");
                                tokCntS++;
                        }
			//SK: fill the data map object for this species
                        gexp[geneName]=data; 
			//SK: fill the dataCnts object with the numbre of measurements for each gene in thsi species.
                        dataCnts[specName]=data.size();
                }
		//SK: ender the data object for this species in the MAPS map object for the data for all  species
                MAPS[specName]=gexp; 
                //cout << specName << "\t" << gexp.size() << endl; 
                //cout << gexp.begin()->second.size() << endl;
                inFileS.close();
        }
        inFile.close(); 
	//SK: at this point the data matrices are read in as string objects, organized by species name and gene name in the MAPS object. 
	//SK: now write header line of the output file 
	oFile << "Name";//write header line of the likes. 
        for(int s=0;s<specOrder.size();s++)
        {
                int exp=dataCnts.find(specOrder[s])->second;
                for(int i=0;i<exp;i++)
                {
                        oFile << "\t" << specOrder[s] << "_Exp" << i;
                }
        }
	oFile << endl; 
	//SK: now being to write out the merged data by orthogroup and duplication level 
	//SK: define the number of species of interest
       	int NSpc=specOrder.size();
        int minSpcCnt=2;
	//SK: obtain orthogroup set
	map<int,MappedOrthogroup*>& ogSet=mor.getMappedOrthogroups();
	//SK: get set of species names and ids for the orthogy information in the mor object
	map<int,string> allSpeciesIDs = mor.getSpeciesIDNameMap();
	//SK: define a map object that contains the same information as allSpeciesIDs from mor, but with reversed mappings, i.e., species name to orthology id, and not the reverse
	map<string,int> allSpeciesIDs_r;
	for(map<int,string>::iterator mIter=allSpeciesIDs.begin();mIter!=allSpeciesIDs.end();mIter++)
	{
		allSpeciesIDs_r[mIter->second]=mIter->first;
	}
	//SK: define mapping of ids in the orthology for ourr subset of species of interest for the analysis
	vector <int> MapS(NSpc,-1); //define bector to map the species of interest in this orthogroup 
	for(int i=0;i<NSpc;i++)
	{
		int m=allSpeciesIDs_r.find(specOrder[i])->second;
		MapS[i]=m;
	}
	//SK: set number of species represetned in the orthology
	int NAllSpecies=0; //set number of annotated species in orthology
	NAllSpecies=mor.getNAllSpecies();
	//SK: loop over all orthogroups 
	for(map<int,MappedOrthogroup*>::iterator oIter=ogSet.begin();oIter!=ogSet.end();oIter++)
	{
		//cout << "Checking Orthogroup " << oIter->first << endl;
		//SK: get orthogroup information
                MappedOrthogroup* grp = oIter->second;
		//leave orthogroup if there are fewer than two members
                if(grp->getOrthoMembers().size()<minSpcCnt)
		{
                        cout << "To few members: Leaving ortho group " << grp->getID() << endl; //check if there are at lest minSpcCnt genes in the orthogroup, if not skip and this orthogroup won't be represented in the output. 
                        continue;
                }
                //else cout << grp->getOrthoMembers().size() << endl;
		//SK: get gene set information if we are proceeding with this orthogroup
                map<int,map<string,string>*>& GS=grp->getGeneSets();
		//SK: get number of duplication levels in this orthogroup
                int DL=GS.size();
                vector <string> SpcList; 
                //SK: define array object to list the gene names associated  with each duplication level of the orthology. This includes the gene names for all species in the orthology not only the species of interest
		vector < vector <string> > GenesInOGID(DL,vector <string>(NAllSpecies));
		//SK: loop over gene set information to fill GenesInOGID
                for(map<int,map<string,string>*>::iterator gsIter=GS.begin();gsIter!=GS.end();gsIter++) //loop over orthogroup information to expreat the list of species and genes that are in the OG
		{
			//SK: define current interer that represetnts teh curent duplication level we're working on 
                        int DLv=gsIter->first;
                        map<string,string>* SET=gsIter->second;
                        for(map<string,string>::iterator igsIter=SET->begin();igsIter!=SET->end();igsIter++)
			{
                                //SpcList.push_back(igsIter->first);
				//SK: get the species id for the species assiciated with the current member gene on this duplication level
				int spc=allSpeciesIDs_r.find(igsIter->first)->second;
				//SK: fill gene name array 
                                GenesInOGID[DLv][spc]=igsIter->second;
				//cout << igsIter->first << "\t" << igsIter->second << endl;
                        }
                }//end of loops
		//SK: have now gotten gene names for the duplication levels of this current orthogroup for all species in this orthology. 
		//SK: now check how mny of out species of interest have data on each orthogroup duplication level and write out information if there are at least 2. 
		//SK: now loop over the duplication levels of this orthogroup
		for(int d=0;d<DL;d++) 
		{
			//SK: variable to define the number of genes in species of interest with data. 	
                        int cov=0; 
			//SK: loop over all species 
                        for(int s=0;s<NSpc;s++)
			{
				//SK: proceed only if there is a non null genename for this species 
				if(!GenesInOGID[d][MapS[s]].empty())
				{
					//SK: get the species gene name for those duplication level (d) of this orthogroup (oIter->first)
					string genename = GenesInOGID[d][MapS[s]];	
					//SK:check if that gene has data for this species data set, and count it if it does.	
					if(MAPS.find(specOrder[s])!=MAPS.end() && MAPS.find(specOrder[s])->second.find(genename)!=MAPS.find(specOrder[s])->second.end()); 
					//check if that gene has data assocaited with it
					{
						cov++;
						//cout << "Found " << genename << " in " << specOrder[s] << endl;
					}
				}
                        }
			//SK: break the loop if this duplication level as fewer than the required minimum number of genes. This means that an orthogroup with multiple duplication levels will have as many entries in the output merged data file as there are duplication levels, starting from no duplication, that have at least the minimum number of genes with data.
			//SK: as genes tend to become sparser with increasing duplication level, the lower duplication levels with a sufficient number of genes are defined in the output and potentially have their own cluster assignment.	
                       if(cov<minSpcCnt)
			{
				continue;
			}
			//SK: vector object to hold the merged data information
			vector <string> entries;
			//SK: define default entry for missing data values, for genes that are not in the orthology or that have no data in a given species.
			string in;
			in.append("<nodata>");
			//SK: keeping it in this format in case there comes a time we want to write the merged data file with <nodata> for missing entries, instead of the 
			//SK: loop over our species of interest again.
                        for(int s=0;s<NSpc;s++) //loop over all species
			{
				//SK: get gene name for this species on this duplication level 
				string genename = GenesInOGID[d][MapS[s]];
				//SK: check if this gene has data and if so add that data to the entries vector for this duplication level .
                                if(MAPS.find(specOrder[s])->second.find(genename)!=MAPS.find(specOrder[s])->second.end()) //if it is a species we are looking at and it is a gene with data do the following
				{
					vector <string> data=MAPS.find(specOrder[s])->second.find(genename)->second;
					for(int i=0;i<data.size();i++)
					{
						entries.push_back(data[i]);
					}
				}
				else //otherwise if it is a species of interest but there is no data for this gene, use the default entry for the mitting gene/data for this species in this like of the meregd data file.
				{
                                        for(int c=0;c<dataCnts.find(specOrder[s])->second;c++)//add as many entries at there are data measuresments
					{
						entries.push_back(in);
                                        }
                                }
                        }
			//SK: a variable to hold the sum of the availalbe data entries for the genes with data on this duplication level 
			double total=0;
			//SK: a variable to hold the value for the number of missing data entries 
			int numEmpty=0;
			//SK: loop over the entries vector and count the number of missing and non-missing entries 
			for(int i=0;i<entries.size();i++)
			{
				if(entries[i]!=in)
				{
					total += atof(entries[i].c_str()); //sum up the total of the non-missing data entries
				}
				else
				{
					numEmpty++; //sum the number of missing data entries. 
				}
			}
			//SK: define the number data entries 
			int dat=entries.size()-numEmpty;
			//SK: dfine the average value of the data entries
			double ave=(double)total/(double)dat; //calculate the average of the data entries.
			//SK: define a double vector object
			vector <double> values;
			//SK: loop over all of the entries and fill the values vector, using the average for the mising values 
			for(int i=0;i<entries.size();i++) //write out the data for this orthogroup duplication level as numbers 
			{
				if(entries[i]!="<nodata>")//if it is not a missing value, than use the data value
				{
					values.push_back(atof(entries[i].c_str()));
				}
				else //if it is a missing value then insert the average value for this orthogroup-duplication level
				{
					values.push_back(ave);
				}
			}
			//SK: Now write out the data for this duplication level using the values learned.  Here the OG id and duplication level define the "name" of this entry in the merged data file
			oFile << "OG" << grp->getID() << "_" << d;
			//SK: set variable for the number of measurements per gene.
			expCntPerGene=values.size();
			for(int i=0;i<values.size();i++)
			{
				oFile << "\t" << values[i];
			}
			oFile << endl;
		}//SK: endloop for this duplication level 
	}//SK: end loop for this orthogroup
	cout << "Exiting non-source species merging function." << endl;
}//SK: finished with function

int
Framework::fgconverter(const char* geneexpFName,const char* outsuffix,int logTrans,int dataStart)
{
	geneMgr.readGeneData(geneexpFName,logTrans,dataStart);
	strcpy(outFName,outsuffix);
        sprintf(pnlGraphFName,"%s.model",outsuffix);
        //expCntPerGene=aCnt;
        vType=CONT;
	BioNetwork bnw;
        bnw.setProteinManager(&protMgr);
        bnw.setGeneManager(&geneMgr);
        bnw.setPPInteractionManager(&ppMgr);
        bnw.setPDInteractionManager(&pdMgr);
        bnw.createNetwork();
        map<int,int> measuredGeneIDs;
        map<string,int> pdIntr;
        map<string,int> ppIntr;
        //Get data for all experiments
        for(int i=0;i<expCntPerGene;i++)
        {
                expIds.push_back(i);
        }
        bnw.getAllNodeIDs(measuredGeneIDs);
        createGraph(measuredGeneIDs, pdIntr, ppIntr);
        phyNw.doTopologicalSort();
        vector<Node*>& topOrder=phyNw.getTopologicalSort();
        writeToFile(topOrder);
        writeToFile_MeanStd(topOrder);
        phyNw.genPNLInputFormat(pnlGraphFName);
	return 0;
}

int
Framework::createGraph(map<int,int>& nodeIds,map<string,int>& pdIntIds, map<string,int>& ppIntIds)
{
        vector<int> expLevels;
        expLevels.push_back(BioNode::LOW);
        expLevels.push_back(BioNode::MEDIUM);
        expLevels.push_back(BioNode::HIGH);
        vector<int> edgeVals;
        edgeVals.push_back(0);
        edgeVals.push_back(1);
        //For every node in the nodeId add a protein and gene node and make the protein a child of the gene node
        for(map<int,int>::iterator aIter=nodeIds.begin();aIter!=nodeIds.end();aIter++)
        {
                Gene* aGene=geneMgr.getGeneNode(aIter->first);
                const char* geneName=aGene->getName();
                phyNw.addNode(geneName,Node::GENE,aIter->first,expLevels);
                Protein* aProt=protMgr.getProteinNode(aIter->first);
                if(aProt!=NULL)
                {
                        const char* proteinName=aProt->getName();
                        phyNw.addNode(proteinName,Node::PROTEIN,aIter->first,expLevels);
                        phyNw.addEdge(geneName,proteinName);
                }
        }
        //For every pdna interaction add a node for the static attribute, another node for dynamic 
        //attribute which is the child of the first. Then make this node and the protein node
        //as the parent of the child node
        for(map<string,int>::iterator aIter=pdIntIds.begin();aIter!=pdIntIds.end();aIter++)
        {
                char sName[256];
                sprintf(sName,"%spg.I",aIter->first.c_str());
                //Here the ids in the source data don't make sense
                phyNw.addNode(sName,Node::STAT_PG,-1,edgeVals);
                char dName[256];
                sprintf(dName,"%spg.N",aIter->first.c_str());
                phyNw.addNode(dName,Node::DYN_PG,-1,edgeVals);
                phyNw.addEdge(sName,dName);

                //Now get the protein node and the gene node and connect them
                int protId=0;
                int geneId=0;
                int srcPos=0;
                int destPos=0;
                char tempStr[256];
                const char* idPtr=aIter->first.c_str();
		while(idPtr[srcPos]!='\0')
                {
                        if(idPtr[srcPos]=='-')
                        {
                                tempStr[destPos]='\0';
                                protId=atoi(tempStr);
                                destPos=0;
                        }
                        else
                        {
                                tempStr[destPos]=idPtr[srcPos];
                                destPos++;
                        }
                        srcPos++;
                }
                tempStr[destPos]='\0';
                geneId=atoi(tempStr);
                Gene* gNode=geneMgr.getGeneNode(geneId);
                Protein* pNode=protMgr.getProteinNode(protId);
                //Make dName as a parent of gNode
                phyNw.addEdge(dName,gNode->getName());
                //Make pNode as a parent of dName 
                phyNw.addEdge(pNode->getName(),dName);
        }
        //For every protein-protein interaction do something similar
        for(map<string,int>::iterator aIter=ppIntIds.begin();aIter!=ppIntIds.end();aIter++)
        {
                char sName[256];
                //Here the ids in the source data don't make sense
                sprintf(sName,"%spp.I",aIter->first.c_str());
                phyNw.addNode(sName,Node::STAT_PP,-1,edgeVals);
                char dName[256];
                sprintf(dName,"%spp.N",aIter->first.c_str());
                phyNw.addNode(dName,Node::DYN_PP,-1,edgeVals);
                phyNw.addEdge(sName,dName);

                //Now get the protein nodes
                int protId1=0;
                int protId2=0;
                int srcPos=0;
                int destPos=0;
                char tempStr[256];
                const char* idPtr=aIter->first.c_str();
		while(idPtr[srcPos]!='\0')
                {
                        if(idPtr[srcPos]=='-')
                        {
                                tempStr[destPos]='\0';
                                protId1=atoi(tempStr);
                                destPos=0;
                        }
                        else
                        {
                                tempStr[destPos]=idPtr[srcPos];
                                destPos++;
                        }
                        srcPos++;
                }
                tempStr[destPos]='\0';
                protId2=atoi(tempStr);

                Protein* pNode1=protMgr.getProteinNode(protId1);
                Protein* pNode2=protMgr.getProteinNode(protId2);
                //Make pNode1 and pNode2 as a parents of dName
                phyNw.addEdge(pNode1->getName(),dName);
                phyNw.addEdge(pNode2->getName(),dName);
        }
}

int
Framework::writeToFile(vector<Node*>& topOrder)
{
        gsl_rng* g_rng=gsl_rng_alloc(gsl_rng_default);
        char exprFName[1024];
        sprintf(exprFName,"%s.data",outFName);
        ofstream oFile(exprFName);
        for(int eId=0;eId<expIds.size();eId++)
        {
                if(eId==167)
                {
                        cout <<"Stop here " << endl;
                }
                int shown=0;
                for(int i=0;i<topOrder.size();i++)
                {
                        //Get the node here
                        Node* anode=topOrder[i];
                        //Depending upon the node type take appropriate action
                        Node::NodeType nType=anode->getNodeType();
                        switch(nType)
                        {
                                case Node::GENE:
                                {
                                        int idInData=anode->getIdInData();
                                        Gene* aGene=geneMgr.getGeneNode(idInData);
                                //      int dval=aGene->getDiscreteExpLevelAt(expIds[eId]);
                                        //double cval=aGene->getRankedExpLevelAt(expIds[eId]);
                                        double cval=aGene->getExpLevelAt(expIds[eId]);
                                        if(isnan(cval))
                                        {
                                                cout <<"Skipping " << aGene->getName() << endl;
                                        }
                                        else
                                        {
                                                if(shown>0)
                                                {
                                                        oFile <<"\t";
                                                }
                                                if(vType==CONT)
                                                {
                                                        oFile << i<<"=["<< cval << "]" ;
                                                }
                                                shown++;
                                        }
                                        break;
                                }
				case Node::PROTEIN:
                                {
                                        int idInData=anode->getIdInData();
                                        Protein* aProtein=protMgr.getProteinNode(idInData);
                                        //int dval=aProtein->getDiscreteExpLevelAt(expIds[eId]);
                                        //double cval=aProtein->getRankedExpLevelAt(expIds[eId]);
                                        double cval=aProtein->getExpLevelAt(expIds[eId]);
                                        if(i>0)
                                        {
                                                oFile <<"\t";
                                        }
                                        if(vType==CONT)
                                        {
                                                oFile << i<<"=["<< cval << "]" ;
                                        }
                                        else
                                        {
                                //              oFile << i<<"=["<< dval << "]" ;
                                        }
                        /*              if((aProtein->getID()==geneID) || (aProtein->getID()==87))
                                        {
                                                cout << " "<< aProtein->getName() << ": " << aProtein->getExpLevelAt(expIds[eId]);
                                        }*/
                                        break;
                                }
                                case Node::STAT_PG:
                                case Node::STAT_PP:
                                {
                                        //Use the multinomial distribution to set a value here
                                        if(i>0)
                                        {
                                                oFile <<"\t";
                                        }
                                        oFile << i<<"=[0|"<<mult_dist[1]
                                                         <<",1|" <<mult_dist[0]<<"]";
                                        break;
                                }
                                case Node::DYN_PP:
                                case Node::DYN_PG:
				{
                                        //Set to 1 with 50-50 chance. It really does not matter
                                /*      if(rNo<0.5)
                                        {
                                                oFile <<" " << 1;
                                        }
                                        else
                                        {
                                                oFile <<" " << 1;
                                        }*/
                                        break;
                                }
                        }

                }
                oFile << endl;
        }
        oFile.close();
	return 0;
}

int
Framework::writeToFile_MeanStd(vector<Node*>& topOrder)
{
        char exprFName[1024];
        sprintf(exprFName,"%s.meansd",outFName);
        ofstream mFile(exprFName);
        for(int i=0;i<topOrder.size();i++)
        {
                Node* anode=topOrder[i];
                int idInData=anode->getIdInData();
                Gene* aGene=geneMgr.getGeneNode(idInData);
                mFile <<aGene->getName() <<"\t" << aGene->getMean() <<"\t" << aGene->getStdev() << endl;
        }
        mFile.close();
        return 0;
}

int
Framework::learnMoE(const char* prefix,const char* clusteringOutputDir,int k)
{
	evManager.setVariableManager(&varManager);
	initTypeL=MotifRegressor::RAND;
	char fName[256];
	sprintf(fName,"%s.model",prefix);
	cout << fName << endl;
	varManager.readVariables(fName);
	sprintf(fName,"%s.data",prefix);
	cout << fName << endl;
	evManager.loadEvidenceFromFile_Continuous(fName);
	//set initialization type for random initialization
	MotifRegressor mlearner;
	mlearner.setInitType(initTypeL);
	//if(initType==MotifRegressor::GMM)
	//{
	//	mlearner.setInitClusterFile(gmmClusterFName);
	//}
	mlearner.setExpertCnt(k);
	mlearner.setVariableManager(&varManager);
	mlearner.setEvidenceManager(&evManager);
	mlearner.setMotifManager(&motifManager);
	char makeDirCmd[1024];
	//SK: set the system varible for creating the output directory
	sprintf(makeDirCmd,"mkdir -p %s",clusteringOutputDir);
	//SK: run the system command to create the output directory
	system(makeDirCmd);
	cout << "Making clustering directory " <<  clusteringOutputDir << endl;
	mlearner.setOutputDir(clusteringOutputDir);
	mlearner.setUntransformedData(untransformedData);
	//mlearner.learnMoE();
	cvfolds=1;
	mlearner.learnMoE_CrossValidation(cvfolds);
	//mlearner.showMoE(outputDir);
        //mlearner.showMoEParameters();
        //mlearner.showClusterAssignment(outputDir);
        mlearner.dispTFsPerCluster();
        //mlearner.showGenatomyModule();
        if(testData.getSize()>0)
        {
                mlearner.predictTestData(testData.getDataSet());
        }
	return 0;
}

int
Framework::genSpeciesClusters(const char* aFName,const char* outDir)
{
	//SK: read in initial assignments
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
                string geneName;
                int clusterid=0;
                while(tok!=NULL)
                {
                        if(tokCnt==0)
                        {
                                geneName.append(tok);
                        }
                        else if(tokCnt==1)
                        {
                                clusterid=atoi(tok);
                        }
                        tok=strtok(NULL,"\t");
                        tokCnt++;
                }
                geneClusterAssignment[geneName]=clusterid;
        }
        inFile.close();
	//SK: write out same assignments
	 int excludeID=0;
	map<string,ofstream*> filePtr;
	for(map<string,int>::iterator spc=speciesList.begin();spc!=speciesList.end();spc++)
	{
                string specName=spc->first;
                if(specName.length()>0)
                {
			char bFName[1024];
                        sprintf(bFName,"%s/%s_initial_clusterassign.txt",(char*)outDir,specName.c_str());
                        ofstream* oFile=new ofstream(bFName);
                        filePtr[specName]=oFile;
                }
        }
	for(map<string,int>::iterator gIter=geneClusterAssignment.begin();gIter!=geneClusterAssignment.end();gIter++)
        {
                MappedOrthogroup* mgrp=mor.getMappedOrthogroup(gIter->first.c_str(),srcSpecies);
                if(mgrp==NULL)
                {
                        cout <<"Null OG for " << gIter->first << endl;
                        continue;
                }
                if((excludedOGList.size()>0) && (excludedOGList.find(mgrp->getID())!=excludedOGList.end()))
                {
                        excludeID++;
                        continue;
                }
                map<string,GeneMap*>& members=mgrp->getOrthoMembers();
                for(map<string,int>::iterator sIter=speciesList.begin();sIter!=speciesList.end();sIter++)
                {
                        GeneMap* spechit=mgrp->getSpeciesHits(sIter->first.c_str());
                        if(spechit==NULL)
                        {
                                continue;
                        }
                        ofstream* file=filePtr[sIter->first];
                        map<string,map<string,STRINTMAP*>*>& specGenes=spechit->getGeneSet();
                        for(map<string,map<string,STRINTMAP*>*>::iterator hIter=specGenes.begin();hIter!=specGenes.end();hIter++)
                        {
				if(hIter->first!="NULL")
				{
                                	(*file)<< hIter->first<<"\t" << gIter->second << endl;
				}
                        }
                }
        }
	return 0;
}

int
Framework::genSpeciesClustersNonSrc(const char* aFName,const char* outDir)
{
	//SK: start by reading in the cluster assignment information for the merged data
	map<string,int> geneClusterAssignment;
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
                string geneName;
                int clusterid=0;
                while(tok!=NULL)
                {
                        if(tokCnt==0)
                        {
                                geneName.append(tok);
                        }
                        else if(tokCnt==1)
                        {
                                clusterid=atoi(tok);
                        }
                        tok=strtok(NULL,"\t");
                        tokCnt++;
                }
                geneClusterAssignment[geneName]=clusterid;
        }
        inFile.close();
	//SK: define map of out file ofstream objects, with one for each species
        map<string,ofstream*> filePtr;
        int excludeID=0;
	//SK: loop over the set of species in the analysis and initialize and initialize the ofstream object for the input cluster assignment file for each species. 
	for(map<string,int>::iterator sIter=speciesList.begin();sIter!=speciesList.end();sIter++)
	{
		 string specName=sIter->first;
                if(specName.length()>0)
                {
                        char bFName[1024];
                        sprintf(bFName,"%s/%s_initial_clusterassign.txt",(char*)outDir,specName.c_str());
                        ofstream* oFile=new ofstream(bFName);
                        filePtr[specName]=oFile;
                }
        }
	//SK: get set of orthogroups 
        map<int,MappedOrthogroup*>& ogSet=mor.getMappedOrthogroups();
	//SK: set minimum number of species to be represented at a given duplicatoin level for a given orthogroup
        int minSpcCnt=2;
	//SK:loop over all orthogroups 
        for(map<int,MappedOrthogroup*>::iterator oIter=ogSet.begin();oIter!=ogSet.end();oIter++)
	{
		//SK: get the orthogroup infofrmation
                MappedOrthogroup* grp = oIter->second;
		//SK: check if there are less than the minimum number of species in this orthogroup
                if(grp->getOrthoMembers().size()<minSpcCnt)
		{
                        //cout << "To few members: Leaving ortho group" << grp->getID() << endl;
                        continue;
                }
                map<string,GeneMap*>& Members=grp->getOrthoMembers();
		//SK: get gene set information 
                map<int,map<string,string>*>& GS=grp->getGeneSets();
		int lastCA=-1;
		int DL_total=GS.size();
		bool first=true;
		int DL_first=-1;
		//SK: if there are more than one partition levels to the gene tree check if the first level with a cluster assignment is not the first _0 level
		if(DL_total>1)
		{
			for(int d=1;d<DL_total;d++)
			{
				char groupid[50];
				sprintf(groupid,"OG%i_%i",grp->getID(),d);
				if(first && geneClusterAssignment.find(groupid)!=geneClusterAssignment.end() && d!=0)
				{
					first=false;
					DL_first=d;
					break;
				}
			}
		}
		//SK: loop over gene set information for each duplication level of the orthogroup
                for(map<int,map<string,string>*>::iterator gsIter=GS.begin();gsIter!=GS.end();gsIter++)
		{
			//SK: define orthogroup level  
                        int DL=gsIter->first;
			//SK: if this orthogroup level is not the first level of this orthogroup represented and there are more than one level to this orthogroup: 
			//SK: this acts to provide a cluster assignment to genes on duplication levels that may be below the first duplication level to be given a cluster assignment from the merged data.
			if(DL_total>1 && DL<DL_first && DL_first!=-1)
			{
				DL=DL_first;
			}
			//SK: define the group id name  as it would appear in the merged clustering output, which will be in the form of  OG<ogid>_<duplication level>
                        char groupid[50];
                        int n = sprintf(groupid,"OG%i_%i",grp->getID(),DL);
                        map<string,string>* SET=gsIter->second;
                        //SK: loop over the gene members of this duplication level of this orthogroup
			for(map<string,string>::iterator igsIter=SET->begin();igsIter!=SET->end();igsIter++)
			{
				//cout << igsIter->first << endl;
				//SK: check if this gene is from a species of interest in our analysis, continue if it is not
				if(speciesList.find(igsIter->first)==speciesList.end() && igsIter->second!="NULL")
				{
					continue;
				}
				//SK: if this is a gene for a species of interest then prepare to write out this information by calling the pointer for the ofstream for the cluster assignment file for this specoes. 
                                ofstream* file=filePtr.find(igsIter->first)->second;
				//SK: check if this group/partition level appears in the clusterassignment information for the merged data, and write out the initial cluster assignment for this gene in this species. 
                                if(geneClusterAssignment.find(groupid)!=geneClusterAssignment.end())
				{
                                        (*file) << igsIter->second << "\t" << geneClusterAssignment.find(groupid)->second << endl;
					lastCA=geneClusterAssignment.find(groupid)->second;
                                }
				//SK: otherwise if this duplication level isn't in the merged data output and we have a gene in a species of interest here, then give that gene the last cluster assignment from the highest duplication level in this orthogroup that has a cluster assignment.
				else if(lastCA!=-1)
				{
					(*file) << igsIter->second << "\t" << lastCA;
				}
                        }//SK: end loop over gene set
                }//SK: end look over duplication levels
        }//SK: end loop over duplication levels
}

int
Framework::setPreClustering(bool in)
{
	preClustering=in;
	return 0;
}

// FIXED COVARIANCE START //
int
Framework::setConstCov(double val)
{
        scMgr.setConstCov(val);
        return 0;
}
// FIXED COVARIANCE END //

int
Framework::setPredictionInputMode(bool in)
{	
	predictionInputMode=in;
	return 0;
}

int
Framework::countOGIDSNonSrc(const char* aFName)
{
        //SK: start by reading in the cluster assignment information for the merged data
        ifstream inFile(aFName);
        char buffer[1024];
 	bool first=true;
        while(inFile.good())
        {       
                inFile.getline(buffer,1023);
                if(strlen(buffer)<=0)
                {       
                        continue;
                }
                char* tok=strtok(buffer,"\t");
                int tokCnt=0;
		if(first)
		{
			first=false;
			continue;
		}
                string Name;
                while(tok!=NULL)
                {       
                        if(tokCnt==0)
                        {       
                                Name.append(tok);
				size_t pt1=Name.find("G");
				size_t pt2=Name.find("_");
				string ogid=Name.substr(pt1+1,pt2-pt1-1);
				int id=atoi(ogid.c_str());
				mergedOgidSet[id]=0;
                        }
                        tok=strtok(NULL,"\t");
                        tokCnt++;
                }
        }
        inFile.close();
	return 0;
}

int
main(int argc, char *argv[])
{
	int opt=0;
	//SK: variable for the file name of the species order list.
	char speciesOrderFName[1024];
	//SK: variable for the OGIDS file name
	char OGIDSFName[1024];
	//SK: variable for the species tree file name             
	char treeFName[1024];
	//SK: variable for the configuration file with the input cluster assignment and data files for each of the species 
	char confFName[1024];
	//SK: variable for the randomization option.
	char randOpt[30];
	//SK:  variable for the name of output directoy. 
	char outputDir[1024];
	//SK: variable for the mode in which the program is supposed to be run.
	char mode[30];
	//SK: variable for the name of the species taht is choosen to be the base species. 
	char sourceSpecies[1024];		
	//SK: variable to define the type of initialization that will be used for the branch lengths of the species tree. 
	char initType[30];
	//SK: the char vaiable to define the branch length information. 
	//SK: Either will be turned into a number or will represent a file name depending on the type of initialization that is requested, which will be defined with the above variable, init type.
	char temp_p_diagonal[1024];
	//SK: the variable containing the directory path of the genetree .tre files
	char geneTreeDir[1024];
	//SK: the variable for the number of clusters to be ised in the analysis, intilized to -1
	int kClusters=-1;
	//SK: the varible to set the number of maximum iterations to be used in the analysis, initialized to the default value in the online code.
	int nMaxIterations=50;
	//SK: The variable to set the number of "folds" of the data use for the crossvalidation analysis.
	int crossValFold=0;
	//SK: The variable to define the threshold for the change in score that will define convergence for the learning algorithm, initilized to the value in the online code. 
	double convThresh=0.5;
	//SK: The variable to define if the output directory can be overwritten, initialized to false so that the arboretum output directory will never be spuriously overwritten. 
	bool overWrite=false;
	//SK: variable for the option to run the second optimization stage in SpeciesClusterManager::dumpAllInferredClusterAssignments(const char* outputDir) or not
	bool secondStage=true;
	bool preClusteringStage=false;
	bool sourceInit=false;
	bool predictionInputMode=false;
	//SK: boolean variables for identifying if the required arguments have been set;
	bool sDefault=false;
	bool eDefault=false;
	bool kDefault=false;
        bool tDefault=false;
	bool cDefault=false;
        bool oDefault=false;
        bool mDefault=false;
        bool iDefault=false;
	bool pDefault=false;
	bool rDefault=false;
        bool bDefault=false;
	bool xDefault=false;    // FIXED COVARIANCE //
	double fixedCov=0.2;    // FIXED COVARIANCE //
	//SK: define argument option
	static struct option long_options[] = {
		{"species-order", 	required_argument, 0,  's' },
		{"ortho-mapping", 	required_argument, 0,  'e' },
		{"k-clusters",    	required_argument, 0,  'k' },
		{"species-tree",   	required_argument, 0,  't' },
		{"clustering-input",	required_argument, 0,  'c' },
		{"randomization",       required_argument, NULL,  'r' },
		{"output",       	required_argument, 0,  'o' },
		{"mode",          	required_argument, 0,  'm' },
		{"base-species",        required_argument, NULL,  'b' },
		{"init",                required_argument, 0,  'i' },
		{"prob",        	required_argument, 0,  'p' },
		{"gene-trees",       	optional_argument, 0,  'g' },
		{"conv-thresh",         optional_argument, 0,  'd' },
		{"num-iter",            optional_argument, 0,  'n' },
		{"validation-fold",	optional_argument, 0,  'v' },
		{"write",               optional_argument, 0,  'w' },
		{"first-clustering",      optional_argument, 0,  'f' },
		{"last-stage",        optional_argument, 0,  'l' },
		{"update-source",		optional_argument, 0,  'u' },
                {"constant covariance", optional_argument, 0,  'x' },   // FIXED COVARIANCE //
		{"help",		optional_argument, 0,  'h' },
		{0,0,0,0}
	};
	int optret='-';
        opterr=1;
        int oldoptind=optind;
        int condCnt=1;
	int long_index=0;
	while ((opt=getopt_long(argc,argv,"s:e:k:t:c:r:o:m:b:i:p:g:d:n:v:w:f:l:u:x:h:",long_options,&long_index)) != -1)  // FIXED COVARIANCE; "x"
	{
		//cout << long_index << endl;
		//cout << opt << optret << endl;
		 if(optret=='?'){
                        cout <<"Option error " << optopt << endl;
                        return -1;
                }
                char c;
                char* my_optarg=NULL;
                c=*(argv[oldoptind]+1);
                if(optind-oldoptind==2){
                        my_optarg=argv[oldoptind+1];
                }
                else{
                        my_optarg=argv[oldoptind]+2;
                }
                //cout << "beginning: " << my_optarg << ", " << c << ", " << optind << ", " << oldoptind << endl;
		switch (opt)
		{
			case 's':
			{
				sDefault=true;
				//SK: if it is, then look at the next string in the argument and copy that string to the SpeciesOrderFName variable.
				strcpy(speciesOrderFName,my_optarg);
				//SK: print out the name of the species order list file name. This is to document in the input argument that has been read
				cout << "Species list file: " << speciesOrderFName << endl;
				break;
			}
                        // FIXED COVARIANCE START //
			case 'x':
                        {
                                fixedCov = atof(my_optarg);
                                xDefault = true;
                                cout << "Fixed covariance is: " << fixedCov << endl;
                                break;
                        }
                        // FIXED COVARIANCE END //
			case 'e':
			{
				eDefault=true;
				//SK: take the next string in the command arguments and set the OGIDS file name variable
				strcpy(OGIDSFName,my_optarg);
				//SK: print out the variable for the orthomapping file name
				cout << "Ortho-mapping file: " << OGIDSFName << endl;
				break;
			}
			case 'k':
			{
				kDefault=true;
				//SK: print out the number of clusters for the output
				cout << "Number of clusters: " << my_optarg << endl;
				//SK: set the kClusters variable
				kClusters=atoi(optarg);
				break;
			}
			case 't':
			{
				tDefault=true;
				//SK: copy the tree file name to the treeFName variable
				strcpy(treeFName,my_optarg);
				//SK: printout the species tree file name
				cout << "Species Tree file: " << treeFName << endl;
				break;
			}
			case 'c':
			{
				cDefault=true;
				//SK: take the next string in the argument set and copy it to the confFName variable
				strcpy(confFName,my_optarg);
				//SK: print out the configuration file name, used for calling the input CA and data files
				cout << "Cluster assignment and data input file: " << confFName << endl;
				break;
			}
			case 'r':
			{ 
				rDefault=true;
				//SK: take the next string in the arguments from the command and set the randOpt variable
				strcpy(randOpt,my_optarg);
				//SK: print out the randomization option that has been set
				cout << "Randomization option: " << randOpt << endl;
				break;
			}
			case 'o':
			{
				oDefault=true;
				//SK: copy next string in the arguments from the command and set the outputDir directory.
				strcpy(outputDir,my_optarg);
				//SK: print out the name of the output directory that has been given in the command
				cout << "Output directory: " << outputDir << endl;
				break;
			}
			case 'm':
			{
				mDefault=true;
				//SK: copy the next string to the mode variable
				strcpy(mode,my_optarg);
				//SK: print out the mode option that has been set.
				cout << "Mode: " << mode << endl;
				break;
			}
			case 'b':
			{
				bDefault=true;
				//SK: take the next string in the arguments from the command and set the sourceSpecies argument
				strcpy(sourceSpecies,my_optarg);
				//SK: print out the source species variable
				cout << "Base species:" << sourceSpecies << endl;
				break;
			}
			case 'i':
			{
				iDefault=true;
				//SK: take the next string in the argument and set the initType variable
				strcpy(initType,my_optarg);
				//SK: print out the setting of the initialization type variable                            
				cout << "Initialization type: " << initType << endl;
				break;
			}
			case 'p':
			{
				pDefault=true;
				//take the next string from the arguments from the command and set it to the temp_o_diagonal variable
				strcpy(temp_p_diagonal,my_optarg);
				//SK: print out this variable. Depending on the init-type variable, this will be a number or a file name
				cout << "p diagonal value: " << temp_p_diagonal << endl;
				break;
			}
			case 'g':	
			{
				//SK: copy the next string in the argument to the geneTreeDirDirectory, which is the name of the directory for the gene trees
				strcpy(geneTreeDir,my_optarg);
				//SK: print out the name of the gene tree directory.
				cout << "Gene tree directory: " << geneTreeDir << endl;
				//SK: set the global directory name variable
				GENETREEDIR=geneTreeDir;
				//SK: set the boolean variable to say that the gene trees should be read in from files
				GENFROMFILE=true;
				break;
			}
			case 'd':
			{
				//SK: set the convergence threshold varibale from the next string in the argument
				convThresh=atof(my_optarg);
				//print out the convergence threshold setting, was initialized to 0.05
				cout << "Converstion threhsold: " << convThresh << endl;
				break;
			}
			case 'n':
			{
				//SK: set the maximum number of iterations with the next string in the argument set
				nMaxIterations=atoi(my_optarg);
				//SK: print out the maximum number of iterations setting, was initialized to 50
				cout << "Max number of iterations: " << nMaxIterations << endl;
				break;
			}
			case 'v':
			{
				//SK: set the number of cross validation folds with the next string in the command
				crossValFold=atoi(my_optarg);
				//SK: print out the number of folds desired for the cross validation work
				cout << "Cross validation parameter: " << crossValFold << endl;
				break;
			} 
			case 'w':
			{
				//SK: take the next string in the command and set the overWrite variable
				if(strcmp(my_optarg,"true")==0)overWrite=true;
				else if(strcmp(my_optarg,"false")==0)overWrite=false;
				// //SK: print out the status of this variable, true or false.
				cout << "Overwrite permission: " << my_optarg << endl;
				break;
			}
			case 'f':
                        {
                                //SK: take the next string in the command and set the overWrite variable
                                if(strcmp(my_optarg,"true")==0)preClusteringStage=true;
                                if(strcmp(my_optarg,"false")==0)preClusteringStage=false;
                                //SK: print out the status of this variable, true or false.
                                cout << "Pre-clustering stage: " << my_optarg << endl;
				break;	
                        }
			case 'l':
                        {
                                //SK: take the next string in the command and set the overWrite variable
                                if(strcmp(my_optarg,"true")==0)secondStage=true;
                                if(strcmp(my_optarg,"false")==0)secondStage=false;
                                //SK: print out the status of this variable, true or false.
                                cout << "Second stage permission: " << my_optarg << endl;
				break;
                        }
			case 'u':
			{
				//SK: take the next string in the command and set the overWrite variable
                                if(strcmp(my_optarg,"true")==0)sourceInit=true;
                                if(strcmp(my_optarg,"false")==0)sourceInit=false;
                                //SK: print out the status of this variable, true or false.
                                cout << "Source species initialization: " << my_optarg << endl;
                                break;
			}
			case 'h':
			{
                                print_usage();
                                return 0;
                        }
			case '?':
			{
        			printf("missing opt: %c\n",opt);
				return 0;
			}
			case ':':
			{
				printf("missing arg: %c\n",my_optarg);
				return 0;
			}
			default: 
			{
				print_usage();
                 		exit(EXIT_FAILURE);
			}
		}//SK: exit switch loop
		 oldoptind=optind;
	}//SK: exit while
	//SK: the following lines check if any of the required input arguments have not been set and prints a statement about what argument is missing and then also the overall usage summary information that is availalbe with the -h or --help information. This is meant to assist the inexperienced user in troublshooting the arguments that may be missing. 
	if(sDefault==false)
	{
		cout << "-s species order list file was not defined" << endl;
		print_usage();
		exit(EXIT_FAILURE);
	}
	if(eDefault==false)
        {
                cout << "-e gene orthology file was not defined" << endl;
		print_usage();
		exit(EXIT_FAILURE);
        }
	if(tDefault==false)
        {
                cout << "-t species tree file was not defined" << endl;
		print_usage();
		exit(EXIT_FAILURE);
        }
	if(sDefault==false)
        {
                cout << "-k number of clusters was not defined" << endl;
		print_usage();
		exit(EXIT_FAILURE);
        }
	if(cDefault==false)
        {
                cout << "-c configuration file was not defined" << endl;
		print_usage();
		exit(EXIT_FAILURE);
        }
	if(rDefault==false)
        {
                cout << "-r randomization option was not defined" << endl;
		exit(EXIT_FAILURE);
        }
	if(oDefault==false)
        {
                cout << "-o output directory path was not defined" << endl;
		print_usage();
		exit(EXIT_FAILURE);
        }
	if(mDefault==false)
        {
                cout << "-m program run mode was not defined" << endl;
		print_usage();
		exit(EXIT_FAILURE);
        }
	if(iDefault==false)
        {
                cout << "-i initialization method  option was not defined" << endl;
		print_usage();
		exit(EXIT_FAILURE);
        }
        if(pDefault==false)
        {
                cout << "-p initial transition probability (p_diagonal) not defined" << endl;
		print_usage();
		exit(EXIT_FAILURE);
        }
	//SK: check if required paramters were set.
	//SK: verify the number of clusters is valid, greater than 1
	if(kClusters<=1)
	{
		//SK: print warning and exit program if it is not
		cout << "Error: invalid number of clusters; exiting program." << endl;
		return 0;
	}
	//SK: check output directory status
	//SK: define stat object
	struct stat statbuf;
	//SK: use stat object
	stat(outputDir, &statbuf);
	//SK: check status of the output directory
	int dirCheck=S_ISDIR(statbuf.st_mode);
	//SK: if the directory doesn't exist enter this if statement and create the output directory as needed
	if(!dirCheck) //SK: if the directory doesn't exist enter this if statement
	{
		//SK: define system command variable
		char outputDirCmd[1024];
		//SK: set the system varible for creating the output directory
		sprintf(outputDirCmd,"mkdir -p %s",outputDir);
		//SK: run the system command to create the output directory
		system(outputDirCmd);
		//SK: state that the output directory is being created.
		cout << "Making output directory " <<  outputDir << endl;
	} 
	//SK: check if the directory exists and if any overwrite is forbidden
	else if(dirCheck && overWrite==false)
	{
		//SK: state that the choosen output directory can't be overwritten and exit the program
		cout << "Error: output directory exists and the (over) write (-w) option is false; exiting program." << endl;
		return 0;
	} 
	//SK: check if the directory exists and if that directory can be overwritten based on the permissions set in the commmand arguments
	else if(dirCheck && overWrite==true)
	{
		//SK: print out warning mesasge that the choosen directory will be written over
		cout << "Note: Now overwriting " << outputDir << "; -w option set to true." << endl;
		//set char string varible for command to remove the directory
		char outputDirCmd1[1024];
		//SK: write to the first commandstring varible, to remove the given directory 
		sprintf(outputDirCmd1,"rm -Rf %s/*",outputDir);
		//SK: run the system rm command
		system(outputDirCmd1);
		//SK: create varible for second command to create empty directory 
		char outputDirCmd2[1024];
		//SK: write command to second variable.
		sprintf(outputDirCmd2,"mkdir -p %s",outputDir);
		//SK: run the system command to create a new, empty version of the directory.
		system(outputDirCmd2);
	}
	//proceed with the analysis, these lines are little changed from the initial version of the code, but I'll make some comments just to be careful
	//SK: create master framework object
	Framework fw;
	//SK: set the option for running the second optimization stage in SpeciesClusterManager::dumpAllInferredClusterAssignments(const char* outputDir)
	fw.setSecondStageOption(secondStage);
	//SK: set the number of iterations for the analysis
	fw.setNMaxIterations(nMaxIterations);
	//SK: set the convergence threshold condition
	fw.setConversionThreshold(convThresh);
	//SK: read in the species tree and number of clusters
	fw.readSpeciesTree(kClusters,treeFName);
	//SK: read in the species order list and mapped orthogroup information to the mor object
	fw.readOrthology(speciesOrderFName,OGIDSFName);
	//SK: seIt the source species name
        fw.setSrcSpecies(sourceSpecies);
	//SK: set option for initializing cluster means from the source species. 
	fw.setSourceInitOption(sourceInit);
	//SK: do preclustering if requested.
        // FIXED COVARIANCE START //
        if (xDefault)
        {
                fw.setConstCov(fixedCov);
        }
        // FIXED COVARIANCE END //
	if(preClusteringStage)
	{
	        cout << "Prep expansion of orthogroup set" << endl;
                char mergedFileNonSource[1024];
                sprintf(mergedFileNonSource,"%s/mergedDataNonSource.txt",outputDir);
                fw.readExpressionDataNonSrc(confFName,mergedFileNonSource);
                fw.countOGIDSNonSrc(mergedFileNonSource);
		fw.setPreClustering(true);
		char mergedFile[1024];
		sprintf(mergedFile,"%s/mergedData.txt",outputDir);
		//fw.readExpressionData(confFName,mergedFile);	// source spc only, 210819
		char outsuffix[1024];
		sprintf(outsuffix,"%s/mergedData",outputDir);
		//fw.fgconverter(mergedFile,outsuffix,1,1);	// source spc only, 210819
		fw.fgconverter(mergedFileNonSource,outsuffix,1,1);
		char clusteringOutputDir[1024];
		sprintf(clusteringOutputDir,"%s/mergedClustering",outputDir);
		fw.learnMoE(outsuffix,clusteringOutputDir,kClusters);
		char clusterFile[1024];
                sprintf(clusterFile,"%s/mergedClustering/fold0/clusterassign.txt",outputDir);
		//fw.genSpeciesClusters(clusterFile,outputDir);	// source spc only, 210819
		fw.genSpeciesClustersNonSrc(clusterFile,outputDir);
	}
	if(strcmp(mode,"prediction")==0)
	{
		predictionInputMode=true;
	}
	if(predictionInputMode)
	{
		cout << "Prepare prediction mode" << endl;
		fw.setPredictionInputMode(true);
		char mergedFile[1024];
		sprintf(mergedFile,"%s/mergedDataAll.txt",outputDir);
		fw.readExpressionDataNonSrc(confFName,mergedFile);
		fw.countOGIDSNonSrc(mergedFile);
		fw.readTransitionMatrices(confFName);
	}
	//SK: check that the program is not being run in cross validation mode
	if(!(strcmp(mode,"crossvalidation")==0))
	{
		//SK:initialize species data and randomization option.
		fw.readSpeciesData(confFName,randOpt,outputDir);
	}
	//SK: check the type of initialization required for the cluster transition probabilities
	if(strcmp(initType,"uniform")==0)
	{
		//SK: if uniform treat the temp_p_diagonal variable as a number and print a comment to that effect
		cout << "Uniform probability value: " << atof(temp_p_diagonal) << endl;
		//SK: set the cluster transition probability with a uniform number
		fw.setClusterTransProb(atof(temp_p_diagonal));
	}
	//SK: if the option is not uniform" but "branchlength", read the initial transition probaility information as a file defining separate transition probabilities 
	else 
	{
		//SK: print out the file name and that this is the initialization method choosen
		cout << "Non-uniform transition probability file:" << temp_p_diagonal << endl;
		 //SK: read in the initial transition matrix probbilities from this file
		fw.setClusterTransProb(temp_p_diagonal);
	}
	//run the program in the selected mode
	//SK: if in learn mode
	if(strcmp(mode,"learn")==0 || strcmp(mode,"prediction")==0)
	{
		//SK: apply the learning algorithm to infer the arboretum model and clusters
		fw.startClustering(outputDir, treeFName);	// repop
	}
	//SK: otherwise, if in crossvalidation mode
	else if(strcmp(mode,"crossvalidation")==0)
	{
		//SK: check that the number of folds to use in the cross validation is greater than 1
		if(crossValFold<1)
		{
			//SK: ifi so print error and exit program
			cout << "Error in main(): The number of cross validation folds is not greater than one; exiting the program." << endl;
			return 0;
		}
		//SK: continue to run the cross validation analysis
		fw.startCrossValidationCheck(outputDir,crossValFold,randOpt,confFName,speciesOrderFName,OGIDSFName,treeFName);	// repop
	}
	//SK: otherwise check if generate mode is requested
	else if(strcmp(mode,"generate")==0)
	{
		//SK: run the generation mode processs
		fw.generateData(outputDir);
	}
	//SK: otherwise assume the program is for the display mode
	else 
	{
		fw.redisplay(outputDir);
	}
	return 0;
}//SK: end of main()
