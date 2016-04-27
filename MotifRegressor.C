
#include <sys/timeb.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include "Variable.H"
#include "Error.H"
#include "VariableManager.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "MotifManager.H"
#include "BFGSWrapperData.H"
#include "BFGSWrapper.H"
#include "ExpertL.H"
#include "Matrix.H"
#include "Kmeans.H"
#include "Distance.H"
#include "ClusterManager.H"
#include "HyperGeomPval.H"
#include "MotifRegressor.H"
#include "GeneExpManager.H"

MotifRegressor::MotifRegressor()
{
	initType=MotifRegressor::RAND;
	foldCnt=1;
	untransformedFName[0]='\0';
}

MotifRegressor::~MotifRegressor()
{
}

int 
MotifRegressor::setVariableManager(VariableManager* aPtr)
{
	varMgr=aPtr;
	return 0;
}

int 
MotifRegressor::setEvidenceManager(EvidenceManager* aPtr)
{
	evMgr=aPtr;
	return 0;
}

int 
MotifRegressor::setMotifManager(MotifManager* aPtr)
{
	motifMgr=aPtr;
}

int
MotifRegressor::setOutputDir(const char* aDirName)
{
	strcpy(outputDir,aDirName);
	return 0;
}

int 
MotifRegressor::setExpertCnt(int expCnt)
{
	expertCnt=expCnt;
	return 0;
}

int
MotifRegressor::setInitType(MotifRegressor::InitType type)
{
	initType=type;
	return 0;
}


int 
MotifRegressor::setInitClusterFile(const char* aFName)
{
	strcpy(clusterFName,aFName);
	return 0;
}

int 
MotifRegressor::learnMoE()
{
	double currScore=0;
	int iter=0;
	bool convergence=false;
	double threshold=0.001;
	r=gsl_rng_alloc(gsl_rng_default);
	//int randseed=16951;
	int randseed=getpid();
	//int randseed=9816;
	cout <<"randseed: " << randseed << endl;
	gsl_rng_set(r,randseed);
	initVarsWithPredictors();
	initPredictorsWithVars();
	if(initType==MotifRegressor::RAND)
	{
		cout <<"Initing experts with rand partition" << endl;
		initExperts();
	}
	else if(initType==MotifRegressor::KMEANS)
	{
		cout <<"Initing experts with kmeans" << endl;
		initExpertsKMeans();
	}
	
	matrixifyData();
	populateBFGSData();
	populateBFGSGammas();
	bfgs.setParamCnt(expertCnt*featIDMatidMap.size());
	bfgs.setFeatureCnt(featIDMatidMap.size());
	bfgs.setStepSize(1e-5);
	bfgs.setTolerance(0.001);
	bfgs.initializeMinimizer();
	while(!convergence && iter<20)
	{
		maximizationStep();
		expectationStep();
		double newScore=getScore();
		if(iter>0)
		{
			double diff=newScore-currScore;
			if(diff<threshold)
			{
				convergence=true;
			}
		}
		currScore=newScore;
		iter++;
		bfgs.reinitializeMinimizer();
	}
	cout <<"Final Score: " << currScore << endl;
		
	return 0;
}


int 
MotifRegressor::learnMoE_CrossValidation(int neededfoldCnt)
{
	foldCnt=neededfoldCnt;
	double threshold=1e-5;
	r=gsl_rng_alloc(gsl_rng_default);
	int randseed=getpid();
	//int randseed=9816;
	cout <<"randseed: " << randseed << endl;
	gsl_rng_set(r,randseed);
	initVarsWithPredictors();
	
	for(int i=0;i<foldCnt;i++)
	{
		split(i);
		initPredictorsWithVars();
		if(initType==MotifRegressor::RAND)
		{
			cout <<"Initing experts with rand partition" << endl;
			initExperts();
		}
		else if(initType==MotifRegressor::KMEANS)
		{
			cout <<"Initing experts with kmeans" << endl;
			initExpertsKMeans();
		}
		else if(initType==MotifRegressor::GMM)
		{
			cout <<"Initing experts with gmm" << endl;
			initExpertsGMM();
		}
		matrixifyData();
		//populateBFGSData();
		//populateBFGSGammas();
		//bfgs.setParamCnt(expertCnt*featIDMatidMap.size());
		//bfgs.setFeatureCnt(featIDMatidMap.size());
		//bfgs.setStepSize(0.5);
		//bfgs.setTolerance(0.1);
		//bfgs.initializeMinimizer();
		bool convergence=false;
		int iter=0;
		double currScore=0;
		while(!convergence && iter<200)
		//while(!convergence && iter<1)
		{
			struct timeval begintime;
			gettimeofday(&begintime,NULL);
			maximizationStep();
			expectationStep();
			double newScore=getScore();
			if(iter>0)
			{
				double diff=(newScore-currScore)/fabs(newScore);
				diff=diff*100;
				if(diff<threshold)
				{
					convergence=true;
				}
				cout <<"Iter: "<< iter <<" DeltaLL " << diff << endl;
			}
			currScore=newScore;
			struct timeval endtime;
			gettimeofday(&endtime,NULL);
			cout << "Time elapsed " << endtime.tv_sec-begintime.tv_sec<< " seconds and " << endtime.tv_usec-begintime.tv_usec << " micro secs" << endl;
			iter++;
			//bfgs.reinitializeMinimizer();
		}
		cout <<"Final Score: " << currScore << endl;
		predictTestExpression();
		dispTFsPerCluster(i);
		double unpenalizedLL=0;
		double testLL=getTestDataLikelihood(unpenalizedLL);
		trainLL[i]=currScore;
		holdoutLL[i]=testLL;
		holdout_unpenLL[i]=unpenalizedLL;
		if(i<foldCnt-1)
		{
			cleanUp();
		}
		
	}	
	char outputFName[1024];
	sprintf(outputFName,"%s/holdoutcc.txt",outputDir);
	ofstream hFile(outputFName);
	VSET& varSet=varMgr->getVariableSet();
	for(map<int,double>::iterator vIter=predictionCorr.begin();vIter!=predictionCorr.end();vIter++)
	{
		hFile << vIter->first <<"\t" << varSet[vIter->first]->getName()  << "\t" << vIter->second  << "\t" << predictionProb[vIter->first]<< endl;
	}
	hFile.close();

	sprintf(outputFName,"%s/likelihood.txt",outputDir);
	ofstream llFile(outputFName);
	for(map<int,double>::iterator tIter=holdoutLL.begin();tIter!=holdoutLL.end();tIter++)
	{
		llFile << trainLL[tIter->first] <<"\t" << tIter->second << "\t" << holdout_unpenLL[tIter->first] << endl;
	}
	llFile.close();
	return 0;
}


int
MotifRegressor::predictTestData(map<string,map<string,int>*>& dataSet)
{
	char outputFName[1024];
	map<string,int>& motifnameIDMap=motifMgr->getMotifNameIDMap();
	sprintf(outputFName,"%s/testdataprofile.txt",outputDir);
	ofstream oFile(outputFName);
	for(map<string,map<string,int>*>::iterator dIter=dataSet.begin();dIter!=dataSet.end();dIter++)
	{
		map<string,int>* motifAssignment=dIter->second;
		map<int,double> motifProf;
		for(map<string,int>::iterator aIter=motifAssignment->begin();aIter!=motifAssignment->end();aIter++)
		{
			if(motifnameIDMap.find(aIter->first)==motifnameIDMap.end())
			{
				continue;
			}
			int mId=motifnameIDMap[aIter->first];
			if(predsWithVars.find(mId)==predsWithVars.end())
			{
				continue;
			}
			motifProf[mId]=aIter->second;
		}
		map<int,double> mWts;
		int eId=0;
		double maxProb=0;
		double sum=0;
		for(map<int,ExpertL*>::iterator eIter=expertSet.begin();eIter!=expertSet.end();eIter++)
		{
			double mwt=eIter->second->getMixtureWeight(&motifProf);
			mWts[eIter->first]=exp(mwt);
			if(mwt>maxProb)
			{
				maxProb=mwt;
				eId=eIter->first;
			}
			sum=sum+exp(mwt);
		}
		ExpertL* e=expertSet[eId];
		Matrix* mean=e->getMean();
		oFile <<dIter->first;
		int dim=evMgr->getNumberOfEvidences();
		for(int i=0;i<dim;i++)
		{
			oFile << "\t" << mean->getValue(0,i);
		}
		oFile << "\t" << exp(maxProb)/sum << endl;
		motifProf.clear();
		mWts.clear();
	}
	VSET& varSet=varMgr->getVariableSet();
	oFile.close();
	return 0;
}

int
MotifRegressor::predictTestExpression()
{
	Distance d;
	VSET& varSet=varMgr->getVariableSet();
	for(map<int,int>::iterator vIter=varsWithPredictors_Test.begin();vIter!=varsWithPredictors_Test.end();vIter++)
	{
		map<int,double>* motifProf=motifMgr->getMotifProfile(varSet[vIter->first]->getName().c_str());
		map<int,double> mWts;
		int eId=0;
		double maxProb=0;
		double sum=0;
		for(map<int,ExpertL*>::iterator eIter=expertSet.begin();eIter!=expertSet.end();eIter++)
		{
			//double mwt=eIter->second->getMixtureWeight(motifProf);
			double mwt=1;
			mWts[eIter->first]=exp(mwt);
			if(mwt>maxProb)
			{
				maxProb=mwt;
				eId=eIter->first;
			}
			sum=sum+exp(mwt);
		}
		ExpertL* e=expertSet[eId];
		Matrix* mean=e->getMean();
		vector<double> predictedExp;
		vector<double> trueExp;
		int dim=evMgr->getNumberOfEvidences();
		for(int i=0;i<dim;i++)
		{
			predictedExp.push_back(mean->getValue(0,i));
			EMAP* evidMap=evMgr->getEvidenceAt(i);
			Evidence* evid=(*evidMap)[vIter->first];
			trueExp.push_back(evid->getEvidVal());
		}
		double cc=d.computeCC(predictedExp,trueExp);
		predictionCorr[vIter->first]=cc;
		double pval=exp(maxProb)/sum;
		predictionProb[vIter->first]=pval;
	}
	return 0;
}


double
MotifRegressor::getTestDataLikelihood(double& unpenalized)
{
	VSET& varSet=varMgr->getVariableSet();
	double ll=0;
	int dim=evMgr->getNumberOfEvidences();
	Matrix* dataPt=new Matrix(1,dim);
	for(map<int,int>::iterator vIter=varsWithPredictors_Test.begin();vIter!=varsWithPredictors_Test.end();vIter++)
	{
		map<int,double> mixOutProbs;
		map<int,double> mixWts;
		for(int i=0;i<dim;i++)
		{
			EMAP* evidenceSet=evMgr->getEvidenceAt(i);
			Evidence* evid=(*evidenceSet)[vIter->first];
			double aval=evid->getEvidVal();
			dataPt->setValue(aval,0,i);
		}
		for(map<int,ExpertL*>::iterator eIter=expertSet.begin();eIter!=expertSet.end();eIter++)
		{
			ExpertL* e=eIter->second;
			double pdf=e->getOutputPDF_Fast(dataPt);
			//double pdf=e->getOutputPDF(dataPt);
			double pp=e->getPrior();
			double lpval=pdf+log(pp);
			mixOutProbs[eIter->first]=lpval;
			mixWts[eIter->first]=lpval;
		}
		normalizeWeights(mixOutProbs);
		/*double sum=0;	
		for(map<int,double>::iterator eIter=mixOutProbs.begin();eIter!=mixOutProbs.end();eIter++)
		{
			sum=sum+(exp(eIter->second));	
		}
		for(map<int,double>::iterator eIter=mixOutProbs.begin();eIter!=mixOutProbs.end();eIter++)
		{
			double val=exp(eIter->second)/sum;
			eIter->second=val;
		}*/
		double expdataLL=0;
		for(map<int,double>::iterator dIter=mixWts.begin();dIter!=mixWts.end();dIter++)
		{
			
			expdataLL=expdataLL+(dIter->second*mixOutProbs[dIter->first]);
		}
		ll=ll+expdataLL;
	}
	delete dataPt;
	int paramCnt=expertSet.size()*2*dim;
	//int paramCnt=expertSet.size()*2;
	double modelComplexity=(paramCnt/2)*log(varsWithPredictors_Test.size());
	unpenalized=ll;
	ll=ll-modelComplexity;
	return ll;
}

int
MotifRegressor::setUntransformedData(const char* dataFName)
{
	strcpy(untransformedFName,dataFName);
	return 0;
}


int
MotifRegressor::showMoE()
{
	char aFName[1024];
	sprintf(aFName,"%s/moemodel.txt",outputDir);
	ofstream oFile(aFName);
	map<int,string>& motifIDInterpret=motifMgr->getMotifIDInterpretation();
	VSET& varSet=varMgr->getVariableSet();
	//assignGenesToExperts();
	assignGenesToExperts_Exclusive();
	for(map<int,ExpertL*>::iterator aIter=expertSet.begin();aIter!=expertSet.end();aIter++)
	{
		aIter->second->sortFeatures();
		//oFile << "Expert" << aIter->first;
		cout << "Expert" << aIter->first;
		vector<int>& sortedFeatures=aIter->second->getSortedFeatures();
		INTDBLMAP& fwts=aIter->second->getFeatureWeights();
		for(int i=0;i<10;i++)
		{
			int mid=sortedFeatures[i];
			if(motifIDInterpret.find(mid)==motifIDInterpret.end())
			{
				cout <<"No motif id " << mid  << " for expert" << aIter->first<< endl;
				continue;
			}
			//oFile <<"\t" << motifIDInterpret[mid].c_str()<<":("<< fwts[mid] <<")";
			cout <<"\t" << motifIDInterpret[mid].c_str()<<":("<< fwts[mid] <<")";
		}
		
		cout << endl;
		//oFile << "Genes for Expert" << aIter->first << endl;
		map<int,int>& geneSet=aIter->second->getGeneSet();
		for(map<int,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			oFile << varSet[gIter->first]->getName() <<"\t" << aIter->first << endl;
		}
	}
	oFile.close();
	return 0;
}

int
MotifRegressor::showMoEParameters(const char* aFName)
{
	ofstream pFile(aFName);
	for(map<int,ExpertL*>::iterator aIter=expertSet.begin();aIter!=expertSet.end();aIter++)
	{
		ExpertL* exp=aIter->second;
		pFile <<"Expert"<<aIter->first<<"_GaussianMean";
		Matrix* m=exp->getMean();
		Matrix* cov=exp->getCovariance();
		int dim=evMgr->getNumberOfEvidences();
		for(int i=0;i<dim;i++)
		{
			pFile << "\t" << m->getValue(0,i);
		}
		pFile << endl;
		pFile <<"Expert"<<aIter->first<<"_GaussianCov";
		for(int i=0;i<dim;i++)
		{
			for(int j=0;j<dim;j++)
			{
				pFile << "\t" << cov->getValue(i,j);
			}
			pFile << endl;
		}
	}
	
	pFile.close();
	return 0;
}


int 
MotifRegressor::showGenatomyModule()
{
	//Expression file
	char expressionFName[1024];
	sprintf(expressionFName,"%s/expr.tab",outputDir);
	ofstream dataFile(expressionFName);

	char modulenetFName[1024];
	sprintf(modulenetFName,"%s/rtee.xml",outputDir);
	ofstream oFile(modulenetFName);
	VSET& varSet=varMgr->getVariableSet();
	map<int,string>& motifIDInterpret=motifMgr->getMotifIDInterpretation();
	oFile <<"<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
	oFile <<"<root>" << endl;

	for(map<int,ExpertL*>::iterator aIter=expertSet.begin();aIter!=expertSet.end();aIter++)
	{
		ExpertL* exp=aIter->second;
		oFile <<"<Module Name=\"Expert"<< aIter->first <<"\">";
		oFile <<"<Samples>" << endl;
		map<int,int>& geneSet=exp->getGeneSet();
		for(map<int,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			if(gIter!=geneSet.begin())
			{
				oFile <<"\t";
			}
			oFile << varSet[gIter->first]->getName();
		}
		oFile<< endl<< "</Samples>" << endl;
		vector<int>& sortedFeatures=aIter->second->getSortedFeatures();
		INTDBLMAP& fwts=aIter->second->getFeatureWeights();
		oFile <<"<Set>" <<endl;
		int dim=evMgr->getNumberOfEvidences();
		for(int d=0;d<dim;d++)
		{
			if(d>0)
			{
				oFile << "\t";
			}
			oFile <<"Expr" <<d;
		}
		oFile << endl <<"</Set>" << endl;
		oFile <<"<Regulators>" <<endl;
		for(int i=0;i<10;i++)
		{
			int mid=sortedFeatures[i];
			if(motifIDInterpret.find(mid)==motifIDInterpret.end())
			{
				cout <<"No motif id " << mid  << " for expert" << aIter->first<< endl;
				continue;
			}
			if(i>0)
			{
				oFile <<"\t";
			}
			oFile << motifIDInterpret[mid].c_str();
		}
		oFile << endl<<"</Regulators>" << endl;
		oFile << "<coefficients>" << endl;
		for(int i=0;i<10;i++)
		{
			int mid=sortedFeatures[i];
			if(i>0)
			{
				oFile <<"\t";
			}
			oFile << fwts[mid];
		}
		oFile<< endl <<"</coefficients>" << endl;
		oFile <<"</Module>"<< endl;
	}
	oFile <<"</root>" << endl;
	oFile.close();

	int did=0;
	dataFile <<"Name";
	for(map<int,ExpertL*>::iterator aIter=expertSet.begin();aIter!=expertSet.end();aIter++)
	{
		ExpertL* exp=aIter->second;
		map<int,int>& geneSet=exp->getGeneSet();
		if(aIter!=expertSet.begin())
		{
			dataFile <<"\tDummy"<< did; 
		}
		for(map<int,int>::iterator gIter=geneSet.begin();gIter!=geneSet.end();gIter++)
		{
			dataFile <<"\t"<< varSet[gIter->first]->getName();
		}
		did++;
	}
	dataFile << endl;
	int dim=evMgr->getNumberOfEvidences();
	for(int d=0;d<dim;d++)
	{
		EMAP* evidSet=evMgr->getEvidenceAt(d);
		dataFile <<"Expr"<<d;
		
		for(map<int,ExpertL*>::iterator aIter=expertSet.begin();aIter!=expertSet.end();aIter++)
		{
			ExpertL* exp=aIter->second;
			map<int,int>& geneset=exp->getGeneSet();
			if(aIter!=expertSet.begin())
			{
				dataFile <<"\t-100";
			}
			for(map<int,int>::iterator gIter=geneset.begin();gIter!=geneset.end();gIter++)
			{
				Evidence* evid=(*evidSet)[gIter->first];
				dataFile << "\t"<< evid->getEvidVal();
			}
		}
		dataFile << endl;
	}
	map<string,INTDBLMAP*>& motifProfiles=motifMgr->getMotifProfileSet();
	//Need to include only those motifs that have a gene with expression value
	for(map<int,int>::iterator pIter=predsWithVars.begin();pIter!=predsWithVars.end();pIter++)
	{
		int mid=pIter->first;
		dataFile << motifIDInterpret[mid].c_str();
		for(map<int,ExpertL*>::iterator aIter=expertSet.begin();aIter!=expertSet.end();aIter++)
		{
			ExpertL* exp=aIter->second;
			map<int,int>& geneset=exp->getGeneSet();
			if(aIter!=expertSet.begin())
			{
				dataFile <<"\t-100";
			}
			for(map<int,int>::iterator gIter=geneset.begin();gIter!=geneset.end();gIter++)
			{
				INTDBLMAP* mProf=motifProfiles[varSet[gIter->first]->getName()];
				if(mProf->find(pIter->first)==mProf->end())
				{
					dataFile << "\t0";
				}
				else
				{
					dataFile<<"\t"<<(*mProf)[pIter->first];
				}	
			}
		}
		dataFile << endl;
	}
	
	dataFile.close();
	return 0;
}

int
MotifRegressor::showClusterAssignment()
{
	assignGenesToExperts();
	//assignGenesToExperts_Exclusive();
	char assignFName[1024];
	sprintf(assignFName,"%s/clusterassign.txt",outputDir);
	ofstream cFile(assignFName);
	VSET& varSet=varMgr->getVariableSet();
	for(map<int,ExpertL*>::iterator aIter=expertSet.begin();aIter!=expertSet.end();aIter++)
	{
		ExpertL* exp=aIter->second;
		map<int,int>& geneset=exp->getGeneSet();
		for(map<int,int>::iterator gIter=geneset.begin();gIter!=geneset.end();gIter++)
		{
			cFile << varSet[gIter->first]->getName() << "\t" << aIter->first << endl;
		}
	}
	cFile.close();
	return 0;
}

int
MotifRegressor::dispTFsPerCluster()
{
	char aFName[1024];
	sprintf(aFName,"%s/tfsforcluster.txt",outputDir);
	ofstream tfFile(aFName);
	sprintf(aFName,"%s/motifsforcluster.txt",outputDir);
	ofstream mFile(aFName);
	sprintf(aFName,"%s/clusterassign.txt",outputDir);
	ofstream oFile(aFName);
	showGammas(outputDir);
	sprintf(aFName,"%s/exprtab.txt",outputDir);
	ofstream eFile(aFName);
	sprintf(aFName,"%s/motifprofile.txt",outputDir);
	ofstream mpFile(aFName);
	sprintf(aFName,"%s/fwt_enr.txt",outputDir);
	ofstream fFile(aFName);
	sprintf(aFName,"%s/moemodel_params.txt",outputDir);
	showMoEParameters(aFName);
	//assignGenesToExperts();
	assignGenesToExperts_Exclusive();
	//assignGenesToExperts_Exclusive_ML();
	map<int,string>& motifIDMap=motifMgr->getMotifIDMap();
	map<string,STRINTMAP*>& motifTFMap=motifMgr->getMotifTFMap();
	VSET& varSet=varMgr->getVariableSet();
	int evidCnt=evMgr->getNumberOfEvidences();
	eFile <<"Gene";
	for(int e=0;e<evidCnt;e++)
	{
		eFile<<"\tExpr"<<e;
	}
	eFile <<endl;
	map<string,int> driverMotifs;
	vector<string> driverMotifsOrder;
	for(map<int,ExpertL*>::iterator aIter=expertSet.begin();aIter!=expertSet.end();aIter++)
	{
		ExpertL* e=aIter->second;
		INTDBLMAP featEnrich;
		getFeatureEnrichment(e,featEnrich);
		//e->sortFeatures();
		e->sortFeatures_Enrichment(featEnrich);
		vector<int>& sortedFeatures=e->getSortedFeatures();
		INTDBLMAP& fwts=e->getFeatureWeights();
		int hitTF=0;
		mFile <<aIter->first;
		for(int i=0;i<sortedFeatures.size();i++)
		{
			int fid=sortedFeatures[i];
			if( (strstr(motifIDMap[fid].c_str(),"MSN")==NULL) 
			&& (strstr(motifIDMap[fid].c_str(),"SKO")==NULL)
			&&(strstr(motifIDMap[fid].c_str(),"HOT")==NULL)
			)
			{
				if(featEnrich[fid]>0.01)
				{
					continue;
				}
			}
			const string& motifname=motifIDMap[fid];
			mFile << "\t" << motifname;
		
			if(motifTFMap.find(motifname)==motifTFMap.end())
			{
				continue;
			}
			if(hitTF==0)
			{
				tfFile <<aIter->first;
				hitTF++;
			}
			STRINTMAP* tflist=motifTFMap[motifname];
			for(STRINTMAP_ITER tIter=tflist->begin();tIter!=tflist->end();tIter++)
			{
				tfFile <<"\t" <<tIter->first;
			}
			if(driverMotifs.find(motifname)==driverMotifs.end())
			{
				driverMotifs[motifname]=fid;
				driverMotifsOrder.push_back(motifname);
			}
			fFile <<motifname<<"\t" << featEnrich[fid] << "\t" << fwts[fid] << endl;
		}
		mFile << endl;
		featEnrich.clear();
		if(hitTF)
		{
			tfFile <<endl;
		}
		map<int,int>& geneset=e->getGeneSet();
		for(map<int,int>::iterator gIter=geneset.begin();gIter!=geneset.end();gIter++)
		{
			oFile << varSet[gIter->first]->getName() << "\t" << aIter->first << endl;
			eFile << varSet[gIter->first]->getName();
			for(int s=0;s<evidCnt;s++)
			{
				EMAP* evidSet=evMgr->getEvidenceAt(s);
				Evidence* evid=(*evidSet)[gIter->first];
				eFile <<"\t" << evid->getEvidVal();
			}
			eFile << endl;
		}
		eFile <<"Dummy"<<aIter->first;
		for(int s=0;s<evidCnt;s++)
		{
			EMAP* evidSet=evMgr->getEvidenceAt(s);
			eFile <<"\t-100";
		}
		eFile << endl;
	}
	mpFile <<"Gene";
	for(int m=0;m<driverMotifsOrder.size();m++)
	{
		mpFile <<"\t"<<driverMotifsOrder[m]<<"-TF";
	}
	mpFile << endl;
	for(map<int,ExpertL*>::iterator aIter=expertSet.begin();aIter!=expertSet.end();aIter++)
	{
		ExpertL* e=aIter->second;
		map<int,int>& geneset=e->getGeneSet();
		for(map<int,int>::iterator gIter=geneset.begin();gIter!=geneset.end();gIter++)
		{
			INTDBLMAP* motifProf=motifMgr->getMotifProfile(varSet[gIter->first]->getName().c_str());
			if(motifProf==NULL)
			{
				continue;
			}
			mpFile << varSet[gIter->first]->getName().c_str();
			for(int m=0;m<driverMotifsOrder.size();m++)
			{
				int mid=driverMotifs[driverMotifsOrder[m]];
				double hit=0;
				if(motifProf->find(mid)!=motifProf->end())
				{
					//hit=(*motifProf)[mid];
					hit=1;
				}
				mpFile <<"\t" << hit;
			}
			mpFile << endl;
		}
		mpFile<<"Dummy"<<aIter->first;
		for(int m=0;m<driverMotifsOrder.size();m++)
		{
			mpFile<<"\t0";
		}
		mpFile << endl;
	}
	tfFile.close();
	oFile.close();
	mFile.close();	
	mpFile.close();	
	eFile.close();
	fFile.close();

	if(strlen(untransformedFName)>0)
	{
		showUntransformedFName();
	}

	return 0;
}

int
MotifRegressor::showUntransformedFName()
{
	GeneExpManager expmgr;
	if(strlen(untransformedFName)==0)
	{
		return 0;
	}
	expmgr.readExpression(untransformedFName);
	char expFile[1024];
	sprintf(expFile,"%s/exprmat_orig.txt",outputDir);
	ofstream oFile(expFile);
	VSET& varSet=varMgr->getVariableSet();
	int samplesize=0;
	vector<string>& headerNames=expmgr.getExpHeaders();
	oFile <<"Gene";
	for(int t=0;t<headerNames.size();t++)
	{
		oFile <<"\t" << headerNames[t];
	}
	oFile << endl;
	for(map<int,ExpertL*>::iterator aIter=expertSet.begin();aIter!=expertSet.end();aIter++)
	{
		ExpertL* e=aIter->second;
		map<int,int>& geneset=e->getGeneSet();
		for(map<int,int>::iterator gIter=geneset.begin();gIter!=geneset.end();gIter++)
		{
			oFile << varSet[gIter->first]->getName();
			vector<double>* expvect=expmgr.getExp(varSet[gIter->first]->getName());
			for(int s=0;s<expvect->size();s++)
			{
				oFile <<"\t" << (*expvect)[s];
			}
			oFile << endl;
			if(samplesize==0)
			{
				samplesize=expvect->size();
			}
		}
		oFile <<"Dummy"<<aIter->first;
		for(int s=0;s<samplesize;s++)
		{
			oFile <<"\t-100";
		}
		oFile << endl;
	}
	oFile.close();
	return 0;
}


int
MotifRegressor::dispTFsPerCluster(int fold)
{
	char dirname[1024];
	sprintf(dirname,"mkdir %s/fold%d",outputDir,fold);
	system(dirname);
	char aFName[1024];
	sprintf(aFName,"%s/fold%d/tfsforcluster.txt",outputDir,fold);
	ofstream tfFile(aFName);
	sprintf(aFName,"%s/fold%d/motifsforcluster.txt",outputDir,fold);
	ofstream mFile(aFName);
	sprintf(aFName,"%s/fold%d/clusterassign.txt",outputDir,fold);
	ofstream oFile(aFName);
	sprintf(aFName,"%s/fold%d/exprtab.txt",outputDir,fold);
	ofstream eFile(aFName);
	sprintf(aFName,"%s/fold%d/motifprofile.txt",outputDir,fold);
	ofstream mpFile(aFName);
	sprintf(aFName,"%s/fold%d/fwt_enr.txt",outputDir,fold);
	ofstream fFile(aFName);
	//assignGenesToExperts();
	assignGenesToExperts_Exclusive();
	//assignGenesToExperts_Exclusive_ML();
	map<int,string>& motifIDMap=motifMgr->getMotifIDMap();
	map<string,STRINTMAP*>& motifTFMap=motifMgr->getMotifTFMap();
	VSET& varSet=varMgr->getVariableSet();
	int evidCnt=evMgr->getNumberOfEvidences();
	eFile <<"Gene";
	for(int e=0;e<evidCnt;e++)
	{
		eFile<<"\tExpr"<<e;
	}
	eFile <<endl;
	map<string,int> driverMotifs;
	vector<string> driverMotifsOrder;
	for(map<int,ExpertL*>::iterator aIter=expertSet.begin();aIter!=expertSet.end();aIter++)
	{
		ExpertL* e=aIter->second;
		INTDBLMAP featEnrich;
		getFeatureEnrichment(e,featEnrich);
		//e->sortFeatures();
		e->sortFeatures_Enrichment(featEnrich);
		vector<int>& sortedFeatures=e->getSortedFeatures();
		INTDBLMAP& fwts=e->getFeatureWeights();
		int hitTF=0;
		mFile <<aIter->first;
		for(int i=0;i<sortedFeatures.size();i++)
		{
			int fid=sortedFeatures[i];
			if(featEnrich[fid]>0.01)
			{
				continue;
			}
			const string& motifname=motifIDMap[fid];
			mFile << "\t" << motifname;
		
			if(motifTFMap.find(motifname)==motifTFMap.end())
			{
				continue;
			}
			if(hitTF==0)
			{
				tfFile <<aIter->first;
				hitTF++;
			}
			STRINTMAP* tflist=motifTFMap[motifname];
			for(STRINTMAP_ITER tIter=tflist->begin();tIter!=tflist->end();tIter++)
			{
				tfFile <<"\t" <<tIter->first;
			}
			if(featEnrich[fid]>0.01)
			{
				continue;
			}
			if(driverMotifs.find(motifname)==driverMotifs.end())
			{
				driverMotifs[motifname]=fid;
				driverMotifsOrder.push_back(motifname);
			}
			fFile <<motifname<<"\t" << featEnrich[fid] << "\t" << fwts[fid] << endl;
		}
		mFile << endl;
		featEnrich.clear();
		if(hitTF)
		{
			tfFile <<endl;
		}
		map<int,int>& geneset=e->getGeneSet();
		for(map<int,int>::iterator gIter=geneset.begin();gIter!=geneset.end();gIter++)
		{
			oFile << varSet[gIter->first]->getName() << "\t" << aIter->first << endl;
			eFile << varSet[gIter->first]->getName();
			for(int s=0;s<evidCnt;s++)
			{
				EMAP* evidSet=evMgr->getEvidenceAt(s);
				Evidence* evid=(*evidSet)[gIter->first];
				eFile <<"\t" << evid->getEvidVal();
			}
			eFile << endl;
		}
		eFile <<"Dummy"<<aIter->first;
		for(int s=0;s<evidCnt;s++)
		{
			EMAP* evidSet=evMgr->getEvidenceAt(s);
			eFile <<"\t-100";
		}
		eFile << endl;
	}
	mpFile <<"Gene";
	for(int m=0;m<driverMotifsOrder.size();m++)
	{
		mpFile <<"\t"<<driverMotifsOrder[m]<<"-TF";
	}
	mpFile << endl;
	for(map<int,ExpertL*>::iterator aIter=expertSet.begin();aIter!=expertSet.end();aIter++)
	{
		ExpertL* e=aIter->second;
		map<int,int>& geneset=e->getGeneSet();
		for(map<int,int>::iterator gIter=geneset.begin();gIter!=geneset.end();gIter++)
		{
			INTDBLMAP* motifProf=motifMgr->getMotifProfile(varSet[gIter->first]->getName().c_str());
			if(motifProf==NULL)
			{
				continue;
			}
			mpFile << varSet[gIter->first]->getName().c_str();
			for(int m=0;m<driverMotifsOrder.size();m++)
			{
				int mid=driverMotifs[driverMotifsOrder[m]];
				double hit=0;
				if(motifProf->find(mid)!=motifProf->end())
				{
					//hit=(*motifProf)[mid];
					hit=1;
				}
				mpFile <<"\t" << hit;
			}
			mpFile << endl;
		}
		mpFile<<"Dummy"<<aIter->first;
		for(int m=0;m<driverMotifsOrder.size();m++)
		{
			mpFile<<"\t0";
		}
		mpFile << endl;
	}
	tfFile.close();
	oFile.close();
	mFile.close();	
	mpFile.close();	
	eFile.close();
	fFile.close();
	return 0;
}

int
MotifRegressor::getFeatureEnrichment(ExpertL* e,INTDBLMAP& featEnrich)
{
	VSET& varSet=varMgr->getVariableSet();
	INTDBLMAP& featWeights=e->getFeatureWeights();
	INTINTMAP featCnts;
	map<int,int>& geneset=e->getGeneSet();
	for(map<int,int>::iterator gIter=geneset.begin();gIter!=geneset.end();gIter++)
	{
		INTDBLMAP* motifProf=motifMgr->getMotifProfile(varSet[gIter->first]->getName().c_str());
		if(motifProf==NULL)
		{
			continue;
		}
		//cout <<"Number of motifs for " << varSet[gIter->first]->getName() <<" "<< motifProf->size() << endl;
		for(INTDBLMAP_ITER dIter=motifProf->begin();dIter!=motifProf->end();dIter++)
		{
			if(featWeights.find(dIter->first)==featWeights.end())
			{
				continue;
			}
			if(featCnts.find(dIter->first)==featCnts.end())
			{
				featCnts[dIter->first]=1;
			}
			else
			{
				featCnts[dIter->first]=featCnts[dIter->first]+1;
			}
		}
	}
	HyperGeomPval hgp;
	for(INTINTMAP_ITER fIter=featCnts.begin();fIter!=featCnts.end();fIter++)
	{
		
		int n1=predsWithVars[fIter->first];
		int n2=varsWithPredictors_Train.size()-n1;
		int k=fIter->second;
		int t=geneset.size();
		double enpval=hgp.getOverRepPval(t,k,n1,n2);
		featEnrich[fIter->first]=enpval;
	}	
	return 0;
}


int
MotifRegressor::readClusterMembership()
{
	ifstream inFile(clusterFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;	
		}
		char* tok=strtok(buffer,"\t");
		string geneName;
		int clusterID;
		int tokCnt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				geneName.append(tok);
			}
			else
			{
				clusterID=atoi(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		map<string,int>* cluster=NULL;
		if(gmmClusters.find(clusterID)==gmmClusters.end())
		{
			cluster=new map<string,int>;
			gmmClusters[clusterID]=cluster;	
		}
		else
		{
			cluster=gmmClusters[clusterID];
		}
		(*cluster)[geneName]=0;
	}
	inFile.close();
	return 0;
}



int
MotifRegressor::showGammas(const char* outputDir)
{
	char outputFName[1024];
	sprintf(outputFName,"%s/gammas.txt",outputDir);
	ofstream oFile(outputFName);
	VSET& varSet=varMgr->getVariableSet();
	for(map<int,map<int,double>*>::iterator dIter=gammas.begin();dIter!=gammas.end();dIter++)
	{
		map<int,double>* gamma_d=dIter->second;
		double maxgamma=0;
		int bestexpert=-1;
		oFile << varSet[dIter->first]->getName();
		for(map<int,double>::iterator gIter=gamma_d->begin();gIter!=gamma_d->end();gIter++)
		{
			oFile << "\t" <<gIter->second;
		}
		oFile << endl;
	}
	oFile.close();
	return 0;
}

int
MotifRegressor::initExpertsKMeans()
{
	ClusterManager cMgr;
	map<int,INTDBLMAP*> kmeansData;
	map<int,int> dataIDMap;
	VSET& varSet=varMgr->getVariableSet();
	int dim=evMgr->getNumberOfEvidences();
	int dId=0;
	for(map<int,int>::iterator vIter=varsWithPredictors_Train.begin();vIter!=varsWithPredictors_Train.end();vIter++)
	{
		map<int,double>* dataPt=new map<int,double>;
		for(int i=0;i<dim;i++)
		{
			EMAP* evidenceSet=evMgr->getEvidenceAt(i);
			Evidence* evid=(*evidenceSet)[vIter->first];
			double aval=evid->getEvidVal();
			(*dataPt)[i]=aval;
		}
		dataIDMap[dId]=vIter->first;
		kmeansData[dId]=dataPt;
		dId++;
	}
	cMgr.setVariableManager(varMgr);
	cMgr.setEvidenceManager(evMgr);
	cMgr.setClusterCnt(expertCnt);
	cMgr.setOutputDir(outputDir);
	vector<INTINTMAP*>& variableSubsets=cMgr.getClusters(kmeansData);
	vector<int> randIndex;
	for(int c=0;c<variableSubsets.size();c++)
	{
		randIndex.clear();
		INTINTMAP* vSet=variableSubsets[c];
		for(INTINTMAP_ITER vIter=vSet->begin();vIter!=vSet->end();vIter++)
		{
			int tempId=vIter->first;
			int vId=dataIDMap[tempId];
			randIndex.push_back(vId);
		}
		ExpertL* expert=new ExpertL;
		//Intialize the potentials with these parameters
		estimateMeanCov(randIndex,0,randIndex.size(),expert);
		//estimateGlobalCov(expert);
		//estimateMean(randIndex,0,randIndex.size(),expert);
		setInitFeatureWeights(expert);
		expertSet[c]=expert;
		//All the data from startInd to endInd have gammas set to 1 for this mixture component
		for(int i=0;i<randIndex.size();i++)
		{
			map<int,double>* gamma_i=new map<int,double>;
			int dataId=randIndex[i];
			gammas[dataId]=gamma_i;
			for(int f=0;f<expertCnt;f++)
			{
				if(c==f)
				{
					(*gamma_i)[f]=1;
				}
				else
				{
					(*gamma_i)[f]=0;
				}
			}
		}
	}
	return 0;
}

int 
MotifRegressor::initExperts()
{
	//partition data into expertCnt number of partitions
	vector<int> randIndex;
	populateRandIntegers(randIndex,varsWithPredictors_Train.size());
	cout <<"Rand integers "<< randIndex.size() << " matindex " << matIdvIdMap.size () << endl;
	//For each partition estimate the mean and covariance
	int testSetSize=varsWithPredictors_Train.size()/expertCnt;
	for(int e=0;e<expertCnt;e++)
	{
		int startInd=e*testSetSize;
		int endInd=(e+1)*testSetSize;
		if(e==expertCnt-1)
		{
			endInd=varsWithPredictors_Train.size();
		}
		ExpertL* expert=new ExpertL;
		//Intialize the potentials with these parameters
		estimateMeanCov(randIndex,startInd,endInd,expert);
		//estimateGlobalCov(expert);
		//estimateMean(randIndex,startInd,endInd,expert);
		setInitFeatureWeights(expert);
		expertSet[e]=expert;
		//All the data from startInd to endInd have gammas set to 1 for this mixture component
		for(int i=startInd;i<endInd;i++)
		{
			map<int,double>* gamma_i=new map<int,double>;
			int dataId=randIndex[i];
			gammas[dataId]=gamma_i;
			for(int f=0;f<expertCnt;f++)
			{
				if(e==f)
				{
					(*gamma_i)[f]=1;
				}
				else
				{
					(*gamma_i)[f]=0;
				}
			}
		}
	}
	randIndex.clear();
	return 0;
}


int 
MotifRegressor::initExpertsGMM()
{
	readClusterMembership();
	int actualExpertCnt=0;
	for(map<int,map<string,int>*>::iterator gIter=gmmClusters.begin();gIter!=gmmClusters.end();gIter++)
	{
		ExpertL* expert=new ExpertL;
		map<string,int>* cluster=gIter->second;
		int addedMembers=0;
		for(map<string,int>::iterator cIter=cluster->begin();cIter!=cluster->end();cIter++)
		{
			if(varsWithPredictors_Train_NameID.find(cIter->first)==varsWithPredictors_Train_NameID.end())
			{
				continue;
			}
			addedMembers++;
			map<int,double>* gamma_i=new map<int,double>;
			int dataId=varsWithPredictors_Train_NameID[cIter->first];
			gammas[dataId]=gamma_i;
			expert->assignGeneToExpert(dataId);
		}
		if(addedMembers==0)
		{
			delete expert;
		}
		else
		{
			setInitFeatureWeights(expert);
			expertSet[actualExpertCnt]=expert;
			actualExpertCnt++;
		}
	}
	expertCnt=actualExpertCnt;
	for(map<int,ExpertL*>::iterator eIter=expertSet.begin();eIter!=expertSet.end();eIter++)
	{
		map<int,int>& memberGenes=eIter->second->getGeneSet();
		for(map<int,int>::iterator mIter=memberGenes.begin();mIter!=memberGenes.end();mIter++)
		{
			map<int,double>* gamma_i=gammas[mIter->first];
			for(map<int,ExpertL*>::iterator fIter=expertSet.begin();fIter!=expertSet.end();fIter++)
			{
				if(fIter->first==eIter->first)
				{
					(*gamma_i)[fIter->first]=1;
				}
				else
				{
					(*gamma_i)[fIter->first]=0;
				}
			}
		}
	}
	int fId=0;
	for(INTINTMAP_ITER fIter=predsWithVars.begin();fIter!=predsWithVars.end();fIter++)
	{
		featMatidIDMap[fId]=fIter->first;
		featIDMatidMap[fIter->first]=fId;
		fId++;
	}
	for(int e=0;e<expertCnt;e++)
	{
		ExpertL* expert=expertSet[e];
		int dim=evMgr->getNumberOfEvidences();
		Matrix* mean=new Matrix(1,dim);
		expert->setMean(mean);
		//Matrix* covariance=new Matrix(dim,dim);
		//expert->setCovariance(covariance);
		estimateMeanCov(e,expert);
	}
	return 0;
}


int 
MotifRegressor::initVarsWithPredictors()
{
	VSET& varSet=varMgr->getVariableSet();
	map<string,INTDBLMAP*>& motifProfiles=motifMgr->getMotifProfileSet();
	int varsWithMotifs=0;
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		if(motifProfiles.find(vIter->second->getName())==motifProfiles.end())
		{
		//	continue;
		}
		else
		{
			varsWithMotifs++;
		}
		varsWithPredictors[vIter->first]=0;
	}
	cout <<"Genes with motifs: " << varsWithMotifs << " Total genes: " << varsWithPredictors.size() << endl;
	return 0;
}

int
MotifRegressor::initPredictorsWithVars()
{
	map<int,string>& motifIDNameMap=motifMgr->getMotifIDMap();
	VSET& varSet=varMgr->getVariableSet();
	map<string,INTDBLMAP*>& motifProfiles=motifMgr->getMotifProfileSet();
	//Need to include only those motifs that have a gene with expression value
	int maxcnt=0;
	int mincnt=1000;
	for(INTINTMAP_ITER vIter=varsWithPredictors_Train.begin();vIter!=varsWithPredictors_Train.end();vIter++)
	{
		Variable* responseRV=varSet[vIter->first];
		if(motifProfiles.find(responseRV->getName())==motifProfiles.end())
		{
			continue;
		}
		INTDBLMAP* profile=motifProfiles[responseRV->getName()];
		cout <<"Number of motifs: "<< responseRV->getName() << "\t" << profile->size() << endl;
		for(INTDBLMAP_ITER dIter=profile->begin();dIter!=profile->end();dIter++)
		{
			if(predsWithVars.find(dIter->first)==predsWithVars.end())
			{
				predsWithVars[dIter->first]=1;
			}
			else
			{
				predsWithVars[dIter->first]=predsWithVars[dIter->first]+1;
			}
		}
	}
	for(INTINTMAP_ITER pIter=predsWithVars.begin();pIter!=predsWithVars.end();pIter++)
	{
		if(maxcnt<pIter->second)
		{
			maxcnt=pIter->second;
		}
		if(mincnt>pIter->second)
		{
			mincnt=pIter->second;
		}
	}
	cout <<"Found instances for " << predsWithVars.size() << " motifs out of total "<< motifIDNameMap.size() << endl;
	cout <<"MaxCnt: "<< maxcnt << " MinCnt: "<< mincnt << endl;
	return 0;
}

int
MotifRegressor::split(int fId)
{
	VSET& varSet=varMgr->getVariableSet();
	varsWithPredictors_Train.clear();
	varsWithPredictors_Test.clear();
	int setSize=(int) (varsWithPredictors.size()/foldCnt);
	int startInd=fId*setSize;
	int endInd=(fId+1)*setSize;
	if(fId==foldCnt-1)
	{
		endInd=varsWithPredictors.size()-1;
	}	
	int dId=0;
	int matid=0;
	for(map<int,int>::iterator vIter=varsWithPredictors.begin();vIter!=varsWithPredictors.end();vIter++)
	{
		if((dId>=startInd && dId<endInd) && (foldCnt>1))
		{
			varsWithPredictors_Test[vIter->first]=vIter->second;
		}	
		else
		{
			varsWithPredictors_Train[vIter->first]=vIter->second;
			Variable* var=varSet[vIter->first];
			varsWithPredictors_Train_NameID[var->getName()]=vIter->first;
			matIdvIdMap[matid]=vIter->first;
			matid++;
		}
		dId++;
	}	
	return 0;
}

int
MotifRegressor::cleanUp()
{
	predsWithVars.clear();
	featIDMatidMap.clear();
	featMatidIDMap.clear();
	matIdvIdMap.clear();
	for(map<int,Matrix*>::iterator eIter=exprProfileSet.begin();eIter!=exprProfileSet.end();eIter++)
	{
		delete eIter->second;
	}
	exprProfileSet.clear();
	for(map<int,map<int,double>*>::iterator gIter=gammas.begin();gIter!=gammas.end();gIter++)
	{
		gIter->second->clear();
		delete gIter->second;
	}
	gammas.clear();
	for(map<int,ExpertL*>::iterator eIter=expertSet.begin();eIter!=expertSet.end();eIter++)
	{
		delete eIter->second;
	}
	expertSet.clear();
	bfgs.reset();
	return 0;
}


int 
MotifRegressor::expectationStep()
{
	VSET& varSet=varMgr->getVariableSet();
	int shownVars=0;
	map<int,int> reallySmallVals;
	for(map<int,int>::iterator vIter=varsWithPredictors_Train.begin();vIter!=varsWithPredictors_Train.end();vIter++)
	{
		Matrix* exprProf=exprProfileSet[vIter->first];
		map<int,double>* motifProf=motifMgr->getMotifProfile(varSet[vIter->first]->getName().c_str());
		map<int,double> mixWts;
		map<int,double> mixOutProbs;
		double wtsum=0;
		Variable* v=varSet[vIter->first];
		for(map<int,ExpertL*>::iterator eIter=expertSet.begin();eIter!=expertSet.end();eIter++)
		{
			ExpertL* e=eIter->second;
			//double wt_unnorm=e->getMixtureWeight(motifProf);
			double wt_unnorm=1;
			double exp_wt_unnorm=exp(wt_unnorm);
			mixWts[eIter->first]=exp_wt_unnorm;
			wtsum=wtsum+exp_wt_unnorm;
			//double pdf=e->getOutputPDF_Fast(exprProf);
			if(strcmp(v->getName().c_str(),"orf19.951")==0|| strcmp(v->getName().c_str(),"orf19.94")==0)
			{
				cout <<"Stop here" << endl;
			}
			//double pdf=e->getOutputPDF(exprProf);
			double pdf=e->getOutputPDF_Fast(exprProf);
			double pp=e->getPrior();
			//double pp=1;
			//mixOutProbs[eIter->first]=pdf*pp;
			mixOutProbs[eIter->first]=pdf+log(pp);
		}
		normalizeWeights(mixOutProbs);
		/*double sum=0;	
		for(map<int,double>::iterator eIter=mixOutProbs.begin();eIter!=mixOutProbs.end();eIter++)
		{
			sum=sum+(exp(eIter->second));	
		}*/
		/*double sum=0;	
		for(map<int,double>::iterator eIter=mixOutProbs.begin();eIter!=mixOutProbs.end();eIter++)
		{
			//double p1=eIter->second;
			double p1=eIter->second;
			//double p2=mixWts[eIter->first]/wtsum;
			double p2=1;
			if(p1==0 || p2==0)
			{
				cout << "p1=" << p1  << " p2=" << p2<< endl;
			}
			mixWts[eIter->first]=p2*p1;
			sum=sum+(p2*p1);
		}*/
		map<int,double>* gamma_i=gammas[vIter->first];
		//for(map<int,double>::iterator eIter=mixWts.begin();eIter!=mixWts.end();eIter++)
		for(map<int,double>::iterator eIter=mixOutProbs.begin();eIter!=mixOutProbs.end();eIter++)
		{
			//double wt_norm=exp(eIter->second)/sum;
			double wt_norm=eIter->second;
			(*gamma_i)[eIter->first]=wt_norm;
			if(isnan(wt_norm))
			{
				//cout << "Gamma nan data " << vIter->first << " expertid: " << eIter->first << " sum: " << sum <<" mixwt: " << eIter->second<< endl;
				cout << "Gamma nan data " << vIter->first << " expertid: " << eIter->first << " mixwt: " << eIter->second<< endl;
			}
		}
		if(strcmp(v->getName().c_str(),"orf19.951")==0|| strcmp(v->getName().c_str(),"orf19.94")==0)
		{
			cout << v->getName();
			for(map<int,double>::iterator aIter=gamma_i->begin();aIter!=gamma_i->end();aIter++)
			{
				cout <<" " << aIter->second;
			}
			cout << endl;
		}
		mixWts.clear();
		mixOutProbs.clear();
		shownVars++;
	}
	for(map<int,int>::iterator vIter=reallySmallVals.begin();vIter!=reallySmallVals.end();vIter++)
	{
		cout << vIter->first <<"\t" << vIter->second << "\tout_of\t" << varsWithPredictors_Train.size() << endl;
	}
	return 0;
}

int 
MotifRegressor::maximizationStep()
{
	for(map<int,ExpertL*>::iterator eIter=expertSet.begin();eIter!=expertSet.end();eIter++)
	{
		ExpertL* expert=eIter->second;
		estimateMeanCov(eIter->first,expert);
	}
	//bfgs.optimize();
	//updateFeatureWeights();
	return 0;
}


int
MotifRegressor::normalizeWeights(map<int,double>& wtVect)
{
	double maxExp=-1000000000;
	for(INTDBLMAP_ITER dIter=wtVect.begin();dIter!=wtVect.end();dIter++)
	{
		double tempval=dIter->second;
		if(dIter->second>maxExp)
		{
			maxExp=dIter->second;
		}
	}
	double pow=-1*maxExp;
	double total=0;
	for(INTDBLMAP_ITER dIter=wtVect.begin();dIter!=wtVect.end();dIter++)
	{
		double tempval=dIter->second;
		total=total+exp(dIter->second+pow);
	}
	total=log(total)-pow;
	//correction factor
	double newtotal=0;
	for(INTDBLMAP_ITER dIter=wtVect.begin();dIter!=wtVect.end();dIter++)
	{
		double tempval=dIter->second;
		double scWt=1e-5+exp(tempval-total);
		dIter->second=scWt;
		newtotal=newtotal+scWt;
	}
	for(INTDBLMAP_ITER dIter=wtVect.begin();dIter!=wtVect.end();dIter++)
	{
		double tempval=dIter->second;
		double scWt=tempval/newtotal;
		dIter->second=scWt;
	}
	return 0;
}

int 
MotifRegressor::populateRandIntegers(vector<int>& randInds,int evidCnt)
{
	double step=1.0/(double)evidCnt;
	map<int,int> usedInit;
	int maxind=0;
	for(int i=0;i<evidCnt;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
			if(rind>maxind)
			{
				maxind=rind;
			}
		}
		if(matIdvIdMap.find(rind)==matIdvIdMap.end())
		{
			cout <<"Did not find " << rind << " matrixidmap " << endl;
		}
		int dataind=matIdvIdMap[rind];
		usedInit[rind]=0;
		randInds.push_back(dataind);
	}
	cout << maxind << endl;
	usedInit.clear();
	return 0;
}

int 
MotifRegressor::estimateMeanCov(vector<int>& randIndex,int startInd,int endInd, ExpertL* e)
{
	int dim=evMgr->getNumberOfEvidences();
	Matrix* mean=new Matrix(1,dim);
	double prior=((double)(endInd-startInd))/((double) randIndex.size());
	for(int i=0;i<dim;i++)
	{
		//Assuming that each evidence is actually a joing assignment to RVs. What we
		//want is a mean and variance for each experiment
		EMAP* evidenceSet=evMgr->getEvidenceAt(i);
		double meanval=0;
		for(int j=startInd;j<endInd;j++)
		{
			int vId=randIndex[j];
			//cout << "Variable " << vId << "\trindex:" << j <<"\tdim:" << i<< endl;
			Evidence* evid=(*evidenceSet)[vId];
			meanval=meanval+evid->getEvidVal();
		}
		meanval=meanval/((double)(endInd-startInd));
		mean->setValue(meanval,0,i);
	}
	Matrix* covariance=new Matrix(dim,dim);
	covariance->setAllValues(0);
	for(int i=0;i<dim;i++)
	{
		double m1=mean->getValue(0,i);
		EMAP* evid1=evMgr->getEvidenceAt(i);
		//for(int j=i;j<dim;j++)
		for(int j=i;j<i+1;j++)
		{
			double m2=mean->getValue(0,j);
			EMAP* evid2=evMgr->getEvidenceAt(j);
			double cov=0;
			if(i==j)
			{
				cov=0.001;
			}
			for(int k=startInd;k<endInd;k++)
			{
				int vId=randIndex[k];
				double diff=((*evid1)[vId]->getEvidVal()-m1)*((*evid2)[vId]->getEvidVal()-m2);
				cov=cov+diff;
			}
			cov=cov/((double)(endInd-startInd)-1);
			if(i==j)
			{
				covariance->setValue(cov,i,j);
			}
			else
			{
				//covariance->setValue(cov,j,i);
				covariance->setValue(0,j,i);
			}
		}
	}
	//covariance->showMatrix();
	e->setMean(mean);
	e->setCovariance(covariance);
	e->setPrior(prior);
	return 0;
}


int 
MotifRegressor::estimateGlobalCov(ExpertL* e)
{
	int dim=evMgr->getNumberOfEvidences();
	Matrix mean(1,dim);
	for(int i=0;i<dim;i++)
	{
		//Assuming that each evidence is actually a joing assignment to RVs. What we
		//want is a mean and variance for each experiment
		EMAP* evidenceSet=evMgr->getEvidenceAt(i);
		double meanval=0;
		for(map<int,int>::iterator vIter=varsWithPredictors_Train.begin();vIter!=varsWithPredictors_Train.end();vIter++)
		{
			int vId=vIter->first;
			Evidence* evid=(*evidenceSet)[vId];
			meanval=meanval+evid->getEvidVal();
		}
		meanval=meanval/((double)(varsWithPredictors_Train.size()));
		mean.setValue(meanval,0,i);
	}
	Matrix* covariance=new Matrix(dim,dim);
	covariance->setAllValues(0);
	for(int i=0;i<dim;i++)
	{
		double m1=mean.getValue(0,i);
		EMAP* evid1=evMgr->getEvidenceAt(i);
		//for(int j=i;j<dim;j++)
		for(int j=i;j<i+1;j++)
		{
			double m2=mean.getValue(0,j);
			EMAP* evid2=evMgr->getEvidenceAt(j);
			double cov=0;
			if(i==j)
			{
				cov=0.001;
			}
			for(map<int,int>::iterator vIter=varsWithPredictors_Train.begin();vIter!=varsWithPredictors_Train.end();vIter++)
			{
				int vId=vIter->first;
				double diff=((*evid1)[vId]->getEvidVal()-m1)*((*evid2)[vId]->getEvidVal()-m2);
				cov=cov+diff;
			}
			cov=cov/((double)(varsWithPredictors_Train.size())-1);
			covariance->setValue(cov,i,j);
			covariance->setValue(cov,j,i);
		}
	}
	e->setCovariance(covariance);
	return 0;

}


int 
MotifRegressor::estimateMean(vector<int>& randIndex,int startInd,int endInd, ExpertL* e)
{
	int dim=evMgr->getNumberOfEvidences();
	Matrix* mean=new Matrix(1,dim);
	for(int i=0;i<dim;i++)
	{
		//Assuming that each evidence is actually a joing assignment to RVs. What we
		//want is a mean and variance for each experiment
		EMAP* evidenceSet=evMgr->getEvidenceAt(i);
		double meanval=0;
		for(int j=startInd;j<endInd;j++)
		{
			int vId=randIndex[j];
			Evidence* evid=(*evidenceSet)[vId];
			meanval=meanval+evid->getEvidVal();
		}
		meanval=meanval/((double)(endInd-startInd));
		mean->setValue(meanval,0,i);
	}
	e->setMean(mean);
	return 0;
}

int 
MotifRegressor::estimateMeanCov(int expertid, ExpertL* e)
{
	int dim=evMgr->getNumberOfEvidences();
	Matrix* mean=e->getMean();
	double prior=0;
	for(int i=0;i<dim;i++)
	{
		EMAP* evidenceSet=evMgr->getEvidenceAt(i);
		double meanval=0;
		double sum=0;
		for(map<int,int>::iterator vIter=varsWithPredictors_Train.begin();vIter!=varsWithPredictors_Train.end();vIter++)
		{	
			if(gammas.find(vIter->first)==gammas.end())
			{
				map<int,double>* gamma_i=new map<int,double>;
				gammas[vIter->first]=gamma_i;
				continue;
			}
			map<int,double>* gamma_i=gammas[vIter->first];
			double g_i=(*gamma_i)[expertid];
			Evidence* evid=(*evidenceSet)[vIter->first];
			meanval=meanval+(g_i*evid->getEvidVal());
			sum=sum+g_i;
		}
		meanval=meanval/sum;
		mean->setValue(meanval,0,i);
		if(i==0)
		{
			prior=(1+sum)/((double)(varsWithPredictors_Train.size()+expertSet.size()));
		}
	}
	Matrix* covariance=e->getCovariance();
	if(covariance==NULL)
	{
		covariance=new Matrix(dim,dim);
	}
	for(int i=0;i<dim;i++)
	{
		double m1=mean->getValue(0,i);
		EMAP* evid1=evMgr->getEvidenceAt(i);
		//for(int j=i;j<dim;j++)
		for(int j=i;j<i+1;j++)
		{
			double m2=mean->getValue(0,j);
			EMAP* evid2=evMgr->getEvidenceAt(j);
			double cov=0;
			if(i==j)
			{
				cov=0.001;
			}
			double sum=0;
			for(map<int,int>::iterator vIter=varsWithPredictors_Train.begin();vIter!=varsWithPredictors_Train.end();vIter++)
			{
				if(gammas.find(vIter->first)==gammas.end())
				{
					continue;
				}
				map<int,double>* gamma_i=gammas[vIter->first];
				double g_i=(*gamma_i)[expertid];
				double diff=g_i*((*evid1)[vIter->first]->getEvidVal()-m1)*((*evid2)[vIter->first]->getEvidVal()-m2);
				cov=cov+diff;
				sum=sum+g_i;
			}
			cov=cov/sum;
			covariance->setValue(cov,i,j);
			covariance->setValue(cov,j,i);
		}
	}
	if(e->getCovariance()==NULL)
	{
		e->setCovariance(covariance);
	}
	else
	{
		e->updateCovariance();
	}
	cout << "Prior: "<< expertid << ": "<<prior<<endl;
	e->setPrior(prior);
	return 0;
}


int 
MotifRegressor::setInitFeatureWeights(ExpertL* e)
{
	cout <<"Expert params"<< endl;
	for(INTINTMAP_ITER mIter=predsWithVars.begin();mIter!=predsWithVars.end();mIter++)
	{
		//double w_p=gsl_ran_flat(r,0,1);
		double w_p=gsl_ran_flat(r,-1,1);

	//	cout<< "Fid: "<< mIter->first << " wt: " << w_p << endl;
		e->setFeatureWeight(mIter->first,w_p);
	}
	return 0;
}

int 
MotifRegressor::populateBFGSData()
{
	map<int,map<int,double>*>& dataMat=bfgs.getData();
	//Not sure if the features and motif ids are ordered and indexed from 0.
	int fId=0;
	for(INTINTMAP_ITER fIter=predsWithVars.begin();fIter!=predsWithVars.end();fIter++)
	{
		featMatidIDMap[fId]=fIter->first;
		featIDMatidMap[fIter->first]=fId;
		fId++;
	}
	VSET& varSet=varMgr->getVariableSet();
	for(INTINTMAP_ITER vIter=varsWithPredictors_Train.begin();vIter!=varsWithPredictors_Train.end();vIter++)
	{
		map<int,double>* dataPt=new map<int,double>;
		for(INTINTMAP_ITER fIter=predsWithVars.begin();fIter!=predsWithVars.end();fIter++)
		{
			int fId=featIDMatidMap[fIter->first];
			double fVal=motifMgr->getMotifVal(varSet[vIter->first]->getName(),fIter->first);
			(*dataPt)[fId]=fVal;
		}
		dataMat[vIter->first]=dataPt;
		if(gammas.find(vIter->first)==gammas.end())
		{
			cout <<"No gammas for " << vIter->first << endl;
		}
	}
	gsl_vector* initx=gsl_vector_alloc(featMatidIDMap.size()*expertCnt);
	for(map<int,ExpertL*>::iterator eIter=expertSet.begin();eIter!=expertSet.end();eIter++)
	{
		ExpertL* e=eIter->second;
		INTDBLMAP& fWts=e->getFeatureWeights();
		for(INTDBLMAP_ITER fIter=fWts.begin();fIter!=fWts.end();fIter++)
		{
			double wt=fIter->second;
			int fid=featIDMatidMap[fIter->first];
			int mid=(featIDMatidMap.size()*eIter->first)+fid;
			gsl_vector_set(initx,mid,wt);
		}
	}
	bfgs.setInitVector(initx);
	return 0;
}

int 
MotifRegressor::populateBFGSGammas()
{
	map<int,map<int,double>*>& bfgs_gammas=bfgs.getGammas();
	for(map<int,map<int,double>*>::iterator aIter=gammas.begin();aIter!=gammas.end();aIter++)
	{
		bfgs_gammas[aIter->first]=aIter->second;
	}
	return 0;
}

int 
MotifRegressor::updateFeatureWeights()
{
	//cout <<"Updating weights " << endl;
	gsl_vector* allparams=bfgs.getParams();
	for(map<int,ExpertL*>::iterator eIter=expertSet.begin();eIter!=expertSet.end();eIter++)
	{
		ExpertL* e=eIter->second;
		for(map<int,int>::iterator fIter=featIDMatidMap.begin();fIter!=featIDMatidMap.end();fIter++)
		{
			int matid=fIter->second;
			int vid=(eIter->first*featIDMatidMap.size())+matid;
			double wt=gsl_vector_get(allparams,vid);
	//		cout << fIter->first<< ":"<< wt << endl;
			e->setFeatureWeight(fIter->first,wt);
		}
	}	
	return 0;
}

double 
MotifRegressor::getScore()
{
	double outputScore=0;
	for(map<int,Matrix*>::iterator dIter=exprProfileSet.begin();dIter!=exprProfileSet.end();dIter++)
	{
		Matrix* m=dIter->second;
		map<int,double>* gamma_i=gammas[dIter->first];
		for(map<int,ExpertL*>::iterator eIter=expertSet.begin();eIter!=expertSet.end();eIter++)
		{
			ExpertL* e=eIter->second;
			double g_i=(*gamma_i)[eIter->first];
			//double pval=e->getOutputPDF_Nocov(m);
			//double pval=e->getOutputPDF(m);
			double pval=e->getOutputPDF_Fast(m);
			//double pval=e->getOutputPDF(m);
			double pp=e->getPrior();
			double lpval=pval+log(pp);
			if(isnan(lpval))
			{
				//cout <<"Found nan! for " << eIter->first << " data " << dIter->first << " pval "<< pval << endl;
			}
				
			//outputScore=outputScore+(g_i*log(e->getOutputPDF_Nocov(m)));
			//outputScore=outputScore+(g_i*log(e->getOutputPDF(m)*pp));
			outputScore=outputScore+(g_i*lpval);
		}
	}
	//cout <<"Score GMM " << outputScore << " MoE " << bfgs.getOptimalFval() << endl;
	cout <<"Score GMM " << outputScore << endl;
	double overallEntropy=0;
	for(map<int,ExpertL*>::iterator eIter=expertSet.begin();eIter!=expertSet.end();eIter++)
	{
		double entropy=eIter->second->getEntropy();
		overallEntropy=overallEntropy+entropy;
	}
	cout <<"Overall entropy " << overallEntropy << endl;
	int paramCnt=expertSet.size()*2*evMgr->getNumberOfEvidences();
	double modelComplexity=(paramCnt/2)*log(exprProfileSet.size());
	double netLL=outputScore-modelComplexity;
	//double netLL=outputScore+bfgs.getOptimalFval();
	return netLL;
}

int
MotifRegressor::matrixifyData()
{
	int dim=evMgr->getNumberOfEvidences();
	for(map<int,int>::iterator vIter=varsWithPredictors_Train.begin();vIter!=varsWithPredictors_Train.end();vIter++)
	{
		Matrix* dataPt=new Matrix(1,dim);
		for(int i=0;i<dim;i++)
		{
			EMAP* evidenceSet=evMgr->getEvidenceAt(i);
			Evidence* evid=(*evidenceSet)[vIter->first];
			double aval=evid->getEvidVal();
			dataPt->setValue(aval,0,i);
		}
		exprProfileSet[vIter->first]=dataPt;
	}
	return 0;
}

int
MotifRegressor::assignGenesToExperts()
{
	for(map<int,ExpertL*>::iterator eIter=expertSet.begin();eIter!=expertSet.end();eIter++)
	{
		ExpertL* exp=eIter->second;
		exp->resetAssignedGenes();
	}
	VSET& varSet=varMgr->getVariableSet();
	for(map<int,map<int,double>*>::iterator dIter=gammas.begin();dIter!=gammas.end();dIter++)
	{
		map<int,double>* gamma_d=dIter->second;
		for(map<int,double>::iterator gIter=gamma_d->begin();gIter!=gamma_d->end();gIter++)
		{
			double pval=gIter->second;
			if(gIter->second>0.8)
			{
				ExpertL* e=expertSet[gIter->first];
				e->assignGeneToExpert(dIter->first);
			}
		}
	}
	return 0;
}


int
MotifRegressor::assignGenesToExperts_Exclusive()
{
	for(map<int,ExpertL*>::iterator eIter=expertSet.begin();eIter!=expertSet.end();eIter++)
	{
		ExpertL* exp=eIter->second;
		exp->resetAssignedGenes();
	}
	VSET& varSet=varMgr->getVariableSet();
	for(map<int,map<int,double>*>::iterator dIter=gammas.begin();dIter!=gammas.end();dIter++)
	{
		map<int,double>* gamma_d=dIter->second;
		double maxgamma=0;
		int bestexpert=-1;
		for(map<int,double>::iterator gIter=gamma_d->begin();gIter!=gamma_d->end();gIter++)
		{
			if(gIter->second>maxgamma)
			{
				maxgamma=gIter->second;
				bestexpert=gIter->first;
			}
		}
		if(bestexpert==-1)
		{
			cout <<"No gene assignment for " << varSet[dIter->first]->getName() << endl;
		}
		ExpertL* e=expertSet[bestexpert];
		e->assignGeneToExpert(dIter->first);
	}
	return 0;
}

int
MotifRegressor::assignGenesToExperts_Exclusive_ML()
{
	for(map<int,ExpertL*>::iterator eIter=expertSet.begin();eIter!=expertSet.end();eIter++)
	{
		ExpertL* exp=eIter->second;
		exp->resetAssignedGenes();
	}
	VSET& varSet=varMgr->getVariableSet();
	for(map<int,map<int,double>*>::iterator dIter=gammas.begin();dIter!=gammas.end();dIter++)
	{
		Matrix* exprProf=exprProfileSet[dIter->first];
		double maxpdf=-1000000;
		int bestexpert=-1;
		for(map<int,ExpertL*>::iterator eIter=expertSet.begin();eIter!=expertSet.end();eIter++)
		{
			ExpertL* e=eIter->second;
			double pdf=e->getOutputPDF(exprProf);
			//double pp=e->getPrior();
			//mixOutProbs[eIter->first]=pdf*pp;
			if(pdf>maxpdf)
			{
				maxpdf=pdf;
				bestexpert=eIter->first;
			}
		}
		if(bestexpert==-1)
		{
			cout <<"No gene assignment for " << varSet[dIter->first]->getName() << endl;
		}
		ExpertL* e=expertSet[bestexpert];
		e->assignGeneToExpert(dIter->first);
	}
	return 0;
}
