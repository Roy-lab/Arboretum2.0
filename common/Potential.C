#include <math.h>
#include <iostream>
#include <stack>
#include "Evidence.H"
#include "CommonTypes.H"
#include "Error.H"
#include "Variable.H"
#include "Potential.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "Rule.H"
#include "Potential.H"

Potential::Potential()
{
	dtree=NULL;
}

Potential::~Potential()
{
}

int 
Potential::setEvidenceManager(EvidenceManager* aPtr)
{
	evMgr=aPtr;
	return 0;
}

int
Potential::setMinLeafSize(int aSize)
{
	minLeafSize=aSize;
	return 0;
}

int 
Potential::setAssocVariable(Variable* var,Potential::VariableRole vRole)
{
	varSet[var->getID()]=var;
	switch(vRole)
	{
		case Potential::FACTOR:
		{
			factorVariables[var->getID()]=0;
			break;
		}
		case Potential::MARKOV_BNKT:
		{
			markovBlnktVariables[var->getID()]=0;
			break;
		}
	}
	return 0;
}

VSET& 
Potential::getAssocVariables()
{
	return varSet;
}

int 
Potential::potZeroInit()
{
	dtree=NULL;
	return 0;
}


int
Potential::calculateConditionalEntropy()
{
	//We are going to use the formular \integral P(X,Mx)logP(X|Mx). 
	//However, we do not know how to estimate P(X,Mx) exactly. P(X|Mx) is given by
	//the regression tree structure. We will approximate P(X,Mx) using two strategies:
	//we will look at the joint configuration of (X,Mx) and treat it as unique and use
	//1/M, where M is the total number of data points in our dataset.
	//In this case the conditional entropy is equal to the sum of the marginal loglikelihood
	//of the child variable weighted by the n_i/M, where n_i is the number of datapoints in 
	//leaf node i. 
	//In the other strategy, we will assume that each configuration Mx which corresponds
	//to the same path in the tree is the same. That is we are partitioning our parent set
	//into P partitions where each partition is a path form the root to a leaf node.
	//We then write P(X,Mx) as P(X|Mx)P(Mx). P(Mx)=P(P_i), where P_i is the path consistent
	//with Mx. This boils down to summing over all the marginal entropies in the leaf nodes
	//weighted by n_i^2/(\sum_j _j).
	//We will begin with this first
	vector<RegressionTree*> allLeafNodes;
	dtree->getLeafNodes(allLeafNodes);
	conditionalEntropy=0;
	jointEntropy=0;
	double totalPaths=evMgr->getNumberOfEvidences();
	for(int i=0;i<allLeafNodes.size();i++)
	{
		RegressionTree* rTree=allLeafNodes[i];
		double aMarginalEntropy=rTree->getMarginalEntropy();
		double pathCnt=(double) rTree->getDataSubset().size();
		double pathProb=pathCnt/totalPaths;
		//conditionalEntropy=aMarginalEntropy* (pathCnt*pathCnt/totalPaths);
		conditionalEntropy=aMarginalEntropy* (pathCnt/totalPaths);
		//jointEntropy=jointEntropy+(pathCnt*pathProb*log(pathProb));
		jointEntropy=jointEntropy+(pathProb*log(pathProb));
	}
	jointEntropy=jointEntropy+conditionalEntropy;
	return 0;
}


double 
Potential::getConditionalEntropy()
{
	return conditionalEntropy;
}

double 
Potential::getJointEntropy()
{
	return jointEntropy;
}

int 
Potential::populateMe(double regval)
{
	lambda=regval;
	dtree=new RegressionTree;
	dtree->setType(RegressionTree::LEAF);
	dtree->setCodingLength(1.0);
	currentLeafNodes.push(dtree);
	for(INTINTMAP_ITER mbIter=markovBlnktVariables.begin();mbIter!=markovBlnktVariables.end();mbIter++)
	{
		dtree->setSubtreeVariable(mbIter->first);
	}
	
	int evidCnt=evMgr->getNumberOfEvidences();
	for(int i=0;i<evidCnt;i++)
	{
		dtree->setDataID(i);
	}
	estimateMarginal(dtree);
	double currentScore=dtree->getMarginalEntropy();
	int splitCnt=1;
	map<int,int> nodeLevel;
	while(!currentLeafNodes.empty())
	{
		double maxGain=0;
		int testVarID=-1;
		double testValue=0;

		RegressionTree* currNode=currentLeafNodes.front();
		INTINTMAP& currSubset=currNode->getDataSubset();
		currentLeafNodes.pop();
		int pNode;
		int branch;
		currNode->getParentInfo(pNode,branch);
		//double penalty=lambda*log(nodeLevel.size()+1);
		//double penalty=lambda*log(splitCnt);
		//cout <<"Splitting child  node of " << pNode <<" branch " << branch << endl;
		INTINTMAP& subtreeVars=currNode->getSubtreeVariables();
		for(INTINTMAP_ITER mbIter=subtreeVars.begin();mbIter!=subtreeVars.end();mbIter++)
		{
			//Partition the training set into sets based on the values of mbIter->first
			//but because this is real, we consider two splits
			//0 corresponds to the partition using testVar<threshold and 1 corresponds to
			//the partition using testVar>threshold
			double gain;
			double splitValue;
			int varId=mbIter->first;
			getPartitions_Cached(mbIter->first,splitValue,currNode->getMarginalEntropy(),currSubset,gain);
			//getPartitions(mbIter->first,splitValue,currNode->getMarginalEntropy(),currSubset,gain);
			if(gain> maxGain)
			{
				maxGain=gain;
				testVarID=mbIter->first;
				testValue=splitValue;
			}
		}
		if(testVarID==-1)
		{
			clearCache();
			//cout << "No good split found " << testVarID << endl;
			continue;
		}
		int testNodeLevel=0;
		if(pNode!=-1)
		{
			int plevel=nodeLevel[pNode];
			testNodeLevel=plevel+1;
		}
		double penalty=2*lambda*log(1+testNodeLevel);
		/*if((maxGain-penalty)<=0)
		{
			clearCache();
			continue;
		}*/
		PARTITION* p=allPartitions[testVarID];
		INTINTMAP* ss1=(*p)[0];
		INTINTMAP* ss2=(*p)[1];
		if( (ss1->size()<minLeafSize) && (ss2->size()<minLeafSize))
		{
			clearCache();
			continue;
		}
		nodeLevel[testVarID]=testNodeLevel;
		splitCnt++;
		Variable* var=varSet[testVarID];
		//cout <<"Gain: "<<maxGain << " for test var "<< testVarID << " " << var->getName().c_str() << " at level " << nodeLevel[testVarID] << endl;
		currNode->setTestVariable(testVarID);
		currNode->setPenalizedScore(maxGain-penalty);
		//Make new leaf nodes using the values of testValue
		currNode->setType(RegressionTree::NONLEAF);
		currNode->split(testValue,allPartitions[testVarID]);
		currNode->setChildParams(allMeans[testVarID],allVariances[testVarID],allMarginalEntropy[testVarID]);
		computeCodingLength(currNode);
		map<int,RegressionTree*>& newLeaves=currNode->getChildren();
		for(map<int,RegressionTree*>::iterator lIter=newLeaves.begin();lIter!=newLeaves.end();lIter++)
		{
			//if((!lIter->second->isPureNode()) && (lIter->second->getDataSubset().size()>10) && (testNodeLevel<1))
			if((!lIter->second->isPureNode()) && (lIter->second->getDataSubset().size()>minLeafSize))
			//if(!lIter->second->isPureNode())
			{
				currentLeafNodes.push(lIter->second);
			}
			lIter->second->setCodingLength(1.0);
		}
		clearCache();
	}
	return 0;
}

int
Potential::clearMe()
{
	if(dtree!=NULL)
	{
		dtree->clear();
		delete dtree;
	}
	dtree=NULL;
	return 0;
}

double
Potential::getJointPotValueFor(INTINTMAP& configMap)
{
	RegressionTree* currNode=dtree;
	double pval=-1;
	while((currNode!=NULL) && (currNode->getType()!=RegressionTree::LEAF))
	{
		int testVar=currNode->getTestVariable();
		int varVal=configMap[testVar];
		RegressionTree* childNode=currNode->getChildAt(varVal);
		currNode=childNode;
	}
	int classVarID=factorVariables.begin()->first;
	int classVarVal=configMap[classVarID];
	if(currNode!=NULL)
	{
		pval=currNode->getMarginalPDF(classVarVal);
	}
	return pval;
}

double 
Potential::getJointPotValueForConf(string& varConf)
{
	cout <<"Not implemented " << endl;
	return 0;
}


int
Potential::generateSample(INTINTMAP& jointConf,int vId, gsl_rng *r)
{
	return 0;
	int sampleVal=-1;
	RegressionTree* currNode=dtree;

	while((currNode!=NULL) && (currNode->getType()!=RegressionTree::LEAF))
	{
		int testVar=currNode->getTestVariable();
		int varVal=jointConf[testVar];
		RegressionTree* childNode=currNode->getChildAt(varVal);
		currNode=childNode;
	}
	if(currNode==NULL)
	{
		return sampleVal;
	}	
	double childMean=currNode->getMean();
	double childVariance=currNode->getVariance();
	double rval=gsl_ran_gaussian(r,sqrt(childVariance));
	//sampleVal=rval+childMean;
	return sampleVal;
}


double 
Potential::predictSample(EMAP* evidMap, int vId)
{
	RegressionTree* currNode=dtree;
	while(currNode->getType()!=RegressionTree::LEAF)
	{
		int testVar=currNode->getTestVariable();
		double varVal=(*evidMap)[testVar]->getEvidVal();
		double testValue=currNode->getTestValue();
		RegressionTree* childNode=NULL;
		if(varVal<=testValue)
		{
			childNode=currNode->getChildAt(0);
		}
		else
		{
			childNode=currNode->getChildAt(1);
		}
		currNode=childNode;
	}
	double predVal=currNode->getMean();
	return predVal;
}

double
Potential::getMSE()
{
	double mse=0;
	queue<RegressionTree*> nodes;
	nodes.push(dtree);
	while(!nodes.empty())
	{
		RegressionTree* cnode=nodes.front();
		nodes.pop();
		map<int,RegressionTree*>& children=cnode->getChildren();
		if(children.size()==0)
		{
			double variance=cnode->getVariance();
			variance=variance*(cnode->getDataSubset().size()-1);
			mse=mse+variance;
		}
		else
		{
			for(map<int,RegressionTree*>::iterator rIter=children.begin();rIter!=children.end();rIter++)
			{
				nodes.push(rIter->second);
			}
		}
	}

	int evidCnt=evMgr->getNumberOfEvidences();
	mse=mse/((double)evidCnt);
	double err2=0;
	int sId=factorVariables.begin()->first;
	for(int i=0;i<evidCnt;i++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt(i);
		double predval=predictSample(evidMap,sId);
		double trueval=(*evidMap)[sId]->getEvidVal();
		double diff=predval-trueval;
		double e=diff*diff;
		err2=err2+e;
	}
	err2=err2/((double) evidCnt);
	cout << "MSE " << mse << " MSE2 " << err2 << endl;

	return mse;
}

int
Potential::dumpPotential(ostream& oFile)
{
	queue<RegressionTree*> nodeList;
	nodeList.push(dtree);
	while(!nodeList.empty())
	{
		RegressionTree* currNode=nodeList.front();
		nodeList.pop();
		map<int,RegressionTree*>& children=currNode->getChildren();
		for(map<int,RegressionTree*>::iterator dIter=children.begin();dIter!=children.end();dIter++)
		{
			oFile << currNode->getTestVariable()<< "\t"<<dIter->first << "\t" << dIter->second->getTestVariable() << endl;		
			nodeList.push(dIter->second);
		}
	}
}


int 
Potential::makeValidJPD()
{
	cout <<"Not implemented yet" << endl;
	return 0;
}

RegressionTree*
Potential::getRegressionTree()
{
	return dtree;
}

int 
Potential::getPartitions(int vId, double& splitValue,double marginalEntropy, INTINTMAP& dataSet, double& infoGain)
{
	//Here we want to partition dataSet into smaller partitions corresponding to the different values of
	//vId. We also want to compute the entropy of the class variable using the partitions. So we store
	//for each value of variable vId, the distribution of the class variable.
	int classVarID=factorVariables.begin()->first;
	vector<int> sortedInd;
	vector<double> sortedValues;
	for(INTINTMAP_ITER dIter=dataSet.begin();dIter!=dataSet.end();dIter++)
	{
		EMAP* evSet=evMgr->getEvidenceAt(dIter->first);
		if(evSet->find(vId)==evSet->end())
		{
			cout <<"No evidence value for " << vId << endl;
			return -1;
		}
		Evidence* evid=(*evSet)[vId];
		double attrVal=evid->getEvidVal();
		sortedInd.push_back(dIter->first);
		sortedValues.push_back(attrVal);
	}
	sortAttrVals(sortedValues,sortedInd);
	int partId=1;
	double maxGain=0;
	int splitId=0;
	double bestEntropyLeft=0;
	double bestMeanLeft=0;
	double bestVarianceLeft=0;
	double bestEntropyRight=0;
	double bestMeanRight=0;
	double bestVarianceRight=0;
	while(partId<sortedValues.size()-1)
	{
		double attrVal=sortedValues[partId];
		/*int leftId=partId+1;
		while(leftId<sortedValues.size() && sortedValues[leftId]==attrVal)
		{
			leftId++;
		}

		partId=leftId-1;
		leftId--;*/
		//Now all datapoints from 0 to partID+leftID are less than or equal to attrVal
		double classEntrLeft=0;
		double meanLeft=0;
		double varianceLeft=0;
		getSubEntropy(sortedInd,0,partId,meanLeft,varianceLeft,classEntrLeft);
		double classEntrRight=0;
		double meanRight=0;
		double varianceRight=0;
		getSubEntropy(sortedInd,partId+1,sortedValues.size()-1,meanRight,varianceRight,classEntrRight);
		//Now compute the new gain
		double weightedEntropy=0;
		double total1=((double)(partId+1));
		double frac1=total1/((double)dataSet.size());
		double frac2=1-frac1;
		weightedEntropy=(frac1*classEntrLeft)+ (frac2*classEntrRight);
		double currInfoGain=marginalEntropy-weightedEntropy;
		if(currInfoGain>maxGain)
		{
			maxGain=currInfoGain;
			splitId=partId;
			splitValue=sortedValues[partId];
			bestEntropyLeft=classEntrLeft;
			bestMeanLeft=meanLeft;
			bestVarianceLeft=varianceLeft;
			bestEntropyRight=classEntrRight;
			bestMeanRight=meanRight;
			bestVarianceRight=varianceRight;
		}
		partId=partId+1;
	}
	infoGain=maxGain;
	PARTITION* partition=new PARTITION;
	INTDBLMAP* mean=new INTDBLMAP; 
	INTDBLMAP* variance=new INTDBLMAP; 
	INTDBLMAP* entropy=new INTDBLMAP;

	for(int i=0;i<sortedInd.size();i++)
	{
		int pseudoVal=-1;
		if(i<=splitId)
		{
			pseudoVal=0;
		}
		else
		{
			pseudoVal=1;
		}
		INTINTMAP* apart=NULL;
		if(partition->find(pseudoVal)==partition->end())
		{
			apart=new INTINTMAP;
			(*partition)[pseudoVal]=apart;
		}
		else
		{
			apart=(*partition)[pseudoVal];
		}
		(*apart)[sortedInd[i]]=0;
	}
	(*mean)[0]=bestMeanLeft;
	(*mean)[1]=bestMeanRight;
	(*variance)[0]=bestVarianceLeft;
	(*variance)[1]=bestVarianceRight;
	(*entropy)[0]=bestEntropyLeft;
	(*entropy)[1]=bestEntropyRight;
	allPartitions[vId]=partition;
	allMeans[vId]=mean;
	allVariances[vId]=variance;
	allMarginalEntropy[vId]=entropy;
	return 0;

}

int 
Potential::getDiscretePartitions(int vId, int value, double marginalEntropy, INTINTMAP& varSet, double& infoGain)
{
	//Here we want to partition dataSet into smaller partitions corresponding to the different values of
	//vId. We also want to compute the entropy of the class variable using the partitions. So we store
	//for each value of variable vId, the distribution of the class variable.
	double maxGain=0;
	int splitId=0;
	double bestEntropyLeft=0;
	double bestMeanLeft=0;
	double bestVarianceLeft=0;
	double bestEntropyRight=0;
	double bestMeanRight=0;
	double bestVarianceRight=0;
	for(INTINTMAP_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		Variable* gene=vMgr->ge
		int attrVal=mMgr.getValue(vId);
		//Now all datapoints from 0 to partID+leftID are less than or equal to attrVal
		double classEntrLeft=0;
		double meanLeft=0;
		double varianceLeft=0;
		getSubEntropy(sortedInd,0,partId,meanLeft,varianceLeft,classEntrLeft);
		double classEntrRight=0;
		double meanRight=0;
		double varianceRight=0;
		getSubEntropy(sortedInd,partId+1,sortedValues.size()-1,meanRight,varianceRight,classEntrRight);
		//Now compute the new gain
		double weightedEntropy=0;
		double total1=((double)(partId+1));
		double frac1=total1/((double)dataSet.size());
		double frac2=1-frac1;
		weightedEntropy=(frac1*classEntrLeft)+ (frac2*classEntrRight);
		double currInfoGain=marginalEntropy-weightedEntropy;
		if(currInfoGain>maxGain)
		{
			maxGain=currInfoGain;
			splitId=partId;
			splitValue=sortedValues[partId];
			bestEntropyLeft=classEntrLeft;
			bestMeanLeft=meanLeft;
			bestVarianceLeft=varianceLeft;
			bestEntropyRight=classEntrRight;
			bestMeanRight=meanRight;
			bestVarianceRight=varianceRight;
		}
		partId=partId+1;
	}
	infoGain=maxGain;
	PARTITION* partition=new PARTITION;
	INTDBLMAP* mean=new INTDBLMAP; 
	INTDBLMAP* variance=new INTDBLMAP; 
	INTDBLMAP* entropy=new INTDBLMAP;

	for(int i=0;i<sortedInd.size();i++)
	{
		int pseudoVal=-1;
		if(i<=splitId)
		{
			pseudoVal=0;
		}
		else
		{
			pseudoVal=1;
		}
		INTINTMAP* apart=NULL;
		if(partition->find(pseudoVal)==partition->end())
		{
			apart=new INTINTMAP;
			(*partition)[pseudoVal]=apart;
		}
		else
		{
			apart=(*partition)[pseudoVal];
		}
		(*apart)[sortedInd[i]]=0;
	}
	(*mean)[0]=bestMeanLeft;
	(*mean)[1]=bestMeanRight;
	(*variance)[0]=bestVarianceLeft;
	(*variance)[1]=bestVarianceRight;
	(*entropy)[0]=bestEntropyLeft;
	(*entropy)[1]=bestEntropyRight;
	allPartitions[vId]=partition;
	allMeans[vId]=mean;
	allVariances[vId]=variance;
	allMarginalEntropy[vId]=entropy;
	return 0;

}


int 
Potential::getPartitions_Cached(int vId, double& splitValue,double marginalEntropy, INTINTMAP& dataSet, double& infoGain)
{
	//Here we want to partition dataSet into smaller partitions corresponding to the different values of
	//vId. We also want to compute the entropy of the class variable using the partitions. So we store
	//for each value of variable vId, the distribution of the class variable.
	int classVarID=factorVariables.begin()->first;
	vector<int> sortedInd;
	vector<double> sortedValues;
	for(INTINTMAP_ITER dIter=dataSet.begin();dIter!=dataSet.end();dIter++)
	{
		EMAP* evSet=evMgr->getEvidenceAt(dIter->first);
		if(evSet->find(vId)==evSet->end())
		{
			cout <<"No evidence value for " << vId << endl;
			return -1;
		}
		Evidence* evid=(*evSet)[vId];
		double attrVal=evid->getEvidVal();
		sortedInd.push_back(dIter->first);
		sortedValues.push_back(attrVal);
	}
	sortAttrVals(sortedValues,sortedInd);
	int partId=1;
	double maxGain=0;
	int splitId=0;
	double bestEntropyLeft=0;
	double bestMeanLeft=0;
	double bestVarianceLeft=0;
	double bestEntropyRight=0;
	double bestMeanRight=0;
	double bestVarianceRight=0;
	bool first=true;
	double sumLeft;
	double sumsqLeft;
	double sumRight;
	double sumsqRight;

	while(partId<sortedValues.size()-2)
	{
		double attrVal=sortedValues[partId];
		/*int leftId=partId+1;
		while(leftId<sortedValues.size() && sortedValues[leftId]==attrVal)
		{
			cout <<"leftid increamenting" << endl;
			leftId++;
		}

		partId=leftId-1;
		leftId--;*/
		//Now all datapoints from 0 to partID+leftID are less than or equal to attrVal
		int dId=sortedInd[partId];
		EMAP* evMap=evMgr->getEvidenceAt(dId);
		double classVarVal=(*evMap)[classVarID]->getEvidVal();

		double classEntrLeft=0;
		double meanLeft=0;
		double varianceLeft=0;
		if(first)
		{
			getSubEntropy(sortedInd,0,partId,meanLeft,varianceLeft,classEntrLeft);
			sumLeft=meanLeft*(partId+1);
			sumsqLeft=(varianceLeft*partId) + (meanLeft*meanLeft*(partId+1));
		}
		else
		{
			double newsumLeft=sumLeft+classVarVal;
			double newsumsqLeft=sumsqLeft+(classVarVal*classVarVal);
			meanLeft=newsumLeft/(partId+1);
			varianceLeft=(((partId+1)*meanLeft*meanLeft) + newsumsqLeft - (2*newsumLeft*meanLeft))/partId;
			/*double meanLeft1=0;
			double varianceLeft1=0;
			getSubEntropy(sortedInd,0,partId,meanLeft1,varianceLeft1,classEntrLeft);
			if(((meanLeft-meanLeft1)!=0) || ((varianceLeft-varianceLeft1)!=0))
			{
				cout <<"LEFT mismatch! ";
				cout <<"Cached m: " << meanLeft << " v: " << varianceLeft << " true m: "<< meanLeft1 << " v: " << varianceLeft1 << endl;
			}*/
			getSubEntropy(meanLeft,varianceLeft,classEntrLeft);
			sumLeft=newsumLeft;
			sumsqLeft=newsumsqLeft;
		}
		double classEntrRight=0;
		double meanRight=0;
		double varianceRight=0;
		double elemCnt=(double) (sortedValues.size()-partId-1);
		if(first)
		{
			getSubEntropy(sortedInd,partId+1,sortedValues.size()-1,meanRight,varianceRight,classEntrRight);
			sumRight=meanRight*(elemCnt);
			sumsqRight=(varianceRight*(elemCnt-1)) + (meanRight*meanRight*elemCnt);
		}
		else
		{
			double newsumRight=sumRight-classVarVal;
			double newsumsqRight=sumsqRight-(classVarVal*classVarVal);
			meanRight=newsumRight/elemCnt;
			varianceRight=((elemCnt*meanRight*meanRight) + (newsumsqRight) - (2*newsumRight*meanRight))/(elemCnt-1);
		/*	double meanRight1=0;
			double varianceRight1=0;

			getSubEntropy(sortedInd,partId+1,sortedValues.size()-1,meanRight1,varianceRight1,classEntrRight);
			if(((meanRight-meanRight1)!=0) || ((varianceRight-varianceRight1)!=0))
			{
				cout <<"RIGHT mismatch! ";
				cout <<"Cached m: " << meanRight << " v: " << varianceRight << " true m: "<< meanRight1 << " v: " << varianceRight1 << endl;
			}*/
			getSubEntropy(meanRight,varianceRight,classEntrRight);
			sumRight=newsumRight;
			sumsqRight=newsumsqRight;
		}
		//Now compute the new gain
		double weightedEntropy=0;
		double total1=((double)(partId+1));
		double frac1=total1/((double)dataSet.size());
		double frac2=1-frac1;
		weightedEntropy=(frac1*classEntrLeft)+ (frac2*classEntrRight);
		double currInfoGain=marginalEntropy-weightedEntropy;
		if(currInfoGain>maxGain)
		{
			maxGain=currInfoGain;
			splitId=partId;
			splitValue=sortedValues[partId];
			bestEntropyLeft=classEntrLeft;
			bestMeanLeft=meanLeft;
			bestVarianceLeft=varianceLeft;
			bestEntropyRight=classEntrRight;
			bestMeanRight=meanRight;
			bestVarianceRight=varianceRight;
		}
		partId=partId+1;
		first=false;
	}
	infoGain=maxGain;
	PARTITION* partition=new PARTITION;
	INTDBLMAP* mean=new INTDBLMAP; 
	INTDBLMAP* variance=new INTDBLMAP; 
	INTDBLMAP* entropy=new INTDBLMAP;

	for(int i=0;i<sortedInd.size();i++)
	{
		int pseudoVal=-1;
		if(i<=splitId)
		{
			pseudoVal=0;
		}
		else
		{
			pseudoVal=1;
		}
		INTINTMAP* apart=NULL;
		if(partition->find(pseudoVal)==partition->end())
		{
			apart=new INTINTMAP;
			(*partition)[pseudoVal]=apart;
		}
		else
		{
			apart=(*partition)[pseudoVal];
		}
		(*apart)[sortedInd[i]]=0;
	}
	(*mean)[0]=bestMeanLeft;
	(*mean)[1]=bestMeanRight;
	(*variance)[0]=bestVarianceLeft;
	(*variance)[1]=bestVarianceRight;
	(*entropy)[0]=bestEntropyLeft;
	(*entropy)[1]=bestEntropyRight;
	allPartitions[vId]=partition;
	allMeans[vId]=mean;
	allVariances[vId]=variance;
	allMarginalEntropy[vId]=entropy;
	return 0;

}

int
Potential::clearCache()
{
	for(map<int,PARTITION*>::iterator pIter=allPartitions.begin();pIter!=allPartitions.end();pIter++)
	{
		PARTITION* sPart=pIter->second;
		for(PARTITION_ITER spIter=sPart->begin();spIter!=sPart->end();spIter++)
		{
			spIter->second->clear();
			delete spIter->second;
		}
		pIter->second->clear();
		delete pIter->second;
	}

	for(map<int,INTDBLMAP*>::iterator aIter=allMeans.begin();aIter!=allMeans.end();aIter++)
	{
		INTDBLMAP* mean=aIter->second;
		INTDBLMAP* variance=allVariances[aIter->first];
		INTDBLMAP* entropy=allMarginalEntropy[aIter->first];
		mean->clear();
		variance->clear();
		entropy->clear();
		delete mean;
		delete variance;
		delete entropy;
	}
	allPartitions.clear();
	allMeans.clear();
	allVariances.clear();
	allMarginalEntropy.clear();
}

int 
Potential::estimateMarginal(RegressionTree* dtree)
{
	int classVarID=factorVariables.begin()->first;
	int evidCnt=evMgr->getNumberOfEvidences();
	vector<int> dataInd;
	for(int i=0;i<evidCnt;i++)
	{
		dataInd.push_back(i);
	}
	double marginalEntropy=0;
	double mean=0;
	double variance=0;
	getSubEntropy(dataInd,0,evidCnt-1,mean,variance,marginalEntropy);
	dtree->setMean(mean);
	dtree->setVariance(variance);
	dtree->setMarginalEntropy(marginalEntropy);
	return 0;
}

int
Potential::sortAttrVals(vector<double>& sortedValues,vector<int>& sortedInd)
{
/*	for(int i=0;i<10;i++)
	{
		cout << i <<" " << sortedValues[i] << endl;
	}*/
	for(int i=0;i<sortedValues.size();i++)
	{
		for(int j=i+1;j<sortedValues.size();j++)
		{
			double aVal=sortedValues[i];
			double bVal=sortedValues[j];
			if(aVal<bVal)
			{
				continue;
			}
			sortedValues[i]=bVal;
			sortedValues[j]=aVal;
			int tempRank=sortedInd[i];
			sortedInd[i]=sortedInd[j];
			sortedInd[j]=tempRank;

		}
	}
	//check for sort
	/*for(int i=0;i<10;i++)
	{
		cout <<i << " " << sortedValues[i] <<" " << sortedValues[sortedInd[i]]<< endl;
	}*/
	return 0;
}

//Assuming that the classVarID follows a Gaussian distribution, we can compute the marginal entropy 
//in closed form
int
Potential::getSubEntropy(vector<int>& sortedInd, int start, int end,double& mean,double& variance,double& marginalEntropy)
{
	marginalEntropy=0;
	int classVarID=factorVariables.begin()->first;
	mean=0;
	for(int i=start;i<=end;i++)
	{
		int dId=sortedInd[i];
		EMAP* evMap=evMgr->getEvidenceAt(dId);
		double classVarVal=(*evMap)[classVarID]->getEvidVal();
		mean=mean+classVarVal;
	}
	mean=mean/((double)(end-start+1));
	variance=0;
	for(int i=start;i<=end;i++)
	{
		int dId=sortedInd[i];
		EMAP* evMap=evMgr->getEvidenceAt(dId);
		double classVarVal=(*evMap)[classVarID]->getEvidVal();
		double diff=mean-classVarVal;
		variance=variance+(diff*diff);
	}
	variance=variance/(end-start);
	double determinant=variance;
	double commFact=1+log(2*PI);
	double n=1;
	marginalEntropy=0.5*((n*commFact) + log(determinant));
	return 0;
}


int
Potential::getSubEntropy(double mean,double variance,double& marginalEntropy)
{
	marginalEntropy=0;
	double determinant=variance;
	double commFact=1+log(2*PI);
	double n=1;
	marginalEntropy=0.5*((n*commFact) + log(determinant));
	return 0;
}


int
Potential::showTree()
{
	stack<RegressionTree*> nodeList;
	nodeList.push(dtree);
	int currLevel=0;
	int nodeID=0;
	int factorVariable=factorVariables.begin()->first;
	while(!nodeList.empty())
	{
		RegressionTree* currNode=nodeList.top();
		nodeList.pop();
		map<int,RegressionTree*>& children=currNode->getChildren();
		int vId=currNode->getTestVariable();
		if(vId==-1)
		{
			vId=factorVariable;
		}
		double testValue=currNode->getTestValue();
		VSET_ITER vIter=varSet.find(vId);
		if(vIter==varSet.end())
		{
			cout <<"No variable with id " << vId << endl;
			return -1;
		}
		Variable* v=vIter->second;
		if(vId!=factorVariable)
		{
			cout << v->getName().c_str() <<" <= " << testValue << endl;
		}
		else
		{
			cout << v->getName().c_str() << endl;
		}
		for(map<int,RegressionTree*>::iterator dIter=children.begin();dIter!=children.end();dIter++)
		{
			RegressionTree* cNode=dIter->second;
			nodeList.push(cNode);
		}
	}
	return 0;
}

int
Potential::prune()
{
	dtree->generateRuleSet();
	vector<Rule*>& ruleSet=dtree->getRuleSet();

	int evidCnt=evMgr->getNumberOfEvidences();
	for(int i=0;i<evidCnt;i++)
	{
		pruneDataSet[i]=0;
	}

	for(int i=0;i<ruleSet.size();i++)
	{
		if(pruneDataSet.size()==0)
		{
			continue;
		}
		Rule* arule=ruleSet[i];
		if(arule==NULL)
		{
			continue;
		}
		//Keep trying to delete conditions from this rule until the entropy does not
		//decrease
		map<int,RegressionTree*>& conditionSet=arule->getAllConditions();
		bool foundDelCondition=true;
		while(foundDelCondition)
		{
			//double currEntropy=arule->getMarginalEntropy();
			//int currCoverage=arule->getCoverage();
			int currCoverage=0;
			double currEntropy=computeEntropyIfDeleted(arule,-1,currCoverage);
			cout <<"Pruning rule " << i << " with coverage: " << currCoverage << endl;
			arule->showRule();
			int oldCoverage=currCoverage;
			double oldComplexity=arule->getRuleComplexity(-1);
			double currComplexity=oldComplexity;
			if(currCoverage==0)
			{
				foundDelCondition=false;
				delete arule;
				ruleSet[i]=NULL;
				arule=NULL;
				continue;
			}
			int toDel=-1;
			for(map<int,RegressionTree*>::iterator lIter=conditionSet.begin();lIter!=conditionSet.end();lIter++)
			{
				if(lIter->second->getChildren().size()==0)
				{
					continue;
				}
				int newCoverage;
				double newEntropy=computeEntropyIfDeleted(arule,lIter->first,newCoverage);
				if(newCoverage<oldCoverage)
				{
					cout <<"Wierd rule. Deleting condition does not increase coverage!"<<endl;
				}
				//We are guaranteed that ruleLen>1
				double ruleLen=(double)conditionSet.size();
				double newComplexity=arule->getRuleComplexity(lIter->first);
				if(((newCoverage*newEntropy)+newComplexity)<=((currCoverage*currEntropy)+currComplexity))
				{
					currEntropy=newEntropy;
					currCoverage=newCoverage;
					toDel=lIter->first;
					currComplexity=newComplexity;
				}
			}
			if(toDel==-1)
			{
				foundDelCondition=false;
				continue;
			}
			map<int,RegressionTree*>::iterator delIter=conditionSet.find(toDel);
			conditionSet.erase(delIter);
			//If the condition set has only one node, it means it is the leaf node and can be deleted.
			if(conditionSet.size()==1)
			{
				foundDelCondition=false;
				delete arule;
				ruleSet[i]=NULL;
				arule=NULL;
				continue;
			}
			arule->setMarginalEntropy(currEntropy);
			arule->setCoverage(currCoverage);
		}
		if(arule!=NULL)
		{
			prunedRuleSet.push_back(arule);
			reducePruneSet(arule);
			//If I am using the Nguyen approach then I need to delete the datapoints
			//already covered by this rule
		}
	}
	cout <<"Number of rules after pruning" << prunedRuleSet.size() << endl;
	removeDuplicateCoverage();
	return 0;
}

int
Potential::getAssocVariables_PostPruning(INTINTMAP& newVarSet)
{
	for(int r=0;r<prunedRuleSet.size();r++)
	{
		Rule* arule=prunedRuleSet[r];
		if(arule==NULL)
		{
			continue;
		}
		map<int,RegressionTree*>& conditionSet=arule->getAllConditions();
		for(map<int,RegressionTree*>::iterator lIter=conditionSet.begin();lIter!=conditionSet.end();lIter++)
		{
			RegressionTree* node=lIter->second;
			if(node->getChildren().size()==0)
			{
				continue;
			}
			int testVarID=node->getTestVariable();
			newVarSet[testVarID]=0;
		}
	}
	return 0;
}


double
Potential::computeEntropyIfDeleted(Rule* arule,int delMe,int& acoverage)
{
	vector<double> filteredVals;
	int evidCnt=evMgr->getNumberOfEvidences();
	double mean=0;
	double var=0;
	int classVarID=factorVariables.begin()->first;
	map<int,RegressionTree*>& conditionSet=arule->getAllConditions();

	//for(int e=0;e<evidCnt;e++)
	for(INTINTMAP_ITER aIter=pruneDataSet.begin();aIter!=pruneDataSet.end();aIter++)
	{
		int e=aIter->first;
		EMAP* evMap=evMgr->getEvidenceAt(e);
		double classVarVal=(*evMap)[classVarID]->getEvidVal();
		bool testTrue=true;
		for(map<int,RegressionTree*>::iterator lIter=conditionSet.begin();lIter!=conditionSet.end();lIter++)
		{
			RegressionTree* rnode=lIter->second;
			if(lIter->first==delMe)
			{
				continue;
			}
			if(rnode->getChildren().size()==0)
			{
				continue;
			}
			int parentBranch=arule->getBranch(lIter->first);
			if(parentBranch==-1)
			{
				cout <<"Invalid branch value " << endl;
				exit(0);
			}
			int testVar=rnode->getTestVariable();
			double testValue=rnode->getTestValue();
			double varVal=(*evMap)[testVar]->getEvidVal();
			if(parentBranch==0)
			{
				if(varVal>testValue)
				{
					testTrue=false;
				}
			}
			else if(parentBranch==1)
			{
				if(varVal<=testValue)
				{
					testTrue=false;
				}
			}
			if(!testTrue)
			{
				break;
			}
		}
		if(testTrue)
		{
			filteredVals.push_back(classVarVal);
			mean=mean+classVarVal;
		}
	}
	//cout <<"Datapoints covered " << filteredVals.size() << endl;
	mean=mean/((double)filteredVals.size());
	for(int i=0;i<filteredVals.size();i++)
	{
		double diff=mean-filteredVals[i];
		var=var+(diff*diff);
	}
	var=var/((double)filteredVals.size()-1);
	double newEntropy=-1;
	getSubEntropy(mean,var,newEntropy);
	acoverage=filteredVals.size();
	filteredVals.clear();
	return newEntropy;
}

int
Potential::reducePruneSet(Rule* r)
{
	INTINTMAP coveredSet;
	map<int,RegressionTree*>& conditionSet=r->getAllConditions();
	for(INTINTMAP_ITER aIter=pruneDataSet.begin();aIter!=pruneDataSet.end();aIter++)
	{
		EMAP* evMap=evMgr->getEvidenceAt(aIter->first);
		bool testTrue=true;
		for(map<int,RegressionTree*>::iterator lIter=conditionSet.begin();lIter!=conditionSet.end();lIter++)
		{
			RegressionTree* rnode=lIter->second;
			if(rnode->getChildren().size()==0)
			{
				continue;
			}
			int parentBranch=r->getBranch(lIter->first);
			if(parentBranch==-1)
			{
				cout <<"Invalid branch value " << endl;
				exit(0);
			}
			int testVar=rnode->getTestVariable();
			double testValue=rnode->getTestValue();
			double varVal=(*evMap)[testVar]->getEvidVal();
			if(parentBranch==0)
			{
				if(varVal>testValue)
				{
					testTrue=false;
				}
			}
			else if(parentBranch==1)
			{
				if(varVal<=testValue)
				{
					testTrue=false;
				}
			}
			if(!testTrue)
			{
				break;
			}
		}
		if(testTrue)
		{
			coveredSet[aIter->first]=0;
			INTINTMAP* coveredByRules=NULL;
			if(dataRuleCoverage.find(aIter->first)==dataRuleCoverage.end())
			{
				coveredByRules=new INTINTMAP;
				dataRuleCoverage[aIter->first]=coveredByRules;
			}
			else
			{
				coveredByRules=dataRuleCoverage[aIter->first];
			}
			(*coveredByRules)[prunedRuleSet.size()-1]=0;
		}
	}

	for(INTINTMAP_ITER dIter=coveredSet.begin();dIter!=coveredSet.end();dIter++)
	{
		INTINTMAP_ITER eIter=pruneDataSet.find(dIter->first);
		if(eIter==pruneDataSet.end())
		{
			cout <<"No datapoint found " << dIter->first <<endl;
			exit(0);
		}
		//pruneDataSet.erase(eIter);
	}
	coveredSet.clear();
	return 0;
}


int 
Potential::computeCodingLength(RegressionTree* nonleafNode)
{
	int evidCnt=evMgr->getNumberOfEvidences();
	int testVarID=nonleafNode->getTestVariable();
	double minVal=0;
	double maxVal=0;
	//for(INTINTMAP_ITER dIter=dataSubset.begin();dIter!=dataSubset.end();dIter++)
	if(minValMap.find(testVarID)==minValMap.end())
	{
		for(int e=0;e<evidCnt;e++)
		{
			EMAP* evMap=evMgr->getEvidenceAt(e);
			double aVal=(*evMap)[testVarID]->getEvidVal();
			if(e==0)
			{
				minVal=aVal;
				maxVal=aVal;
			}
			else
			{
				if(aVal<minVal)
				{
					minVal=aVal;
				}
				if(aVal>maxVal)
				{
					maxVal=aVal;
				}
			}
		}
		minValMap[testVarID]=minVal;
		maxValMap[testVarID]=maxVal;
	}
	else
	{
		minVal=minValMap[testVarID];
		maxVal=maxValMap[testVarID];
	}
	double codingLen=1.0+log((double)varSet.size())/log(2.0);
	double l=maxVal-minVal;
	codingLen=codingLen+(log(l/0.1)/log(2.0));
	nonleafNode->setCodingLength(codingLen);
	return 0;
}


int 
Potential::removeDuplicateCoverage()
{
	int overlappedData=0;
	for(map<int,INTINTMAP*>::iterator dIter=dataRuleCoverage.begin();dIter!=dataRuleCoverage.end();dIter++)
	{
		INTINTMAP* covRuleSet=dIter->second;
		if(covRuleSet->size()==1)
		{
			continue;
		}
		overlappedData++;
		//If there are multiple rules, then keep the one with the least MDL
		double minEntropy=0;
		double minComplexity=0;
		double coverage=0;
		int bestRule=-1;
		for(INTINTMAP_ITER rIter=covRuleSet->begin();rIter!=covRuleSet->end();rIter++)
		{
			Rule* rule=prunedRuleSet[rIter->first];
			if(rIter==covRuleSet->begin())
			{
				minEntropy=rule->getMarginalEntropy();
				minComplexity=rule->getRuleComplexity(-1);
				coverage=rule->getCoverage()-1;
				bestRule=rIter->first;
			}
			else
			{
				int newCoverage=rule->getCoverage()-1;
				double newComplexity=rule->getRuleComplexity(-1);
				double newEntropy=rule->getMarginalEntropy();
				if((newComplexity+(newCoverage*newEntropy))
				  <(minComplexity+(coverage*minEntropy)))
				{
					minEntropy=newEntropy;
					coverage=newCoverage;
					minComplexity=newComplexity;
					bestRule=rIter->first;
				}
			}
		}
		for(INTINTMAP_ITER rIter=covRuleSet->begin();rIter!=covRuleSet->end();rIter++)
		{
			if(rIter->first==bestRule)
			{
				continue;
			}
			Rule* arule=prunedRuleSet[rIter->first];
			int newCov=arule->getCoverage()-1;
			arule->setCoverage(newCov);
		}
		covRuleSet->clear();
		(*covRuleSet)[bestRule]=0;
	}
	cout <<"Number of datapoints overlapped "<<  overlappedData << endl;
	//Now finally get rid of any rules that do not have any coverage;
	int blankrule=0;
	for(int r=0;r<prunedRuleSet.size();r++)
	{
		if(prunedRuleSet[r]->getCoverage()>0)
		{
			continue;
		}
		blankrule++;
		delete prunedRuleSet[r];
		prunedRuleSet[r]=NULL;
	}
	cout <<"Blank rules " << blankrule << endl;
	return 0;
}
