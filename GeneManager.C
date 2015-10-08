
#include "Path.H"
#include "Gene.H"
#include "GeneManager.H"

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>

GeneManager::GeneManager()
{
	maxExp=0;
	minExp=1e50;
	geneIds=0;
}

GeneManager::~GeneManager()
{
}

int
GeneManager::readGeneData(const char* inFName,int logtrans)
{
	ifstream inFile(inFName);
	char buffer[EXPRFILE_LEN];
	totalMeasurements=0;
	while(inFile.good())
	{
		inFile.getline(buffer,EXPRFILE_LEN-1);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		addGeneNode(buffer,logtrans);
	}
	inFile.close();
	return 0;
}

int
GeneManager::readGeneData(const char* inFName,int logtrans, int datastart)
{
	ifstream inFile(inFName);
	char* buffer=NULL;
	int bufflen=0;
	totalMeasurements=0;
	int linecnt=0;
	string buffstr;
	while(inFile.good())
	{
		getline(inFile,buffstr);
		if(buffstr.length()<=0)
		{
			continue;
		}
		if((buffstr.length()+1)>=bufflen)
		{	
			bufflen=buffstr.length()+1;
			if(buffer!=NULL)
			{
				delete[] buffer;
			}
			buffer=new char[bufflen];
		}
		strcpy(buffer,buffstr.c_str());
		if(linecnt<datastart)
		{
			linecnt++;
			continue;
		}
		addGeneNode(buffer,logtrans);
	}
	inFile.close();
	return 0;
}

//The parameter is a single line of the gene expression file
int 
GeneManager::addGeneNode(const char* geneData,int logtrans)
{
	int len=strlen(geneData)+1;
	char* buffer=new char[len];
	strcpy(buffer,geneData);
	//The data is tab-delimited
	int tokCnt=0;
	char* tok=strtok(buffer,"\t ");
	Gene* aGene=NULL;
	while(tok!=NULL)
	{
		if(tokCnt==0)
		{
			aGene=new Gene;
			aGene->setName(tok);
			aGene->setType(BioNode::GENE);
			aGene->setID(geneIds);
			geneSet.push_back(aGene);
			string aKey(tok);
			geneMap[aKey]=aGene;
			geneIds++;
		}
		else
		{
			double aVal=atof(tok);
			if(logtrans==1)
			{
				aVal=pow(2,aVal);
			}
			if(aVal>maxExp)
			{
				maxExp=aVal;
			}
			if(aVal<minExp)
			{
				minExp=aVal;
			}
			aGene->addExpLevel(aVal);
		}
		tokCnt++;
		tok=strtok(NULL,"\t ");
	}
	if(totalMeasurements==0)
	{
		totalMeasurements=tokCnt-1;
	}
	else 
	{
		if(tokCnt-1 !=totalMeasurements)
		{
			cout <<"Found only " << tokCnt-1 << " out of total " << totalMeasurements  << " of total measurements for " << aGene->getName() << endl;
		}
	}

	return 0;
}

Gene* 
GeneManager::getGeneNode(int id)
{
	if((id<0) || (id>=geneSet.size()))
	{
		cout <<"Bad gene ID given to GeneManager " << endl;
		return NULL;
	}
	return geneSet[id];
}

Gene* 
GeneManager::getGeneWithName(const char* aName)
{
	string aKey(aName);
	if(geneMap.find(aKey)==geneMap.end())
	{
		return NULL;
	}
	return geneMap[aKey];
}

int 
GeneManager::getTotalNumberOfGenes()
{
	return geneSet.size();
}


int 
GeneManager::discretizeGeneExpr()
{
	for(int i=0;i<geneSet.size();i++)
	{
		//geneSet[i]->discretizeLevels();
		//geneSet[i]->discretizeEqualBins();
		geneSet[i]->discretizeEqualBins_Global(minExp,maxExp);
	}
	return 0;
}

int
GeneManager::rankOrderGeneExpr()
{
	for(int i=0;i<geneSet.size();i++)
	{
		//geneSet[i]->rankOrder();
		geneSet[i]->rankOrder_Store();
	}
	return 0;
}

int
GeneManager::scaleGeneExpr()
{
	for(int i=0;i<geneSet.size();i++)
	{
		geneSet[i]->scaleExpression();
	}
	return 0;
}

int
GeneManager::standardizeGeneExpr()
{
	for(int i=0;i<geneSet.size();i++)
	{
		geneSet[i]->standardizeData();
	}
	return 0;
}


int
GeneManager::zeroMeanGeneExpr()
{
	for(int i=0;i<geneSet.size();i++)
	{
		geneSet[i]->zeroMeanData();
	}
	return 0;
}
	
	
int 
GeneManager::assocProtWithCodingGene(const char* protName, int id)
{
	if((id<0) || (id>=geneSet.size()))
	{
		cout<< "Bad Gene Id " << endl;
		return -1;
	}
	Gene* gene=geneSet[id];
	gene->setCodedProteinName(protName);
	return 0;
}
