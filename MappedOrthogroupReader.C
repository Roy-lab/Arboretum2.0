#include <iostream>

#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <unistd.h>
#include "GeneMap.H"
#include "MappedOrthogroup.H"
#include "MappedOrthogroupReader.H"

MappedOrthogroupReader::MappedOrthogroupReader()
{
}

MappedOrthogroupReader::~MappedOrthogroupReader()
{
}

//This function reads the orthogroup mapping produced in any one of the PROPER, FILLED and UNFILLED formats.
//We are essentialy interested in the first two columns
int 
MappedOrthogroupReader::readFile(const char* aFName)
{
	ifstream inFile(aFName);
	char* buffer=NULL;
	int bufflen=0;
	string buffstr;
	int lineCnt=0;
	int oldogid=-2;
	MappedOrthogroup* ogrp=NULL;
	while(inFile.good())
	{
		getline(inFile,buffstr);
		if(lineCnt==0)
		{
			lineCnt++;
			continue;
		}
		if(bufflen<=buffstr.length())
		{
			if(buffer!=NULL)
			{
				delete[] buffer;
			}
			bufflen=buffstr.length()+1;
			buffer=new char[bufflen];
		}
		strcpy(buffer,buffstr.c_str());
		if(strstr(buffer,"orf19.993")!=NULL)
		{
			cout <<"Stop here" << endl;
		}
		//Get the orthogroup and the genes
		int ogid=-1;
		string geneSet;
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				char* pos=strchr(tok,'_');
				if(pos==NULL)
				{
					cout <<"Bad OG name " << tok<<  endl;
					exit(0);
				}
				*pos='\0';
				ogid=atoi(tok+2);
				if(oldogid!=ogid)
				{
					ogrp=new MappedOrthogroup;
					ogrp->setID(ogid);
					orthogroupSet[ogrp->getID()]=ogrp;
					oldogid=ogid;
				}
				ogrp->incrCnt();
				
			}
			else if(tokCnt==1)
			{
				addMembers(tok,ogrp);
			}
			tokCnt++;
			tok=strtok(NULL,"\t");
		}
		lineCnt++;
	}
	inFile.close();
	generateGeneOrthoMap();	
	return 0;
}

int
MappedOrthogroupReader::readSpeciesMapping(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	int lineCnt=0;
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string speciesName(buffer);
		speciesIDNameMap[lineCnt]=speciesName;
		lineCnt++;
	}
	inFile.close();
	return 0;
}




int 
MappedOrthogroupReader::getMappedOrthogroupID(const char* geneName, const char* species)
{
	string speciesKey(species);
	if(geneOrthoMap.find(speciesKey)==geneOrthoMap.end())
	{
		return -1;
	}
	string geneKey(geneName);
	STRINTMAP* geneSet=geneOrthoMap[speciesKey]; 
	if(geneSet->find(geneKey)==geneSet->end())
	{
		return -1;
	}
	int orthoID=(*geneSet)[geneKey];
	return orthoID;
}


map<int,MappedOrthogroup*>& 
MappedOrthogroupReader::getMappedOrthogroups()
{
	return orthogroupSet;
}

MappedOrthogroup*
MappedOrthogroupReader::getMappedOrthogroup(const char* geneName, const char* species)
{
	int orthoID=getMappedOrthogroupID(geneName,species);
	MappedOrthogroup* og=orthogroupSet[orthoID];
	return og;
}


STRINTMAP* 
MappedOrthogroupReader::getOrtholog(const char* srcSpecName, const char* geneName, const char* targetSpecName)
{
	int ogid=getMappedOrthogroupID(geneName,srcSpecName);
	if(ogid==-1)
	{
		return NULL;
	}
	MappedOrthogroup* mgrp=orthogroupSet[ogid];
	STRINTMAP* orthohits=mgrp->getSpeciesHitsForGene(srcSpecName,targetSpecName,geneName);
	return orthohits;
}

int
MappedOrthogroupReader::addMembers(char* abuffer, MappedOrthogroup* ogrp)
{
	char* tok=abuffer;
	int tokCnt=0;
	map<string,string> specGeneMap;
	while(tok!=NULL)
	{
		char* end=strchr(tok,',');
		if(end!=NULL)
		{
			*end='\0';
		}
		if(speciesIDNameMap.find(tokCnt)==speciesIDNameMap.end())
		{
			cout <<"No species with id " << tokCnt << endl;
			exit(0);
		}
		string speciesName=speciesIDNameMap[tokCnt];
		string geneName(tok);
		if(strcmp(geneName.c_str(),"NONE")!=0)
		{
			specGeneMap[speciesName]=geneName;
		}
		tokCnt++;
		if(end!=NULL)
		{
			tok=end+1;
		}	
		else
		{
			tok=end;
		}
	}
	ogrp->setMembers(specGeneMap);
	specGeneMap.clear();
	return 0;
}

int 
MappedOrthogroupReader::generateGeneOrthoMap()
{
	for(map<int,MappedOrthogroup*>::iterator oIter=orthogroupSet.begin();oIter!=orthogroupSet.end();oIter++)
	{
		MappedOrthogroup* og=oIter->second;
		map<string,GeneMap*>& members=og->getOrthoMembers();
		for(map<string,GeneMap*>::iterator sIter=members.begin();sIter!=members.end();sIter++)
		{
			map<string,map<string,STRINTMAP*>*>& genes=sIter->second->getGeneSet();
			STRINTMAP* allgenesForSpecies=NULL;
			if(geneOrthoMap.find(sIter->first)==geneOrthoMap.end())
			{
				allgenesForSpecies=new STRINTMAP;
				geneOrthoMap[sIter->first]=allgenesForSpecies;
			}
			else
			{
				allgenesForSpecies=geneOrthoMap[sIter->first];
			}
			for(map<string,map<string,STRINTMAP*>*>::iterator gIter=genes.begin();gIter!=genes.end();gIter++)
			{
				if(strcmp(gIter->first.c_str(),"orf19.993")==0)
				{
					cout <<"Found " << gIter->first.c_str() << endl;
				}
				(*allgenesForSpecies)[gIter->first]=oIter->first; 
			}
		}
	}
	return 0;
}
