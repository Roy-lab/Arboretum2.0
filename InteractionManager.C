#include <fstream>
#include <iostream>
#include <string.h>

#include "Interaction.H"
#include "InteractionManager.H"

InteractionManager::InteractionManager()
{
}

InteractionManager::~InteractionManager()
{
}

int
InteractionManager::setDirectionality(bool dirFlag)
{
	directionality=dirFlag;
	return 0;
}

int
InteractionManager::setDelimiter(const char* aDel)
{
	strcpy(delim,aDel);
	return 0;
}

int
InteractionManager::readInteractions(const char* fileName)
{
	ifstream inFile(fileName);
	char buffer[1024];
	int totalIntr=0;
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		int tokCnt=0;
		char* tok=strtok(buffer,delim);
		char firstGene[GENENAME_LEN];
		char secondGene[GENENAME_LEN];
		string actiontype;
		while(tok!=NULL)
		{
			switch(tokCnt)
			{
				case 0:
				{
					strcpy(firstGene,tok);
					break;
				}
				case 1:
				{
					actiontype.append(tok);
					break;
				}
				case 2:
				{
					strcpy(secondGene,tok);
					break;
				}
			}
			tok=strtok(NULL,delim);
			tokCnt++;
		}
		string aKey;
		createKey(firstGene,secondGene,aKey);
		totalIntr++;
		if(interactionSet.find(aKey)==interactionSet.end())
		{
			Interaction* intr=new Interaction;
			intr->setMembers(firstGene,secondGene);
			intr->setProperty("actiontype",actiontype.c_str());
			interactionSet[aKey]=intr;
			intrIDtoString[interactionSet.size()-1]=aKey;
		}
	}
	inFile.close();

	cout << "Added " << interactionSet.size() << " out of a total " << totalIntr <<  " interactions"<< endl;
	return 0;
}

INTR_MAP&
InteractionManager::getInteractions()
{
	return interactionSet;
}

bool
InteractionManager::isInteraction(const char* firstOrf, const char* secondOrf)
{
	string aKey;
	createKey(firstOrf,secondOrf,aKey);
	if(interactionSet.find(aKey)==interactionSet.end())
	{
		return false;
	}
	return true;
}

int
InteractionManager::addInteraction(const char* firstOrf, const char* secondOrf)
{
	string aKey;
	createKey(firstOrf,secondOrf,aKey);
	if(interactionSet.find(aKey)==interactionSet.end())
	{
		Interaction* intr=new Interaction;
		intr->setMembers(firstOrf,secondOrf);
		interactionSet[aKey]=intr;
	}
	return 0;
}


int
InteractionManager::getNumberOfIntr()
{
	return interactionSet.size();
}

Interaction*
InteractionManager::getInteractionAt(int intrId)
{
	map<int,string>::iterator aIter=intrIDtoString.find(intrId);
	if(aIter==intrIDtoString.end())
	{
		return NULL;
	}
	return interactionSet[aIter->second];
}

int
InteractionManager::dumpInteractions(const char* ofName)
{
	ofstream oFile(ofName);
	for(INTR_MAP_ITER imIter=interactionSet.begin();imIter!=interactionSet.end();imIter++)
	{
		oFile <<imIter->first.c_str() <<  endl;
	}
	oFile.close();
	return 0;
}

int
InteractionManager::createKey(const char* firstGene, const char* secondGene, string& aKey)
{
	char keyChar[1024];
	if(directionality)
	{
		/*int len1=strlen(firstGene);
		int len2=strlen(secondGene);
		int compLen=len1;
		if(len2<len1)
		{
			compLen=len2;
		}
		int i=0;
		int smallId=0;
		while((i<compLen) && (smallId==0))
		{
			if(firstGene[i]<secondGene[i])
			{
				smallId=1;
			}
			else if(firstGene[i]>secondGene[i])
			{
				smallId=2;
			}
			i++;
		}
		if(smallId==1)
		{
			sprintf(keyChar,"%s-%s",firstGene,secondGene);
		}
		else if(smallId==2)
		{
			sprintf(keyChar,"%s-%s",secondGene,firstGene);
		}
		else
		{
			//cout << "First and second genes are the same " << firstGene << endl;
			sprintf(keyChar,"%s-%s",secondGene,firstGene);
		}*/
		if(strcmp(firstGene,secondGene)<=0)
		{
			sprintf(keyChar,"%s-%s",firstGene,secondGene);
		}
		else
		{
			sprintf(keyChar,"%s-%s",secondGene,firstGene);
		}
	}
	else
	{
		sprintf(keyChar,"%s-%s",firstGene,secondGene);
	}
	aKey.append(keyChar);
	
	return 0;
}


