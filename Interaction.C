#include <string.h>
#include "Interaction.H"

Interaction::Interaction()
{
	occurCnt=1;
}

Interaction::~Interaction()
{
}

int 
Interaction::setMembers(const char* mem1, const char* mem2)
{
	strcpy(firstMember,mem1);
	strcpy(secondMember,mem2);
	return 0;
}

const char*
Interaction::getFirstMember()
{
	return firstMember;
}

const char*
Interaction::getSecondMember()
{
	return secondMember;
}

int 
Interaction::setProperty(const char* name, const char* value)
{
	string nameStr(name);
	string valueStr(value);
	properties[nameStr]=valueStr;
	return 0;
}


const char* 
Interaction::getProperty(const char* name)
{
	string nameStr(name);
	map<string,string>::iterator aIter=properties.find(nameStr);
	if(aIter==properties.end())
	{
		return NULL;
	}
	return aIter->second.c_str();
}

