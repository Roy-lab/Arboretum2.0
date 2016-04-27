#include <iostream>
#include <string.h>
#include "Path.H"
#include "Gene.H"

Gene::Gene()
{
}

Gene::~Gene()
{
}

int 
Gene::setCodedProteinName(const char* protName)
{
	strcpy(assocProtName,protName);
	return 0;
}

const char*
Gene::getCodedProteinName()
{
	return assocProtName;
}
