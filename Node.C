#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "Node.H"

Node::Node()
{
	valueStr[0]='\0';
}

Node::~Node()
{
}

int 
Node::setName(const char* aName)
{
	strcpy(nodeName,aName);
	return 0;
}

const char* 
Node::getName()
{
	return nodeName;
}

int 
Node::setNodeType(Node::NodeType aType)
{
	ntype=aType;
	return 0;
}

int 
Node::setIdInData(int id)
{
	nodeId=id;
	return 0;
}

int 
Node::getIdInData()
{
	return nodeId;
}

int 
Node::setTopSortId(int id)
{
	topsortId=id;
}

int 
Node::getTopSortId()
{
	return topsortId;
}


Node::NodeType 
Node::getNodeType()
{
	return ntype;
}

int 
Node::setValues(vector<int>& valSet)
{
	for(int i=0;i<valSet.size();i++)
	{
		values.push_back(valSet[i]);
		if(i==0)
		{
			sprintf(valueStr,"%d",valSet[i]);
		}
		else
		{
			char tempStr[5];
			sprintf(tempStr,",%d",valSet[i]);
			strcat(valueStr,tempStr);
		}
	}
	return 0;
}

vector<int>& 
Node::getValues()
{
	return values;
}

const char*
Node::getValueStr()
{
	return valueStr;
}

int 
Node::addChild(Node* aNode)
{
	children.push_back(aNode);
	return 0;
}

int 
Node::addParent(Node* aNode)
{
	parents.push_back(aNode);
	return 0;
}

vector<Node*>&
Node::getChildren()
{
	return children;
}
	
vector<Node*>&
Node::getParents()
{
	return parents;
}

int
Node::getNeighbourCnt()
{
	return parents.size()+children.size();
}
