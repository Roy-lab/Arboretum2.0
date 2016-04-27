#ifndef _BIONODE
#define _BIONODE
#include <vector>
#include <map>
#include <string>

using namespace std;

#include "CommonDefs.H"
typedef vector<double> EXPRLEVELS;

class Path;
class BioNode;

typedef map<string,BioNode*> BNODE_MAP;
typedef map<string,BioNode*>::iterator BNODE_MAP_ITER;

class BioNode
{
	public:
		BioNode();
		~BioNode();
		int setName(const char*);
		int setID(int);

		const char* getName();
		int getID();
		
		typedef enum
		{
			PROTEIN,
			GENE
		}BioNodeType;

		int setType(BioNodeType);
		BioNodeType getType();

		int addExpLevel(double);
		double getExpLevelAt(int);
		double getRankedExpLevelAt(int);
	
		typedef enum
		{
			LOW,
			MEDIUM,
			HIGH
		}BioNodeLevel;

		typedef enum
		{
			IN_NEIGHBOUR,
			OUT_NEIGHBOUR
		}NeighbourType;

		
		int discretizeLevels();
		int rankOrder();
		int rankOrder_Store();
		int discretizeEqualBins();
		int discretizeEqualBins_Global(double,double);
		int scaleExpression();
		int standardizeData();
		int zeroMeanData();
		int getDiscreteExpLevelAt(int);
		int addNeighbour(BioNode*,NeighbourType);
		BNODE_MAP& getOutNeighbours();
		BNODE_MAP& getInNeighbours();
			
		int findReachableNodes();
		int getReachableNodeCnt();
		map<string,Path*>& getPaths();
		int dumpPaths();
		int showInNeighbours();
		int showPathTo(const char*);
		int calculateMean();
		int calculateStdev(double);
		double getMean();
		double getStdev();
		double getMaxMin(double&,double&);	
	protected:
		int buildReachableMap();
		int buildNodePathAssoc();
		
		int sort(INTVECT&);

		char name[GENENAME_LEN];
		int nId;
		EXPRLEVELS expLevels;
		INTVECT discreteExpLevels;
		EXPRLEVELS rankedExpLevels;
		
		BNODE_MAP inNeighbours;
		BNODE_MAP outNeighbours;
		STRINTMAP reachableNeighbours;

		//The set of reachable nodes from this node using a particular path length.
		//The reachable nodes maybe a proteins or genes. E.g. this gene will be connected
		//to another gene, if it encodes a tf which acts on this gene
		//This gene will be connected to a protein if there is a protein which is not a tf.
		map<int,STRVECT*> neighbourSets; 

		//The number of paths will be equal to the number of reachable neighbours
		map<string,Path*> pathSet;
		map<string,string> nodePathAssoc;
		BioNodeType nodeType;
		double expMin;
		double mean;
		double stdev;
};
#endif
