#ifndef _GENEEXPMGR_
#define _GENEEXPMGR_

#include <map>
#include <string>
#include <vector>

using namespace std;

class GeneExpManager
{
	public:
	GeneExpManager();
	~GeneExpManager();
	int readExpression(const char* aFName);
	vector<double>* getExp(const string&);
	map<string,vector<double>*>& getGeneSet();
	
	private:
	map<string,vector<double>*> exprSet;
};
#endif
