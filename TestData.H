#ifndef _TEST_DATA_
#define _TEST_DATA_

#include <map>
#include <string>
using namespace std;
class TestData
{
	public:
	TestData();
	~TestData();
	int readData(const char*);
	int getSize();
	map<string,map<string,int>*>& getDataSet();
	private:
	map<string,map<string,int>*> dataSet;
};

#endif
