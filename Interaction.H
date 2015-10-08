#ifndef _INTERACTIONS
#define _INTERACTIONS

#include <map>
#include <string>

using namespace std;

#include "CommonDefs.H"

class Interaction
{
	public:
		Interaction();
		~Interaction();
		int setMembers(const char*, const char*);
		const char* getFirstMember();
		const char* getSecondMember();
		int setProperty(const char*, const char*);
		const char* getProperty(const char*);
	private:	
		map<string,string> properties;
		char firstMember[GENENAME_LEN];
		char secondMember[GENENAME_LEN];
		int occurCnt;
};
#endif
