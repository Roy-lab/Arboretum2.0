#ifndef _INTERACTIONMANAGER
#define _INTERACTIONMANAGER

#include <map>
#include <string>
#include "Interaction.H"

using namespace std;

class Interaction;
typedef map<string, Interaction*> INTR_MAP;
typedef map<string, Interaction*>::iterator INTR_MAP_ITER;

class InteractionManager 
{
	public:
		InteractionManager();
		~InteractionManager();
		int setDirectionality(bool);
		int setDelimiter(const char*);
		int readInteractions(const char*);
		INTR_MAP& getInteractions();
		
		bool isInteraction(const char*, const char*);
		int addInteraction(const char*,const char*);
		int getNumberOfIntr();
		Interaction* getInteractionAt(int);
		int dumpInteractions(const char*);
	private:
		int createKey(const char*, const char* ,string&);
		INTR_MAP interactionSet;
		map<int,string> intrIDtoString;
		bool directionality;
		char delim[5];
};
#endif
