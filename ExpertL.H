#ifndef _EXPERTL_
#define _EXPERTL_
#include <vector>
using namespace std;
#include "CommonTypes.H"
#define PI 3.1472
class Matrix;
class ExpertL
{
	public:
		ExpertL();
		~ExpertL();
		int setMean(Matrix*);
		int setCovariance(Matrix*);
		int setFeatureWeight(int,double);
		double getOutputPDF(Matrix* y);
		double getOutputPDF_Fast(Matrix* y);
		double getOutputPDF_Nocov(Matrix* y);
		double getMixtureWeight(INTDBLMAP* x);
		Matrix* getMean();
		Matrix* getCovariance();
		int updateCovariance();
		INTDBLMAP& getFeatureWeights();
		int sortFeatures();
		int sortFeatures_Enrichment(INTDBLMAP&);
		vector<int>& getSortedFeatures();
		int assignGeneToExpert(int);
		map<int,int>& getGeneSet();
		int resetAssignedGenes();
		double getEntropy();
		int setPrior(double);
		double getPrior();
	private:
		 Matrix* covariance;
		 Matrix* invCovariance;
		 Matrix* mean;
		 double normFactor;
		 INTDBLMAP featureWeights;
		 vector<int> sortedFeatureIDs;
		 map<int,int> geneSet;
		double prior;
};
#endif

