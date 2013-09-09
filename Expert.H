#ifndef _EXPERT_
#define _EXPERT_
#include <vector>
using namespace std;
#include "gsl/gsl_randist.h"
#include "CommonTypes.H"
#define PI 3.1472
class Matrix;
class Expert
{
	public:
		Expert();
		~Expert();
		int setMean(Matrix*);
		int setCovariance(Matrix*);
		double getOutputPDF(vector<double>* y);
		double getOutputPDF_Nocov(vector<double>* y);
		int setPrior(double);
		double getPrior();
		Matrix* getMean();
		Matrix* getCovariance();
		int updateCovariance();
		int assignGeneToExpert(const char*);
		map<string,int>& getGeneSet();
		int resetAssignedGenes();
		double getEntropy();
		int generateSample(gsl_rng* r, vector<double>&);
		int resetClip();
		int getClip();
		inline double getMinPDF()
		{
			return minpdf;
		}
	private:
		 Matrix* covariance;
		 Matrix* invCovariance;
		 Matrix* mean;
		 double normFactor;
		 map<string,int> geneSet;
		double priorProb;
		int clipCnt;
		double minpdf;
};
#endif

