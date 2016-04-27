#include <iostream>
#include <fstream>
using namespace std;

#include "gsl/gsl_randist.h"
#include "Error.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "BFGSWrapperData.H"
#include "BFGSWrapper.H"


BFGSWrapper::BFGSWrapper()
{
	initx=NULL;
	optimizer=NULL;
	optimizer_func=NULL;
}

BFGSWrapper::~BFGSWrapper()
{
	if(initx!=NULL)
	{
		gsl_vector_free(initx);
	}
	if(optimizer!=NULL)
	{
		gsl_multimin_fdfminimizer_free(optimizer);
	}
}

int
BFGSWrapper::setParamCnt(int aCnt)
{
	paramCnt=aCnt;
	return 0;
}


int
BFGSWrapper::setFeatureCnt(int aCnt)
{
	featureCnt=aCnt;
	return 0;
}

int
BFGSWrapper::setStepSize(double asize)
{
	stepSize=asize;
	return 0;
}

int
BFGSWrapper::setTolerance(double tol)
{
	tolerance=tol;
	return 0;
}

gsl_vector*
BFGSWrapper::getInitVector()
{
	return initx;
}

int
BFGSWrapper::setInitVector(gsl_vector* x)
{
	initx=x;
	return 0;
}
	
map<int,map<int,double>*>& 
BFGSWrapper::getData()
{
	return data;
}

map<int,map<int,double>*>& 
BFGSWrapper::getGammas()
{
	return gammas;
}

int 
BFGSWrapper::initializeMinimizer()
{
	bwdata.data=&data;
	bwdata.gammas=&gammas;
	const gsl_multimin_fdfminimizer_type* T=gsl_multimin_fdfminimizer_vector_bfgs;
	//const gsl_multimin_fdfminimizer_type* T=gsl_multimin_fdfminimizer_conjugate_fr;
	optimizer=gsl_multimin_fdfminimizer_alloc(T,paramCnt);
	optimizer_func=new gsl_multimin_function_fdf;
	optimizer_func->n=paramCnt;
	optimizer_func->f=&likelihood_function;
	optimizer_func->df=&likelihood_gradient;
	optimizer_func->fdf=&likelihood_function_gradient;
	optimizer_func->params=(void*)(&bwdata);
	//gsl_vector* initx=getStartingPoint();
	gsl_multimin_fdfminimizer_set(optimizer,optimizer_func,initx,stepSize,tolerance);
	return 0;
}

int
BFGSWrapper::reinitializeMinimizer()
{
	gsl_multimin_fdfminimizer_restart(optimizer);
	//gsl_multimin_fdfminimizer_set(optimizer,optimizer_func,optimizer->x,stepSize,tolerance);
	return 0;
}

int
BFGSWrapper::optimize()
{
	int iter=0;
	int status;
	do
	{
		iter++;
		status=gsl_multimin_fdfminimizer_iterate(optimizer);
		cout <<"Curr Multimin_Iter " << iter << " function_val " << optimizer->f << endl;
		if(status)
		{
			cout <<"Breaking: gsl_errno " << status << endl;
			status=gsl_multimin_test_gradient(optimizer->gradient,1e-3);
			break;
		}
		status=gsl_multimin_test_gradient(optimizer->gradient,1e-3);
		if(status==GSL_SUCCESS)
		{
			cout <<"Minima found "<< iter << endl;
			gsl_vector* wt=optimizer->x;
			for(int i=0;i<paramCnt;i++)
			{
		       		cout << i << " "<< gsl_vector_get(wt,i) << endl;
			}
			cout <<"function val " << optimizer->f << endl;
		}
	} while(status==GSL_CONTINUE && iter<100);
	return 0;
}


gsl_vector* 
BFGSWrapper::getParams()
{
	gsl_vector* params=optimizer->x;
	return params;
}

double 
BFGSWrapper::getOptimalFval()
{
	double f=optimizer->f;
	return f;
}

int
BFGSWrapper::reset()
{
	gammas.clear();
	for(map<int,map<int,double>*>::iterator dIter=data.begin();dIter!=data.end();dIter++)
	{
		dIter->second->clear();
		delete dIter->second;
	}
	data.clear();
	if(initx!=NULL)
	{
		gsl_vector_free(initx);
	}
	if(optimizer!=NULL)
	{
		gsl_multimin_fdfminimizer_free(optimizer);
	}
	return 0;
}

double 
likelihood_function(const gsl_vector* x, void* params)
{
	//cout <<"In likelihood function " << endl;
	//Compute the softmax probability
	double ans=0;
	double den=0;
	//The vector looks as follows
	//The first k dimensions have the gammas
	//The remaining dimensions have the motifs
	double gvals[MAXCOMP];
	BFGSWrapperData bwdata=*((BFGSWrapperData*) params);
	map<int,map<int,double>*>& data=*(bwdata.data);
	map<int,map<int,double>*>& gammas=*(bwdata.gammas);
	//cout <<"Analyzing data of size " << data.size() << " with gammas " << gammas.size() << endl;
	double lambda=1.5;
	double priorcontr=0;
	for(map<int,map<int,double>*>::iterator dIter=data.begin();dIter!=data.end();dIter++)
	{
		map<int,double>* dataPt=dIter->second;
		map<int,double>* gamma_dataPt=gammas[dIter->first];
		if(gammas.find(dIter->first)==gammas.end())
		{
			cout <<"No gammas for " << dIter->first << endl;
			exit(0);
		}
		if(gamma_dataPt==NULL)
		{
			cout <<"Gamma null for " << dIter->first << endl;
			exit(0);
		}
		int featureCnt=dataPt->size();
		//cout <<"Analyzing datapoint " << dIter->first << endl;
		for(map<int,double>::iterator gIter=gamma_dataPt->begin();gIter!=gamma_dataPt->end();gIter++)
		{
			double sum=0;
			for(map<int,double>::iterator fIter=dataPt->begin();fIter!=dataPt->end();fIter++)
			{
				double val=gsl_vector_get(x,(gIter->first*featureCnt)+fIter->first);
				if(dIter==data.begin())
				{
					double wt=gsl_vector_get(x,(gIter->first*featureCnt)+fIter->first);
					priorcontr=priorcontr+(lambda*wt);
				}
				double dataPtVal=fIter->second;
				sum=sum+(val*dataPtVal);
				if(isnan(val))
				{
				//	cout <<"val is nan " << val << " for mixture " << gIter->first << " featureid " << fIter->first << endl;
				}
			}	
			gvals[gIter->first]=exp(sum);
			den=den+gvals[gIter->first];
		}
		for(map<int,double>::iterator gIter=gamma_dataPt->begin();gIter!=gamma_dataPt->end();gIter++)
		{
			double gamma_i=gIter->second;
			if(dIter==data.begin())
			{
	//			cout <<" Gamma_"  << gIter->first << "="<< gamma_i;
			}
			double g_i=gvals[gIter->first]/den;
			ans=ans+(gamma_i*log(g_i));
			if(isnan(ans))
			{
			//	cout <<"ans is nan gamma_i: " << gamma_i<< " wt_i:" << g_i 
				//<< " expert_id:" <<  gIter->first<< " datapt:" << dIter->first<< endl;
			}
		}
		if(dIter==data.begin())
		{
	//		cout << endl;
		}
	}
	ans=ans-priorcontr;
	//cout <<"Current fval " << -1*ans << endl;
	return -1.0*ans;
	//return ans;
}

//Store the gradient of the function in g using the value of x
/*void 
likelihood_gradient(const gsl_vector* x, void* params,gsl_vector* g)
{
	//cout <<"In likelihood_gradient function " << endl;
	//Iterate over the complete data to get the common parts first
	map<int,map<int,double>*> sum1;
	map<int,map<int,double>*> sum2;
	BFGSWrapperData bwdata=*((BFGSWrapperData*) params);
	map<int,map<int,double>*>& data=*(bwdata.data);
	map<int,map<int,double>*>& gammas=*(bwdata.gammas);

	for(map<int,map<int,double>*>::iterator dIter=data.begin();dIter!=data.end();dIter++)
	{
		map<int,double>* dataPt=dIter->second;
		map<int,double>* gamma_dataPt=gammas[dIter->first];	
		int featureCnt=dataPt->size();
		//Data point specific sum
		map<int,double> dataPt_sum;
		double dataPt_sum_experts=0;
		for(map<int,double>::iterator gIter=gamma_dataPt->begin();gIter!=gamma_dataPt->end();gIter++)
		{
			int expert_i=gIter->first;
			double gval_k=gIter->second;
			map<int,double>* sum1_expert_i=NULL;
			if(sum1.find(gIter->first)==sum1.end())
			{
				sum1_expert_i=new map<int,double>;
				sum1[gIter->first]=sum1_expert_i;
			}
			else
			{
				sum1_expert_i=sum1[gIter->first];
			}
			double sum=0;
			for(int f=0;f<featureCnt;f++)
			{
				double wx=gsl_vector_get(x,(expert_i*featureCnt)+f);
				double dval=(*dataPt)[f];
				if(sum1_expert_i->find(f)==sum1_expert_i->end())
				{
					//(*sum1_expert_i)[f]=dval*gval_k;
					(*sum1_expert_i)[f]=dval;
				}
				else
				{
					//(*sum1_expert_i)[f]=(*sum1_expert_i)[f]+(dval*gval_k);
					(*sum1_expert_i)[f]=(*sum1_expert_i)[f]+(dval);
				}
				sum=sum+(wx*dval);
			}
			double xi=exp(sum);
			dataPt_sum[expert_i]=xi;
			dataPt_sum_experts=dataPt_sum_experts+xi;
		}
		double overall_sum=0;
		for(map<int,double>::iterator gIter=gamma_dataPt->begin();gIter!=gamma_dataPt->end();gIter++)
		{
			int expert_i=gIter->first;
			double gval_k=gIter->second;
			dataPt_sum[expert_i]=dataPt_sum[expert_i]/dataPt_sum_experts;
			overall_sum=overall_sum+(gval_k*dataPt_sum[expert_i]);
		}
		dataPt_sum.clear();
		for(map<int,double>::iterator gIter=gamma_dataPt->begin();gIter!=gamma_dataPt->end();gIter++)
		{
			int expert_i=gIter->first;
			map<int,double>* sum2_expert_i=NULL;
			if(sum2.find(gIter->first)==sum2.end())
			{
				sum2_expert_i=new map<int,double>;
				sum2[gIter->first]=sum2_expert_i;
			}
			else
			{
				sum2_expert_i=sum2[gIter->first];
			}
			for(int f=0;f<featureCnt;f++)
			{
				double dval=(*dataPt)[f];
				double expdval=overall_sum*dval;
				
				if(sum2_expert_i->find(f)==sum2_expert_i->end())
				{
					(*sum2_expert_i)[f]=expdval;
				}
				else
				{
					(*sum2_expert_i)[f]=(*sum2_expert_i)[f]+expdval;
				}
			}
			
		}
	}
	int paramIter=0;
	for(map<int,map<int,double>*>::iterator gIter=sum1.begin();gIter!=sum1.end();gIter++)
	{
		map<int,double>* sum1_expert=gIter->second;
		map<int,double>* sum2_expert=sum2[gIter->first];
		for(map<int,double>::iterator fIter=sum1_expert->begin();fIter!=sum1_expert->end();fIter++)
		{
			double val=(*sum2_expert)[fIter->first]-fIter->second;
			//double val=fIter->second-(*sum2_expert)[fIter->first];
			gsl_vector_set(g,paramIter,val);
			if(paramIter<10)
			{
			//	cout << "Param " << paramIter<< " before " << gsl_vector_get(x,paramIter) << " after " << val << endl;
			}
			paramIter++;
		}
	}
	for(map<int,map<int,double>*>::iterator gIter=sum1.begin();gIter!=sum1.end();gIter++)
	{
		map<int,double>* sum1_expert=gIter->second;
		map<int,double>* sum2_expert=sum2[gIter->first];
		sum1_expert->clear();
		sum2_expert->clear();
		delete sum1_expert;
		delete sum2_expert;
	}
	sum1.clear();
	sum2.clear();
		
}*/


//Store the gradient of the function in g using the value of x
void 
likelihood_gradient(const gsl_vector* x, void* params,gsl_vector* g)
{
	//cout <<"In likelihood_gradient function " << endl;
	//Iterate over the complete data to get the common parts first
	map<int,map<int,double>*> sum1;
	map<int,map<int,double>*> sum2;
	BFGSWrapperData bwdata=*((BFGSWrapperData*) params);
	map<int,map<int,double>*>& data=*(bwdata.data);
	map<int,map<int,double>*>& gammas=*(bwdata.gammas);

	for(map<int,map<int,double>*>::iterator dIter=data.begin();dIter!=data.end();dIter++)
	{
		map<int,double>* dataPt=dIter->second;
		map<int,double>* gamma_dataPt=gammas[dIter->first];	
		//Data point specific sum
		map<int,double> dataPt_sum;
		double dataPt_sum_experts=0;
		int featureCnt=dataPt->size();
		for(map<int,double>::iterator gIter=gamma_dataPt->begin();gIter!=gamma_dataPt->end();gIter++)
		{
			int expert_i=gIter->first;
			double sum=0;
			for(int f=0;f<featureCnt;f++)
			{
				double wx=gsl_vector_get(x,(expert_i*featureCnt)+f);
				double dval=(*dataPt)[f];
				sum=sum+(wx*dval);
			}
			double xi=exp(sum);
			double gval_di=gIter->second;
			dataPt_sum[expert_i]=xi*gval_di;
			dataPt_sum_experts=dataPt_sum_experts+xi;
		}
		double overall_sum=0;
		for(map<int,double>::iterator gIter=gamma_dataPt->begin();gIter!=gamma_dataPt->end();gIter++)
		{
			int expert_i=gIter->first;
			dataPt_sum[expert_i]=dataPt_sum[expert_i]/dataPt_sum_experts;
			overall_sum=overall_sum+dataPt_sum[expert_i];
		}
		dataPt_sum.clear();
		for(map<int,double>::iterator gIter=gamma_dataPt->begin();gIter!=gamma_dataPt->end();gIter++)
		{
			int expert_i=gIter->first;
			map<int,double>* sum2_expert_i=NULL;
			if(sum2.find(gIter->first)==sum2.end())
			{
				sum2_expert_i=new map<int,double>;
				sum2[gIter->first]=sum2_expert_i;
			}
			else
			{
				sum2_expert_i=sum2[gIter->first];
			}
			double gval=gIter->second;
			for(int f=0;f<featureCnt;f++)
			{
				double dval=(*dataPt)[f];
				double expdval=(overall_sum-gval)*dval;
				
				if(sum2_expert_i->find(f)==sum2_expert_i->end())
				{
					(*sum2_expert_i)[f]=expdval;
				}
				else
				{
					(*sum2_expert_i)[f]=(*sum2_expert_i)[f]+expdval;
				}
			}
			
		}
	}
	int paramIter=0;
	double lambda=1.5;
	for(map<int,map<int,double>*>::iterator gIter=sum2.begin();gIter!=sum2.end();gIter++)
	{
		map<int,double>* sum2_expert=sum2[gIter->first];
		int expert_i=gIter->first;
		int featureCnt=sum2_expert->size();
		for(map<int,double>::iterator fIter=sum2_expert->begin();fIter!=sum2_expert->end();fIter++)
		{
			double val=fIter->second;
			double wx=gsl_vector_get(x,(expert_i*featureCnt)+fIter->first);
			val=(lambda*(wx/(fabs(wx))))-val;
			//double val=fIter->second-(*sum2_expert)[fIter->first];
			gsl_vector_set(g,paramIter,val);
			if(paramIter<10)
			{
			//	cout << "Param " << paramIter<< " before " << gsl_vector_get(x,paramIter) << " after " << val << endl;
			}
			paramIter++;
		}
	}
	for(map<int,map<int,double>*>::iterator gIter=sum1.begin();gIter!=sum1.end();gIter++)
	{
		map<int,double>* sum1_expert=gIter->second;
		map<int,double>* sum2_expert=sum2[gIter->first];
		sum1_expert->clear();
		sum2_expert->clear();
		delete sum1_expert;
		delete sum2_expert;
	}
	sum1.clear();
	sum2.clear();
		
}


void
likelihood_function_gradient(const gsl_vector* x, void* params, double* f, gsl_vector* g)
{
	*f=likelihood_function(x,params);
	likelihood_gradient(x,params,g);
}



//Get the starting point for the parameters that we are optimizing
gsl_vector* 
BFGSWrapper::getStartingPoint()
{
	gsl_vector* initx=gsl_vector_alloc(paramCnt);
	gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
	for(int p=0;p<paramCnt;p++)
	{
		double w_p=gsl_ran_flat(r,0,1);
		gsl_vector_set(initx,p,w_p);
	}
	return initx;
}
