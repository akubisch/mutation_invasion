  //////////////////////////////////////////////////////////////////////
//																		//
//		       		  SOME BASIC PROCEDURES (GSL & co.)					//
//						-----------------------------					//
//								Alex Kubisch							//
//						  Evolutionary Ecology Group					//
//			  			   University of Wuerzburg						//
//																		//
//	    	  	  licensed under GNU GPL Version 2, June 1991			//
//																		//
  //////////////////////////////////////////////////////////////////////

const gsl_rng *gBaseRand;

//________________________________________________________________________________________
//------------------------------------------------------Initialize Random Number Generator

void specify_rng(unsigned long randSeed)
{
	gBaseRand = gsl_rng_alloc(gsl_rng_rand);
	
	srand(randSeed);
	unsigned long r = random();
	gsl_rng_set(gBaseRand, r);
}

//________________________________________________________________________________________
//-------------------------------------------------------------------------Simplifications

//-------------------------------------------------Simplify Random Drawing between 0 and 1

double ran()
{
	return gsl_rng_uniform(gBaseRand);
}

//---------------------------------------------------------------Simplify Gaussian Randoms

double gauss(double sd)
{
	return gsl_ran_gaussian(gBaseRand,sd);
}

//------------------------------------------------------------------Simplify Cauchy Randoms

double cauchy(double d)
{
	return gsl_ran_cauchy(gBaseRand,d);
}

//-----------------------------------------------------------------Simplify Poisson Random

int poisson(double sd)
{
	return gsl_ran_poisson(gBaseRand,sd);
}

//---------------------------------------------------------------Simplify Lognormal Random

double lognorm(double zeta, double sigma)
{
	double var;											//variance of resulting randoms
	double mu,s_d;										//mean and sd of lognormal distr.
														//to be calculated by mean and
														//sigma of resulting randoms
	var = sigma*sigma;
	s_d = sqrt(log((var/(2*zeta))+1));
	mu = log(zeta)-0.5*(s_d*s_d);
	return gsl_ran_lognormal(gBaseRand,mu,s_d);
}

//---------------------------------------------------------------Simplify mean calculation

double mean(double data[], int n)
{
	return gsl_stats_mean(data,1,n);
}

//--------------------------------------------------------Simplify vector mean calculation

double mean(vector<double> vec)
{
	int n = int(vec.size());
	double data[n];
	for(int i = 0; i<n; i++) {
		data[i] = double(vec.at(i));
	}
	return gsl_stats_mean(data,1,n);
}

//-----------------------------------------------------Simplify integer median calculation

double median(vector<int> vec, int n)
{
	int data[n];
	for(int i = 0; i<n; i++) {
		data[i] = vec.at(i);
	}
	gsl_sort_int(data,1,n);
	double medi = gsl_stats_int_median_from_sorted_data(data,1,n);
	return(medi);
}

//---------------------------------------------------Simplify integer quantile calculation

double q005(vector<int> vec, int n)
{
	int data[n];
	for(int i = 0; i<n; i++) {
		data[i] = vec.at(i);
	}
	gsl_sort_int(data,1,n);
	double q = gsl_stats_int_quantile_from_sorted_data(data,1,n,0.05);
	return(q);
}

double q095(vector<int> vec, int n)
{
	int data[n];
	for(int i = 0; i<n; i++) {
		data[i] = vec.at(i);
	}
	gsl_sort_int(data,1,n);
	double q = gsl_stats_int_quantile_from_sorted_data(data,1,n,0.95);
	return(q);
}

//----------------------------------------------------Simplify minimum and maximum finding

double min_int(vector<int> vec, int n)
{
	double data[n];
	for(int i = 0; i<n; i++) {
		data[i] = double(vec.at(i));
	}
	return(round(gsl_stats_min(data,1,n)));
}

double max_int(vector<int> vec, int n)
{
	double data[n];
	for(int i = 0; i<n; i++) {
		data[i] = double(vec.at(i));
	}
	return(round(gsl_stats_max(data,1,n)));
}

//-------------------------------------------------------------Simplify standard deviation

double sd(double data[], int n)
{
	return gsl_stats_sd(data,1,n);
}

double var(vector<double> vec)
{
	int n = int(vec.size());
	double data[n];
	for(int i = 0; i<n; i++) {
		data[i] = double(vec.at(i));
	}
	return gsl_stats_variance(data,1,n);
}

//---------------------------------------------------------------------Round double to int

int rnd(double a) 
{
	return round(a);
}
