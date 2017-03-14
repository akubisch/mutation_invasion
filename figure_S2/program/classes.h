const int XMAX  = 250;     //max. x-dimension of the world
const int YMAX = 20;     //max. y-dimension of the world
const float VARIANCE = 0.2;   //variance used for many calculations
const float VARIANCE_ADAP = 0.5;

const int MAXF = 1000;
const int MAXM = 1000;
const int DELTA_ZERO_GROWTH = 10;       //individual fitness is decreased to zero, when
//the mismatch between trait and environment
//exceeds this value
const int TIME_OUT_INT = 100;
const int INIT_AREA = 50;
const int BURN_IN_TIME = 10000;
const float PROB_DEATH = 0.9;

class tInd
{
public:
								double ld[2]; //one locus with 2 dispersal alleles
								double la[2]; //one locus with 2 adaptation alleles
								double le[2]; //one locus with 2 evolvability alleles
								double ln[2]; //one locus with 2 neutral alleles
};

class tPatch
{
public:
								tPatch();
								vector<tInd> males,females; //population
								vector<tInd> malimmis,femimmis; //immigrants
								float temp; //temperature in this patch
};

//Konstruktor für tCell

tPatch::tPatch() {
								temp = 0;
								males.clear();
								males.reserve(MAXM);
								females.clear();
								females.reserve(MAXF);
								malimmis.clear();
								malimmis.reserve(MAXM);
								femimmis.clear();
								malimmis.reserve(MAXF);
}

class tParameters
{
public:
								int capacity; //habitat capacity
								double lambda_0; //per capita growth rate
								double sigma; //environmental fluctuations
								double disp_mort;
								double ext_prob; //extinction rate
								double mut_prob; //mutation probability (for dispersal)
								int tmax; //max. no. of generations to be simulated
								bool ld_evol,la_evol; //evolution of ld/la?

								double niche; //niche widths of the species

								bool ddd; //density-dependent dispersal yes/no

								double temp_max; //maximum temperature (x_max = temp_grad) == steepness of gradient (x0 = 0)

								int replications; //number of replicates per experiment
};

class tAnalysis
{
public:
								double dens[XMAX][YMAX]; //population density of patches
								double sum_adap_x[XMAX];
								unsigned long cnt_adap_x[XMAX];
};
