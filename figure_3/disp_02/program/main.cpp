#include <iostream>
#include <cstdlib>                                                      //standard C library
#include <ctime>                                                        //access system time library
#include <fstream>                                                      //file streaming library
#include <string>                                                       //string library included
#include <sstream>                                                      //string streaming for reading numbers from
//infiles
#include <vector>
#include <cmath>                                                        //standard math library

#include <gsl/gsl_rng.h>                                        //gsl random number generator
#include <gsl/gsl_randist.h>                            //gsl random distributions
#include <gsl/gsl_statistics.h>                         //gsl some statistical methods
#include <gsl/gsl_statistics_int.h>                     //gsl some integer statistical methods
#include <gsl/gsl_sort_int.h>                           //gsl sort integer arrays
#include <gsl/gsl_math.h>                                       //provides additional standard math functions

using namespace std;

#include "procedures.h"                                         //procedure simplifications
#include "classes.h"                                            //class definitions

//Variables

vector<vector<tPatch> >                 world;                          //simulated landscape
tParameters par;                                                        //parameters
tAnalysis ana;
double maxtemp=0;
double mintemp=1000;
double initial_range;

int real_border,pred_border;                                    //realized and predicted interspe. range border

fstream parinfile("data/para.in", ios::in);                     //parameter infile
fstream adapoutfile("data/adaptation.txt", ios::out);
fstream dispoutfile("data/dispersal.txt", ios::out);
fstream evoloutfile("data/evolvability.txt", ios::out);
fstream densoutfile("data/densities.txt", ios::out);
fstream occoutfile("data/occupancy.txt", ios::out);
fstream divoutfile("data/diversity.txt", ios::out);
fstream neutdivoutfile("data/neutral_diversity.txt", ios::out);
fstream mutdivoutfile("data/mut_diversity.txt", ios::out);
fstream dispdivoutfile("data/disp_diversity.txt", ios::out);

///////////////////////////////////////////////////////////////////////////////////////
//---------------------------------SET PARAMETERS------------------------------------//
///////////////////////////////////////////////////////////////////////////////////////

void set_parameters() {                                 //read in standard simulation parameters

        string buffer;
        istringstream is;

        getline(parinfile,buffer);
        getline(parinfile,buffer);
        if (buffer=="yes") {par.ddd = true; } else {par.ddd = false; }; //density-dependent emigration?

        getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> par.niche; //species' niche widths

        getline(parinfile,buffer);
        getline(parinfile,buffer);
        if (buffer=="1") {par.ld_evol = true; } else {par.ld_evol = false; } //evolving dispersal?
        getline(parinfile,buffer);
        getline(parinfile,buffer);
        if (buffer=="1") {par.la_evol = true; } else {par.la_evol = false; } //evolving adaptation?
        getline(parinfile,buffer);
        getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> par.capacity;                       //habitat capacity
        getline(parinfile,buffer);
        getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> par.lambda_0;                       //per capita growth rate
        getline(parinfile,buffer);
        getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> par.sigma;                          //environmental stochasticity
        getline(parinfile,buffer);
        getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> par.disp_mort;                          //environmental stochasticity
        getline(parinfile,buffer);
        getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> par.ext_prob;                       //extinction probability
        getline(parinfile,buffer);
        getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> par.mut_prob;                       //mutation probability
        getline(parinfile,buffer);
        getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> par.temp_max;                       //max. temperature in gradient (≃ steepness)
        getline(parinfile,buffer);
        getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> par.tmax;                                   //simulated time
        getline(parinfile,buffer);
        getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> par.replications;                   //no of replications
}

///////////////////////////////////////////////////////////////////////////////////////
//-----------------------------INITIALIZATION FUNCTIONS------------------------------//
///////////////////////////////////////////////////////////////////////////////////////

void apply_gradient() {
        for (int x = 0; x < XMAX; x++)
                for (int y = 0; y < YMAX; y++) {
                        world[x][y].temp = par.temp_max-(par.temp_max/double(XMAX)) * double(x);
                }
}

void init_individual(tInd *newind, int x, int y) {
        newind->la[0] = world[x][y].temp + gauss(VARIANCE);
        newind->la[1] = world[x][y].temp + gauss(VARIANCE);
        newind->ld[0] = ran();
        newind->ld[1] = ran();
        newind->le[0] = 4 + gauss(1);
        newind->le[1] = 4 + gauss(1);
        newind->ln[0] = newind->la[0];
        newind->ln[1] = newind->la[1];
}

void init_world() {

        tInd newind;

        world.resize(XMAX);
        for (int x=0; x < XMAX; x++) {
                world[x].resize(YMAX);
        }

        apply_gradient();

        for (int x = 0; x < INIT_AREA; x++)
                for (int y = 0; y < YMAX; y++) {

                        for (int i = 0; i < (par.capacity/2); i++) {
                                init_individual(&newind,x,y);
                                world[x][y].females.push_back(newind);
                        }
                        for (int i = 0; i < (par.capacity/2); i++) {
                                init_individual(&newind,x,y);
                                world[x][y].males.push_back(newind);
                        }
                }

}

///////////////////////////////////////////////////////////////////////////////////////
//-------------------------------SIMULATION FUNCTIONS--------------------------------//
///////////////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________________
//-Dispersal---------------------------------------------------------------------------

double disp_prob(int x, int y, tInd ind) {

        return 0.2;
}

void find_target_patch(int xs, int ys, int t, int *xt, int *yt, bool burn) {
        double dist;                      //drawn dispersal distance
        double phi;

        dist = 1;                                 //gives the possibility to implement a kernel

        phi = 2 * M_PI * ran();           //direction is uniformly drawn (0°-360°)

        double xstart = double(xs) + ran() - 0.5; //simulating area to area dispersal				<<== IF KERNEL IS USED!!!
        double ystart = double(ys) + ran() - 0.5; //(avoid artefacts by drawing starting position <<== IF KERNEL IS USED!!!
        //for dispersal randomly inside the patch)              <<== IF KERNEL IS USED!!!

        int xtt = rnd(xstart + dist*cos(phi));
        int ytt = rnd(ystart + dist*sin(phi));

        if (burn) {
                if (xtt>INIT_AREA) {xtt = INIT_AREA-1; }
        }

        *xt = xtt;
        *yt = ytt;
}

void dispersal(int x, int y, int xt, int yt, unsigned int i, int sex) {

        if (yt<0) {yt = YMAX; }; if (yt>=YMAX) {yt = 0; };  //building a tube

        if ((xt>=0)&&(xt<XMAX)) {

                double mort = par.disp_mort;

                if (sex==0) { //female dispersal
                        if (ran()>mort) {                     //a tube with absorbing x-conditions
                                world[xt][yt].femimmis.push_back(world[x][y].females.at(i)); //immigration
                        }
                }
                if (sex==1) { //male dispersal
                        if (ran()>mort) {
                                world[xt][yt].malimmis.push_back(world[x][y].males.at(i));
                        }
                }
        }
}

void dispersal_loop(int t, bool burn) {
        int xt,yt;

        for (int x = 0; x < XMAX; x++)
                for (int y = 0; y < YMAX; y++) {

                        int fpop = world[x][y].females.size();
                        int mpop = world[x][y].males.size();

                        if ((fpop+mpop)>0) { //dispersal only, when there are beetles in the patch

                                for (unsigned int f = 0; f < world[x][y].females.size(); f++) {
                                        double d = disp_prob(x,y,world[x][y].females.at(f));

                                        if (ran()<d) {
                                                find_target_patch(x,y,t,&xt,&yt,burn);
                                                dispersal(x,y,xt,yt,f,0); //dispersal
                                        } else {
                                                world[x][y].femimmis.push_back(world[x][y].females.at(f));
                                        }

                                }
                                for (unsigned int m = 0; m < world[x][y].males.size(); m++) {
                                        double d = disp_prob(x,y,world[x][y].males.at(m));
                                        if (ran()<d) {
                                                find_target_patch(x,y,t,&xt,&yt,burn);
                                                dispersal(x,y,xt,yt,m,1);
                                        } else {
                                                world[x][y].malimmis.push_back(world[x][y].males.at(m));
                                        }
                                }

                        }

                }

        for (int x = 0; x < XMAX; x++)
                for (int y = 0; y < YMAX; y++) {
                        world[x][y].females = world[x][y].femimmis;
                        world[x][y].males       = world[x][y].malimmis;
                        world[x][y].femimmis.clear();
                        world[x][y].malimmis.clear();
                }

}

//_____________________________________________________________________________________
//-Reproduction------------------------------------------------------------------------

void genetics(tInd *newborn, tInd mother, tInd father) {

        if (ran()>0.5) {newborn->la[0] = mother.la[0]; } else {newborn->la[0] = mother.la[1]; }
        if (ran()>0.5) {newborn->la[1] = father.la[0]; } else {newborn->la[1] = father.la[1]; }
        if (ran()>0.5) {newborn->ld[0] = mother.ld[0]; } else {newborn->ld[0] = mother.ld[1]; }
        if (ran()>0.5) {newborn->ld[1] = father.ld[0]; } else {newborn->ld[1] = father.ld[1]; }
        if (ran()>0.5) {newborn->le[0] = mother.le[0]; } else {newborn->le[0] = mother.le[1]; }
        if (ran()>0.5) {newborn->le[1] = father.le[0]; } else {newborn->le[1] = father.le[1]; }
        if (ran()>0.5) {newborn->ln[0] = mother.ln[0]; } else {newborn->ln[0] = mother.ln[1]; }
        if (ran()>0.5) {newborn->ln[1] = father.ln[0]; } else {newborn->ln[1] = father.ln[1]; }

        double expo = mean(newborn->le,2);
        double mut_la = pow(10,-expo);

        if (ran()<mut_la) {
          newborn->la[0] += gauss(VARIANCE);
          if (ran()<PROB_DEATH) {
              newborn->la[0] = 1000000;
          }
        }
        if (ran()<mut_la) {
          newborn->la[1] += gauss(VARIANCE);
          if (ran()<PROB_DEATH) {
              newborn->la[1] = 1000000;
          }
        }

        if (ran()<mut_la) {newborn->ld[0] += gauss(VARIANCE); }
        if (ran()<mut_la) {newborn->ld[1] += gauss(VARIANCE); }

        if (ran()<mut_la) {newborn->le[0] += gauss(VARIANCE); }
        if (ran()<mut_la) {newborn->le[1] += gauss(VARIANCE); }

        if (ran()<mut_la) {newborn->ln[0] += gauss(VARIANCE); }
        if (ran()<mut_la) {newborn->ln[1] += gauss(VARIANCE); }

}

double adaptation(int x, int y, tInd ind) {

        double diff = world[x][y].temp - mean(ind.la,2); //calculate difference between individual trait & temperature

        double surv = (1/(par.niche*sqrt(2*M_PI))*exp(-0.5*pow((diff/par.niche),2)))/
                      (1/(par.niche*sqrt(2*M_PI)));

        return(surv);
}

void reproduction(int x, int y, double lambda, double survival) {

        vector<tInd> girls,boys;
        tInd newborn;

        girls.clear(); boys.clear();

        for (unsigned int f = 0; f < world[x][y].females.size(); f++) {

                int kids = poisson(2.0*lambda*survival); //children numbers are poisson distributed
                int partner = floor(ran() * world[x][y].males.size()); //partner is searched at random

                for (int child = 0; child < kids; child++) {
                        float sex = ran(); //gender is randomly drawn

                        genetics(&newborn,
                                 world[x][y].females.at(f),
                                 world[x][y].males.at(partner));
                        double surv_prob = adaptation(x,y,newborn);

                        ana.sum_adap_x[x] += surv_prob;
                        ana.cnt_adap_x[x]++;

                        if (ran()<surv_prob) {
                                if (sex>=0.5) {girls.push_back(newborn); } else {boys.push_back(newborn); }
                        }

                }

        }

        world[x][y].females = girls;
        world[x][y].males = boys;
}

void reproduction_loop() {

        for (int x = 0; x < XMAX; x++)
                for (int y = 0; y < YMAX; y++) {

                        double lambda = lognorm(par.lambda_0,par.sigma);

                        unsigned long sum_pops; //the population sizes of all species are summed up
                        sum_pops = 0; //for calculating the logistic growth (i.e. all species
                        //compete for the same resources, determined by K)

                        sum_pops += world[x][y].females.size()+
                                    world[x][y].males.size();

                        double a = (par.lambda_0-1)/double(par.capacity);
                        double survival = 1/(1+a*double(sum_pops)); //assuming logistic growth

                        if ((world[x][y].females.size()>0)&&(world[x][y].males.size()>0)) {
                                reproduction(x,y,lambda,survival);
                        } else {
                                world[x][y].females.clear();
                                world[x][y].males.clear();
                        }

                }
}

//_____________________________________________________________________________________
//-Extinction--------------------------------------------------------------------------

void extinction() {
        for (int x = 0; x < XMAX; x++)
                for (int y = 0; y < YMAX; y++) {
                        if (ran()<par.ext_prob) {
                                world[x][y].females.clear();
                                world[x][y].males.clear();
                        }
                }
}

///////////////////////////////////////////////////////////////////////////////////////
//------------------------------ANALYZING FUNCTIONS----------------------------------//
///////////////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________________
//-Analysis----------------------------------------------------------------------------

void analyze(int t) {

}

void calculate_densities() {
        //double alldens = 0;
        //int cnt = 0;
        for (int x = 0; x < XMAX; x++)
                for (int y = 0; y < YMAX; y++) {
                        unsigned long sum_pops;
                        sum_pops = 0;
                        sum_pops += world[x][y].females.size()+world[x][y].males.size();
                        ana.dens[x][y] = double(sum_pops)/double(par.capacity);
                        //if (ana.dens[x][y]>0) {alldens+=ana.dens[x][y]; cnt++;}
                }
        //cout << alldens / double(cnt) << endl;
}

//_____________________________________________________________________________________
//-Write out Analysis Results

void write_dispersal(int t) {

        for (int x=0; x<XMAX; x++) {
                double ldsum = 0;
                unsigned long ldcnt = 0;
                for (int y=0; y<YMAX; y++) {
                        for (size_t f=0; f<world[x][y].females.size(); f++) {
                                ldsum += world[x][y].females.at(f).ld[0]+
                                         world[x][y].females.at(f).ld[1];
                                ldcnt += 2;
                        }
                        for (size_t m=0; m<world[x][y].males.size(); m++) {
                                ldsum += world[x][y].males.at(m).ld[0]+
                                         world[x][y].males.at(m).ld[1];
                                ldcnt += 2;
                        }
                }
                dispoutfile << ldsum/double(ldcnt) << "  ";
        }
        dispoutfile << endl;
}

void write_adaptation(int t) {
        for (int x=0; x<XMAX; x++) {
                double adap = ana.sum_adap_x[x] / double(ana.cnt_adap_x[x]);
                adapoutfile << adap << "  ";
                ana.sum_adap_x[x] = 0; ana.cnt_adap_x[x] = 0;
        }
        adapoutfile << endl;
}

void write_evolvability(int t) {
        for (int x=0; x<XMAX; x++) {
                double lesum = 0;
                unsigned long lecnt = 0;
                for (int y=0; y<YMAX; y++) {
                        for (size_t f=0; f<world[x][y].females.size(); f++) {
                                lesum += world[x][y].females.at(f).le[0]+
                                         world[x][y].females.at(f).le[1];
                                lecnt += 2;
                        }
                        for (size_t m=0; m<world[x][y].males.size(); m++) {
                                lesum += world[x][y].males.at(m).le[0]+
                                         world[x][y].males.at(m).le[1];
                                lecnt += 2;
                        }
                }
                evoloutfile << lesum/double(lecnt) << "  ";
        }
        evoloutfile << endl;
}

void write_densities(int t) {
        for (int x=0; x<XMAX; ++x) {

                double density = 0;
                int cnt = 0;

                for (int y=0; y<YMAX; ++y) {

                        unsigned long sum_pops = world[x][y].females.size()+world[x][y].males.size();

                        if (sum_pops>0) {
                                density += double(sum_pops)/double(par.capacity);
                                cnt++;
                        }

                }

                densoutfile << density/double(cnt) << "  ";

        }

        densoutfile << endl;
}

void write_occupancies(int t) {
        for (int x=0; x<XMAX; ++x) {
                int occ = 0;
                for (int y=0; y<YMAX; ++y) {
                        unsigned long sum_pops = world[x][y].females.size()+world[x][y].males.size();
                        if (sum_pops>0) {
                                occ++;
                        }
                }
                occoutfile << double(occ)/double(YMAX) << "  ";
        }
        occoutfile << endl;
}

void write_diversity(int t) {

        for (int x=0; x<XMAX; ++x) {

                vector<double> mismatch;
                vector<double> diversity;

                for (int y=0; y<YMAX; ++y) {
                        mismatch.clear();
                        for (size_t f=0; f<world[x][y].females.size(); ++f) {
                                mismatch.push_back(world[x][y].temp - world[x][y].females[f].la[0]);
                                mismatch.push_back(world[x][y].temp - world[x][y].females[f].la[1]);
                        }
                        for (size_t m=0; m<world[x][y].males.size(); ++m) {
                                mismatch.push_back(world[x][y].temp - world[x][y].males[m].la[0]);
                                mismatch.push_back(world[x][y].temp - world[x][y].males[m].la[1]);
                        }
                        diversity.push_back(var(mismatch));
                }

                divoutfile << mean(diversity) << "  ";
        }
        divoutfile << endl;
}

void write_neutral_diversity(int t) {

        for (int x=0; x<XMAX; ++x) {

                vector<double> neutr_allele;
                vector<double> diversity;

                for (int y=0; y<YMAX; ++y) {
                        neutr_allele.clear();
                        for (size_t f=0; f<world[x][y].females.size(); ++f) {
                                neutr_allele.push_back(world[x][y].females[f].ln[0]);
                                neutr_allele.push_back(world[x][y].females[f].ln[1]);
                        }
                        for (size_t m=0; m<world[x][y].males.size(); ++m) {
                                neutr_allele.push_back(world[x][y].males[m].ln[0]);
                                neutr_allele.push_back(world[x][y].males[m].ln[1]);
                        }
                        diversity.push_back(var(neutr_allele));
                }

                neutdivoutfile << mean(diversity) << "  ";
        }
        neutdivoutfile << endl;
}

void write_mut_diversity(int t) {

    for (int x=0; x<XMAX; ++x) {

        vector<double> alleles; alleles.clear();

        for (int y=0; y<YMAX; ++y) {

          for (size_t f=0; f<world[x][y].females.size(); ++f) {
            alleles.push_back(world[x][y].females[f].le[0]);
            alleles.push_back(world[x][y].females[f].le[1]);
          }
          for (size_t m=0; m<world[x][y].males.size(); ++m) {
            alleles.push_back(world[x][y].males[m].le[0]);
            alleles.push_back(world[x][y].males[m].le[1]);
          }

        }

        mutdivoutfile << var(alleles) << "  ";
    }
    mutdivoutfile << endl;
}

void write_disp_diversity(int t) {

    for (int x=0; x<XMAX; ++x) {

        vector<double> alleles; alleles.clear();

        for (int y=0; y<YMAX; ++y) {

          for (size_t f=0; f<world[x][y].females.size(); ++f) {
            alleles.push_back(world[x][y].females[f].ld[0]);
            alleles.push_back(world[x][y].females[f].ld[1]);
          }
          for (size_t m=0; m<world[x][y].males.size(); ++m) {
            alleles.push_back(world[x][y].males[m].ld[0]);
            alleles.push_back(world[x][y].males[m].ld[1]);
          }

        }

        dispdivoutfile << var(alleles) << "  ";
    }
    dispdivoutfile << endl;
}

///////////////////////////////////////////////////////////////////////////////////////
//----------------------------------MAIN FUNCTIONS-----------------------------------//
///////////////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________________
//-Initialize--------------------------------------------------------------------------

void initialize() {

        set_parameters();

        init_world();

        for (int x=0; x<XMAX; x++) {
                ana.sum_adap_x[x] = 0; ana.cnt_adap_x[x] = 0;
        }


}


void simulate() {

        for (int t = 0; t < par.tmax+BURN_IN_TIME; t++) {

                bool burn;
                if (t<BURN_IN_TIME) {burn = 1; } else {burn = 0; }

                cout << "t: " << t << endl;
                calculate_densities();

                dispersal_loop(t,burn);
                if (((t+1) % TIME_OUT_INT)==0) {
                        write_adaptation(t);
                }
                reproduction_loop();
                extinction();
                analyze(t);

                if (((t+1) % TIME_OUT_INT)==0) {
                        write_dispersal(t);
                        write_evolvability(t);
                        write_densities(t);
                        write_occupancies(t);
                        write_diversity(t);
                        write_neutral_diversity(t);
                        write_mut_diversity(t);
                        write_disp_diversity(t);
                }

        }
}

int main() {

        cout.setf(ios_base::scientific);  //exponential data output
        //cout.precision(5);			//precision of data output

        specify_rng(time(NULL));          //random randomization?

        initialize();

        simulate();

        densoutfile.close();
        adapoutfile.close();
        dispoutfile.close();
        evoloutfile.close();
        occoutfile.close();
        divoutfile.close();
        mutdivoutfile.close();
        dispdivoutfile.close();

        return 0;
}
