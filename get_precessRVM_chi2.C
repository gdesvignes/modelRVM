#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include "RVMnest.h"
#include <complex>

#define DEG_TO_RAD (M_PI/180.0)

using namespace std;

void get_precessRVM_chi2(MNStruct *par) {

    double Ln, PA;
    complex <double> L, arg;
    
    double alpha = par->alpha;
    double rmsQ = par->rmsQ[0];
    double rmsU = par->rmsU[0];
    double Ltot = 0.0;

    par-> chi = 0.0;
    //cout << alpha / DEG_TO_RAD << " " << beta  / DEG_TO_RAD<< " "<< phi0  / DEG_TO_RAD << " " << rmsQ << endl;
    for (unsigned int j = 0; j < par->n_epoch; j++) {
	totnbpts += par->npts[j];
	for(unsigned int i = 0; i < par->npts[j]; i++) {
	    PA = get_RVM(alpha, par->beta[j], par->phi0[j], par->psi0[j], par->phase[j][i]);

	    Ln = (par->Q[j][i]*cos(2*PA) / (par->rmsQ[j]*par->rmsQ[j]) + par->U[j][i]*sin(2*PA) / (par->rmsU[j]*par->rmsU[j])) 
		/ (cos(2*PA)*cos(2*PA) / (par->rmsQ[j]*par->rmsQ[j]) + sin(2*PA)*sin(2*PA) / (par->rmsU[j]*par->rmsU[j]));
	    
	    arg = complex<double>(0., 2 * PA);
	    L =  Ln * exp (arg);
	    par->chi += (par->Q[j][i]-real(L))*(par->Q[j][i]-real(L))/(par->rmsQ[j]*par->rmsQ[j])
		+ (par->U[j][i]-imag(L))*(par->U[j][i]-imag(L))/(par->rmsU[j]*par->rmsU[j]);
	    Ltot += Ln;
	}	
	//cout << par->chi << endl;
	
	
        // shift the reference P.A. by 90 degrees and ensure that PA0 lies on -pi/2 -> pi/2
        if (Ltot < 0.0) {
            par->psi0[j] += M_PI /2.;
        }
        par->psi0[j] = atan( tan(par->psi0[j]) );


	if (par->do_plot) {
	  
	  stringstream s;
	  s.str("");
	  s << (int)par->epoch[j]<< "-prof.log";
	  string result = s.str();
	  ofstream myf;
	  myf.open(result.c_str());

	  double Lv, Lve;
	  for(unsigned int i = 0; i < par->nbin[j]; i++)
	    {
	      PA = get_RVM(alpha, par->beta[j], par->phi0[j], par->psi0[j], par->phase[j][i]);
	      myf << par->phase[j][i] << " "<< par->I[j][i] << " "<< par->L[j][i] << " "<< par->V[j][i]<< " " <<  PA * 180./M_PI << endl;
	    }
	  myf.close();

	  s.str("");
	  

	  s << (int)par->epoch[j]<< "-PA.log";
	  result = s.str();
	  myf.open(result.c_str());
	  for(unsigned int i = 0; i < par->npts[j]; i++) {
	    PA = 0.5 * atan2(par->U[j][i], par->Q[j][i]) * 180./ M_PI;
	    Lv = sqrt(pow(par->U[j][i],2) + pow(par->Q[j][i],2));
	    Lve = 28.65 * par->rmsI[j]/Lv;
	    myf << par->phase[j][i]*180./M_PI << " "<< PA << " " << Lve << endl;
	  }
	  myf.close();


	  cout << "Final Chi2 : " << par->chi<< endl;
	  cout << "Tot nbpts: " << par->npts[j] << endl;
	}
	//exit(0);
    }
}
