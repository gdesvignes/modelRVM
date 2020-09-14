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

void get_precessmRVM_chi2(MNStruct *par) {

    double Ln, PA;
    complex <double> L, arg;
    int totnbpts=0;
    double alpha = par->alpha;
    double prev_chi = 0.0;
    double rms2Q, rms2U;
    double beta;

    int debug = 0;
    if (debug) {
	//alpha = 100 * DEG_TO_RAD;
	//par->beta[0] = 3 * DEG_TO_RAD;
	//par->phi0[0] = 180 * DEG_TO_RAD;
	//par->psi0[0] = 45 * DEG_TO_RAD;
	//par->beta[1] = -5 * DEG_TO_RAD;
        //par->phi0[1] = 180 * DEG_TO_RAD;
        //par->psi0[1] = 45 * DEG_TO_RAD;
	par->do_plot=1;
    }

    par->chi = 0.0;
    par->logdetN = 0.0; 
    //cout << alpha / DEG_TO_RAD << " " << par->beta[0]  / DEG_TO_RAD<< " "<< par->phi0[0]  / DEG_TO_RAD << " " << par->psi0[0] << endl;
    for (unsigned int j = 0; j < par->n_epoch; j++) {

	rms2Q = par->rmsQ[j] * par->rmsQ[j];
	rms2U = par->rmsU[j] * par->rmsU[j];

	if (par->have_efac) {
	    rms2Q *= par->efac[0]*par->efac[0];
	    rms2U *= par->efac[0]*par->efac[0];
	}
	
	totnbpts += par->npts[j];
	par->Ltot[j] = 0;

	// Get beta from the model
	beta = par->beta0 + par->betaamp * sin(2*M_PI* par->pperiod / exp(par->epoch[j] - par->mjd0)  * (par->epoch[j] - par->mjd0) + par->pphase);

	for(unsigned int i = 0; i < par->npts[j]; i++) {
	    PA = get_RVM(alpha, beta, par->phi0[j], par->psi0[j], par->phase[j][i]);

	    Ln = (par->Q[j][i]*cos(2*PA) / rms2Q + par->U[j][i]*sin(2*PA) / rms2U) 
		/ (cos(2*PA)*cos(2*PA) / rms2Q + sin(2*PA)*sin(2*PA) / rms2U);
	    
	    arg = complex<double>(0., 2 * PA);
	    L =  Ln * exp (arg);
	    par->chi += (par->Q[j][i]-real(L))*(par->Q[j][i]-real(L)) / rms2Q
		+ (par->U[j][i]-imag(L))*(par->U[j][i]-imag(L)) / rms2U;
	    par->Ltot[j] += Ln;
	    par->logdetN += log(rms2Q) + log(rms2U);
	}	
	//cout << par->chi << endl;
	
	
        // shift the reference P.A. by 90 degrees and ensure that PA0 lies on -pi/2 -> pi/2
        //if (Ltot < 0.0) {
        //    par->psi0[j] += M_PI /2.;
        //}
        par->psi0[j] = atan( tan(par->psi0[j]) );

	if (par->do_plot) {
	  
	  stringstream s;
	  s.str("");
	  //s << par->epoch[j]<< "-prof.log";
	  s << "Profile_" << j << "-prof.log";
	  string result = s.str();
	  ofstream myf;
	  myf.open(result.c_str());
	  myf << "# " << par->epoch[j] << " " << endl;

	  double Lv, Lve;
	  for (unsigned int i = 0; i < par->nbin[j]; i++) {
	      PA = get_RVM(alpha, beta, par->phi0[j], par->psi0[j], (i+.5)*360./par->nbin[j] * DEG_TO_RAD);
	      myf << (i+.5)*360./par->nbin[j] << " "<< par->I[j][i] << " "<< par->L[j][i] << " "<< par->V[j][i]<< " " <<  PA * 180./M_PI << endl;
	  }
	  myf.close();

	  s.str("");
	  

	  //s << par->epoch[j]<< "-PA.log";
	  s << "Profile_" << j << "-PA.log";
	  result = s.str();
	  myf.open(result.c_str());
	  for(unsigned int i = 0; i < par->npts[j]; i++) {
	    PA = 0.5 * atan2(par->U[j][i], par->Q[j][i]) * 180./ M_PI;
	    Lv = sqrt(pow(par->U[j][i],2) + pow(par->Q[j][i],2));
	    Lve = 28.65 * par->rmsI[j]/Lv;
	    myf << par->phase[j][i]*180./M_PI << " "<< PA << " " << Lve << endl;
	  }
	  myf.close();


	  cout << "Beta : " << beta / DEG_TO_RAD  << endl;
	  cout << "Chi2 : " << par->chi-prev_chi<< endl;
	  cout << "Nbpts: " << par->npts[j] << endl;
	  prev_chi = par->chi;
	}
    }
    if (debug) exit(0);
}
