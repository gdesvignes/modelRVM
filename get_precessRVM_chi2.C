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

  double Ln, PA,PA2;
    complex <double> L, arg;
    int totnbpts=0;
    double prev_chi = 0.0;
    double rms2Q, rms2U, alpha_tmp;
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
	for(unsigned int i = 0; i < par->npts[j]; i++) {
	    if (par->pmodel==0) // pmodel=0 forced precession, pmodel=1 free precession 
		PA2 = 2*get_RVM(par->alpha, par->beta[j], par->phi0[j], par->psi0[j], par->phase[j][i]);
	    else if (par->pmodel==1) // in this case, par->alpha is really zeta
		PA2 = 2*get_RVM(par->alpha-par->beta[j], par->beta[j], par->phi0[j], par->psi0[j], par->phase[j][i]);
	    else if (par->pmodel==2) {
	      alpha_tmp = par->alpha + par->alpha1*(par->epoch[j]-par->epoch[0]) + par->alpha2*pow(par->epoch[j]-par->epoch[0], 2);
	      PA2 = 2*get_RVM(alpha_tmp, par->beta[j], par->phi0[j], par->psi0[j], par->phase[j][i]);
	    }
	    else if (par->pmodel==3) {
	      alpha_tmp = par->alpha + par->alpha1*(par->epoch[j]-par->epoch[0]) + par->alpha2*pow(par->epoch[j]-par->epoch[0], 2);
	      PA2 = 2*get_RVM(alpha_tmp - par->beta[j], par->beta[j], par->phi0[j], par->psi0[j], par->phase[j][i]);
	    }

	    Ln = (par->Q[j][i]*cos(PA2) / rms2Q + par->U[j][i]*sin(PA2) / rms2U) 
		/ (cos(PA2)*cos(PA2) / rms2Q + sin(PA2)*sin(PA2) / rms2U);
	    
	    arg = complex<double>(0., PA2);
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
	    if (par->pmodel==0) // pmodel=0 forced precession, pmodel=1 free precession 
	      PA = get_RVM(par->alpha, par->beta[j], par->phi0[j], par->psi0[j], (i+.5)*360./par->nbin[j] * M_PI/180.);
	    else // in this case, par->alpha is really zeta
	      PA = get_RVM(par->alpha-par->beta[j], par->beta[j], par->phi0[j], par->psi0[j], (i+.5)*360/par->nbin[j] * M_PI/180.);
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


	  
	  //cout << "Chi2 : " << par->chi-prev_chi<< endl;
	  //cout << "Nbpts: " << par->npts[j] << endl;
	  prev_chi = par->chi;
	}
    }
}
