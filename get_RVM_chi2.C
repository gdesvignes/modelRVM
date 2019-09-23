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

void get_RVM_chi2(MNStruct *par) {

        double Ln, PA;
	complex <double> L, arg;

	double alpha = par->alpha;
	double beta = par->beta;
	double phi0 = par->phi0[0];
	double psi0 = par->psi0;
	double rmsQ = par->rmsQ[0];
	double rmsU = par->rmsU[0];
	double Ltot = 0.0;

	//alpha = 100 * DEG_TO_RAD;
	//beta = -16 * DEG_TO_RAD;
	//phi0 = 90 * DEG_TO_RAD;
	//psi0 = 10 * DEG_TO_RAD;

	par-> chi = 0.0;
	//cout << alpha / DEG_TO_RAD << " " << beta  / DEG_TO_RAD<< " "<< phi0  / DEG_TO_RAD << " " << rmsQ << endl;

	for(unsigned int i = 0; i < par->npts[0]; i++) {
		PA = get_RVM(alpha, beta, phi0, psi0, par->phase[0][i]);
		//cout << alpha / DEG_TO_RAD << " " << beta  / DEG_TO_RAD<< " "<< phi0  / DEG_TO_RAD << " " << rmsQ << " " << PA << endl;
		Ln = (par->Q[0][i]*cos(2*PA) / (rmsQ*rmsQ) + par->U[0][i]*sin(2*PA) / (rmsU*rmsU)) 
			/ (cos(2*PA)*cos(2*PA) / (rmsQ*rmsQ) + sin(2*PA)*sin(2*PA) / (rmsU*rmsU));

		arg = complex<double>(0., 2 * PA);
		L =  Ln * exp (arg);
		par->chi += (par->Q[0][i]-real(L))*(par->Q[0][i]-real(L))/(rmsQ*rmsQ)
			+ (par->U[0][i]-imag(L))*(par->U[0][i]-imag(L))/(rmsU*rmsU);
		Ltot += Ln;
	}	
	//cout << par->chi << endl;


        // shift the reference P.A. by 90 degrees and ensure that PA0 lies on -pi/2 -> pi/2
        if (Ltot < 0.0) {
            psi0 += M_PI /2.;
        }
        par->psi0 = atan( tan(psi0) );


	if (par->do_plot) {
	  
	  stringstream s;
	  s.str("");
	  s << (int)par->epoch[0]<< "-prof.log";
	  string result = s.str();
	  ofstream myf;
	  myf.open(result.c_str());

	  double Lv, Lve;
	  for(unsigned int i = 0; i < par->nbin[0]; i++)
	    {
	      PA = get_RVM(alpha, beta, phi0, psi0, (i+.5)/par->nbin[0]*M_PI*2);
	      myf << i*360./par->nbin[0] << " "<< par->I[0][i] << " "<< par->L[0][i] << " "<< par->V[0][i]<< " " <<  PA * 180./M_PI << endl;
	    }
	  myf.close();

	  s.str("");
	  

	  s << (int)par->epoch[0]<< "-PA.log";
	  result = s.str();
	  myf.open(result.c_str());
	  for(unsigned int i = 0; i < par->npts[0]; i++) {
	    PA = 0.5 * atan2(par->U[0][i], par->Q[0][i]) * 180./ M_PI;
	    Lv = sqrt(pow(par->U[0][i],2) + pow(par->Q[0][i],2));
	    Lve = 28.65 * par->rmsI[0]/Lv;
	    myf << par->phase[0][i]*180./M_PI << " "<< PA << " " << Lve << endl;
	  }
	  myf.close();


	  cout << "Final Chi2 : " << par->chi<< endl;
	  cout << "Tot nbpts: " << par->npts[0] << endl;
	}
	//exit(0);
	    
}