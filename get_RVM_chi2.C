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
	double Qu,Uu;

	// Off centred dipole stuff
	double phas = par->phas;
	double Minc = par->Minc;
	double ita = par->ita;
	double eps = par->eps;

	//eps = 0.0;
	par->logdetN = 0.0;
	par-> chi = 0.0;
	double phi_off;
	
	Qu = rmsQ*rmsQ*par->efac[0]*par->efac[0];
	Uu = rmsU*rmsU*par->efac[0]*par->efac[0];

	for(unsigned int i = 0; i < par->npts[0]; i++) {

	  // Offset in phase between the two RVM inflexion points. The offset is set at phase 180 deg.
	  if (par->have_aberr_offset && par->phase[0][i] > 180.*DEG_TO_RAD && par->phase[0][i] < 360.*DEG_TO_RAD)
	    phi_off = par->phi_aberr_offset[0];
	  else phi_off = 0.0;

	  
	  if (par->have_offset_dipole)
	    PA = get_offRVM(alpha, beta, phi0 + phi_off, psi0, par->phase[0][i], phas, Minc, ita, eps);
	  else
	    PA = get_RVM(alpha, beta, phi0 + phi_off, psi0, par->phase[0][i]);

	  //cout << alpha / DEG_TO_RAD << " " << beta  / DEG_TO_RAD<< " "<< phi0  / DEG_TO_RAD << " " << rmsQ << " " << PA << endl;
	  Ln = (par->Q[0][i]*cos(2*PA) / Qu + par->U[0][i]*sin(2*PA) / Uu) 
	    / (cos(2*PA)*cos(2*PA) / Qu + sin(2*PA)*sin(2*PA) / Uu);
	  
	  arg = complex<double>(0., 2 * PA);
	  L =  Ln * exp (arg);
	  par->chi += (par->Q[0][i]-real(L))*(par->Q[0][i]-real(L))/Qu
	    + (par->U[0][i]-imag(L))*(par->U[0][i]-imag(L))/Uu;
	  par->logdetN += log(Uu) + log(Qu);
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
	      if (par->have_aberr_offset && (i+.5)/par->nbin[0]*M_PI*2 > 180.*DEG_TO_RAD && (i+.5)/par->nbin[0]*M_PI*2 < 360.*DEG_TO_RAD)
		phi_off = par->phi_aberr_offset[0];
	      else phi_off = 0.0;
	      PA = get_RVM(alpha, beta, phi0+phi_off, psi0, (i+.5)/par->nbin[0]*M_PI*2);
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
	    
}
