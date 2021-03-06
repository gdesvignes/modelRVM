#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <complex>
#include <float.h>

#include "RVMnest.h"

#define DEG_TO_RAD	(M_PI/180.0)
#define RAD_TO_DEG	(180.0/M_PI)
#define VERBOSE 0

using namespace std;

/******************************************** loglikelihood routine ****************************************************/

// Now an example, sample an egg box likelihood

// Input arguments
// ndim                                                 = dimensionality (total number of free parameters) of the problem
// npars                                                = total number of free plus derived parameters
// context                                              void pointer, any additional information
//
// Input/Output arguments
// Cube[npars]                                          = on entry has the ndim parameters in unit-hypercube
//                                                      on exit, the physical parameters plus copy any derived parameters you want to store with the free parameters
//       
// Output arguments
// lnew                                                 = loglikelihood

void get_globalRVM_chi2(MNStruct *par) {

	char label[16];	
	int ie=0;
	int ix = 0;
	int totnbpts=0;
	double Qu,Uu;
        double Ln, PA2;
	double cosPA2, sinPA2;
	complex <double> L, arg;
	double phi, eta, xsi;
	double M, cslam, snlam, cseta, sneta, csi, sni, csphi, snphi, csdel, sndel;
	double lambda, beta;
	double dt;
	double rmsQ;
        double rmsU;
	double prev_chi = 0.0;
	stringstream s;
	ofstream myf;
	par->chi = 0.0;
	par->logdetN = 0.0;

	par->omega = par->prate / 365.25 * DEG_TO_RAD;

	for (unsigned int j = 0; j < par->n_epoch; j++) {
	    par->Ltot[j] = 0.0;

	    totnbpts += par->npts[j];
	    dt = par->epoch[j] - par->epoch[0];

	    // -- compute beta and psi first --
	    phi = par->phase0 + par->omega * dt; // in rad
	    csdel = cos(par->delta);
	    sndel = sin(par->delta);

	    csi = cos(par->inc);
	    sni = sin(par->inc);
	    csphi = cos(phi);
	    snphi = sin(phi);

	    cslam = csdel * csi - sndel*sni*csphi;
	    snlam = sqrt(1.0-cslam*cslam);

	    //  as 0 < lambda < 180, sin(lambda) >0
	    //  hence we can just compute acos(cos(lambda))

	    lambda = acos(cslam);
	    beta = M_PI - par->alpha - lambda;

            cseta = sndel*snphi/snlam;
            sneta = (cslam*csi-csdel)/sni/snlam;

            eta = atan2(sneta,cseta);
            //psi0 = eta;

	    // Take into account Psi Jumps
	    for(unsigned l=0; l<par->njump; l++) {
	      if (par->psi_jump_MJD[2*l] < par->epoch[j] && par->epoch[j] < par->psi_jump_MJD[2*l+1]) eta += par->psi_jumps[l] ;
	    }

	    if (VERBOSE)
	      {
		printf("delta=%lf  inc=%lf  phi=%lf\n", par->delta * 180./M_PI, par->inc * 180./M_PI, phi);
		printf("%2d  epoch=%lf  lambda=%lf  beta=%lf  eta=%lf\n", j, par->epoch[j], lambda * 180./M_PI, beta * 180./M_PI, eta * 180./M_PI);
	      }

	    rmsQ = par->rmsQ[j];
	    rmsU = par->rmsU[j];

	    xsi = par->alpha + beta;
	    
	    ie=0;
	    Qu = rmsQ*rmsQ;
	    Uu = rmsU*rmsU;
	    if (par->have_efac) {
	      Qu *= par->efac[ie]*par->efac[ie];
	      Uu *= par->efac[ie]*par->efac[ie];
	    }

	    //if (par->epoch[j] < 57700) ix ++;

	    for(unsigned int i = 0; i < par->npts[j]; i++) {
		double phi_off = 0.0;
		if (par->epoch[j] < 57300 && par->have_aberr_offset && par->phase[j][i] > 180.*DEG_TO_RAD && par->phase[j][i] < 360.*DEG_TO_RAD)
		    phi_off = par->phi_aberr_offset[ix];
		PA2 = 2*get_RVM(par->alpha, beta, par->phi0[j] + phi_off, par->psi00 + eta, par->phase[j][i]);

		cosPA2 = cos(PA2);
		sinPA2 = sin(PA2);

		Ln = (par->Q[j][i]*cosPA2 / (Qu) + par->U[j][i]*sinPA2 / (Uu)) 
		  / (cosPA2*cosPA2 / (Qu) + sinPA2*sinPA2 / (Uu));

		//arg = 1i * PA2;
		arg = std::complex<double>(0., PA2);
		L =  Ln * exp (arg);
		par->chi += (par->Q[j][i]-real(L))*(par->Q[j][i]-real(L))/(Qu)
			+ (par->U[j][i]-imag(L))*(par->U[j][i]-imag(L))/(Uu);

		par->logdetN += log(Uu) + log(Qu);
		par->Ltot[j] += Ln;
	      }

	    if (par->do_plot)
	      {
		cout << "File #"<<j << " red. chi2= " << (par->chi-prev_chi)/par->npts[j] << endl;
		prev_chi = par->chi;

		s.str("");
                s << (int)par->epoch[j]<< "-prof.log";
                string result = s.str();
                myf.open(result.c_str());
		double PA, Lv, Lve;
		for(unsigned int i = 0; i < par->nbin[j]; i++)
		  { 
		    double phi_off = 0.0;
		    if (par->epoch[j] < 57300 && par->have_aberr_offset && (i+.5)/par->nbin[j]*M_PI*2 > 180.*DEG_TO_RAD && (i+.5)/par->nbin[j]*M_PI*2 < 360.*DEG_TO_RAD)
		      phi_off = par->phi_aberr_offset[ix];
		    PA2 = get_RVM(par->alpha, beta, par->phi0[j] + phi_off, par->psi00 + eta, (i+.5)/par->nbin[j]*M_PI*2);
		    myf << i*360./par->nbin[j] << " "<< par->I[j][i] << " "<< par->L[j][i] << " "<< par->V[j][i]<< " " <<  PA2 * 180./M_PI << endl;
		  }
		myf.close();

		s.str("");
                s << (int)par->epoch[j]<< "-PA.log";
                result = s.str();
                myf.open(result.c_str());
                for(unsigned int i = 0; i < par->npts[j]; i++)
                  {
		    
		    PA = 0.5 * atan2(par->U[j][i], par->Q[j][i]) * 180./ M_PI;
		    Lv = sqrt(pow(par->U[j][i],2) + pow(par->Q[j][i],2));
		    Lve = 28.65 * par->rmsI[j]/Lv;

                    myf << par->phase[j][i]*180./M_PI << " "<< PA << " " << Lve << endl;
                  }
                myf.close();
	      }
	    
	    // Increment for aberration 
	    if (par->epoch[j] < 57300) ix++;
	}

	//if (par->Ltot < 0.0) {
	//  par->psi00 += M_PI /2.;
	//}

	// Print summary
	if (par->do_plot)
	  {
	    cout << "Final Chi2 : " << par->chi<< endl;
	    cout << "Tot nbpts: " << totnbpts << endl;
	  }
}
