#include <config.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex>



#include "RVMnest.h"
//#endif

#define DEG_TO_RAD	(M_PI/180.0)
#define RAD_TO_DEG	(180.0/M_PI)

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

void precessLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context) {

	MNStruct *par = ((MNStruct *)context);

	char label[16];	
	double Qu,Uu;
        double chi = 0.0, Ln, PA2;
	double Ltot = 0.;
	double cosPA2, sinPA2;
	complex <double> L, arg;
	double phi, eta, xsi;
	double M, cslam, snlam, cseta, sneta, csi, sni, csphi, snphi, csdel, sndel;
	double lambda, beta;
	double dt;

	int npar=0;
	par->alpha = Cube[0] * (par->r_alpha[1]*DEG_TO_RAD - par->r_alpha[0]*DEG_TO_RAD) + par->r_alpha[0]*DEG_TO_RAD;
	npar++;
	for (unsigned int j = 0; j < par->n_epoch; j++) {
	    par->beta[j] = Cube[npar] * (par->r_beta[1]*DEG_TO_RAD - par->r_beta[0]*DEG_TO_RAD) + par->r_beta[0]*DEG_TO_RAD;
	    npar++;
	    par->phi0[j] = Cube[npar] * (par->r_phi0[1]*DEG_TO_RAD - par->r_phi0[0]*DEG_TO_RAD) + par->r_phi0[0]*DEG_TO_RAD;
            npar++;
	    par->psi0[j] = Cube[npar] * (par->r_psi0[1]*DEG_TO_RAD - par->r_psi0[0]*DEG_TO_RAD) + par->r_psi0[0]*DEG_TO_RAD;
            npar++;
	    
	}

#if 0
	// EFACs
	if (par->have_efac) {
	  for (unsigned int j = 0; j < 2; j++) {
	    par->efac[j] = Cube[npar] * (par->r_efac[1]-par->r_efac[0]) + par->r_efac[0];
	    npar++;
	  }
	} else {
	  for (unsigned int j = 0; j < 2; j++) par->efac[j] = 1.;
	}
#endif
	// Compute chi**2
	get_precessRVM_chi2(par);

	npar = 0;
	Cube[npar] = par->alpha * RAD_TO_DEG;
	npar = 1;
        for (unsigned int j = 0; j < par->n_epoch; j++) {
	    Cube[npar] = par->beta[j] * RAD_TO_DEG; npar++;
	    Cube[npar] = par->phi0[j] * RAD_TO_DEG; npar++;
	    Cube[npar] = par->psi0[j] * RAD_TO_DEG; npar++;
        }

	if (par->have_efac) {
	    for (unsigned int j = 0; j < 2; j++) {
	        Cube[npar] = par->efac[j];
		npar++;
	    }
	}

        //lnew = -1.*par->chi/2 - 0.5*par->logdetN;
	lnew = -1.*par->chi/2.;
}
