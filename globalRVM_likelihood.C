#include <config.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex>

#ifdef HAVE_MINUIT
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/FCNBase.h"
#include "globalRVMKJ_Fcn.h"
#else
#include "RVMnest.h"
#endif

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

void globalRVMLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context) {

	MNStruct *par = ((MNStruct *)context);

	char label[16];	
	double Qu,Uu;
        double chi = 0.0, Ln, PA2;
	double Ltot = 0., logdetN = 0.;
	double cosPA2, sinPA2;
	complex <double> L, arg;
	double phi, eta, xsi;
	double M, cslam, snlam, cseta, sneta, csi, sni, csphi, snphi, csdel, sndel;
	double lambda, beta;
	double dt;

	int npar = 4;
	par->alpha = Cube[0] * (par->r_alpha[1]*DEG_TO_RAD - par->r_alpha[0]*DEG_TO_RAD) + par->r_alpha[0]*DEG_TO_RAD;
	par->delta = Cube[1] * (par->r_delta[1]*DEG_TO_RAD - par->r_delta[0]*DEG_TO_RAD) + par->r_delta[0]*DEG_TO_RAD;
	par->phase0 = Cube[2] * (par->r_Phi0[1]*DEG_TO_RAD - par->r_Phi0[0]*DEG_TO_RAD) + par->r_Phi0[0]*DEG_TO_RAD;
	par->psi00 = Cube[3] * M_PI - M_PI_2;

	// RVM center Phase of the profiles - Can be marginalized
	if (!par->margin_phi0) {
	    for (unsigned int j = 0; j < par->n_epoch; j++) {
		par->phi0[j] = Cube[npar] * (par->r_phi0[1]*DEG_TO_RAD - par->r_phi0[0]*DEG_TO_RAD) + par->r_phi0[0]*DEG_TO_RAD;
		npar++;
	    }
	} else {
	    for (unsigned int j = 0; j < par->n_epoch; j++) par->phi0[j] = 5.0 * DEG_TO_RAD;
	}

	// Precession rate
	if (!par->prate_fixed) {
	    par->prate = Cube[npar] * (par->r_prate[1]-par->r_prate[0]) + par->r_prate[0];
	    npar++;
	}

	// Inclination angle
	if (!par->inc_fixed) {
 	    par->inc = Cube[npar] * (par->r_inc[1]*DEG_TO_RAD - par->r_inc[0]*DEG_TO_RAD) + par->r_inc[0]*DEG_TO_RAD;
	    npar++;
	}
	

	// EFACs
	if (par->have_efac) {
	    for (unsigned int j = 0; j < par->n_epoch; j++) {
	        par->efac[j] = Cube[npar] * (par->r_efac[1]-par->r_efac[0]) + par->r_efac[0];
		npar++;
	    }
	} else {
	    for (unsigned int j = 0; j < par->n_epoch; j++) par->efac[j] = 1.;
	}

	// Psi0 jumps between different datasets
	//par->psi_jumps = (double *) malloc(par->njump * sizeof(double));
	if (par->njump) {
	  if(!par->psi_jump_fixed)
	    {
	      for (unsigned  l=0; l<par->njump; l++) {
		par->psi_jumps[l] = Cube[npar] * M_PI/2 + M_PI/4.;
		npar++;
	      }
	    }
	  //else
	  //{
	  //  for (unsigned  l=0; l<par->njump; l++) par->psi_jumps[l] = M_PI/2.;
	  //}
	}

	double rmsQ;
	double rmsU;

	vector<double> vrmsQ, vrmsU;
	vector< vector<double> > vQ, vU, vphase;
	vQ.resize(par->n_epoch); vU.resize(par->n_epoch); vphase.resize(par->n_epoch);

	for(unsigned int j = 0; j < par->n_epoch; j++) {
	  vrmsQ.push_back(par->rmsQ[j]);
	  vrmsU.push_back(par->rmsU[j]);
	  for(unsigned int i = 0; i < par->npts[j]; i++) {
	    vQ[j].push_back(par->Q[j][i]);
	    vU[j].push_back(par->U[j][i]);
	    vphase[j].push_back(par->phase[j][i]);
	  }
	}

#ifdef HAVE_MINUIT
	// Marginalize over Phi0
	if (par->margin_phi0) {
	    //RVM_Fcn fFCN();
	    RVM_Fcn fFCN(vphase, vQ, vU, vrmsQ, vrmsU);
	    fFCN.set_params(par);

	    MnUserParameters upar;
	    for(unsigned int j = 0; j < par->n_epoch; j++) {
		sprintf(label, "phi0_%d", j);
		upar.Add(label, par->phi0[j], 5.);
	    }

	    MnMinimize migrad(fFCN, upar);
	    //cout << upar << endl;
	    FunctionMinimum min = migrad();
	    //cout<<"minimum: "<<min<<endl;

	    for(unsigned int j = 0; j < par->n_epoch; j++) {
		sprintf(label, "phi0_%d", j);
		par->phi0[j] = min.UserState().Value(label);
		//printf("phi0[%02d] = %lf\n",j, par->phi0[j]);
	    }
	}
#endif

	get_globalRVM_chi2(par);

	// shift the reference P.A. by 90 degrees and ensure that PA0 lies on -pi/2 -> pi/2
	if (par->Ltot < 0.0) {
	    par->psi00 += M_PI /2.;
	}
	par->psi00 = atan( tan(par->psi00) );
	
	npar = 0;
	Cube[0] = par->alpha * RAD_TO_DEG;
	Cube[1] = par->delta * RAD_TO_DEG;
	Cube[2] = par->phase0 * RAD_TO_DEG;
	Cube[3] = par->psi00 * RAD_TO_DEG;
	npar = 4;

	// phi0
	if (!par->margin_phi0) {
	    for (unsigned int j = 0; j < par->n_epoch; j++) {
		Cube[npar] = par->phi0[j] * RAD_TO_DEG;
		npar++;
	    }
	}

	if (!par->prate_fixed) {
	    Cube[npar] = par->prate;
	    npar++;
	}
	if (!par->inc_fixed) {
	    Cube[npar] = par->inc * RAD_TO_DEG;
	    npar++;
	}

	if (par->have_efac) {
	    for (unsigned int j = 0; j < par->n_epoch; j++) {
	        Cube[npar] = par->efac[j];
		npar++;
	    }
	}

        if (par->njump) {
          for (unsigned  l=0; l<par->njump; l++) {
	    Cube[npar] = par->psi_jumps[l] * RAD_TO_DEG;
            npar++;
          }
        }

        lnew = -1.*par->chi/2 - 0.5*par->logdetN;
	//free(par->psi_jumps);
}
