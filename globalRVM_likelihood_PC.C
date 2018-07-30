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
#endif

#include "RVMnest.h"


#include "globalRVM_likelihood_PC.h"

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


void prior (double cube[], double theta[], int nDims,  void *context)
{
  int ipar = 0;

  theta[ipar] = cube[ipar] * (sp->r_alpha[1] - sp->r_alpha[0]) + sp->r_alpha[0];
  ipar++;
  theta[ipar] = cube[ipar] * (sp->r_delta[1] - sp->r_delta[0]) + sp->r_delta[0];
  ipar++;
  theta[ipar] = cube[ipar] * (sp->r_Phi0[1] - sp->r_Phi0[0]) + sp->r_Phi0[0];
  ipar++;
  
  if (sp->sin_psi) {
    theta[ipar] = cube[ipar] * 2 - 1.;
    ipar++;
    theta[ipar] = cube[ipar] * 2 - 1;
    ipar++;
  } else {
    theta[ipar] = cube[ipar] * 180;
    ipar++;
  }

  
  // RVM center Phase of the profiles - Can be marginalized                                                                 
  if (!sp->margin_phi0) {
    for (unsigned int j = 0; j < sp->n_epoch; j++) {
      theta[ipar] = cube[ipar] * (sp->r_phi0[1] - sp->r_phi0[0]) + sp->r_phi0[0];
      ipar++;
    }
  }

  // include an offset between the PA curves of the MP and IP                                                               
  if (sp->have_aberr_offset) {
    for (unsigned int j = 0; j < sp->nfiles_aberr; j++) {
      theta[ipar] = (cube[ipar] *  12 - 6);
      ipar++;
    }
  }

  // Precession rate                                                                                                        
  if (!sp->prate_fixed) {
    theta[ipar] = cube[ipar] * (sp->r_prate[1] - sp->r_prate[0]) + sp->r_prate[0];
    ipar++;
  }

  // Inclination angle                                                                                                      
  if (!sp->inc_fixed) {
    theta[ipar] = cube[ipar] * (sp->r_inc[1] - sp->r_inc[0]) + sp->r_inc[0];
    ipar++;
  }

  // EFAC
  if (sp->have_efac) {
    for (unsigned int j=0; j < 1; j++) {
      theta[ipar] = cube[ipar] * (sp->r_efac[1] - sp->r_efac[0]) + sp->r_efac[0];
      ipar++;
    }
  }

  // Psi0 jumps between different datasets                                                                                  
  if (sp->njump) {
    if(!sp->psi_jump_fixed) {
      for (unsigned  l=0; l<sp->njump; l++) {
	theta[ipar] = cube[ipar] *  (sp->r_psi_jump[1] - sp->r_psi_jump[0]) + sp->r_psi_jump[0];
	ipar++;
      }
    }
  }

}


double globalRVMLogLike_PC(double theta[], int nDims, double phi[], int nDerived, void *context) {

  //void globalRVMLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context) {

  //	MNStruct *par = ((MNStruct *)context);

	char label[16];	
	int ipar=0;
	double Qu,Uu;
	double lnew;

	sp->alpha = theta[ipar] * DEG_TO_RAD; ipar++;
	sp->delta = theta[ipar] * DEG_TO_RAD; ipar++;
	sp->phase0 = theta[ipar] * DEG_TO_RAD; ipar++;
	if (sp->sin_psi) {
	  sp->psi00 = atan(theta[ipar+1]/theta[ipar]);
	  ipar += 2;
	} else {
	  sp->psi00 = theta[ipar] * DEG_TO_RAD; ipar++;
	}

	// RVM center Phase of the profiles - Can be marginalized
	if (!sp->margin_phi0) {
	    for (unsigned int j = 0; j < sp->n_epoch; j++) {
		sp->phi0[j] = theta[ipar] * DEG_TO_RAD;
		ipar++;
	    }
	} else {
	    for (unsigned int j = 0; j < sp->n_epoch; j++) sp->phi0[j] = 5.0 * DEG_TO_RAD;
	}

	// include an offset between the PA curves of the MP and IP
	if (sp->have_aberr_offset) {
	    for (unsigned int j = 0; j < sp->nfiles_aberr; j++) {
		sp->phi_aberr_offset[j] = theta[ipar] * DEG_TO_RAD;
		ipar++;
	    }
	}

	// Precession rate
	if (!sp->prate_fixed) {
	    sp->prate = theta[ipar];
	    ipar++;
	}

	// Inclination angle
	if (!sp->inc_fixed) {
 	    sp->inc = theta[ipar] * DEG_TO_RAD;
	    ipar++;
	}
	
	// EFACs
	if (sp->have_efac) {
	    for (unsigned int j = 0; j <3; j++) {
	        sp->efac[j] = theta[ipar];
		ipar++;
	    }
	} else {
	    for (unsigned int j = 0; j < 3; j++) sp->efac[j] = 1.;
	}

	// Psi0 jumps between different datasets
	if (sp->njump) {
	  if(!sp->psi_jump_fixed) {
	    for (unsigned  l=0; l<sp->njump; l++) {
	      sp->psi_jumps[l] = theta[ipar] * DEG_TO_RAD;
	      ipar++;
	    }
	  }
	}

	double rmsQ;
	double rmsU;

	vector<double> vrmsQ, vrmsU;
	vector< vector<double> > vQ, vU, vphase;
	vQ.resize(sp->n_epoch); vU.resize(sp->n_epoch); vphase.resize(sp->n_epoch);

	for(unsigned int j = 0; j < sp->n_epoch; j++) {
	  vrmsQ.push_back(sp->rmsQ[j]);
	  vrmsU.push_back(sp->rmsU[j]);
	  for(unsigned int i = 0; i < sp->npts[j]; i++) {
	    vQ[j].push_back(sp->Q[j][i]);
	    vU[j].push_back(sp->U[j][i]);
	    vphase[j].push_back(sp->phase[j][i]);
	  }
	}

#ifdef HAVE_MINUIT
	// Marginalize over Phi0
	if (sp->margin_phi0) {
	    //RVM_Fcn fFCN();
	    RVM_Fcn fFCN(vphase, vQ, vU, vrmsQ, vrmsU);
	    fFCN.set_params(sp);

	    MnUserParameters upar;
	    for(unsigned int j = 0; j < sp->n_epoch; j++) {
		sprintf(label, "phi0_%d", j);
		upar.Add(label, sp->phi0[j], 5.);
	    }

	    MnMinimize migrad(fFCN, upar);
	    //cout << upar << endl;
	    FunctionMinimum min = migrad();
	    //cout<<"minimum: "<<min<<endl;

	    for(unsigned int j = 0; j < sp->n_epoch; j++) {
		sprintf(label, "phi0_%d", j);
		sp->phi0[j] = min.UserState().Value(label);
		//printf("phi0[%02d] = %lf\n",j, par->phi0[j]);
	    }
	}
#endif

	get_globalRVM_chi2(sp);

	// shift the reference P.A. by 90 degrees and ensure that PA0 lies on -pi/2 -> pi/2
	//if (sp->Ltot < 0.0) {
	  //sp->psi00 += M_PI /2.;
	    //}
	//sp->psi00 = atan( tan(sp->psi00) );
	//theta[3] = sp->psi00 / DEG_TO_RAD;

        lnew = -1.*sp->chi/2 - 0.5*sp->logdetN;
	return lnew;
}
