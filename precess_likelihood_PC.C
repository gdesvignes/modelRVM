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
#include "precess_margin.h"
#endif

#include "precess_likelihood_PC.h"
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

void precessprior (double cube[], double theta[], int nDims, void *context) {
    int ipar = 0;
     
    theta[ipar] = cube[ipar] * (sp->r_alpha[1] - sp->r_alpha[0]) + sp->r_alpha[0];
    ipar++;

    if (sp->pmodel==2 || sp->pmodel==3) {
      theta[ipar] = cube[ipar] * (sp->r_alpha1[1] - sp->r_alpha1[0]) + sp->r_alpha1[0];
      ipar++;
      theta[ipar] = cube[ipar] * (sp->r_alpha2[1] - sp->r_alpha2[0]) + sp->r_alpha2[0];
      ipar++;
    }
    
    for (unsigned int j = 0; j < sp->n_epoch; j++) {
	//if (sp->epoch[j]<58570) theta[ipar]= cube[ipar] *(30-0) + 0;
	//else theta[ipar] = cube[ipar] * (0 + 30) - 30;
	//ipar++;
	theta[ipar] = cube[ipar] * (sp->r_beta[1] - sp->r_beta[0]) + sp->r_beta[0]; ipar++;
	if (!sp->margin_phi0) {theta[ipar] = cube[ipar] * (sp->r_phi0[1] - sp->r_phi0[0]) + sp->r_phi0[0]; ipar++; }
	if (!sp->margin_psi0) {theta[ipar] = cube[ipar] * (sp->r_psi0[1] - sp->r_psi0[0]) + sp->r_psi0[0]; ipar++; }
     }

    // EFACs                                                                                                                                        
    if (sp->have_efac) {
	for (unsigned int j = 0; j < 1; j++) {
            theta[ipar] = cube[ipar] * (sp->r_efac[1]-sp->r_efac[0]) + sp->r_efac[0]; ipar++;
	}
    }

    //theta[0] = 100;
    //theta[1] = 3 ;
    //theta[2] = 180 ;
    //theta[3] = 45;
}

double precessLogLike_PC(double theta[], int nDims, double phi[], int nDerived, void *context) {

    //MNStruct *par = ((MNStruct *)context);
    
    char label[16];	
    double Qu,Uu;
    double chi = 0.0, Ln, PA2;
    double Ltot = 0.;
    double cosPA2, sinPA2;
    complex <double> L, arg;
    double M, cslam, snlam, cseta, sneta, csi, sni, csphi, snphi, csdel, sndel;
    double lambda, beta;
    double dt;
    
    int ipar=0;
    sp->alpha = theta[ipar] * DEG_TO_RAD; ipar++;
    if (sp->pmodel==2 || sp->pmodel==3) {
      sp->alpha1 = theta[ipar] * DEG_TO_RAD; ipar++;
      sp->alpha2 = theta[ipar] * DEG_TO_RAD; ipar++;
    }
    for (unsigned int j = 0; j < sp->n_epoch; j++) {
	sp->beta[j] = theta[ipar] * DEG_TO_RAD; ipar++;
	if (sp->margin_phi0) {
	  sp->phi0[j] = (sp->r_phi0[0] + sp->r_phi0[1])/2. * DEG_TO_RAD;
	} else {
	  sp->phi0[j] = theta[ipar] * DEG_TO_RAD; ipar++;
	}
	if (sp->margin_psi0) {
	  sp->psi0[j] = (sp->r_psi0[0] + sp->r_psi0[1])/2. * DEG_TO_RAD;
	} else {
	  sp->psi0[j] = theta[ipar] * DEG_TO_RAD; ipar++;
	}
    }

    // EFACs
    if (sp->have_efac) {
        for (unsigned int j = 0; j < 1; j++) {
            sp->efac[0] = theta[ipar]; ipar++;
        }
    }

#ifdef HAVE_MINUIT
    // Init arrays to be passed to Minuit
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
    
        if (sp->margin_phi0 || sp->margin_psi0) {
            //RVM_Fcn fFCN();                                                                                                                                         
            Pmargin fFCN(vphase, vQ, vU, vrmsQ, vrmsU);
            fFCN.set_params(sp);

            MnUserParameters upar;
	    if (sp->margin_phi0) {
	      for(unsigned int j = 0; j < sp->n_epoch; j++) {
		sprintf(label, "phi0_%d", j);
                upar.Add(label, sp->phi0[j], 1., sp->r_phi0[0] * DEG_TO_RAD, sp->r_phi0[1] * DEG_TO_RAD);
	      }	    }
	    if (sp->margin_psi0) {
	      for(unsigned int j = 0; j < sp->n_epoch; j++) {
                sprintf(label, "psi0_%d", j);
                upar.Add(label, sp->psi0[j], 1., sp->r_psi0[0] * DEG_TO_RAD, sp->r_psi0[1] * DEG_TO_RAD);
              }
            }

	    MnMinimize migrad(fFCN, upar);
	    //cout << migrad.NFcn() << endl;                                                                                                                                   
	    FunctionMinimum min = migrad(2., 1e-2);
	    //cout << min.NFcn() << endl;  
            //cout<<"minimum: "<<min<<endl;
	    //migrad.PrintResults();
	    
	    if (sp->margin_phi0) {
	      for(unsigned int j = 0; j < sp->n_epoch; j++) {
		sprintf(label, "phi0_%d", j);
	        sp->phi0[j] = min.UserState().Value(label);
	      }
	    }
	    
	    if (sp->margin_psi0) {
	      for(unsigned int j = 0; j < sp->n_epoch; j++) {
		sprintf(label, "psi0_%d", j);
                sp->psi0[j] = min.UserState().Value(label);
	      }
            }
        }
#endif

    
    
    // Compute chi**2
    get_precessRVM_chi2(sp);

    for (unsigned int j = 0; j < sp->n_epoch; j++) phi[j] = sp->Ltot[j];
    
    //lnew = -1.*par->chi/2 - 0.5*par->logdetN;
    return -1.*sp->chi/2. - 0.5*sp->logdetN;
}
