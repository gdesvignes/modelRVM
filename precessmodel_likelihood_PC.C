#include <config.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex>

#include "precessm_likelihood_PC.h"
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

void precessmprior (double cube[], double theta[], int nDims, void *context) {
    int ipar = 0;

    // Alpha
    //theta[ipar] = cube[ipar] * (sp->r_alpha[1] - sp->r_alpha[0]) + sp->r_alpha[0];
    //ipar++;

    // alpha
    theta[ipar] = cube[ipar] * 180; ipar++;

    // First coeff
    theta[ipar] = cube[ipar] * 2 - 1; ipar++;

    // Second coeff
    theta[ipar] = cube[ipar] * 2 - 1; ipar++;

    // Third coeff
    //theta[ipar] = cube[ipar] * 2 - 1; ipar++;


    for (unsigned int j = 0; j < sp->n_epoch; j++) {
	theta[ipar] = cube[ipar] * (sp->r_beta[1] - sp->r_beta[0]) + sp->r_beta[0]; ipar++;
	theta[ipar] = cube[ipar] * (sp->r_phi0[1] - sp->r_phi0[0]) + sp->r_phi0[0]; ipar++;
	theta[ipar] = cube[ipar] * (sp->r_psi0[1] - sp->r_psi0[0]) + sp->r_psi0[0]; ipar++;
    }
    
    // EFAC
    if (sp->have_efac) {
        for (unsigned int j = 0; j < 1; j++) {
            theta[ipar] = cube[ipar] * (sp->r_efac[1]-sp->r_efac[0]) + sp->r_efac[0]; ipar++;
        }
    }

}

double precessmLogLike_PC(double theta[], int nDims, double phi[], int nDerived, void *context) {

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
    sp->mjd0 = theta[ipar]; ipar++;
    sp->pperiod = theta[ipar]; ipar++;
    //sp->pphase = theta[ipar]; ipar++;

    for (unsigned int j = 0; j < sp->n_epoch; j++) {
	sp->beta[j] = theta[ipar] * DEG_TO_RAD; ipar++;
	sp->phi0[j] = theta[ipar] * DEG_TO_RAD; ipar++;
	sp->psi0[j] = theta[ipar] * DEG_TO_RAD; ipar++;
    }

    // EFAC
    if (sp->have_efac) {
	for (unsigned int j = 0; j < 1; j++) {
            sp->efac[0] = theta[ipar]; ipar++;
        }
    }
    
    // Compute chi**2
    get_precessmRVM_chi2(sp);

    for (unsigned int j = 0; j < sp->n_epoch; j++) phi[j] = sp->Ltot[j];
    
    //lnew = -1.*par->chi/2 - 0.5*par->logdetN;
    return -1.*sp->chi/2. - 0.5*sp->logdetN;
}
