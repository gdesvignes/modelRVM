#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "RVMnest.h"
#include <complex>

#define DEG_TO_RAD (M_PI/180.0)

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

void RVMLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context) {

	MNStruct *par = ((MNStruct *)context);


	if (par->r_alpha==NULL) par->alpha = Cube[0] * M_PI;
	else par->alpha = Cube[0] * (par->r_alpha[1]*DEG_TO_RAD - par->r_alpha[0]*DEG_TO_RAD) + par->r_alpha[0]*DEG_TO_RAD;

	if (par->r_beta==NULL) par->beta = Cube[1] * M_PI - M_PI/2.;
	else par->beta = Cube[1] * (par->r_beta[1]*DEG_TO_RAD - par->r_beta[0]*DEG_TO_RAD) + par->r_beta[0]*DEG_TO_RAD;
	
	if (par->r_phi0==NULL) par->phi0[0] = Cube[2] * M_PI;
	else par->phi0[0] = Cube[2] * (par->r_phi0[1]*DEG_TO_RAD - par->r_phi0[0]*DEG_TO_RAD) + par->r_phi0[0]*DEG_TO_RAD;

	if (par->r_psi0==NULL) par->psi0 = Cube[3] * M_PI;
	else par->psi0 = Cube[3] * (par->r_psi0[1]*DEG_TO_RAD - par->r_psi0[0]*DEG_TO_RAD) + par->r_psi0[0]*DEG_TO_RAD;

	if (par->have_efac) par->efac[0] = Cube[4] * (par->r_efac[1]-par->r_efac[0]) + par->r_efac[0];
	else par->efac[0] = 1.0;

        get_RVM_chi2(par);

	Cube[0] = par->alpha / DEG_TO_RAD;
	Cube[1] = par->beta / DEG_TO_RAD;
	Cube[2] = par->phi0[0] / DEG_TO_RAD;
	Cube[3] = par->psi0 / DEG_TO_RAD;
	if (par->have_efac) Cube[4] = par->efac[0];

        lnew = -par->chi/2 - 0.5*par->logdetN;;
}
