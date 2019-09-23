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

	int ipar = 0;
	if (par->r_alpha==NULL) par->alpha = Cube[ipar] * M_PI;
	else par->alpha = Cube[ipar] * (par->r_alpha[1]*DEG_TO_RAD - par->r_alpha[0]*DEG_TO_RAD) + par->r_alpha[0]*DEG_TO_RAD;
	ipar++;

	if (par->r_beta==NULL) par->beta = Cube[ipar] * M_PI - M_PI/2.;
	else par->beta = Cube[ipar] * (par->r_beta[1]*DEG_TO_RAD - par->r_beta[0]*DEG_TO_RAD) + par->r_beta[0]*DEG_TO_RAD;
	ipar++;

	if (par->r_phi0==NULL) par->phi0[0] = Cube[ipar] * M_PI;
	else par->phi0[0] = Cube[ipar] * (par->r_phi0[1]*DEG_TO_RAD - par->r_phi0[0]*DEG_TO_RAD) + par->r_phi0[0]*DEG_TO_RAD;
	ipar++;

	if (par->r_psi0==NULL) par->psi0 = Cube[ipar] * M_PI;
	else par->psi0 = Cube[ipar] * (par->r_psi0[1]*DEG_TO_RAD - par->r_psi0[0]*DEG_TO_RAD) + par->r_psi0[0]*DEG_TO_RAD;
	ipar++;

	if (par->have_efac) {par->efac[0] = Cube[ipar] * (par->r_efac[1]-par->r_efac[0]) + par->r_efac[0]; ipar++;}
	else par->efac[0] = 1.0;

	if  (par->have_aberr_offset) {par->phi_aberr_offset[0] = (Cube[ipar] * 40-20)*DEG_TO_RAD ;ipar++;}
	else par->phi_aberr_offset[0] = 0.0;

	if  (par->have_offset_dipole) {
	  par->phas = Cube[ipar] * M_PI; ipar++;
	  par->Minc = Cube[ipar] * M_PI; ipar++;
	  par->ita = Cube[ipar] * 1; ipar++; // emission altitude, h = (1+ita)*R,  up to 500 R above the magnetosphere
	  par->eps = Cube[ipar]; ipar++; // ok, between 0 and 1
	} else {par->eps=0;}

        get_RVM_chi2(par);

	Cube[0] = par->alpha / DEG_TO_RAD;
	Cube[1] = par->beta / DEG_TO_RAD;
	Cube[2] = par->phi0[0] / DEG_TO_RAD;
	Cube[3] = par->psi0 / DEG_TO_RAD;
	ipar=4;
	if (par->have_efac) {Cube[ipar] = par->efac[0];ipar++;}
	if (par->have_aberr_offset) {Cube[ipar] = par->phi_aberr_offset[0] / DEG_TO_RAD;ipar++;}

	if  (par->have_offset_dipole) {
	  Cube[ipar] = par->phas / DEG_TO_RAD; ipar++;
	  Cube[ipar] = par->Minc / DEG_TO_RAD; ipar++;
	  Cube[ipar] = par->ita; ipar++;
	  Cube[ipar] = par->eps; ipar++;
	}

        lnew = -par->chi/2 - 0.5*par->logdetN;;
}
