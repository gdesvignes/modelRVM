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

        double chi = 0.0, Ln, PA;
	complex <double> L, arg;

	double alpha, beta, phi0, psi0;
	if (par->r_alpha==NULL) alpha = Cube[0] * M_PI;
	else alpha = Cube[0] * (par->r_alpha[1]*DEG_TO_RAD - par->r_alpha[0]*DEG_TO_RAD) + par->r_alpha[0]*DEG_TO_RAD;

	if (par->r_beta==NULL) beta = Cube[1] * M_PI - M_PI/2.;
	else beta = Cube[1] * (par->r_beta[1]*DEG_TO_RAD - par->r_beta[0]*DEG_TO_RAD) + par->r_beta[0]*DEG_TO_RAD;
	
	if (par->r_phi0==NULL) phi0 = Cube[2] * M_PI;
	else phi0 = Cube[2] * (par->r_phi0[1]*DEG_TO_RAD - par->r_phi0[0]*DEG_TO_RAD) + par->r_phi0[0]*DEG_TO_RAD;

	if (par->r_psi0==NULL) psi0 = Cube[3] * M_PI;
	else psi0 = Cube[3] * (par->r_psi0[1]*DEG_TO_RAD - par->r_psi0[0]*DEG_TO_RAD) + par->r_psi0[0]*DEG_TO_RAD;

	double rmsQ = par->rmsQ[0];
	double rmsU = par->rmsU[0];
	double Ltot = 0.0;

	for(unsigned int i = 0; i < par->npts[0]; i++) {
		PA = get_RVM(alpha, beta, phi0, psi0, par->phase[0][i]);
		Ln = (par->Q[0][i]*cos(2*PA) / (rmsQ*rmsQ) + par->U[0][i]*sin(2*PA) / (rmsU*rmsU)) 
			/ (cos(2*PA)*cos(2*PA) / (rmsQ*rmsQ) + sin(2*PA)*sin(2*PA) / (rmsU*rmsU));

		arg = complex<double>(0., 2 * PA);
		L =  Ln * exp (arg);
		chi += (par->Q[0][i]-real(L))*(par->Q[0][i]-real(L))/(rmsQ*rmsQ)
			+ (par->U[0][i]-imag(L))*(par->U[0][i]-imag(L))/(rmsU*rmsU);
		Ltot += Ln;
	}	

        // shift the reference P.A. by 90 degrees and ensure that PA0 lies on -pi/2 -> pi/2
        if (Ltot < 0.0) {
            psi0 += M_PI /2.;
        }
        psi0 = atan( tan(psi0) );

	Cube[0] = alpha / DEG_TO_RAD;
	Cube[1] = beta / DEG_TO_RAD;
	Cube[2] = phi0 / DEG_TO_RAD;
	Cube[3] = psi0 / DEG_TO_RAD;

        lnew = -chi/2;
}
