#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "RVMnest.h"
#include <complex>

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

	double alpha = Cube[0];
	double xsi = Cube[1];
	double phi0 = Cube[2];
	double psi0 = Cube[3];

	//double *Q = &par->Q[0];
	//double *U = &par->U[0];
	double rmsQ = par->rmsQ[0];
	double rmsU = par->rmsU[0];

	alpha *= M_PI;
	xsi *= M_PI;
	//phi0 *=  2*M_PI;
	phi0 =  phi0*M_PI + M_PI/2. ;
	psi0 = psi0*M_PI - M_PI/2.;
/*
	alpha =  98.9205  * M_PI / 180;
	xsi =  -5.7105  * M_PI / 180;
	phi0 =  18.9991  * M_PI / 180;
	psi0 =  -81.1533  * M_PI / 180;

	alpha =  81.0831  * M_PI / 180;
	xsi =  12.1269  * M_PI / 180;
	phi0 =  198.998  * M_PI / 180;
	psi0 =  -81.1622  * M_PI / 180;
*/
	//alpha =  80  * M_PI / 180;
	//xsi =  90  * M_PI / 180;
	//phi0 =  180  * M_PI / 180;
	//psi0 =  20  * M_PI / 180;

	double Ltot = 0.0;

	for(unsigned int i = 0; i < par->npts[0]; i++) {
		//printf("%lf %lf %lf %lf  %lf %lf  %lf\n", alpha, xsi, phi0, psi0, par->Q[0][i], par->U[0][i], par->x[0][i]);
		PA = get_RVM(alpha, xsi, phi0, psi0, par->phase[0][i]);
		Ln = (par->Q[0][i]*cos(2*PA) / (rmsQ*rmsQ) + par->U[0][i]*sin(2*PA) / (rmsU*rmsU)) 
			/ (cos(2*PA)*cos(2*PA) / (rmsQ*rmsQ) + sin(2*PA)*sin(2*PA) / (rmsU*rmsU));

		arg = complex<double>(0., 2 * PA);
		L =  Ln * exp (arg);
		chi += (par->Q[0][i]-real(L))*(par->Q[0][i]-real(L))/(rmsQ*rmsQ)
			+ (par->U[0][i]-imag(L))*(par->U[0][i]-imag(L))/(rmsU*rmsU);
		Ltot += Ln;
	        //printf("%d %lf %lf\n", i,par->x[0][i], PA);
	}	

        // shift the reference P.A. by 90 degrees and ensure that PA0 lies on -pi/2 -> pi/2
        if (Ltot < 0.0) {
            psi0 += M_PI /2.;
        }
        psi0 = atan( tan(psi0) );

	
	Cube[0] *= 180;
	Cube[1] *= 180;
	//Cube[2] *= 2*180.;
	Cube[2] = Cube[2]*180. + 90.;
	Cube[3] = psi0 * 180. / M_PI;

        lnew = -chi/2;
}
