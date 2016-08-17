#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/FCNBase.h"


//#include "RVMnest.h"
#include "globalRVMKJ_Fcn.h"
#include <complex>

#define DEG_TO_RAD	(M_PI/180.0)
#define RAD_TO_DEG	(180.0/M_PI)


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
	    for (unsigned int j = 0; j < par->n_epoch; j++) par->phi0[j] = 5.0 * M_PI/180.;
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

	double *psi_jumps;
	psi_jumps = (double *) malloc(par->njump * sizeof(double));
	if (par->njump) {
	  for (unsigned  l=0; l<par->njump; l++) {
	    psi_jumps[l] = Cube[npar] * 180.;
	    npar++;
	  }
	}

	double rmsQ;
	double rmsU;

	// Test check
#if 0
	par->alpha = 99.0484141926 * M_PI/180.;
	par->delta = 117.5365848921 * M_PI/180.;
	par->phase0 = 230.8024672783 * M_PI/180.;
	par->psi00 = 53.97486341145 * M_PI/180.;
	par->phi0[0] = 18.96421428197 * M_PI/180.;
	par->phi0[1] = 19.39861029502 * M_PI/180.;
	par->phi0[2] = 19.60041431936 * M_PI/180.;
	par->phi0[3] = 19.92786958068 * M_PI/180.;
	par->phi0[4] = 20.19144456874 * M_PI/180.;
	par->phi0[5] = 20.61378549792 * M_PI/180.;
	par->phi0[6] = 20.82955408938 * M_PI/180.;
	par->phi0[7] = 21.17505866282 * M_PI/180.;
	par->phi0[8] = 21.32304445259 * M_PI/180.;
	par->phi0[9] = 22.20055547289 * M_PI/180.;
	par->phi0[10] = 22.34049877568 * M_PI/180.;
	par->phi0[11] = 22.96195389531 * M_PI/180.;
	par->phi0[12] = 24.17391776309 * M_PI/180.;
	par->margin_phi0 = 1;
#endif
	//sni = sin(inc);
	//M = par->omdot * (1. - par->ecc*par->ecc) * r13 / (pow(TSUN,r23)) * pow((par->pb/(2*M_PI)),r53);
	//M = pow(M, 1.5);
	//exit(0);

	par->omega = par->prate / 365.25 * DEG_TO_RAD;

	/*
	par->mc = pow((par->massfn*M*M),r13) / sini;
	par->mp = M-par->mc;

	par->omega = pow((2.0*M_PI/par->pb),r53) * pow(TSUN,r23) / (1.-par->ecc*par->ecc);
	par->omega = par->omega * par->mc*(4.*par->mp+3.0*par->mc)/2.0/pow(M,r43);   // rad/sec
	par->omega = par->omega * 86400.; // rad/day
	*/

	// Marginalize over Phi0
	if (par->margin_phi0) {
	    //RVM_Fcn fFCN();
	    RVM_Fcn fFCN(par->phase, par->Q, par->U, par->rmsQ, par->rmsU);
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
	//exit(0);

	for (unsigned int j = 0; j < par->n_epoch; j++) {

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

	    //printf("delta=%lf  inc=%lf  phi=%lf\n", par->delta * 180./M_PI, par->inc * 180./M_PI, phi);
	    //printf("%2d  epoch=%lf  lambda=%lf  beta=%lf  eta=%lf\n", j, par->epoch[j], lambda * 180./M_PI, beta * 180./M_PI, eta * 180./M_PI);

	    rmsQ = par->rmsQ[j];
	    rmsU = par->rmsU[j];

	    xsi = par->alpha + beta;
	    
	    Qu = rmsQ*rmsQ*par->efac[j]*par->efac[j];
	    Uu = rmsU*rmsU*par->efac[j]*par->efac[j];

	    for(unsigned int i = 0; i < par->npts[j]; i++) {
		//printf("%lf %lf %lf %lf  %lf %lf  %lf\n", alpha, xsi, phi0, psi0, par->Q[0][i], par->U[0][i], par->x[0][i]);
		PA2 = 2*get_RVM(par->alpha, xsi, par->phi0[j], par->psi00 + eta, par->phase[j][i]);

		// Take into account Psi Jumps
		for(unsigned l=0; l<par->njump; l++) {
		  if (par->epoch[j] > par->psi_jump_MJD[l]) PA2 += psi_jumps[l] ;
		}

		cosPA2 = cos(PA2);
		sinPA2 = sin(PA2);

		Ln = (par->Q[j][i]*cosPA2 / (Qu) + par->U[j][i]*sinPA2 / (Uu)) 
			/ (cosPA2*cosPA2 / (Qu) + sinPA2*sinPA2 / (Uu));

		//arg = 1i * PA2;
		arg = std::complex<double>(0., PA2);
		L =  Ln * exp (arg);
		chi += (par->Q[j][i]-real(L))*(par->Q[j][i]-real(L))/(Qu)
			+ (par->U[j][i]-imag(L))*(par->U[j][i]-imag(L))/(Uu);

		logdetN += log(Uu) + log(Qu);
		
		Ltot += Ln;
	        //printf("%d %lf %lf\n", i,par->phase[0][i], PA);
	    }	
	}
	//printf("chi = %lf\n", chi);
	//exit(0);

	// shift the reference P.A. by 90 degrees and ensure that PA0 lies on -pi/2 -> pi/2
	if (Ltot < 0.0) {
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
	    Cube[npar] = psi_jumps[l];
            npar++;
          }
        }

	//printf("%d ", npar);
	//for (unsigned int j = 0; j < npar; j++) printf("%lf ", Cube[j]);
	//printf("\n");
        lnew = -chi/2 - 0.5*logdetN;
	free(psi_jumps);
}
