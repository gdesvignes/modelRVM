#include "RVMnest.h"
#include <stdlib.h>
#include <math.h>

#include <iostream>

using namespace std;

MNStruct* init_struct(int nfiles, vector< vector<double> > phase, vector< vector<double> > Q, vector< vector<double> > U, vector<double> rmsQ, vector<double> rmsU) {

	// Init struct
	MNStruct* MNS = (MNStruct*)malloc(sizeof(MNStruct));
	
	// Malloc
	MNS->npts = (int *)malloc(nfiles * sizeof(int));
	MNS->efac = (double *)malloc(nfiles * sizeof(double));
	MNS->epoch = (double *)malloc(nfiles * sizeof(double));
	MNS->phi0 = (double *)malloc(nfiles * sizeof(double));

	MNS->rmsQ = (double *)malloc(nfiles * sizeof(double));
	MNS->rmsU = (double *)malloc(nfiles * sizeof(double));

	MNS->phase = (double **)malloc(nfiles * sizeof(double *));
	MNS->Q = (double **)malloc(nfiles * sizeof(double *));
	MNS->U = (double **)malloc(nfiles * sizeof(double *));

	MNS->n_epoch = nfiles;
	for (int i=0; i<nfiles; i++) {
		MNS->npts[i] = phase[i].size(); // TBD 
		cout << i << " "<<  nfiles << " "<<  MNS->npts[i] <<endl;
		MNS->phase[i] = (double *)malloc(MNS->npts[i] * sizeof(double));
		MNS->Q[i] = (double *)malloc(MNS->npts[i] * sizeof(double));
		MNS->U[i] = (double *)malloc(MNS->npts[i] * sizeof(double));
		MNS->rmsQ[i] = rmsQ[i];
		MNS->rmsU[i] = rmsU[i];
	}
		
	cout << "Malloc finished" << endl;

	// 1906 params
	MNS->a1 = 1.420; // lt-s
	MNS->pb = 0.1659930; // days
	MNS->ecc = 0.08530;
	MNS->omdot = 7.5841; // deg/yr
	MNS->massfn = 0.111568;

	// -- Tune the parameters --
	MNS->pb *= 86400.0; // seconds
	MNS->omdot *=  M_PI/180.0 / ( 86400. * 365.245 ); // rad/seconds
	
	// Assign data
	for (int i=0; i<nfiles; i++) {
	    for (int j=0; j<MNS->npts[i]; j++) {
		MNS->phase[i][j] = phase[i][j]; // TBD
		MNS->Q[i][j] = Q[i][j];
		MNS->U[i][j] = U[i][j];
	    }
	}

	cout << "Params finished" << endl;

	return MNS;
}

double get_RVM(const double &alpha, const double &xsi, const double &phi0, const double &psi0, const double &x) {

    //printf("%lf %lf %lf %lf %lf\n", alpha, xsi, phi0, psi0, x);

    double k1 = sin(alpha)*sin(x-phi0);
    double k2 = sin(xsi)*cos(alpha)-cos(xsi)*sin(alpha)*cos(x-phi0);
    double PA = -atan(k1/k2) + psi0;
    return  PA;
}

