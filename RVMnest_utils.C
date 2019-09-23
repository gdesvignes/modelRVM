#include "RVMnest.h"
#include <stdlib.h>
#include <math.h>

#include <iostream>

using namespace std;

MNStruct* init_struct(int nfiles, vector< vector<double> > phase, vector< vector<double> > I, vector< vector<double> > Q, vector< vector<double> > U, vector< vector<double> > L, vector< vector<double> > V, vector<double> rmsI, vector<double> rmsQ, vector<double> rmsU, vector<int> &nbin, int njump) {

	// Init struct
	MNStruct* MNS = (MNStruct*)malloc(sizeof(MNStruct));
	
	// Malloc
	MNS->npts = (int *)malloc(nfiles * sizeof(int));
	MNS->efac = (double *)malloc(nfiles * sizeof(double));
	MNS->epoch = (double *)malloc(nfiles * sizeof(double));
	MNS->beta = (double *)malloc(nfiles * sizeof(double));
	MNS->phi_aberr_offset = (double *)malloc(nfiles * sizeof(double));
	MNS->phi0 = (double *)malloc(nfiles * sizeof(double));
	MNS->psi0 = (double *)malloc(nfiles * sizeof(double));

	MNS->Q = (double **)malloc(nfiles * sizeof(double *));
	MNS->U = (double **)malloc(nfiles * sizeof(double *));
	MNS->phase = (double **)malloc(nfiles * sizeof(double *));
	MNS->rmsQ = (double *)malloc(nfiles * sizeof(double));
	MNS->rmsU = (double *)malloc(nfiles * sizeof(double));

	MNS->I = (double **)malloc(nfiles * sizeof(double *));
	MNS->L = (double **)malloc(nfiles * sizeof(double *));
	MNS->V = (double **)malloc(nfiles * sizeof(double *));
	MNS->rmsI = (double *)malloc(nfiles * sizeof(double));
	MNS->nbin = (int *)malloc(nfiles * sizeof(int));

	MNS->psi_jumps = (double *) malloc(njump * sizeof(double));
	
	
	MNS->n_epoch = nfiles;
	for (int i=0; i<nfiles; i++) {
	  MNS->efac[i] = 1.;
		MNS->npts[i] = phase[i].size(); // TBD 
		MNS->rmsQ[i] = rmsQ[i];
		MNS->rmsU[i] = rmsU[i];
		MNS->rmsI[i] = rmsI[i];
		MNS->nbin[i] = nbin[i];
		MNS->phase[i] = (double *)malloc(MNS->npts[i] * sizeof(double));
		MNS->Q[i] = (double *)malloc(MNS->npts[i] * sizeof(double));
		MNS->U[i] = (double *)malloc(MNS->npts[i] * sizeof(double));

		MNS->I[i] = (double *)malloc(MNS->nbin[i] * sizeof(double));
		MNS->L[i] = (double *)malloc(MNS->nbin[i] * sizeof(double));
		MNS->V[i] = (double *)malloc(MNS->nbin[i] * sizeof(double));


		for (int j=0; j<MNS->npts[i]; j++) {
		  MNS->phase[i][j] = phase[i][j];
		  MNS->Q[i][j] = Q[i][j];
		  MNS->U[i][j] = U[i][j];

		}
		for (int j=0; j<MNS->nbin[i]; j++) {
		  MNS->I[i][j] = I[i][j];
		  MNS->L[i][j] = L[i][j];
		  MNS->V[i][j] = V[i][j];
		}
       	}
	//cout << "Malloc finished" << endl;
	//cout << "Params finished" << endl;
	return MNS;
}

double get_RVM(const double &alpha, const double &beta, const double &phi0, const double &psi0, const double &x) {

    //printf("%lf %lf %lf %lf %lf\n", alpha, beta, phi0, psi0, x);

    double k1 = sin(alpha)*sin(x-phi0);
    double k2 = sin(alpha+beta)*cos(alpha)-cos(alpha+beta)*sin(alpha)*cos(x-phi0);
    double PA = -atan(k1/k2) + psi0;
    return  PA;
}

double get_offRVM(const double &alpha, const double &beta, const double &phi0, const double &psi0, const double &x, const double &phas, const double &Minc, const double &ita, const double &eps) {
  // phas is beta, phase de rotation du pulsar
  // Minc is delta, angle entre l'axe de rotation et OM
  // ita: altitude d'emission dans la magnetosphere
  // Zita: angle d'inclinaison de la ligne de visee
  // eps: distance par rapport au centre de l'etoile [0,1]
  double zita = alpha + beta;
  double k1 = sin(alpha) * (1+ita-eps*cos(Minc)*cos(zita)) * sin(phas + x-phi0) + eps*sin(Minc)* (cos(alpha)*cos(zita)*sin(x-phi0) - sin(alpha)*sin(phas)*sin(zita));
  double k2 = (1+ita)*(cos(alpha)*sin(zita)-sin(alpha)*cos(zita)*cos(phas + x-phi0)) + eps*(sin(alpha)*cos(Minc)*cos(phas + x-phi0) - cos(alpha)*sin(Minc)*cos(x-phi0));

  double PA = -atan(k1/k2) + psi0;
  return PA;

}
