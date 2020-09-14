#ifndef RVM_H_INCLUDED
#define RVM_H_INCLUDED

#include<vector>

#define NB_PTS_PROF 2048.0
#define TSUN 4.925490947e-06
#define r13 (1./3.)
#define r23 (2./3.)
#define r43 (4./3.)
#define r53 (5./3.)

typedef struct {
    
    int n_epoch;    /* Number of epochs */
    int ndims;
    double *epoch;
    int *npts;
    int *nbin;
    int sampler;
    int prate_fixed;
    int inc_fixed;
    int have_efac;
    int have_aberr_offset;
    int nfiles_aberr;
    int margin_phi0;
    int psi_jump_fixed;
    int do_plot;
    int sin_psi;
    
    double alpha;   /*  par[0] */
    double *beta;
    double delta;   /*  par[1] */
    double phase0;  /* */
    double inc;             /* Orbital inclination */
    double omega;
    double prate;
    double psi00;   /*  */
    double *phi0;  /* */
    double *psi0;
    double *phi_aberr_offset;
  
    // Offset dipole
    int have_offset_dipole;
    double phas;
    double Minc;
    double ita;
    double eps;

    // Free precession model
    double pperiod;
    double pphase;
    double pfact;
    double betaamp;
    double beta0;
    double mjd0;

    double **phase;
    double **I;
    double **Q;
    double **U;
    double **V;
    double **L;
    
    double *rmsI;
    double *rmsQ;
    double *rmsU;
    double *efac;
    
    double *r_alpha;
    double *r_beta;
    double *r_delta;
    double *r_Phi0;
    double *r_phi0;
    double *r_psi0;
    double *r_inc;
    double *r_prate;
    double *r_efac;
    int njump;
    double *r_psi_jump;
    double *psi_jumps;
    double *psi_jump_MJD;
    
    double chi;
    double *Ltot;
    double logdetN;
} MNStruct;

void readPA(char filename[127], std::vector<double> &  x, std::vector<double> & Q, std::vector<double> & U);
MNStruct* init_struct(int nfiles, std::vector<std::vector<double> > x, std::vector<std::vector<double> > I, std::vector<std::vector<double> > Q, std::vector<std::vector<double> > U, std::vector<std::vector<double> > L, std::vector<std::vector<double> > V, std::vector<double> rmsI, std::vector<double> rmsQ, std::vector<double> rmsU, std::vector<int>& nbin, int njump);
void RVMLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void globalRVMLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void precessLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void get_precessRVM_chi2(MNStruct *par);
void get_precessmRVM_chi2(MNStruct *par);
void get_globalRVM_chi2(MNStruct *par);
void get_RVM_chi2(MNStruct *par);
double get_RVM(const double &al, const double &be, const double &ph0, const double &ps0, const double &x);
double get_offRVM(const double &alpha, const double &beta, const double &phi0, const double &psi0, const double &x, const double &phas, const double &Minc, const double &ita, const double &eps);

#endif
