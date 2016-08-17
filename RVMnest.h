#include<vector>

#define NB_PTS_PROF 2048.0
#define TSUN 4.925490947e-06
#define r13 (1./3.)
#define r23 (2./3.)
#define r43 (4./3.)
#define r53 (5./3.)

typedef struct {

        int n_epoch;    /* Number of epochs */

        double *epoch;
        int *npts;

	int prate_fixed;
	int inc_fixed;
	int have_efac;
	int margin_phi0;

        double alpha;   /*  par[0] */
        //double beta;
        //double xsi;
        double delta;   /*  par[1] */
        //double lambda;   /* */
        double phase0;  /* */
        double inc;             /* Orbital inclination */
        double omega;
	double prate;
        double psi00;   /*  */
        double *phi0;  /* */
        double psi0;

  std::vector<std::vector<double> > phase;
  std::vector<std::vector<double> >  Q;
  std::vector<std::vector<double> >  U;

  std::vector<double> rmsQ;
  std::vector<double> rmsU;
	double *efac;

	double *r_alpha;
	double *r_delta;
	double *r_Phi0;
	double *r_phi0;
	double *r_inc;
	double *r_prate;
	double *r_efac;
  int njump;
  double *psi_jump_MJD;

} MNStruct;

void readPA(char filename[127], std::vector<double> &  x, std::vector<double> & Q, std::vector<double> & U);
MNStruct* init_struct(int nfiles, std::vector<std::vector<double> > x, std::vector<std::vector<double> > Q, std::vector<std::vector<double> > U, std::vector<double> rmsQ, std::vector<double> rmsU);
void RVMLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void globalRVMLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
double get_RVM(const double &al, const double &be, const double &ph0, const double &ps0, const double &x);
