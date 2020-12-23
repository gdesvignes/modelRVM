#include <complex>
#include "RVMnest.h"
using namespace std;

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/FCNBase.h"
using namespace ROOT::Minuit2;

class Pmargin : public FCNBase {

public:


  Pmargin(vector<vector<double> >  phase, vector<vector<double> >  Q, vector<vector<double> >  U, vector<double>  Q_err, vector<double> U_err)
  		: fPhase(phase), polQ(Q), polU(U), varQ(Q_err), varU(U_err), theErrorDef(1.) {}

  virtual ~Pmargin() {}

  virtual double operator()(const vector<double>& par) const {

    double chi2 = 0.;
    double Ln, PA,PA2;
    complex <double> L, arg;
    int totnbpts=0, ipar=0;
    double prev_chi = 0.0;
    double rms2Q, rms2U;

    // -- FIT PARAMS --
    if (params->margin_phi0) {
      for(unsigned int i = 0; i < params->n_epoch; i++) {
	params->phi0[i] = par[ipar]; ipar++;
      }
    }
    
    if (params->margin_psi0) {
      for(unsigned int i = 0; i < params->n_epoch; i++) {
	params->psi0[i] = par[ipar]; ipar++;
      }
    }
    

    //cout << params->margin_psi0<< " "<< params->n_epoch << " " << params->alpha<< " "<< params->beta[0]<< endl;
    for(unsigned int j = 0; j < params->n_epoch; j++) {
      rms2Q = params->rmsQ[j] * params->rmsQ[j];
      rms2U = params->rmsU[j] * params->rmsU[j];

      if (params->have_efac) {
            rms2Q *= params->efac[0]*params->efac[0];
            rms2U *= params->efac[0]*params->efac[0];
        }

        totnbpts += params->npts[j];
        params->Ltot[j] = 0;
        for(unsigned int i = 0; i < params->npts[j]; i++) {
            if (params->pmodel==0) // pmodel=0 forced precession, pmodel=1 free precession                                                                               
                PA2 = 2*RVMz(params->alpha, params->alpha+params->beta[j], params->phi0[j], params->psi0[j], params->phase[j][i]);
            else // in this case, params->alpha is really zeta                                                                                                           
                PA2 = 2*RVMz(params->alpha-params->beta[j], params->alpha+params->beta[j], params->phi0[j], params->psi0[j], params->phase[j][i]);

            Ln = (params->Q[j][i]*cos(PA2) / rms2Q + params->U[j][i]*sin(PA2) / rms2U)
                / (cos(PA2)*cos(PA2) / rms2Q + sin(PA2)*sin(PA2) / rms2U);

            arg = complex<double>(0., PA2);
            L =  Ln * exp (arg);
            chi2 += (params->Q[j][i]-real(L))*(params->Q[j][i]-real(L)) / rms2Q
	      + (params->U[j][i]-imag(L))*(params->U[j][i]-imag(L)) / rms2U;
        }
    }
    //cout << chi2 << endl;
    return chi2;
  }

  virtual double Up() const {return theErrorDef;}

  void setErrorDef(double def) {theErrorDef = def;}

  void set_params(MNStruct *p) {
      params = p;
  }

  void printQU() {
    for(int i=0; i<polQ[0].size(); ++i)
      cout << polQ[0][i] << " " << polU[0][i] << endl;
  }

private:

  // Data passed to the function
  vector<vector<double> > fPhase;
  vector<vector<double> > polQ;
  vector<vector<double> > polU;
  vector<double> varQ;
  vector<double> varU;

  double theErrorDef;

  MNStruct *params;

  double RVMz(const double &alpha, const double &xsi, const double &phi0, const double &psi0, const double &x) const {

      double k1 = sin(alpha)*sin(x-phi0);
      double k2 = sin(xsi)*cos(alpha)-cos(xsi)*sin(alpha)*cos(x-phi0);
      return -atan(k1/k2) + psi0;
  }

};
