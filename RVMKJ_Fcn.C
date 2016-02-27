#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/MnPlot.h"

using namespace std;
using namespace ROOT::Minuit2;

#include "stdio.h"
#include <complex.h> 
#include <math.h> 


class RVM_Fcn : public FCNBase {


public:

  RVM_Fcn(const vector<double>& phase, const vector<double>& Q, const vector<double>& U, const double Q_err, const double U_err)
  		: fPhase(phase), polQ(Q), polU(U), rmsQ(Q_err), rmsU(U_err), theErrorDef(1.) {}

  virtual ~RVM_Fcn() {}

  virtual double operator()(const vector<double>& par) const {

    double RVM_val, delta;
    double chi2 = 0.;
    double Ln, Ltot;
    double PA;

    //printf("alpha=%lf beta=%lf\n", par[0], par[1]);
    complex<double> L;

    for(unsigned int i = 0; i < polQ.size(); i++) {
      PA = get_RVM(par, fPhase[i]);
      Ln = (polQ[i]*cos(2*PA) / (rmsQ*rmsQ) + polU[i]*sin(2*PA) / (rmsU*rmsU)) / (cos(2*PA)*cos(2*PA) / (rmsQ*rmsQ) + sin(2*PA)*sin(2*PA) / (rmsU*rmsU)); 
      //L =  par[4+i] * exp ( complex<double >(0., 2.*get_RVM(par, fPhase[i])));
      //L =  par[4+i] * exp (2. *1.i*get_RVM(par, fPhase[i]));
      //L =  Ln * cexp (2. *I*get_RVM(par, fPhase[i]));
      L =  Ln * exp ( complex<double >(0., 2.*get_RVM(par, fPhase[i])));

      Ltot += Ln;
      //printf("%d %lf %lf\n", i, (polQ[i]-real(L))*(polQ[i]-real(L))/(rmsQ*rmsQ), (polU[i]-imag(L))*(polU[i]-imag(L))/(rmsU*rmsU));
      chi2 += (polQ[i]-real(L))*(polQ[i]-real(L))/(rmsQ*rmsQ)  + (polU[i]-imag(L))*(polU[i]-imag(L))/(rmsU*rmsU);
    }
    return chi2;
  }

  virtual double Up() const {return theErrorDef;}

  void setErrorDef(double def) {theErrorDef = def;}

private:
  vector<double> fPhase;
  vector<double> polQ;
  vector<double> polU;
  double theErrorDef;

  double rmsQ, rmsU;

  double get_RVM(const vector<double>& par, const double &x) const {

    double alpha, beta, xsi, phi0, psi0;

    alpha = par[0]*M_PI/180.;
    beta = par[1]*M_PI/180.;
    phi0 = par[2]*M_PI/180.;
    psi0 = par[3]*M_PI/180.;
    xsi = alpha + beta;

    //alpha = par[0];
    //xsi = par[1];
    //phi0 = par[2];
    //psi0 = par[3];
    
    double k1 = sin(alpha)*sin(x-phi0);
    double k2 = sin(xsi)*cos(alpha)-cos(xsi)*sin(alpha)*cos(x-phi0);
    double PA = -atan(k1/k2) + psi0;
    return  PA;
  }
};
