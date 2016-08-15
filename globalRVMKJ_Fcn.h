#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/FCNBase.h"

using namespace std;
using namespace ROOT::Minuit2;

#include "RVMnest.h"
//#include "globalRVM_Fcn.h"
#include <complex> 


class RVM_Fcn : public FCNBase {


public:

  //RVM_Fcn() : theErrorDef(1.) {}

  RVM_Fcn(vector<vector<double> >  phase, vector<vector<double> >  Q, vector<vector<double> >  U, vector<double>  Q_err, vector<double> U_err)
  		: fPhase(phase), polQ(Q), polU(U), varQ(Q_err), varU(U_err), theErrorDef(1.) {}

  virtual ~RVM_Fcn() {}

  virtual double operator()(const vector<double>& par) const {

      double chi2 = 0.;
      double phi, csdel, sndel, csi, sni, csphi, snphi, cslam, snlam;
      double eta, cseta, sneta, beta, lambda, xsi;
      double Qu, Uu, dt;
      double PA2, cosPA2, sinPA2, Ln;
      //double logdetN = 0;
      complex<double> L;
      complex<double> arg;

      // -- FIT PARAMS --
      for(unsigned int i = 0; i < params->n_epoch; i++) {
          params->phi0[i] = par[i];
      }

      for(unsigned int i = 0; i < params->n_epoch; i++) {
          dt = params->epoch[i] - params->epoch[0];

          // -- compute beta and psi first --
          phi = params->phase0 + params->omega * dt; // in rad

          csdel = cos(params->delta);
          sndel = sin(params->delta);

          csi = cos(params->inc);
          sni = sin(params->inc);

          csphi = cos(phi);
          snphi = sin(phi);

          cslam = csdel * csi - sndel*sni*csphi;
          snlam = sqrt(1.0-cslam*cslam);

          //  as 0 < lambda < 180, sin(lambda) >0 
          // hence we can just compute acos(cos(lambda))

          lambda = acos(cslam);
          beta = M_PI - params->alpha - lambda;

          cseta = sndel*snphi/snlam;
          sneta = (cslam*csi-csdel)/sni/snlam;

          eta = atan2(sneta,cseta);

	  xsi = params->alpha + beta;


	  Qu = varQ[i] * params->efac[i]*params->efac[i];
	  Uu = varU[i] * params->efac[i]*params->efac[i];
	 
	  //cout << dt << " " << params->alpha * 180./M_PI << " " << beta * 180./M_PI << endl; 

          for(unsigned int j = 0; j < params->npts[i]; j++) {
              PA2 = 2 * get_RVM(params->alpha, xsi, params->phi0[i], params->psi00 + eta, params->phase[i][j]);

	      cosPA2 = cos(PA2);
	      sinPA2 = sin(PA2);

	      Ln = (polQ[i][j]*cosPA2 / Qu + polU[i][j]*sinPA2 / Uu) / (cosPA2*cosPA2 / Qu + sinPA2*sinPA2 / Uu);

	      arg = std::complex<double>(0., PA2);
	      L =  Ln * exp (arg);

	      chi2 += (polQ[i][j]-real(L))*(polQ[i][j]-real(L))/Qu +
	      		(polU[i][j]-imag(L))*(polU[i][j]-imag(L))/Uu;

	      //logdetN += log(Uu) + log(Qu);

	      //printf("i=%d j=%d Ln=%lf ipar=%lf varQ=%lf varU=%lf \n", i, j, Ln, par[6 + fPhase.size() + tot_param + j], varQ[i], varU[i]);


          }
      }
      return chi2;
  }

  virtual double Up() const {return theErrorDef;}

  void setErrorDef(double def) {theErrorDef = def;}

  void set_params(MNStruct *p) {
      params = p;
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

  double get_RVM(const double &alpha, const double &xsi, const double &phi0, const double &psi0, const double &x) const {

      double k1 = sin(alpha)*sin(x-phi0);
      double k2 = sin(xsi)*cos(alpha)-cos(xsi)*sin(alpha)*cos(x-phi0);
      return -atan(k1/k2) + psi0;
  }
};
