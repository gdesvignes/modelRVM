#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

double get_RVM_model(const double &x, const double &alpha, const double &beta, const double &phi0, const double &psi0);

//float naghi(float pa, float lintrue, float int_rms);

void write_results( char *filename, double mjd, int nbpts, float rms, float* totI,  float* totQ,  float* totU,  float* totV, vector<int> phase_bin, vector<double> Q, vector<double>U, double al, double be, double phi, double psi){

  string str(filename);
  string ofn = str + "-rvm.txt";

  ofstream myfile;
  myfile.open (ofn.c_str(), ofstream::out);

  float PA, PA_model;
  double L, sigPA;

  myfile << "# " << mjd << endl;

  for (int i=0; i< nbpts; i++) {

      L = sqrt(totQ[i]*totQ[i]+totU[i]*totU[i]);


      for(int j=0; j<phase_bin.size(); j++) {
	//printf("i=%d j=%d\n", i, phase_bin[j]);
	if (phase_bin[j]==i) {
	    PA = 1/2. * atan2(U[j],Q[j]) * 180. / M_PI;
	    sigPA = 28.65 * rms/L;
	    //if (L>10*rms) sigPA = 28.65 * rms/L;
	    //else sigPA = naghi(PA, L, rms);
	    break;
	}
	else {
	    PA = 0;
	    sigPA = 0.;
	}	
      }
      if (PA < 0.) PA+= 180;

      PA_model = get_RVM_model((i+.5)/2048.*2*M_PI, al*M_PI/180., be*M_PI/180., phi*M_PI/180., psi*M_PI/180.);

      myfile << i << " " << totI[i] << " " << L << " " << totV[i] << " " << PA << " " << sigPA << " " << PA_model * 180./M_PI + 180<< " " << endl;
  }
  myfile.close();


}

double get_RVM_model(const double &x, const double &alpha, const double &beta, const double &phi0, const double &psi0) {
      double xsi;

      xsi = alpha + beta;
      double k1 = sin(alpha)*sin(x-phi0);
      double k2 = sin(xsi)*cos(alpha)-cos(xsi)*sin(alpha)*cos(x-phi0);
      return -atan(k1/k2) + psi0;
  }

