/**********************************************************************
 *                                                                    *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MnContours.h"

#include "Pulsar/Archive.h"
#include "Pulsar/Integration.h"
#include "Pulsar/Profile.h"
//#include "Pulsar/Interpreter.h"
#include "Pulsar/FaradayRotation.h"
#include "Pulsar/PolnProfileStats.h"

#include "RVMKJ_Fcn.C"
#include "write_results.h"

#include "getopt.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"

#define NB_PTS_PROF 1024.0
#define R2D (180.0/M_PI)

using namespace std;
using namespace ROOT::Minuit2;
using namespace Pulsar;


int main(int argc, char **argv) {

  int nfree = 0;
  double threshold=3;
  double alpha=80.0, beta=6.0, phi0=18.0, psi0=30.0;
  bool alpha_fixed=false, beta_fixed=false, phi0_fixed=false, psi0_fixed=false;
  char filename[128]="test.pa";
  char label[16];

  static struct option long_opts[] = {
      {"thresh",  1, NULL, 't'},
      {"alpha",  1, NULL, 'a'},
      {"beta",   1, NULL, 'b'},
      {"phi0",   1, NULL, 'p'},
      {"psi0",   1, NULL, 'z'},
      {"Alpha",  0, NULL, 'A'},
      {"Beta",   0, NULL, 'B'},
      {"Phi0",   0, NULL, 'P'},
      {"Psi0",   0, NULL, 'Z'},
      {"help",   0, NULL, 'h'},
      {"file",   0, NULL, 'f'},
      {0,0,0,0}
  };

  int opt, opti;
  while ((opt=getopt_long(argc,argv,"t:a:b:p:z:ABPZhf:",long_opts,&opti))!=-1) {
      switch (opt) {
          case 't':
              threshold = atof(optarg);
              break;
          case 'a':
              alpha = atof(optarg);
              break;
          case 'b':
              beta = atof(optarg);
              break;
          case 'p':
              phi0 = atof(optarg);
              break;
          case 'z':
              psi0 = atof(optarg);
              break;
          case 'A':
              alpha_fixed = true;
              break;
          case 'B':
              beta_fixed = true;
              break;
          case 'P':
              phi0_fixed = true;
              break;
          case 'Z':
              psi0_fixed = true;
              break;
          case 'f':
              strncpy(filename, optarg, 127);
              break;
          case 'h':
          default:
              //usage();
              exit(0);
              break;
      }
  }



  vector<int> phase_bin;
  vector<double> phase;
  vector<double> totI, L, Q, U, V;

  Reference::To< Archive > archive = Archive::load( filename );
  if( !archive ) return -1;

  // correct PA to infinite frequency
  FaradayRotation xform;
  xform.set_reference_wavelength( 0 );
  xform.set_measure( archive->get_rotation_measure() );
  xform.execute( archive );


  // Scrunch and Remove baseline
  archive->tscrunch();
  archive->fscrunch();
  if( archive->get_state() != Signal::Stokes) archive->convert_state(Signal::Stokes);
  archive->remove_baseline();

  // Get Data
  Pulsar::Integration* integration = archive->get_Integration(0);
  Pulsar::PolnProfileStats stats;
  stats.set_profile(integration->new_PolnProfile(0));
  Estimate<double> rmsI = sqrt( stats.get_baseline_variance(0) );
  Estimate<double> rmsQ = sqrt( stats.get_baseline_variance(1) );
  Estimate<double> rmsU = sqrt( stats.get_baseline_variance(2) );

  //printf("Rms I = %lf Rms Q = %lf  rmsU= %lf\n", rmsI.get_value(), rmsQ.get_value(), rmsU.get_value());

  double max_L=0.;
  double ph;

  /* 1906 - specific */
  double ex1l = .013867128125, ex1h=.470;
  double ex2l = .542178828125, ex2h=.986787109375;

  double mjd = (double)integration->get_epoch().intday() + integration->get_epoch().fracday();

#if 0
if (mjd < 54700) {ex1l = .0111;}
if ( mjd > 53700) {ex1h = .480;}
if ( mjd > 54300) {ex1h = .490398828125;}
if  ( mjd > 54150 && mjd < 54200) { ex1h = .480468828125;}
//if  ( mjd > 54200 && mjd < 54300) { ex1l = .008298828125; ex1h = .489298828125; ex2l = .533898828125;}
if  ( mjd > 54200 && mjd < 54300) { ex1l = .008298828125; ex1h = .49298828125; ex2l = .533898828125; ex2h = .99; }
if ( mjd > 54400 && mjd < 54480) {ex2l = .527298828125; ex1h = .498099499499; ex2h = .9888; ex1l = 0.0138;}
//if ( mjd > 54400 && mjd < 54480) {ex2l = .526298828125; ex1h = .500198828125; ex2h = .989; ex1l = 0.014;}
if ( mjd > 55000) {ex1h = .494; ex2l=0.549;}
//if ( mjd > 55000) {ex1h = .505;}

#endif


  for (int ibin=0; ibin<archive->get_nbin(); ibin++) {
      ph = ibin/(double) archive->get_nbin();
#if 0
      if ((ex1l <= ph && ph <= ex1h) || (ex2l <= ph && ph <= ex2h)) continue;
#endif
      //if ((ex1l <= ph && ph <= ex1h) || (ex2l <= ph && ph <= ex2h) || (ex3l <= ph && ph <= ex3h)) continue;
      if (integration->get_Profile(0,0)->get_amps()[ibin] > threshold * rmsI.get_value()) {
        phase.push_back((ibin+.5)*(2*M_PI/NB_PTS_PROF));
        phase_bin.push_back(ibin);
        totI.push_back(integration->get_Profile(0,0)->get_amps()[ibin]);
        Q.push_back(integration->get_Profile(1,0)->get_amps()[ibin]);
        U.push_back(integration->get_Profile(2,0)->get_amps()[ibin]);
        V.push_back(integration->get_Profile(3,0)->get_amps()[ibin]);
	L.push_back( sqrt(U.back()*U.back() + Q.back()*Q.back()));
        //cout << totI[ibin] << " " << Q[ibin]<< " " << U[ibin]<<endl;
        cout << ibin << " " <<totI.back() << " " << 0.5 * atan2(U.back(),  Q.back()) <<endl;

	if (L.back() > max_L) max_L = L.back();
      }
  }

	
  RVM_Fcn fFCN(phase, Q, U, rmsQ.get_value(), rmsU.get_value());

  MnUserParameters upar;
  upar.Add("alpha", alpha, 10.);
  upar.Add("beta", beta , 40.);
  upar.Add("phi0", phi0 , 90.);
  upar.Add("psi0", psi0 , 90.);
  nfree = 4;

  if (alpha_fixed) {upar.Fix("alpha"); nfree--;}
  if (beta_fixed) {upar.Fix("beta"); nfree--;}
  if (phi0_fixed) {upar.Fix("phi0"); nfree--;}
  if (psi0_fixed) {upar.Fix("psi0"); nfree--;}

  cout<<"initial parameters: "<<upar<<endl;

  cout<<"start migrad "<<endl;
  MnMinimize migrad(fFCN, upar);
  FunctionMinimum min = migrad();
  /*if(!min.IsValid()) {
    //try with higher strategy
    cout<<"FM is invalid, try with strategy = 2."<<endl;
    MnMigrad migrad2(fFCN, min.UserState(), MnStrategy(2));
    min = migrad2();
  } */
  cout<<"minimum: "<<min<<endl;
  cout<<"Number of points: "<<phase.size()<<endl;
  cout<<"chi**2: " << min.Fval()<<endl;
  cout<<"R chi**2: " << min.Fval()/(double)((int)Q.size() - nfree - 1) <<endl;

  
  cout<<"start Minos"<<endl;
  fFCN.setErrorDef(1.0);
  fFCN.setErrorDef(4.);
  MnMinos Minos(fFCN, min, MnStrategy(2));
  //fFCN.setErrorDef(4.);
  pair<double, double> e0, e1, e2, e3;
  if (!alpha_fixed) e0 = Minos(0);
  if (!beta_fixed) e1 = Minos(1);
  if (!phi0_fixed) e2 = Minos(2);
  if (!psi0_fixed) e3 = Minos(3);

  
  if (!alpha_fixed) cout<<"par0: "<<min.UserState().Value("alpha")<<" + "<<e0.first<<" "<<e0.second<<endl;
  if (!beta_fixed) cout<<"par1: "<<min.UserState().Value("beta")<<" + "<<e1.first<<" "<<e1.second<<endl;
  if (!phi0_fixed) cout<<"par2: "<<min.UserState().Value("phi0")<<" + "<<e2.first<<" "<<e2.second<<endl;
  if (!psi0_fixed) cout<<"par3: "<<min.UserState().Value("psi0")<<" + "<<e3.first<<" "<<e3.second<<endl;

  /*{
    MnScan scan(fFCN, upar);
    cout<<"scan parameters: "<<scan.Parameters()<<endl;
    MnPlot plot;
    for(unsigned int i = 0; i < upar.VariableParameters(); i++) {
      vector<pair<double, double> > xy = scan.Scan(i);
//       vector<pair<double, double> > xy = scan.scan(0);
      plot(xy);
    }
    cout<<scan.Parameters()<<endl;
  }*/

  {
    // create contours factory with FCN and minimum
    MnContours contours(fFCN, min);
    //fFCN.setErrorDef(5.99);
    vector<pair<double,double> > cont4 = contours(0, 1, 20);
    // plot the contours
    MnPlot plot;
    //cont4.insert(cont4.end(), cont.begin(), cont.end());
    plot(min.UserState().Value("alpha"), min.UserState().Value("beta"), cont4);
  }

  cout <<"Simple> "<< mjd<< " " <<min.UserState().Value("alpha")<<" "<<e0.first<<" "<<e0.second <<" "<<min.UserState().Value("beta")<< " "<<e1.first<<" "<<e1.second<<" " <<min.UserState().Value("phi0")<<" "<<e2.first<<" "<<e2.second <<" "<<min.UserState().Value("psi0")<<" "<<e3.first<<" "<<e3.second<< " "<<min.Fval()/(phase.size() - 4 - 1)<<endl;

  printf("Latex> %.1f \& \$%.2f\_\{%.2f\}\^\{+%.2f\}\$ \& \$%.2f\_\{%.2f\}\^\{+%.2f\}\$ \& \$%.2f\_\{%.2f\}\^\{+%.2f\}\$ \& \$%.2f\_\{%.2f\}\^\{+%.2f\}\$ \& %.1f \& %d  \n", mjd,
	min.UserState().Value("alpha"), e0.first, e0.second,
	min.UserState().Value("beta"), e1.first, e1.second,
	min.UserState().Value("phi0"), e2.first, e2.second,
	min.UserState().Value("psi0"), e3.first, e3.second,
	min.Fval(), phase.size());

  cout << "chi**2: " << min.Fval() << endl;
  printf("chi**2: %.14f\n", min.Fval());

  write_results(filename, mjd, archive->get_nbin(), rmsI.get_value(), integration->get_Profile(0,0)->get_amps(), integration->get_Profile(1,0)->get_amps(), integration->get_Profile(2,0)->get_amps(), integration->get_Profile(3,0)->get_amps(), phase_bin, Q, U, min.UserState().Value("alpha"), min.UserState().Value("beta"), min.UserState().Value("phi0"), min.UserState().Value("psi0"));

  return 0;
}
