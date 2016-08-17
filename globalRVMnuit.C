/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
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

#include "globalRVMKJ_Fcn.h"

#include "getopt.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"

#define NB_PTS_PROF 2048.0
#define R2D (180.0/M_PI)

using namespace std;
using namespace ROOT::Minuit2;
using namespace Pulsar;

int main(int argc, char **argv) {

  int ndat = 0,nfree = 0;
  double threshold=3;
  double alpha=100.0, delta=-6.0, phase0=18.0, psi00=30.0, inc=45.2, prate= 2.26;
  bool alpha_fixed=false, delta_fixed=false, phase0_fixed=false, inc_fixed=false, psi00_fixed=false, phi0_fixed=false;
  bool prate_fixed=false;
  char filename[128]="test.pa";
  char label[16];
  double phi0[20];
phi0[0]=18.9077;                                                                 
phi0[1]=19.1344;                                                                 
phi0[2]=19.4348;                                                                 
phi0[3]=20.0993;                                                                 
phi0[4]=20.0993;                                                                 
phi0[5]=20.2774;                                                                 
phi0[6]=20.6447;                                                                 
phi0[7]=21.47;                                                                   
phi0[8]=20.9732;                                                                 
phi0[9]=21.8713;                                                                 
phi0[10]=23.6219;                                                                
phi0[11]=23.1055;                                                                
phi0[12]=24.1256;


  static struct option long_opts[] = {
      {"thresh",  1, NULL, 't'},
      {"alpha",  1, NULL, 'a'},
      {"delta",  1, NULL, 'd'},
      {"phase0", 1, NULL, 'p'},
      {"inc",    1, NULL, 'i'},
      {"psi00",  1, NULL, 'z'},
      {"rate",   1, NULL, 'r'},
      {"Alpha",  0, NULL, 'A'},
      {"Delta",  0, NULL, 'D'},
      {"Phase0", 0, NULL, 'P'},
      {"Inc",    0, NULL, 'I'},
      {"Psi00",  0, NULL, 'Z'},
      {"Phi0",   0, NULL, 'B'},
      {"Rate",   0, NULL, 'B'},
      {"help",   0, NULL, 'h'},
      {"file",   0, NULL, 'f'},
      {0,0,0,0}
  };

  int opt, opti;
  while ((opt=getopt_long(argc,argv,"t:a:d:p:i:z:r:ADPIZBRhf:",long_opts,&opti))!=-1) {
      switch (opt) {
          case 't':
              threshold = atof(optarg);
              break;
          case 'a':
              alpha = atof(optarg);
              break;
          case 'd':
              delta = atof(optarg);
              break;
          case 'p':
              phase0 = atof(optarg);
	      break;
	  case 'i':
	      inc = atof(optarg);
              break;
          case 'z':
              psi00 = atof(optarg);
              break;
          case 'r':
              prate = atof(optarg);
              break;
          case 'A':
              alpha_fixed = true;
              break;
          case 'D':
              delta_fixed = true;
              break;
          case 'P':
              phase0_fixed = true;
              break;
          case 'I':
              inc_fixed = true;
              break;
          case 'Z':
              psi00_fixed = true;
              break;
          case 'B':
              phi0_fixed = true;
              break;
          case 'R':
              prate_fixed = true;
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

  // -- PSR Parameters --
  MNStruct *params;
  params = (MNStruct *) malloc(sizeof(MNStruct));

  vector< vector<double> > phase;
  vector< vector<double> > totI, L, Q, U;
  vector<double> varI, varQ, varU;

  // Read list of filenames
  int nfiles = 0;
  char file[128], line[128];
  FILE *pfi;
  vector< string > filenames;
  if( (pfi=fopen(filename,"r"))==NULL) {printf("Can't open %s\n",filename);exit(0);}
  while (fgets(line,128,pfi) &&  !feof(pfi)) {
      if (strncmp(&line[0],"#", 1)==0) continue; // Skip comments
      sscanf(line, "%s", file);
      filenames.push_back(file);
      nfiles++;
  }    
  // Resize vectors
  phase.resize(nfiles);
  totI.resize(nfiles);
  L.resize(nfiles);
  Q.resize(nfiles);
  U.resize(nfiles);


  // Now loop over all files
  for (unsigned ifile=0; ifile<filenames.size(); ifile++) {
      Reference::To< Archive > archive = Archive::load( filenames[ifile] );
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
      varI.push_back( stats.get_baseline_variance(0).get_value() );
      varQ.push_back( stats.get_baseline_variance(1).get_value() );
      varU.push_back( stats.get_baseline_variance(2).get_value() );

      double ph=0.;
      double ex1l = 0.065, ex1h=0.495;
      double ex2l = 0.62, ex2h=1.;
      double ex3l = 0.0, ex3h=0.03;

      for (unsigned ibin=0; ibin<archive->get_nbin(); ibin++) {
	  ph = ibin/(double) archive->get_nbin();

	  // Exclude phase range
	  if ((ex1l <= ph && ph <= ex1h) || (ex2l <= ph && ph <= ex2h) || (ex3l <= ph && ph <= ex3h)) continue;
	  if (integration->get_Profile(0,0)->get_amps()[ibin] > threshold * sqrt(varI[ifile])) {
	    phase[ifile].push_back((ibin+.5)*(2*M_PI/NB_PTS_PROF));
	    totI[ifile].push_back(integration->get_Profile(0,0)->get_amps()[ibin]);
	    Q[ifile].push_back(integration->get_Profile(1,0)->get_amps()[ibin]);
	    U[ifile].push_back(integration->get_Profile(2,0)->get_amps()[ibin]);
	    L[ifile].push_back( sqrt(U[ifile].back()*U[ifile].back() + Q[ifile].back()*Q[ifile].back()));
	    ndat++;
	  }
      }
      params->epoch[ifile] = (double)integration->get_epoch().intday() + integration->get_epoch().fracday();
      cout << params->epoch[ifile] << " "<< ndat<< " "<< Q[ifile].size() << endl;
  }    

  params->n_epoch = nfiles;

  RVM_Fcn fFCN(phase, Q, U, varQ, varU);
  fFCN.set_params(params);

  MnUserParameters upar;
  upar.Add("alpha", alpha, 10.); nfree++;
  upar.Add("delta", delta , 40.); nfree++;
  upar.Add("phase0", phase0 , 90.); nfree++;
  upar.Add("inc",   inc , 1.); nfree++;
  upar.Add("psi00", psi00 , 90.); nfree++;
  upar.Add("prec_rate", prate , 5.); nfree++;
  for(unsigned int i = 0; i < Q.size(); i++) {
      sprintf(label, "phi0_%d", i);
      upar.Add(label, phi0[i], 20.);
      nfree++;
  }
#if 0   
  for (unsigned i=0; i< Q.size(); i++) {
    for (unsigned j=0; j< Q[i].size(); j++) {
      sprintf(label, "L%d_%d", i, j);
      upar.Add(label, L[i][j], 10.0);
      //printf("i=%d j=%d\n", i, j);fflush(stdout);
    }  
  }    
#endif

  if (alpha_fixed) {upar.Fix("alpha"); nfree--;}
  if (delta_fixed) {upar.Fix("delta"); nfree--;}
  if (inc_fixed)   {upar.Fix("inc"); nfree--;}
  if (phase0_fixed) {upar.Fix("phase0"); nfree--;}
  if (psi00_fixed)  {upar.Fix("psi00"); nfree--;}
  if (prate_fixed)  {upar.Fix("prec_rate"); nfree--;}
  if (phi0_fixed) {
      for(unsigned int i = 0; i < Q.size(); i++) {
          sprintf(label, "phi0_%d", i);
          upar.Fix(label); nfree--;
      }	  
  }    

  cout<<"initial parameters: "<<upar<<endl;
  fflush(stdout);

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
  cout<<"Number of free parameters: "<< nfree <<endl;
  cout<<"R chi**2:" << min.Fval()/(2*ndat - nfree - 1) <<endl;

  
  cout<<"start Minos"<<endl;
  fFCN.setErrorDef(1.0);
  MnMinos Minos(fFCN, min, MnStrategy(2));
  //fFCN.setErrorDef(4.);
  pair<double, double> e0, e1, e2, e3, e4, e5;
  if (!alpha_fixed) e0 = Minos(0);
  if (!delta_fixed) e1 = Minos(1);
  if (!phase0_fixed) e2 = Minos(2);
  if (!inc_fixed) e3 = Minos(3);
  if (!psi00_fixed) e4 = Minos(4);
  if (!prate_fixed) e5 = Minos(5);

  if (!alpha_fixed)  cout<<"    alpha: "<<min.UserState().Value("alpha")<<" +/- "<<e0.first<<" "<<e0.second<<endl;
  if (!delta_fixed)  cout<<"    delta: "<<min.UserState().Value("delta")<<" +/- "<<e1.first<<" "<<e1.second<<endl;
  if (!phase0_fixed) cout<<"   phase0: "<<min.UserState().Value("phase0")<<" +/- "<<e2.first<<" "<<e2.second<<endl;
  if (!inc_fixed)    cout<<"      inc: "<<min.UserState().Value("inc")<<" +/- "<<e3.first<<" "<<e3.second<<endl;
  if (!psi00_fixed)  cout<<"    psi00: "<<min.UserState().Value("psi00")<<" +/- "<<e4.first<<" "<<e4.second<<endl;
  if (!prate_fixed)  cout<<"prec_rate: "<<min.UserState().Value("prec_rate")<<" +/- "<<e5.first<<" "<<e5.second<<endl;

  
  fFCN.setErrorDef(4.);
  if (!alpha_fixed) e0 = Minos(0);
  if (!delta_fixed) e1 = Minos(1);
  if (!phase0_fixed) e2 = Minos(2);
  if (!inc_fixed) e3 = Minos(3);
  if (!psi00_fixed) e4 = Minos(4);
  if (!prate_fixed) e5 = Minos(5);

  if (!alpha_fixed)  cout<<"    alpha: "<<min.UserState().Value("alpha")<<" +/- "<<e0.first<<" "<<e0.second<<endl;
  if (!delta_fixed)  cout<<"    delta: "<<min.UserState().Value("delta")<<" +/- "<<e1.first<<" "<<e1.second<<endl;
  if (!phase0_fixed) cout<<"   phase0: "<<min.UserState().Value("phase0")<<" +/- "<<e2.first<<" "<<e2.second<<endl;
  if (!inc_fixed)    cout<<"      inc: "<<min.UserState().Value("inc")<<" +/- "<<e3.first<<" "<<e3.second<<endl;
  if (!psi00_fixed)  cout<<"    psi00: "<<min.UserState().Value("psi00")<<" +/- "<<e4.first<<" "<<e4.second<<endl;
  if (!prate_fixed)  cout<<"prec_rate: "<<min.UserState().Value("prec_rate")<<" +/- "<<e5.first<<" "<<e5.second<<endl;

  /*{
    // create contours factory with FCN and minimum
    MnContours contours(fFCN, min);
    //fFCN.setErrorDef(5.99);
    vector<pair<double,double> > cont4 = contours(0, 1, 20);
    // plot the contours
    MnPlot plot;
    //cont4.insert(cont4.end(), cont.begin(), cont.end());
    plot(min.UserState().Value("alpha"), min.UserState().Value("beta"), cont4);
  }*/

 // cerr <<integration->get_epoch() << " "<< min.UserState().Value("alpha")<<" "<<e0.first<<" "<<e0.second <<" "<<min.UserState().Value("beta")<< " "<<e1.first<<" "<<e1.second<<" " <<min.UserState().Value("phi0")<<" "<<e2.first<<" "<<e2.second <<" "<<min.UserState().Value("psi0")<<" "<<e3.first<<" "<<e3.second<< " "<<min.Fval()/(phase.size() - 4 - 1)<<endl;

  return 0;
}
