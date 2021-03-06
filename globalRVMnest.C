#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include "multinest.h"
#include "Parameters.h"
#include "RVMnest.h"
#include "read_results.h"

#include <mpi.h>

// psrchive stuff
#include "Pulsar/Archive.h"
#include "Pulsar/Integration.h"
#include "Pulsar/Profile.h"
#include "Pulsar/FaradayRotation.h"
#include "Pulsar/PolnProfileStats.h"

using namespace std;
using namespace Pulsar;
#ifdef HAVE_POLYCHORD
#include "interfaces.hpp"
#include "globalRVM_likelihood_PC.h"
#endif
MNStruct *sp;

/************************************************* dumper routine ******************************************************/

// The dumper routine will be called every updInt*10 iterations
// MultiNest doesn not need to the user to do anything. User can use the arguments in whichever way he/she wants
//
//
// Arguments:
//
// nSamples 						= total number of samples in posterior distribution
// nlive 						= total number of live points
// nPar 						= total number of parameters (free + derived)
// physLive[1][nlive * (nPar + 1)] 			= 2D array containing the last set of live points (physical parameters plus derived parameters) along with their loglikelihood values
// posterior[1][nSamples * (nPar + 2)] 			= posterior distribution containing nSamples points. Each sample has nPar parameters (physical + derived) along with the their loglike value & posterior probability
// paramConstr[1][4*nPar]:
// paramConstr[0][0] to paramConstr[0][nPar - 1] 	= mean values of the parameters
// paramConstr[0][nPar] to paramConstr[0][2*nPar - 1] 	= standard deviation of the parameters
// paramConstr[0][nPar*2] to paramConstr[0][3*nPar - 1] = best-fit (maxlike) parameters
// paramConstr[0][nPar*4] to paramConstr[0][4*nPar - 1] = MAP (maximum-a-posteriori) parameters
// maxLogLike						= maximum loglikelihood value
// logZ							= log evidence value from the default (non-INS) mode
// INSlogZ						= log evidence value from the INS mode
// logZerr						= error on log evidence value
// context						void pointer, any additional information

void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context)
{
	// convert the 2D Fortran arrays to C++ arrays
	
	
	// the posterior distribution
	// postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
	
	int i, j;
	
	double postdist[nSamples][nPar + 2];
	for( i = 0; i < nPar + 2; i++ )
		for( j = 0; j < nSamples; j++ )
			postdist[j][i] = posterior[0][i * nSamples + j];
	
	
	
	// last set of live points
	// pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
	
	double pLivePts[nlive][nPar + 1];
	for( i = 0; i < nPar + 1; i++ )
		for( j = 0; j < nlive; j++ )
			pLivePts[j][i] = physLive[0][i * nlive + j];
}

/***********************************************************************************************************************/




/************************************************** Main program *******************************************************/



int main(int argc, char *argv[])
{
	// set the MultiNest sampling parameters
	
	int IS = 0;	 // IS =1 for evidence comparison				// do Nested Importance Sampling?
	int mmodal = 0;					// do mode separation?
	int ceff = 1;					// run in constant efficiency mode?
	int nlive = 2000;				// number of live points
	double efr = 0.05;				// set the required efficiency
	double tol = 0.5;				// tol, defines the stopping criteria
	int ndims = 4;					// dimensionality (no. of free parameters)
	int nPar = 4;					// total no. of parameters including free & derived parameters
	int nClsPar = 2;				// no. of parameters to do mode separation on
	int updInt = 5000;				// after how many iterations feedback is required & the output files should be updated
							// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored
	int maxModes = 10;				// expected max no. of modes (used only for memory allocation)
	
	char filename[128];
	char root[160] = "chains/globalRVMnest-";		// root for output files
	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock
	int fb = 1;					// need feedback on standard output?
	int resume = 1;					// resume from a previous job?
	int outfile = 1;				// write output files?
	int initMPI = 0;				// initialize MPI routines?, relevant only if compiling with MPI
							// set it to F if you want your main program to handle MPI initialization
	double logZero = -1E90;				// points with loglike < logZero will be ignored by MultiNest
	int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
							// has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
	void *context = 0;				// not required by MultiNest, any additional information user wants to pass

	char confname[160];				// root for output files
	int inc_fixed=1;
	int prate_fixed=1;
	int have_efac=1;
	int have_aberr_offset=0;
	int nfiles_aberr = 0;
	double threshold=1.8;
	int margin_phi0=0;
	int nfiles = argc - 1;
	int psi_jump_fixed=1;
	int sampler = 0;
	int sin_psi = 0;
	int ascii_output =1;
	int totnpts = 0;
	int rank, size;
	MPI_Comm world_comm;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_dup(MPI_COMM_WORLD, &world_comm);
	
	// Read Params from config files
	param p;
	p.numfiles = nfiles;
	if (rank==0)  cout << "Reading parameters" << endl; 
	int rv = readParameters(&p, "config.txt");

	if (rv == EXIT_SUCCESS) {
	  sampler = p.sampler;
	  IS = p.IS;
	  nlive = p.nlive;
	  ceff = p.ceff;
	  efr = p.efr;
	  strcpy(root, p.basename);
	  
	  inc_fixed = p.inc_fixed;
	  prate_fixed = p.prate_fixed;
	  psi_jump_fixed = p.psi_jump_fixed;
	  have_efac = p.have_efac;
	  have_aberr_offset = p.have_aberr_offset;
	  threshold = p.threshold;
	  margin_phi0 = p.margin_phi0;
	  sin_psi = p.sin_psi;
	  
	    // Copy config file
	  sprintf(confname,"%s.config", root);
	  ifstream  src("config.txt", ios::binary);
	  ofstream  dst(confname,   ios::binary);
	  dst << src.rdbuf();
	} else {
	  for (int i=0; i<nfiles; i++)  p.n_phs_exclude[i] = 0;
	}


	vector <int> nbin;
	vector< vector<double> > phase, I, Q, U, L, V;
	vector <double> RMS_I, RMS_Q, RMS_U, MJD;
	phase.resize(nfiles);	
	I.resize(nfiles);	
	Q.resize(nfiles);	
	U.resize(nfiles);	
	L.resize(nfiles);
	V.resize(nfiles);

	// Process and read files
	//stringstream s;
	//string fn="file";
	for (unsigned i=0; i<nfiles; i++) {

	  Reference::To< Archive > archive = Archive::load( argv[i+1] );
	  if (rank==0) cout << "Opening " << argv[i+1] << endl;
	  if( !archive ) return -1;

	  nbin.push_back(archive->get_nbin());
	  //cout << "  nbin : " << nbin.back()<<endl;
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
	  
	  double mjd = (double)integration->get_epoch().intday() + integration->get_epoch().fracday();

	  double max_L=0.;
	  double ph=0.;
	  bool skip_bin = false;

	  // TODO: Specific to 1906 analysis
	  if (mjd < 57300) nfiles_aberr++;

	  if (rank==0) {
	    for (int nphs=0 ; nphs < p.n_phs_exclude[i]; nphs++)
	      cout << "File #" << i << " : will exclude phase " << p.phs_exclude[i][nphs*2] << " to " <<  p.phs_exclude[i][1 + nphs*2] << endl;
	  }

	  // 
	  stringstream str;
	  ofstream myf;
	  str << integration->get_epoch().intday();
	  if (ascii_output) {
	    string fnf = str.str();
	    myf.open(fnf.c_str());
	    myf << scientific;
	  }

	  for (int ibin=0; ibin<archive->get_nbin(); ibin++) {
	    if (ascii_output) {
	      myf << ibin << " ";
	      for (int ik=0; ik < 4; ik++)
		myf << integration->get_Profile(ik,0)->get_amps()[ibin]<< " "; 
	      myf << endl;
	    }

            ph = ibin/(double) archive->get_nbin();
	
	    // Add all points to I vector
	    I[i].push_back(integration->get_Profile(0,0)->get_amps()[ibin]);
	    L[i].push_back( sqrt( pow(integration->get_Profile(1,0)->get_amps()[ibin], 2) + pow(integration->get_Profile(2,0)->get_amps()[ibin], 2) ));
	    V[i].push_back(integration->get_Profile(3,0)->get_amps()[ibin]);

            // Exclude phase range
            skip_bin = false;
	    for (int nphs=0; nphs < p.n_phs_exclude[i]; nphs++)
	      if (p.phs_exclude[i][nphs*2] <= ph && ph <= p.phs_exclude[i][1 + nphs*2]) skip_bin = true;
	    
	    if (skip_bin) continue;
            
	    if (integration->get_Profile(0,0)->get_amps()[ibin] > threshold * rmsI.get_value()) {
	      phase[i].push_back((ibin+.5)*(2*M_PI/(double) archive->get_nbin()));
	      Q[i].push_back(integration->get_Profile(1,0)->get_amps()[ibin]);
	      U[i].push_back(integration->get_Profile(2,0)->get_amps()[ibin]);

	    }
	  }
	  if (ascii_output) myf.close();

	  if (rank==0)  cout << "Number of data points " << Q[i].size() << endl;
	  totnpts += Q[i].size();
	  MJD.push_back(integration->get_epoch().in_days());
	  RMS_I.push_back(rmsI.get_value());
	  RMS_Q.push_back(rmsQ.get_value());
	  RMS_U.push_back(rmsU.get_value());
	}
	
	// Init struct
	context = init_struct(nfiles, phase , I, Q, U, L, V, RMS_I, RMS_Q, RMS_U, nbin, p.njump);

	MNStruct *par = ((MNStruct *)context);

	par->sampler = sampler;
	par->inc_fixed = inc_fixed;
	par->prate_fixed = prate_fixed;
	par->psi_jump_fixed = psi_jump_fixed;
	par->inc = 43.7 * M_PI / 180.;
	par->prate = 2.234;
	par->have_efac = have_efac;
	par->have_aberr_offset = have_aberr_offset;
	par->margin_phi0 = margin_phi0;
	par->sin_psi = sin_psi;
	par->nfiles_aberr = nfiles_aberr;
	for(unsigned i = 0; i < nfiles; i++) par->epoch[i] = MJD[i];
		
	ndims = nPar = 4;
	if (!par->margin_phi0) {ndims+=nfiles; nPar+=nfiles;}  // nfiles for phi0
	//ndims+=nfiles; nPar+=nfiles;  // nfiles for phi0
	if (!par->prate_fixed) {ndims+=1; nPar+=1;}
	if (!par->inc_fixed) {ndims+=1; nPar+=1;}
	if (par->have_efac) {ndims+=1; nPar+=1;}
	if (par->have_aberr_offset) {ndims+=nfiles_aberr; nPar+=nfiles_aberr;}
	if (par->sin_psi) {ndims+=1; nPar+=1;}

	par->njump = 0;
	// Copy the range of parameters
	if (rv == EXIT_SUCCESS) {
	  par->inc = p.inc * M_PI / 180.;
	  par->prate = p.prate;
	  par->r_alpha = p.alpha;
	  par->r_delta = p.delta;
	  par->r_Phi0 = p.Phi0;
	  par->r_phi0 = p.phi0;
	  par->r_inc = p.r_inc;
	  par->r_prate = p.r_prate;
	  par->r_efac = p.efac;
	  
	  if (p.njump) {
	      par->njump = p.njump;
	      par->r_psi_jump = p.r_psi_jump;
	      par->psi_jumps = p.psi_jumps;
	      par->psi_jump_MJD = p.psi_jump_MJD;
	      if (!par->psi_jump_fixed) {
		  ndims += p.njump;
		  nPar += p.njump;
	      }
	      for(int i = 0; i < par->njump; i++) par->psi_jumps[i] *= M_PI/180;
	  }
	}
	par->do_plot = 0;

	int pWrap[ndims];
	for(int i = 0; i < ndims; i++) pWrap[i] = 0;
	if (par->r_delta[0]==0 && par->r_delta[0]==180.) pWrap[1] = 1;
	if (par->r_Phi0[0]==0 && par->r_Phi0[0]==360.) pWrap[2] = 1;
        if (!par->sin_psi) pWrap[3] = 1;
	
	if (rank==0) {
	  // calling MultiNest
	  cout << endl << " -- Parameters -- " << endl;
	  cout << "Inc fixed   = " << inc_fixed << endl;
	  if (prate_fixed) cout << "Inclination angle (deg) = " << par->inc * 180./M_PI<< endl;
	  cout << "Prate fixed = " << prate_fixed << endl;
	  if (prate_fixed) cout << "Prate (deg/yr) = " << par->prate << endl;
	  cout << "Have EFAC   = " << have_efac << endl;
	  cout << "Threshold   = " << threshold << endl;
	  cout << "Margin_phi0 = " << margin_phi0 << endl;
	  cout << "nlive = " << nlive << endl;
	  cout << "ndims = " << ndims << endl;
	  cout << "sin_psi = " << sin_psi << endl;
	  cout << "Will model "<< p.njump << " Psi0 jumps "<< endl;
	  for(int i = 0; i < par->njump; i++) cout << "Introduced a Psi0 offset between MJDs "<< par->psi_jump_MJD[2*i] << " "<< par->psi_jump_MJD[2*i+1]<< endl;
	  cout << "Psi_jump fixed = " << psi_jump_fixed << endl;
	  cout << "Assuming reading " << nfiles << " files" << endl;
	  cout << "Basefilename " << root << endl;
	  cout << "Total data points " << totnpts << endl;
	  cout << endl;
	}

	if (sampler==0) {
	  sprintf(filename,"%s/chainsMN-", root);
	  nested::run(IS, mmodal, ceff, p.nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, filename, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, globalRVMLogLike, dumper, context);
	}
	else if (sampler==1) {
#ifdef HAVE_POLYCHORD
	  Settings settings;
          settings.nDims         = ndims;
          settings.nDerived      = 1; // Nderived is the sign of Ltot. Will be used to correct Psi0, e.g if Ltot<0, psi0+=Pi/2
          settings.nlive         = 500;
          settings.num_repeats   = settings.nDims*4;
          settings.do_clustering = false;
          settings.precision_criterion = 1e-3;
          settings.base_dir.assign(root);
          settings.file_root     = "chainsPC";
          settings.write_resume  = true;
          settings.read_resume   = true;
          settings.write_live    = true;
          settings.write_dead    = false;
          settings.write_stats   = true;
          settings.equals        = true;
          settings.posteriors    = true;
          settings.cluster_posteriors = false;
          settings.feedback      = 2;
          //settings.update_files  = settings.nlive;
          settings.boost_posterior= 5.0;

          //setup_loglikelihood();                                                                             
          sp = par;
          run_polychord(globalRVMLogLike_PC, prior, settings, world_comm) ;
#else
	  cerr << "PolyChord library not detected during configure. Aborting! "<< endl;
	  return(-1);
#endif
	}
	

	if (rank == 0) {
	  // Read results from stats file
	  read_stats(root, nPar, par);
	  get_globalRVM_chi2(par);
	}

	//MPI_Finalize();

	return(0);
}

/***********************************************************************************************************************/
