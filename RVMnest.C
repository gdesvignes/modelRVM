#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "multinest.h"
#include "Parameters.h"
#include "RVMnest.h"
#include "read_results.h"

// psrchive stuff
#include "Pulsar/Archive.h"
#include "Pulsar/Integration.h"
#include "Pulsar/Profile.h"
#include "Pulsar/FaradayRotation.h"
#include "Pulsar/PolnProfileStats.h"

using namespace std;
using namespace Pulsar;
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
	int IS = 0;					// do Nested Importance Sampling?
	int mmodal = 0;					// do mode separation?
	int ceff = 1;					// run in constant efficiency mode?
	int nlive = 2000;				// number of live points
	double efr = 0.05;				// set the required efficiency
	double tol = 0.5;				// tol, defines the stopping criteria
	int ndims = 4;					// dimensionality (no. of free parameters)
	int nPar = 4;					// total no. of parameters including free & derived parameters
	int nClsPar = 2;				// no. of parameters to do mode separation on
	int updInt = 1000;				// after how many iterations feedback is required & the output files should be updated
							// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored
	int maxModes = 100;				// expected max no. of modes (used only for memory allocation)
	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(int i = 0; i < ndims; i++) pWrap[i] = 1;
	char filename[128];
	char root[100];
	char tmproot[100] = "chains";		        // root for output files
	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock
	int fb = 1;					// need feedback on standard output?
	int resume = 1;					// resume from a previous job?
	int outfile = 1;				// write output files?
	int initMPI = 1;				// initialize MPI routines?, relevant only if compiling with MPI
							// set it to F if you want your main program to handle MPI initialization
	double logZero = -1E90;				// points with loglike < logZero will be ignored by MultiNest
	int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
							// has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
	void *context = 0;				// not required by MultiNest, any additional information user wants to pass
	double threshold=1.0;
	int have_efac=0;
	int nfiles = 1;
	vector< vector<double> > phase, I, Q, U, L, V;
	vector <int> nbin;
	vector <double> RMS_I, RMS_Q, RMS_U;

	// Read Params from config files
	param p;
	p.numfiles = 1;
        cout << "Reading parameters" << endl;
	int rv = readsimpleParameters(&p, "config.txt");
	if (rv == EXIT_SUCCESS) {
	  IS = p.IS;
	  nlive = p.nlive;
	  ceff = p.ceff;
	  efr = p.efr;
	  strcpy(tmproot, p.basename);
	  threshold = p.threshold;
	  have_efac = p.have_efac;
	}

	strcpy(filename, argv[1]);
	phase.resize(nfiles);
        I.resize(nfiles);
        Q.resize(nfiles);
        U.resize(nfiles);
        L.resize(nfiles);
	V.resize(nfiles);
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

	bool skip_bin = false;
	double max_L=0.;
	double ph;
	double ex1l, ex1h, ex2l, ex2h;
	nbin.push_back(archive->get_nbin());

	mkdir(tmproot, 0700);
	sprintf(root, "%s/RVM-", tmproot);

	for (int ibin=0; ibin<archive->get_nbin(); ibin++) {
	    ph = ibin/(double) archive->get_nbin();

	    // Exclude phase range
            skip_bin = false;
	    for (int nphs=0; nphs < p.n_phs_exclude[0]; nphs++)
              if (p.phs_exclude[0][nphs*2] <= ph && ph <= p.phs_exclude[0][1 + nphs*2]) skip_bin = true;
            if (skip_bin) continue;

	    if (integration->get_Profile(0,0)->get_amps()[ibin] > threshold * rmsI.get_value()) {
		phase[0].push_back((ibin+.5)*(2*M_PI/(double) archive->get_nbin()));
		I[0].push_back(integration->get_Profile(0,0)->get_amps()[ibin]);
		Q[0].push_back(integration->get_Profile(1,0)->get_amps()[ibin]);
		U[0].push_back(integration->get_Profile(2,0)->get_amps()[ibin]);
		L[0].push_back( sqrt(U[0].back()*U[0].back() + Q[0].back()*Q[0].back()));
		V[0].push_back(integration->get_Profile(3,0)->get_amps()[ibin]);
		cout << ibin << " " <<I[0].back() << endl;
	    }
	}

	RMS_I.push_back(rmsI.get_value());
	RMS_Q.push_back(rmsQ.get_value());
	RMS_U.push_back(rmsU.get_value());

	// Init struct
	context = init_struct(nfiles, phase , I, Q, U, L, V, RMS_I, RMS_Q, RMS_U, nbin, 0);
	MNStruct *par = ((MNStruct *)context);
	par->do_plot = 0;
	par->have_efac = have_efac;
	par->epoch[0] = (double)integration->get_epoch().intday() + integration->get_epoch().fracday();

	if (par->have_efac) {ndims+=1; nPar+=1;}
	// Copy the range of parameters
        if (rv == EXIT_SUCCESS) {
	  par->r_alpha = p.alpha;
	  par->r_beta = p.beta;
	  par->r_phi0 = p.phi0;
	  par->r_psi0 = p.psi0;
	  par->r_efac = p.efac;
	} else {
	  par->r_alpha = NULL;
	  par->r_beta = NULL;
	  par->r_phi0 = NULL;
	  par->r_psi0 = NULL;
	}
	
	// calling MultiNest
	nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI,
	logZero, maxiter, RVMLogLike, dumper, context);


	read_statsRVM(root, nPar, par);

	get_RVM_chi2(par);
}

/***********************************************************************************************************************/
