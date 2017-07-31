#include <iostream>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <limits>
#include "RVMnest.h"

using namespace std;

int read_stats(char *root, int npar, MNStruct *p)
{

  int ipar=0,lpar=npar+1;
  string line;
  char filename[256];
  double likelihood=numeric_limits<double>::min();
  double val, val_err;

  cout << "Reading results " << endl;

  if (p->sampler == 0) {
    sprintf(filename, "%s/chainsMN-phys_live.points", root);
    lpar--;
  }
  else
  sprintf(filename, "%s/chainsPC_phys_live.txt", root);

  ifstream stats(filename);

  // Ask to output the result
  p->do_plot = 1;

  if (!stats.is_open())
    {
      cerr << "Error opening file "<<filename << endl;
      return(-1);
    }

  while (getline(stats, line)) 
    {
      istringstream iss(line);
      vector<double> cols{istream_iterator<double>{iss},
	  istream_iterator<double>{}};

      ipar = 0;
      if (cols[lpar] > likelihood) {
	likelihood = cols[lpar];

	p->alpha = cols[ipar] * M_PI/180.; ipar++;
	p->delta = cols[ipar] * M_PI/180.; ipar++;
        p->phase0 = cols[ipar] * M_PI/180.; ipar++;
        p->psi00 = cols[ipar] * M_PI/180.; ipar++;

	// Set RVM phase
	if(!p->margin_phi0)
	  {
	    for (unsigned int j = 0; j < p->n_epoch; j++)
	      {
		p->phi0[j] = cols[ipar] * M_PI/180.; ipar++;
	      }
	  }

	if (p->have_aberr_offset)
	  {
	    for (unsigned int j = 0; j < p->n_epoch; j++)
	      {
		p->phi_aberr_offset[j] = cols[ipar] * M_PI/180.; ipar++;
	      }
	  }
	
	// Set precession rate                                                                                                                               
	if (!p->prate_fixed)
	  {
	    p->prate = cols[ipar]; ipar++;
	  }
	
	// Set inclination angle                                                                                                                             
	if (!p->inc_fixed)
	  {
	    p->inc = cols[ipar] * M_PI/180.; ipar++;
	  }
	
	// Phase Jump
	if (p->njump && !p->psi_jump_fixed)
	  {
	    for (unsigned  j=0; j<p->njump; j++)
	      {
		p->psi_jumps[j] = cols[ipar] * M_PI/180.; ipar++;
	      }
	  }

      }
    }
  
  cout << "Finished reading stats file "<< filename << endl;
  cout << "Alpha = " << p->alpha * 180./M_PI << endl
       << "Delta = " << p->delta * 180./M_PI << endl
       << "Phase = " << p->phase0 * 180./M_PI << endl
       << "Psi00 = " << p->psi00 * 180./M_PI << endl
       << "Inc   = " << p->inc * 180./M_PI << endl
       << "Prate = " << p->prate << endl;
  for (unsigned int j = 0; j < p->n_epoch; j++)
    {
      cout << "Phi  = " << p->phi0[j] * 180./M_PI << endl;
    }    
  cout << "Likelihood = " << likelihood << endl;
  
  return(0);
}
