#include "string.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include "RVMnest.h"

using namespace std;

void read_stats(char *root, int npar, MNStruct *p)
{

  int i, ipar,ioff, nmodes;
  string line, a,b,c;
  char filename[256];

  double val, val_err;

  sprintf(filename, "%sstats.dat", root);

  ifstream stats(filename);

  // Ask to output the result
  p->do_plot = 1;

  if (stats.is_open())
    {
      //getline(stats, line);

      for(i=0; i<6; i++) getline(stats, line);
      // Read the number of modes
      istringstream iss(line);
      iss >> a; iss>>b; iss >> c; iss>> nmodes;
      cout << "Found Nmodes: " << nmodes << endl;
      
      for(i=0; i<7; i++) getline(stats, line);

      for(i=0; i<npar; i++) getline(stats, line);
      
      // Skip 3 lines to go to ML results
      for(i=0; i<3; i++) getline(stats, line);

      for(i=0; i<npar; i++) 
	{
	  getline(stats, line);
	  istringstream iss(line);
	  iss >> ipar; iss >> val;
	  cout << ipar << " " << val << endl;
	  
	  if (i==0) p->alpha = val * M_PI/180.;
	  if (i==1) p->delta = val * M_PI/180.;
	  if (i==2) p->phase0 = val * M_PI/180.;
	  if (i==3) p->psi00 = val * M_PI/180.;
	  
	  // Set RVM phase
	  if(!p->margin_phi0) 
	    {
	      for (unsigned int j = 0; j < p->n_epoch; j++)
		{
		  if (i==4+j) p->phi0[j] = val * M_PI/180.; 
		}
	    }
	  
	  // Set precession rate
	  if (!p->prate_fixed) {
	    if (i==4 + p->n_epoch)
	      p->prate = val;
	  }
	  
	  // EFACS are not supported yet !!!
	  
	  // Set inclination angle
	  if (!p->inc_fixed) 
	    {
	      if (!p->prate_fixed) ioff=1;
	      else ioff=0;
	      if (i==4 + p->n_epoch + ioff) p->inc = val * M_PI/180.;
	    }
	  
	  if (p->njump) 
	    {
	      if (!p->prate_fixed) ioff=1;
	      else ioff=0;
	      if (!p->inc_fixed) ioff+=1;
	      for (unsigned  j=0; j<p->njump; j++)
		{
		  if (i == 4+p->n_epoch + ioff + j) 
		    {
		      cout <<"Test " << 4+p->n_epoch + ioff + j << " "<< val << endl;
		      p->psi_jumps[j] = val * M_PI/180.;
		    }
		}
	    }
	}
    }
  else
    {
      cout << "Unable to open file "<< filename << endl;
    }
  cout << "Finished reading stats file "<< filename << endl;
}

