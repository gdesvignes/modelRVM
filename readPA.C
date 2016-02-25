#include "string.h"
#include "stdio.h"
#include "stdlib.h"

#include <iostream>
#include <vector>
#include <math.h>

#include "RVMnest.h"


using namespace std;

void readPA(char filename[127], vector<double> &  x, vector<double> & Q, vector<double> & U)
{
    FILE *pfi;
    char line[256];
    if( (pfi=fopen(filename,"r"))==NULL) {printf("Can't open %s\n",filename);exit(0);}

    int a, b, x_tmp;
    float polI, polQ, polU, polV;
    double y_tmp, e_tmp;
    fgets(line,256,pfi); // -- Skip first line --
    cout << "Reading file from psrchive" << endl;
    while(fgets(line,256,pfi) &&  !feof(pfi)) {
      sscanf(line,"%d %d %d %f %f %f %f %lf %lf",&a,&b,&x_tmp,&polI,&polQ,&polU,&polV,&y_tmp,&e_tmp);
      if (e_tmp > 0.0) {
        x.push_back(x_tmp*(2*M_PI/NB_PTS_PROF));
        // //if (y_tmp> 70.) y_tmp-=180.;
        Q.push_back(polQ);
        U.push_back(polU);
      }
    }

    fclose(pfi);
    cout << "Finished" << endl;
}

