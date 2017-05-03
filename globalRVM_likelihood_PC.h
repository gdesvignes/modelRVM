
#include "RVMnest.h"

extern MNStruct *sp;

double globalRVMLogLike_PC(double theta[], int nDims, double phi[], int nDerived, void *context);
void prior (double cube[], double theta[], int nDims, void * context);
void setup_loglikelihood();
