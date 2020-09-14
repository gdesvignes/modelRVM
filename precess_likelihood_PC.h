#pragma once
#include "RVMnest.h"

extern MNStruct *sp;

double precessLogLike_PC(double theta[], int nDims, double phi[], int nDerived, void *context);
void precessprior (double cube[], double theta[], int nDims, void * context);
void setup_loglikelihood();
