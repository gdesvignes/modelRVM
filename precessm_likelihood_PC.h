#pragma once
#include "RVMnest.h"

extern MNStruct *sp;

double precessmLogLike_PC(double theta[], int nDims, double phi[], int nDerived, void *context);
void precessmprior (double cube[], double theta[], int nDims, void * context);
void setup_loglikelihood();
