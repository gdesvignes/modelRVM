
typedef struct param {
  // Sampler
  int sampler;

    // Multinest
    int IS;
    int nlive;
    int ceff;
    double efr;
    char *basename;
    
    // config
    int have_efac;  
    double threshold;
    int margin_phi0;  
    int have_aberr_offset;
    int sin_psi;

    // Params
    double *alpha;
    double *beta;
    double *delta;
    double *Phi0;
    double *phi0;
    double *psi0;
    double *efac;
    int inc_fixed;
    double inc;
    double *r_inc;
    int prate_fixed;
    double prate;
    double *r_prate;
    int njump;
    int psi_jump_fixed;
    double *psi_jumps;
    double *r_psi_jump;
    double *psi_jump_MJD;
    
    // Data
    int numfiles;
    double **phs_exclude;
    int *n_phs_exclude;
    
} param;

int readParameters(param *p, char *paramFile);
int readsimpleParameters(param *p, char *paramFile);
