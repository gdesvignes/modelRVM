
typedef struct param {
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

    // Params
    double *alpha;
    double *delta;
    double *Phi0;
    double *phi0;
    double *efac;
    int inc_fixed;
    double *inc;
    int prate_fixed;
    double *prate;

    // Data
    int numfiles;
    double **phs_exclude;
    int *n_phs_exclude;

} param;

int readParameters(param *p, char *paramFile);
