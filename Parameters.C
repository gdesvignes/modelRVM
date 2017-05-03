#include <stdio.h>
#include <stdlib.h>
#include <glib.h>

#include "Parameters.h"

int readParameters(param *p, char *paramFile){
    GKeyFile* gkf; /* Notice we declared a pointer */

    char string[32];

    gkf = g_key_file_new();
    gsize length;

    if (!g_key_file_load_from_file(gkf, paramFile, G_KEY_FILE_NONE, NULL)){
        fprintf (stderr, "Could not read config file %s\n", paramFile);
        return EXIT_FAILURE;
    }
    p->sampler = g_key_file_get_integer(gkf,"sampler","sampler",NULL);

    p->IS = g_key_file_get_integer(gkf,"multinest","IS",NULL);
    p->nlive = g_key_file_get_integer(gkf,"multinest","nlive",NULL);
    p->ceff = g_key_file_get_integer(gkf,"multinest","ceff",NULL);
    p->efr = g_key_file_get_double(gkf,"multinest","efr",NULL);
    if (g_key_file_has_key(gkf,"multinest","basename",NULL))
        p->basename = g_key_file_get_string(gkf,"multinest","basename",NULL);

    p->threshold = g_key_file_get_double(gkf,"config","threshold",NULL);
    p->have_efac = g_key_file_get_integer(gkf,"config","have_efac",NULL);
    p->margin_phi0 = g_key_file_get_integer(gkf,"config","margin_phi0",NULL);
    p->have_aberr_offset = g_key_file_get_integer(gkf,"config","have_aberr_offset",NULL);

    p->alpha = g_key_file_get_double_list(gkf,"params","alpha", &length, NULL);
    p->delta = g_key_file_get_double_list(gkf,"params","delta", &length, NULL);
    p->Phi0 = g_key_file_get_double_list(gkf,"params","Phi0", &length, NULL);
    p->phi0 = g_key_file_get_double_list(gkf,"params","phi0", &length, NULL);
    p->inc_fixed = g_key_file_get_integer(gkf,"params","inc_fixed", NULL);
    p->inc = g_key_file_get_double(gkf,"params","inc",NULL); // Read inclination in deg
    p->r_inc = g_key_file_get_double_list(gkf,"params","r_inc", &length, NULL);
    p->prate_fixed = g_key_file_get_integer(gkf,"params","prate_fixed", NULL);
    p->prate = g_key_file_get_double(gkf,"params","prate",NULL);
    p->r_prate = g_key_file_get_double_list(gkf,"params","r_prate", &length, NULL);
    p->efac = g_key_file_get_double_list(gkf,"params","efac", &length, NULL);
    p->psi_jump_fixed = g_key_file_get_integer(gkf,"params","psi_jump_fixed", NULL);
    p->psi_jumps = g_key_file_get_double_list(gkf,"params","psi_jumps", &length, NULL);
    p->r_psi_jump =  g_key_file_get_double_list(gkf,"params","r_psi_jump", &length, NULL);
    p->psi_jump_MJD = g_key_file_get_double_list(gkf,"params","psi_jump_MJD", &length, NULL);
    p->njump = length;

    // Read the range of phases to be excluded
    p->phs_exclude = (double **) malloc( p->numfiles * sizeof(double *));
    p->n_phs_exclude = (int *) malloc( p->numfiles * sizeof(int));
    for (int i=0; i < p->numfiles; i++) {
	sprintf(string, "phs_exclude_%d", i);
	p->phs_exclude[i] = g_key_file_get_double_list(gkf,"data", string, &length, NULL);
	p->n_phs_exclude[i] = length / 2;
    }

    fprintf (stderr, "Finished reading config file %s\n", paramFile);

    g_key_file_free (gkf);

    return EXIT_SUCCESS;
}
