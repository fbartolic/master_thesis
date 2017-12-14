#include "binary_c.h"
#include <stdio.h>

int main(int argc, const char *argv[])
{
    int n = 0;
    int nbytes = 0;
    char * buffer = NULL;
    char * argstring = "binary_c M_1 1.0 M_2 3.0 metallicity 0.02 orbital_period 10.0 max_evolution_time 100.0 eccentricity 0.0";

    struct libbinary_c_store_t * store = NULL;
    struct libbinary_c_stardata_t * stardata = NULL;

    while(n < 10000)
    {
        binary_c_new_system(&stardata,
                            NULL,
                            NULL,
                            &store,
                            &argstring,
                            -1);
        
        stardata->preferences->internal_buffering = 2;
        stardata->preferences->internal_buffering_compression = 0;
        stardata->preferences->batchmode = BATCHMODE_LIBRARY;
        binary_c_evolve_for_dt(stardata,stardata->model.max_evolution_time);
        binary_c_buffer_info(stardata, &buffer,&nbytes);

        /* Free memory */
        binary_c_free_memory(&stardata, TRUE, TRUE, TRUE);

        /*
         * The memory which was pointed to is freed,
         * but in order that the store is recreated the pointers
         * must be set to NULL.
         */
        store = NULL;
        stardata = NULL;
        //buffer = NULL; 

        n++;
    }


    return 0;
}
