/* 
    This code sets up and evolves a system using the binary_c API.
    It is to be compiled as a shared library which can then be used
    in Python via the ctypes module. The output variables which are 
    saved at each timestep are defined in the function
    log_every_timestep.c in /src/logging.
*/

#include "binary_c.h"
#include <stdio.h>
#include <string.h>

void evolve_binary(char *, char *);

void evolve_binary(char * argstring, char * pythonbuffer)
{
    struct libbinary_c_stardata_t * stardata = NULL;
    struct libbinary_c_store_t * store = NULL;
    
    binary_c_new_system(&stardata,
                NULL,
                NULL,
                &store,
                &argstring,
                -1);

    char * temp = NULL;
    int nbytes = 0;

    stardata->preferences->internal_buffering = 2;
    stardata->preferences->internal_buffering_compression = 0;
    stardata->preferences->batchmode = BATCHMODE_LIBRARY;
    binary_c_evolve_for_dt(stardata, stardata->model.max_evolution_time);
    binary_c_buffer_info(stardata, &temp, &nbytes);
    
    // Copy output to string buffer allocated in python code
    strcpy(pythonbuffer, temp); 
    
    //* Free memory */
    binary_c_free_memory(&stardata, TRUE, TRUE, TRUE);
    stardata = NULL;
    store = NULL;
    temp = NULL;
}
