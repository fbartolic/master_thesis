/*
 * Subroutine to do all the logging at the end of a timestep.
 *
 * To put your own logging in here, choose a string (e.g. "FRED") and put
 * your log function in #defines, e.g.
 * 
 * #define FRED
 * 
 *     PRINTF("FRED %g %g %g\n",
 *            stardata->model.time,
 *            stardata->star[0].stellar_type,
 *            stardata->common.separation);
 *
 * #endif // FRED
 * 
 * You should use the PRINTF() binary_c macro instead
 * of printf(). This allows binary_c to (perhaps) compress
 * the output which is then sent to stdout (and binary_grid).
 *
 * Whether compression is faster, or not, depends on many 
 * factors. RGI's basic experiments suggest it is not, but
 * your experience on your hardware may differ.
 */

#include "../binary_c.h"

#include "buffering/printf_stack.h"

void log_every_timestep(struct stardata_t * RESTRICT stardata)
{
    /* 
     * No logging when RLOF-interpolating 
     * (negative timesteps are bad!).
     */
    if(stardata->model.dtm<-TINY) return; 
    
    /* timestep and dtp */
    Boolean first_timestep=IS_ZERO(stardata->model.time);
    if(first_timestep==TRUE) stardata->model.prev_log_t=0.0;
    double dt = stardata->model.time - stardata->model.prev_log_t;
    double dtp MAYBE_UNUSED = dt * stardata->model.probability;
    stardata->model.prev_log_t = stardata->model.time;
    
    /* FRAN'S OUTPUT CODE */
    PRINTF("%30.12e,%g,%g,%g,%g,%g,%g,%g,%d,%d,%d,%d,%d,",
             stardata->model.time,
             stardata->star[0].mass,
             stardata->star[1].mass,
             stardata->star[0].radius,
             stardata->star[1].radius,
             stardata->common.eccentricity,
             stardata->common.orbital_period,
             stardata->common.separation,
             stardata->star[0].stellar_type,
             stardata->star[1].stellar_type,
             stardata->model.in_RLOF,
             stardata->model.coel,
             stardata->model.merged
             );
    /* END FRAN's OUTPUT CODE */

    /* single star condition */
    Boolean single MAYBE_UNUSED = system_is_single(stardata); 
    Boolean binary MAYBE_UNUSED = NOT(single);

    if(0)PRINTF("TESTLOG %30.12e %g %g %d %g %g %g\n",
           stardata->model.time, // 1
           stardata->star[0].mass,//2
           stardata->star[0].core_mass,
           stardata->star[0].stellar_type,//4
           stardata->star[1].mass,
           stardata->common.separation ,//6
           stardata->common.eccentricity//7
        );
    
    if(0)PRINTF("TESTLOG t=%30.12e M1=%g Mc1=%g k1=%d M2=%g a=%g e=%g R=%g RL=%g dt=%g\n",
           stardata->model.time, // 1
           stardata->star[0].mass,//2
           stardata->star[0].core_mass,
           stardata->star[0].stellar_type,//4
           stardata->star[1].mass,
           stardata->common.separation ,//6

                  stardata->common.eccentricity//7
           ,stardata->star[0].radius//8
           ,stardata->star[0].roche_radius //9
           ,stardata->model.dt
        ); 

    //if(stardata->star[0].stellar_type>=6)
    if(0)printf("HRDHRD %d %g %g\n",
               stardata->star[0].stellar_type,
               TEFF(0),
               stardata->star[0].luminosity);

#if defined GAIA && defined STELLAR_COLOURS
    /*
     * Colours for Gaia
     */
    Star_number k;
    STARLOOP(k)
    {
        SETstar(k);

        /* Kurucz colours (Hurley) */
        double M_v,U_minus_B,B_minus_V,V_minus_I,M_b,M_u,M_i,M_bol;
        binary_colours(stardata,
                       stardata->store,
                       k,
                       &M_v,
                       &U_minus_B,
                       &B_minus_V,
                       &V_minus_I,
                       &M_b,
                       &M_u,
                       &M_i,
                       &M_bol,
                       NULL,
                       NULL,
                       NULL,
                       NULL);

        /* star information */
        char specstring[10]; // spectral type string
        spectral_type_string(specstring,stardata,star);
        PRINTF("t=%g dt=%g P=%g sn=%d st=%d M=%g logL=%g Mbol=%g logTeff=%g sp=%s\n",
               stardata->model.time,
               dt,
               stardata->model.probability,
               k,
               star->stellar_type,
               star->mass,
               log10(star->luminosity),
               M_bol,
               log10(TEFFSTAR(star)),
               specstring);

        /* Kurucz information */
        PRINTF("GAIA%d KURUCZ  U=% 8g B=% 8g V=% 8g R=******** I=% 8g\n",
               k,
               M_u,
               M_b,
               M_v,
               M_i);
               
        /* Eldridge colours */
        double magnitudes[NUMBER_OF_STELLAR_MAGNITUDES];
        char *cstring[NUMBER_OF_STELLAR_MAGNITUDES] = STELLAR_COLOUR_STRINGS;
        int c;
        
        /* 2012 */
        eldridge2012_colours(stardata,star,magnitudes);
        PRINTF("GAIA%d JJE2012 ",k);
        for(c=0;c<8;c++)
        {
            PRINTF("%s=% 8g ",
                   cstring[c],
                   magnitudes[c]);
        }
        PRINTF("\n");

        /* 2015 */
        eldridge2015_colours(stardata,star,magnitudes);
        gaia_colours(stardata,magnitudes,GAIA_CONVERSION_DEFAULT_METHOD);

        PRINTF("GAIA%d JJE2015 ",k);
        for(c=0;c<NUMBER_OF_STELLAR_MAGNITUDES;c++)
        {
            PRINTF("%s=% 8g ",
                   cstring[c],
                   magnitudes[c]);
        }
        PRINTF("\n");
    }
    

#endif // GAIA and STELLAR_COLOURS

    
//#ifdef NUCSYN 
//    nucsyn_binary_yield(stardata,FALSE);
//    RESET_NTP(0); 
//    RESET_NTP(1);
//#ifdef NUCSYN_GCE
//    /* GCE populations every timestep */
//    nucsyn_gce_log(stardata,FALSE);
//#endif
//#endif
  
#define MIMF
#ifdef MIMF
    if(stardata->star[0].stellar_type==TPAGB)
    {
        PRINTF("MCWD %g %g\n",
               stardata->star[0].pms_mass,
               stardata->star[0].core_mass);
    }
#endif // MIMF

#ifdef TPAGBTRACKS
    if(first_timestep)
    {
        PRINTF("TPAGBSTAR\n");
    }
    else if(stardata->star[0].stellar_type==TPAGB)
    {
        PRINTF("TPAGBSTAR %g %g\n",
               stardata->star[0].core_mass,
               stardata->star[0].mass);
    }
#endif // TPAGBTRACKS

 

#ifdef BLUE_STRAGGLER_PROJECT
    if(stardata->star[0].blue_straggler==TRUE ||
       stardata->star[1].blue_straggler==TRUE)
    {
        if(single)
        {
            PRINTF("BSS0 %g %g\n",
                   dtp,
                   stardata->model.time);
        }
        else
        {
            /* which star is the blue straggler? */
            int kbs= stardata->star[0].mass > stardata->star[0].pms_mass ? 0 : 1;
            struct star_t *bs = &(stardata->star[kbs]);
            struct star_t *comp = &(stardata->star[Other_star(kbs)]);

            PRINTF("BSS1 %g %g %g %g %g %d\n",
                   dtp,
                   stardata->model.time,
                   stardata->common.orbital_period*365.25,
                   bs->mass - bs->pms_mass, // mass accreted
                   comp->mass, // companion mass
                   comp->stellar_type // companion stellar type
                );
        }
    }
#endif //BLUE_STRAGGLER_PROJECT

#ifdef BLUE_STRAGGLER_BARIUM_PROJECT

    


    // start counting at OLD_TIME
#define OLD_TIME 10e3

    if(stardata->model.time > OLD_TIME)
    { 
	/* do not count back beyond OLD_TIME */
	double ddtp = MIN(stardata->model.time - OLD_TIME,dt) * 
	    stardata->model.probability;
	
	if(stardata->model.bss==TRUE)
	{
	    /* blue straggler statistics */

	    /* which star is the blue straggler? */
	    int kbs= stardata->star[0].mass > stardata->star[0].pms_mass ? 0 : 1;
	    struct star_t *bs = &(stardata->star[kbs]);

	    /* look for old, blue straggler dwarfs */
	    if(ON_MAIN_SEQUENCE(bs->stellar_type))
	    {
		//struct star_t *comp = &(stardata->star[Other_star(kbs)]);
		double *X = nucsyn_observed_surface_abundances(bs);
	
		if(0)
		    PRINTF("Ba=%g/%g, Fe=%g/%g\n",
			   nucsyn_elemental_abundance("Ba",X,stardata,stardata->store),
			   nucsyn_elemental_abundance("Ba",stardata->common.Xsolar,stardata,stardata->store),
			   nucsyn_elemental_abundance("Fe",X,stardata,stardata->store),
			   nucsyn_elemental_abundance("Fe",stardata->common.Xsolar,stardata,stardata->store)
			);
	       
	    
//	    PRINTF("BSSBa ddtp=%g t=%g P/days=%g dm=%g [C/Fe]=%g [Ba/Fe]=%g\n",
		PRINTF("BSSBa%d %g %g %g %g %g %g\n",
		       single ? 0 : 1,
		       ddtp,
		       stardata->model.time,
		       stardata->common.orbital_period*365.25,
		       bs->mass - bs->pms_mass, // mass accreted
		       nucsyn_elemental_square_bracket("C",
						       "Fe",
						       X,
						       stardata->common.Xsolar,
						       stardata),
		       nucsyn_elemental_square_bracket("Ba",
						       "Fe",
						       X,
						       stardata->common.Xsolar,
						       stardata)
		    );
	    }
	}
    
	/*
	 * MS dwarf statistics:
	 * count all MS stars, also MS-MS binaries
	 * but only log the primary (brightest) star
	 */
	struct star_t * star = NULL;
	if(ON_MAIN_SEQUENCE(stardata->star[0].stellar_type))
	{ 
	    // star 1 is a MS star
	    int j;
	    if(ON_MAIN_SEQUENCE(stardata->star[1].stellar_type))
	    { 
		// MS MS binary : log the "primary" == brightest star
	        j = stardata->star[0].luminosity > stardata->star[1].luminosity ? 0 : 1;
	    }
	    else
	    {
		// log only the MS star
		j = 1;
	    }
	    star = &(stardata->star[j]);
	}
	else if(ON_MAIN_SEQUENCE(stardata->star[1].stellar_type))
	{
	    // star 1 is not MS, but star 2 is
	    star = &(stardata->star[1]);
	}

	// proceed if we have a star to log
	if(star!=NULL)
	{
	    double *X = nucsyn_observed_surface_abundances(star);
	    PRINTF("MS%d %g %g %g %g %g %g\n",
		   single ? 0 : 1,
		   ddtp,
		   stardata->model.time,
		   stardata->common.orbital_period*365.25,
		   star->mass - star->pms_mass, // mass accreted
		   nucsyn_elemental_square_bracket("C",
						   "Fe",
						   X,
						   stardata->common.Xsolar,
						   stardata),
		   nucsyn_elemental_square_bracket("Ba",
						   "Fe",
						   X,
						   stardata->common.Xsolar,
						   stardata)
		);
	}
    }
#endif // BLUE_STRAGGLER_BARIUM_PROJECT

#ifdef WTTS_LOG
    wtts_log(stardata,dt);
#endif // WTTS_LOG

#if (defined(NUCSYN) && defined(YONG_PROJECT))
    {
        /*
         * Project to plot the magnesium ratios
         */
#define YONGOUT(K)                                      \
        if(stardata->star[(K)].stellar_type>1 &&        \
           stardata->star[(K)].stellar_type<7)          \
        {                                               \
            PRINTF("YONG%c %g %g\n",                    \
                   single?'s':'b',                      \
                   stardata->star[(K)].Xenv[XMg26]/    \
                   stardata->star[(K)].Xenv[XMg24],    \
                   dtp);                                \
        }

        YONGOUT(0);
        YONGOUT(1);
    }
#endif//YONG_PROJECT

#ifdef CANBERRA_PROJECT
    {
        /* output MS star */
        if(first_timestep==TRUE)
        {
            if(stardata->common.canberra_done_first_timestep==FALSE)
            {
                stardata->common.canberra_done_first_timestep=TRUE;

                // choose non-zero mass star
                int k=IS_ZERO(stardata->star[0].mass) ? 1 : 0;
                SETstar(k);
                PRINTF("CANBERRA_TEST%d %g %g %d",
                       1-single,
                       stardata->model.probability,
                       star->mass,
                       spectral_type(stardata,star));

                if(binary==TRUE)
                {
                    PRINTF(" %g %g %g sp2=%d\n",
                           stardata->star[1].mass/star->mass,
                           stardata->common.orbital_period*365.25,
                           stardata->common.eccentricity,
                           spectral_type(stardata,&stardata->star[1]));
                }
                else
                {
                    PRINTF("\n");
                }

                PRINTF("CANBERRA_ZAMS_SPECTRAL_TYPE %g %d %g %g\n",
                       star->mass,
                       spectral_type(stardata,star),
                       TEFFSTAR(star), 
                       long_spectral_type(stardata,star));
            }
        }
        else
        {
            stardata->common.canberra_done_first_timestep=FALSE;
        }

        int more_massive = 
            stardata->star[0].mass>stardata->star[1].mass ? 0 : 1;
        int more_massive_spectype = spectral_type(&stardata->star[more_massive]);

        /* MS (G/K) type binaries */
        if(binary==TRUE && 
           ON_MAIN_SEQUENCE(stardata->star[more_massive].stellar_type)&&
           (more_massive_spectype==SPECTRAL_TYPE_G ||
            more_massive_spectype==SPECTRAL_TYPE_K))
        {
            /* MS q distribution */

            /* select only G/K binaries */

            /* q = M1/M2 or M2/M1, whichever is less */
            double q = stardata->star[Other_star(more_massive)].mass/stardata->star[more_massive].mass;
            PRINTF("CANBERRAMSMSQ %g %g (from M1=%g M2=%g)\n",q,dtp,stardata->star[0].mass,stardata->star[1].mass);
        }


        /* WD system logging */
        if(WHITE_DWARF(stardata->star[0].stellar_type) ||
           WHITE_DWARF(stardata->star[1].stellar_type))
        {
            /* always log the WD first: find out which star it is */
            int k_wd = WHITE_DWARF(stardata->star[0].stellar_type) ? 0 : 1;
          
            if(single)
            {
                /* single WD : log its dtp */
                PRINTF("CANBERRA0 %g %g %d %d %g\n",
                       dtp,
                       stardata->model.time,
                       stardata->star[k_wd].stellar_type,
                       spectral_type(stardata,&stardata->star[k_wd]),
                       TEFF(k_wd)
                    );
            }
            else
            {
                /* binary with a white dwarf */
                PRINTF("CANBERRA1 %g %g %d %d %g %d %d %g %g %g\n",
                       dtp,//1
                       stardata->model.time, //2
                       stardata->star[k_wd].stellar_type,//3
                       spectral_type(stardata,&stardata->star[k_wd]),//4
                       TEFF(k_wd),//5
                       stardata->star[Other_star(k_wd)].stellar_type,//6
                       spectral_type(stardata,&stardata->star[Other_star(k_wd)]),//7
                       TEFF(Other_star(k_wd)),//8
                       stardata->common.orbital_period*365.25,//9
                       stardata->common.eccentricity//10
                    );
            }
        }
    }
#endif


#ifdef CARLO_OUTPUT //CAB
#ifndef NUCSYN
    PRINTF("testCarlo t= %1.5e  P= %1.5e  a= %1.5e  %i  m1= %1.5e  m2= %1.5e \n",
           stardata->model.time,
           stardata->common.orbital_period*365.25,
           stardata->common.separation,
           stardata->star[0].stellar_type,
           stardata->star[0].mass,
           stardata->star[1].mass
        );
#endif
//#ifdef NUCSYN
//    struct star_t * star = &(stardata->star[0]);
//      
//    Abundance * Xenv=nucsyn_observed_surface_abundances(star);
//    /***** t     type M1 Menv ntp  C  N  O  F  Ne  ******/
//    int k;
//  
//    STARLOOP(k)
//    {
//        PRINTF("CAB %1.4e %i %i %1.5f %1.5f %g %1.2f %1.2f %1.2f %1.2f %1.2f \n",
//               stardata->model.time,
//               k,
//               stardata->star[k].stellar_type,
//               stardata->star[k].mass,
//               stardata->star[k].core_mass,
//               stardata->star[k].num_thermal_pulses,
//               log10(stardata->star[k].Xenv[XC12]+stardata->star[k].Xenv[XC13]),
//               log10(stardata->star[k].Xenv[XN14]+stardata->star[k].Xenv[XN15]), 
//               log10(stardata->star[k].Xenv[XO16]+stardata->star[k].Xenv[XO17]+stardata->star[k].Xenv[XO18]), 
//               log10(stardata->star[k].Xenv[XF19]), 
//               log10(stardata->star[k].Xenv[XNe22]+stardata->star[k].Xenv[XNe20]+stardata->star[k].Xenv[XNe21])
//            );
//    }
//#endif // NUCSYN
#endif // CARLO_OUTPUT

  

#ifdef TESTPLOTXXX
    if(stardata->star[0].stellar_type<10)
    {
        if(stardata->model.in_RLOF && stardata->star[0].stellar_type<=2)
        {
            PRINTF("PLOT %g %g %g %g %g\n",       
                   stardata->star[0].mass,
                   stardata->star[0].radius,
                   stardata->star[0].roche_radius,
                   stardata->model.time,
                   -stardata->common.dM_RLOF_lose/stardata->model.the_timestep
                );  
        }
        else
        {
            PRINTF("PLOT %g %g %g %g %g\n",       
                   stardata->star[0].mass,
                   stardata->star[0].radius,
                   stardata->star[0].roche_radius,
                   stardata->model.time,0.0
                );  
        }}
#endif //TESTTPLOTXX

#ifdef NUCSYN

    if(0)
    {
        int k;
     
        STARLOOP(k)
        {
            if(stardata->star[k].stellar_type==TPAGB)
            {
                PRINTF("TPAGBPOP %d %g %g %g %g\n",
                       k,
                       stardata->star[k].pms_mass,
                       stardata->model.probability*dt,
                       stardata->star[Other_star(k)].pms_mass,
                       stardata->common.zams_separation
                    );
            }
        }
    }


    if(0)
    {
        // MCa=(MOo+MOa)*16/12-MCo
        struct star_t *star = &(stardata->star[0]);
        if(first_timestep==TRUE)

            PRINTF("carbon required %g\n",
                   (star->mass * star->Xenv[XO16]*12.0/16.0
                    -
                    star->mass * star->Xenv[XC12])
                );
    
        double vorb = 1e-5*sqrt ( GRAVITATIONAL_CONSTANT * 
                                  (stardata->star[0].mass+stardata->star[1].mass)*M_SUN/
                                  (stardata->common.separation*R_SUN));
    
        if(stardata->star[0].stellar_type==COWD &&
           stardata->star[1].roche_radius > stardata->star[1].radius)
        {
            PRINTF("vorb %g %g %g\n",
                   stardata->model.time,
                   vorb,

                   1e-5*sqrt(GRAVITATIONAL_CONSTANT * M_SUN * stardata->star[0].mass/(stardata->star[0].radius * R_SUN)) 

                );
        }
    }
#endif//NUCSYN

#ifdef TEST_NATURE
  
    {
        int k;
   
        if(dt>-TINY)
        {
            int spec[NUMBER_OF_STARS];
            double dtp=dt*stardata->model.probability;
            char *cc[]=SPECTRAL_TYPE_STRINGS;

            STARLOOP(k)
            {
                if(stardata->star[k].stellar_type<6)
                {
                    spec[k]=spectral_type(stardata,
                                          &stardata->star[k]);
                    //      PRINTF("SPEC TYPE %s -> stellar type %d\n",
                    //cc[spec[k]],stardata->star[k].stellar_type);
                }
                else
                {
                    spec[k]=-1; // ignore!
                }
            }

            /* select only systems with an O star */
            if((spec[1]==SPECTRAL_TYPE_O)||(spec[2]==SPECTRAL_TYPE_O))
            {
                /* single- or binary-star statistics */
                STARLOOP(k)
                {
                    if(spec[k]==SPECTRAL_TYPE_O)
                    {
                        PRINTF("OTYPE%d %d %g %g\n",
                               k,
                               stardata->model.sgl,
                               stardata->model.time,
                               dtp);
                    }
                }
    
                /* binary stars only */
                if(stardata->model.sgl==FALSE)
                {
                    if((spec[1]==SPECTRAL_TYPE_O)&&(spec[2]==SPECTRAL_TYPE_O))
                    {
                        /* O-O binary */
                        PRINTF("binOO %g %g\n",
                               stardata->model.time,
                               dtp);
                    }
                    else
                    {
                        /* O-else binary */
                        PRINTF("binOX %g %g\n",
                               stardata->model.time,
                               dtp);
                    }
                }
            }
        }
    }
#endif //TEST_NATURE

#if defined(SUPERNOVA_COMPANION_LOG)&&defined(LOG_SUPERNOVAE)
    if(stardata->star[0].went_sn_last_time!=0)
    {
        log_supernova_sc(stardata);
    }
#endif
  
    Dprint("t=%g m1=%g m1_0=%g m2=%g m2_0=%g\n",
           stardata->model.time,
           stardata->star[0].mass,stardata->star[0].phase_start_mass,
           stardata->star[1].mass,stardata->star[1].phase_start_mass);
    if(0)
    {
        kelvin_helmholtz_time(stardata);
        if(stardata->star[0].stellar_type<10)
        {
            double tdyn1=5.05e-5*sqrt(POW3(stardata->star[0].radius)/stardata->star[0].mass);
            double tdyn2=5.05e-5*sqrt(POW3(stardata->star[1].radius)/stardata->star[1].mass);

            PRINTF("TIMESCALES %g 1: %d dyn=%g kh=%g :2: %d dyn=%g kh=%g \n", 
                   stardata->model.time,
                   stardata->star[0].stellar_type,
                   tdyn1,
                   stardata->star[0].tkh,
                   stardata->star[1].stellar_type,
                   tdyn2,
                   stardata->star[1].tkh);
        }
    }



#ifdef FABIAN_IMF_LOG
    {

        Boolean fabian_log = FALSE;

        // force logging at time 0.0
        if (FEQUAL(stardata->model.time, 0.0) && (!stardata->common.fabian_logged_zero)) {
            /* first timestep */
            stardata->model.next_fabian_imf_log_time = 
                stardata->model.fabian_imf_log_timestep;
            fabian_log = TRUE;
            stardata->common.fabian_logged_zero = TRUE;
        }
        else
        {
            //stardata->common.fabian_logged_zero=FALSE;
        }
        if((MORE_OR_EQUAL(stardata->model.time, 
                          stardata->model.next_fabian_imf_log_time))) {
            // set next logging point and force logging
            stardata->model.next_fabian_imf_log_time = stardata->model.next_fabian_imf_log_time + stardata->model.fabian_imf_log_timestep;
            fabian_log = TRUE;
        }

        Dprint("Fabian IMF log t=%g next_fabian_imf_log_time=%g log_timestep=%g : do log? %d\n",
               stardata->model.time,
               stardata->model.next_fabian_imf_log_time,
               stardata->model.fabian_imf_log_timestep,
               fabian_log);

        if (fabian_log) 
        {
            // do logging if star is on MS
            if ((ON_MAIN_SEQUENCE(stardata->star[0].stellar_type)) || 
                (ON_MAIN_SEQUENCE(stardata->star[1].stellar_type)) )
            {
                PRINTF("FABIAN %16.12e %16.12e %16.12e %16.12e %16.12e %d %d %16.12e %16.12e %16.12e %16.12e %16.12e %16.12e %g %g %g\n", 
                       stardata->model.time,
                       stardata->star[0].mass,
                       stardata->star[1].mass,
                       stardata->star[0].luminosity,
                       stardata->star[1].luminosity,
                       stardata->star[0].stellar_type,
                       stardata->star[1].stellar_type,
                       stardata->common.separation,
                       stardata->common.eccentricity,
                       stardata->common.metallicity,
                       stardata->star[0].age,
                       stardata->star[1].age,
                       stardata->model.probability,
                       stardata->model.dt,
                       TEFF(0),
                       TEFF(1)
                    );
                //PRINTF("TEFF %g %g\n",TEFF(0),TEFF(1));
            }
        }
    }
  
#endif

#ifdef PERIOD_ECC_PLOT
    {
    
        if((stardata->star[0].stellar_type<10)&&
           (stardata->star[1].stellar_type<10)&&
           (stardata->model.sgl==FALSE))
        {
            PRINTF("PE %d %g %g\n",
                   stardata->star[0].stellar_type,
                   stardata->common.orbital_period*365.25,
                   stardata->common.eccentricity);
        }
    
        if(stardata->star[0].stellar_type>6)
        { 
            Exit_binary_c(NORMAL_EXIT,"Post AGB star : exit\n");
        }
    }
#endif

//#ifdef SINGLE_STAR_LIFETIMES
//    if(stardata->star[0].stellar_type < HERTZSPRUNG_GAP)
//    {
//        stardata->common.done_single_star_lifetime = FALSE;
//    }
//    else if(stardata->common.done_single_star_lifetime==FALSE &&
//            LATER_THAN_WHITE_DWARF(stardata->star[0].stellar_type))
//    {
//        PRINTF("SINGLE_STAR_LIFETIME %g %g\n",
//               stardata->star[0].pms_mass,
//               stardata->model.time);
//        fflush(stdout);
//        stardata->common.done_single_star_lifetime = TRUE;
//    }
//#endif

#ifdef FILE_LOG
    {
#ifdef PRE_MAIN_SEQUENCE
	if(stardata->preferences->pre_main_sequence)
	{
	    
            if(first_timestep) stardata->common.pms[0]=stardata->common.pms[1]=TRUE;
	    
	    Star_number k;
	    STARLOOP(k)
	    {
		if(stardata->common.pms[k]==TRUE &&
                   stardata->star[k].PMS==FALSE)
                { 
                    Append_logstring(LOG_PREMAINSEQUENCE," ZAMS %d",k);
                    stardata->common.pms[k]=FALSE;
                }
            }
	}
#endif // PRE_MAIN_SEQUENCE
    }
#endif //FILE_LOG

#ifdef NUCSYN
    //nucsyn_short_log(stardata);
#endif

#ifdef NANCHECKS
    Dprint("execute not-a-number checks\n");
#define NCHECK(A,B,C) if(isnan((A))!=0){PRINTF("NaN %s star %d\n",(B),(C));fflush(stdout);fflush(stderr);kill(0,SIGSEGV);Exit_binary_c(EXIT_NAN,"NaN detected");}
    {
        int k;
        STARLOOP(k)
        {
            SETstar(k);
            NCHECK(star->mass,"Mass (evolution)",k);
            NCHECK(star->radius,"Radius (evolution)",k);
            NCHECK(star->luminosity,"Luminosity (evolution)",k);
            NCHECK(star->angular_momentum,"jspin",k);
        }
    }
    if(stardata->model.sgl == FALSE)
    {
        NANCHECK(stardata->common.separation);
        NANCHECK(stardata->common.eccentricity);
    }
    NANCHECK(stardata->model.time);
#endif

#ifdef CONSERVATIVE_SYSTEM
    {
        /* conservative systems should not lose mass */
        double dm=
            fabs(stardata->star[0].mass+stardata->star[1].mass - 
                 stardata->star[0].pms_mass - stardata->star[1].pms_mass);
        Dprint("Conservative system check : %g\n",dm); 

        if(dm>0.1)
        {
            Exit_binary_c(OUT_OF_RANGE,
                          "System is non-conservative when it should be conservative! M1=%g M2=%g : total = %g : should be %g\n",
                          stardata->star[0].mass,
                          stardata->star[1].mass,
                          stardata->star[0].mass+stardata->star[1].mass,
                          stardata->star[0].pms_mass+stardata->star[1].pms_mass);
            
        }
    }

#endif

#ifdef NUCSYN_LOWZ_SUPERNOVAE

    if((stardata->common.nucsyn_metallicity>-TINY)&&
       (stardata->common.nucsyn_metallicity<NUCSYN_LOWZ_SNE_THRESHOLD))
    {
        /* Use zero-Z supernovae for both stars */
        double *X;
        
        X=CALLOC(ISOTOPE_ARRAY_SIZE,sizeof(double));
        stardata->model.time=0.0;
        nucsyn_binary_yield(stardata,TRUE);
        stardata->model.time=4.000; // Assume output at 4Myr
#ifdef FILE_LOG
        output_string_to_log(stardata,"PopIII stars\n");
#endif

        // explode star 1
        nucsyn_lowz_yields(X,stardata->star[0].mass);
        nucsyn_calc_yields(stardata,&(stardata->star[0]),1.0,X,0,X,1,TRUE,NUCSYN_SOURCE_POPIII);
        stardata->star[0].mass=0.0;
        stardata->star[0].stellar_type=MASSLESS_REMNANT;
        // explode star 2
        nucsyn_lowz_yields(X,stardata->star[1].mass);
        nucsyn_calc_yields(stardata,&(stardata->star[1]),1.0,X,0,X,2,TRUE,NUCSYN_SOURCE_POPIII);
        stardata->star[1].mass=0.0;
        stardata->star[1].stellar_type=MASSLESS_REMNANT;

        /* Force output at 4Myr */
        nucsyn_binary_yield(stardata,TRUE);
        stardata->model.time=4.001; 
        end_of_iteration(stardata);
        nucsyn_binary_yield(stardata,TRUE);
        stardata->model.time=stardata->model.max_evolution_time;
        nucsyn_binary_yield(stardata,TRUE);

#ifdef FILE_LOG
        close_log_files(&(model->log_fp));
#endif
        SAFE_FREE(X);
        return(0);
    }
#endif //NUCSYN_LOWZ_SUPERNOVAE 


#ifdef NUCSYN_NETWORK_TEST
    nucsyn_network_test(stardata);
#endif

#if defined ADAPTIVE_RLOF && defined ADAPTIVE_RLOF_LOG
    /* new RLOF log */
    if(//(stardata->preferences->RLOF_method!=0)&&
        ((stardata->star[0].stellar_type<10)||
         (stardata->star[1].stellar_type<10)))
    {
        /* find RLOFing star */
        int krl=0;

        if((stardata->model.time>103.0) && (stardata->model.time<105.0))
        {
            krl=1;
        }
        else
        {
            if(stardata->star[0].radius>stardata->star[0].roche_radius*0.99)
            {
                krl=1;
            }
            else if(stardata->star[1].radius>stardata->star[1].roche_radius*0.99)
            {
                krl=2;
            }
        }

        if(krl>0)
        {
            /* set thermal timescales and log */
            kelvin_helmholtz_time(stardata);
            struct star_t * donor = &(stardata->star[krl]);
            struct star_t * accretor = &(stardata->star[Other_star(krl)]);
            PRINTF("RLOFLOG %20.12e %g %g %g %g %g %g %g %g %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                   stardata->model.time, // 1
                   stardata->common.mdot_RLOF, // 2 what we got
                   stardata->common.mdot_RLOF_adaptive, // 3 what the adaptive wants
                   donor->mass/donor->tkh, // 4 thermal rate
                   donor->radius, // 5 
                   donor->roche_radius, // 6
                   donor->radius/donor->roche_radius, // 7
                   stardata->star[0].mass, //8
                   stardata->star[1].mass, // 9
                   krl, // 10 (number of rlofing star)
                   donor->rdot*1e-6, // 11 da/dt (t=Myr)
                   donor->roldot, // 12 dRL/dt (t=yr)
                   stardata->common.separation, //13 (Rsun)
                   donor->mass, // 14
                   stardata->model.dt, // 15 dt (years)
                   stardata->model.dtm*1e6,   // 16 dtm (years)
                   accretor->radius, // 17
                   accretor->roche_radius, // 18
                   donor->nth, // 19
                   MAX(0.0,donor->radius/donor->roche_radius), // 20
                   donor->radius, //21
                   accretor->radius, // 22
                   donor->luminosity, //23
                   accretor->luminosity, // 24
                   stardata->common.orbital_period*365.25, // 25
                   stardata->common.orbital_angular_momentum // 26

                );
        }
        else
        {
            //PRINTF("RLOFLOG * * * * * * * * * * * * * * * * * * * * * * *\n");
        }
        stardata->common.mdot_RLOF=0.0;
        stardata->common.mdot_RLOF_adaptive=0.0;

    }

    {
      

        double dM_RLOF_losedt;
        double dt=(1e6*stardata->model.dtm);
        if(stardata->model.time<TINY)
        {
            dM_RLOF_losedt=0.0;
        }
        else if(dt>TINY)
        {
            dM_RLOF_losedt=(stardata->common.RLOF_log_m1was-stardata->star[0].mass)/dt;
        }
        else
        {
            dM_RLOF_losedt=0.0;
        }
      
     
        if(dt>TINY)
        {
            PRINTF("MDOT %10.10e ",stardata->model.time); // 1
            if(stardata->model.in_RLOF==TRUE)
            {
                PRINTF("%g %g",
                       stardata->common.mdot_RLOF_adaptive, // 2
                       stardata->common.mdot_RLOF_H02 /* 3 */);
            }
            else
            {
                PRINTF("* *");
            }
            /* thermal and dynamical rates (years) */
            double tkh1=1e7*POW2(stardata->star[0].mass)/(stardata->star[0].radius * stardata->star[0].luminosity);
            double tdyn1=5.05e-5*sqrt(POW3(stardata->star[0].radius)/stardata->star[0].mass);

            /* omega crit (s^-1) */
            /* code omega = (seconds per year) * v (cm/s) / r(cm) */
            /* code omega units = per year */
            /* so to convert to cgs divide by YEAR_LENGTH_IN_SECONDS */
            calculate_rotation_variables(&(stardata->star[0]));
            calculate_rotation_variables(&(stardata->star[1]));

          
            //        4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39
            PRINTF(" %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                   dM_RLOF_losedt,//4
                   stardata->star[0].mass,//5
                   stardata->star[1].mass,
                   stardata->star[0].radius / MAX(1e-14,stardata->star[0].roche_radius),
                   stardata->star[1].radius / MAX(1e-14,stardata->star[1].roche_radius),
                   HRdiagram(0),//9,10  - macro (two elements)
                   HRdiagram(1),//11,12
                   stardata->common.separation,//13
                   stardata->star[0].v_eq,// 14, equatorial velocity km/s
                   stardata->star[1].v_eq,// 15
                   stardata->star[0].v_eq_ratio,//16
                   stardata->star[1].v_eq_ratio,//17
                   stardata->star[0].phase_start_mass,//18
                   stardata->star[1].phase_start_mass,//19
                   stardata->star[0].radius,//20
                   stardata->star[1].radius,//21
                   stardata->star[0].mass/tdyn1,//22
                   stardata->star[0].mass/tkh1,//23
                   stardata->star[0].roche_radius,//24
                   stardata->star[1].roche_radius,//25
                   stardata->common.used_rate,//26
                   stardata->common.runaway_rlof,//27
                   stardata->common.thermal_cap,//28
                   stardata->star[0].stellar_type,//29
                   stardata->star[1].stellar_type,//30
                   stardata->star[0].angular_momentum,//31
                   /* wind enhancement factors */
                   stardata->star[0].rotfac,//32
                   stardata->star[1].rotfac,//33
                   stardata->star[0].v_sph, //34
                   stardata->star[1].v_sph,//35
                   stardata->star[0].v_sph_ratio, //36
                   stardata->star[1].v_sph_ratio,//37
                   stardata->star[0].omega,//38
                   stardata->star[1].omega,//39
                   stardata->star[0].omega_ratio,//40
                   stardata->star[1].omega_ratio,//41
                   stardata->star[0].radius, //42
                   stardata->star[1].radius, //43
                   stardata->star[0].luminosity, //44
                   stardata->star[1].luminosity //45
                );

            // reset RLOF parameters
            stardata->common.used_rate = RLOF_RATE_UNDEF;
            stardata->common.runaway_rlof=FALSE;
            stardata->common.thermal_cap=FALSE;
          
        }
        stardata->common.RLOF_log_m1was=stardata->star[0].mass;
    }
#endif //  ADAPTIVE_RLOF && ADAPTIVE_RLOF_LOG



#if DEBUG==1 && defined NUCSYN
    {
        Abundance * Xenv=nucsyn_observed_surface_abundances(&(stardata->star[1]));
        if(stardata->star[1].stellar_type<10)
        {
            Dprint("AB%d %g %g %d %g %g %g\n",
                   2,
                   stardata->model.time,
                   stardata->star[1].mass,
                   stardata->star[1].stellar_type,
                   Xenv[XC12],
                   Xenv[XC13],
                   Xenv[XN14]);
        }
    }
#endif

#ifdef NS_NS_BIRTH_LOG
    if(stardata->common.nsns==FALSE)
    {
        if((stardata->star[0].stellar_type==NEUTRON_STAR)&&
           (stardata->star[1].stellar_type==NEUTRON_STAR))
        {
            PRINTF("NASCENT_NSNS at t=%g p=%g\n",
                   stardata->model.time,
                   stardata->model.probability);
            stardata->common.nsns=TRUE;
        }
    }

    /* first timestep: restore nsns birth boolean to FALSE */
    if(first_timestep==TRUE)
    {
        stardata->common.nsns=FALSE;
    }
#endif
 
#ifdef FINAL_MASSES_LOG
    /* if star 1 is post-AGB then stop */
    if(stardata->star[0].stellar_type>TPAGB)
    {
        stardata->model.max_evolution_time=      stardata->model.time ;
        return;
    }
#endif
#ifdef LOG_COMENV_RSTARS
    STARLOOP(k)
    {
        log_r_star(&(stardata->star[k]),stardata);
    }
#endif
#ifdef LOG_HERBERT
    {
        struct star_t * star=&(stardata->star[0]);

        if((star->stellar_type<10)&&(star->stellar_type!=6))
        {

            double mch=0.0,mche=0.0,lhe=0.0;

            /* H-burning shell */
            if((star->stellar_type>=HG)&&(star->stellar_type<=TPAGB))
            {
                mch=star->core_mass;
            }
            /* He-burning shell */
            if((star->stellar_type>=EAGB)&&(star->stellar_type<=HeGB))
            {
                mche=star->core_mass;
            }
            /* He luminosity : guess */
            if((star->stellar_type==CHeB)||((star->stellar_type>=HeMS)&&(star->stellar_type<=HeGB)))
            {
                lhe=star->luminosity;
            } 
      
            PRINTF("HERBERT 0 %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                   stardata->model.time,
                   log10(MAX(1e-20,star->radius)),
                   log10(MAX(1e-20,TEFF(0))), // (4)
                   log10(star->luminosity),
                   star->mass, // (6)
                   mch,
                   mche, // (8)
                   log10(MAX(1e-20,lhe)),
                   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, // convective boundaries 10-21
                   mch,mche, // masses of max energy generation
                   0.0, // (24) opacity
                   stardata->common.metallicity, // metallicity
                   stardata->model.dt, // 26 timestep
                   star->menv, // convective envelope
                   0.0,0.0,0.0,0.0,0.0, // H,He,C,N,O
                   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, // 33-44 burning shell boundaries
                   0.0,0.0,0.0,0.0,0.0,0.0,0.0, // 45-51 LH,C,th,nu,pp,CNO,pp/H
                   0.0,0.0,0.0,0.0,0.0,0.0,0.0); // 51-58 core H,He,C,N,O,T,density
        }}
#endif

#ifdef NUCSYN
#ifdef LOG_BARIUM_STARS
    log_barium_stars(stardata);
#endif
#endif 


#ifdef NUCSYN
#ifdef RS_LOG_ACCREITON_LOG
    if((stardata->star[0].stellar_type==TPAGB)&&
       (stardata->star[1].stellar_type<=MAIN_SEQUENCE))
    {
        /* Choose to output (log) abundance or square bracket */ 

#define ACCABUND(A) (stardata->star[1].Xacc[(A)])
        

        PRINTF("ACCLAYER t=%g Macc2=%g ISOTOPES ",
               stardata->model.time,
               stardata->star[1].dmacc         
            );

        Isotope a;
        Ordered_isotope_loop(a)
        {
            PRINTF("%g ",ACCABUND(a));
        }
        PRINTF("\n");
        
    }

#endif
#endif


#ifdef TEST_AXEL_RGBF_FIX
    if(stardata->star[0].stellar_type==GIANT_BRANCH)
    {
        PRINTF("RADIUS %g %g %g %g\n",
               stardata->model.time,
               (1000.0*pow((1130.0*stardata->star[0].luminosity/POW2(stardata->star[0].radius)),0.25)),
               stardata->star[0].radius,
               stardata->star[0].luminosity);
    }
#endif


    Dprint("ITERATE t=%g sep=%6.8lf intpol=%d dtm=%g maxt=%g\n",
           stardata->model.time,
           stardata->common.separation,
           stardata->model.intpol,
           stardata->model.dtm,
           stardata->model.max_evolution_time);

#ifdef NUCSYN_GIANT_ABUNDS_LOG
    giant_abundance_log(stardata);
#endif //NUCSYN_GIANT_ABUNDS_LOG

#ifdef ANG_MOM_CHECKS 
    determine_roche_lobe_radii(stardata); // Update RLOBes 
    PRINTF("ANGMOM %g %d %g %g %d %g %g %g %g %g %g %g %g %g %g %g\n",stardata->model.time, 
           stardata->star[0].stellar_type, 
           MAX(0,stardata->star[0].angular_momentum*M_SUN*R_SUN*R_SUN/YEAR_LENGTH_IN_SECONDS), 
           MAX(0,2.0*PI*stardata->star[0].omega/YEAR_LENGTH_IN_SECONDS), 
           stardata->star[1].stellar_type, 
           MAX(0,stardata->star[1].angular_momentum*M_SUN*R_SUN*R_SUN/YEAR_LENGTH_IN_SECONDS), 
           MAX(0,2.0*PI*stardata->star[1].omega/YEAR_LENGTH_IN_SECONDS), 
           stardata->common.separation, 
           stardata->common.orbital_period, 
           stardata->star[0].mass, 
           stardata->star[1].mass, 
           stardata->star[0].radius/stardata->star[0].roche_radius, 
           stardata->star[1].radius/stardata->star[1].roche_radius, 
           stardata->star[0].dmr, // dmr=wind 
           stardata->star[1].dmr, 
           stardata->star[0].dmr+stardata->common.dM_RLOF_transfer /* wind + RLOF */);  
#endif/* ANG_MOM_CHECKS */ 



#ifdef STELLAR_COLOURS

    if(0){
    double colours[NUMBER_OF_STARS][9]; /* colour data for each star */
    double l[NUMBER_OF_STARS][5]; /* luminosities */ 

//        if(stardata->star[0].stellar_type<=1)
    if(stardata->model.time<TINY)
    {
        /* Obtain colours for star k */
        Star_number k;
        STARLOOP(k)
        {
            binary_colours(stardata,
                           stardata->store,
                           k,
                           colours[k],
                           colours[k]+1,
                           colours[k]+2,
                           colours[k]+3,
                           colours[k]+4,
                           colours[k]+5,
                           colours[k]+6,
                           colours[k]+7,
                           l[k],
                           l[k]+1,
                           l[k]+2,
                           l[k]+3);
                
            PRINTF("STELLAR_COLOURS# Mu Mb Mv Mi  U-B B-V V-I Mbol (star %d)\n",k);
            PRINTF("STELLAR_COLOURS %g %g %g %g %g %g %g %g\n",
                   colours[k][5],colours[k][4],colours[k][0],colours[k][6],
                   colours[k][1],colours[k][2],colours[k][3],colours[k][7]);
                
            PRINTF("\n\n");
                

        }

        /* calculate M_x for the binary */
        if(stardata->star[1].mass<=stardata->star[0].mass)
        {
            PRINTF("BINARY_BINARY_COLOURS %g %g %g %g %g %g %g\n",
                   6.34482 - 2.5 * log10(l[1][0]+l[2][0]),
                   6.07549 - 2.5 * log10(l[1][1]+l[2][1]),
                   5.35129 - 2.5 * log10(l[1][2]+l[2][2]),
                   4.60777 - 2.5 * log10(l[1][3]+l[2][3]),
                   stardata->star[0].mass,
                   stardata->star[1].mass,
                   stardata->common.separation
                );
        }
        /* skip the rest of the code */
        stardata->model.time=16000;
            
    }
    }
#endif


#ifdef HRDIAG
    log_hr(stardata);
#endif

    //#define STAR2_MONITOR
#ifdef STAR2_MONITOR
    if(stardata->star[1].stellar_type<10)
    {
        PRINTF("STAR2 m=%g (accreted %g) kw=%d Xenv=%g ",
               stardata->star[1].mass,
               stardata->star[1].mass-stardata->star[1].zams_mass,
               stardata->star[1].stellar_type,
               stardata->star[1].Xenv[XC12]);
        if(stardata->star[1].dmacc>0.0)
        {
            PRINTF(" dmacc=%g Xacc=%g",
                   stardata->star[1].dmacc,
                   stardata->star[1].Xacc[XC12]);
        }
        else
        {
            PRINTF(" no acc layer ");
        }
        PRINTF("\n");
    }
#endif
#ifdef DETAILED_LOG
    detailed_log(stardata);
#endif /* DETAILED_LOG */

#ifdef XRAY_BINARIES
    // derivative[DERIVATIVE_STELLAR_MASS_WIND_GAIN] is msun /yr so multiply by M_SUN / YEAR_LENGTH_IN_SECONDS
    // to get into CGS
    log_xray_binary(stardata,2,1,stardata->star[0].derivative[DERIVATIVE_STELLAR_MASS_WIND_GAIN]);
    log_xray_binary(stardata,1,2,stardata->star[1].derivative[DERIVATIVE_STELLAR_MASS_WIND_GAIN]);
#endif 

#ifdef STELLAR_TYPE_LOG
    stellar_type_log(stardata);
#endif
        
           
#if defined NUCSYN && defined NUCSYN_AND_HRD
    if(stardata->star[0].stellar_type<10)
    {
        struct star_t * star = &(stardata->star[0]);
        PRINTF("HRDCHEM prob=%g ZAMSmass=%g logTeff=%g logL=%g H1=%g He4=%g Li7=%g C12=%g C13=%g N14=%g O16=%g Na23=%g Al26=%g CO=%g Ba=%g Pb=%g\n",
               stardata->model.probability,
               star->pms_mass,
               log10(TEFF(0)),
               log10(star->luminosity),
               star->Xenv[XH1],
               star->Xenv[XHe4],
               star->Xenv[XLi7],
               star->Xenv[XC12],
               star->Xenv[XC13],
               star->Xenv[XN14],
               star->Xenv[XO16],
               star->Xenv[XNa23],
               star->Xenv[XAl26],
               nucsyn_elemental_abundance("C",star->Xenv,stardata,stardata->store)/
               nucsyn_elemental_abundance("O",star->Xenv,stardata,stardata->store),
               nucsyn_elemental_abundance("Ba",star->Xenv,stardata,stardata->store),
               nucsyn_elemental_abundance("Pb",star->Xenv,stardata,stardata->store)
            );
    }
#endif

    //#define TEST_JJES_RATES
#ifdef TEST_JJES_RATES

    /* test JJE's wind loss rates */

    for(stardata->star[0].Xenv[XH1]=0.7;
        stardata->star[0].Xenv[XH1]>0;
        stardata->star[0].Xenv[XH1]-=0.1)
    {
        stardata->common.metallicity=0.003;
        stardata->star[0].Xenv[XO16]=0.2;
        stardata->star[0].Xenv[XC12]=0.2;
        stardata->star[0].Xenv[XH1]=0.005;
        stardata->star[0].Xenv[XHe4]=1.0-
            stardata->star[0].Xenv[XH1]-
            stardata->star[0].Xenv[XC12]-
            stardata->star[0].Xenv[XO16]-
            stardata->common.metallicity;
            
        for(stardata->star[0].luminosity=100;
            stardata->star[0].luminosity<=10000;
            stardata->star[0].luminosity+=100)
        {
            for(teff=4000;
                teff<=40000;
                teff+=1000)
            {
                /*      stardata->star[0].mass=20;
                        stardata->star[0].luminosity=10000;
                        teff=12603.9004;
                */
                stardata->star[0].radius=sqrt(1130.0*stardata->star[0].luminosity)/POW2(teff/1e3);

                PRINTF("%g %g %g %g\n",
                       stardata->star[0].luminosity,
                       1000.0*pow((1130.0*stardata->star[0].luminosity/
                                   POW2(stardata->star[0].radius)),
                                  0.25),
                       stardata->star[0].Xenv[XH1],
                       calc_stellar_wind_mass_loss(1,
                                                   stardata->star[0].luminosity,
                                                   stardata->star[0].radius,
                                                   20,
                                                   1,
                                                   1e6*stardata->star[0].radius,
                                                   stardata->common.metallicity,
                                                   0.0,
                                                   stardata,
                                                   1)
                    );
            }
            PRINTF("\n");
        }
        Exit_binary_c(NORMAL_EXIT,"Test JJE rates exit\n");
    }
    Exit_binary_c(NORMAL_EXIT,"Test JJE rates exit\n");
#endif

#ifdef SELMA
    STARLOOP(k)
    {
        if((stardata->star[k].radius>=stardata->star[k].roche_radius)
           &&
           (stardata->star[Other_star(k)].radius<stardata->star[Other_star(k)].roche_radius)
           && 
           (stardata->star[Other_star(k)].stellar_type==MAIN_SEQUENCE)
           &&
           (stardata->star[k].stellar_type==MAIN_SEQUENCE)
            )
        {
            PRINTF("SELMA %d %g %g %g %g %g %g %g\n",
                   k,
                   stardata->model.time,
                   stardata->star[k].radius,
                   stardata->star[Other_star(k)].radius,
                   stardata->star[k].mass,
                   stardata->star[Other_star(k)].mass,
                   stardata->common.separation,
                   stardata->common.orbital_period
                );
        }
    }
#endif //SELMA

#ifdef LTEST
        
    if(stardata->star[0].stellar_type<10)
    {
        PRINTF("LTEST %g %g %g %g %g %g %g %g %g %g %g %g\n",
               1e6*stardata->model.time,
               log10(stardata->star[0].luminosity),
               log10(1000.0*pow((1130.0*stardata->star[0].luminosity/
                                 POW2(MAX(stardata->star[0].radius,TINY))),
                                0.25)),
               stardata->star[0].mass,
               stardata->star[0].Xenv[XH1],
               stardata->star[0].Xenv[XHe4],
               stardata->star[0].Xenv[XC12],
               stardata->star[0].Xenv[XN14],
               stardata->star[0].Xenv[XO16],
               ((stardata->star[0].Xenv[XC12]/12.0+stardata->star[0].Xenv[XC13]/13.0+stardata->star[0].Xenv[XO16]/16.0)/(1e-14+stardata->star[0].Xenv[XHe4]/4.0))     ,
               log10(stardata->star[0].radius),
               stardata->star[0].core_mass
                   
            );
    }
#endif


#if defined CGI_EXTENSIONS && defined CGI_EXTENSIONS_HRD
    /*
     * HRD for cgi (web) output
     */
    if((stardata->model.log_fp!=NULL)&&
       (stardata->model.time>1e-5))
    {
      fprintf(stardata->model.log_fp,
	      "CGIHRD %d %g %g %d %g %g\n",
	      stardata->star[0].stellar_type,
	      TEFF(0),
	      stardata->star[0].luminosity,
	      stardata->star[1].stellar_type,
	      TEFF(1),
	      stardata->star[1].luminosity);
    }
#endif
    
#if defined NUCSYN && defined CGI_EXTENSIONS \
  && defined CGI_EXTENSIONS_ABUNDS && defined FILE_LOG

    /* Extra stuff to the logfile for cgi (web) output */
    if((stardata->model.log_fp!=NULL)&&
       (stardata->model.time>1e-5))
    {
        fprintf(stardata->model.log_fp,"ABUNDS %10.10e ",stardata->model.time);
        Star_number j;
        double *X;
	STARLOOP(j)
        {
            X=nucsyn_observed_surface_abundances(&(stardata->star[j]));
            fprintf(stardata->model.log_fp,
                    " %3.3e %3.3e %3.3e %3.3e %3.3e %3.3e %3.3e",
                    X[XH1],X[XHe4],X[XC12],X[XN14],X[XO16],
                    X[XFe56],
                    nucsyn_elemental_abundance("Ba",X,stardata,stardata->store));
        }
        fprintf(stardata->model.log_fp,"\n");
    }
#endif/* NUCSYN */

    //#define NUCSYN_ABUNDANCE_LOG
#ifdef NUCSYN_ABUNDANCE_LOG

    STARLOOP(k)
    {
        if(stardata->star[k].stellar_type<10)
        {
            Abundance * X=nucsyn_observed_surface_abundances(&(stardata->star[k]));
            PRINTF("ABUNDS %g %d %d ",
                   stardata->model.time,
                   k,
                   stardata->star[k].stellar_type);
            Isotope ii;
            Ordered_isotope_loop(ii)
            {
                PRINTF("%g ",X[ii]);
            }
            PRINTF("\n");
        }
    }
#endif

#ifdef NUCSYN_CEMP_LOGGING
    nucsyn_cemp_log(stardata);
#endif

#ifdef SDB_CHECKS
    sdb_check(stardata);
#endif
    /* and some logging */
#ifdef NUCSYN_J_LOG
    nucsyn_j_log(stardata);
#endif
#ifdef NUCSYN_SHORT_LOG
    nucsyn_short_log(stardata);
#endif
#ifdef NUCSYN_LONG_LOG
    nucsyn_long_log(stardata);
#endif
#ifdef LOG_JJE
    log_jje(stardata);
#endif
#ifdef NUCSYN_LOG_JL    
    log_jl(stardata);
#endif


#ifdef LOG_RSTARS
    /*
     * Log R-stars
     */
    double teff,dt=stardata->model.dtm;
    double min_thick_disk_age=thick_disk_age(stardata->common.metallicity);
        
    if(stardata->model.time<TINY)
    {
        stardata->model.rstar_stardata->model.prev_log_time=0.0;
    }
    dt=stardata->model.time-stardata->model.rstar_stardata->model.prev_log_time;
    stardata->model.rstar_stardata->model.prev_log_time=stardata->model.time;
    Boolean GKgiant=FALSE;
    /*
     * Convert C-rich CHeB stars into R-stars
     */
#ifdef NUCSYN
    STARLOOP(k)
    {
        if((stardata->star[k].stellar_type==CHeB)&&
           ((stardata->star[k].Xenv[XC12]/12.0)/
            (stardata->star[k].Xenv[XO16]/16.0)>1.0)
           &&(stardata->star[k].rstar==0)
            )
        {
            // should be an r-star!
            stardata->star[k].rstar=4;
        }
    }
#endif

    /*
     * Star must be older than the minimum age of the thick disk
     */
    if((stardata->model.time>min_thick_disk_age)&&
       (dt>0.0)&&
       (stardata->model.time>0.0))
    {

        STARLOOP(k)
        {
            teff = 1000.0*pow((1130.0*stardata->star[k].luminosity/
                               POW2(stardata->star[k].radius)),
                              0.25);


            if((stardata->star[k].rstar>0)
               &&(stardata->star[k].stellar_type==CHeB)
                )
            {
                /* calculate radial velocity parameter K1 */
                double v=radial_velocity_K(stardata,PI/2.0,k);
                    
                PRINTF("RSTAR%d %g %g %g %g %g %g %g %g %g\n",
                       stardata->star[k].rstar,
                       stardata->model.probability*dt,
                       //stardata->model.time,
                       stardata->star[k].mass,
                       stardata->star[k].luminosity,
                       stardata->star[k].radius,
                       // binary information: M2,p(years),a,ecc,v_periastron (i.e. MAX radial velocity, km/s)
                       stardata->star[Other_star(k)].mass,
                       stardata->common.orbital_period,
                       stardata->common.separation,
                       stardata->common.eccentricity,
                       v*1e-5 
                    );
            }

#ifdef NUCSYN
            if(((stardata->star[k].Xenv[XC12]/12.0)/(stardata->star[k].Xenv[XO16]/16.0)>1.0)
               &&((stardata->star[k].stellar_type==GIANT_BRANCH)||
                  (stardata->star[k].stellar_type==EAGB)||
                  (stardata->star[k].stellar_type==TPAGB))
               &&(teff<3800)
                )
            {
                // giant C stars i.e. N stars (AGB-AGB binaries counted twice! oops)
                PRINTF("NSTAR %g %d\n",
                       stardata->model.probability*dt,
                       stardata->star[k].stellar_type);
            }
#endif
            /*
             * Detect K+G giants (numbers of Jaschek and Jaschek)
             */
            if((teff>=3800.0)&&
               (teff<=5850.0)&&
               ((stardata->star[k].stellar_type==GIANT_BRANCH)||
                (stardata->star[k].stellar_type==EAGB)||
                (stardata->star[k].stellar_type==TPAGB))
                )
            {
                GKgiant=TRUE;
            }
        } //k loop

        // check if one star or the other is a G/K giant
        if(GKgiant==TRUE)
        {
            PRINTF("GKgiant %g\n",stardata->model.probability*dt);
        }

        /* Choose how to select clump stars */

#define CLUMP_CHeB
//#define CLUMP_by_luminosity

#ifdef CLUMP_CHeB
        /* clump stars are any CHeB stars */
        if((stardata->star[0].stellar_type==CHeB)||
           (stardata->star[1].stellar_type==CHeB))
#endif

#ifdef CLUMP_by_luminosity
            /* choose clump stars by luminosity */
            double lclump_min=pow(10.0,1.8);
        double lclump_max=pow(10.0,2.8);

        if(((stardata->star[0].stellar_type<=CHeB)
            &&(stardata->star[0].luminosity>lclump_min)&&(stardata->star[0].luminosity<lclump_max)
               )
           ||
           ((stardata->star[1].stellar_type<=CHeB)
            &&(stardata->star[1].luminosity>lclump_min)&&(stardata->star[1].luminosity<lclump_max)
               )
            )
#endif
        {
            PRINTF("CHEB %g\n",stardata->model.probability*dt);
        }
    }
#endif

#ifdef NS_BH_AIC_LOG
    STARLOOP(k)
    {
        SETstar(k);
        star->prev_luminosity=star->luminosity;
        star->prev_radius=star->radius;
        star->prev_mass=star->mass;
        star->prev_core_mass=star->core_mass;
        star->prev_stellar_type=star->stellar_type;
    }
#endif

    //#define COMENV_LECTURE1
#ifdef COMENV_LECTURE1
    /* for comenv lecture */

    if(first_timestep==TRUE)
    {
        stardata->common.star1_agb=FALSE;
        /*      PRINTF("BINSTART %g %g %g\n",
                stardata->model.probability,
                stardata->common.separation,
                stardata->common.orbital_period);
        */
    }
    /*
      else if(FEQUAL(stardata->model.max_evolution_time,
      stardata->model.time))
      {
      PRINTF("BINSTOP %g %g %g\n",
      stardata->model.probability,
      stardata->common.separation,
      stardata->common.orbital_period);
      }
    */

    if(stardata->star[0].stellar_type==TPAGB)
    {
        if(stardata->common.star1_agb==FALSE)
        {
            /* start AGB */
            //PRINTF("Start AGB at t=%g\n",stardata->model.time);
            stardata->common.star1_agb=TRUE;
            stardata->common.maxrad=0.0;
            stardata->common.minrad=stardata->star[0].radius;
        }
        stardata->common.maxrad=MAX(stardata->common.maxrad,stardata->star[0].radius);
        //PRINTF("stardata->common.maxrad=%g\n",stardata->common.maxrad);
    }
    else if((stardata->common.star1_agb==TRUE)&&
            (stardata->star[0].stellar_type>=COWD))
    {
        /* post-AGB */
        PRINTF("TPAGBEND %g %g %g %g\n",
               stardata->star[0].pms_mass,
               stardata->star[0].mass,
               stardata->common.minrad,stardata->common.maxrad);
        stardata->common.star1_agb=FALSE;
    }

#endif

#ifdef COMENV_LECTURE2
    
#endif

#ifdef ADAPTIVE_RLOG
    /*
     * RLOF stuff for lectures
     */
    if(stardata->star[0].stellar_type<10)
    {
        PRINTF("RLOF %g %g %g %g %g\n",
               stardata->model.time,
               stardata->star[0].mass,
               stardata->star[1].mass,
               stardata->star[1].v_eq,
               stardata->star[1].v_crit_eq
             
            );
    }
#endif
   
    //#define TYL_CE
#ifdef TYL_CE
    /*
     * Logging for Tyl: 
     * output lambda ion (the common envelope binding energy factor)
     * as a function of various stellar properties during the TPAGB
     */
    if(stardata->star[0].stellar_type==TPAGB)
    {
        struct star_t * star = &(stardata->star[0]);
        double rzams=rzamsf(star->phase_start_mass,
                            stardata->common.main_sequence_parameters);
        double lambda=celamf(star->stellar_type,
                             star->mass,
                             star->luminosity,
                             star->radius,
                             rzams,
                             1.0, // assume convfrac = 1
                             stardata->preferences->lambda_ionisation,
                             stardata);
                                          
        PRINTF("TYLCE %g %g %g %g %g %g\n",
               stardata->model.time-star->model_time_first_pulse,
               star->mass,
               star->radius,
               star->num_thermal_pulses,
               lambda,
               star->luminosity
            );
    }
#endif

#ifdef DETAILED_COMPACT_LOG
    output_to_detailed_compact_logfile(stardata);
#endif

#ifdef ANTI_TZ
    /* model to debunk TZ object myth */
    {
        int k;
        STARLOOP(k)
        {
            SETstar(k);

            if(star->stellar_type == EAGB &&
               star->started_EAGB == FALSE)
            {
                star->started_EAGB = TRUE;
                PRINTF("ANTI_TZ EAGB %g %g\n",
                       star->mass,
                       stardata->model.probability);
            }
            else if(star->TZ_object==TRUE)
            {
                PRINTF("ANTI_TZ TZ %d %g %g %g\n",
                       star->TZ_channel,
                       star->TZ_mass,
                       star->TZ_NS_mass,
                       stardata->model.probability);
                star->TZ_object=FALSE;
            }
        }
    }
#endif // ANTI_TZ

#ifdef WIND_ENERGY_LOGGING
    /*
     * Calculate energy outflow in winds
     */
    //if(binary==TRUE)
    { 
        /* orbital velocity in cgs */
        double vorb = sqrt (GRAVITATIONAL_CONSTANT * 
                            (stardata->star[0].mass+stardata->star[1].mass)*M_SUN/
                            (stardata->common.separation*R_SUN));
         
        int k;
        STARLOOP(k)
        {
            SETstar(k);
            if(star->stellar_type < HeWD)
            {

                /* vwind and mdot, in cgs */
                double vwind = star->vwind;
                double mdot = fabs(star->derivative[DERIVATIVE_STELLAR_MASS_WIND_LOSS]) * M_SUN / YEAR_LENGTH_IN_SECONDS;
                double dt_seconds = dt *1e6 * YEAR_LENGTH_IN_SECONDS;
 
                /* output only if there is a wind */
                if(mdot>TINY && star->vwind>TINY)
                {
                    double Lwind = mdot * (POW2(vorb) + POW2(vwind));
                    star->Ewind += Lwind * dt_seconds;
                    star->Elum += star->luminosity * L_SUN * dt_seconds;
                    PRINTF("EWIND%d t=%g Myr, stellar_type=%d, mdot=%g (Msun/y) vorb=%g (km/s), vwind=%g (km/s), wind power=%g (erg s-1) enh.fac=%g Ewind=%g Elum=%g\n",
                           k,
                           stardata->model.time,
                           star->stellar_type,
                           mdot * YEAR_LENGTH_IN_SECONDS/M_SUN,
                           vorb*1e-5, 
                           vwind*1e-5,
                           Lwind,
                           (mdot>TINY && vwind>TINY) ? Lwind : 0,
                           star->Ewind,
                           star->Elum);
                }
            }  
        }
    }

#endif // WIND_ENERGY_LOGGING

#ifdef AMANDA_HELIUM_CORE_SHIFT_LOGGING

    {
        /* log the core mass when the stellar type changes */
        if(first_timestep) stardata->common.amanda_st = -1;
        if(stardata->common.amanda_st != stardata->star[0].stellar_type)
        {
            PRINTF("CORE %d %g\n",
                   stardata->star[0].stellar_type,
                   stardata->star[0].core_mass);
            stardata->common.amanda_st = stardata->star[0].stellar_type;
        }
    }

#endif // AMANDA_HELIUM_CORE_SHIFT_LOGGING


#ifdef NS_BH_AIC_LOG

    /* TODO : needs updating */

    int this_star=star->starnum;
    int other_star=Other_star(star->starnum);
    if((*stellar_type==BLACK_HOLE)&&(stellar_typein==NEUTRON_STAR))
    {
        //double r2=stardata->star[other_star].radius;

        double m1,r1,l1,m2,r2,l2;
        Stellar_type st1,st2;

        if(stardata->star[other_star].stellar_type == MASSLESS_REMNANT)
        {
            /* there is a problem when there is a merger and the other
               star is reported as a massless remnant when really it was 
               something else (e.g. BLACK_HOLE) : in this case use the prev_
               variables
            */
            r1=stardata->star[this_star].prev_radius;
            l1=stardata->star[this_star].prev_luminosity;
            m1=stardata->star[this_star].prev_mass;
            st1=stardata->star[this_star].prev_stellar_type;
            r2=stardata->star[other_star].prev_radius;
            l2=stardata->star[other_star].prev_luminosity;
            m2=stardata->star[other_star].prev_mass;
            st2=stardata->star[other_star].prev_stellar_type;
        }
        else
        {
            r1=stardata->star[this_star].radius;
            l1=stardata->star[this_star].luminosity;
            m1=stardata->star[this_star].mass;
            st1=stardata->star[this_star].stellar_type;           
            r2=stardata->star[other_star].radius;
            l2=stardata->star[other_star].luminosity;
            m2=stardata->star[other_star].mass;
            st2=stardata->star[other_star].stellar_type;
        }

        //double angle=2.0*r2/stardata->common.separation;
        //double solid_angle=4.0*PI*(PI*POW2(r2))/(4.0*PI*POW2(stardata->common.separation));
        PRINTF("AIC NS->BH sgl=%d comenv=%d p=%12.12e a=%12.12e e=%12.12e t=%12.12e star1 m=%12.12e r=%12.12e  L=%12.12e type=%d star2 m=%12.12e r=%12.12e L=%12.12e type=%d\n",
               stardata->model.sgl,
               stardata->model.comenv_type,
               stardata->model.probability,
               stardata->common.separation,
               stardata->common.eccentricity,
               stardata->model.time,
               m1,
               r1,
               l1,
               st1,
               m2,
               r2,
               l2,
               st2
            );

        /* we don't care about the subsequent evolution */
        stardata->model.max_evolution_time = stardata->model.time;

    }
#endif

#ifdef MASSERON
    {
        struct star_t *star = &(stardata->star[0]);
        if(star->stellar_type<=TPAGB &&
           IS_NOT_ZERO(stardata->model.time))
        {
            /* fix spin rate on the main sequence */
            if(star->stellar_type<HERTZSPRUNG_GAP)
            {
                star->v_eq = star->vrot0;
                star->omega=45.35*star->v_eq/star->radius;
                star->angular_momentum = star->omega*moment_of_inertia(star,star->radius);
            }

            char spectype[10];
            spectral_type_string(spectype,stardata,star);
            
#define MASSERON_COLUMNS stardata->model.time,  \
                star->stellar_type,             \
                star->mass,                     \
                star->v_eq,                     \
                star->luminosity,               \
                star->radius,                   \
                logg(star),                     \
                TEFFSTAR(star),                 \
                spectype
            
            PRINTF("MASSERON %g %d %g %g %g %g %g %g %s\n",
                   MASSERON_COLUMNS
                );
        }
    }
#endif //MASSERON


#if defined NUCSYN && defined CN_THICK_DISC
    {
        /*
         * CN logging in thick-disc giants
         */
        int k;

        /* 
         * Selection:
         * choose only giants 
         * with age 4-10Gyr i.e. single star mass 0.93-1.25Msun
         * if Z=Z_SMC=0.004
         */

#define THICK_DISC_STAR (\
            star->stellar_type>=GIANT_BRANCH &&         \
            star->stellar_type<HeMS &&                  \
            stardata->model.time > stardata->preferences->thick_disc_end_age && \
            stardata->model.time < stardata->preferences->thick_disc_start_age && \
            IN_RANGE(loggv[k],                                          \
                     stardata->preferences->thick_disc_logg_min,        \
                     stardata->preferences->thick_disc_logg_max)        \
            )

        Boolean out=FALSE;
        double loggv[NUMBER_OF_STARS+1];
        STARLOOP(k)
        {            
            SETstar(k);
            loggv[k] = logg(star);
            if(THICK_DISC_STAR) out = TRUE;
            /*
              printf("STAR %d, type %d : age = %g cf %g .. %g, logg = %g vs %g .. %g\n",
                    k,
                    star->stellar_type,
                    stardata->model.time,
                    stardata->preferences->thick_disc_end_age,
                    stardata->preferences->thick_disc_start_age,
                    loggv[k],
                    stardata->preferences->thick_disc_logg_min,
                    stardata->preferences->thick_disc_logg_max 
                 );
            */
            if(star->blue_straggler) star->was_blue_straggler=TRUE;
        }

        if(out==TRUE)
        {
            PRINTF("CND %g %g %g %g %g %g %g ",
                   stardata->model.time,
                   dtp,

                   /* initial conditions */
                   stardata->common.zams_period,
                   stardata->common.zams_separation,
                   stardata->common.zams_eccentricity,
                   stardata->star[0].pms_mass,
                   stardata->star[1].pms_mass);

            if(SYSTEM_IS_BINARY)
            {
                /* current binary conditions */
                PRINTF("%g %g %g %g %g ",
                       stardata->common.orbital_period,
                       stardata->common.separation,
                       stardata->common.eccentricity,
                       stardata->star[0].mass,
                       stardata->star[1].mass
                    );
            }
            else
            {
                /* currently single : only give mass */
                PRINTF("* * * %g * ",
                       stardata->star[(stardata->star[0].stellar_type!=MASSLESS_REMNANT ? 0 : 1)].mass
                    );
            }

            STARLOOP(k)
            {
                SETstar(k);

                /*
                 * choose only hydrogen-rich 'giants' which are 
                 * between 8 and 10 Gyr old
                 */
                if(THICK_DISC_STAR)
                {
                    double * X = nucsyn_observed_surface_abundances(star);
                    PRINTF("%d %g %g %g %g %g %d ",
                           star->stellar_type,
                           1e-5*radial_velocity_K(stardata,PI/2.0,k),   // convert from cm/s to km/s
                           star->mass,
                           nucsyn_elemental_square_bracket(
                               "C","N",
                               X,
                               stardata->common.Xsolar,
                               stardata),
                           star->core_mass,
                           loggv[k],
                           star->was_blue_straggler
                        );
                }
                else
                {
                    PRINTF("* * * * * * * ");
                }
            }
            PRINTF("\n");
        }
    }     
#endif //CN_THICK_DISC

    if(0)
    printf("HRDHRD %g %g %g %g %g %g\n",
           stardata->model.time,
           HRdiagram(0),
           HRdiagram(1),
           MAX(stardata->star[0].derivative[DERIVATIVE_STELLAR_MASS_RLOF_GAIN],
               stardata->star[1].derivative[DERIVATIVE_STELLAR_MASS_RLOF_GAIN]));
/*    
    if(stardata->star[0].radius > stardata->star[0].roche_radius)
    {
        struct star_t * star = &(stardata->star[0]);
        printf("DTT %g %g %g %g %g %g\n",
               stardata->model.time,
               stardata->common.separation,
               star->radius,
               star->roche_radius,
               star->mass,
               -star->derivative[DERIVATIVE_STELLAR_MASS_RLOF_LOSS]
            );
        // ~0.002/Myr
    }
*/

    if(0)
    {
        printf("POO # t a M1 R1/RL1 M2 R2/RL2 loss trans adot Jdotsys\n");
        printf("POO %20.12g %20.12g %20.12g %g %20.12g %g %g %g %g %g %g\n",
               stardata->model.time,
               stardata->common.separation,
               stardata->star[0].mass,
               stardata->star[0].radius/stardata->star[0].roche_radius,
               stardata->star[1].mass, // 5 
               stardata->star[1].radius/stardata->star[1].roche_radius,
               stardata->star[1].derivative[DERIVATIVE_STELLAR_MASS_RLOF_LOSS], // 7 
               stardata->star[1].derivative[DERIVATIVE_STELLAR_MASS_RLOF_TRANSFER],
               stardata->star[0].derivative[DERIVATIVE_STELLAR_MASS_RLOF_GAIN],
               stardata->model.derivative[DERIVATIVE_ORBIT_SEMI_MAJOR_AXIS],
               stardata->model.derivative[DERIVATIVE_ORBIT_ANGMOM]
            );
    }    

#ifdef BUFFERED_STACK
    /*
     * Clear the binary_c printf buffer here, if required
     */
    if(stardata->preferences->internal_buffering ==
       INTERNAL_BUFFERING_PRINT)
    {
        Clear_printf_buffer;
    }
#endif//BUFFERED_STACK


}

