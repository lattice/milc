/************************* control.c *******************************/
/* MIMD version 7 */
/* Main procedure for SU3 with dynamical staggered fermions        */
/* general quark action, general gauge action */

/* This version combines code for the PHI algorithm (approriate for 4
   flavors) and the R algorithm for "epsilon squared" updating of 
   1 to 4 flavors.  Compilation should occur with PHI_ALGORITHM defined
   for the former and not defined for the latter.  It also contains code
   for the hybrid Monte Carlo algorithm, for which HMC_ALGORITHM and
   PHI_ALGORITHM should be defined.  (Actually, the
   changes to control.c are minimal and the real differences will appear
   in update.c */

#define CONTROL
#include "ks_imp_includes.h"	/* definitions files and prototypes */
#include "lattice_qdp.h"
#ifdef HAVE_QUDA
#include <quda_milc_interface.h>
#endif

#ifdef HAVE_QIO
#include <qio.h>
#include "../include/io_scidac.h"
#endif

EXTERN gauge_header start_lat_hdr;	/* Input gauge field header */

int
main( int argc, char **argv )
{
  int meascount,traj_done;
  int prompt;
  int s_iters, avs_iters, avbcorr_iters;
  double dtime, dclock();
  
  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();
  /* set up */
  prompt = setup();

//  restore_random_state_scidac_to_site("randsave", F_OFFSET(site_prn));
//  restore_color_vector_scidac_to_site("xxx1save", F_OFFSET(xxx1),1);
//  restore_color_vector_scidac_to_site("xxx2save", F_OFFSET(xxx2),1);

  /* loop over input sets */
  while( readin(prompt) == 0) {
    
    /* perform warmup trajectories */
    dtime = -dclock();
    for( traj_done=0; traj_done < warms; traj_done++ ){
      update();
    }
    node0_printf("WARMUPS COMPLETED\n"); fflush(stdout);
    
    /* perform measuring trajectories, reunitarizing and measuring 	*/
    meascount=0;		/* number of measurements 		*/
    avs_iters = avbcorr_iters = 0;
    for( traj_done=0; traj_done < trajecs; traj_done++ ){ 
      
      /* do the trajectories */
      s_iters=update();
      
      /* measure every "propinterval" trajectories */
      if( (traj_done%propinterval)==(propinterval-1) ){
	
	/* call gauge_variable fermion_variable measuring routines */
	/* results are printed in output file */
	rephase(OFF);
	g_measure( );
	rephase(ON);

	restore_fermion_links_from_site(fn_links, prec_pbp);

	/* Measure pbp, etc */
#ifdef ONEMASS
	f_meas_imp( npbp_reps_in, prec_pbp, 
		    F_OFFSET(phi),F_OFFSET(xxx),mass, 0, fn_links);
#else
	f_meas_imp( npbp_reps_in, prec_pbp, 
		    F_OFFSET(phi1), F_OFFSET(xxx1), mass1, 0, fn_links);
	f_meas_imp( npbp_reps_in, prec_pbp, 
		    F_OFFSET(phi2), F_OFFSET(xxx2), mass2, 0, fn_links);
#endif

	/* Measure derivatives wrto chemical potential */
#ifdef D_CHEM_POT
#ifdef ONEMASS
	Deriv_O6( npbp_reps_in, prec_pbp, F_OFFSET(phi1), 
		  F_OFFSET(xxx1), F_OFFSET(xxx2), mass, fn_links);
#else
	Deriv_O6( npbp_reps_in, prec_pbp, F_OFFSET(phi1), 
		  F_OFFSET(xxx1), F_OFFSET(xxx2), mass1, fn_links);
	Deriv_O6( npbp_reps_in, prec_pbp, F_OFFSET(phi1), 
		  F_OFFSET(xxx1), F_OFFSET(xxx2), mass2, fn_links);
#endif
#endif

	avs_iters += s_iters;
	++meascount;
	fflush(stdout);
      }
    }	/* end loop over trajectories */
    
    node0_printf("RUNNING COMPLETED\n"); fflush(stdout);
    if(meascount>0)  {
      node0_printf("average cg iters for step= %e\n",
		   (double)avs_iters/meascount);
    }
    
    dtime += dclock();
    if(this_node==0){
      printf("Time = %e seconds\n",dtime);
      printf("total_iters = %d\n",total_iters);
#ifdef HISQ_SVD_COUNTER
      printf("hisq_svd_counter = %d\n",hisq_svd_counter);
#endif
#ifdef HYPISQ_SVD_COUNTER
      printf("hypisq_svd_counter = %d\n",hypisq_svd_counter);
#endif
      
#ifdef HISQ_FORCE_FILTER_COUNTER
      printf("hisq_force_filter_counter = %d\n",hisq_force_filter_counter);
#endif
#ifdef HYPISQ_FORCE_FILTER_COUNTER
      printf("hypisq_force_filter_counter = %d\n",hypisq_force_filter_counter);
#endif
    }
    fflush(stdout);
    
    /* save lattice if requested */
    if( saveflag != FORGET ){
      rephase( OFF );
      save_lattice( saveflag, savefile, stringLFN );
      rephase( ON );

    /* Destroy fermion links (created in readin() */

#if FERM_ACTION == HISQ
    destroy_fermion_links_hisq(fn_links);
#elif FERM_ACTION == HYPISQ
    destroy_fermion_links_hypisq(fn_links);
#else
    destroy_fermion_links(fn_links);
#endif
    fn_links = NULL;


#ifdef HAVE_QIO
//       save_random_state_scidac_from_site("randsave", "Dummy file XML",
//        "Random number state", QIO_SINGLEFILE, F_OFFSET(site_prn));
//       save_color_vector_scidac_from_site("xxx1save", "Dummy file XML",
//        "xxx vector", QIO_SINGLEFILE, F_OFFSET(xxx1),1);
//       save_color_vector_scidac_from_site("xxx2save", "Dummy file XML",
//        "xxx vector", QIO_SINGLEFILE, F_OFFSET(xxx2),1);
#endif

    }
  }

#ifdef HAVE_QUDA
  qudaFinalize();
#endif

  normal_exit(0);
  return 0;
}
