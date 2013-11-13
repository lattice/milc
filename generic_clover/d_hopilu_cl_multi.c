/******* d_hopilu_cl_multi.c - multi mass multisource hopping
	 expansion for heavy clover fermions ****/
/* MIMD version 7 */
/* Memory stingy version
   "r" overwrites src on even sites
   "p" overwrites src on odd sites
   3/29/00 EVENFIRST is the rule now. CD.
   */

/* Modifications:

   8/02/01 Uses dslash_w_site_special - CD
   4/25/98 Initial version based on UMH's d_cgilu_cl_lean.c CD
   */

/* Requires qic wilson vector temporaries wv1 and wv2 */

/* ------------------------------------------------------------
   The matrix to be inverted is written in block even-odd form as


   M = ( R_o     -K D_oe )
       ( -K D_eo  R_e    )

   where R_o and R_e are 1 + K Clov_c/u_0^3 i sigma_mu,nu F_mu,nu
   are the site-diagonal hermitian clover matrices on odd and even
   sites.

   The ILU decomposition results in M = L A U


   M = ( 1            0 ) ( R_o  0   ) ( 1  -K/R_o D_oe )
       ( -K D_eo/R_o  1 ) ( 0    M_e ) ( 0     1        )

   where

         M_e = R_e - K^2 D_eo/R_o D_oe

   acts only on even sites. We invert M_e using a hopping
   expansion as follows:

   M_e = R_e ( 1 - K^2/R_e D_eo/R_o D_oe ) = R_e ( 1 - H )

   so

   M_e^(-1) = ( 1 + H + H^2 + ... ) R_e^(-1)

 ------------------------------------------------------------ */

#include "generic_clover_includes.h"

/*#define CGTIME*/          /* Uncomment if you want timing info */


int hopilu_cl_multi(       /* Return value is number of iterations taken */
    Real *kappas[],        /* kappas[isrc][ikappa] */
    int nkappa[],          /* nkappa[isrc] */
    wilson_vector* src[]   /* src[isrc][siteindex] src OVERWRITTEN! */
    wilson_vector** dest[] /* dest[isrc][ikappa][siteindex] ans and guess */
    quark_invert_control **qic[], /* *qic[isrc][ikappa] controls inversion */
    void *dmp              /* parameters defining the Dirac matrix */
    int nsrc,              /* Number of sources */
    )
{
  /* Unpack required members of the structures */
  int MinHOP = qic->min;      /* minimum number of iterations */
  int MaxHOP = qic->max;      /* maximum number of iterations */
  Real RsdCG = qic->resid * qic->resid;      /* desired residual - 
				normalized as (r*r)/(src_e*src_e */
  Real RRsdCG = qic->relresid * qic->relresid; /* desired relative residual - */
  int flag = qic->start_flag;   /* ignored - no guessing with hopping expansion! */
  /* End of unpacking required members of structures */

  wilson_vector *tmp, *mp;
  
  dirac_clover_param *dcp 
    = (dirac_clover_param *)dmp; /* Cast pass-through pointer */
  
  Real Clov_c = dcp->Clov_c;   /* Perturbative clover coeff */
  Real U0 = dcp->U0;           /* Tadpole correction to Clov_c */
  
  int N_iter;
  int ikappa,isrc;
  register int i;
  register site *s;
  Real *size_src;
  double sr;
  Real **Ksq;
  su3_matrix *r;
  Real CKU0 = Kappa*Clov_c/(U0*U0*U0);
  double dtime;
  msg_tag *tage[8],*tago[8];
  int is_startedo, is_startede;
  char myname[] = "hopilu_cl_multi";

  Ksq = (Real **)malloc(sizeof(Real*)*nsrc);
  if(Ksq == NULL){
    printf("%s(%d): Can't malloc Ksq\n",myname,this_node);
    terminate(1);
  }
  for(isrc = 0; isrc < nsrc; isrc++){
    Ksq[isrc] = (Real *)malloc(sizeof(Real)*nkappa[isrc]);
    if(Ksq[isrc] == NULL){
      printf("%s(%d): Can't malloc Ksq[%d]\n",isrc,myname,this_node);
      terminate(1);
    }
  }

  for(isrc = 0; isrc < nsrc; isrc++){
    for(ikappa = 0; ikappa < nkappa[isrc]; ikappa++){
      Ksq[isrc][ikappa] = kappas[isrc][ikappa]*kappas[isrc][ikappa];  

  tmp = (su3_matrix *)malloc(sizeof(su3_matrix)*sites_on_node);
  if(tmp == NULL){
    printf("%s(%d): Can't malloc tmp\n",myname,this_node);
    terminate(1);
  }

  mp = (su3_matrix *)malloc(sizeof(su3_matrix)*sites_on_node);
  if(mp == NULL){
    printf("%s(%d): Can't malloc mp\n",myname,this_node);
    terminate(1);
  }


  is_startedo = is_startede = 0;
  
  if(even_sites_on_node!=odd_sites_on_node){
    printf("Need same number of even and odd sites on each node\n");
    terminate(1);
  }
  
  /* Compute R_e and R_o and put in "clov" and "clov_diag" */
  compute_clov(gen_clov, CKU0);
  
  /* Invert BOTH R_o and R_e in place in "clov" and "clov_diag" */
  compute_clovinv(gen_clov, EVENANDODD);
  
  /* "r" is the residual vector, which is a pointer to src since the
     source is overwritten to save space. */
  r = src;
  
  /** if(this_node==0)printf("HOPILU: p=%d\n",p); **/
  
  /* HOP_ILU: */
  
  /* Start Inversion */
  
#ifdef CGTIME
  dtime = -dclock();
#endif

  /* ---------  src = L^(-1)*src  ------------- */

  /* ( src_o ) = ( 1           0 ) ( src_o )
     ( src_e )   ( K D_eo/R_o  1 ) ( src_e )

     */
  
  /* src_e = srce_e + K D_eo/R_o srce_o */
  /* (leaving src_o = src_o)   */

  /* mp_o = 1/R_o srce_o */
  mult_this_ldu_field(gen_clov, src, mp, ODD);
  /* mp_e = D_eo/R_o srce_o */
  dslash_w_field_special(mp, mp, PLUS, EVEN, tage, is_startede);
  is_startede = 1;

  /* src_e = srce_e + K D_eo/R_o srce_o */
  FOREVENSITES(i,s)
    scalar_mult_add_wvec( (wilson_vector *)F_PT(s,src), 
			  mp+i, Kappa, 
			  (wilson_vector *)F_PT(s,src) );

  /* r and p overwrite the source (bad for
     dynamical fermions but good for quenched calcs) */

  /* --------- dest_o = src_o ---------- */

  /* set dest = src on odd sites, whatever else you do
     (even if we restart with a nonzero solution vector, the end of the
     subroutine rebuilds the odd component from the even one. The (trivial)
     solution of the odd component of the equation is dest = src, before
     we rotate back to the basis in which  M is not checkerboard-diagonal) */
  FORODDSITES(i,s) {
    copy_wvec( (wilson_vector *)F_PT(s,src),
	       (wilson_vector *)F_PT(s,dest) );
  }
  
  /* --------- dest_e <- R_e^(-1) srce_e ---------- */
  mult_this_ldu_site(gen_clov, src, dest, EVEN);
  
  /* ------------ r_e <- dest_e -------------- */
  FOREVENSITES(i,s) {
    copy_wvec( (wilson_vector *)F_PT(s,dest),
	       (wilson_vector *)F_PT(s,r) );
  }

  sr=0.0;
  FOREVENSITES(i,s) {
    sr += (double)magsq_wvec( (wilson_vector *)F_PT(s,r) );
  }

  /* --------- size_src = sqrt(r_e*r_e) ---------- */
  g_doublesum(&sr);
  size_src = (Real)sqrt(sr);
  
  /* node0_printf("beginning inversion--size_src=%e\n",
     (double)size_src); */

  /* --------- Beginning of hop iterations --------- */

  for( N_iter = 1; N_iter < MaxHOP; N_iter = N_iter + 1) {
    
    /* ---- r_e = H r_e =  K^2/R_e D_eo/R_o D_oe r_e ----- */
    /* mp_o = D_oe r_e */
    dslash_w_field_special(r, mp, PLUS, ODD, tago, is_startedo);
    is_startedo = 1;
    /* tmp_o = 1/R_o D_oe r_e */
    mult_this_ldu_field(gen_clov, mp, tmp, ODD);
    /* mp_e = D_eo/R_o D_oe r_e */
    dslash_w_field_special(tmp, mp, PLUS, EVEN, tage, is_startede);
    is_startede = 1;
    /* r_e = 1/R_e D_eo/R_o D_oe r_e */
    mult_this_ldu_field(gen_clov, mp, r, EVEN);
    /* r_e = K^2/R_e D_eo/R_o D_oe r_e */
    sr=0.0;
    FOREVENSITES(i,s){
      scalar_mult_wvec((wilson_vector *)F_PT(s,r), Ksq, 
		       (wilson_vector *)F_PT(s,r));
      sr += (double)magsq_wvec( (wilson_vector *)F_PT(s,r) );
    }
    g_doublesum(&sr);
    qic->size_r = (Real)sqrt(sr)/size_src;
    
    /* node0_printf("N_iter = %d size_r = %g\n",N_iter,qic->size_r); */
    
    /* --------  dest_e <- dest_e + r_e -------*/
    FOREVENSITES(i,s){
      add_wilson_vector((wilson_vector *)F_PT(s,dest),
	       (wilson_vector *)F_PT(s,r), 
	       (wilson_vector *)F_PT(s,dest));
    }
    
    /* Stopping criterion */
    if(( N_iter > MinHOP ) && 
       ( RsdHOP > qic->size_r)) break;
  }
  /* ------------------------------------ */
  /* --------- End of iterations --------- */

  qic->final_iters = N_iter;
  qic->final_restart = 1;

#ifdef CGTIME
  dtime += dclock();
#endif
  if(this_node==0){
    if(N_iter==0)
      printf("HOPILU: NO ITERATIONS TAKEN size_r= %.2e\n",qic->size_r);
#ifdef CGTIME
    else
      printf("HOPILU: time = %.2e size_r= %.2e iters= %d MF = %.1f\n",
	     dtime,qic->size_r,N_iter,
	     (double)3744*N_iter*even_sites_on_node/(dtime*1e6));
#endif
    fflush(stdout);
  }
  
  /** if( (qic->size_r) > RsdHOP ) {
    if(this_node==0)printf(" HOP_ILU: Not Converged\n");
    } **/
  
  /* --------- dest_o =  dest_o + K/R_o D_oe dest_e --------- */
  /* mp_o = D_oe * dest_e */
  dslash_w_field_special(dest, mp, PLUS, ODD, tago, is_startedo);

  /* mp_o = dest_o + K D_oe * dest_e (remember dest_o = original src_o still)*/
  FORODDSITES(i,s) {
    scalar_mult_add_wvec( (wilson_vector *)F_PT(s,dest), 
			  mp+i, Kappa, mp+i );
  }
  /* dest_o = 1/R_o dest_o + K/R_o D_oe * dest_e */
  mult_this_ldu_field(gen_clov, mp, dest, ODD);
  
  for( i=XUP; i <= TUP; i++) {
    if(is_startedo)cleanup_gather(tago[i]);
    if(is_startedo)cleanup_gather(tago[OPP_DIR(i)]);
    if(is_startede)cleanup_gather(tage[i]);
    if(is_startede)cleanup_gather(tage[OPP_DIR(i)]);
  }
  is_startede = is_startedo = 0;

  /* Clean up */
  free(mp);
  free(tmp);
  for(isrc = 0; isrc < nsrc; isrc++)
    free(Ksq[isrc]);
  free(Ksq);

  cleanup_dslash_wtemps();
  return(N_iter);

} /* d_hopilu_cl_lean2.c */

