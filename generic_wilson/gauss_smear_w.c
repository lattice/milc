/************************** gauss_smear_w.c *****************************/
/* MIMD version 7 */
/*
 * Create a Gauss-smeared source, Chroma style
 * 
 * CD 4/07 Stolen from Chroma. MILC version.
 */

#include "generic_wilson_includes.h"
#include <string.h>

static wilson_vector *wtmp[8] ;

/*------------------------------------------------------------*/
static void 
malloc_kg_temps(){
  int dir;
  
  for(dir=0;dir<8;dir++)wtmp[dir] = NULL;

  FORALLUPDIRBUT(TUP,dir){
    wtmp[dir] =(wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector));
    if(wtmp[dir] == NULL){
      printf("node %d can't malloc wtmp[%d]\n",this_node,dir);
      terminate(1);
    }
    memset(wtmp[dir],'\0',sites_on_node*sizeof(wilson_vector));
    
    wtmp[OPP_DIR(dir)] =(wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector));
    if(wtmp[OPP_DIR(dir)] == NULL){
      printf("node %d can't malloc wtmp[%d]\n",this_node,OPP_DIR(dir));
      terminate(1);
    }
    memset(wtmp[OPP_DIR(dir)],'\0',sites_on_node*sizeof(wilson_vector));
  }
}

/*------------------------------------------------------------*/
static void 
cleanup_kg_temps(){
  int i ;
  for(i=0;i<8;i++){
    if(wtmp[i] != NULL){
      free(wtmp[i]); 
      wtmp[i] = NULL;
    }
  }
}

/*------------------------------------------------------------*/
/* Double forward parallel transport the quick and dirty way.
   Result in dest */

static void 
forward2(int dir, wilson_vector *dest, wilson_vector *src, 
	 su3_matrix *t_links, int t0)
{
  int i;
  site *s;
  msg_tag *tag;
  wilson_vector *tmp = create_wv_field();
  
  /* start parallel transport of src from up dir */
  tag = start_gather_field( src, sizeof(wilson_vector),
			    dir, EVENANDODD, gen_pt[dir] );
  wait_gather(tag);
  
  /* tmp <- U(up,dir) shift(up,dir)(src) */
  FORALLSITES(i,s){
    if(t0 == ALL_T_SLICES || s->t == t0){
      mult_mat_wilson_vec( t_links + 4*i + dir,  
			   (wilson_vector * )(gen_pt[dir][i]), 
			   tmp + i ); 
    }
  }

  cleanup_gather(tag);

  /* start parallel transport of tmp from up dir */
  tag = start_gather_field( tmp, sizeof(wilson_vector),
			    dir, EVENANDODD, gen_pt[dir] );
  wait_gather(tag);

  /* dest <- U(up,dir) shift(up,2dir)(src) */
  FORALLSITES(i,s){
    if(t0 == ALL_T_SLICES || s->t == t0){
      mult_mat_wilson_vec( t_links + 4*i + dir,  
			   (wilson_vector * )(gen_pt[dir][i]), 
			   dest + i ); 
    }
  }

  cleanup_gather(tag);
  destroy_wv_field(tmp);
}

/*------------------------------------------------------------*/
/* Double backward parallel transport the quick and dirty way.
   Result in dest */

static void 
backward2(int dir, wilson_vector *dest, wilson_vector *src, 
	 su3_matrix *t_links, int t0)
{
  int i;
  site *s;
  msg_tag *tag;
  wilson_vector *tmp = create_wv_field();

  /* prepare parallel transport of psi from down dir */

  /* dest <- U^dagger(down,dir) src */
  FORALLSITES(i,s){
    /* Work only on the specified time slice(s) */
    if(t0 == ALL_T_SLICES || s->t == t0){
      mult_adj_mat_wilson_vec( t_links +  4*i + dir, src + i, 
			       dest + i );
    }
  }
  
  /* gen_pt_array <- shift(down,dir)(dest) */
  tag = start_gather_field(dest, 
			   sizeof(wilson_vector), OPP_DIR(dir),
			   EVENANDODD, gen_pt[OPP_DIR(dir)] );
  wait_gather(tag);
  
  /* tmp <- U^dagger(down,dir) gen_pt_array */
  FORALLSITES(i,s){
    if(t0 == ALL_T_SLICES || s->t == t0){
      mult_adj_mat_wilson_vec( t_links + 4*i + dir,  
			       (wilson_vector * )(gen_pt[OPP_DIR(dir)][i]), 
			       tmp + i ); 
    }
  }
  
  cleanup_gather(tag);

  /* gen_pt_array <- shift(down,dir)(tmp) */
  tag = start_gather_field(tmp, 
			   sizeof(wilson_vector), OPP_DIR(dir),
			   EVENANDODD, gen_pt[OPP_DIR(dir)] );

  wait_gather(tag);

  /* dest <- gen_pt_array */
  FORALLSITES(i,s){
    if(t0 == ALL_T_SLICES || s->t == t0){
      dest[i] = *((wilson_vector *)gen_pt[OPP_DIR(dir)][i]);
    }
  }

  cleanup_gather(tag);
  
  destroy_wv_field(tmp);
}

/*------------------------------------------------------------*/
/* Compute chi <- msq * psi - Lapl_3d psi
   where Lapl_3d psi(r) = -6 psi + sum_{dir=1}^3 [psi(r+dir) + psi(r-dir)] */

static void 
klein_gord_wv_field_stride1(wilson_vector *psi, wilson_vector *chi, 
			    su3_matrix *t_links, Real msq, int t0)
{
  Real ftmp = 6 + msq;  /* for 3D */
  int i, dir;
  site *s;
  msg_tag *tag[8];

  malloc_kg_temps();

  /* chi = psi * ftmp; */
  FORALLSITES(i,s){
    if(t0 == ALL_T_SLICES || s->t == t0)
      scalar_mult_wvec(psi + i, ftmp, chi + i);
  }

  /* start parallel transport of psi from up dir */
  FORALLUPDIRBUT(TUP,dir){
    tag[dir]=start_gather_field( psi, sizeof(wilson_vector),
				 dir, EVENANDODD, gen_pt[dir] );
  }

  /* prepare parallel transport of psi from down dir */
  FORALLSITES(i,s){
    /* Work only on the specified time slice */
    if(t0 == ALL_T_SLICES || s->t == t0){
      FORALLUPDIRBUT(TUP,dir){
	mult_adj_mat_wilson_vec( t_links +  4*i + dir, psi + i, 
				 wtmp[OPP_DIR(dir)] + i );
      }
    }
  }
  
  FORALLUPDIRBUT(TUP,dir){
    tag[OPP_DIR(dir)] = 
      start_gather_field(wtmp[OPP_DIR(dir)], 
			 sizeof(wilson_vector), OPP_DIR(dir),
			 EVENANDODD, gen_pt[OPP_DIR(dir)] );
  }
  
  FORALLUPDIRBUT(TUP,dir){
    wait_gather(tag[dir]);
  }
  
  /* chi <- chi - sum_dir U(up,dir) shift(up,dir)(psi) */
  FORALLSITES(i,s){
    if(t0 == ALL_T_SLICES || s->t == t0){
      FORALLUPDIRBUT(TUP,dir){
	mult_mat_wilson_vec( t_links + 4*i + dir,  
			     (wilson_vector * )(gen_pt[dir][i]), 
			     wtmp[dir] + i ); 
      }
      sub_wilson_vector( chi + i, wtmp[XUP] + i, chi + i);
      sub_wilson_vector( chi + i, wtmp[YUP] + i, chi + i);
      sub_wilson_vector( chi + i, wtmp[ZUP] + i, chi + i);
    }
  }
  
  FORALLUPDIRBUT(TUP,dir){
    cleanup_gather(tag[dir]);
  }
  
  FORALLUPDIRBUT(TUP,dir){
    wait_gather(tag[OPP_DIR(dir)]);
  }
  
  /* chi <- chi - sum_dir U(down,dir) shift(down,dir)(psi) */
  FORALLSITES(i,s){
    if(t0 == ALL_T_SLICES || s->t == t0){
      sub_wilson_vector( chi + i,
			 (wilson_vector *)(gen_pt[XDOWN][i]), chi + i);
      sub_wilson_vector( chi + i,
			 (wilson_vector *)(gen_pt[YDOWN][i]), chi + i);
      sub_wilson_vector( chi + i,
			 (wilson_vector *)(gen_pt[ZDOWN][i]), chi + i);
    }
  }
  
  FORALLUPDIRBUT(TUP,dir){
    cleanup_gather(tag[OPP_DIR(dir)]);
  }
  cleanup_kg_temps();
}

/*------------------------------------------------------------*/
/* For staggered fermions we compute the Laplacian on sites displaced
   by 2 lattice units */
/* Compute chi <- msq * psi - Lapl_3d psi
   where Lapl_3d psi(r) = -6 psi + sum_{dir=1}^3 [psi(r+2*dir) + psi(r-2*dir)] */

static void 
klein_gord_wv_field_stride2(wilson_vector *psi, wilson_vector *chi, 
			    su3_matrix *t_links, Real msq, int t0)
{
  Real ftmp = 6 + msq;  /* for 3D */
  int i, dir;
  site *s;

  malloc_kg_temps();

  /* chi = psi * ftmp; */
  FORALLSITES(i,s){
    if(t0 == ALL_T_SLICES || s->t == t0)
      scalar_mult_wvec(psi + i, ftmp, chi + i);
  }

  /* do 2-link parallel transport of psi in all dirs */
  FORALLUPDIRBUT(TUP,dir){
    forward2(dir, wtmp[dir], psi, t_links, t0);
    backward2(dir, wtmp[OPP_DIR(dir)], psi, t_links, t0);
  }
  
  /* chi <- chi - sum_dir U(up,dir) shift(up,dir)(psi) -
     sum_dir shift(down,dir) U^\dagger(down,dir)(psi) */
  FORALLSITES(i,s){
    if(t0 == ALL_T_SLICES || s->t == t0){
      sub_wilson_vector( chi + i, wtmp[XUP] + i, chi + i);
      sub_wilson_vector( chi + i, wtmp[YUP] + i, chi + i);
      sub_wilson_vector( chi + i, wtmp[ZUP] + i, chi + i);
      sub_wilson_vector( chi + i, wtmp[XDOWN] + i, chi + i);
      sub_wilson_vector( chi + i, wtmp[YDOWN] + i, chi + i);
      sub_wilson_vector( chi + i, wtmp[ZDOWN] + i, chi + i);
    }
  }

  cleanup_kg_temps();
}

/*------------------------------------------------------------*/

/* Computes 
   src <- exp[-width^2/4 * Lapl_3d] src 
   by approximating exp(a) as (1 + a/iters)^iters 
   and Lap_3d is the discrete three dimensional Laplacian
*/

void gauss_smear_wv_field(wilson_vector *src, su3_matrix *t_links,
			  int stride, Real width, int iters, int t0)
{
  wilson_vector *tmp;
  Real ftmp = -(width*width)/(4*iters);
  Real ftmpinv = 1. / ftmp;
  int i, j;
  site *s;

  if(t_links == NULL){
    printf("gauss_smear_field(%d): NULL t_links\n",this_node);
    terminate(1);
  }

  tmp = (wilson_vector *)malloc(sizeof(wilson_vector)*sites_on_node);
  if(tmp == NULL){
    printf("gauss_smear_site(%d): No room for temporary source\n",this_node);
    terminate(1);
  }
  
  /* We want (1 + ftmp * Lapl ) = (Lapl + 1/ftmp)*ftmp */

  for(j = 0; j < iters; j++)
    {
      FORALLSITES(i,s){
	/* tmp = src * ftmp; */
	if(t0 == ALL_T_SLICES || s->t == t0)
	  scalar_mult_wvec(src+i, ftmp, tmp+i);
      }
      if(stride == 1)
	klein_gord_wv_field_stride1(tmp, src, t_links, ftmpinv, t0);
      else
	klein_gord_wv_field_stride2(tmp, src, t_links, ftmpinv, t0);
    }

  free(tmp);
}

/*------------------------------------------------------------*/

void gauss_smear_wv_site(field_offset src, su3_matrix *t_links, 
			 int stride, Real width, int iters, int t0)
{
  wilson_vector *srctmp;
  int i;
  site *s;
  
  /* Copy source to temporary field */
  srctmp = (wilson_vector *)malloc(sizeof(wilson_vector)*sites_on_node);
  if(srctmp == NULL){
    printf("gauss_smear_site(%d): No room for temporary source\n",this_node);
    terminate(1);
  }
  
  FORALLSITES(i,s){
    srctmp[i] = *((wilson_vector *)F_PT(s,src));
  }

  /* Smear in temporary field */
  gauss_smear_wv_field(srctmp, t_links, stride, width, iters, t0);
    

  FORALLSITES(i,s){
    *((wilson_vector *)F_PT(s,src)) = srctmp[i];
  }

  free(srctmp);
}
