/******* dslash_fn_u1.c - dslash for improved KS fermions with U(1) gauge field  ****/
/* MIMD version 7 */
/* Kogut-Susskind fermions -- improved.  This version for "fat plus
   Naik" quark action.  Connection to nearest neighbors stored in
   fatlink and to third nearest neighbors in longlink */

/* This version overlaps computation and gathers from negative
   directions, and has an extra lattice loop devoted to exclusively to
   sub_four_vectors (traditional algorithm) */

/* Jim Hetrick, Kari Rummukainen, Doug Toussaint, Steven Gottlieb */
/* C. DeTar 9/29/01 Standardized prefetching and synced versions */

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#define LOOPEND
#include "../include/loopend.h"
#define FETCH_UP 1

#define INDEX_3RD(dir) (dir - 8)      /* this gives the 'normal' direction */


#define FORSOMEPARITYNOTZ(i,s,choice) \
{ register int loopend;  \
loopend= (choice)==EVEN ? even_sites_on_node : sites_on_node ; \
for( i=((choice)==ODD ? even_sites_on_node : 0 ), s= &(lattice[i]); \
i<loopend; i++,s++)if(s->z==0)
#define END_LOOP }

#define FORALLUP_3_DIR(dir) for(dir=X3UP; dir<=T3UP; dir++)

#define FORALLUP_3_DIRBUT(direction,dir) \
   FORALLUP_3_DIR(dir)if(dir != direction)

#define FORALLMYUP_3_DIR(dir) FORALLUP_3_DIRBUT(ZUP, dir)
//#define FORALLMYUP_3_DIR(dir) FORALLUP_3_DIR(dir)

#define FORALLMYUPDIR(dir) FORALLUPDIRBUT(ZUP, dir)
//#define FORALLMYUPDIR(dir) FORALLUPDIR(dir)

#define FORMYSITESANDPARITY(i,s,choice) FORSOMEPARITYNOTZ(i,s,choice)
//#define FORMYSITESANDPARITY(i,s,choice) FORSOMEPARITY(i,s,choice)

/* Temporary work space for dslash_fn_field_special */ 
static su3_vector *temp[9] ;
/* Flag indicating if temp is allocated               */
static int temp_not_allocated=1 ;

static void 
cleanup_one_gather_set(msg_tag *tags[])
{
  int i;

  FORALLMYUPDIR(i){
    cleanup_gather( tags[i] );
    cleanup_gather( tags[OPP_DIR(i)] );
  }

  FORALLMYUP_3_DIR(i){
    cleanup_gather( tags[i] );
    cleanup_gather( tags[OPP_3_DIR(i)] );
  }
}

void cleanup_gathers(msg_tag *tags1[], msg_tag *tags2[])
{
  cleanup_one_gather_set(tags1);
  cleanup_one_gather_set(tags2);
}

void cleanup_dslash_temps(){
  register int i ;
  if(!temp_not_allocated)
    for(i=0;i<9;i++) {
      free(temp[i]) ; 
    }
  temp_not_allocated=1 ;
}

static void allocate_dslash_temps(){
  FORALLMYUPDIR(dir){ 
    temp[dir]  =(complex *)malloc(sites_on_node*sizeof(complex));
    temp[dir+4]=(complex *)malloc(sites_on_node*sizeof(complex));
  }
  temp[8]=(complex *)malloc(sites_on_node*sizeof(complex));
  temp_not_allocated = 0 ;
}


/* D_slash routine - sets dest. on each site equal to sum of
   sources parallel transported to site, with minus sign for transport
   from negative directions.  Use "fatlinks" for one link transport,
   "longlinks" for three link transport. */

void dslash_fn_u1_field( complex *src, complex *dest, int parity,
			 ferm_links_t *fn) {
  msg_tag *tag[16];
    
   dslash_fn_field_special(src, dest, parity, tag, 1, fn);
   cleanup_one_gather_set(tag);
}

/* Special dslash for use by congrad.  Uses restart_gather_field() when
  possible. Next to last argument is an array of message tags, to be set
  if this is the first use, otherwise reused. If start=1,use
  start_gather_field, otherwise use restart_gather_field. 
  The calling program must clean up the gathers and temps! */
void dslash_fn_u1_field_special(complex *src, complex *dest,
				int parity, msg_tag **tag, int start,
				ferm_links_t *fn){
  register int i;
  register site *s;
  register int dir,otherparity=0;
  register su3_matrix *fat4, *long4;
  su3_matrix *t_fatlink;
  su3_matrix *t_longlink;
  
  /* allocate temporary work space only if not already allocated */
  if(temp_not_allocated)
    allocate_dslash_temps();
  
  /* load fatlinks and longlinks */
  if(!fn->fl.valid){
    printf("dslash_fn_field_special: invalid fn links!\n");
    terminate(1);
  }
  t_longlink = fn->fl.lng;
  t_fatlink = fn->fl.fat;

  switch(parity)
    {
    case EVEN:	otherparity=ODD; break;
    case ODD:	otherparity=EVEN; break;
    case EVENANDODD:	otherparity=EVENANDODD; break;
    }
  
  /* Start gathers from positive directions */
  /* And start the 3-step gather too */
  FORALLMYUPDIR(dir){ 
    if(start==1)
      {
	tag[dir] = start_gather_field( src, sizeof(complex), 
					   dir, parity,gen_pt[dir] );
	tag[DIR3(dir)] = start_gather_field(src, sizeof(complex),
						DIR3(dir),parity, 
						gen_pt[DIR3(dir)] );
      }
    else
      {
	restart_gather_field( src, sizeof(complex), 
				  dir, parity,gen_pt[dir], tag[dir]);
	restart_gather_field(src, sizeof(complex), DIR3(dir), parity, 
				 gen_pt[DIR3(dir)], tag[DIR3(dir)]);
      }
  }
  
  /* Multiply by adjoint matrix at other sites */
  /* Use fat link for single link transport */
  FORMYSITESANDPARITY( i, s, otherparity ){
    fat4 = &(t_fatlink[4*i]);
    long4 = &(t_longlink[4*i]);

    CMULJ_(fat4[0], src[i], temp[0][i]);
    CMULJ_(fat4[1], src[i], temp[1][i]);
    CMULJ_(fat4[2], src[i], temp[2][i]);
    CMULJ_(fat4[3], src[i], temp[3][i]);

    /* multiply by 3-link matrices too */

    CMULJ_(long4[0], src[i], temp[4][i]);
    CMULJ_(long4[1], src[i], temp[5][i]);
    CMULJ_(long4[2], src[i], temp[6][i]);
    CMULJ_(long4[3], src[i], temp[7][i]);

  } END_LOOP
      
  /* Start gathers from negative directions */
  FORALLMYUPDIR(dir){ 
    if (start==1) tag[OPP_DIR(dir)] = start_gather_field( temp[dir],
	  sizeof(complex), OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)] );
    else restart_gather_field( temp[dir], sizeof(complex), 
	   OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)], tag[OPP_DIR(dir)] );
  }

  /* Start 3-neighbour gathers from negative directions */
  FORALLMYUP_3_DIR(dir){
    if (start==1) tag[OPP_3_DIR(dir)]=start_gather_field(
		 temp[INDEX_3RD(dir)+4], sizeof(complex), 
		 OPP_3_DIR( dir), parity, gen_pt[OPP_3_DIR(dir)] );
    else restart_gather_field(temp[INDEX_3RD(dir)+4], 
	      sizeof(complex), OPP_3_DIR( dir),parity, 
	      gen_pt[OPP_3_DIR(dir)], tag[OPP_3_DIR(dir)] );
  }

  /* Wait gathers from positive directions, multiply by matrix and
     accumulate */
  /* wait for the 3-neighbours from positive directions, multiply */
  FORALLMYUPDIR(dir){ 
    wait_gather(tag[dir]);
    wait_gather(tag[DIR3(dir)]);
  }
  
  FORMYSITESANDPARITY(i,s,parity){
    fat4 = &(t_fatlink[4*i]);
    long4 = &(t_longlink[4*i]);

    mult_su3_mat_vec_sum_4dir( fat4,
	       (complex *)gen_pt[XUP][i], (complex *)gen_pt[YUP][i],
	       (complex *)gen_pt[ZUP][i], (complex *)gen_pt[TUP][i],
	       &(dest[i]) );
    
    mult_su3_mat_vec_sum_4dir( long4,
	    (complex *)gen_pt[X3UP][i], (complex *)gen_pt[Y3UP][i],
	    (complex *)gen_pt[Z3UP][i], (complex *)gen_pt[T3UP][i],
	    &(temp[8][i]));
  } END_LOOP
   
  /* Wait gathers from negative directions, accumulate (negative) */
  /* and the same for the negative 3-rd neighbours */
  FORALLMYUPDIR(dir){ 
    wait_gather(tag[OPP_DIR(dir)]);
  }
  FORALLMYUP_3_DIR(dir){
    wait_gather(tag[OPP_3_DIR(dir)]);
  }
  
  FORMYSITESANDPARITY(i,s,parity){
    
    sub_four_su3_vecs( &(dest[i]),
		       (complex *)(gen_pt[XDOWN][i]),
		       (complex *)(gen_pt[YDOWN][i]),
		       (complex *)(gen_pt[ZDOWN][i]),
		       (complex *)(gen_pt[TDOWN][i]) );
    sub_four_su3_vecs( &(temp[8][i]), 
		       (complex *)(gen_pt[X3DOWN][i]),
		       (complex *)(gen_pt[Y3DOWN][i]),
		       (complex *)(gen_pt[Z3DOWN][i]),
		       (complex *)(gen_pt[T3DOWN][i]) );
    /* Now need to add these things together */
    add_su3_vector(&(dest[i]), &(temp[8][i]),&(dest[i]));
  } END_LOOP 
      
}


