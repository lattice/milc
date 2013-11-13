#include "generic_includes.h"

#include <quda.h>
#include <quda_milc_interface.h>
#include "../include/generic_quda.h"


void update_u_gpu(Real eps){

  int i,dir;
  site *s;
  int j;

  int dim[4] = {nx, ny, nz, nt};

  initialize_quda();

  printf("EPS: %lf\n", eps);

  Real *momentum = (Real*)malloc(sites_on_node*4*sizeof(anti_hermitmat));
  Real *gauge = (Real*)malloc(sites_on_node*4*sizeof(su3_matrix));
  // Populate gauge and momentum fields
  FORALLSITES(i,s){
    for(dir=XUP; dir<=TUP; ++dir){
      for(j=0; j<18; ++j){
        gauge[(4*i + dir)*18 + j] = ((Real*)(&(s->link[dir])))[j];
      }
      for(j=0; j<10; ++j){
        momentum[(4*i + dir)*10 + j] = ((Real*)(&(s->mom[dir])))[j];
      }
    } // dir
  }

  

//  qudaUpdateU(PRECISION, eps, momentum, gauge); 

  // Copy updated gauge field back to site structure
  FORALLSITES(i,s){
    for(dir=XUP; dir<=TUP; ++dir){
      for(j=0; j<18; ++j){
       ((Real*)(&(s->link[dir])))[j] = gauge[(4*i + dir)*18 + j];
      }
    }
  }
  free(momentum);
  free(gauge);

  return;
}
