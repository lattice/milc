#include <stdio.h>

/*------------------------------------------------------------------*/
/* Compute the hypercube coordinate relative to an offset.  We assume
   that all lattice dimensions are even, as they should be for
   staggered fermions! */
static short *
hyp_coord(int s[], int r0[]){
  static short h[4];
  h[0] = (s[0] - r0[0]) & 0x1;
  h[1] = (s[1] - r0[1]) & 0x1;
  h[2] = (s[2] - r0[2]) & 0x1;
  h[3] = (s[3] - r0[3]) & 0x1;
  return h;
}

/* Compute the parity of the site s relative to offset r0. */
/* Return 0x0 for even and 0x1 for odd */
static short 
hyp_parity_bit(int s[], int r0[]){
  short p;
  int s_parity = (s[0] + s[1] + s[2] + s[3]) % 2;
  if((r0[0] + r0[1] + r0[2] + r0[3]) % 2 == 0)
    p = 0;
  else p = 1;
  if(p == s_parity)
    return 0x0;
  else
    return 0x1;
}

/*------------------------------------------------------------------*/
static float
spin_sign(int spin, int r0[], int s[]){
  /* Compute (-)^[spin dot (x-r0)] epsilon(x-r0)^spin */
  /* Same as prod_\mu [eta_\mu(x-r0) zeta_\mu(x-r0)]^(s_\mu) */
  int j, mask;
  float sign = 1.;
  short *h = hyp_coord(s, r0);
  short hp = hyp_parity_bit(s, r0);
  
  /* For each nonzero gamma_mu bit in "spin",
     a factor of (-)^(x[mu]-r0[mu]) epsilon(x) */
  mask = 1;
  for(j = 0; j < 3; j++){
    if( (spin & mask) && (hp ^ h[j]) ) sign = -sign;
    mask <<= 1;
  }

  return sign;
}

int main(){
  int r0[4], s[4], spin;
  scanf("%d %d %d %d %d %d %d %d %d", &spin, &r0[0], &r0[1], &r0[2], &r0[3], &s[0], &s[1], &s[2], &s[3]);
  printf("%g\n", spin_sign(spin, r0, s));
  return 1;
}
