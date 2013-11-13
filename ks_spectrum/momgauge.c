/* ************************************************************	*/
/*								*/
/* 			       MOMGAUGE.C	   		*/
/*								*/
/* Generate U(1) fields A_\mu(p) in momentum space: 		*/
/*  S= 1/2 \sum(k) [A0(k)^2 \vec{k'}^2 + \vec{A}^2 k'^2]	*/
/*     k'^2 = \sum [4 (sin k/2)^2] 				*/
/*  \eta(k)_\mu = gaussrand()					*/
/*  A0(k) = [\vec{k'}^2/2]^{-1/2} * \eta{k}			*/
/*  \vec{A}(k) = [k'^2/2]^{-1/2} * \eta{k}			*/
/*								*/
/* Last updated on 08.08.07					*/
/*								*/
/* ************************************************************	*/

#include "include_u1g.h"

void momgauge(void)
{

  complex Am,Ap;
  Real lk[4],mm,mp,r1,r2;
  Real lkssq,lktsq,lksq;
  register site *s;
  register int i,dir;

  /* initialize */
  FORALLSITES(i,s){
      FORALLUPDIR(dir){
	  s->u1gf[dir]=cmplx(0.0,0.0);
	  s->u1link[dir]=cmplx(1.0,0.0);
	 }
    }

  /* A0 */
  FORALLSITES(i,s){
      /* lattice mom: \vec{k}^2, k1,k2,k3 and squares */
      lkssq=0.0;
      for(dir=XUP;dir<=ZUP;dir++){
          lk[dir]=2.0*(Real)sin((double)(s->mom[dir])/2.0);
	  lkssq+=sqr(lk[dir]);
	 }
      if(lkssq==0.0)				/* \vec{k}==0 */
	{
	if(latin[i]!=junk_id)
	  {
	  s->u1gf[TUP]=cmplx(0.0,0.0);
	  lattice[latin[i]].u1gf[TUP]=conjg(&(s->u1gf[TUP]));
	  }
	}
      else					/* \vec{k}!=0 */
	{
	if(latin[i]!=junk_id)
	  {
          s->u1gf[TUP].real=
		(Real)(gaussian_rand_no(&(s->site_prn))/sqrt(s->hp*lkssq/2.0));
          s->u1gf[TUP].imag=
		(Real)(gaussian_rand_no(&(s->site_prn))/sqrt(s->hp*lkssq/2.0));
	  lattice[latin[i]].u1gf[TUP]=conjg(&(s->u1gf[TUP]));
	  }
	}
     } /* FORALLSITES-ends */

  /* \vec{A} */
  FORALLSITES(i,s){
      /* lattice mom: k^2, k0,k1,k2,k3 and squares */
      lkssq=lktsq=lksq=0.0;
      for(dir=XUP;dir<=TUP;dir++){
	  lk[dir]=2.0*(Real)sin((double)(s->mom[dir])/2.0);
	  if(dir!=TUP) lkssq+=sqr(lk[dir]);
	  if(dir==TUP) lktsq+=sqr(lk[dir]);
	 }
      lksq=lkssq+lktsq;

      if(lkssq==0.0 && lktsq==0.0)		/* \vec{k}=0 & k0=0 */
	{
	if(latin[i]!=junk_id)
	  {
	  for(dir=XUP;dir<=ZUP;dir++){
	      s->u1gf[dir]=cmplx(0.0,0.0);
	      lattice[latin[i]].u1gf[dir]=conjg(&(s->u1gf[dir]));
	     }
	  }
	}
      else if(lkssq==0.0 && lktsq!=0.0)		/* \vec{k}=0 & k0!=0 */
	{
	if(latin[i]!=junk_id)
	  {
	  for(dir=XUP;dir<=ZUP;dir++){
	      s->u1gf[dir].real=
		(Real)(gaussian_rand_no(&(s->site_prn))/sqrt(s->hp*lksq/2.0));
	      s->u1gf[dir].imag=
		(Real)(gaussian_rand_no(&(s->site_prn))/sqrt(s->hp*lksq/2.0));
	      lattice[latin[i]].u1gf[dir]=conjg(&(s->u1gf[dir]));
	     }
	  }
	}
      else if(lkssq!=0.0) 			/* \vec{k}!=0 & k0=any */
	{
	if(lk[YUP]==0.0 && lk[ZUP]==0.0)	/* k1!=0 & k2=k3=0 */
	  {
	  if(latin[i]!=junk_id)
	    {
	    for(dir=YUP;dir<=ZUP;dir++){
                s->u1gf[dir].real=
                  (Real)(gaussian_rand_no(&(s->site_prn))/sqrt(s->hp*lksq/2.0));
                s->u1gf[dir].imag=
                  (Real)(gaussian_rand_no(&(s->site_prn))/sqrt(s->hp*lksq/2.0));
	        lattice[latin[i]].u1gf[dir]=conjg(&(s->u1gf[dir]));
               }
            s->u1gf[XUP].real=-(lk[YUP]*s->u1gf[YUP].real+
                                lk[ZUP]*s->u1gf[ZUP].real)/lk[XUP];
            s->u1gf[XUP].imag=-(lk[YUP]*s->u1gf[YUP].imag+
                                lk[ZUP]*s->u1gf[ZUP].imag)/lk[XUP];
	    lattice[latin[i]].u1gf[XUP]=conjg(&(s->u1gf[XUP]));
	    }
	  }
	else if(lk[ZUP]==0.0 && lk[XUP]==0.0)	/* k2!=0 & k3=k1=0 */
	  {
	  if(latin[i]!=junk_id)
	    {
	    for(dir=XUP;dir<=ZUP;dir++)if(dir!=YUP){
                s->u1gf[dir].real=
                  (Real)(gaussian_rand_no(&(s->site_prn))/sqrt(s->hp*lksq/2.0));
                s->u1gf[dir].imag=
                  (Real)(gaussian_rand_no(&(s->site_prn))/sqrt(s->hp*lksq/2.0));
	        lattice[latin[i]].u1gf[dir]=conjg(&(s->u1gf[dir]));
               }
            s->u1gf[YUP].real=-(lk[ZUP]*s->u1gf[ZUP].real+
                                lk[XUP]*s->u1gf[XUP].real)/lk[YUP];
            s->u1gf[YUP].imag=-(lk[ZUP]*s->u1gf[ZUP].imag+
                                lk[XUP]*s->u1gf[XUP].imag)/lk[YUP];
	    lattice[latin[i]].u1gf[YUP]=conjg(&(s->u1gf[YUP]));
	    }
	  }
	else if(lk[XUP]==0.0 && lk[YUP]==0.0)	/* k3!=0 & k1=k2=0 */
	  {
	  if(latin[i]!=junk_id)
	    {
	    for(dir=XUP;dir<=YUP;dir++){
                s->u1gf[dir].real=
                  (Real)(gaussian_rand_no(&(s->site_prn))/sqrt(s->hp*lksq/2.0));
                s->u1gf[dir].imag=
                  (Real)(gaussian_rand_no(&(s->site_prn))/sqrt(s->hp*lksq/2.0));
	        lattice[latin[i]].u1gf[dir]=conjg(&(s->u1gf[dir]));
               }
            s->u1gf[ZUP].real=-(lk[XUP]*s->u1gf[XUP].real+
                                lk[YUP]*s->u1gf[YUP].real)/lk[ZUP];
            s->u1gf[ZUP].imag=-(lk[XUP]*s->u1gf[XUP].imag+
                                lk[YUP]*s->u1gf[YUP].imag)/lk[ZUP];
	    lattice[latin[i]].u1gf[ZUP]=conjg(&(s->u1gf[ZUP]));
	    }
	  }
	else if(lk[XUP]!=0.0 && lk[YUP]!=0.0 && lk[ZUP]==0.0)
	  {					/* k1,k2!=0, k3=0 */
	  if(latin[i]!=junk_id)
	    {
	    mm=lksq;
            mp=lksq*(1+((sqr(lk[XUP])+sqr(lk[ZUP]))/sqr(lk[YUP])));
            r1=lk[XUP]/((Real)sqrt(sqr(lk[XUP])+sqr(lk[ZUP])));
            r2=lk[ZUP]/((Real)sqrt(sqr(lk[XUP])+sqr(lk[ZUP])));
            Am.real=(Real)(gaussian_rand_no(&(s->site_prn))/sqrt(s->hp*mm/2.0));
            Am.imag=(Real)(gaussian_rand_no(&(s->site_prn))/sqrt(s->hp*mm/2.0));
            Ap.real=(Real)(gaussian_rand_no(&(s->site_prn))/sqrt(s->hp*mp/2.0));
            Ap.imag=(Real)(gaussian_rand_no(&(s->site_prn))/sqrt(s->hp*mp/2.0));

            s->u1gf[XUP].real=r2*Am.real+r1*Ap.real;
            s->u1gf[XUP].imag=r2*Am.imag+r1*Ap.imag;
            s->u1gf[ZUP].real=-r1*Am.real+r2*Ap.real;
            s->u1gf[ZUP].imag=-r1*Am.imag+r2*Ap.imag;
            s->u1gf[YUP].real=-(lk[XUP]*s->u1gf[XUP].real+
                                lk[ZUP]*s->u1gf[ZUP].real)/lk[YUP];
            s->u1gf[YUP].imag=-(lk[XUP]*s->u1gf[XUP].imag+
                                lk[ZUP]*s->u1gf[ZUP].imag)/lk[YUP];
	    for(dir=XUP;dir<=ZUP;dir++){
		lattice[latin[i]].u1gf[dir]=conjg(&(s->u1gf[dir]));
	       }
	    }
	  }
	else					/* k3!=0 & k1 &/or k2!=0 */
	  {
	  if(latin[i]!=junk_id)
	    {
	    mm=lksq;
            mp=lksq*(1+((sqr(lk[XUP])+sqr(lk[YUP]))/sqr(lk[ZUP])));
            r1=lk[XUP]/((Real)sqrt(sqr(lk[XUP])+sqr(lk[YUP])));
            r2=lk[YUP]/((Real)sqrt(sqr(lk[XUP])+sqr(lk[YUP])));
            Am.real=(Real)(gaussian_rand_no(&(s->site_prn))/sqrt(s->hp*mm/2.0));
            Am.imag=(Real)(gaussian_rand_no(&(s->site_prn))/sqrt(s->hp*mm/2.0));
            Ap.real=(Real)(gaussian_rand_no(&(s->site_prn))/sqrt(s->hp*mp/2.0));
            Ap.imag=(Real)(gaussian_rand_no(&(s->site_prn))/sqrt(s->hp*mp/2.0));

            s->u1gf[XUP].real=r2*Am.real+r1*Ap.real;
            s->u1gf[XUP].imag=r2*Am.imag+r1*Ap.imag;
            s->u1gf[YUP].real=-r1*Am.real+r2*Ap.real;
            s->u1gf[YUP].imag=-r1*Am.imag+r2*Ap.imag;
            s->u1gf[ZUP].real=-(lk[XUP]*s->u1gf[XUP].real+
                                lk[YUP]*s->u1gf[YUP].real)/lk[ZUP];
            s->u1gf[ZUP].imag=-(lk[XUP]*s->u1gf[XUP].imag+
                                lk[YUP]*s->u1gf[YUP].imag)/lk[ZUP];
	    for(dir=XUP;dir<=ZUP;dir++){
		lattice[latin[i]].u1gf[dir]=conjg(&(s->u1gf[dir]));
	       }
	    }
	  }
	} /* lkssq!=0 ends */

     } /* FORALLSITES-ends */

  FORALLSITES(i,s)if(i==latin[i]){
      for(dir=XUP;dir<=TUP;dir++) s->u1gf[dir].imag=0.0;
     }

} /* end of momgauge() */

/* ************************************************************	*/

