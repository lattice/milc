#define _inline_sse_mult_su3_mat_hwvec(aa,bb,cc) \
{ \
__asm__ __volatile__ ("movlps %0, %%xmm7 \n\t" \
                      "movaps %%xmm7, %%xmm0 \n\t" \
                      "unpcklps %%xmm0, %%xmm0 \n\t" \
                      "movhlps %%xmm0, %%xmm3 \n\t" \
                      "movlhps %%xmm0, %%xmm0 \n\t" \
                      "movhps %1, %%xmm7 \n\t" \
                      "movaps %%xmm7, %%xmm1 \n\t" \
                      "unpckhps %%xmm1, %%xmm1 \n\t" \
                      "movhlps %%xmm1, %%xmm4 \n\t" \
                      "movlhps %%xmm1, %%xmm1 \n\t" \
                      "movlps %2, %%xmm10 \n\t" \
                      "movaps %%xmm10, %%xmm2 \n\t" \
                      "unpcklps %%xmm2, %%xmm2 \n\t" \
                      "movhlps %%xmm2, %%xmm5 \n\t" \
                      "movlhps %%xmm2, %%xmm2 \n\t" \
                      "movlps %3, %%xmm12 \n\t" \
                      "movhps %4, %%xmm12 \n\t" \
                      "mulps %%xmm12, %%xmm0 \n\t" \
                      "movlps %5, %%xmm13" \
                      : \
                      : \
                      "m" ((aa)->e[0][0]), \
                      "m" ((aa)->e[0][1]), \
                      "m" ((aa)->e[0][2]), \
                      "m" ((bb)->h[0].c[0]), \
                      "m" ((bb)->h[1].c[0]), \
                      "m" ((bb)->h[0].c[1])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm10", \
                      "%xmm12", \
                      "%xmm13", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("movhps %0, %%xmm13 \n\t" \
                      "mulps %%xmm13, %%xmm1 \n\t" \
                      "addps %%xmm1, %%xmm0 \n\t" \
                      "movlps %1, %%xmm14 \n\t" \
                      "movhps %2, %%xmm14 \n\t" \
                      "mulps %%xmm14, %%xmm2 \n\t" \
                      "addps %%xmm2, %%xmm0 \n\t" \
                      "movlhps %%xmm3, %%xmm3 \n\t" \
                      "movlhps %%xmm4, %%xmm4 \n\t" \
                      "movlhps %%xmm5, %%xmm5 \n\t" \
                      "mulps %%xmm12, %%xmm3 \n\t" \
                      "mulps %%xmm13, %%xmm4 \n\t" \
                      "mulps %%xmm14, %%xmm5 \n\t" \
                      "addps %%xmm4, %%xmm3 \n\t" \
                      "addps %%xmm5, %%xmm3 \n\t" \
                      "xorps %3, %%xmm3 \n\t" \
                      "shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      : \
                      : \
                      "m" ((bb)->h[1].c[1]), \
                      "m" ((bb)->h[0].c[2]), \
                      "m" ((bb)->h[1].c[2]), \
                      "m" (_sse_sgn24)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm12", \
                      "%xmm13", \
                      "%xmm14", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm0, %0 \n\t" \
                      "movhps %%xmm0, %1 \n\t" \
                      : \
                      "=m" ((cc)->h[0].c[0]), \
                      "=m" ((cc)->h[1].c[0])\
                      : \
                      : \
                      "%xmm0", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm8 \n\t" \
                      "movaps %%xmm8, %%xmm1 \n\t" \
                      "unpcklps %%xmm1, %%xmm1 \n\t" \
                      "movhlps %%xmm1, %%xmm4 \n\t" \
                      "movlhps %%xmm1, %%xmm1 \n\t" \
                      "movhps %1, %%xmm8 \n\t" \
                      "movaps %%xmm8, %%xmm2 \n\t" \
                      "unpckhps %%xmm2, %%xmm2 \n\t" \
                      "movhlps %%xmm2, %%xmm5 \n\t" \
                      "movlhps %%xmm2, %%xmm2 \n\t" \
                      "movhps %2, %%xmm10 \n\t" \
                      "movaps %%xmm10, %%xmm3 \n\t" \
                      "unpckhps %%xmm3, %%xmm3 \n\t" \
                      "movhlps %%xmm3, %%xmm6 \n\t" \
                      "movlhps %%xmm3, %%xmm3 \n\t" \
                      "mulps %%xmm12, %%xmm1 \n\t" \
                      "mulps %%xmm13, %%xmm2 \n\t" \
                      "mulps %%xmm14, %%xmm3 \n\t" \
                      "addps %%xmm2, %%xmm1 \n\t" \
                      "addps %%xmm3, %%xmm1 \n\t" \
                      "movlhps %%xmm4, %%xmm4 \n\t" \
                      "movlhps %%xmm5, %%xmm5 \n\t" \
                      "movlhps %%xmm6, %%xmm6 \n\t" \
                      "mulps %%xmm12, %%xmm4 \n\t" \
                      "mulps %%xmm13, %%xmm5 \n\t" \
                      "mulps %%xmm14, %%xmm6 \n\t" \
                      "addps %%xmm5, %%xmm4 \n\t" \
                      "addps %%xmm6, %%xmm4 \n\t" \
                      "xorps %3, %%xmm4 \n\t" \
                      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      : \
                      : \
                      "m" ((aa)->e[1][0]), \
                      "m" ((aa)->e[1][1]), \
                      "m" ((aa)->e[1][2]), \
                      "m" (_sse_sgn24)\
                      : \
                      "%xmm1", \
                      "%xmm10", \
                      "%xmm12", \
                      "%xmm13", \
                      "%xmm14", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm8", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm1, %0 \n\t" \
                      "movhps %%xmm1, %1 \n\t" \
                      : \
                      "=m" ((cc)->h[0].c[1]), \
                      "=m" ((cc)->h[1].c[1])\
                      : \
                      : \
                      "%xmm1", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm3 \n\t" \
                      "unpcklps %%xmm3, %%xmm3 \n\t" \
                      "movhlps %%xmm3, %%xmm6 \n\t" \
                      "movlhps %%xmm3, %%xmm3 \n\t" \
                      "movhps %1, %%xmm4 \n\t" \
                      "unpckhps %%xmm4, %%xmm4 \n\t" \
                      "movhlps %%xmm4, %%xmm7 \n\t" \
                      "movlhps %%xmm4, %%xmm4 \n\t" \
                      "movlps %2, %%xmm11 \n\t" \
                      "movaps %%xmm11, %%xmm5 \n\t" \
                      "unpcklps %%xmm5, %%xmm5 \n\t" \
                      "movhlps %%xmm5, %%xmm10 \n\t" \
                      "movlhps %%xmm5, %%xmm5 \n\t" \
                      "mulps %%xmm12, %%xmm3 \n\t" \
                      "mulps %%xmm13, %%xmm4 \n\t" \
                      "mulps %%xmm14, %%xmm5 \n\t" \
                      "addps %%xmm4, %%xmm3 \n\t" \
                      "addps %%xmm5, %%xmm3 \n\t" \
                      "movlhps %%xmm6, %%xmm6 \n\t" \
                      "movlhps %%xmm7, %%xmm7 \n\t" \
                      "movlhps %%xmm10, %%xmm10 \n\t" \
                      "mulps %%xmm6, %%xmm12 \n\t" \
                      "mulps %%xmm7, %%xmm13 \n\t" \
                      "mulps %%xmm10, %%xmm14 \n\t" \
                      "addps %%xmm13, %%xmm12 \n\t" \
                      "addps %%xmm14, %%xmm12 \n\t" \
                      "xorps %3, %%xmm12 \n\t" \
                      "shufps $0xb1, %%xmm12, %%xmm12 \n\t" \
                      "addps %%xmm12, %%xmm3 \n\t" \
                      : \
                      : \
                      "m" ((aa)->e[2][0]), \
                      "m" ((aa)->e[2][1]), \
                      "m" ((aa)->e[2][2]), \
                      "m" (_sse_sgn24)\
                      : \
                      "%xmm10", \
                      "%xmm11", \
                      "%xmm12", \
                      "%xmm13", \
                      "%xmm14", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm3, %0 \n\t" \
                      "movhps %%xmm3, %1 \n\t" \
                      : \
                      "=m" ((cc)->h[0].c[2]), \
                      "=m" ((cc)->h[1].c[2])\
                      : \
                      : \
                      "%xmm3", \
                      "memory"); \
}
