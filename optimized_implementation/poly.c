#include <stdint.h>
#include <stdlib.h>
#include "reduce.h"
#include "consts.h"
#include "rounding.h"
#include "symmetric.h"
#include <string.h>
#include <time.h>
#include "packing.h"
#include "sign.h"
#include <stdio.h>

extern void ntt1_asm(__m256i *a, const __m256i *qdata);
extern void ntt3_asm(__m256i *a, const __m256i *qdata);
extern void ntt5_asm(__m256i *a, const __m256i *qdata);
extern void ntt5plus_asm(__m256i *a, const __m256i *qdata);
extern void poly_base_mul1_asm(__m256i *c, const __m256i *a, const __m256i *b, const __m256i *qdata);
extern void poly_base_mul3_asm(__m256i *c, const __m256i *a, const __m256i *b, const __m256i *qdata);
extern void poly_base_mul5_asm(__m256i *c, const __m256i *a, const __m256i *b, const __m256i *qdata);
extern void poly_base_mul5plus_asm(__m256i *c, const __m256i *a, const __m256i *b, const __m256i *qdata);
extern void invntt1_asm(__m256i *a, const __m256i *qdata);
extern void invntt3_asm(__m256i *a, const __m256i *qdata);
extern void invntt5_asm(__m256i *a, const __m256i *qdata);
extern void invntt5plus_asm(__m256i *a, const __m256i *qdata);

void ntt_avx_asm(poly *Out, poly *A)
{
    if(Out!=A) memcpy(Out,A,sizeof(int32_t)*N);
    
#if NIMS_TRI_NTT_MODE == 1
    ntt1_asm(Out->vec, qdata.vec);
#elif NIMS_TRI_NTT_MODE == 3
    ntt3_asm(Out->vec, qdata.vec);
#elif NIMS_TRI_NTT_MODE == 5
    ntt5_asm(Out->vec, qdata.vec);    
#elif NIMS_TRI_NTT_MODE == 55
    ntt5plus_asm(Out->vec, qdata.vec);    
#endif
}

void poly_base_mul_avx_asm(poly *c, poly *a, poly *b)
{
#if NIMS_TRI_NTT_MODE == 1
    poly_base_mul1_asm(c->vec, a->vec, b->vec, qdata.vec);
#elif NIMS_TRI_NTT_MODE == 3
    poly_base_mul3_asm(c->vec, a->vec, b->vec, qdata.vec);
#elif NIMS_TRI_NTT_MODE == 5
    poly_base_mul5_asm(c->vec, a->vec, b->vec, qdata.vec);
#elif NIMS_TRI_NTT_MODE == 55
    poly_base_mul5plus_asm(c->vec, a->vec, b->vec, qdata.vec);
#endif
}

void invntt_avx_asm(poly *Out, poly *A)
{
    if(Out!=A) memcpy(Out,A,sizeof(int32_t)*N);
    
#if NIMS_TRI_NTT_MODE == 1
    invntt1_asm(Out->vec, qdata.vec);
#elif NIMS_TRI_NTT_MODE == 3
    invntt3_asm(Out->vec, qdata.vec);
#elif NIMS_TRI_NTT_MODE == 5
    invntt5_asm(Out->vec, qdata.vec);
#elif NIMS_TRI_NTT_MODE == 55
    invntt5plus_asm(Out->vec, qdata.vec);  
#endif
}
const uint8_t idxlut[256][8] = {
    {0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0},
    {1, 0, 0, 0, 0, 0, 0, 0},
    {0, 1, 0, 0, 0, 0, 0, 0},
    {2, 0, 0, 0, 0, 0, 0, 0},
    {0, 2, 0, 0, 0, 0, 0, 0},
    {1, 2, 0, 0, 0, 0, 0, 0},
    {0, 1, 2, 0, 0, 0, 0, 0},
    {3, 0, 0, 0, 0, 0, 0, 0},
    {0, 3, 0, 0, 0, 0, 0, 0},
    {1, 3, 0, 0, 0, 0, 0, 0},
    {0, 1, 3, 0, 0, 0, 0, 0},
    {2, 3, 0, 0, 0, 0, 0, 0},
    {0, 2, 3, 0, 0, 0, 0, 0},
    {1, 2, 3, 0, 0, 0, 0, 0},
    {0, 1, 2, 3, 0, 0, 0, 0},
    {4, 0, 0, 0, 0, 0, 0, 0},
    {0, 4, 0, 0, 0, 0, 0, 0},
    {1, 4, 0, 0, 0, 0, 0, 0},
    {0, 1, 4, 0, 0, 0, 0, 0},
    {2, 4, 0, 0, 0, 0, 0, 0},
    {0, 2, 4, 0, 0, 0, 0, 0},
    {1, 2, 4, 0, 0, 0, 0, 0},
    {0, 1, 2, 4, 0, 0, 0, 0},
    {3, 4, 0, 0, 0, 0, 0, 0},
    {0, 3, 4, 0, 0, 0, 0, 0},
    {1, 3, 4, 0, 0, 0, 0, 0},
    {0, 1, 3, 4, 0, 0, 0, 0},
    {2, 3, 4, 0, 0, 0, 0, 0},
    {0, 2, 3, 4, 0, 0, 0, 0},
    {1, 2, 3, 4, 0, 0, 0, 0},
    {0, 1, 2, 3, 4, 0, 0, 0},
    {5, 0, 0, 0, 0, 0, 0, 0},
    {0, 5, 0, 0, 0, 0, 0, 0},
    {1, 5, 0, 0, 0, 0, 0, 0},
    {0, 1, 5, 0, 0, 0, 0, 0},
    {2, 5, 0, 0, 0, 0, 0, 0},
    {0, 2, 5, 0, 0, 0, 0, 0},
    {1, 2, 5, 0, 0, 0, 0, 0},
    {0, 1, 2, 5, 0, 0, 0, 0},
    {3, 5, 0, 0, 0, 0, 0, 0},
    {0, 3, 5, 0, 0, 0, 0, 0},
    {1, 3, 5, 0, 0, 0, 0, 0},
    {0, 1, 3, 5, 0, 0, 0, 0},
    {2, 3, 5, 0, 0, 0, 0, 0},
    {0, 2, 3, 5, 0, 0, 0, 0},
    {1, 2, 3, 5, 0, 0, 0, 0},
    {0, 1, 2, 3, 5, 0, 0, 0},
    {4, 5, 0, 0, 0, 0, 0, 0},
    {0, 4, 5, 0, 0, 0, 0, 0},
    {1, 4, 5, 0, 0, 0, 0, 0},
    {0, 1, 4, 5, 0, 0, 0, 0},
    {2, 4, 5, 0, 0, 0, 0, 0},
    {0, 2, 4, 5, 0, 0, 0, 0},
    {1, 2, 4, 5, 0, 0, 0, 0},
    {0, 1, 2, 4, 5, 0, 0, 0},
    {3, 4, 5, 0, 0, 0, 0, 0},
    {0, 3, 4, 5, 0, 0, 0, 0},
    {1, 3, 4, 5, 0, 0, 0, 0},
    {0, 1, 3, 4, 5, 0, 0, 0},
    {2, 3, 4, 5, 0, 0, 0, 0},
    {0, 2, 3, 4, 5, 0, 0, 0},
    {1, 2, 3, 4, 5, 0, 0, 0},
    {0, 1, 2, 3, 4, 5, 0, 0},
    {6, 0, 0, 0, 0, 0, 0, 0},
    {0, 6, 0, 0, 0, 0, 0, 0},
    {1, 6, 0, 0, 0, 0, 0, 0},
    {0, 1, 6, 0, 0, 0, 0, 0},
    {2, 6, 0, 0, 0, 0, 0, 0},
    {0, 2, 6, 0, 0, 0, 0, 0},
    {1, 2, 6, 0, 0, 0, 0, 0},
    {0, 1, 2, 6, 0, 0, 0, 0},
    {3, 6, 0, 0, 0, 0, 0, 0},
    {0, 3, 6, 0, 0, 0, 0, 0},
    {1, 3, 6, 0, 0, 0, 0, 0},
    {0, 1, 3, 6, 0, 0, 0, 0},
    {2, 3, 6, 0, 0, 0, 0, 0},
    {0, 2, 3, 6, 0, 0, 0, 0},
    {1, 2, 3, 6, 0, 0, 0, 0},
    {0, 1, 2, 3, 6, 0, 0, 0},
    {4, 6, 0, 0, 0, 0, 0, 0},
    {0, 4, 6, 0, 0, 0, 0, 0},
    {1, 4, 6, 0, 0, 0, 0, 0},
    {0, 1, 4, 6, 0, 0, 0, 0},
    {2, 4, 6, 0, 0, 0, 0, 0},
    {0, 2, 4, 6, 0, 0, 0, 0},
    {1, 2, 4, 6, 0, 0, 0, 0},
    {0, 1, 2, 4, 6, 0, 0, 0},
    {3, 4, 6, 0, 0, 0, 0, 0},
    {0, 3, 4, 6, 0, 0, 0, 0},
    {1, 3, 4, 6, 0, 0, 0, 0},
    {0, 1, 3, 4, 6, 0, 0, 0},
    {2, 3, 4, 6, 0, 0, 0, 0},
    {0, 2, 3, 4, 6, 0, 0, 0},
    {1, 2, 3, 4, 6, 0, 0, 0},
    {0, 1, 2, 3, 4, 6, 0, 0},
    {5, 6, 0, 0, 0, 0, 0, 0},
    {0, 5, 6, 0, 0, 0, 0, 0},
    {1, 5, 6, 0, 0, 0, 0, 0},
    {0, 1, 5, 6, 0, 0, 0, 0},
    {2, 5, 6, 0, 0, 0, 0, 0},
    {0, 2, 5, 6, 0, 0, 0, 0},
    {1, 2, 5, 6, 0, 0, 0, 0},
    {0, 1, 2, 5, 6, 0, 0, 0},
    {3, 5, 6, 0, 0, 0, 0, 0},
    {0, 3, 5, 6, 0, 0, 0, 0},
    {1, 3, 5, 6, 0, 0, 0, 0},
    {0, 1, 3, 5, 6, 0, 0, 0},
    {2, 3, 5, 6, 0, 0, 0, 0},
    {0, 2, 3, 5, 6, 0, 0, 0},
    {1, 2, 3, 5, 6, 0, 0, 0},
    {0, 1, 2, 3, 5, 6, 0, 0},
    {4, 5, 6, 0, 0, 0, 0, 0},
    {0, 4, 5, 6, 0, 0, 0, 0},
    {1, 4, 5, 6, 0, 0, 0, 0},
    {0, 1, 4, 5, 6, 0, 0, 0},
    {2, 4, 5, 6, 0, 0, 0, 0},
    {0, 2, 4, 5, 6, 0, 0, 0},
    {1, 2, 4, 5, 6, 0, 0, 0},
    {0, 1, 2, 4, 5, 6, 0, 0},
    {3, 4, 5, 6, 0, 0, 0, 0},
    {0, 3, 4, 5, 6, 0, 0, 0},
    {1, 3, 4, 5, 6, 0, 0, 0},
    {0, 1, 3, 4, 5, 6, 0, 0},
    {2, 3, 4, 5, 6, 0, 0, 0},
    {0, 2, 3, 4, 5, 6, 0, 0},
    {1, 2, 3, 4, 5, 6, 0, 0},
    {0, 1, 2, 3, 4, 5, 6, 0},
    {7, 0, 0, 0, 0, 0, 0, 0},
    {0, 7, 0, 0, 0, 0, 0, 0},
    {1, 7, 0, 0, 0, 0, 0, 0},
    {0, 1, 7, 0, 0, 0, 0, 0},
    {2, 7, 0, 0, 0, 0, 0, 0},
    {0, 2, 7, 0, 0, 0, 0, 0},
    {1, 2, 7, 0, 0, 0, 0, 0},
    {0, 1, 2, 7, 0, 0, 0, 0},
    {3, 7, 0, 0, 0, 0, 0, 0},
    {0, 3, 7, 0, 0, 0, 0, 0},
    {1, 3, 7, 0, 0, 0, 0, 0},
    {0, 1, 3, 7, 0, 0, 0, 0},
    {2, 3, 7, 0, 0, 0, 0, 0},
    {0, 2, 3, 7, 0, 0, 0, 0},
    {1, 2, 3, 7, 0, 0, 0, 0},
    {0, 1, 2, 3, 7, 0, 0, 0},
    {4, 7, 0, 0, 0, 0, 0, 0},
    {0, 4, 7, 0, 0, 0, 0, 0},
    {1, 4, 7, 0, 0, 0, 0, 0},
    {0, 1, 4, 7, 0, 0, 0, 0},
    {2, 4, 7, 0, 0, 0, 0, 0},
    {0, 2, 4, 7, 0, 0, 0, 0},
    {1, 2, 4, 7, 0, 0, 0, 0},
    {0, 1, 2, 4, 7, 0, 0, 0},
    {3, 4, 7, 0, 0, 0, 0, 0},
    {0, 3, 4, 7, 0, 0, 0, 0},
    {1, 3, 4, 7, 0, 0, 0, 0},
    {0, 1, 3, 4, 7, 0, 0, 0},
    {2, 3, 4, 7, 0, 0, 0, 0},
    {0, 2, 3, 4, 7, 0, 0, 0},
    {1, 2, 3, 4, 7, 0, 0, 0},
    {0, 1, 2, 3, 4, 7, 0, 0},
    {5, 7, 0, 0, 0, 0, 0, 0},
    {0, 5, 7, 0, 0, 0, 0, 0},
    {1, 5, 7, 0, 0, 0, 0, 0},
    {0, 1, 5, 7, 0, 0, 0, 0},
    {2, 5, 7, 0, 0, 0, 0, 0},
    {0, 2, 5, 7, 0, 0, 0, 0},
    {1, 2, 5, 7, 0, 0, 0, 0},
    {0, 1, 2, 5, 7, 0, 0, 0},
    {3, 5, 7, 0, 0, 0, 0, 0},
    {0, 3, 5, 7, 0, 0, 0, 0},
    {1, 3, 5, 7, 0, 0, 0, 0},
    {0, 1, 3, 5, 7, 0, 0, 0},
    {2, 3, 5, 7, 0, 0, 0, 0},
    {0, 2, 3, 5, 7, 0, 0, 0},
    {1, 2, 3, 5, 7, 0, 0, 0},
    {0, 1, 2, 3, 5, 7, 0, 0},
    {4, 5, 7, 0, 0, 0, 0, 0},
    {0, 4, 5, 7, 0, 0, 0, 0},
    {1, 4, 5, 7, 0, 0, 0, 0},
    {0, 1, 4, 5, 7, 0, 0, 0},
    {2, 4, 5, 7, 0, 0, 0, 0},
    {0, 2, 4, 5, 7, 0, 0, 0},
    {1, 2, 4, 5, 7, 0, 0, 0},
    {0, 1, 2, 4, 5, 7, 0, 0},
    {3, 4, 5, 7, 0, 0, 0, 0},
    {0, 3, 4, 5, 7, 0, 0, 0},
    {1, 3, 4, 5, 7, 0, 0, 0},
    {0, 1, 3, 4, 5, 7, 0, 0},
    {2, 3, 4, 5, 7, 0, 0, 0},
    {0, 2, 3, 4, 5, 7, 0, 0},
    {1, 2, 3, 4, 5, 7, 0, 0},
    {0, 1, 2, 3, 4, 5, 7, 0},
    {6, 7, 0, 0, 0, 0, 0, 0},
    {0, 6, 7, 0, 0, 0, 0, 0},
    {1, 6, 7, 0, 0, 0, 0, 0},
    {0, 1, 6, 7, 0, 0, 0, 0},
    {2, 6, 7, 0, 0, 0, 0, 0},
    {0, 2, 6, 7, 0, 0, 0, 0},
    {1, 2, 6, 7, 0, 0, 0, 0},
    {0, 1, 2, 6, 7, 0, 0, 0},
    {3, 6, 7, 0, 0, 0, 0, 0},
    {0, 3, 6, 7, 0, 0, 0, 0},
    {1, 3, 6, 7, 0, 0, 0, 0},
    {0, 1, 3, 6, 7, 0, 0, 0},
    {2, 3, 6, 7, 0, 0, 0, 0},
    {0, 2, 3, 6, 7, 0, 0, 0},
    {1, 2, 3, 6, 7, 0, 0, 0},
    {0, 1, 2, 3, 6, 7, 0, 0},
    {4, 6, 7, 0, 0, 0, 0, 0},
    {0, 4, 6, 7, 0, 0, 0, 0},
    {1, 4, 6, 7, 0, 0, 0, 0},
    {0, 1, 4, 6, 7, 0, 0, 0},
    {2, 4, 6, 7, 0, 0, 0, 0},
    {0, 2, 4, 6, 7, 0, 0, 0},
    {1, 2, 4, 6, 7, 0, 0, 0},
    {0, 1, 2, 4, 6, 7, 0, 0},
    {3, 4, 6, 7, 0, 0, 0, 0},
    {0, 3, 4, 6, 7, 0, 0, 0},
    {1, 3, 4, 6, 7, 0, 0, 0},
    {0, 1, 3, 4, 6, 7, 0, 0},
    {2, 3, 4, 6, 7, 0, 0, 0},
    {0, 2, 3, 4, 6, 7, 0, 0},
    {1, 2, 3, 4, 6, 7, 0, 0},
    {0, 1, 2, 3, 4, 6, 7, 0},
    {5, 6, 7, 0, 0, 0, 0, 0},
    {0, 5, 6, 7, 0, 0, 0, 0},
    {1, 5, 6, 7, 0, 0, 0, 0},
    {0, 1, 5, 6, 7, 0, 0, 0},
    {2, 5, 6, 7, 0, 0, 0, 0},
    {0, 2, 5, 6, 7, 0, 0, 0},
    {1, 2, 5, 6, 7, 0, 0, 0},
    {0, 1, 2, 5, 6, 7, 0, 0},
    {3, 5, 6, 7, 0, 0, 0, 0},
    {0, 3, 5, 6, 7, 0, 0, 0},
    {1, 3, 5, 6, 7, 0, 0, 0},
    {0, 1, 3, 5, 6, 7, 0, 0},
    {2, 3, 5, 6, 7, 0, 0, 0},
    {0, 2, 3, 5, 6, 7, 0, 0},
    {1, 2, 3, 5, 6, 7, 0, 0},
    {0, 1, 2, 3, 5, 6, 7, 0},
    {4, 5, 6, 7, 0, 0, 0, 0},
    {0, 4, 5, 6, 7, 0, 0, 0},
    {1, 4, 5, 6, 7, 0, 0, 0},
    {0, 1, 4, 5, 6, 7, 0, 0},
    {2, 4, 5, 6, 7, 0, 0, 0},
    {0, 2, 4, 5, 6, 7, 0, 0},
    {1, 2, 4, 5, 6, 7, 0, 0},
    {0, 1, 2, 4, 5, 6, 7, 0},
    {3, 4, 5, 6, 7, 0, 0, 0},
    {0, 3, 4, 5, 6, 7, 0, 0},
    {1, 3, 4, 5, 6, 7, 0, 0},
    {0, 1, 3, 4, 5, 6, 7, 0},
    {2, 3, 4, 5, 6, 7, 0, 0},
    {0, 2, 3, 4, 5, 6, 7, 0},
    {1, 2, 3, 4, 5, 6, 7, 0},
    {0, 1, 2, 3, 4, 5, 6, 7}};

void poly_caddq(poly *a)
{

  unsigned int i;
  __m256i f, g;
  const __m256i q = _mm256_set1_epi32(Q);
  const __m256i zero = _mm256_setzero_si256();

  for (i = 0; i < N / 8; i++)
  {
    f = _mm256_load_si256(&a->vec[i]);
    g = _mm256_blendv_epi32(zero, q, f);
    f = _mm256_add_epi32(f, g);
    _mm256_store_si256(&a->vec[i], f);
  }
}

void poly_add(poly *c, const poly *a, const poly *b)
{
  unsigned int i;
  __m256i f, g;

  for (i = 0; i < N / 8; i++)
  {
    f = _mm256_load_si256(&a->vec[i]);
    g = _mm256_load_si256(&b->vec[i]);
    f = _mm256_add_epi32(f, g);
    _mm256_store_si256(&c->vec[i], f);
  }
}

void poly_modadd(poly *c, poly *a, poly *b)
{
  unsigned int i;

  for (i = 0; i < N; i++)
    c->coeffs[i] = mod_add(a->coeffs[i], b->coeffs[i]);
}

void poly_sub(poly *c, poly *a, poly *b)
{
  unsigned int i;
  __m256i f, g;

  for (i = 0; i < N / 8; i++)
  {
    f = _mm256_load_si256(&a->vec[i]);
    g = _mm256_load_si256(&b->vec[i]);
    f = _mm256_sub_epi32(f, g);
    _mm256_store_si256(&c->vec[i], f);
  }
}

void poly_modsub(poly *c, poly *a, poly *b)
{
  unsigned int i;
  for (i = 0; i < N; i++)
    c->coeffs[i] = mod_sub(a->coeffs[i], b->coeffs[i]);
}

void poly_shiftl(poly *a)
{
  unsigned int i;
  __m256i f;

  for (i = 0; i < N / 8; i++)
  {
    f = _mm256_load_si256(&a->vec[i]);
    f = _mm256_slli_epi32(f, D);
    _mm256_store_si256(&a->vec[i], f);
  }
}

void poly_power2round(poly *a1, poly *a0, const poly *a)
{
  power2round_avx(a1->vec, a0->vec, a->vec);
}

void poly_decompose(poly *a1, poly *a0, const poly *a)
{
  decompose_avx(a1->vec, a0->vec, a->vec);
}
unsigned int poly_make_hint(poly *h, const poly *a0, const poly *a1)
{
  unsigned int i, s = 0;

  for (i = 0; i < N; i++)
  {
    h->coeffs[i] = make_hint(a0->coeffs[i], a1->coeffs[i]);
    s += h->coeffs[i];
  }

  return s;
}

void poly_use_hint(poly *b, const poly *a, const poly *h)
{
  use_hint_avx(b->vec, a->vec, h->vec);
}

int poly_chknorm(poly *a, int32_t B)
{
  unsigned int i;
  int r;
  __m256i f, t;
  const __m256i bound = _mm256_set1_epi32(B - 1);

  if (B > (Q - 1) / 8)
    return 1;

  t = _mm256_setzero_si256();
  for (i = 0; i < N / 8; i++)
  {
    f = _mm256_load_si256(&a->vec[i]);
    f = _mm256_abs_epi32(f);
    f = _mm256_cmpgt_epi32(f, bound);
    t = _mm256_or_si256(t, f);
  }

  r = 1 - _mm256_testz_si256(t, t);
  return r;
}

#if N == 2048
#define POLY_UNIFORM_NBLOCKS ((3*N + STREAM128_BLOCKBYTES - 1)/STREAM128_BLOCKBYTES)
#else
#define POLY_UNIFORM_NBLOCKS ((6*N + STREAM128_BLOCKBYTES - 1)/STREAM128_BLOCKBYTES)
#endif
#define POLY_UNIFORM_BUFLEN (POLY_UNIFORM_NBLOCKS*STREAM128_BLOCKBYTES)

static unsigned int rej_eta(int32_t *a, unsigned int len, const uint8_t *buf, unsigned int buflen)
{
  unsigned int ctr, pos;
  uint32_t t0, t1;
  uint32_t t2, t3;

  ctr = pos = 0;
  while (ctr < len && pos < buflen)
  {

    t0 = buf[pos] & 0x03;
    t1 = (buf[pos] >> 2) & 0x03;
    t2 = (buf[pos] >> 4) & 0x03;
    t3 = buf[pos++] >> 6;
    if (t0 < 3)
    {
      a[ctr++] = 1 - t0;
    }
    if (t1 < 3 && ctr < len)
    {
      a[ctr++] = 1 - t1;
    }
    if (t2 < 3 && ctr < len)
    {
      a[ctr++] = 1 - t2;
    }
    if (t3 < 3 && ctr < len)
    {
      a[ctr++] = 1 - t3;
    }
  }
  return ctr;
}

#if N == 1152
#define POLY_UNIFORM_ETA_NBLOCKS ((3* 136 + STREAM256_BLOCKBYTES - 1)/STREAM256_BLOCKBYTES)
#elif N == 1536
#define POLY_UNIFORM_ETA_NBLOCKS ((4* 136 + STREAM256_BLOCKBYTES - 1)/STREAM256_BLOCKBYTES)
#else
#define POLY_UNIFORM_ETA_NBLOCKS ((6 * 136 + STREAM256_BLOCKBYTES - 1)/STREAM256_BLOCKBYTES)
#endif
#define POLY_UNIFORM_ETA_BUFLEN (POLY_UNIFORM_ETA_NBLOCKS*STREAM256_BLOCKBYTES)

static void merge_alternate_m256(__m256i *g1, __m256i *g2, __m256i *result)
{
  __m256i g1_low = _mm256_loadu_si256(&g1[0]);
  __m256i g1_high = _mm256_loadu_si256(&g1[1]);
  __m256i g2_low = _mm256_loadu_si256(&g2[0]);
  __m256i g2_high = _mm256_loadu_si256(&g2[1]);

  __m256i result1 = _mm256_unpacklo_epi16(g1_low, g2_low);
  __m256i result2 = _mm256_unpackhi_epi16(g1_low, g2_low);
  __m256i result3 = _mm256_unpacklo_epi16(g1_high, g2_high);
  __m256i result4 = _mm256_unpackhi_epi16(g1_high, g2_high);

  _mm256_storeu_si256(&result[0], result1);
  _mm256_storeu_si256(&result[1], result2);
  _mm256_storeu_si256(&result[2], result3);
  _mm256_storeu_si256(&result[3], result4);
}
static unsigned int rej_eta_avx(int32_t *restrict r, const uint8_t buf[POLY_UNIFORM_ETA_BUFLEN])
{

  unsigned int ctr, pos;
  uint32_t good;
  __m256i f0, f1, f2, f3, f4;
  __m128i g0, g1;
  const __m256i masked = _mm256_set1_epi8(3);
  const __m256i eta = _mm256_set1_epi8(1);
  const __m256i bound = _mm256_set1_epi8(3);

  ctr = pos = 0;
  while (ctr <= N - 8 && pos <= POLY_UNIFORM_ETA_BUFLEN - 16)
  {
    f0 = _mm256_cvtepu8_epi16(_mm_loadu_si128((__m128i *)&buf[pos]));

    f1 = _mm256_slli_epi16(f0, 6);

    f2 = _mm256_slli_epi16(f0, 2);
    f3 = _mm256_srli_epi16(f0, 4);
    f2 = _mm256_or_si256(f2, f3);
    f2 = _mm256_and_si256(f2, masked);

    f0 = _mm256_or_si256(f0, f1);
    f0 = _mm256_and_si256(f0, masked);
    merge_alternate_m256(&f0, &f2, &f4);

    f1 = _mm256_sub_epi8(f4, bound);
    f0 = _mm256_sub_epi8(eta, f0);
    f2 = _mm256_sub_epi8(eta, f2);
    merge_alternate_m256(&f0, &f2, &f4);
    good = _mm256_movemask_epi8(f1);

    g0 = _mm256_castsi256_si128(f4);
    g1 = _mm_loadl_epi64((__m128i *)&idxlut[good & 0xFF]);
    g1 = _mm_shuffle_epi8(g0, g1);
    f1 = _mm256_cvtepi8_epi32(g1);
    _mm256_storeu_si256((__m256i *)&r[ctr], f1);
    ctr += _mm_popcnt_u32(good & 0xFF);
    good >>= 8;
    pos += 2;

    if(ctr > N - 8) break;
    g0 = _mm_bsrli_si128(g0,8);
    g1 = _mm_loadl_epi64((__m128i *)&idxlut[good & 0xFF]);
    g1 = _mm_shuffle_epi8(g0, g1);
    f1 = _mm256_cvtepi8_epi32(g1);
    _mm256_storeu_si256((__m256i *)&r[ctr], f1);
    ctr += _mm_popcnt_u32(good & 0xFF);
    pos += 2;
  }

  uint32_t t0, t1, t2, t3;
  while (ctr < N && pos < POLY_UNIFORM_ETA_BUFLEN)
  {

    t0 = buf[pos] & 0x03;
    t1 = (buf[pos] >> 2) & 0x03;
    t2 = (buf[pos] >> 4) & 0x03;
    t3 = buf[pos++] >> 6;
    if (t0 < 3)
    {
      r[ctr++] = 1 - t0;
    }
    if (t1 < 3 && ctr < N)
    {
      r[ctr++] = 1 - t1;
    }
    if (t2 < 3 && ctr < N)
    {
      r[ctr++] = 1 - t2;
    }
    if (t3 < 3 && ctr < N)
    {
      r[ctr++] = 1 - t3;
    }
  }
  return ctr;
}
void poly_uniform_eta(poly *a, const uint8_t seed[CRHBYTES], uint16_t nonce)
{
  unsigned int ctr;
  uint8_t buf[POLY_UNIFORM_ETA_NBLOCKS * STREAM256_BLOCKBYTES];
  stream256_state state;

  stream256_init(&state, seed, nonce);
  stream256_squeezeblocks(buf, POLY_UNIFORM_ETA_NBLOCKS, &state);

  ctr = rej_eta_avx(a->coeffs, buf);
  while (ctr < N)
  {
    stream256_squeezeblocks(buf, 1, &state);
    ctr += rej_eta(a->coeffs + ctr, N - ctr, buf, STREAM256_BLOCKBYTES);
  }
}

#define POLY_UNIFORM_GAMMA1_NBLOCKS ((POLYZ_PACKEDBYTES + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES)

void poly_uniform_gamma1(poly *a, const uint8_t seed[CRHBYTES], uint16_t nonce)
{
  uint8_t buf[POLY_UNIFORM_GAMMA1_NBLOCKS * STREAM256_BLOCKBYTES + 14];
  stream256_state state;

  stream256_init(&state, seed, nonce);
  stream256_squeezeblocks(buf, POLY_UNIFORM_GAMMA1_NBLOCKS, &state);
  polyz_unpack(a, buf);
}

// sample in ball
void poly_challenge(poly *c, const uint8_t seed[SEEDBYTES])
{
  unsigned int i, b, pos;
  uint32_t signs;
  uint8_t buf[SHAKE256_RATE];
  keccak_state state;

  shake256_init(&state);
  shake256_absorb(&state, seed, SEEDBYTES);
  shake256_finalize(&state);
  shake256_squeezeblocks(buf, 1, &state);

  signs = 0;
  for (i = 0; i < 4; i++)
    signs |= (uint32_t)buf[i] << 8 * i;
  pos = 4;

  for (i = 0; i < N; i++)
    c->coeffs[i] = 0;
#if N == 2304
  for (i = N - TAU; i < N; i++)
  {
    do
    {
      if (pos >= SHAKE256_RATE)
      {
        shake256_squeezeblocks(buf, 1, &state);
        pos = 0;
      }

      b = (uint32_t)buf[pos++] << 4;
      b |= (buf[pos++] & 0xF);
    } while (b > i);

    c->coeffs[i] = c->coeffs[b];
    c->coeffs[b] = 1 - 2 * (signs & 1);
    signs >>= 1;
  }
#else
  for (i = N - TAU; i < N; i++)
  {
    do
    {
      if (pos >= SHAKE256_RATE)
      {
        shake256_squeezeblocks(buf, 1, &state);
        pos = 0;
      }

      b = (uint32_t)buf[pos++] << 3;
      b |= (buf[pos++] & 0x7);
    } while (b > i);

    c->coeffs[i] = c->coeffs[b];
    c->coeffs[b] = 1 - 2 * (signs & 1);
    signs >>= 1;
  }
#endif
}

uint8_t convToIdx(uint16_t* res, const uint8_t res_length, const int32_t* op,
    const size_t op_length) {
    uint8_t index = 0, b = 0;
    uint8_t index_arr[2] = { 0, res_length - 1 }; // 0 for positive, 1 for
    // negative
    for (size_t i = 0; i < op_length; ++i) {


        index = ((op[i] & 0x80000000) >> 31) & 0x00000001;
        b = (~(uint64_t)(uint32_t)op[i] + 1) >> 63;
        res[index_arr[index]] ^= (~b + 1) & (res[index_arr[index]] ^ i);
        index_arr[index] += op[i];
    }

    return index_arr[0];
}


static void poly_sparse_add(poly_sparse *res, const poly *op1, const uint16_t deg)
{
  __m256i v_res, v_op1;
  for (size_t i = 0; i < N / 8; ++i)
  {
    v_res = _mm256_loadu_si256((__m256i *)&res->coeffs[deg + 8 * i]);
    v_op1 = _mm256_loadu_si256((__m256i *)&op1->coeffs[8 * i]);
    v_res = _mm256_add_epi32(v_res, v_op1);
    _mm256_storeu_si256((__m256i *)&res->coeffs[deg + 8 * i], v_res);
  }
}

static void poly_sparse_sub(poly_sparse *res, const poly *op1, const uint16_t deg)
{
  __m256i v_res, v_op1;
  for (size_t i = 0; i < N / 8; ++i)
  {
    v_res = _mm256_loadu_si256((__m256i *)&res->coeffs[deg + 8 * i]);
    v_op1 = _mm256_loadu_si256((__m256i *)&op1->coeffs[8 * i]);
    v_res = _mm256_sub_epi32(v_res, v_op1);
    _mm256_storeu_si256((__m256i *)&res->coeffs[deg + 8 * i], v_res);
  }
}

static void poly_sparse_reduce(poly *res, poly_sparse *c)
{
  uint32_t i;
  for (i = N + (N >> 1) - 1; i < 2 * N - 1; i++)
  {
    c->coeffs[i - (N >> 1)] = (c->coeffs[i - (N >> 1)] + c->coeffs[i]);
    c->coeffs[i - N] = (c->coeffs[i - N] - c->coeffs[i]);
  }
  for (i = N; i < N + (N >> 1) - 1; i++)
  {
    c->coeffs[i - (N >> 1)] = (c->coeffs[i - (N >> 1)] + c->coeffs[i]);
    c->coeffs[i - N] = (c->coeffs[i - N] - c->coeffs[i]);
  }

  for (i = 0; i < N / 8; i++)
  {
    res->vec[i] = c->vec[i];
  }
}
void poly_mult_add(poly *res, const poly *op1, const uint16_t *op2, const  uint8_t neg_start) {
    poly_sparse temp = {0};

    for (size_t j = 0; j < neg_start; ++j) {
        poly_sparse_add(&temp, op1, op2[j]);
    }

    for (size_t j = neg_start; j < TAU; ++j) {
        poly_sparse_sub(&temp, op1, op2[j]);
    }
    
    poly_sparse_reduce(res, &temp);
    
}
void poly_reduce(poly *a)
{
  unsigned int i;

  __m256i f, g;
  const __m256i q = _mm256_set1_epi32(Q);
  const __m256i off = _mm256_set1_epi32(1 << 22);

  for (i = 0; i < N / 8; i++)
  {
    f = _mm256_load_si256(&a->vec[i]);
    g = _mm256_add_epi32(f, off);
    g = _mm256_srai_epi32(g, 23);
    g = _mm256_mullo_epi32(g, q);
    f = _mm256_sub_epi32(f, g);
    _mm256_store_si256(&a->vec[i], f);
  }
}

static unsigned int rej_uniform_avx(int32_t * restrict r, const uint8_t buf[POLY_UNIFORM_BUFLEN+8])
{
  unsigned int ctr, pos;
  uint32_t good;
  __m256i d, tmp;
  const __m256i bound = _mm256_set1_epi32(Q);
  #if NIMS_TRI_NTT_MODE == 5
	const __m256i local_mask  = _mm256_set1_epi32(0x7FFFFF);
  #else
	const __m256i local_mask  = _mm256_set1_epi32(0xFFFFFF);
  #endif
  
  const __m256i idx8  = _mm256_set_epi8(-1,15,14,13,-1,12,11,10,
                                        -1, 9, 8, 7,-1, 6, 5, 4,
                                        -1,11,10, 9,-1, 8, 7, 6,
                                        -1, 5, 4, 3,-1, 2, 1, 0);

  ctr = pos = 0;
  while(pos <= POLY_UNIFORM_BUFLEN - 24) {
    d = _mm256_loadu_si256((__m256i *)&buf[pos]);
    d = _mm256_permute4x64_epi64(d, 0x94);
    d = _mm256_shuffle_epi8(d, idx8);
    d = _mm256_and_si256(d, local_mask);
    pos += 24;

    tmp = _mm256_sub_epi32(d, bound);
    good = _mm256_movemask_ps((__m256)tmp);
    tmp = _mm256_cvtepu8_epi32(_mm_loadl_epi64((__m128i *)&idxlut[good]));
    d = _mm256_permutevar8x32_epi32(d, tmp);

    _mm256_storeu_si256((__m256i *)&r[ctr], d);
    ctr += _mm_popcnt_u32(good);

    if(ctr > N - 8) break;
  }

  uint32_t t;

#if NIMS_TRI_NTT_MODE == 5
	while(ctr < N && pos <= POLY_UNIFORM_BUFLEN - 3) {
		t  = buf[pos++];
		t |= (uint32_t)buf[pos++] << 8;
		t |= (uint32_t)buf[pos++] << 16;
		t &= 0x7FFFFF;

		if(t < Q)
		r[ctr++] = t;
	}

	return ctr;
#else
	while(ctr < (N) && pos <= POLY_UNIFORM_BUFLEN - 3) 
	{
		t  = buf[pos++];
		t |= (uint32_t)buf[pos++] << 8;
		t |= (uint32_t)buf[pos++] << 16;
		t &= 0xFFFFFF;

		if(t < Q)
		r[ctr++] = t;
	}
  	return ctr;
#endif
}

static unsigned int rej_uniform(int32_t *a,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen)
{
  unsigned int ctr, pos;
  uint32_t t;
  
#if NIMS_TRI_NTT_MODE == 5
	ctr = pos = 0;
	while (ctr < len && pos + 3 <= buflen) {
		t = buf[pos++];
		t |= (uint32_t)buf[pos++] << 8;
		t |= (uint32_t)buf[pos++] << 16;
		t &= 0x7FFFFF;	// Q -> 23-bit
		if (t < Q && t != 0)
			a[ctr++] = t;
	}
#else
  ctr = pos = 0;
  while (ctr < len && pos + 3 <= buflen)
  {
    t = buf[pos++];
    t |= (uint32_t)buf[pos++] << 8;
    t |= (uint32_t)buf[pos++] << 16;
    t &= 0xFFFFFF;

    if (t < Q)
      a[ctr++] = t;
  }
#endif

  return ctr;
}
void poly_uniform(poly *a, const uint8_t seed[SEEDBYTES], uint16_t nonce)
{
  unsigned int ctr; 
  uint8_t buf[POLY_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES + 2] = {0};
  stream128_state state;
  stream128_init(&state, seed, nonce);
  stream128_squeezeblocks(buf, POLY_UNIFORM_NBLOCKS, &state);

  ctr = rej_uniform_avx(a->coeffs, buf);

  while (ctr < N)
  {
    stream128_squeezeblocks(buf, 1, &state);
    ctr += rej_uniform(a->coeffs + ctr, N - ctr, buf, STREAM128_BLOCKBYTES);
  }
}
