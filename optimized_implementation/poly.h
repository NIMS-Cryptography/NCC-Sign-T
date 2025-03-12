#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include <immintrin.h>
#include "params.h"

typedef union _poly
{
    int32_t coeffs[N];
    __m256i vec[(N + 7) / 8 + 4];
} poly;

typedef union _poly2 {
    int32_t coeffs[N << 1];
    __m256i vec[((N+7)/8+4) << 1];
} poly_sparse;
extern const uint8_t idxlut[256][8];
void ntt_avx_asm(poly *Out, poly *A);
void poly_base_mul_avx_asm(poly *c, poly *a, poly *b);
void invntt_avx_asm(poly *Out, poly *A);
uint8_t convToIdx(uint16_t* res, const uint8_t res_length, const int32_t* op,
    const size_t op_length);
void poly_mult_add(poly *res, const poly *op1, const uint16_t *op2, const  uint8_t neg_start);

void poly_modadd(poly *c, poly *a, poly *b);
void poly_modsub(poly *c, poly *a, poly *b);

void poly_reduce(poly *a);
void poly_caddq(poly *a);

void poly_add(poly *c, const poly *a, const poly *b);
void poly_sub(poly *c, poly *a, poly *b);
void poly_shiftl(poly *a);

void poly_power2round(poly *a1, poly *a0, const poly *a);
void poly_decompose(poly *a1, poly *a0, const poly *a);
unsigned int poly_make_hint(poly *h, const poly *a0, const poly *a1);
void poly_use_hint(poly *b, const poly *a, const poly *h);
int poly_chknorm(poly *a, int32_t B);
void poly_uniform(poly *a,
                      const uint8_t seed[SEEDBYTES],
                      uint16_t nonce);
void poly_uniform_eta(poly *a,
                      const uint8_t seed[CRHBYTES],
                      uint16_t nonce);
void poly_uniform_gamma1(poly *a,
                         const uint8_t seed[CRHBYTES],
                         uint16_t nonce);
void poly_challenge(poly *c, const uint8_t seed[SEEDBYTES]);

#define _mm256_blendv_epi32(a, b, mask)                          \
    _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(a), \
                                         _mm256_castsi256_ps(b), \
                                         _mm256_castsi256_ps(mask)));


#endif
