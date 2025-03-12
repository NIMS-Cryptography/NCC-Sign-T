#include <stdio.h>
#include <stdint.h>
#include "params.h"
#include "reduce.h"

int32_t caddq(int32_t a)
{
    a += (a >> 31) & Q;
    return a;
}

int32_t csubq(int32_t a)
{
    a -= Q;
    a += ((int32_t)a >> 31) & Q;
    return a;
}

int32_t freeze(int32_t a)
{
    a = reduce32(a);
    a = caddq(a);
    return a;
}

int32_t mod_add(int32_t a, int32_t b)
{
    int32_t t;
    t = (a + b);
    t = t - Q;
    t += (t >> 31) & Q;
    t += (t >> 31) & Q;
    t += (t >> 31) & Q;

    return (uint32_t)t;
}

int32_t mod_sub(int32_t a, int32_t b)
{
    int32_t t;
    t = a - b;
    t = t - Q;
    t += (t >> 31) & Q;
    t += (t >> 31) & Q;
    t += (t >> 31) & Q;

    return (uint32_t)t;
}
/*************************************************
 * Name:        reduce32
 *
 * Description: For finite field element a with a <= 2^{31} - 2^{22} - 1,
 *              compute r \equiv a (mod Q) such that -6283009 <= r <= 6283007.
 *
 * Arguments:   - int32_t: finite field element a
 *
 * Returns r.
 **************************************************/
int32_t reduce32(int32_t a)
{
    int32_t t;

    t = (a + (1 << 22)) >> 23;
    t = a - t * Q;
    // t = a % Q;
    return t;
}
#define _mm256_blendv_epi32(a, b, mask)                          \
    _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(a), \
                                         _mm256_castsi256_ps(b), \
                                         _mm256_castsi256_ps(mask)));

__m256i reduce32_avx(__m256i a)
{
    __m256i t, temp;
    __m256i Qx4 = _mm256_set1_epi64x(Q);
    __m256i _22 = _mm256_set_epi32(0, 0, 0, 0x400000, 0, 0x400000, 0, 0x400000);
    __m256i t1 = _mm256_set1_epi64x((int64_t)0x7fffff00000000);

    temp = _mm256_add_epi32(a, _22);
    t = _mm256_slli_epi64(temp, 32);
    temp = _mm256_blendv_epi32(temp, t1, t);
    t = _mm256_srli_epi64(temp, 23);

    t = _mm256_sub_epi32(a, _mm256_mullo_epi32(t, Qx4));

    return t;
}

__m256i reduce32_avx_4(__m256i a)
{
    __m256i t, temp;
    __m256i Qx4 = _mm256_set1_epi64x(Q);
    __m256i _22 = _mm256_set_epi32(0, 0x400000, 0, 0x400000, 0, 0x400000, 0, 0x400000);
    __m256i t1 = _mm256_set1_epi64x((int64_t)0x7fffff00000000);

    temp = _mm256_add_epi32(a, _22);
    t = _mm256_slli_epi64(temp, 32);
    temp = _mm256_blendv_epi32(temp, t1, t);
    t = _mm256_srli_epi64(temp, 23);

    t = _mm256_sub_epi32(a, _mm256_mullo_epi32(t, Qx4));

    return t;
}
