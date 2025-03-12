#include <stdint.h>
#include "params.h"
#include "rounding.h"
#include "poly.h"
void power2round_avx(__m256i *a1, __m256i *a0, const __m256i *a)
{
  unsigned int i;
  __m256i f, f0, f1;
  const __m256i mask = _mm256_set1_epi32(-(1 << D));
  const __m256i half = _mm256_set1_epi32((1 << (D - 1)) - 1);

  for (i = 0; i < N / 8; ++i)
  {
    f = _mm256_load_si256(&a[i]);
    f1 = _mm256_add_epi32(f, half);
    f0 = _mm256_and_si256(f1, mask);
    f1 = _mm256_srli_epi32(f1, D);
    f0 = _mm256_sub_epi32(f, f0);
    _mm256_store_si256(&a1[i], f1);
    _mm256_store_si256(&a0[i], f0);
  }
}


static __m256i mulhi_epu32(__m256i a, __m256i b)
{

  __m256i mul1 = _mm256_mul_epu32(a, b);
  __m256i a_shifted = _mm256_srli_epi64(a, 32);
  __m256i b_shifted = _mm256_srli_epi64(b, 32);
  __m256i mul2 = _mm256_mul_epu32(a_shifted, b_shifted);

  __m256i result1 = _mm256_srli_epi64(mul1, 32);
  __m256i result2 = _mm256_srli_epi64(mul2, 32);

  __m256i result = _mm256_blend_epi32(result1, _mm256_slli_epi64(result2, 32), 0xAA);

  return result;
}

void decompose_avx(__m256i *a1, __m256i *a0, const __m256i *a)
{

  #if N == 2304
     unsigned int i;
  __m256i f,f0,f1;
  const __m256i q = _mm256_set1_epi32(Q);
  const __m256i hq = _mm256_srli_epi32(q,1);
  const __m256i v = _mm256_set1_epi32(1022);
  const __m256i alpha = _mm256_set1_epi32(2*GAMMA2);
  const __m256i off = _mm256_set1_epi32(127);
  const __m256i shift = _mm256_set1_epi32(512);
  const __m256i mask = _mm256_set1_epi32(15);

  for(i=0;i<N/8;i++) {
    f = _mm256_load_si256(&a[i]);
    f1 = _mm256_add_epi32(f,off);
    f1 = _mm256_srli_epi32(f1,7);
    f1 = _mm256_mulhi_epu16(f1,v);
    f1 = _mm256_mulhrs_epi16(f1,shift);
    f1 = _mm256_and_si256(f1,mask);
    f0 = _mm256_mullo_epi32(f1,alpha);
    f0 = _mm256_sub_epi32(f,f0);
    f = _mm256_cmpgt_epi32(f0,hq);
    f = _mm256_and_si256(f,q);
    f0 = _mm256_sub_epi32(f0,f);
    _mm256_store_si256(&a1[i],f1);
    _mm256_store_si256(&a0[i],f0);
  }
 
#elif N == 1536
    unsigned int i;
    __m256i f,f0,f1;
    const __m256i q = _mm256_set1_epi32(Q);
    const __m256i hq = _mm256_srli_epi32(q,1);
    const __m256i v = _mm256_set1_epi32(1047489);
    const __m256i alpha = _mm256_set1_epi32(2*GAMMA2);
    const __m256i off = _mm256_set1_epi32(7);
    const __m256i shift = _mm256_set1_epi32(4096);
    const __m256i mask = _mm256_set1_epi32(31);

  for(i=0;i<N/8;i++) {
    f = _mm256_load_si256(&a[i]);
    f1 = _mm256_add_epi32(f,off);
    f1 = _mm256_srli_epi32(f1,3);
    f1 = mulhi_epu32(f1,v);
    f1 = _mm256_mulhrs_epi16(f1,shift);
    f1 = _mm256_and_si256(f1,mask);
    f0 = _mm256_mullo_epi32(f1,alpha);
    f0 = _mm256_sub_epi32(f,f0);
    f = _mm256_cmpgt_epi32(f0,hq);
    f = _mm256_and_si256(f,q);
    f0 = _mm256_sub_epi32(f0,f);
    _mm256_store_si256(&a1[i],f1);
    _mm256_store_si256(&a0[i],f0);
  }

#elif N == 1152
     unsigned int i;
    __m256i f,f0,f1;
    const __m256i q = _mm256_set1_epi32(Q);
    const __m256i hq = _mm256_srli_epi32(q,1);
    const __m256i v = _mm256_set1_epi32(4187849);
    const __m256i alpha = _mm256_set1_epi32(2*GAMMA2);
    const __m256i off = _mm256_set1_epi32(1);
    const __m256i shift = _mm256_set1_epi32(256);
    const __m256i mask = _mm256_set1_epi32(31);

  for(i=0;i<N/8;i++) {
    f = _mm256_load_si256(&a[i]);
    f1 = _mm256_add_epi32(f,off);
    f1 = _mm256_srli_epi32(f1,1);
    f1 = mulhi_epu32(f1,v);
    f1 = _mm256_mulhrs_epi16(f1,shift);
    f1 = _mm256_and_si256(f1,mask);
    f0 = _mm256_mullo_epi32(f1,alpha);
    f0 = _mm256_sub_epi32(f,f0);
    f = _mm256_cmpgt_epi32(f0,hq);
    f = _mm256_and_si256(f,q);
    f0 = _mm256_sub_epi32(f0,f);
    _mm256_store_si256(&a1[i],f1);
    _mm256_store_si256(&a0[i],f0);
  }
#elif N == 2048
	unsigned int i;
  	__m256i f,f0,f1;
	const __m256i q = _mm256_set1_epi32(Q);
	const __m256i hq = _mm256_srli_epi32(q,1);
	const __m256i v = _mm256_set1_epi32(1025);
	const __m256i alpha = _mm256_set1_epi32(2*GAMMA2);
	const __m256i off = _mm256_set1_epi32(127);
	const __m256i shift = _mm256_set1_epi32(1024);
	const __m256i mask = _mm256_set1_epi32(31);

	for(i=0;i<N/8;i++) {
		f = _mm256_load_si256(&a[i]);
		f1 = _mm256_add_epi32(f,off);
		f1 = _mm256_srli_epi32(f1,7);
		f1 = _mm256_mulhi_epu16(f1,v);
		f1 = _mm256_mulhrs_epi16(f1,shift);
		f1 = _mm256_and_si256(f1,mask);
		f0 = _mm256_mullo_epi32(f1,alpha);
		f0 = _mm256_sub_epi32(f,f0);
		f = _mm256_cmpgt_epi32(f0,hq);
		f = _mm256_and_si256(f,q);
		f0 = _mm256_sub_epi32(f0,f);
		_mm256_store_si256(&a1[i],f1);
		_mm256_store_si256(&a0[i],f0);
	}
#endif
}

unsigned int make_hint(int32_t a0, int32_t a1)
{
  if (a0 > GAMMA2 || a0 < -GAMMA2 || (a0 == -GAMMA2 && a1 != 0))
    return 1;

  return 0;
}

void use_hint_avx(__m256i *b, const __m256i *a, const __m256i *restrict hint)
{
  unsigned int i;
  __m256i a0[N / 8];
  __m256i f, g, h, t;
  const __m256i zero = _mm256_setzero_si256();
#if GAMMA2 == (Q - 1) / 32
  const __m256i mask = _mm256_set1_epi32(15);
#elif GAMMA2 == (Q - 1) / 64
  const __m256i max = _mm256_set1_epi32(31);
#endif

  decompose_avx(b, a0, a);
  for (i = 0; i < N / 8; i++)
  {
    f = _mm256_load_si256(&a0[i]);
    g = _mm256_load_si256(&b[i]);
    h = _mm256_load_si256(&hint[i]);
    t = _mm256_blendv_epi32(zero, h, f);
    t = _mm256_slli_epi32(t, 1);
    h = _mm256_sub_epi32(h, t);
    g = _mm256_add_epi32(g, h);

#if GAMMA2 == (Q - 1) / 32
    g = _mm256_and_si256(g, mask);
#elif GAMMA2 == (Q - 1) / 64
    g = _mm256_blendv_epi32(g, max, g);
    f = _mm256_cmpgt_epi32(g, max);
    g = _mm256_blendv_epi32(g, zero, f);
#endif

    _mm256_store_si256(&b[i], g);
  }
}
