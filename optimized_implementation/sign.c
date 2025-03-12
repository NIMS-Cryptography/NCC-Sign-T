#include <stdint.h>
#include "params.h"
#include "sign.h"
#include "packing.h"
#include "poly.h"
#include "randombytes.h"
#include "symmetric.h"
#include "fips202.h"
#include "stdio.h"
#include <stdlib.h>
#include "consts.h"

int crypto_sign_keypair(uint8_t *pk, uint8_t *sk)
{
	uint8_t zeta[SEEDBYTES];
	uint8_t seedbuf[3 * SEEDBYTES];
	uint8_t tr[SEEDBYTES];
	const uint8_t *xi_1, *xi_2, *key;
	poly mat_8way;
	poly s1, s2, t1, t0;
	poly temp_result, temp_s1;

	randombytes(zeta, SEEDBYTES);
	randombytes(seedbuf, SEEDBYTES);
	shake256(seedbuf, 3 * SEEDBYTES, seedbuf, SEEDBYTES);

	xi_1 = seedbuf;
	xi_2 = seedbuf + SEEDBYTES;
	key = seedbuf + 2 * SEEDBYTES;

	poly_uniform(&mat_8way, zeta, 0);
	poly_uniform_eta(&s1, xi_1, 0);
	poly_uniform_eta(&s2, xi_2, 0);

	ntt_avx_asm(&temp_s1, &s1);
	poly_base_mul_avx_asm(&temp_result, &temp_s1, &mat_8way);
	invntt_avx_asm(&t1, &temp_result);

	poly_caddq(&t1);
	poly_add(&t1, &t1, &s2);
	poly_caddq(&t1);
	poly_power2round(&t1, &t0, &t1);
	pack_pk(pk, zeta, &t1);

	shake256(tr, SEEDBYTES, pk, CRYPTO_PUBLICKEYBYTES);

	pack_sk(sk, zeta, tr, key, &t0, &s1, &s2);

	return 0;
}

int crypto_sign_signature(uint8_t *sig,
						  size_t *siglen,
						  const uint8_t *m,
						  size_t mlen,
						  const uint8_t *sk)
{
	unsigned int n;
	uint8_t seedbuf[3 * SEEDBYTES + 2 * CRHBYTES];
	uint8_t *zeta, *tr, *key, *mu, *rho;
	uint16_t nonce = 0;
	poly s1, y, z, t0, s2, w1, w0, h, cp, mat_8way,temp_result;
	keccak_state state;
	uint8_t neg_start = 0;
	uint16_t sparse[TAU] = {0};
	zeta = seedbuf;
	tr = zeta + SEEDBYTES;
	key = tr + SEEDBYTES;
	mu = key + SEEDBYTES;
	rho = mu + CRHBYTES;

	unpack_sk(zeta, tr, key, &t0, &s1, &s2, sk);
	
	shake256_init(&state);
	shake256_absorb(&state, tr, SEEDBYTES);
	shake256_absorb(&state, m, mlen);
	shake256_finalize(&state);
	shake256_squeeze(mu, CRHBYTES, &state);

#ifdef NIMS_RANDOMIZED_SIGNING
	randombytes(rho, CRHBYTES);
#else
	shake256(rho, CRHBYTES, key, SEEDBYTES + CRHBYTES);
#endif
	poly_uniform(&mat_8way,zeta,0);

rej:
	poly_uniform_gamma1(&y, rho, nonce++);
	z = y;
	ntt_avx_asm(&z, &z);
	poly_base_mul_avx_asm(&temp_result,&z, &mat_8way);
	invntt_avx_asm(&w1,&temp_result);

	poly_caddq(&w1);
	poly_decompose(&w1, &w0, &w1);
	
	polyw1_pack(sig, &w1);

	shake256_init(&state);
	shake256_absorb(&state, mu, CRHBYTES);
	shake256_absorb(&state, sig, POLYW1_PACKEDBYTES);
	shake256_finalize(&state);
	shake256_squeeze(sig, SEEDBYTES, &state);

	poly_challenge(&cp, sig);
	neg_start = convToIdx(sparse,TAU,cp.coeffs,N);

	poly_mult_add(&h, &s2, sparse, neg_start);
	poly_sub(&w0, &w0, &h);
  	poly_reduce(&w0);
	if(poly_chknorm(&w0, GAMMA2 - BETA))	goto rej;

	poly_mult_add(&z, &s1, sparse, neg_start);
	poly_add(&z, &z, &y);
  	poly_reduce(&z);

	if (poly_chknorm(&z, GAMMA1 - BETA))	goto rej;
	poly_mult_add(&h, &t0, sparse, neg_start);
  	poly_reduce(&h);

	if (poly_chknorm(&h, GAMMA2))	goto rej;

	poly_add(&w0, &w0, &h);
	n = poly_make_hint(&h, &w0, &w1);
	if (n > OMEGA)	goto rej;

	pack_sig(sig, sig, &z, &h);
	*siglen = CRYPTO_BYTES;
	return 0;
}

int crypto_sign(uint8_t *sm,
				size_t *smlen,
				const uint8_t *m,
				size_t mlen,
				const uint8_t *sk)
{
	size_t i;

	for (i = 0; i < mlen; i++) sm[CRYPTO_BYTES + mlen - 1 - i] = m[mlen - 1 - i];
	crypto_sign_signature(sm, smlen, sm + CRYPTO_BYTES, mlen, sk);
	*smlen += mlen;
	return 0;
}

int crypto_sign_verify(const uint8_t *sig,
					   size_t siglen,
					   const uint8_t *m,
					   size_t mlen,
					   const uint8_t *pk)
{
	unsigned int i;
	uint8_t buf[POLYW1_PACKEDBYTES];
	uint8_t zeta[SEEDBYTES];
	uint8_t mu[CRHBYTES];
	uint8_t c[SEEDBYTES];
	uint8_t c2[SEEDBYTES];

	poly cp, z, t1, w1, h;
	poly mat_8way, temp_w1, temp_t1, temp_result;
	keccak_state state;

	if (siglen != CRYPTO_BYTES)	return -1;
	unpack_pk(zeta, &t1, pk);

	if (unpack_sig(c, &z, &h, sig))	return -1;

	if (poly_chknorm(&z, GAMMA1 - BETA))	return -1;

	shake256(mu, SEEDBYTES, pk, CRYPTO_PUBLICKEYBYTES);
	shake256_init(&state);
	shake256_absorb(&state, mu, SEEDBYTES);
	shake256_absorb(&state, m, mlen);
	shake256_finalize(&state);
	shake256_squeeze(mu, CRHBYTES, &state);
	poly_challenge(&cp, c);
	poly_uniform(&mat_8way, zeta, 0);

	ntt_avx_asm(&z, &z);
	poly_base_mul_avx_asm(&temp_w1, &z, &mat_8way);

	ntt_avx_asm(&cp, &cp);
	poly_shiftl(&t1);

	ntt_avx_asm(&t1, &t1);
	poly_base_mul_avx_asm(&temp_t1, &cp, &t1);

	poly_sub(&temp_result, &temp_w1, &temp_t1);
	invntt_avx_asm(&w1, &temp_result);

	poly_caddq(&w1);
	poly_use_hint(&w1, &w1, &h);
	polyw1_pack(buf, &w1);
	shake256_init(&state);
	shake256_absorb(&state, mu, CRHBYTES);
	shake256_absorb(&state, buf, POLYW1_PACKEDBYTES);
	shake256_finalize(&state);
	shake256_squeeze(c2, SEEDBYTES, &state);

	for (i = 0; i < SEEDBYTES; i++) if (c[i] != c2[i]) return -1;

	return 0;
}

int crypto_sign_open(uint8_t *m,
					 size_t *mlen,
					 const uint8_t *sm,
					 size_t smlen,
					 const uint8_t *pk)
{
	size_t i;

	if (smlen < CRYPTO_BYTES)
		goto badsig;

	*mlen = smlen - CRYPTO_BYTES;
	if (crypto_sign_verify(sm, CRYPTO_BYTES, sm + CRYPTO_BYTES, *mlen, pk))
		goto badsig;
	else
	{
		for (i = 0; i < *mlen; i++)
			m[i] = sm[CRYPTO_BYTES + i];
		return 0;
	}

badsig:
	/* Signature verification failed */
	*mlen = -1;
	for (i = 0; i < smlen; i++)
		m[i] = 0;

	return -1;
}
