#include <stdint.h>
#include "../sign.h"
#include "../poly.h"
#include "../params.h"
#include "cpucycles.h"
#include "speed_print.h"
#include "stdio.h"

#define NTESTS 10000

uint64_t t[NTESTS];

uint64_t t_mul = 0, t_notmul = 0, t_overhead;
int main(void)
{
  unsigned int i;
  size_t smlen;
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];
  uint8_t sm[CRYPTO_BYTES + CRHBYTES];
  uint8_t seed[CRHBYTES];
  poly a;
  __m256i b;
  int32_t c[N_avx * 3];

#ifdef SHAKE_ASSEMBLY
  printf("shake assembly mode \n");
#endif

  t_overhead = 19;

  for (i = 0; i < NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_uniform_eta(&a, seed, 0);
  }
  print_results("poly_uniform_eta:", t, NTESTS);

  for (i = 0; i < NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_uniform_gamma1(&a, seed, 0);
  }
  print_results("poly_uniform_gamma1:", t, NTESTS);

  for (i = 0; i < NTESTS; i++)
  {
    t[i] = cpucycles();
    ntt_avx(&b, c);
  }
  print_results("poly_ntt:", t, NTESTS);

  for (i = 0; i < NTESTS; i++)
  {
    t[i] = cpucycles();
    invntt_tomont_avx(c, &b);
  }
  print_results("poly_invntt:", t, NTESTS);

  for (i = 0; i < NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_challenge(&a, seed);
  }
  print_results("poly_challenge:", t, NTESTS);

  for (i = 0; i < NTESTS; i++)
  {
    t[i] = cpucycles();
    crypto_sign_keypair(pk, sk);
  }
  print_results("Keypair:", t, NTESTS);

  for (i = 0; i < NTESTS; i++)
  {
    t[i] = cpucycles();
    crypto_sign(sm, &smlen, sm, CRHBYTES, sk);
  }
  print_results("Sign:", t, NTESTS);

  for (i = 0; i < NTESTS; i++)
  {
    t[i] = cpucycles();
    crypto_sign_verify(sm, CRYPTO_BYTES, sm, CRHBYTES, pk);
  }
  print_results("Verify:", t, NTESTS);

  return 0;
}
