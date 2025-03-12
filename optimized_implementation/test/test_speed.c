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
  poly mat;
  poly eta, gamma, a, b, c;
  uint8_t seed[3 * SEEDBYTES];

  t_overhead = 19;

  for (i = 0; i < NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_uniform(&mat, seed, 0);
  }
  print_results("nims", t, NTESTS);

  for (i = 0; i < NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_uniform_eta(&eta, seed, 0);
  }
  print_results("poly_uniform_eta:", t, NTESTS);

  for (i = 0; i < NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_uniform_gamma1(&gamma, seed, 0);
  }
  print_results("poly_uniform_gamma1:", t, NTESTS);

  for (i = 0; i < NTESTS; i++)
  {
    t[i] = cpucycles();
    ntt_avx_asm(&b, &a);
  }
  print_results("poly_ntt:", t, NTESTS);

  for (i = 0; i < NTESTS; i++)
  {
    t[i] = cpucycles();
    invntt_avx_asm(&c, &b);
  }
  print_results("poly_invntt:", t, NTESTS);

  for (i = 0; i < NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_base_mul_avx_asm(&c, &b, &a);
  }
  print_results("poly_base_mul:", t, NTESTS);

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
