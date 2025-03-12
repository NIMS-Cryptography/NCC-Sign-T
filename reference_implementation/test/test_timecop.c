#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include "../randombytes.h"
#include "../sign.h"
#include "poison.h"

//! please ref to "https://www.post-apocalyptic-crypto.org/timecop/crypto_sign/index.html"
//! also one need to add "-g" to CFLAGS in Makefile
// valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./test/test_timecop1
#define MLEN 59

static void timcop_doit(void)
{
	size_t mlen, smlen;
	uint8_t m[MLEN] = {0};
	uint8_t sm[MLEN + CRYPTO_BYTES] = {0};
	uint8_t m2[MLEN + CRYPTO_BYTES] = {0};
	uint8_t pk[CRYPTO_PUBLICKEYBYTES] = {0};
	uint8_t sk[CRYPTO_SECRETKEYBYTES] = {0};

	randombytes(m, MLEN);

	poison(sk, CRYPTO_SECRETKEYBYTES);
	crypto_sign(sm, &smlen, m, mlen, sk);
	unpoison(pk, CRYPTO_PUBLICKEYBYTES);
	unpoison(sm, smlen);
	crypto_sign_open(m2, &mlen, sm, smlen, pk);
}

int main(void)
{
	timcop_doit();
	return 0;
}
