#include <stdint.h>
#include "symmetric.h"
#include "aes.h"

void NIMS_aes256ctr_init(aes256ctr_ctx *state, const uint8_t key[32], uint64_t nonce)
{
  aes256ctr_init(state, key, nonce);
}
