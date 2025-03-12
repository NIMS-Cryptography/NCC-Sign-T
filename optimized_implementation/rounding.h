#ifndef ROUNDING_H
#define ROUNDING_H

#include <stdint.h>
#include "params.h"

void power2round_avx(__m256i *a1, __m256i *a0, const __m256i *a);


unsigned int make_hint(int32_t a0, int32_t a1);

void decompose_avx(__m256i *a1, __m256i *a0, const __m256i *a);
void use_hint_avx(__m256i *b, const __m256i *a, const __m256i *restrict hint);

#endif
