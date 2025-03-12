#include "rounding.h"

#include <stdint.h>

#include "params.h"
int32_t power2round(int32_t *a0, int32_t a) {
	int32_t a1;
	a1 = (a + (1 << (D - 1)) - 1) >> D;
	*a0 = a - (a1 << D);
	return a1;
}

int32_t decompose_origin(int32_t *a0, int32_t a) {
	int32_t a1;

	a = a % Q;
	if (a < 0)
		a += Q;
	*a0 = a % (2 * GAMMA2);
	if (*a0 > GAMMA2)
		*a0 -= (2 * GAMMA2);
	if (a - (*a0) == (Q - 1)) {
		a1 = 0;
		*a0 = *a0 - 1;
	} else {
		a1 = (a - (*a0)) / (2 * GAMMA2);
	}
	return a1;
}

int32_t decompose(int32_t *a0, int32_t a) {
#if N == 2304
	int32_t a1;
	a1 = (a + 127) >> 7;
	a1 = (a1 * 1022 + (1UL << 21)) >> 22;
	a1 &= 15;

#elif N == 1536
	int64_t a1;
	a1 = (a + 7) >> 3;
	a1 = ((a1 * 1047489) + (1UL << 34)) >> 35;
	a1 &= 31;

#elif N == 1152
	int64_t a1;
	a1 = (a + 1) >> 1;
	a1 = ((a1 * 4187849) + (1UL << 38)) >> 39;
	a1 &= 31;

#elif N == 2048
	int32_t a1;
	a1 = (a + 127) >> 7;
	a1 = ((a1 * 1025) + (1UL << 20)) >> 21;
	a1 &= 31;

#endif
	*a0 = a - a1 * 2 * GAMMA2;
	*a0 -= (((Q - 1) / 2 - *a0) >> 31) & Q;
	return a1;
}

unsigned int make_hint(int32_t a0, int32_t a1) {
	if (a0 > GAMMA2 || a0 < -GAMMA2 || (a0 == -GAMMA2 && a1 != 0)) return 1;

	return 0;
}

int32_t use_hint(int32_t a, unsigned int hint) {
	int32_t a0, a1;

	a1 = decompose(&a0, a);

	if (hint == 0) return a1;

#if N == 2304
	if (a0 > 0)
		return (a1 == 15) ? 0 : a1 + 1;
	else
		return (a1 == 0) ? 15 : a1 - 1;
#else
	if (a0 > 0)
		return (a1 == 31) ? 0 : a1 + 1;
	else
		return (a1 == 0) ? 31 : a1 - 1;
#endif
}
