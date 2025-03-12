#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "../params.h"
#include "../randombytes.h"
#include "../poly.h"
#include "../packing.h"

#define NTESTS 10000
void polyw1_unpack(poly *r, const uint8_t *a)
{
    unsigned int i;
#if N == 1024
    for (i = 0; i < 128; ++i)
    { // r 5 a 8
        r->coeffs[8 * i + 0] = (a[5 * i + 0] >> 0) & 0x1F;
        r->coeffs[8 * i + 1] = ((a[5 * i + 0] >> 5) | ((uint32_t)a[5 * i + 1] << 3)) & 0x1F;
        r->coeffs[8 * i + 2] = (a[5 * i + 1] >> 2) & 0x1F;
        r->coeffs[8 * i + 3] = ((a[5 * i + 1] >> 7) | ((uint32_t)a[5 * i + 2] << 1)) & 0x1F;
        r->coeffs[8 * i + 4] = ((a[5 * i + 2] >> 4) | ((uint32_t)a[5 * i + 3] << 4)) & 0x1F;
        r->coeffs[8 * i + 5] = (a[5 * i + 3] >> 1) & 0x1F;
        r->coeffs[8 * i + 6] = ((a[5 * i + 3] >> 6) | ((uint32_t)a[5 * i + 4] << 2)) & 0x1F;
        r->coeffs[8 * i + 7] = (a[5 * i + 4] >> 3) & 0x1F;
    }

#elif N == 1458
    for (i = 0; i < 182; ++i)
    { // r 5 a 8
        r->coeffs[8 * i + 0] = (a[5 * i + 0] >> 0) & 0x1F;
        r->coeffs[8 * i + 1] = ((a[5 * i + 0] >> 5) | ((uint32_t)a[5 * i + 1] << 3)) & 0x1F;
        r->coeffs[8 * i + 2] = (a[5 * i + 1] >> 2) & 0x1F;
        r->coeffs[8 * i + 3] = ((a[5 * i + 1] >> 7) | ((uint32_t)a[5 * i + 2] << 1)) & 0x1F;
        r->coeffs[8 * i + 4] = ((a[5 * i + 2] >> 4) | ((uint32_t)a[5 * i + 3] << 4)) & 0x1F;
        r->coeffs[8 * i + 5] = (a[5 * i + 3] >> 1) & 0x1F;
        r->coeffs[8 * i + 6] = ((a[5 * i + 3] >> 6) | ((uint32_t)a[5 * i + 4] << 2)) & 0x1F;
        r->coeffs[8 * i + 7] = (a[5 * i + 4] >> 3) & 0x1F;
    }
    r->coeffs[1456] = (a[910] >> 0) & 0x1F;
    r->coeffs[1457] = ((a[910] >> 5) | ((uint32_t)a[911] << 3)) & 0x1F;

#elif N == 1944
    for (i = 0; i < 972; ++i)
    { // r 4 a 8
        r->coeffs[2 * i + 0] = (a[i] >> 0) & 0xF;
        r->coeffs[2 * i + 1] = (a[i] >> 4) & 0xF;
    }

#endif
}

int main(void)
{
    unsigned int i, j;
    uint8_t seed[CRHBYTES];
    uint16_t nonce = 0;
    poly a, b, c, d, e, f, g, h;
    poly tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
    uint8_t p_eta[POLYETA_PACKEDBYTES];
    uint8_t p_t1[POLYT1_PACKEDBYTES];
    uint8_t p_t0[POLYT0_PACKEDBYTES];
    uint8_t p_z[POLYZ_PACKEDBYTES];
    uint8_t p_w1[POLYW1_PACKEDBYTES];
    //	randombytes(seed, sizeof(seed));
    for (int iter = 0; iter < CRHBYTES; iter++)
    {
        seed[iter] = iter;
    }

    for (j = 0; j < 100000; j++)
    {
        poly_uniform(&a, seed, nonce++);
        poly_uniform(&b, seed, nonce++);
        poly_uniform(&c, seed, nonce++);
        poly_uniform(&d, seed, nonce++);
        poly_uniform(&e, seed, nonce++);
        poly_uniform(&f, seed, nonce++);
        poly_uniform(&g, seed, nonce++);
        poly_uniform(&h, seed, nonce++);

        for (i = 0; i < N; i++)
        {
            a.coeffs[i] &= (1 << 23) - 1;
        }
        for (i = 0; i < N; i++)
        {
            b.coeffs[i] &= (1 << 23) - 1;
        }
        for (i = 0; i < N; i++)
        {
            c.coeffs[i] &= (1 << 23) - 1;
        }
        for (i = 0; i < N; i++)
        {
            d.coeffs[i] &= (1 << 23) - 1;
        }
        for (i = 0; i < N; i++)
        {
            e.coeffs[i] &= (1 << 23) - 1;
        }
        for (i = 0; i < N; i++)
        {
            f.coeffs[i] &= (1 << 23) - 1;
        }
        for (i = 0; i < N; i++)
        {
            g.coeffs[i] &= (1 << 23) - 1;
        }
        for (i = 0; i < N; i++)
        {
            h.coeffs[i] &= (1 << 23) - 1;
        }
        // poly_mul(&tmp1, &a, &b);
        // poly_add(&d, &c, &tmp1);
        // poly_mul(&g, &e, &d);
        // poly_mul(&tmp2, &e, &tmp1);
        // poly_mul(&tmp3, &e, &c);
        // poly_add(&h, &tmp2, &tmp3);
        // 교환법칙
        poly_mul(&tmp1, &a, &b);
        poly_mul(&tmp2, &tmp1, &c);

        poly_mul(&tmp3, &a, &c);
        poly_mul(&tmp4, &tmp3, &b);

        // poly_add(&d, &tmp1, &tmp2);

        // poly_add(&tmp3, &b, &c);
        // poly_mul(&e, &a, &tmp3);

        for (int iter = 0; iter < N; iter++)
        {
            printf("%08x  %08x\n", tmp2.coeffs[iter] & 0x7FFFFF, tmp4.coeffs[iter] & 0x7FFFFF);
        }

        // 분배법칙
        poly_mul(&tmp1, &a, &b);
        poly_mul(&tmp2, &a, &c);
        poly_add(&d, &tmp1, &tmp2);

        poly_add(&tmp3, &b, &c);
        poly_mul(&e, &a, &tmp3);

        // for(int iter=0;iter<N;iter++){
        //     printf("%08x  %08x\n", d.coeffs[iter] &0x7FFFFF, e.coeffs[iter] & 0x7FFFFF);

        // }

        poly_add(&d, &c, &a);

        poly_mul(&g, &e, &d);
        poly_mul(&tmp2, &e, &tmp1);
        poly_mul(&tmp3, &e, &c);
        poly_add(&h, &tmp2, &tmp3);

        ///        poly_mul(&c,&a,&b);
        //       poly_mul_schoolbook(&c_sc,&a,&b);
        //     for(i=0;i<N;i++){
        //         if(g.coeffs[i]!=h.coeffs[i]) {
        //             printf("다름다름");
        //             exit(0);
        //         }
        //     }
        //     printf(" 다 같다\n");

        // }

        printf("asdfasdf %d\n", N);
        printf("polyeta.......");
        // eta
        for (j = 0; j < 100000; j++)
        {
            poly_uniform_eta(&a, seed, nonce++);
            polyeta_pack(p_eta, &a);
            polyeta_unpack(&b, p_eta);
            for (i = 0; i < N; i++)
            {
                if (a.coeffs[i] != b.coeffs[i])
                {
                    printf("\npack_eta error\n");
                    return -1;
                }
            }
        }

        printf("good~\n");
        // t1
        printf("polyt1.......");

        for (j = 0; j < 100000; j++)
        {
            poly_uniform(&a, seed, nonce++);
            for (i = 0; i < N; i++)
            {
                a.coeffs[i] &= ((1 << (23 - D)) - 1);
            }
            polyt1_pack(p_t1, &a);
            polyt1_unpack(&b, p_t1);
            for (i = 0; i < N; i++)
            {
                if (a.coeffs[i] != b.coeffs[i])
                {
                    printf("\n%d , pack_t1 error\n", i);
                    printf("a[%d]:%x\n", i, a.coeffs[i]);
                    printf("b[%d]:%x\n", i, b.coeffs[i]);
                    return -1;
                }
            }
        }

        printf("good~\n");
        // t0
        printf("polyt0.......");

        for (j = 0; j < 100000; j++)
        {
            poly_uniform(&a, seed, nonce++);
            for (i = 0; i < N; i++)
            {
                a.coeffs[i] &= ((1 << D) - 1);
                a.coeffs[i] = (1 << (D - 1)) - a.coeffs[i];
            }
            polyt0_pack(p_t0, &a);
            polyt0_unpack(&b, p_t0);
            for (i = 0; i < N; i++)
            {
                if (a.coeffs[i] != b.coeffs[i])
                {
                    printf("\n%d , pack_t0 error\n", i);
                    printf("a[%d]:%x\n", i, a.coeffs[i]);
                    printf("b[%d]:%x\n", i, b.coeffs[i]);
                    return -1;
                }
            }
        }

        printf("good~\n");
        // z
        printf("polyz.......");

        for (j = 0; j < 100000; j++)
        {
            poly_uniform(&a, seed, nonce++);
            for (i = 0; i < N; i++)
            {
                a.coeffs[i] &= (GAMMA1 << 1) - 1;
                a.coeffs[i] = GAMMA1 - a.coeffs[i];
            }
            polyz_pack(p_z, &a);
            polyz_unpack(&b, p_z);
            for (i = 0; i < N; i++)
            {
                if (a.coeffs[i] != b.coeffs[i])
                {
                    printf("\n%d , pack_z error\n", i);
                    printf("a[%d]:%x\n", i, a.coeffs[i]);
                    printf("b[%d]:%x\n", i, b.coeffs[i]);
                    return -1;
                }
            }
        }
        printf("good~\n");
        // w1
        printf("polyw1.......");
#if N == 1944
#define X ((1 << 4) - 1)
#else
#define X ((1 << 5) - 1)
#endif
        for (j = 0; j < 100000; j++)
        {
            poly_uniform(&a, seed, nonce++);
            for (i = 0; i < N; i++)
            {
                a.coeffs[i] &= X;
            }
            polyw1_pack(p_w1, &a);
            polyw1_unpack(&b, p_w1);
            for (i = 0; i < N; i++)
            {
                if (a.coeffs[i] != b.coeffs[i])
                {
                    printf("\n%d , pack_w1 error\n", i);
                    printf("a[%d]:%x\n", i, a.coeffs[i]);
                    printf("b[%d]:%x\n", i, b.coeffs[i]);
                    return -1;
                }
            }
        }

        printf("good~\n");

        return 0;
    }
}