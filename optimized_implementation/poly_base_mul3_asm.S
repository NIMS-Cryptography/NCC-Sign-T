#include "consts.h"

.section .text

.macro shuffle3 s0,s1,s2,r0, r1, r2
vpunpckldq  %ymm\s1,   %ymm\s0,  %ymm12
vpunpckhdq  %ymm\s1,   %ymm\s0,  %ymm13
vpunpckldq  %ymm\s1,   %ymm\s2,  %ymm14
vpunpckhdq  %ymm\s1,   %ymm\s2,  %ymm15

vpunpcklqdq %ymm14,   %ymm12,  %ymm\r0
vpunpckhqdq %ymm14,   %ymm12,  %ymm\r1
vpunpcklqdq %ymm15,   %ymm13,  %ymm\r2
.endm

   # rdi : c
   # rsi : a
   # rdx : b
   # rcx : zetas
.global poly_base_mul3_asm
poly_base_mul3_asm:
   # initial settng
   xor   %rax,  %rax
   vmovdqa      _8XQINV*4(%rcx),%ymm0
   vmovdqa      _8XQ*4(%rcx),%ymm1
   
   .p2align 5 # for Cache Hit
   __loop:
   vmovdqu  0*4(%rsi),%xmm4   # a:[ 0,  1,  2, ?]
   vmovdqu  3*4(%rsi),%xmm5   # a:[ 3,  4,  5, ?]
   vmovdqu  6*4(%rsi),%xmm6   # a:[ 6,  7,  8, ?]
   vmovdqu  9*4(%rsi),%xmm7   # a:[ 9, 10, 11, ?]
   vmovdqu  12*4(%rsi),%xmm8  # a:[12, 13, 14, ?]
   vmovdqu  15*4(%rsi),%xmm9  # a:[15, 16, 17, ?]
   vmovdqu  0*4(%rdx),%xmm10  # b:[ 0,  1,  2, ?]
   vmovdqu  3*4(%rdx),%xmm11  # b:[ 3,  4,  5, ?]
   vmovdqu  6*4(%rdx),%xmm12  # b:[ 6,  7,  8, ?]
   vmovdqu  9*4(%rdx),%xmm13  # b:[ 9, 10, 11, ?]
   vmovdqu  12*4(%rdx),%xmm14 # b:[12, 13, 14, ?]
   vmovdqu  15*4(%rdx),%xmm15 # b:[15, 16, 17, ?]

   # packing 6 coeff
   vinserti128  $1, %xmm7,  %ymm4,  %ymm4 # a:[ 0,  1,  2, ?] [ 9, 10, 11, ?]
   vinserti128  $1, %xmm8,  %ymm5,  %ymm5 # a:[ 3,  4,  5, ?] [12, 13, 14, ?]
   vinserti128  $1, %xmm9,  %ymm6,  %ymm6 # a:[ 6,  7,  8, ?] [15, 16, 17, ?]
   vinserti128  $1, %xmm13, %ymm10, %ymm7 # b:[ 0,  1,  2, ?] [ 9, 10, 11, ?]
   vinserti128  $1, %xmm14, %ymm11, %ymm8 # b:[ 3,  4,  5, ?] [12, 13, 14, ?]
   vinserti128  $1, %xmm15, %ymm12, %ymm9 # b:[ 6,  7,  8, ?] [15, 16, 17, ?]

   # a:ymm4 : [ 0,  3,  6, ?] [ 9, 12, 15, ?] b:ymm7 : [ 0,  3,  6, ?] [ 9, 12, 15, ?]
   # a:ymm5 : [ 1,  4,  7, ?] [10, 13, 16, ?] b:ymm8 : [ 1,  4,  7, ?] [10, 13, 16, ?]
   # a:ymm6 : [ 2,  5,  8, ?] [11, 14, 17, ?] b:ymm9 : [ 2,  5,  8, ?] [11, 14, 17, ?]
   shuffle3 4, 5, 6, 4, 5, 6
   shuffle3 7, 8, 9, 7, 8, 9

   # odd coeff -----------------------------------------------------------------------
   vpmuldq     %ymm6,%ymm8,%ymm10   # C[0] : A[2]*B[1]
   vpmuldq     %ymm5,%ymm9,%ymm3    # C[0] : A[1]*B[2]
   vpmuldq     %ymm6,%ymm9,%ymm11   # C[1] : A[2]* B[2]
   vpmuldq     %ymm6,%ymm7,%ymm13   # C[2] : A[2]*B[0]
   vpmuldq     %ymm5,%ymm8,%ymm14   # C[2] : A[1]*B[1]
   vpmuldq     %ymm4,%ymm9,%ymm15   # C[2] : A[0]*B[2]

   vpaddq      %ymm10,%ymm3, %ymm10 # C[0] : A[2]*B[1] + A[1]*B[2]
   vpaddq      %ymm13,%ymm14,%ymm14 # C[2] : A[2]*B[0] + A[1]*B[1]
   vpaddq      %ymm14,%ymm15,%ymm12 # C[2] : A[0]*B[2] + A[1]*B[1] + A[2]*B[0] 
   
   vpmuldq     %ymm0,%ymm10,%ymm13
   vpmuldq     %ymm0,%ymm11,%ymm14
   vpmuldq     %ymm0,%ymm12,%ymm15
   vpmuldq     %ymm1,%ymm13,%ymm13
   vpmuldq     %ymm1,%ymm14,%ymm14
   vpmuldq     %ymm1,%ymm15,%ymm15

   vpsubq      %ymm13,%ymm10,%ymm10
   vpsubq      %ymm14,%ymm11,%ymm11
   vpsubq      %ymm15,%ymm12,%ymm12

   vmovdqu     (_ZETAS+512)*4(%rcx),%ymm2   
   vpsrlq      $32,%ymm10,%ymm10
   vpmuldq     %ymm2,%ymm10, %ymm13 # C[0] : C[0]*ZETA 
   vpmuldq     %ymm4,%ymm7,  %ymm15 # C[0] : A[0]*B[0]
   vpaddq      %ymm13,%ymm15,%ymm10 # C[0] : C[0]*ZETA + A[0]*B[0]
   vpmuldq     %ymm0,%ymm10,%ymm13
   vpmuldq     %ymm1,%ymm13,%ymm13
   vpsubq      %ymm13,%ymm10,%ymm10  

   vpsrlq      $32,%ymm11,%ymm11
   vpmuldq     %ymm2,%ymm11, %ymm3  # C[1] : C[1]*ZETA
   vpmuldq     %ymm4,%ymm8,  %ymm14 # C[1] : A[0]*B[1]
   vpmuldq     %ymm5,%ymm7,  %ymm15 # C[1] : A[1]*B[0]
   vpaddq      %ymm3,%ymm14, %ymm14 # C[1] : C[1]*ZETA + A[0]*B[1]
   vpaddq      %ymm15,%ymm14,%ymm11 # C[1] : C[1]*ZETA + A[0]*B[1] + A[1]*B[0]
   vpmuldq     %ymm0,%ymm11,%ymm13
   vpmuldq     %ymm1,%ymm13,%ymm13
   vpsubq      %ymm13,%ymm11,%ymm11

   vpsrlq      $32,%ymm4,%ymm4
   vpsrlq      $32,%ymm5,%ymm5
   vmovshdup   %ymm6,%ymm6
   vpsrlq      $32,%ymm7,%ymm7
   vpsrlq      $32,%ymm8,%ymm8
   vmovshdup   %ymm9,%ymm9

   vpmuldq     %ymm6,%ymm8,%ymm2    # C[0] : A[2]*B[1] 
   vpmuldq     %ymm5,%ymm9,%ymm3    # C[0] : A[1]*B[2]
   vpmuldq     %ymm5,%ymm8,%ymm14   # C[2] : A[1]*B[1]
   vpmuldq     %ymm6,%ymm7,%ymm13   # C[2] : A[2]*B[0]
   vpmuldq     %ymm4,%ymm9,%ymm15   # C[2] : A[0]*B[2]
   vpmuldq     %ymm6,%ymm9,%ymm9    # C[1] : A[2]* B[2]

   vpaddq      %ymm2,%ymm3,%ymm2    # C[0] : A[2]*B[1] + A[1]*B[2]
   vpaddq      %ymm13,%ymm14,%ymm14 # C[2] : A[2]*B[0] + A[1]*B[1]
   vpaddq      %ymm14,%ymm15,%ymm6  # C[2] : A[0]*B[2] + A[1]*B[1] + A[2]*B[0] 

   vpmuldq     %ymm0,%ymm2,%ymm13
   vpmuldq     %ymm0,%ymm9,%ymm14
   vpmuldq     %ymm0,%ymm6,%ymm15
   vpmuldq     %ymm1,%ymm13,%ymm13
   vpmuldq     %ymm1,%ymm14,%ymm14
   vpmuldq     %ymm1,%ymm15,%ymm15

   vpsubq      %ymm13,%ymm2,%ymm2
   vpsubq      %ymm14,%ymm9,%ymm9
   vpsubq      %ymm15,%ymm6,%ymm6

   vmovdqu     (_ZETAS+512)*4(%rcx),%ymm3
   vpsrlq      $32,%ymm3,%ymm3
   vpsrlq      $32,%ymm2,%ymm2
   vpmuldq     %ymm3,%ymm2, %ymm13  # C[0] : C[0]*ZETA 
   vpmuldq     %ymm4,%ymm7,  %ymm15 # C[0] : A[0]*B[0]
   vpaddq      %ymm13,%ymm15,%ymm2  # C[0] : C[0]*ZETA  + A[0]*B[0]
   vpmuldq     %ymm0,%ymm2,%ymm13
   vpmuldq     %ymm1,%ymm13,%ymm13
   vpsubq      %ymm13,%ymm2,%ymm2
   
   vpsrlq      $32,%ymm9,%ymm9
   vpmuldq     %ymm3,%ymm9, %ymm3   # C[1] : C[1]*ZETA
   vpmuldq     %ymm4,%ymm8,  %ymm14 # C[1] : A[0]*B[1]
   vpmuldq     %ymm5,%ymm7,  %ymm15 # C[1] : A[1]*B[0]
   vpaddq      %ymm3,%ymm14, %ymm14 # C[1] : C[1]*ZETA + A[0]*B[1]
   vpaddq      %ymm15,%ymm14,%ymm9  # C[1] : C[1]*ZETA + A[0]*B[1] + A[1]*B[0]
   vpmuldq     %ymm0,%ymm9,%ymm13
   vpmuldq     %ymm1,%ymm13,%ymm13
   vpsubq      %ymm13,%ymm9,%ymm9
     
   vpsrlq      $32,%ymm10,%ymm10
   vpsrlq      $32,%ymm11,%ymm11
   vmovshdup   %ymm12,%ymm12
   vpblendd    $0xAA,%ymm2,%ymm10,%ymm2
   vpblendd    $0xAA,%ymm9,%ymm11,%ymm9
   vpblendd    $0xAA,%ymm6,%ymm12,%ymm6
   # c:ymm4 : [ 0,  1,  2, ?] [ 9, 10, 11, ?] 
   # c:ymm5 : [ 3,  4,  5, ?] [12, 13, 14, ?] 
   # c:ymm6 : [ 6,  7,  8, ?] [15, 16, 17, ?] 
   shuffle3 2, 9, 6, 4, 5, 6

   # unpacking 6 coeff
   vextracti128 $1, %ymm4, %xmm7 # c:[ 0,  1,  2, ?] [ 9, 10, 11, ?]
   vextracti128 $1, %ymm5, %xmm8 # c:[ 3,  4,  5, ?] [12, 13, 14, ?]
   vextracti128 $1, %ymm6, %xmm9 # c:[ 6,  7,  8, ?] [15, 16, 17, ?]

   vmovdqu     6*4(%rsi),%xmm14
   vpblendd    $0xF8,%xmm14,%xmm6,%xmm6
   vmovdqu     15*4(%rsi),%xmm15
   vpblendd    $0xF8,%xmm15,%xmm9,%xmm9
   vmovdqu     %xmm4,   0*4(%rdi) # c:[ 0,  1,  2, ?] 
   vmovdqu     %xmm5,   3*4(%rdi) # c:[ 3,  4,  5, ?] 
   vmovdqu     %xmm6,   6*4(%rdi) # c:[ 6,  7,  8, ?]
   vmovdqu     %xmm7,   9*4(%rdi) # c:[ 9, 10, 11, ?]
   vmovdqu     %xmm8,  12*4(%rdi) # c:[12, 13, 14, ?]
   vmovdqu     %xmm9,  15*4(%rdi) # c:[15, 16, 17, ?]

   add     $72,   %rsi
   add     $72,   %rdi
   add     $72,   %rdx
   add     $32,   %rcx
   add     $1,    %rax
   cmp     $85,   %rax
   jb      __loop
   
   # last 6 coeff ********************************************************************
   vmovdqu     0*4(%rsi),%xmm4  # a:[ 0,  1,  2, ?]
   vmovdqu     3*4(%rsi),%xmm5  # a:[ 3,  4,  5, ?]
   vmovdqu     6*4(%rsi),%xmm6  # a:[ 6,  7,  8, ?]
   vmovdqu     0*4(%rdx),%xmm10 # b:[ 0,  1,  2, ?]
   vmovdqu     3*4(%rdx),%xmm11 # b:[ 3,  4,  5, ?]
   vmovdqu     6*4(%rdx),%xmm12 # b:[ 6,  7,  8, ?]

   # packing 6 coeff
   vinserti128  $1, %xmm7,  %ymm4,  %ymm4 # a:[ 0,  1,  2, ?] [ 9, 10, 11, ?]
   vinserti128  $1, %xmm8,  %ymm5,  %ymm5 # a:[ 3,  4,  5, ?] [12, 13, 14, ?]
   vinserti128  $1, %xmm9,  %ymm6,  %ymm6 # a:[ 6,  7,  8, ?] [15, 16, 17, ?]
   vinserti128  $1, %xmm13, %ymm10, %ymm7 # b:[ 0,  1,  2, ?] [ 9, 10, 11, ?]
   vinserti128  $1, %xmm14, %ymm11, %ymm8 # b:[ 3,  4,  5, ?] [12, 13, 14, ?]
   vinserti128  $1, %xmm15, %ymm12, %ymm9 # b:[ 6,  7,  8, ?] [15, 16, 17, ?]

   # a:ymm4 : [ 0,  3,  6, ?] [ 9, 12, 15, ?] b:ymm7 : [ 0,  3,  6, ?] [ 9, 12, 15, ?]
   # a:ymm5 : [ 1,  4,  7, ?] [10, 13, 16, ?] b:ymm8 : [ 1,  4,  7, ?] [10, 13, 16, ?]
   # a:ymm6 : [ 2,  5,  8, ?] [11, 14, 17, ?] b:ymm9 : [ 2,  5,  8, ?] [11, 14, 17, ?]
   shuffle3 4, 5, 6, 4, 5, 6
   shuffle3 7, 8, 9, 7, 8, 9

   # odd coeff -----------------------------------------------------------------------
   vpmuldq     %ymm6,%ymm8,%ymm10   # C[0] : A[2]*B[1]
   vpmuldq     %ymm5,%ymm9,%ymm3    # C[0] : A[1]*B[2]
   vpmuldq     %ymm6,%ymm9,%ymm11   # C[1] : A[2]* B[2]
   vpmuldq     %ymm6,%ymm7,%ymm13   # C[2] : A[2]*B[0]
   vpmuldq     %ymm5,%ymm8,%ymm14   # C[2] : A[1]*B[1]
   vpmuldq     %ymm4,%ymm9,%ymm15   # C[2] : A[0]*B[2]

   vpaddq      %ymm10,%ymm3, %ymm10 # C[0] : A[2]*B[1] + A[1]*B[2]
   vpaddq      %ymm13,%ymm14,%ymm14 # C[2] : A[2]*B[0] + A[1]*B[1]
   vpaddq      %ymm14,%ymm15,%ymm12 # C[2] : A[0]*B[2] + A[1]*B[1] + A[2]*B[0] 
   
   vpmuldq     %ymm0,%ymm10,%ymm13
   vpmuldq     %ymm0,%ymm11,%ymm14
   vpmuldq     %ymm0,%ymm12,%ymm15
   vpmuldq     %ymm1,%ymm13,%ymm13
   vpmuldq     %ymm1,%ymm14,%ymm14
   vpmuldq     %ymm1,%ymm15,%ymm15

   vpsubq      %ymm13,%ymm10,%ymm10
   vpsubq      %ymm14,%ymm11,%ymm11
   vpsubq      %ymm15,%ymm12,%ymm12

   vmovdqu     (_ZETAS+512)*4(%rcx),%ymm2   
   vpsrlq      $32,%ymm10,%ymm10
   vpmuldq     %ymm2,%ymm10, %ymm13 # C[0] : C[0]*ZETA 
   vpmuldq     %ymm4,%ymm7,  %ymm15 # C[0] : A[0]*B[0]
   vpaddq      %ymm13,%ymm15,%ymm10 # C[0] : C[0]*ZETA  + A[0]*B[0]
   vpmuldq     %ymm0,%ymm10,%ymm13
   vpmuldq     %ymm1,%ymm13,%ymm13
   vpsubq      %ymm13,%ymm10,%ymm10  

   vpsrlq      $32,%ymm11,%ymm11
   vpmuldq     %ymm2,%ymm11, %ymm3  # C[1] : C[1]*ZETA
   vpmuldq     %ymm4,%ymm8,  %ymm14 # C[1] : A[0]*B[1]
   vpmuldq     %ymm5,%ymm7,  %ymm15 # C[1] : A[1]*B[0]
   vpaddq      %ymm3,%ymm14, %ymm14 # C[1] : C[1]*ZETA + A[0]*B[1]
   vpaddq      %ymm15,%ymm14,%ymm11 # C[1] : C[1]*ZETA + A[0]*B[1] + A[1]*B[0]
   vpmuldq     %ymm0,%ymm11,%ymm13
   vpmuldq     %ymm1,%ymm13,%ymm13
   vpsubq      %ymm13,%ymm11,%ymm11

   vpsrlq      $32,%ymm4,%ymm4
   vpsrlq      $32,%ymm5,%ymm5
   vmovshdup   %ymm6,%ymm6
   vpsrlq      $32,%ymm7,%ymm7
   vpsrlq      $32,%ymm8,%ymm8
   vmovshdup   %ymm9,%ymm9

   vpmuldq     %ymm6,%ymm8,%ymm2    # C[0] : A[2]*B[1] 
   vpmuldq     %ymm5,%ymm9,%ymm3    # C[0] : A[1]*B[2]
   vpmuldq     %ymm5,%ymm8,%ymm14   # C[2] : A[1]*B[1]
   vpmuldq     %ymm6,%ymm7,%ymm13   # C[2] : A[2]*B[0]
   vpmuldq     %ymm4,%ymm9,%ymm15   # C[2] : A[0]*B[2]
   vpmuldq     %ymm6,%ymm9,%ymm9    # C[1] : A[2]* B[2]

   vpaddq      %ymm2,%ymm3,%ymm2    # C[0] : A[2]*B[1] + A[1]*B[2]
   vpaddq      %ymm13,%ymm14,%ymm14 # C[2] : A[2]*B[0] + A[1]*B[1]
   vpaddq      %ymm14,%ymm15,%ymm6  # C[2] : A[0]*B[2] + A[1]*B[1] + A[2]*B[0] 

   vpmuldq     %ymm0,%ymm2,%ymm13
   vpmuldq     %ymm0,%ymm9,%ymm14
   vpmuldq     %ymm0,%ymm6,%ymm15
   vpmuldq     %ymm1,%ymm13,%ymm13
   vpmuldq     %ymm1,%ymm14,%ymm14
   vpmuldq     %ymm1,%ymm15,%ymm15

   vpsubq      %ymm13,%ymm2,%ymm2
   vpsubq      %ymm14,%ymm9,%ymm9
   vpsubq      %ymm15,%ymm6,%ymm6

   vmovdqu     (_ZETAS+512)*4(%rcx),%ymm3
   vpsrlq      $32,%ymm3,%ymm3
   vpsrlq      $32,%ymm2,%ymm2
   vpmuldq     %ymm3,%ymm2, %ymm13  # C[0] : C[0]*ZETA 
   vpmuldq     %ymm4,%ymm7,  %ymm15 # C[0] : A[0]*B[0]
   vpaddq      %ymm13,%ymm15,%ymm2  # C[0] : C[0]*ZETA  + A[0]*B[0]
   vpmuldq     %ymm0,%ymm2,%ymm13
   vpmuldq     %ymm1,%ymm13,%ymm13
   vpsubq      %ymm13,%ymm2,%ymm2
   
   vpsrlq      $32,%ymm9,%ymm9
   vpmuldq     %ymm3,%ymm9, %ymm3   # C[1] : C[1]*ZETA
   vpmuldq     %ymm4,%ymm8,  %ymm14 # C[1] : A[0]*B[1]
   vpmuldq     %ymm5,%ymm7,  %ymm15 # C[1] : A[1]*B[0]
   vpaddq      %ymm3,%ymm14, %ymm14 # C[1] : C[1]*ZETA + A[0]*B[1]
   vpaddq      %ymm15,%ymm14,%ymm9  # C[1] : C[1]*ZETA + A[0]*B[1] + A[1]*B[0]
   vpmuldq     %ymm0,%ymm9,%ymm13
   vpmuldq     %ymm1,%ymm13,%ymm13
   vpsubq      %ymm13,%ymm9,%ymm9
     
   vpsrlq      $32,%ymm10,%ymm10
   vpsrlq      $32,%ymm11,%ymm11
   vmovshdup   %ymm12,%ymm12
   vpblendd    $0xAA,%ymm2,%ymm10,%ymm2
   vpblendd    $0xAA,%ymm9,%ymm11,%ymm9
   vpblendd    $0xAA,%ymm6,%ymm12,%ymm6
   # c:ymm4 : [ 0,  1,  2, ?] [ 9, 10, 11, ?] 
   # c:ymm5 : [ 3,  4,  5, ?] [12, 13, 14, ?] 
   # c:ymm6 : [ 6,  7,  8, ?] [15, 16, 17, ?] 
   shuffle3 2, 9, 6, 4, 5, 6

   # unpacking 6 coeff
   vextracti128 $1, %ymm4, %xmm7 # c:[ 0,  1,  2, ?] [ 9, 10, 11, ?]
   vextracti128 $1, %ymm5, %xmm8 # c:[ 3,  4,  5, ?] [12, 13, 14, ?]
   vextracti128 $1, %ymm6, %xmm9 # c:[ 6,  7,  8, ?] [15, 16, 17, ?]

   vmovdqu     %xmm4,   0*4(%rdi) # c:[ 0,  1,  2, ?] 
   vmovdqu     %xmm5,   3*4(%rdi) # c:[ 3,  4,  5, ?] 

   ret