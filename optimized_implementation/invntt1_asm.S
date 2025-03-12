#include "consts.h"

.section .text

.macro rdx2_butterfly l,h,zl0=1,zl1=1,zh0=2,zh1=2
   vpsubd      %ymm\h,%ymm\l,%ymm12
   vpaddd      %ymm\h,%ymm\l,%ymm\l
   vpmuldq     %ymm\zl0,%ymm12,%ymm13
   vmovshdup   %ymm12,%ymm\h
   vpmuldq     %ymm\zl1,%ymm\h,%ymm14
   vpmuldq     %ymm\zh0,%ymm12,%ymm12
   vpmuldq     %ymm\zh1,%ymm\h,%ymm\h
   vpmuldq     %ymm0,%ymm13,%ymm13
   vpmuldq     %ymm0,%ymm14,%ymm14
   vpsubd      %ymm13,%ymm12,%ymm12
   vpsubd      %ymm14,%ymm\h,%ymm\h
   vmovshdup   %ymm12,%ymm12
   vpblendd    $0xAA,%ymm\h,%ymm12,%ymm\h
.endm

.macro mont_reduce h,zl0=1,zl1=1,zh0=2,zh1=2
   vpmuldq     %ymm\zl0,%ymm\h,%ymm13
   vmovshdup   %ymm\h,%ymm12
   vpmuldq     %ymm\zl1,%ymm12,%ymm14
   vpmuldq     %ymm\zh0,%ymm\h,%ymm\h
   vpmuldq     %ymm\zh1,%ymm12,%ymm12
   vpmuldq     %ymm0,%ymm13,%ymm13
   vpmuldq     %ymm0,%ymm14,%ymm14
   vmovshdup   %ymm\h,%ymm\h
   vpblendd    $0xAA,%ymm12,%ymm\h,%ymm\h
   vmovshdup   %ymm13,%ymm13
   vpblendd    $0xAA,%ymm14,%ymm13,%ymm13
   vpsubd      %ymm13,%ymm\h,%ymm\h
.endm

.macro shuffle3 s0,s1,s2,r0, r1, r2
   vpunpckldq  %ymm\s1,   %ymm\s0,  %ymm12
   vpunpckhdq  %ymm\s1,   %ymm\s0,  %ymm13
   vpunpckldq  %ymm\s1,   %ymm\s2,  %ymm14
   vpunpckhdq  %ymm\s1,   %ymm\s2,  %ymm15
   vpunpcklqdq %ymm14,   %ymm12,  %ymm\r0
   vpunpckhqdq %ymm14,   %ymm12,  %ymm\r1
   vpunpcklqdq %ymm15,   %ymm13,  %ymm\r2
.endm

.global invntt1_asm
invntt1_asm:
   # initial settng
   vmovdqa _8XQ*4(%rsi),%ymm0
   mov     %rdi,   %r8
   mov     %rsi,   %rcx
   xor     %rax,   %rax
   xor     %rdx,   %rdx

   .p2align 5  # for Cache Hit
   __loop3_out:
   vmovdqu		  0*4(%rdi),%xmm4 
   vmovdqu	 	  3*4(%rdi),%xmm5 
   vmovdqu		  6*4(%rdi),%xmm6
   vmovdqu		  9*4(%rdi),%xmm7
   vmovdqu	 	 12*4(%rdi),%xmm8
   vmovdqu		 15*4(%rdi),%xmm9

   # packing 6 coeff
   # ymm4 : [0 1 2, ?,  9, 10, 11, ?]
   # ymm5 : [3 4 5, ?, 12, 13, 14, ?]
   # ymm6 : [6 7 8, ?, 15, 16, 17, ?] 
   vinserti128  $1, %xmm7, %ymm4, %ymm4
   vinserti128  $1, %xmm8, %ymm5, %ymm5
   vinserti128  $1, %xmm9, %ymm6, %ymm6
   
   # radix-3 : layer 8 ***********************************************************
   # ymm4 : [0 3 6, ?,  9, 12, 15, ?]
   # ymm5 : [1 4 7, ?, 10, 13, 16, ?]
   # ymm6 : [2 5 8, ?, 11, 14, 17, ?] 
   shuffle3 4, 5, 6, 4, 5, 6

   vpsubd      %ymm4, %ymm5, %ymm3     # t1 = Out[j] - Out[j + len];
   vpsubd      %ymm4, %ymm6, %ymm10    # t2 = Out[j + 2 * len] - Out[j];
   vpbroadcastd   _8XWmont_QINV*4(%rcx),%ymm1
   vpbroadcastd   _8XWmont*4(%rcx),%ymm2
   mont_reduce 3                       # t1 = montgomery_reduce((int64_t)Wmont * t1);

   vpaddd      %ymm4, %ymm5, %ymm4     # Out[j] = Out[j] + Out[j + len];
   vpaddd      %ymm4, %ymm6, %ymm4     # Out[j] = Out[j] + Out[j + len] + Out[j + 2 * len];
   vpaddd      %ymm3, %ymm10, %ymm11   # t2 = t2 + t1;
   vpsubd      %ymm3, %ymm6, %ymm3     # t3 = Out[j + 2 * len] - t1;
   vpsubd      %ymm5, %ymm3, %ymm3     # t3 = Out[j + 2 * len] - Out[j + len] - t1;

   vmovdqu     (_INVZETAS_QINV+0)*4(%rsi),%ymm1
   vmovdqu     (_INVZETAS+0)*4(%rsi),%ymm2
   vpsrlq		$32,%ymm1,%ymm10
   vmovshdup	%ymm2,%ymm15
   mont_reduce 11,1,10,2,15 # t1 = montgomery_reduce((int64_t)zeta1 * Out[j + len]);

   vmovdqu     (_INVZETAS_QINV+8)*4(%rsi),%ymm1
   vmovdqu     (_INVZETAS+8)*4(%rsi),%ymm2
   vpsrlq		$32,%ymm1,%ymm10
   vmovshdup	%ymm2,%ymm15
   mont_reduce 3,1,10,2,15 # # t2 = montgomery_reduce((int64_t)zeta2 * Out[j + 2*len]);

   # ymm4 : [0 1 2, ?,  9, 10, 11, ?]
   # ymm5 : [3 4 5, ?, 12, 13, 14, ?]
   # ymm6 : [6 7 8, ?, 15, 16, 17, ?] 
   shuffle3 4, 11, 3, 4, 5, 6

   # radix-3 : layer 7 ***********************************************************
   vpsubd      %ymm4, %ymm5, %ymm3     # t1 = Out[j] - Out[j + len];
   vpsubd      %ymm4, %ymm6, %ymm10    # t2 = Out[j + 2 * len] - Out[j];
   vpbroadcastd   _8XWmont_QINV*4(%rcx),%ymm1
   vpbroadcastd   _8XWmont*4(%rcx),%ymm2
   mont_reduce 3                       # t1 = montgomery_reduce((int64_t)Wmont * t1);

   vpaddd      %ymm4, %ymm5, %ymm4     # Out[j] = Out[j] + Out[j + len];
   vpaddd      %ymm4, %ymm6, %ymm4     # Out[j] = Out[j] + Out[j + len] + Out[j + 2 * len];
   vpaddd      %ymm3, %ymm10, %ymm11   # t2 = t2 + t1;
   vpsubd      %ymm3, %ymm6, %ymm3     # t3 = Out[j + 2 * len] - t1;
   vpsubd      %ymm5, %ymm3, %ymm3     # t3 = Out[j + 2 * len] - Out[j + len] - t1;

   vmovdqu     (_INVZETAS_QINV+16)*4(%rsi),%ymm1
   vmovdqu     (_INVZETAS+16)*4(%rsi),%ymm2
   vpsrlq		$32,%ymm1,%ymm10
   vmovshdup	%ymm2,%ymm15
   mont_reduce 11,1,10,2,15 # t1 = montgomery_reduce((int64_t)zeta1 * Out[j + len]);

   vmovdqu     (_INVZETAS_QINV+24)*4(%rsi),%ymm1
   vmovdqu     (_INVZETAS+24)*4(%rsi),%ymm2
   vpsrlq		$32,%ymm1,%ymm10
   vmovshdup	%ymm2,%ymm15
   mont_reduce 3,1,10,2,15 # # t2 = montgomery_reduce((int64_t)zeta2 * Out[j + 2*len]);

   vextracti128 $1, %ymm4,    %xmm7
   vextracti128 $1, %ymm11,   %xmm8
   vextracti128 $1, %ymm3,    %xmm9

   vmovdqu		6*4(%rdi),%xmm14
   vpblendd    $0xF8,%xmm14,%xmm3,%xmm3
   vmovdqu		15*4(%rdi),%xmm15
   vpblendd    $0xF8,%xmm15,%xmm9,%xmm9
   vmovdqu		 %xmm4,   0*4(%rdi) 
   vmovdqu	 	 %xmm11,  3*4(%rdi) 
   vmovdqu		 %xmm3,   6*4(%rdi)
   vmovdqu		 %xmm7,   9*4(%rdi)
   vmovdqu	 	 %xmm8,  12*4(%rdi)
   vmovdqu		 %xmm9,  15*4(%rdi)

   add     $128,  %rsi
   add     $72,   %rdi
   add     $18,   %rdx
   cmp     $1152, %rdx
   jb  __loop3_out

   mov     %rcx,   %rsi
   mov     %r8,    %rdi   
   xor     %rdx,   %rdx
   .p2align 5  # for Cache Hit
   __loop2_out:
   
   vmovdqu		 0*4(%rdi),%ymm4
   vmovdqu	 	 1*4(%rdi),%ymm5
   vmovdqu		 9*4(%rdi),%ymm6
   vmovdqu	 	10*4(%rdi),%ymm7
   vmovdqu		18*4(%rdi),%ymm8
   vmovdqu     19*4(%rdi),%ymm9
   vmovdqu		27*4(%rdi),%ymm10
   vmovdqu	 	28*4(%rdi),%ymm11

   # radix-2 : layer 6
   vpbroadcastd   (_INVZETAS_QINV+2048+0)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+2048+0)*4(%rsi),%ymm2
   rdx2_butterfly   4,6
   rdx2_butterfly   5,7
   
   vpbroadcastd   (_INVZETAS_QINV+2048+1)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+2048+1)*4(%rsi),%ymm2
   rdx2_butterfly   8,10
   rdx2_butterfly   9,11

   # radix-2 : layer 5
   vpbroadcastd   (_INVZETAS_QINV+2048+2)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+2048+2)*4(%rsi),%ymm2

   rdx2_butterfly	4,8
   rdx2_butterfly	5,9
   rdx2_butterfly	6,10
   rdx2_butterfly	7,11

   vmovdqu     %ymm4 ,  0*4(%rdi)
   vmovdqu     %ymm5 ,  1*4(%rdi)
   vmovdqu     %ymm6 ,  9*4(%rdi)
   vmovdqu     %ymm7 , 10*4(%rdi)
   vmovdqu     %ymm8 , 18*4(%rdi)
   vmovdqu     %ymm9 , 19*4(%rdi)
   vmovdqu     %ymm10, 27*4(%rdi)
   vmovdqu     %ymm11, 28*4(%rdi)
   
   add     $144,   %rdi
   add     $12,    %rsi
   add     $36,    %rdx
   cmp     $1152,  %rdx
   jb  __loop2_out

   mov     %rcx,   %rsi
   mov     %r8,    %rdi   
   xor     %rdx,   %rdx
   xor     %rax,   %rax 

   .p2align 5 # for Cache Hit
   __loop1_out:
   __loop1_in:
   vmovdqu		  0*4(%rdi),%ymm4
   vmovdqu		 36*4(%rdi),%ymm5
   vmovdqu		 72*4(%rdi),%ymm6
   vmovdqu	 	108*4(%rdi),%ymm7
   vmovdqu		144*4(%rdi),%ymm8
   vmovdqu		180*4(%rdi),%ymm9
   vmovdqu		216*4(%rdi),%ymm10
   vmovdqu	 	252*4(%rdi),%ymm11
   
   # radix-2 : layer 4
   vpbroadcastd   (_INVZETAS_QINV+2144+0)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+2144+0)*4(%rsi),%ymm2
   rdx2_butterfly   4,5

   vpbroadcastd   (_INVZETAS_QINV+2144+1)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+2144+1)*4(%rsi),%ymm2
   rdx2_butterfly   6,7

   vpbroadcastd   (_INVZETAS_QINV+2144+2)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+2144+2)*4(%rsi),%ymm2
   rdx2_butterfly   8,9

   vpbroadcastd   (_INVZETAS_QINV+2144+3)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+2144+3)*4(%rsi),%ymm2
   rdx2_butterfly   10,11

   # radix-2 : layer 3
   vpbroadcastd   (_INVZETAS_QINV+2144+4)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+2144+4)*4(%rsi),%ymm2
   rdx2_butterfly   4,6
   rdx2_butterfly   5,7
   
   vpbroadcastd   (_INVZETAS_QINV+2144+5)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+2144+5)*4(%rsi),%ymm2
   rdx2_butterfly   8,10
   rdx2_butterfly   9,11

   # radix-2 : layer 2
   vpbroadcastd   (_INVZETAS_QINV+2144+6)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+2144+6)*4(%rsi),%ymm2
   rdx2_butterfly	4,8
   rdx2_butterfly	5,9
   rdx2_butterfly	6,10
   rdx2_butterfly	7,11
   
   vmovdqu     %ymm4 ,   0*4(%rdi)
   vmovdqu     %ymm5 ,  36*4(%rdi)
   vmovdqu     %ymm6 ,  72*4(%rdi)
   vmovdqu     %ymm7 , 108*4(%rdi)
   vmovdqu     %ymm8 , 144*4(%rdi)
   vmovdqu     %ymm9 , 180*4(%rdi)
   vmovdqu     %ymm10, 216*4(%rdi)
   vmovdqu     %ymm11, 252*4(%rdi)

   add     $32,    %rdi
   add     $32,    %rax
   cmp     $128,   %rax
   jb      __loop1_in

   xor     %rax,   %rax
   sub     $128,   %rdi

   vmovdqu		 32*4(%rdi),%xmm4
   vmovdqu	 	 68*4(%rdi),%xmm5
   vmovdqu		104*4(%rdi),%xmm6
   vmovdqu	 	140*4(%rdi),%xmm7
   vmovdqu		176*4(%rdi),%xmm8
   vmovdqu		212*4(%rdi),%xmm9
   vmovdqu		248*4(%rdi),%xmm10
   vmovdqu	 	284*4(%rdi),%xmm11
   
   # radix-2 : layer 4
   vpbroadcastd   (_INVZETAS_QINV+2144+0)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+2144+0)*4(%rsi),%ymm2
   rdx2_butterfly   4,5

   vpbroadcastd   (_INVZETAS_QINV+2144+1)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+2144+1)*4(%rsi),%ymm2
   rdx2_butterfly   6,7

   vpbroadcastd   (_INVZETAS_QINV+2144+2)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+2144+2)*4(%rsi),%ymm2
   rdx2_butterfly   8,9

   vpbroadcastd   (_INVZETAS_QINV+2144+3)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+2144+3)*4(%rsi),%ymm2
   rdx2_butterfly   10,11

   # radix-2 : layer 3
   vpbroadcastd   (_INVZETAS_QINV+2144+4)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+2144+4)*4(%rsi),%ymm2
   rdx2_butterfly   4,6
   rdx2_butterfly   5,7
   
   vpbroadcastd   (_INVZETAS_QINV+2144+5)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+2144+5)*4(%rsi),%ymm2
   rdx2_butterfly   8,10
   rdx2_butterfly   9,11

   # radix-2 : layer 2
   vpbroadcastd   (_INVZETAS_QINV+2144+6)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+2144+6)*4(%rsi),%ymm2
   rdx2_butterfly	4,8
   rdx2_butterfly	5,9
   rdx2_butterfly	6,10
   rdx2_butterfly	7,11

   vmovdqu     %xmm4 ,  32*4(%rdi)
   vmovdqu     %xmm5 ,  68*4(%rdi)
   vmovdqu     %xmm6 , 104*4(%rdi)
   vmovdqu     %xmm7 , 140*4(%rdi)
   vmovdqu     %xmm8 , 176*4(%rdi)
   vmovdqu     %xmm9 , 212*4(%rdi)
   vmovdqu     %xmm10, 248*4(%rdi)
   vmovdqu     %xmm11, 284*4(%rdi)

   add     $1152,  %rdi
   add     $28,    %rsi
   add     $1,     %rdx
   cmp     $4,     %rdx
   jb  __loop1_out

   mov     %rcx,   %rsi
   mov     %r8,    %rdi      
   xor     %rax,   %rax 
   .p2align 5 # for Cache Hit
   vpbroadcastd   (_8XF2_QINV)*4(%rsi),%ymm3
   vpbroadcastd   (_8XF2)*4(%rsi),%ymm15
   __loop0:
   vmovdqa         4*0(%rdi),%ymm4
   vmovdqa        4*144(%rdi),%ymm5
   vmovdqa        4*288(%rdi),%ymm6
   vmovdqa        4*432(%rdi),%ymm7
   vmovdqa        4*576(%rdi),%ymm8
   vmovdqa        4*720(%rdi),%ymm9
   vmovdqa        4*864(%rdi),%ymm10
   vmovdqa       4*1008(%rdi),%ymm11

   # radix-2 : layer 1
   vpbroadcastd   (_INVZETAS_QINV+2172+0)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+2172+0)*4(%rsi),%ymm2
   rdx2_butterfly   4,6
   rdx2_butterfly   5,7
   
   vpbroadcastd   (_INVZETAS_QINV+2172+1)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+2172+1)*4(%rsi),%ymm2
   rdx2_butterfly   8,10
   rdx2_butterfly   9,11

   # trinomial : layer 0
   vpbroadcastd    (_INVZETAS_QINV+2172+2)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+2172+2)*4(%rsi),%ymm2

   rdx2_butterfly   4,8 
   rdx2_butterfly   5,9 
   rdx2_butterfly   6,10
   rdx2_butterfly   7,11

   vpbroadcastd   (_8XF1_QINV)*4(%rsi),%ymm1
   vpbroadcastd   (_8XF1)*4(%rsi),%ymm2
   vpsubd      %ymm8,%ymm4,%ymm4
   vpsubd      %ymm9,%ymm5,%ymm5
   vpsubd      %ymm10,%ymm6,%ymm6
   vpsubd      %ymm11,%ymm7,%ymm7

   mont_reduce 4
   mont_reduce 8, 3, 3, 15, 15
   mont_reduce 5
   mont_reduce 9, 3, 3, 15, 15
   mont_reduce 6
   mont_reduce 10, 3, 3, 15, 15
   mont_reduce 7
   mont_reduce 11, 3, 3, 15, 15

   vmovdqa     %ymm4 ,    4*0(%rdi)
   vmovdqa     %ymm5 ,  4*144(%rdi)
   vmovdqa     %ymm6 ,  4*288(%rdi)
   vmovdqa     %ymm7 ,  4*432(%rdi)
   vmovdqa     %ymm8 ,  4*576(%rdi)
   vmovdqa     %ymm9 ,  4*720(%rdi)
   vmovdqa     %ymm10,  4*864(%rdi)
   vmovdqa     %ymm11, 4*1008(%rdi)

   add     $32,    %rdi
   add     $32,    %rax
   cmp     $576,   %rax
   jb      __loop0

   ret