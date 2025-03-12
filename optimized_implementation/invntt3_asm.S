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

.global invntt3_asm
invntt3_asm:
   # initial settng
   vmovdqa _8XQ*4(%rsi),%ymm0
   mov     %rsi,   %rcx
   mov     %rdi,   %r8
   xor     %rax,   %rax

   .p2align 5 # for Cache Hit
   # radix-2 : layer 8
   __loop3_out:
   vmovdqu         4*0(%rdi),%ymm4 #1
   vmovdqu         4*3(%rdi),%ymm5 #1
   vmovdqu         4*3(%rdi),%ymm3 #1
   vpbroadcastd   (_INVZETAS_QINV+0)*4(%rsi),%ymm1  #1
   vpbroadcastd   (_INVZETAS+0)*4(%rsi),%ymm2  #1
   rdx2_butterfly   4,5 #1
   vmovdqu         4*6(%rdi),%ymm6 #2
   vmovdqu         4*9(%rdi),%ymm7 #2
   vmovdqu         4*9(%rdi),%ymm8 #2
   vpblendd $0xF8, %ymm3, %ymm5, %ymm5 #1
   vmovdqu     %ymm4 ,  4*0(%rdi) #1
   vmovdqu     %ymm5 ,  4*3(%rdi) #1

   vpbroadcastd   (_INVZETAS_QINV+1)*4(%rsi),%ymm1 #2
   vpbroadcastd   (_INVZETAS+1)*4(%rsi),%ymm2 #2
   rdx2_butterfly   6,7 #2
   vmovdqu         4*12(%rdi),%ymm9 #3
   vmovdqu         4*15(%rdi),%ymm10 #3
   vmovdqu         4*15(%rdi),%ymm11 #3
   vpblendd $0xF8, %ymm8, %ymm7, %ymm7 #2
   vmovdqu     %ymm6 ,  4*6(%rdi) #2
   vmovdqu     %ymm7 ,  4*9(%rdi) #2

   vpbroadcastd   (_INVZETAS_QINV+2)*4(%rsi),%ymm1 #3
   vpbroadcastd   (_INVZETAS+2)*4(%rsi),%ymm2 #3
   rdx2_butterfly   9,10 #3
   vmovdqu          4*18(%rdi),%ymm4 #4
   vmovdqu         4*21(%rdi),%ymm5 #4
   vmovdqu         4*21(%rdi),%ymm3 #4
   vpblendd $0xF8, %ymm11, %ymm10, %ymm10 #3
   vmovdqu     %ymm9 ,  4*12(%rdi) #3
   vmovdqu     %ymm10 ,  4*15(%rdi) #3

   vpbroadcastd   (_INVZETAS_QINV+3)*4(%rsi),%ymm1  #4
   vpbroadcastd   (_INVZETAS+3)*4(%rsi),%ymm2  #4
   rdx2_butterfly   4,5 #4
   vmovdqu         4*24(%rdi),%ymm6 #5
   vmovdqu         4*27(%rdi),%ymm7 #5
   vmovdqu         4*27(%rdi),%ymm8 #5
   vpblendd $0xF8, %ymm3, %ymm5, %ymm5 #4
   vmovdqu     %ymm4 ,  4*18(%rdi) #4
   vmovdqu     %ymm5 ,  4*21(%rdi) #4
   
   vpbroadcastd   (_INVZETAS_QINV+4)*4(%rsi),%ymm1 #5
   vpbroadcastd   (_INVZETAS+4)*4(%rsi),%ymm2 #5
   rdx2_butterfly   6,7 #5
   vmovdqu         4*30(%rdi),%ymm9 #6
   vmovdqu         4*33(%rdi),%ymm10 #6
   vmovdqu         4*33(%rdi),%ymm11 #6
   vpblendd $0xF8, %ymm8, %ymm7, %ymm7 #5
   vmovdqu     %ymm6 ,  4*24(%rdi) #5
   vmovdqu     %ymm7 ,  4*27(%rdi) #5

   vpbroadcastd   (_INVZETAS_QINV+5)*4(%rsi),%ymm1 #6
   vpbroadcastd   (_INVZETAS+5)*4(%rsi),%ymm2 #6
   rdx2_butterfly   9,10 #6
   
   add     $24,    %rsi
   vpblendd $0xF8, %ymm11, %ymm10, %ymm10 #6
   add     $36,    %rax
   vmovdqu     %ymm9 ,  4*30(%rdi) #6
   vmovdqu     %ymm10 ,  4*33(%rdi) #6
   add     $144,   %rdi
   cmp     $1536, %rax
   jb      __loop3_out

   mov      %r8,   %rdi
   mov     %rcx,   %rsi
   xor     %rax,   %rax
  
   .p2align 5 # for Cache Hit
   __loop2_out:
   vmovdqu          4*0(%rdi),%ymm4 #1
   vmovdqu          4*6(%rdi),%ymm5 #1
   vmovdqu         4*12(%rdi),%ymm6 #1
   vmovdqu         4*18(%rdi),%ymm7 #1
   vmovdqu         4*18(%rdi),%ymm3 #1
   vmovdqu         4*24(%rdi),%ymm8 #2
   vmovdqu         4*30(%rdi),%ymm9 #2
   vmovdqu         4*36(%rdi),%ymm10 #2
   vmovdqu         4*42(%rdi),%ymm11 #2
   vmovdqu         4*42(%rdi),%ymm15 #2
   
   # radix-2 : layer 7
   vpbroadcastd   (_INVZETAS_QINV+256+0)*4(%rsi),%ymm1 #1
   vpbroadcastd   (_INVZETAS+256+0)*4(%rsi),%ymm2 #1
   rdx2_butterfly   4,5 #1

   vpbroadcastd   (_INVZETAS_QINV+256+1)*4(%rsi),%ymm1 #1
   vpbroadcastd   (_INVZETAS+256+1)*4(%rsi),%ymm2 #1
   rdx2_butterfly   6,7 #1

   vpbroadcastd   (_INVZETAS_QINV+256+2)*4(%rsi),%ymm1 #2
   vpbroadcastd   (_INVZETAS+256+2)*4(%rsi),%ymm2 #2
   rdx2_butterfly   8,9 #1

   vpbroadcastd   (_INVZETAS_QINV+256+3)*4(%rsi),%ymm1 #2
   vpbroadcastd   (_INVZETAS+256+3)*4(%rsi),%ymm2 #2
   rdx2_butterfly   10,11 #2

   # radix-2 : layer 6
   vpbroadcastd   (_INVZETAS_QINV+256+4)*4(%rsi),%ymm1 #1
   vpbroadcastd   (_INVZETAS+256+4)*4(%rsi),%ymm2 #1
   rdx2_butterfly   4,6 #1
   rdx2_butterfly   5,7 #1

   vpbroadcastd   (_INVZETAS_QINV+256+5)*4(%rsi),%ymm1 #2
   vpbroadcastd   (_INVZETAS+256+5)*4(%rsi),%ymm2 #2
   rdx2_butterfly   8,10 #2
   rdx2_butterfly   9,11 #2

   vpblendd $0xC0, %ymm3, %ymm7, %ymm7 #1
   vpblendd $0xC0, %ymm15, %ymm11, %ymm11 #2
   vmovdqu     %ymm4 ,    4*0(%rdi) #1
   vmovdqu     %ymm5 ,    4*6(%rdi) #1
   vmovdqu     %ymm6 ,   4*12(%rdi) #1
   vmovdqu     %ymm7 ,   4*18(%rdi) #1
   vmovdqu     %ymm8 ,   4*24(%rdi) #2
   vmovdqu     %ymm9 ,   4*30(%rdi) #2
   vmovdqu     %ymm10 ,  4*36(%rdi) #2
   vmovdqu     %ymm11 ,  4*42(%rdi) #2

   add     $24,   %rsi
   add     $192,  %rdi
   add     $48,   %rax
   cmp     $1536, %rax
   jb      __loop2_out

   mov      %r8,   %rdi
   mov     %rcx,   %rsi
   xor     %rdx,   %rdx 
   xor     %rax,   %rax

   .p2align 5 # for Cache Hit
   __loop1_out:
   xor     %rax,  %rax
   vpbroadcastd    (_INVZETAS_QINV+448+6)*4(%rsi),%ymm3
   vpbroadcastd   (_INVZETAS+448+6)*4(%rsi),%ymm15
   __loop1_in:
   vmovdqu         4*0(%rdi),%ymm4
   vmovdqu        4*24(%rdi),%ymm5
   vmovdqu        4*48(%rdi),%ymm6
   vmovdqu        4*72(%rdi),%ymm7
   vmovdqu        4*96(%rdi),%ymm8
   vmovdqu       4*120(%rdi),%ymm9
   vmovdqu       4*144(%rdi),%ymm10
   vmovdqu       4*168(%rdi),%ymm11

   # radix-2 : layer 5
   vpbroadcastd   (_INVZETAS_QINV+448+0)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+448+0)*4(%rsi),%ymm2
   rdx2_butterfly   4,5

   vpbroadcastd   (_INVZETAS_QINV+448+1)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+448+1)*4(%rsi),%ymm2
   rdx2_butterfly   6,7

   vpbroadcastd   (_INVZETAS_QINV+448+2)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+448+2)*4(%rsi),%ymm2
   rdx2_butterfly   8,9

   vpbroadcastd   (_INVZETAS_QINV+448+3)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+448+3)*4(%rsi),%ymm2
   rdx2_butterfly   10,11

   # radix-2 : layer 4
   vpbroadcastd   (_INVZETAS_QINV+448+4)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+448+4)*4(%rsi),%ymm2
   rdx2_butterfly   4,6
   rdx2_butterfly   5,7

   vpbroadcastd   (_INVZETAS_QINV+448+5)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+448+5)*4(%rsi),%ymm2
   rdx2_butterfly   8,10
   rdx2_butterfly   9,11

   # radix-2 : layer 3
   rdx2_butterfly   4,8 ,3, 3, 15, 15
   rdx2_butterfly   5,9 ,3, 3, 15, 15
   rdx2_butterfly   6,10,3, 3, 15, 15
   rdx2_butterfly   7,11,3, 3, 15, 15

   vmovdqu     %ymm4 ,   4*0(%rdi)
   vmovdqu     %ymm5 ,  4*24(%rdi)
   vmovdqu     %ymm6 ,  4*48(%rdi)
   vmovdqu     %ymm7 ,  4*72(%rdi)
   vmovdqu     %ymm8 ,  4*96(%rdi)
   vmovdqu     %ymm9 , 4*120(%rdi)
   vmovdqu     %ymm10, 4*144(%rdi)
   vmovdqu     %ymm11, 4*168(%rdi)

   add     $32,   %rdi
   add     $8,    %rax
   cmp     $24,   %rax
   jb      __loop1_in

   add     $28,    %rsi
   add     $672,  %rdi
   add     $1,    %rdx
   cmp     $8,    %rdx
   jb      __loop1_out

   mov      %r8,   %rdi
   mov     %rcx,   %rsi
   xor     %rax,   %rax
   .p2align 5 # for Cache Hit
   vpbroadcastd   (_8XF2_QINV)*4(%rsi),%ymm3
   vpbroadcastd   (_8XF2)*4(%rsi),%ymm15
   __loop0:
   vmovdqa          4*0(%rdi),%ymm4
   vmovdqa        4*192(%rdi),%ymm5
   vmovdqa        4*384(%rdi),%ymm6
   vmovdqa        4*576(%rdi),%ymm7
   vmovdqa        4*768(%rdi),%ymm8
   vmovdqa        4*960(%rdi),%ymm9
   vmovdqa       4*1152(%rdi),%ymm10
   vmovdqa       4*1344(%rdi),%ymm11

  # radix-2 : layer 2 : 7424965 | 714BC5
   vpbroadcastd   (_INVZETAS_QINV+504+0)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+504+0)*4(%rsi),%ymm2
   rdx2_butterfly   4,5

   vpbroadcastd   (_INVZETAS_QINV+504+1)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+504+1)*4(%rsi),%ymm2
   rdx2_butterfly   6,7

   vpbroadcastd   (_INVZETAS_QINV+504+2)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+504+2)*4(%rsi),%ymm2
   rdx2_butterfly   8,9

   vpbroadcastd   (_INVZETAS_QINV+504+3)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+504+3)*4(%rsi),%ymm2
   rdx2_butterfly   10,11

   # radix-2 : layer 1
   vpbroadcastd   (_INVZETAS_QINV+504+4)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+504+4)*4(%rsi),%ymm2
   rdx2_butterfly   4,6
   rdx2_butterfly   5,7
 
   vpbroadcastd   (_INVZETAS_QINV+504+5)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+504+5)*4(%rsi),%ymm2
   rdx2_butterfly   8,10
   rdx2_butterfly   9,11

   # trinomial : layer 0
   vpbroadcastd   (_INVZETAS_QINV+504+6)*4(%rsi),%ymm1
   vpbroadcastd   (_INVZETAS+504+6)*4(%rsi),%ymm2

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
   vmovdqa     %ymm5 ,  4*192(%rdi)
   vmovdqa     %ymm6 ,  4*384(%rdi)
   vmovdqa     %ymm7 ,  4*576(%rdi)
   vmovdqa     %ymm8 ,  4*768(%rdi)
   vmovdqa     %ymm9 ,  4*960(%rdi)
   vmovdqa     %ymm10, 4*1152(%rdi)
   vmovdqa     %ymm11, 4*1344(%rdi)

   add     $32,    %rdi
   add     $8,     %rax
   cmp     $192,   %rax
   jb      __loop0

   ret