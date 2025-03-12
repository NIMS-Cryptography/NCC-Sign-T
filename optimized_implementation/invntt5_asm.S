#include "consts.h"

.macro shuffle8 r0,r1,r2,r3
vperm2i128	$0x20,%ymm\r1,%ymm\r0,%ymm\r2
vperm2i128	$0x31,%ymm\r1,%ymm\r0,%ymm\r3
.endm

.macro shuffle4 r0,r1,r2,r3
vpunpcklqdq	%ymm\r1,%ymm\r0,%ymm\r2
vpunpckhqdq	%ymm\r1,%ymm\r0,%ymm\r3
.endm

.macro shuffle2 r0,r1,r2,r3
vmovsldup	%ymm\r1,%ymm\r2
vpblendd	$0xAA,%ymm\r2,%ymm\r0,%ymm\r2
vpsrlq		$32,%ymm\r0,%ymm\r0
vpblendd	$0xAA,%ymm\r1,%ymm\r0,%ymm\r3
.endm


.macro forwardshuffle r0,r1
    vpunpckldq 	%ymm\r1,%ymm\r0,%ymm1
    vpunpckhdq 	%ymm\r1,%ymm\r0,%ymm2
    vpunpcklqdq %ymm2, %ymm1,%ymm12
    vpunpckhqdq %ymm2, %ymm1,%ymm13
    vpermd %ymm12, %ymm3, %ymm\r0
    vpermd %ymm13, %ymm3, %ymm\r1
.endm




.macro butterfly l,h,zl0=1,zl1=1,zh0=2,zh1=2
vpsubd		%ymm\h,%ymm\l,%ymm12 // a-b = t
vpaddd		%ymm\h,%ymm\l,%ymm\l // a + b = a

vpmuldq		%ymm\zl0,%ymm12,%ymm13 // t * zqinv = r0
vmovshdup	%ymm12,%ymm\h // t >> 32 = b
vpmuldq		%ymm\zl1,%ymm\h,%ymm14 //b*zqinv = r1

vpmuldq		%ymm\zh0,%ymm12,%ymm12 // t*zeta = t
vpmuldq		%ymm\zh1,%ymm\h,%ymm\h // b*zeta = b

vpmuldq		%ymm0,%ymm13,%ymm13 // r0 * q = r0
vpmuldq		%ymm0,%ymm14,%ymm14 // r1 * q = r1

vpsubd		%ymm13,%ymm12,%ymm12 // t-r0 = t
vpsubd		%ymm14,%ymm\h,%ymm\h // b-r1 = b

vmovshdup	%ymm12,%ymm12
vpblendd	$0xAA,%ymm\h,%ymm12,%ymm\h
.endm


.macro butterfly0 l,h,zl0=1,zl1=1,zh0=2,zh1=2
vpsubd		%ymm\h,%ymm\l,%ymm12 // a-b = t
vpaddd		%ymm\h,%ymm\l,%ymm\l // a + b = a

vpmuldq		%ymm\zl0,%ymm12,%ymm13 // t * zqinv = r0
vmovshdup	%ymm12,%ymm\h // t >> 32 = b
vmovshdup	%ymm\zl1,%ymm\zl1  // zqinv >> 32
vpmuldq		%ymm\zl1,%ymm\h,%ymm14 //b*zqinv = r1

vpmuldq		%ymm\zh0,%ymm12,%ymm12 // t*zeta = t
vmovshdup	%ymm\zh1,%ymm\zh1  // zqinv >> 32
vpmuldq		%ymm\zh1,%ymm\h,%ymm\h // b*zeta = b

vpmuldq		%ymm0,%ymm13,%ymm13 // r0 * q = r0
vpmuldq		%ymm0,%ymm14,%ymm14 // r1 * q = r1

vpsubd		%ymm13,%ymm12,%ymm12 // t-r0 = t
vpsubd		%ymm14,%ymm\h,%ymm\h // b-r1 = b

vmovshdup	%ymm12,%ymm12
vpblendd	$0xAA,%ymm\h,%ymm12,%ymm\h
.endm


.macro reduce h
vpmuldq      %ymm1,%ymm\h,%ymm13  // dqiv * b = r0
vmovshdup    %ymm\h,%ymm12 // b >> 32 = t
vpmuldq      %ymm1,%ymm12,%ymm14 // t*dqiv = r1

vpmuldq      %ymm2,%ymm12,%ymm12 //t*div = t
vpmuldq      %ymm2,%ymm\h,%ymm\h // b*div = b

vpmuldq      %ymm0,%ymm13,%ymm13 //r0 *q = r0
vpmuldq      %ymm0,%ymm14,%ymm14 //r1 *q = r1

vpsubd      %ymm14,%ymm12,%ymm12 //t - r1 = t
vpsubd      %ymm13,%ymm\h,%ymm\h //b - r0 = b
vmovshdup   %ymm\h,%ymm\h

vpblendd    $0xAA,%ymm12,%ymm\h,%ymm\h 
.endm



.macro levels0t2 off

vmovdqu	(_INVZETAS_QINV+\off)*4(%rsi),%ymm1
vmovdqu	(_INVZETAS+\off)*4(%rsi),%ymm2
butterfly0 4,5
shuffle2 4,5,3,5

vmovdqu	(_INVZETAS_QINV+\off+8)*4(%rsi),%ymm1
vmovdqu	(_INVZETAS+\off+8)*4(%rsi),%ymm2
butterfly 3,5
shuffle4 3,5,4,5

vmovdqu	(_INVZETAS_QINV+\off+16)*4(%rsi),%ymm1
vmovdqu	(_INVZETAS+\off+16)*4(%rsi),%ymm2
butterfly 4,5
shuffle8 4,5,3,5



vmovdqu	(_INVZETAS_QINV+\off+24)*4(%rsi),%ymm1
vmovdqu	(_INVZETAS+\off+24)*4(%rsi),%ymm2
butterfly0 6,7
shuffle2 6,7,4,7
vmovdqu	(_INVZETAS_QINV+\off+32)*4(%rsi),%ymm1
vmovdqu	(_INVZETAS+\off+32)*4(%rsi),%ymm2
butterfly 4,7
shuffle4 4,7,6,7
vmovdqu	(_INVZETAS_QINV+\off+40)*4(%rsi),%ymm1
vmovdqu	(_INVZETAS+\off+40)*4(%rsi),%ymm2
butterfly 6,7
shuffle8 6,7,4,7



vmovdqu	(_INVZETAS_QINV+\off+48)*4(%rsi),%ymm1
vmovdqu	(_INVZETAS+\off+48)*4(%rsi),%ymm2
butterfly0 8,9
shuffle2 8,9,6,9
vmovdqu	(_INVZETAS_QINV+\off+56)*4(%rsi),%ymm1
vmovdqu	(_INVZETAS+\off+56)*4(%rsi),%ymm2
butterfly 6,9
shuffle4 6,9,8,9
vmovdqu	(_INVZETAS_QINV+\off+64)*4(%rsi),%ymm1
vmovdqu	(_INVZETAS+\off+64)*4(%rsi),%ymm2
butterfly 8,9
shuffle8 8,9,6,9



vmovdqu	(_INVZETAS_QINV+\off+72)*4(%rsi),%ymm1
vmovdqu	(_INVZETAS+\off+72)*4(%rsi),%ymm2
butterfly0 10,11
shuffle2 10,11,8,11

vmovdqu	(_INVZETAS_QINV+\off+80)*4(%rsi),%ymm1
vmovdqu	(_INVZETAS+\off+80)*4(%rsi),%ymm2
butterfly 8,11
shuffle4 8,11,10,11

vmovdqu	(_INVZETAS_QINV+\off+88)*4(%rsi),%ymm1
vmovdqu	(_INVZETAS+\off+88)*4(%rsi),%ymm2
butterfly 10,11
shuffle8 10,11,8,11

.endm



.macro levels0t5 off, off2
vmovdqa		  0+ 256*\off2(%rdi),%ymm4   
vmovdqa		32+256*\off2(%rdi),%ymm5   
vmovdqa		64+256*\off2(%rdi),%ymm6   
vmovdqa	 	96+256*\off2(%rdi),%ymm7   
vmovdqa	    128+256*\off2(%rdi),%ymm8   
vmovdqa		160+256*\off2(%rdi),%ymm9   
vmovdqa		192+256*\off2(%rdi),%ymm10  
vmovdqa	 	224+256*\off2(%rdi),%ymm11  

/* level 0 */

vmovdqu	(_8idx)*4(%rsi),%ymm3
forwardshuffle 4,5
forwardshuffle 6,7
forwardshuffle 8,9
forwardshuffle 10,11

levels0t2 \off
vpbroadcastd	(_INVZETAS_QINV+\off + 96)*4(%rsi),%ymm1
vpbroadcastd	(_INVZETAS+\off + 96)*4(%rsi),%ymm2
butterfly 3,5

vpbroadcastd	(_INVZETAS_QINV+\off + 97)*4(%rsi),%ymm1
vpbroadcastd	(_INVZETAS+\off + 97)*4(%rsi),%ymm2
butterfly 4,7

vpbroadcastd	(_INVZETAS_QINV+\off + 98)*4(%rsi),%ymm1
vpbroadcastd	(_INVZETAS+\off + 98)*4(%rsi),%ymm2
butterfly 6,9

vpbroadcastd	(_INVZETAS_QINV+\off + 99)*4(%rsi),%ymm1
vpbroadcastd	(_INVZETAS+\off + 99)*4(%rsi),%ymm2
butterfly 8,11

vpbroadcastd	(_INVZETAS_QINV+\off + 100)*4(%rsi),%ymm1
vpbroadcastd	(_INVZETAS+\off + 100)*4(%rsi),%ymm2
butterfly 3,4
butterfly 5,7

vpbroadcastd	(_INVZETAS_QINV+\off + 101)*4(%rsi),%ymm1
vpbroadcastd	(_INVZETAS+\off + 101)*4(%rsi),%ymm2
butterfly 6,8
butterfly 9,11

vpbroadcastd	(_INVZETAS_QINV+\off + 102)*4(%rsi),%ymm1
vpbroadcastd	(_INVZETAS+\off + 102)*4(%rsi),%ymm2
butterfly 3,6
butterfly 5,9
butterfly 4,8
butterfly 7,11


vmovdqa		%ymm3,256*\off2+  0(%rdi)
vmovdqa		%ymm5,256*\off2+ 32(%rdi)
vmovdqa		%ymm4,256*\off2+ 64(%rdi)
vmovdqa		%ymm7,256*\off2+ 96(%rdi)
vmovdqa		%ymm6,256*\off2+128(%rdi)
vmovdqa		%ymm9,256*\off2+160(%rdi)
vmovdqa		%ymm8,256*\off2+192(%rdi)
vmovdqa		%ymm11,256*\off2+224(%rdi)
.endm




.macro levels2t4 off , off2
/* level 0 */


vmovdqa		  0+ 32*\off+2048*\off2(%rdi),%ymm4   
vmovdqa		256+32*\off+2048*\off2(%rdi),%ymm5   
vmovdqa		512+32*\off+2048*\off2(%rdi),%ymm6   
vmovdqa	 	768+32*\off+2048*\off2(%rdi),%ymm7   
vmovdqa	    1024+32*\off+2048*\off2(%rdi),%ymm8   
vmovdqa		1280+32*\off+2048*\off2(%rdi),%ymm9   
vmovdqa		1536+32*\off+2048*\off2(%rdi),%ymm10  
vmovdqa	 	1792+32*\off+2048*\off2(%rdi),%ymm11  

vpbroadcastd	(_INVZETAS_QINV+3296+ 7*\off2)*4(%rsi),%ymm1
vpbroadcastd	(_INVZETAS+ 3296 +7*\off2)*4(%rsi),%ymm2

butterfly	4,5
/* level 1 */
vpbroadcastd	(_INVZETAS_QINV+3297+7*\off2)*4(%rsi),%ymm1
vpbroadcastd	(_INVZETAS+3297+7*\off2)*4(%rsi),%ymm2
butterfly	6,7

vpbroadcastd	(_INVZETAS_QINV+3298+7*\off2)*4(%rsi),%ymm1
vpbroadcastd	(_INVZETAS+3298+7*\off2)*4(%rsi),%ymm2
butterfly	8,9

vpbroadcastd	(_INVZETAS_QINV+3299+7*\off2)*4(%rsi),%ymm1
vpbroadcastd	(_INVZETAS+3299+7*\off2)*4(%rsi),%ymm2
butterfly	10,11

vpbroadcastd	(_INVZETAS_QINV+3300+7*\off2)*4(%rsi),%ymm1
vpbroadcastd	(_INVZETAS+3300+7*\off2)*4(%rsi),%ymm2
butterfly	4,6
butterfly	5,7

vpbroadcastd	(_INVZETAS_QINV+3301+7*\off2)*4(%rsi),%ymm1
vpbroadcastd	(_INVZETAS+3301+7*\off2)*4(%rsi),%ymm2
butterfly	8,10
butterfly	9,11


vpbroadcastd	(_INVZETAS_QINV+3302+7*\off2)*4(%rsi),%ymm1
vpbroadcastd	(_INVZETAS+3302+7*\off2)*4(%rsi),%ymm2
butterfly	4,8
butterfly	5,9
butterfly	6,10
butterfly	7,11

vmovdqa		%ymm4,  0+32*\off+2048*\off2(%rdi)
vmovdqa		%ymm5,256+32*\off+2048*\off2(%rdi)
vmovdqa		%ymm6,512+32*\off+2048*\off2(%rdi)
vmovdqa		%ymm7,768+32*\off+2048*\off2(%rdi)
vmovdqa		%ymm8,1024+32*\off+2048*\off2(%rdi)
vmovdqa		%ymm9,1280+32*\off+2048*\off2(%rdi)
vmovdqa		%ymm10,1536+32*\off+2048*\off2(%rdi)
vmovdqa		%ymm11,1792+32*\off+2048*\off2(%rdi)
.endm



.macro levels0t1 off

vmovdqa		  0+64*\off(%rdi),%ymm4   
vmovdqa		32+64*\off(%rdi),%ymm5   
vmovdqa		2048+64*\off(%rdi),%ymm6   
vmovdqa	 	2080+64*\off(%rdi),%ymm7   
vmovdqa		4096+64*\off(%rdi),%ymm8   
vmovdqa		4128+64*\off(%rdi),%ymm9   
vmovdqa		6144+64*\off(%rdi),%ymm10  
vmovdqa	 	6176+64*\off(%rdi),%ymm11  

vpbroadcastd	(_INVZETAS_QINV+3324)*4(%rsi),%ymm1
vpbroadcastd	(_INVZETAS+3324)*4(%rsi),%ymm2
butterfly 4,6
butterfly 5,7


vpbroadcastd	(_INVZETAS_QINV+3325)*4(%rsi),%ymm1
vpbroadcastd	(_INVZETAS+3325)*4(%rsi),%ymm2
butterfly	8,10
butterfly	9,11

vpbroadcastd	(_INVZETAS_QINV+3326)*4(%rsi),%ymm1
vpbroadcastd	(_INVZETAS+3326)*4(%rsi),%ymm2
butterfly	4,8
butterfly	5,9
butterfly	6,10
butterfly	7,11


vpbroadcastd	(_8XF1_QINV)*4(%rsi),%ymm1
vpbroadcastd	(_8XF1)*4(%rsi),%ymm2
vpsubd %ymm8,%ymm4,%ymm4
vpsubd %ymm9,%ymm5,%ymm5
vpsubd %ymm10,%ymm6,%ymm6
vpsubd %ymm11,%ymm7,%ymm7
reduce 4
reduce 5
reduce 6
reduce 7
vpbroadcastd	(_8XF2_QINV)*4(%rsi),%ymm1
vpbroadcastd	(_8XF2)*4(%rsi),%ymm2
reduce 8
reduce 9
reduce 10
reduce 11

vmovdqa		%ymm4,  0+64*\off(%rdi)
vmovdqa		%ymm5,32+64*\off(%rdi)
vmovdqa		%ymm6,2048+64*\off(%rdi)
vmovdqa		%ymm7,2080+64*\off(%rdi)
vmovdqa		%ymm8,4096+64*\off(%rdi)
vmovdqa		%ymm9,4128+64*\off(%rdi)
vmovdqa		%ymm10,6144+64*\off(%rdi)
vmovdqa		%ymm11,6176+64*\off(%rdi)
.endm



.text
.global invntt5_asm
invntt5_asm:
vmovdqa		_8XQ*4(%rsi),%ymm0

levels0t5	0, 0
levels0t5	103, 1
levels0t5	206, 2
levels0t5	309, 3
levels0t5	412, 4
levels0t5	515, 5
levels0t5	618, 6
levels0t5	721, 7

levels0t5	824, 8
levels0t5	927, 9
levels0t5	1030, 10
levels0t5	1133, 11
levels0t5	1236, 12
levels0t5	1339, 13
levels0t5	1442, 14
levels0t5	1545, 15

levels0t5	1648, 16
levels0t5	1751, 17
levels0t5	1854, 18
levels0t5	1957, 19
levels0t5	2060, 20
levels0t5	2163, 21
levels0t5	2266, 22
levels0t5	2369, 23


levels0t5	2472, 24
levels0t5	2575, 25
levels0t5	2678, 26
levels0t5	2781, 27
levels0t5	2884, 28
levels0t5	2987, 29
levels0t5	3090, 30
levels0t5	3193, 31



levels2t4	0,0
levels2t4	1,0
levels2t4	2,0
levels2t4	3,0
levels2t4	4,0
levels2t4	5,0
levels2t4	6,0
levels2t4	7,0

levels2t4	0,1
levels2t4	1,1
levels2t4	2,1
levels2t4	3,1
levels2t4	4,1
levels2t4	5,1
levels2t4	6,1
levels2t4	7,1

levels2t4	0,2
levels2t4	1,2
levels2t4	2,2
levels2t4	3,2
levels2t4	4,2
levels2t4	5,2
levels2t4	6,2
levels2t4	7,2


levels2t4	0,3
levels2t4	1,3
levels2t4	2,3
levels2t4	3,3
levels2t4	4,3
levels2t4	5,3
levels2t4	6,3
levels2t4	7,3


levels0t1	0
levels0t1	1
levels0t1	2
levels0t1	3
levels0t1	4
levels0t1	5
levels0t1	6
levels0t1	7
levels0t1	8
levels0t1	9
levels0t1	10
levels0t1	11
levels0t1	12
levels0t1	13
levels0t1	14
levels0t1	15
levels0t1	16
levels0t1	17
levels0t1	18
levels0t1	19
levels0t1	20
levels0t1	21
levels0t1	22
levels0t1	23
levels0t1	24
levels0t1	25
levels0t1	26
levels0t1	27
levels0t1	28
levels0t1	29
levels0t1	30
levels0t1	31
ret

