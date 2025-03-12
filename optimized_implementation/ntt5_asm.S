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



.macro shuffle_inv r0,r1
    vpunpckldq 	%ymm\r1,%ymm\r0,%ymm1
    vpunpckhdq 	%ymm\r1,%ymm\r0,%ymm2

    vperm2i128 	$0x20,%ymm2,%ymm1,%ymm\r0
    vperm2i128 	$0x31,%ymm2,%ymm1,%ymm\r1
.endm

.macro tri_butterfly l,h,zl0=1,zl1=1,zh0=2,zh1=2
vmovaps     %ymm\h,%ymm3
vpmuldq      %ymm\zl0,%ymm\h,%ymm13
vmovshdup   %ymm\h,%ymm12
vpmuldq      %ymm\zl1,%ymm12,%ymm14

vpmuldq      %ymm\zh0,%ymm\h,%ymm\h
vpmuldq      %ymm\zh1,%ymm12,%ymm12

vpmuldq      %ymm0,%ymm13,%ymm13
vpmuldq      %ymm0,%ymm14,%ymm14

vmovshdup   %ymm\h,%ymm\h
vpblendd   $0xAA,%ymm12,%ymm\h,%ymm\h

vpsubd      %ymm\h,%ymm\l,%ymm12
vpaddd      %ymm\h,%ymm\l,%ymm\l

vmovshdup   %ymm13,%ymm13
vpblendd   $0xAA,%ymm14,%ymm13,%ymm13

vpaddd      %ymm13,%ymm12,%ymm\h
vpaddd      %ymm\h,%ymm3,%ymm\h
vpsubd      %ymm13,%ymm\l,%ymm\l
.endm

.macro butterfly l,h,zl0=1,zl1=1,zh0=2,zh1=2
vpmuldq		%ymm\zl0,%ymm\h,%ymm13   // zinvq * b = r0
vmovshdup	%ymm\h,%ymm12            // b >>32  = h
vpmuldq		%ymm\zl1,%ymm12,%ymm14 //zinvq * h = r1

vpmuldq		%ymm\zh0,%ymm\h,%ymm\h //zeta * b = b
vpmuldq		%ymm\zh1,%ymm12,%ymm12 //zeta * h = h

vpmuldq		%ymm0,%ymm13,%ymm13 // q * r0  = r0
vpmuldq		%ymm0,%ymm14,%ymm14 // r1 * q = r1

vpsubd		%ymm12,%ymm14,%ymm14
vpsubd		%ymm\h,%ymm13,%ymm13

vmovshdup	%ymm13,%ymm13
vpblendd	$0xAA,%ymm14,%ymm13,%ymm14

vpaddd		%ymm14,%ymm\l,%ymm\h
vpsubd		%ymm14,%ymm\l,%ymm\l
.endm


.macro butterfly7 l,h,zl0=1,zl1=1,zh0=2,zh1=2
vpmuldq		%ymm\zl0,%ymm\h,%ymm13   // zinvq * b = r0
vmovshdup	%ymm\h,%ymm12            // b >>32  = h

vmovshdup	%ymm\zl1,%ymm\zl1            // zinvq >>32  = zinvq
vpmuldq		%ymm\zl1,%ymm12,%ymm14 //zinvq * h = r1

vpmuldq		%ymm\zh0,%ymm\h,%ymm\h //zeta * b = b

vmovshdup	%ymm\zh1,%ymm\zh1            // zinvq >>32  = zinvq
vpmuldq		%ymm\zh1,%ymm12,%ymm12 //zeta * h = h

vpmuldq		%ymm0,%ymm13,%ymm13 // q * r0  = r0
vpmuldq		%ymm0,%ymm14,%ymm14 // r1 * q = r1

vpsubd		%ymm12,%ymm14,%ymm14
vpsubd		%ymm\h,%ymm13,%ymm13

vmovshdup	%ymm13,%ymm13
vpblendd	$0xAA,%ymm14,%ymm13,%ymm14

vpaddd		%ymm14,%ymm\l,%ymm\h
vpsubd		%ymm14,%ymm\l,%ymm\l
.endm

.macro levels0t1 off
/* level 0 */
vpbroadcastd	(_ZETAS_QINV+1)*4(%rsi),%ymm1
vpbroadcastd	(_ZETAS+1)*4(%rsi),%ymm2

vmovdqa		  0+64*\off(%rdi),%ymm4   
vmovdqa		32+64*\off(%rdi),%ymm5   
vmovdqa		2048+64*\off(%rdi),%ymm6   
vmovdqa	 	2080+64*\off(%rdi),%ymm7   
vmovdqa		4096+64*\off(%rdi),%ymm8   
vmovdqa		4128+64*\off(%rdi),%ymm9   
vmovdqa		6144+64*\off(%rdi),%ymm10  
vmovdqa	 	6176+64*\off(%rdi),%ymm11  

tri_butterfly	4,8
tri_butterfly	5,9
tri_butterfly	6,10
tri_butterfly	7,11

/* level 1 */
vpbroadcastd	(_ZETAS_QINV+2)*4(%rsi),%ymm1
vpbroadcastd	(_ZETAS+2)*4(%rsi),%ymm2
butterfly	4,6
butterfly	5,7

vpbroadcastd	(_ZETAS_QINV+3)*4(%rsi),%ymm1
vpbroadcastd	(_ZETAS+3)*4(%rsi),%ymm2
butterfly	8,10
butterfly	9,11

vmovdqa		%ymm4,  0+64*\off(%rdi)
vmovdqa		%ymm5,32+64*\off(%rdi)
vmovdqa		%ymm6,2048+64*\off(%rdi)
vmovdqa		%ymm7,2080+64*\off(%rdi)
vmovdqa		%ymm8,4096+64*\off(%rdi)
vmovdqa		%ymm9,4128+64*\off(%rdi)
vmovdqa		%ymm10,6144+64*\off(%rdi)
vmovdqa		%ymm11,6176+64*\off(%rdi)
.endm



.macro levels2t4 off , off2
/* level 0 */
vpbroadcastd	(_ZETAS_QINV+4+ 7*\off2)*4(%rsi),%ymm1
vpbroadcastd	(_ZETAS+4 +7*\off2)*4(%rsi),%ymm2

vmovdqa		  0+ 32*\off+2048*\off2(%rdi),%ymm4   
vmovdqa		256+32*\off+2048*\off2(%rdi),%ymm5   
vmovdqa		512+32*\off+2048*\off2(%rdi),%ymm6   
vmovdqa	 	768+32*\off+2048*\off2(%rdi),%ymm7   
vmovdqa	    1024+32*\off+2048*\off2(%rdi),%ymm8   
vmovdqa		1280+32*\off+2048*\off2(%rdi),%ymm9   
vmovdqa		1536+32*\off+2048*\off2(%rdi),%ymm10  
vmovdqa	 	1792+32*\off+2048*\off2(%rdi),%ymm11  

butterfly	4,8
butterfly	5,9
butterfly	6,10
butterfly	7,11

/* level 1 */
vpbroadcastd	(_ZETAS_QINV+5+7*\off2)*4(%rsi),%ymm1
vpbroadcastd	(_ZETAS+5+7*\off2)*4(%rsi),%ymm2
butterfly	4,6
butterfly	5,7

vpbroadcastd	(_ZETAS_QINV+6+7*\off2)*4(%rsi),%ymm1
vpbroadcastd	(_ZETAS+6+7*\off2)*4(%rsi),%ymm2
butterfly	8,10
butterfly	9,11

vpbroadcastd	(_ZETAS_QINV+7+7*\off2)*4(%rsi),%ymm1
vpbroadcastd	(_ZETAS+7+7*\off2)*4(%rsi),%ymm2
butterfly	4,5

vpbroadcastd	(_ZETAS_QINV+8+7*\off2)*4(%rsi),%ymm1
vpbroadcastd	(_ZETAS+8+7*\off2)*4(%rsi),%ymm2
butterfly	6,7

vpbroadcastd	(_ZETAS_QINV+9+7*\off2)*4(%rsi),%ymm1
vpbroadcastd	(_ZETAS+9+7*\off2)*4(%rsi),%ymm2
butterfly	8,9


vpbroadcastd	(_ZETAS_QINV+10+7*\off2)*4(%rsi),%ymm1
vpbroadcastd	(_ZETAS+10+7*\off2)*4(%rsi),%ymm2
butterfly	10,11

vmovdqa		%ymm4,  0+32*\off+2048*\off2(%rdi)
vmovdqa		%ymm5,256+32*\off+2048*\off2(%rdi)
vmovdqa		%ymm6,512+32*\off+2048*\off2(%rdi)
vmovdqa		%ymm7,768+32*\off+2048*\off2(%rdi)
vmovdqa		%ymm8,1024+32*\off+2048*\off2(%rdi)
vmovdqa		%ymm9,1280+32*\off+2048*\off2(%rdi)
vmovdqa		%ymm10,1536+32*\off+2048*\off2(%rdi)
vmovdqa		%ymm11,1792+32*\off+2048*\off2(%rdi)
.endm


.macro levels8t10 off 

vmovdqu	(_ZETAS_QINV+\off)*4(%rsi),%ymm1
vmovdqu	(_ZETAS+\off)*4(%rsi),%ymm2
shuffle8 4,5,3,5
butterfly 3,5

vmovdqu	(_ZETAS_QINV+\off+8)*4(%rsi),%ymm1
vmovdqu	(_ZETAS+\off+8)*4(%rsi),%ymm2
shuffle4 3,5,4,5
butterfly 4,5
vmovdqu	(_ZETAS_QINV+\off+16)*4(%rsi),%ymm1
vmovdqu	(_ZETAS+\off+16)*4(%rsi),%ymm2
shuffle2 4,5,3,5
butterfly7 3,5
shuffle_inv 3,5


vmovdqu	(_ZETAS_QINV+\off+24)*4(%rsi),%ymm1
vmovdqu	(_ZETAS+\off+24)*4(%rsi),%ymm2
shuffle8 6,7,4,7
butterfly 4,7

vmovdqu	(_ZETAS_QINV+\off+32)*4(%rsi),%ymm1
vmovdqu	(_ZETAS+\off+32)*4(%rsi),%ymm2
shuffle4 4,7,6,7
butterfly 6,7
vmovdqu	(_ZETAS_QINV+\off+40)*4(%rsi),%ymm1
vmovdqu	(_ZETAS+\off+40)*4(%rsi),%ymm2
shuffle2 6,7,4,7
butterfly7 4,7
shuffle_inv 4,7


vmovdqu	(_ZETAS_QINV+\off+48)*4(%rsi),%ymm1
vmovdqu	(_ZETAS+\off+48)*4(%rsi),%ymm2
shuffle8 8,9,6,9
butterfly 6,9

vmovdqu	(_ZETAS_QINV+\off+56)*4(%rsi),%ymm1
vmovdqu	(_ZETAS+\off+56)*4(%rsi),%ymm2
shuffle4 6,9,8,9
butterfly 8,9
vmovdqu	(_ZETAS_QINV+\off+64)*4(%rsi),%ymm1
vmovdqu	(_ZETAS+\off+64)*4(%rsi),%ymm2
shuffle2 8,9,6,9
butterfly7 6,9
shuffle_inv 6,9

vmovdqu	(_ZETAS_QINV+\off+72)*4(%rsi),%ymm1
vmovdqu	(_ZETAS+\off+72)*4(%rsi),%ymm2
shuffle8 10,11,8,11
butterfly 8,11

vmovdqu	(_ZETAS_QINV+\off+80)*4(%rsi),%ymm1
vmovdqu	(_ZETAS+\off+80)*4(%rsi),%ymm2
shuffle4 8,11,10,11
butterfly 10,11
vmovdqu	(_ZETAS_QINV+\off+88)*4(%rsi),%ymm1
vmovdqu	(_ZETAS+\off+88)*4(%rsi),%ymm2
shuffle2 10,11,8,11
butterfly7 8,11
shuffle_inv 8,11

.endm



.macro levels5t10 off , off2 // 32 , 0
vpbroadcastd	(_ZETAS_QINV+\off)*4(%rsi),%ymm1
vpbroadcastd	(_ZETAS+\off)*4(%rsi),%ymm2

vmovdqa		  0+ 256*\off2(%rdi),%ymm4   
vmovdqa		32+256*\off2(%rdi),%ymm5   
vmovdqa		64+256*\off2(%rdi),%ymm6   
vmovdqa	 	96+256*\off2(%rdi),%ymm7   
vmovdqa	    128+256*\off2(%rdi),%ymm8   
vmovdqa		160+256*\off2(%rdi),%ymm9   
vmovdqa		192+256*\off2(%rdi),%ymm10  
vmovdqa	 	224+256*\off2(%rdi),%ymm11  


butterfly	4,8
butterfly	5,9
butterfly	6,10
butterfly	7,11

vpbroadcastd	(_ZETAS_QINV+\off+1)*4(%rsi),%ymm1
vpbroadcastd	(_ZETAS+\off+1)*4(%rsi),%ymm2


butterfly	4,6
butterfly	5,7

vpbroadcastd	(_ZETAS_QINV+\off+2)*4(%rsi),%ymm1
vpbroadcastd	(_ZETAS+\off+2)*4(%rsi),%ymm2


butterfly	8,10
butterfly	9,11

vpbroadcastd	(_ZETAS_QINV+\off+3)*4(%rsi),%ymm1
vpbroadcastd	(_ZETAS+\off+3)*4(%rsi),%ymm2

butterfly	4,5

vpbroadcastd	(_ZETAS_QINV+\off+4)*4(%rsi),%ymm1
vpbroadcastd	(_ZETAS+\off+4)*4(%rsi),%ymm2

butterfly	6,7

vpbroadcastd	(_ZETAS_QINV+\off+5)*4(%rsi),%ymm1
vpbroadcastd	(_ZETAS+\off+5)*4(%rsi),%ymm2

butterfly	8,9

vpbroadcastd	(_ZETAS_QINV+\off+6)*4(%rsi),%ymm1
vpbroadcastd	(_ZETAS+\off+6)*4(%rsi),%ymm2

butterfly	10,11

levels8t10 \off+7



vmovdqa		%ymm3,  0+256*\off2(%rdi)
vmovdqa		%ymm5,32+256*\off2(%rdi)
vmovdqa		%ymm4,64+256*\off2(%rdi)
vmovdqa		%ymm7,96+256*\off2(%rdi)
vmovdqa		%ymm6,128+256*\off2(%rdi)
vmovdqa		%ymm9,160+256*\off2(%rdi)
vmovdqa		%ymm8,192+256*\off2(%rdi)
vmovdqa		%ymm11,224+256*\off2(%rdi)
.endm






.text
.global ntt5_asm
ntt5_asm:
vmovdqa		_8XQ*4(%rsi),%ymm0

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

levels5t10 32,0
levels5t10 135,1
levels5t10 238,2
levels5t10 341,3
levels5t10 444,4
levels5t10 547,5
levels5t10 650,6
levels5t10 753,7


levels5t10 856,8
levels5t10 959,9
levels5t10 1062,10
levels5t10 1165,11
levels5t10 1268,12
levels5t10 1371,13
levels5t10 1474,14
levels5t10 1577,15


levels5t10 1680,16
levels5t10 1783,17
levels5t10 1886,18
levels5t10 1989,19
levels5t10 2092,20
levels5t10 2195,21
levels5t10 2298,22
levels5t10 2401,23

levels5t10 2504,24
levels5t10 2607,25
levels5t10 2710,26
levels5t10 2813,27
levels5t10 2916,28
levels5t10 3019,29
levels5t10 3122,30
levels5t10 3225,31



ret