.section .text /* text section */
.syntax unified /* Unified Assembly Syntax - UAL */
.thumb /* Thumb instruction set */

#include "macros.i"

// void asm_ntt_1(int32_t * Out, int32_t* zeta);
.global asm_ntt_1
.type asm_ntt_1,%function
.align 2
asm_ntt_1:
    push.w {r4-r11, r14}
    ptr_p     .req R0
  	ptr_zeta  .req R1
  	zeta      .req R1
  	qinv      .req R2
  	q         .req R3
  	cntr      .req R4
  	pol4      .req R4
  	pol0      .req R5
  	pol1      .req R6
  	pol2      .req R7
  	pol3      .req R8
  	temp_h    .req R9
  	temp_l    .req R10
  	zeta0     .req R11
  	zeta1     .req R12
  	zeta2     .req R14
  	pol5     .req R11
  	pol6     .req R12
  	pol7     .req R14

    movw qinv, #:lower16:0x7a29f27f
    movt qinv, #:upper16:0x7a29f27f
    movw q, #0x3281
    movt q, #0x80



//stage 1 - 3
	.equ distance, 144 //144*4*4
  	.equ strincr, 4

  	vldm ptr_zeta!, {s2-s8}
  	vmov s0, ptr_zeta

  	add.w temp_l, ptr_p, #distance*strincr //  8*144 = 1152
  	vmov s9, temp_l

  1:
    .rept 2
    ldr.w pol0, [ptr_p]
    ldr.w pol1, [ptr_p, #1*distance*strincr]
    ldr.w pol2, [ptr_p, #2*distance*strincr]
    ldr.w pol3, [ptr_p, #3*distance*strincr]
    ldr.w pol4, [ptr_p, #4*distance*strincr]
    ldr.w pol5, [ptr_p, #5*distance*strincr]
    ldr.w pol6, [ptr_p, #6*distance*strincr]
    ldr.w pol7, [ptr_p, #7*distance*strincr]

    _first_3_layer_CT_32 pol0, pol1, pol2, pol3, pol4, pol5, pol6, pol7, s2, s3, s4, s5, s6, s7, s8, zeta, qinv, q, temp_h, temp_l


    str.w pol1, [ptr_p, #1*distance*strincr]
    str.w pol2, [ptr_p, #2*distance*strincr]
    str.w pol3, [ptr_p, #3*distance*strincr]
    str.w pol4, [ptr_p, #4*distance*strincr]
    str.w pol5, [ptr_p, #5*distance*strincr]
    str.w pol6, [ptr_p, #6*distance*strincr]
    str.w pol7, [ptr_p, #7*distance*strincr]
    str.w pol0, [ptr_p], #strincr
    .endr

	vmov temp_l, s9
	cmp.w ptr_p, temp_l
    bne 1b

    sub ptr_p, #distance*strincr

//stage 4 - 5
	.equ distance, 36	//36*4*4
	add.w temp_l, ptr_p, #4608
  	vmov s9, temp_l

  1:
	add.w temp_l, ptr_p, #36*strincr
    vmov s10, temp_l
    vmov ptr_zeta, s0
    ldr.w zeta1, [ptr_zeta, #4]
    ldr.w zeta2, [ptr_zeta, #8]
    ldr.w zeta0, [ptr_zeta], #12
    vmov s0, ptr_zeta

	  2:
		.rept 2
		ldr.w pol0, [ptr_p]
    	ldr.w pol1, [ptr_p, #1*distance*strincr]
    	ldr.w pol2, [ptr_p, #2*distance*strincr]
    	ldr.w pol3, [ptr_p, #3*distance*strincr]

    	_2_layer_CT_32 pol0, pol1, pol2, pol3, zeta0, zeta1, zeta2, qinv, q, temp_h, temp_l

		str.w pol1, [ptr_p, #1*distance*strincr]
    	str.w pol2, [ptr_p, #2*distance*strincr]
    	str.w pol3, [ptr_p, #3*distance*strincr]
    	str.w pol0, [ptr_p], #strincr
		.endr

		vmov temp_l, s10
		cmp.w ptr_p, temp_l
    	bne 2b

	add.w ptr_p, #108*strincr //next ring index-0 : 36 + 108 = 144
    vmov temp_l, s9
    cmp.w ptr_p, temp_l
    bne 1b

//stage 6 - 7

	sub ptr_p, #4608

	.equ distance, 9	//9*4*4

	add.w temp_l, ptr_p, #4608
  	vmov s9, temp_l

  1:
	add.w temp_l, ptr_p, #9*strincr
    vmov s10, temp_l
    vmov ptr_zeta, s0
    ldr.w zeta1, [ptr_zeta, #4]
    ldr.w zeta2, [ptr_zeta, #8]
    ldr.w zeta0, [ptr_zeta], #12
    vmov s0, ptr_zeta

	  2:
		.rept 3
		ldr.w pol0, [ptr_p]
    	ldr.w pol1, [ptr_p, #1*distance*strincr]
    	ldr.w pol2, [ptr_p, #2*distance*strincr]
    	ldr.w pol3, [ptr_p, #3*distance*strincr]

    	_2_layer_CT_32 pol0, pol1, pol2, pol3, zeta0, zeta1, zeta2, qinv, q, temp_h, temp_l

		str.w pol1, [ptr_p, #1*distance*strincr]
    	str.w pol2, [ptr_p, #2*distance*strincr]
    	str.w pol3, [ptr_p, #3*distance*strincr]
    	str.w pol0, [ptr_p], #strincr
		.endr

		vmov temp_l, s10
		cmp.w ptr_p, temp_l
    	bne 2b

	add.w ptr_p, #27*strincr //next ring index-0 : 9 + 27 = 36
    vmov temp_l, s9
    cmp.w ptr_p, temp_l
    bne 1b

    pop.w {r4-r11, pc}


//extern void asm_ntt_1_radix3(int32_t * Out, int32_t* zeta);
.global asm_ntt_1_radix3
.type asm_ntt_1_radix3,%function
.align 2
asm_ntt_1_radix3:
    push.w {r4-r11, r14}
    vpush.w {s16-s24}
    ptr_p     .req R0
  	ptr_zeta  .req R1
  	qinv      .req R2
  	q         .req R3
  	cntr      .req R4
  	pol0      .req R5
  	pol1      .req R6
  	pol2      .req R7
  	pol3      .req R8
  	temp_h    .req R9
  	temp_l    .req R10
  	zeta0     .req R11
  	zeta1     .req R12
  	wmont     .req R14


    movw qinv, #:lower16:0x7a29f27f
    movt qinv, #:upper16:0x7a29f27f
    movw q, #0x3281
    movt q, #0x80
    movw wmont, #0x571f // 1857311
    movt wmont, #0x1c // 1857311


	movw cntr, #128

1:
	vldm ptr_p, {s0-s8} //2-8
	vldm ptr_zeta!, {s9-s16}

	_2_layer_CT_32_radix3 pol0, pol1, pol2, pol3, s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, zeta0, zeta1, qinv, q, wmont, temp_l, temp_h

	vmov.w pol0, s1
	vmov.w pol1, s2
	vmov.w pol2, s3
    str.w pol0, [ptr_p, #4]
    str.w pol1, [ptr_p, #8]
    str.w pol2, [ptr_p, #12]

    vmov.w pol0, s4
	vmov.w pol1, s5
	vmov.w pol2, s6
    str.w pol0, [ptr_p, #16]
    str.w pol1, [ptr_p, #20]
    str.w pol2, [ptr_p, #24]

    vmov.w pol0, s7
	vmov.w pol1, s8
	vmov.w pol2, s0
    str.w pol0, [ptr_p, #28]
    str.w pol1, [ptr_p, #32]
    str.w pol2, [ptr_p], #36

	subs cntr, #1
    bne.w 1b

	vpop.w {s16-s24}
    pop.w {r4-r11, pc}


// void asm_ntt_3(int32_t * Out, int32_t* zeta);
.global asm_ntt_3
.type asm_ntt_3,%function
.align 2
asm_ntt_3:
    push.w {r4-r11, r14}
    ptr_p     .req R0
  	ptr_zeta  .req R1
  	zeta      .req R1
  	qinv      .req R2
  	q         .req R3
  	cntr      .req R4
  	pol4      .req R4
  	pol0      .req R5
  	pol1      .req R6
  	pol2      .req R7
  	pol3      .req R8
  	temp_h    .req R9
  	temp_l    .req R10
  	zeta0     .req R11
  	zeta1     .req R12
  	zeta2     .req R14
  	pol5     .req R11
  	pol6     .req R12
  	pol7     .req R14

	movw qinv, #:lower16:0x83fc21ff
    movt qinv, #:upper16:0x83fc21ff
    movw q, #0x2201
    movt q, #0x80

//stage 1 - 3
	.equ distance, 192 //1536/8
  	.equ strincr, 4

  	vldm ptr_zeta!, {s2-s8}
  	vmov s0, ptr_zeta

  	add.w temp_l, ptr_p, #distance*strincr // 8*192 = 1536
  	vmov s9, temp_l

	add.w temp_l, ptr_p, #4*distance*strincr
	vmov s14, temp_l
  1:
    .rept 2
    vmov temp_l, s14
    ldr.w pol0, [ptr_p]
    ldr.w pol1, [ptr_p, #1*distance*strincr]
    ldr.w pol2, [ptr_p, #2*distance*strincr]
    ldr.w pol3, [ptr_p, #3*distance*strincr]
    ldr.w pol4, [temp_l]
    ldr.w pol5, [temp_l, #1*distance*strincr]
    ldr.w pol6, [temp_l, #2*distance*strincr]
    ldr.w pol7, [temp_l, #3*distance*strincr]

    _first_3_layer_CT_32 pol0, pol1, pol2, pol3, pol4, pol5, pol6, pol7, s2, s3, s4, s5, s6, s7, s8, zeta, qinv, q, temp_h, temp_l

	vmov temp_l, s14
    str.w pol1, [ptr_p, #1*distance*strincr]
    str.w pol2, [ptr_p, #2*distance*strincr]
    str.w pol3, [ptr_p, #3*distance*strincr]
    str.w pol0, [ptr_p], #strincr
    str.w pol5, [temp_l, #1*distance*strincr]
    str.w pol6, [temp_l, #2*distance*strincr]
    str.w pol7, [temp_l, #3*distance*strincr]
    str.w pol4, [temp_l], #strincr

	vmov s14, temp_l
    .endr

	vmov temp_l, s9
	cmp.w ptr_p, temp_l
    bne 1b

    sub ptr_p, #distance*strincr


//stage 4 - 6
	.equ distance, 24	//24*8 = 192
	add.w temp_l, ptr_p, #6144
  	vmov s9, temp_l

  1:
	add.w temp_l, ptr_p, #distance*strincr
    vmov s10, temp_l
    vmov ptr_zeta, s0
    vldm ptr_zeta!, {s2-s8}
    vmov s0, ptr_zeta

	  2:
		.rept 2
		ldr.w pol0, [ptr_p]
    	ldr.w pol1, [ptr_p, #1*distance*strincr]
    	ldr.w pol2, [ptr_p, #2*distance*strincr]
    	ldr.w pol3, [ptr_p, #3*distance*strincr]
    	ldr.w pol4, [ptr_p, #4*distance*strincr]
    	ldr.w pol5, [ptr_p, #5*distance*strincr]
    	ldr.w pol6, [ptr_p, #6*distance*strincr]
    	ldr.w pol7, [ptr_p, #7*distance*strincr]

    	_3_layer_CT_32 pol0, pol1, pol2, pol3, pol4, pol5, pol6, pol7, s2, s3, s4, s5, s6, s7, s8, zeta, qinv, q, temp_h, temp_l


    	str.w pol1, [ptr_p, #1*distance*strincr]
    	str.w pol2, [ptr_p, #2*distance*strincr]
    	str.w pol3, [ptr_p, #3*distance*strincr]
    	str.w pol4, [ptr_p, #4*distance*strincr]
    	str.w pol5, [ptr_p, #5*distance*strincr]
    	str.w pol6, [ptr_p, #6*distance*strincr]
    	str.w pol7, [ptr_p, #7*distance*strincr]
    	str.w pol0, [ptr_p], #strincr
		.endr

		vmov temp_l, s10
		cmp.w ptr_p, temp_l
    	bne 2b

	add.w ptr_p, #168*strincr //next ring index-0 : 24 + 168 = 192
    vmov temp_l, s9
    cmp.w ptr_p, temp_l
    bne 1b

    sub ptr_p, #6144

//stage 7 - 9
	.equ distance, 3	//3*8 = 24
	add.w temp_l, ptr_p, #6144
  	vmov s9, temp_l

  1:
	add.w temp_l, ptr_p, #distance*strincr
    vmov s10, temp_l
    vmov ptr_zeta, s0
    vldm ptr_zeta!, {s2-s8}
    vmov s0, ptr_zeta

	  2:
		ldr.w pol0, [ptr_p]
    	ldr.w pol1, [ptr_p, #1*distance*strincr]
    	ldr.w pol2, [ptr_p, #2*distance*strincr]
    	ldr.w pol3, [ptr_p, #3*distance*strincr]
    	ldr.w pol4, [ptr_p, #4*distance*strincr]
    	ldr.w pol5, [ptr_p, #5*distance*strincr]
    	ldr.w pol6, [ptr_p, #6*distance*strincr]
    	ldr.w pol7, [ptr_p, #7*distance*strincr]

    	 _3_layer_CT_32 pol0, pol1, pol2, pol3, pol4, pol5, pol6, pol7, s2, s3, s4, s5, s6, s7, s8, zeta, qinv, q, temp_h, temp_l


    	str.w pol1, [ptr_p, #1*distance*strincr]
    	str.w pol2, [ptr_p, #2*distance*strincr]
    	str.w pol3, [ptr_p, #3*distance*strincr]
    	str.w pol4, [ptr_p, #4*distance*strincr]
    	str.w pol5, [ptr_p, #5*distance*strincr]
    	str.w pol6, [ptr_p, #6*distance*strincr]
    	str.w pol7, [ptr_p, #7*distance*strincr]
    	str.w pol0, [ptr_p], #strincr

		vmov temp_l, s10
		cmp.w ptr_p, temp_l
    	bne 2b

	add.w ptr_p, #21*strincr //next ring index-0 : 3 + 21 = 24
    vmov temp_l, s9
    cmp.w ptr_p, temp_l
    bne 1b

    pop.w {r4-r11, pc}



// void asm_ntt_55(int32_t * Out, int32_t* zeta);
.global asm_ntt_55
.type asm_ntt_55,%function
.align 2
asm_ntt_55:
    push.w {r4-r11, r14}
    ptr_p     .req R0
  	ptr_zeta  .req R1
  	zeta      .req R1
  	qinv      .req R2
  	q         .req R3
  	cntr      .req R4
  	pol4      .req R4
  	pol0      .req R5
  	pol1      .req R6
  	pol2      .req R7
  	pol3      .req R8
  	temp_h    .req R9
  	temp_l    .req R10
  	zeta0     .req R11
  	zeta1     .req R12
  	zeta2     .req R14
  	pol5     .req R11
  	pol6     .req R12
  	pol7     .req R14

	movw qinv, #:lower16:0xf0803fff
    movt qinv, #:upper16:0xf0803fff
    movw q, #0x4001
    movt q, #0x80

//stage 1 - 3
	.equ distance, 288 //2304/8
  	.equ strincr, 4

  	vldm ptr_zeta!, {s2-s8}
  	vmov s0, ptr_zeta

  	add.w temp_l, ptr_p, #distance*strincr //8*288 = 1536
  	vmov s9, temp_l

	add.w temp_l, ptr_p, #4*distance*strincr
	vmov s14, temp_l
  1:
    .rept 2
    vmov temp_l, s14
    ldr.w pol0, [ptr_p]
    ldr.w pol1, [ptr_p, #1*distance*strincr]
    ldr.w pol2, [ptr_p, #2*distance*strincr]
    ldr.w pol3, [ptr_p, #3*distance*strincr]
    ldr.w pol4, [temp_l]
    ldr.w pol5, [temp_l, #1*distance*strincr]
    ldr.w pol6, [temp_l, #2*distance*strincr]
    ldr.w pol7, [temp_l, #3*distance*strincr]

    _first_3_layer_CT_32 pol0, pol1, pol2, pol3, pol4, pol5, pol6, pol7, s2, s3, s4, s5, s6, s7, s8, zeta, qinv, q, temp_h, temp_l

	vmov temp_l, s14
    str.w pol1, [ptr_p, #1*distance*strincr]
    str.w pol2, [ptr_p, #2*distance*strincr]
    str.w pol3, [ptr_p, #3*distance*strincr]
    str.w pol0, [ptr_p], #strincr
    str.w pol5, [temp_l, #1*distance*strincr]
    str.w pol6, [temp_l, #2*distance*strincr]
    str.w pol7, [temp_l, #3*distance*strincr]
    str.w pol4, [temp_l], #strincr

	vmov s14, temp_l
    .endr

	vmov temp_l, s9
	cmp.w ptr_p, temp_l
    bne 1b

    sub ptr_p, #distance*strincr


//stage 4 - 6
	.equ distance, 36	//288/8 = 36
	add.w temp_l, ptr_p, #9216
  	vmov s9, temp_l

  1:
	add.w temp_l, ptr_p, #distance*strincr
    vmov s10, temp_l
    vmov ptr_zeta, s0
    vldm ptr_zeta!, {s2-s8}
    vmov s0, ptr_zeta

	  2:
		.rept 2
		ldr.w pol0, [ptr_p]
    	ldr.w pol1, [ptr_p, #1*distance*strincr]
    	ldr.w pol2, [ptr_p, #2*distance*strincr]
    	ldr.w pol3, [ptr_p, #3*distance*strincr]
    	ldr.w pol4, [ptr_p, #4*distance*strincr]
    	ldr.w pol5, [ptr_p, #5*distance*strincr]
    	ldr.w pol6, [ptr_p, #6*distance*strincr]
    	ldr.w pol7, [ptr_p, #7*distance*strincr]

    	_3_layer_CT_32 pol0, pol1, pol2, pol3, pol4, pol5, pol6, pol7, s2, s3, s4, s5, s6, s7, s8, zeta, qinv, q, temp_h, temp_l


    	str.w pol1, [ptr_p, #1*distance*strincr]
    	str.w pol2, [ptr_p, #2*distance*strincr]
    	str.w pol3, [ptr_p, #3*distance*strincr]
    	str.w pol4, [ptr_p, #4*distance*strincr]
    	str.w pol5, [ptr_p, #5*distance*strincr]
    	str.w pol6, [ptr_p, #6*distance*strincr]
    	str.w pol7, [ptr_p, #7*distance*strincr]
    	str.w pol0, [ptr_p], #strincr
		.endr

		vmov temp_l, s10
		cmp.w ptr_p, temp_l
    	bne 2b

	add.w ptr_p, #252*strincr //next ring index-0 : 36 + 252 = 288
    vmov temp_l, s9
    cmp.w ptr_p, temp_l
    bne 1b

    sub ptr_p, #9216

//stage 7 - 8

	.equ distance, 9	//36/4 = 9

	add.w temp_l, ptr_p, #9216
  	vmov s9, temp_l

  1:
	add.w temp_l, ptr_p, #distance*strincr
    vmov s10, temp_l
    vmov ptr_zeta, s0
    ldr.w zeta1, [ptr_zeta, #4]
    ldr.w zeta2, [ptr_zeta, #8]
    ldr.w zeta0, [ptr_zeta], #12
    vmov s0, ptr_zeta

	  2:
		.rept 3
		ldr.w pol0, [ptr_p]
    	ldr.w pol1, [ptr_p, #1*distance*strincr]
    	ldr.w pol2, [ptr_p, #2*distance*strincr]
    	ldr.w pol3, [ptr_p, #3*distance*strincr]

    	_2_layer_CT_32 pol0, pol1, pol2, pol3, zeta0, zeta1, zeta2, qinv, q, temp_h, temp_l

		str.w pol1, [ptr_p, #1*distance*strincr]
    	str.w pol2, [ptr_p, #2*distance*strincr]
    	str.w pol3, [ptr_p, #3*distance*strincr]
    	str.w pol0, [ptr_p], #strincr
		.endr

		vmov temp_l, s10
		cmp.w ptr_p, temp_l
    	bne 2b

	add.w ptr_p, #27*strincr //next ring index-0 : 9 + 27 = 36
    vmov temp_l, s9
    cmp.w ptr_p, temp_l
    bne 1b

    pop.w {r4-r11, pc}

//extern void asm_ntt_55_radix3(int32_t * Out, int32_t* zeta);
.global asm_ntt_55_radix3
.type asm_ntt_55_radix3,%function
.align 2
asm_ntt_55_radix3:
    push.w {r4-r11, r14}
    vpush.w {s16-s24}
    ptr_p     .req R0
  	ptr_zeta  .req R1
  	qinv      .req R2
  	q         .req R3
  	cntr      .req R4
  	pol0      .req R5
  	pol1      .req R6
  	pol2      .req R7
  	pol3      .req R8
  	temp_h    .req R9
  	temp_l    .req R10
  	zeta0     .req R11
  	zeta1     .req R12
  	wmont     .req R14


    movw qinv, #:lower16:0xf0803fff
    movt qinv, #:upper16:0xf0803fff
    movw q, #0x4001
    movt q, #0x80
    movw wmont, #0xdc5d // 6675549
    movt wmont, #0x65 // 6675549


	movw cntr, #256

1:
	vldm ptr_p, {s0-s8} //2-8
	vldm ptr_zeta!, {s9-s16}

	_2_layer_CT_32_radix3 pol0, pol1, pol2, pol3, s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, zeta0, zeta1, qinv, q, wmont, temp_l, temp_h

	vmov.w pol0, s1
	vmov.w pol1, s2
	vmov.w pol2, s3
    str.w pol0, [ptr_p, #4]
    str.w pol1, [ptr_p, #8]
    str.w pol2, [ptr_p, #12]

    vmov.w pol0, s4
	vmov.w pol1, s5
	vmov.w pol2, s6
    str.w pol0, [ptr_p, #16]
    str.w pol1, [ptr_p, #20]
    str.w pol2, [ptr_p, #24]

    vmov.w pol0, s7
	vmov.w pol1, s8
	vmov.w pol2, s0
    str.w pol0, [ptr_p, #28]
    str.w pol1, [ptr_p, #32]
    str.w pol2, [ptr_p], #36

	subs cntr, #1
    bne.w 1b

	vpop.w {s16-s24}
    pop.w {r4-r11, pc}

// void asm_ntt_55(int32_t * Out, int32_t* zeta);
.global asm_ntt_5
.type asm_ntt_5,%function
.align 2
asm_ntt_5:
    push.w {r4-r11, r14}
    ptr_p     .req R0
  	ptr_zeta  .req R1
  	zeta      .req R1
  	qinv      .req R2
  	q         .req R3
  	cntr      .req R4
  	pol4      .req R4
  	pol0      .req R5
  	pol1      .req R6
  	pol2      .req R7
  	pol3      .req R8
  	temp_h    .req R9
  	temp_l    .req R10
  	zeta0     .req R11
  	zeta1     .req R12
  	zeta2     .req R14
  	pol5     .req R11
  	pol6     .req R12
  	pol7     .req R14

	movw qinv, #:lower16:0xfc7fdfff
    movt qinv, #:upper16:0xfc7fdfff
    movw q, #0xe001
    movt q, #0x7f

//stage 1 - 3
	.equ distance, 256 //2304/8
  	.equ strincr, 4

  	vldm ptr_zeta!, {s2-s8}
  	vmov s0, ptr_zeta

  	add.w temp_l, ptr_p, #distance*strincr //8*288 = 1536
  	vmov s9, temp_l

	add.w temp_l, ptr_p, #4*distance*strincr
	vmov s14, temp_l
  1:
    .rept 2
    vmov temp_l, s14
    ldr.w pol0, [ptr_p]
    ldr.w pol1, [ptr_p, #1*distance*strincr]
    ldr.w pol2, [ptr_p, #2*distance*strincr]
    ldr.w pol3, [ptr_p, #3*distance*strincr]
    ldr.w pol4, [temp_l]
    ldr.w pol5, [temp_l, #1*distance*strincr]
    ldr.w pol6, [temp_l, #2*distance*strincr]
    ldr.w pol7, [temp_l, #3*distance*strincr]

    _first_3_layer_CT_32 pol0, pol1, pol2, pol3, pol4, pol5, pol6, pol7, s2, s3, s4, s5, s6, s7, s8, zeta, qinv, q, temp_h, temp_l

	vmov temp_l, s14
    str.w pol1, [ptr_p, #1*distance*strincr]
    str.w pol2, [ptr_p, #2*distance*strincr]
    str.w pol3, [ptr_p, #3*distance*strincr]
    str.w pol0, [ptr_p], #strincr
    str.w pol5, [temp_l, #1*distance*strincr]
    str.w pol6, [temp_l, #2*distance*strincr]
    str.w pol7, [temp_l, #3*distance*strincr]
    str.w pol4, [temp_l], #strincr

	vmov s14, temp_l
    .endr

	vmov temp_l, s9
	cmp.w ptr_p, temp_l
    bne 1b

    sub ptr_p, #distance*strincr

	//stage 4 - 6
	.equ distance, 32	//256/8 = 32
	add.w temp_l, ptr_p, #8192
  	vmov s9, temp_l

  1:
	add.w temp_l, ptr_p, #distance*strincr
    vmov s10, temp_l
    vmov ptr_zeta, s0
    vldm ptr_zeta!, {s2-s8}
    vmov s0, ptr_zeta

	  2:
		.rept 2
		ldr.w pol0, [ptr_p]
    	ldr.w pol1, [ptr_p, #1*distance*strincr]
    	ldr.w pol2, [ptr_p, #2*distance*strincr]
    	ldr.w pol3, [ptr_p, #3*distance*strincr]
    	ldr.w pol4, [ptr_p, #4*distance*strincr]
    	ldr.w pol5, [ptr_p, #5*distance*strincr]
    	ldr.w pol6, [ptr_p, #6*distance*strincr]
    	ldr.w pol7, [ptr_p, #7*distance*strincr]

    	_3_layer_CT_32 pol0, pol1, pol2, pol3, pol4, pol5, pol6, pol7, s2, s3, s4, s5, s6, s7, s8, zeta, qinv, q, temp_h, temp_l


    	str.w pol1, [ptr_p, #1*distance*strincr]
    	str.w pol2, [ptr_p, #2*distance*strincr]
    	str.w pol3, [ptr_p, #3*distance*strincr]
    	str.w pol4, [ptr_p, #4*distance*strincr]
    	str.w pol5, [ptr_p, #5*distance*strincr]
    	str.w pol6, [ptr_p, #6*distance*strincr]
    	str.w pol7, [ptr_p, #7*distance*strincr]
    	str.w pol0, [ptr_p], #strincr
		.endr

		vmov temp_l, s10
		cmp.w ptr_p, temp_l
    	bne 2b

	add.w ptr_p, #224*strincr //next ring index-0 : 36 + 252 = 288
    vmov temp_l, s9
    cmp.w ptr_p, temp_l
    bne 1b

    sub ptr_p, #8192

	//stage 7 - 9
	.equ distance, 4	//32/8 = 4
	add.w temp_l, ptr_p, #8192
  	vmov s9, temp_l

  1:
	add.w temp_l, ptr_p, #distance*strincr
    vmov s10, temp_l
    vmov ptr_zeta, s0
    vldm ptr_zeta!, {s2-s8}
    vmov s0, ptr_zeta

	  2:
		.rept 2
		ldr.w pol0, [ptr_p]
    	ldr.w pol1, [ptr_p, #1*distance*strincr]
    	ldr.w pol2, [ptr_p, #2*distance*strincr]
    	ldr.w pol3, [ptr_p, #3*distance*strincr]
    	ldr.w pol4, [ptr_p, #4*distance*strincr]
    	ldr.w pol5, [ptr_p, #5*distance*strincr]
    	ldr.w pol6, [ptr_p, #6*distance*strincr]
    	ldr.w pol7, [ptr_p, #7*distance*strincr]

    	_3_layer_CT_32 pol0, pol1, pol2, pol3, pol4, pol5, pol6, pol7, s2, s3, s4, s5, s6, s7, s8, zeta, qinv, q, temp_h, temp_l


    	str.w pol1, [ptr_p, #1*distance*strincr]
    	str.w pol2, [ptr_p, #2*distance*strincr]
    	str.w pol3, [ptr_p, #3*distance*strincr]
    	str.w pol4, [ptr_p, #4*distance*strincr]
    	str.w pol5, [ptr_p, #5*distance*strincr]
    	str.w pol6, [ptr_p, #6*distance*strincr]
    	str.w pol7, [ptr_p, #7*distance*strincr]
    	str.w pol0, [ptr_p], #strincr
		.endr

		vmov temp_l, s10
		cmp.w ptr_p, temp_l
    	bne 2b

	add.w ptr_p, #28*strincr //next ring index-0 : 4 + 28 = 32
    vmov temp_l, s9
    cmp.w ptr_p, temp_l
    bne 1b

    sub ptr_p, #8192

    //stage 10 - 11

	.equ distance, 1	//36/4 = 9

	add.w temp_l, ptr_p, #8192
  	vmov s9, temp_l

  1:

    vmov ptr_zeta, s0
    ldr.w zeta1, [ptr_zeta, #4]
    ldr.w zeta2, [ptr_zeta, #8]
    ldr.w zeta0, [ptr_zeta], #12
    vmov s0, ptr_zeta

	ldr.w pol0, [ptr_p]
    ldr.w pol1, [ptr_p, #1*distance*strincr]
    ldr.w pol2, [ptr_p, #2*distance*strincr]
    ldr.w pol3, [ptr_p, #3*distance*strincr]

   	_2_layer_CT_32 pol0, pol1, pol2, pol3, zeta0, zeta1, zeta2, qinv, q, temp_h, temp_l

	str.w pol1, [ptr_p, #1*distance*strincr]
   	str.w pol2, [ptr_p, #2*distance*strincr]
   	str.w pol3, [ptr_p, #3*distance*strincr]
   	str.w pol0, [ptr_p], #strincr

	add.w ptr_p, #3*strincr
    vmov temp_l, s9
    cmp.w ptr_p, temp_l
    bne 1b

    pop.w {r4-r11, pc}

.end /* End of file */
