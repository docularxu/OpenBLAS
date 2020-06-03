/*******************************************************************************
Copyright (c) 2015, The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*******************************************************************************/

#define ASSEMBLER
#include "common.h"

/*                   X0          X1          X2          s0        X3        x4       x5           x6            x7*/
/*int CNAME(BLASLONG bm,BLASLONG bn,BLASLONG bk,FLOAT alpha0,FLOAT* ba,FLOAT* bb,FLOAT* C,BLASLONG ldc, BLASLONG offset) */

#define origM		x0
#define origN		x1
#define origK		x2
#define origPA		x3
#define origPB		x4
#define pC		x5
#define LDC		x6
#define offset		x7
#define counterL	x8
#define counterI	x9
#define counterJ	x10
#define pB		x11
#define pCRow0		x12
#define pCRow1		x13
#define pCRow2		x14
#define pCRow3		x15
#define pA		x16
#define alpha		x17
#define temp		x18
#define tempOffset	x19
#define tempK		x20
#define pCRow4		x21
#define pCRow5		x22
#define pCRow6		x23
#define pCRow7		x24

#define alpha0		d10
#define alphaV0		v10.d[0]

#define A_PRE_SIZE	2560
#define B_PRE_SIZE	448
#define C_PRE_SIZE	128

// 00 origM
// 01 origN
// 02 origK
// 03 origPA
// 04 origPB
// 05 pC
// 06 origLDC -> LDC
// 07 offset
// 08 counterL
// 09 counterI
// 10 counterJ
// 11 pB
// 12 pCRow0
// 13 pCRow1
// 14 pCRow2
// 15 pA
// 16 temp
// 17 tempOffset
// 18 must save tempK
// 19 must save
// 20 must save
// 21 must save
// 22 must save
// 23 must save
// 24 must save
// 25 must save
// 26 must save
// 27 must save
// 28 must save
// 29 frame
// 30 link
// 31 sp

//v00 ALPHA -> pA0_0, pA0_1
//v01 pA0_2, pA0_3
//v02 pA0_4, pA0_5
//v03 pA0_6, pA0_7
//v04 pA1_0, pA1_1
//v05 pA1_2, pA1_3
//v06 pA1_4, pA1_5
//v07 pA1_6, pA1_7
//v08 must save pB0_0
//v09 must save pB0_1
//v10 must save pB0_2 --> ALPHA0
//v11 must save pB0_3
//v12 must save pB1_0
//v13 must save pB1_1
//v14 must save pB1_2
//v15 must save pB1_3
//v16 must save C00, C01
//v17 must save C02, C03
//v18 C04, C05
//v19 C06, C07
//v20 C10, C11
//v21 C12, C13
//v22 C14, C15
//v23 C16, C17
//v24 C20, C21
//v25 C22, C23
//v26 C24, C25
//v27 C26, C27
//v28 C30, C31
//v29 C32, C33
//v30 C34, C35
//v31 C36, C37

/*******************************************************************************
* Macro definitions
*******************************************************************************/

/* Input:
 *   pA: address of loading from *a
 *   pB: address of loading from *b
 *   p0.d: all true predicate register
 * Vector registers planning:
 *   z8 ~ z31  : c
 *   z0 ~ z2   : a
 *   z4 ~ z7   : b
 *   z3        : free, used in loop unrolling
 * Output:
 *   z8 ~ z31          : updated c
 *   pA: pointing to next elements of *a
 *   pB: pointing to next elements of *b
 */
.macro KERNEL12x8_I
	/* load a */
	ldr  z0, [pA]                        // load a(  0:3, l )
	ldr  z1, [pA, #1, MUL VL]            // load a(  4:7, l )
	ldr  z2, [pA, #2, MUL VL]            // load a(  8:11,l )
	/* load b */
	ld1rqd  {z4.d}, p0/z, [pB]           // load b( l,0:1 )

	/* Input:
	 *        z0, z1, z2    : a( 0:11,l )
	 *        z4, z5?, z6?, z7?: b( l,0:7  )
	 *        z3            : free
	 * Output:
	 *        z3, z1, z0      : a( 0:11,l+1 )
	 *        z4, z5, z6, z7? : b( l+1,0:7  )
	 *        z2              : free
	 * Note:
	 *   fmul implies z8 ~ z31 can be any initial value, no need to zero out.
	 *        That's different with using fmla.
	 */
	fmul  z8.d,  z0.d, z4.d[0]           //
	fmul  z9.d,  z1.d, z4.d[0]           //
	fmul  z10.d, z2.d, z4.d[0]           //
	ldr  z3, [pA, #3, MUL VL]            // load a(  0:3,l+1 )

	fmul  z11.d, z0.d, z4.d[1]           //
	fmul  z12.d, z1.d, z4.d[1]           //
	ld1rqd  {z5.d}, p0/z, [pB, #16]      // load b( l,2:3 )

	fmul  z13.d, z2.d, z4.d[1]           //
	fmul  z14.d, z0.d, z5.d[0]           //
	ld1rqd  {z4.d}, p0/z, [pB, #64]      // load b( l+1,0:1 )

	fmul  z17.d, z0.d, z5.d[1]           //
	fmul  z15.d, z1.d, z5.d[0]           //
	ld1rqd  {z6.d}, p0/z, [pB, #32]      // load b( l,4:5 )

	fmul  z18.d, z1.d, z5.d[1]           //
	fmul  z16.d, z2.d, z5.d[0]           //
	ld1rqd  {z7.d}, p0/z, [pB, #48]      // load b( l,6:7 )

	fmul  z19.d, z2.d, z5.d[1]           //
	fmul  z20.d, z0.d, z6.d[0]           //
	fmul  z23.d, z0.d, z6.d[1]           //
	ld1rqd  {z5.d}, p0/z, [pB, #80]      // load b( l+1,2:3 )

	fmul  z26.d, z0.d, z7.d[0]           //
	fmul  z29.d, z0.d, z7.d[1]           //
	fmul  z21.d, z1.d, z6.d[0]           //
	ldr  z0, [pA, #5, MUL VL]            // load a(  8:11,l+1 )

	fmul  z22.d, z2.d, z6.d[0]           //
	fmul  z24.d, z1.d, z6.d[1]           //
	fmul  z25.d, z2.d, z6.d[1]           //
	fmul  z27.d, z1.d, z7.d[0]           //
	ld1rqd  {z6.d}, p0/z, [pB, #96]      // load b( l+1,4:5 )

	fmul  z28.d, z2.d, z7.d[0]           //
	fmul  z30.d, z1.d, z7.d[1]           //
	fmul  z31.d, z2.d, z7.d[1]           //
	ldr  z1, [pA, #4, MUL VL]            // load a(  4:7, l )
.endm

.macro KERNEL12x8_M2
	/* Input:
	 *        z3, z1, z0      : a( 0:11,l )
	 *        z4, z5, z6, z7? : b( l,0:7  )
	 *        z2              : free
	 * Output:
	 *        z0, z1, z2      : a( 0:11,l )
	 *        z4, z5, z6, z7  : b( l,0:7  )
	 *        z3              : free
	 */
	fmla  z8.d,  z3.d, z4.d[0]           //
	fmla  z9.d,  z1.d, z4.d[0]           //
	fmla  z10.d, z0.d, z4.d[0]           //
	ldr  z2, [pA, #8, MUL VL]            // load a(  8:11,l+2 )

	/* adjust pB address */
	add pB, pB, #128                     //
	fmla  z11.d, z3.d, z4.d[1]           //
	fmla  z12.d, z1.d, z4.d[1]           //
	fmla  z13.d, z0.d, z4.d[1]           //
	fmla  z14.d, z3.d, z5.d[0]           //
	ld1rqd  {z4.d}, p0/z, [pB, #0]       // load b( l+2,0:1 )

	fmla  z15.d, z1.d, z5.d[0]           //
	fmla  z16.d, z0.d, z5.d[0]           //
	fmla  z17.d, z3.d, z5.d[1]           //
	ld1rqd  {z7.d}, p0/z, [pB, #-16]     // load b( l+1,6:7 )

	fmla  z18.d, z1.d, z5.d[1]           //
	fmla  z19.d, z0.d, z5.d[1]           //
	fmla  z20.d, z3.d, z6.d[0]           //
	ld1rqd  {z5.d}, p0/z, [pB, #16]      // load b( l+2,2:3 )

	fmla  z21.d, z1.d, z6.d[0]           //
	fmla  z22.d, z0.d, z6.d[0]           //
	fmla  z23.d, z3.d, z6.d[1]           //
	fmla  z24.d, z1.d, z6.d[1]           //
	fmla  z25.d, z0.d, z6.d[1]           //
	ld1rqd  {z6.d}, p0/z, [pB, #32]      // load b( l+2,4:5 )

	fmla  z28.d, z0.d, z7.d[0]           //
	fmla  z31.d, z0.d, z7.d[1]           //
	fmla  z26.d, z3.d, z7.d[0]           //
	ldr  z0, [pA, #6, MUL VL]            // load a(  0:3,l+2 )

	fmla  z30.d, z1.d, z7.d[1]           //
	fmla  z27.d, z1.d, z7.d[0]           //
	fmla  z29.d, z3.d, z7.d[1]           //
	ldr  z1, [pA, #7, MUL VL]            // load a(  4:7,l+2 )

	/* note: z7 be loaded in next iteration (M1) */
.endm

.macro KERNEL12x8_M1
	/* Input:
	 *        z0, z1, z2     : a( 0:11,l+2 )
	 *        z4, z5, z6, z7?: b( l+2,0:7  )
	 *        z3             : free
	 * Output:
	 *        z3, z1, z0      : a( 0:11,l+1 )
	 *        z4, z5, z6, z7? : b( l+1,0:7  )
	 *        z2              : free
	 */
	fmla  z8.d,  z0.d, z4.d[0]           //
	fmla  z9.d,  z1.d, z4.d[0]           //
	fmla  z10.d, z2.d, z4.d[0]           //
	/* load a for next iteration */
	ldr  z3, [pA, #9, MUL VL]            // load a(  0:3, l+3 )

	fmla  z11.d, z0.d, z4.d[1]           //
	fmla  z12.d, z1.d, z4.d[1]           //
	fmla  z13.d, z2.d, z4.d[1]           //
	fmla  z14.d, z0.d, z5.d[0]           //
	/* load b for next iteration */
	ld1rqd  {z4.d}, p0/z, [pB, #64]      // load b( l+3,0:1 )

	fmla  z15.d, z1.d, z5.d[0]           //
	fmla  z16.d, z2.d, z5.d[0]           //
	fmla  z17.d, z0.d, z5.d[1]           //
	ld1rqd  {z7.d}, p0/z, [pB, #48]      // load b( l+2,6:7 )

	fmla  z18.d, z1.d, z5.d[1]           //
	fmla  z19.d, z2.d, z5.d[1]           //
	fmla  z20.d, z0.d, z6.d[0]           //
	ld1rqd  {z5.d}, p0/z, [pB, #80]     // load b( l+3,2:3 )

	fmla  z23.d, z0.d, z6.d[1]           //
	fmla  z29.d, z0.d, z7.d[1]           //
	fmla  z26.d, z0.d, z7.d[0]           //
	fmla  z21.d, z1.d, z6.d[0]           //
	ldr  z0, [pA, #11, MUL VL]           // load a(  8:11,l+3 )

	fmla  z27.d, z1.d, z7.d[0]           //
	fmla  z30.d, z1.d, z7.d[1]           //
	fmla  z24.d, z1.d, z6.d[1]           //
	fmla  z22.d, z2.d, z6.d[0]           //
	ldr  z1, [pA, #10, MUL VL]           // load a(  4:7, l+3 )

	fmla  z25.d, z2.d, z6.d[1]           //
	fmla  z28.d, z2.d, z7.d[0]           //
	fmla  z31.d, z2.d, z7.d[1]           //
	ld1rqd  {z6.d}, p0/z, [pB, #96]      // load b( l+3,4:5 )

	/* adjust pA address */
	add pA, pA, #192                     // 6 * (MUL VL)

	/* note: z7 be loaded in next iteration (M2 or E) */
.endm

.macro KERNEL12x8_E
	/* Iteration E */
	/* Input:
	 *        z3, z1, z0      : a( 0:11,l )
	 *        z4, z5, z6, z7? : b( l,0:7  )
	 *        z2              : free
	 * Output:
	 *   z8 ~ z31         : updated c
	 *   z0 ~ z7          : free
	 *   pA: pointing to next elements of *a
	 *   pB: pointing to next elements of *b
	 */ 
	fmla  z8.d,  z3.d, z4.d[0]           //
	fmla  z9.d,  z1.d, z4.d[0]           //
	fmla  z10.d, z0.d, z4.d[0]           //
	ld1rqd  {z7.d}, p0/z, [pB, #112]     // load b( l+3,6:7 )

	fmla  z11.d, z3.d, z4.d[1]           //
	fmla  z12.d, z1.d, z4.d[1]           //
	fmla  z13.d, z0.d, z4.d[1]           //

	/* adjust pA address */
	add pA, pA, #192                     // 6 * (MUL VL)
	fmla  z14.d, z3.d, z5.d[0]           //
	fmla  z15.d, z1.d, z5.d[0]           //
	fmla  z16.d, z0.d, z5.d[0]           //

	fmla  z17.d, z3.d, z5.d[1]           //
	fmla  z18.d, z1.d, z5.d[1]           //
	fmla  z19.d, z0.d, z5.d[1]           //

	/* adjust pB address */
	add pB, pB, #128                     //
	fmla  z20.d, z3.d, z6.d[0]           //
	fmla  z21.d, z1.d, z6.d[0]           //
	fmla  z22.d, z0.d, z6.d[0]           //

	fmla  z23.d, z3.d, z6.d[1]           //
	fmla  z24.d, z1.d, z6.d[1]           //
	fmla  z25.d, z0.d, z6.d[1]           //

	fmla  z26.d, z3.d, z7.d[0]           //
	fmla  z27.d, z1.d, z7.d[0]           //
	fmla  z28.d, z0.d, z7.d[0]           //

	fmla  z29.d, z3.d, z7.d[1]           //
	fmla  z30.d, z1.d, z7.d[1]           //
	fmla  z31.d, z0.d, z7.d[1]           //
.endm

/*******************************************************************************
* End of macro definitions
*******************************************************************************/

	PROLOGUE

	.align 5
	add	sp, sp, #-(11 * 16)
	stp	d8, d9, [sp, #(0 * 16)]
	stp	d10, d11, [sp, #(1 * 16)]
	stp	d12, d13, [sp, #(2 * 16)]
	stp	d14, d15, [sp, #(3 * 16)]
	stp	d16, d17, [sp, #(4 * 16)]
	stp	x18, x19, [sp, #(5 * 16)]
	stp	x20, x21, [sp, #(6 * 16)]
	stp	x22, x23, [sp, #(7 * 16)]
	stp	x24, x25, [sp, #(8 * 16)]
	stp	x26, x27, [sp, #(9 * 16)]
	str	x28, [sp, #(10 * 16)]

	prfm	PLDL1KEEP, [origPB]
	prfm	PLDL1KEEP, [origPA]

	/* transfer alpha to scalar register 'alpha' */
	fmov	alpha, d0                      // alpha
	ptrue	p0.d, all

	lsl	LDC, LDC, #3			// ldc = ldc * 8

#if !defined(LEFT)
	neg	tempOffset, offset
#endif
	mov	pB, origPB

	mov	counterJ, origN
	asr 	counterJ, counterJ, #3		// J = J / 8
	cmp 	counterJ, #0
	ble	.Ldtrmm_kernel_L4_BEGIN

/******************************************************************************/

.Ldtrmm_kernel_L8_BEGIN:
	mov	pCRow0, pC
	add	pCRow1, pCRow0, LDC
	add	pCRow2, pCRow1, LDC
	add	pCRow3, pCRow2, LDC
	add	pCRow4, pCRow3, LDC
	add	pCRow5, pCRow4, LDC
	add	pCRow6, pCRow5, LDC
	add	pCRow7, pCRow6, LDC

	add	pC, pCRow7, LDC

#if defined(LEFT)
	mov	tempOffset, offset
#endif
	mov	pA, origPA			// pA = start of A array

.Ldtrmm_kernel_L8_M12_BEGIN:

	mov	counterI, origM
	subs	counterI, counterI, #12		// M12
	blt	.Ldtrmm_kernel_L8_M8_BEGIN	// if counterI < 12

	.align 5
.Ldtrmm_kernel_L8_M12_20:

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
						// 12 * sizeof(double) = (8+4)*8 = 2^6 + 2^5
						//  , where 12 is M12
	add	pA, pA, tempOffset, lsl #6
	add	pA, pA, tempOffset, lsl #5
						// #6 = 2^6 = 8 * sizeof( double ), where 8 is L8
	add	pB, pB, tempOffset, lsl #6	// pB moves in its k dimension, L4
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset        // tempK is the number of non-zero rows/cols.
#elif defined(LEFT)                             //  defined(LEFT) &&  defined(TRANSA)
	add	tempK, tempOffset, #12          // M12
#else                                           // !defined(LEFT) && !defined(TRANSA)
	add	tempK, tempOffset, #8           // L8
#endif

	asr 	counterL , tempK, #3		// L = K / 8, unroll at 8
	cmp	counterL , #2			// is there at least 4 to do?
	blt	.Ldtrmm_kernel_L8_M12_32

	KERNEL12x8_I				// do one in the K
	KERNEL12x8_M2				// do another in the K
	KERNEL12x8_M1
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_M2

	subs	counterL, counterL, #2		// subtract 2
	ble	.Ldtrmm_kernel_L8_M12_22a

	.align 5
.Ldtrmm_kernel_L8_M12_22:

	KERNEL12x8_M1
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_M2

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L8_M12_22

	.align 5
.Ldtrmm_kernel_L8_M12_22a:

	KERNEL12x8_M1
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_E

	b	 .Ldtrmm_kernel_L8_M12_44

	.align 5
.Ldtrmm_kernel_L8_M12_32:

	tst	counterL, #1
	ble	.Ldtrmm_kernel_L8_M12_40

	KERNEL12x8_I
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_E

	b	.Ldtrmm_kernel_L8_M12_44

.Ldtrmm_kernel_L8_M12_40:

	INIT12x8

.Ldtrmm_kernel_L8_M12_44:

	ands	counterL , tempK, #7
	ble	.Ldtrmm_kernel_L8_M12_100

	.align 5
.Ldtrmm_kernel_L4_M16_46:

	KERNEL12x8_SUB

	subs	counterL, counterL, #1
	bne	.Ldtrmm_kernel_L8_M12_46

.Ldtrmm_kernel_L8_M12_100:

	SAVE12x8

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)                               //  defined(LEFT) && defined(TRANSA)
	sub	tempK, tempK, #12               // M12
#else                                           // !defined(LEFT) && !defined(TRANSA)
	sub	tempK, tempK, #8                // L8
#endif
						// 12 * sizeof(double) = (8+4)*8 = 2^6 + 2^5
						//  , where 12 is M12
	add	pA, pA, tempK, lsl #6
	add	pA, pA, tempK, lsl #5
	add	pB, pB, tempK, lsl #6           // #6 = 2^6 = 8 * sizeof(double), where 8 is L8
#endif

#if defined(LEFT)
	add	tempOffset, tempOffset, #12     // M12
#endif
/* TODO: adjust prfm
 * 	prfm	PLDL1KEEP, [pA]
 * 	prfm	PLDL1KEEP, [pA, #64]
 * 	prfm	PLDL1KEEP, [origPB]
 */
.Ldtrmm_kernel_L8_M12_END:
	subs	counterI, counterI, #12		// 12 rows in each cycle
	bge	.Ldtrmm_kernel_L8_M12_20	// greater than or equal

.Ldtrmm_kernel_L8_M8_BEGIN:
	adds	counterI, counterI, #12		// counterI = origM % 12
	beq	.Ldtrmm_kernel_L8_END           // if zero, no need to run L8_M8/M4/M2/M1

	tst	counterI, #8                    // 8 is from M8, to test wether L8_M8 need to run
	ble	.Ldtrmm_kernel_L8_M4_BEGIN      // if zero, skip L8_M8, go test next bit (#4) in L8_M4

.Ldtrmm_kernel_L8_M8_20:

	INIT8x8

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #6            // #6 = 2^6 = 8 * sizeof( double ), where 8 is L8/M8
	add	pB, pB, temp                    // pB is moved 'tempOffset' pos. in its k dimension L8
	add	pA, pA, temp                    // pA is moved 'tempOffset' pos. in its k dimension M8
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #8           // 8 is M8
#else
	add	tempK, tempOffset, #8           // 8 is L8
#endif

	asr 	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Ldtrmm_kernel_L8_M8_40

.Ldtrmm_kernel_L8_M8_22:

/* TODO: it worths to write below kernels,
 *          unroll them, and interleave fops & memops
 */
	KERNEL8x8_SUB
	KERNEL8x8_SUB
	KERNEL8x8_SUB
	KERNEL8x8_SUB

	KERNEL8x8_SUB
	KERNEL8x8_SUB
	KERNEL8x8_SUB
	KERNEL8x8_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L8_M8_22

.Ldtrmm_kernel_L8_M8_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L8_M8_100

.Ldtrmm_kernel_L8_M8_42:

	KERNEL8x8_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L8_M8_42

.Ldtrmm_kernel_L8_M8_100:

	SAVE8x8

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #8                // 8 is M8
#else
	sub	tempK, tempK, #8                // 8 is L8
#endif
	lsl	temp, tempK, #6                 // #6 = 2^6 = 8 * sizeof(double), where 8 is M8/L8
	add	pA, pA, temp
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #8      // 8 is M8
#endif

.Ldtrmm_kernel_L8_M8_END:


.Ldtrmm_kernel_L8_M4_BEGIN:
						// counterI = origM % 12
	tst	counterI , #7                   // 7 is (M8-1)
	ble	.Ldtrmm_kernel_L8_END           // if zero, no need to run L8_M4/M2/M1

	tst	counterI, #4	        	// 4 is M4, to test whether L8_M4 need to run
	ble	.Ldtrmm_kernel_L8_M2_BEGIN      // if zero, skip L8_M4, go test next bit (#2) in L8_M2

.Ldtrmm_kernel_L8_M4_20:

	INIT4x8

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
						// #6 = 2^6 = 8 * sizeof( double ), where 8 is L8
	add	pB, pB, tempOffset, lsl #6      // pB is moved 'tempOffset' pos. in its k dimension L8
						// #5 = 2^5 = 4 * sizeof( double ), where 4 is M4
	add	pA, pA, tempOffset, lsl #5	// pA is moved 'tempOffset' pos. in its k dimension M4
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #4           // 4 is M4
#else
	add	tempK, tempOffset, #8           // 8 is L8
#endif
	asr 	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Ldtrmm_kernel_L8_M4_40

.Ldtrmm_kernel_L8_M4_22:

/* TODO: to be unrolled, for better fops/memops interleaving */
	KERNEL4x8_SUB
	KERNEL4x8_SUB
	KERNEL4x8_SUB
	KERNEL4x8_SUB

	KERNEL4x8_SUB
	KERNEL4x8_SUB
	KERNEL4x8_SUB
	KERNEL4x8_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L8_M4_22


.Ldtrmm_kernel_L8_M4_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L8_M4_100

.Ldtrmm_kernel_L8_M4_42:

	KERNEL4x8_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L8_M4_42

.Ldtrmm_kernel_L8_M4_100:

	SAVE4x8

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #4                // 4 is M4
#else
	sub	tempK, tempK, #8                // 8 is L8
#endif
	add	pA, pA, tempK, lsl #5           // #5 = 2^5 = 4 * sizeof(double), where 4 is M4
	add	pB, pB, tempK, lsl #6           // #6 = 2^6 = 8 * sizeof(double), where 8 is L8
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #4      // 4 is M4
#endif

.Ldtrmm_kernel_L8_M4_END:

.Ldtrmm_kernel_L8_M2_BEGIN:
						// counterI = origM % 12
	tst	counterI , #3                   // 3 is (M4-1)
	ble	.Ldtrmm_kernel_L8_END           // if zero, no need to run L8_M2/M1

	tst	counterI, #2			// 2 is M2
	ble	.Ldtrmm_kernel_L8_M1_BEGIN      // if zero, skip L8_M2, go test next bit (#1) for L8_M1

.Ldtrmm_kernel_L8_M2_20:

	INIT2x8

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
						// #4 = 2^4 = 2 * sizeof( double ), where 2 is M2
	add	pA, pA, tempOffset, lsl #4      // pA is moved M2 in its k dimension
						// #6 = 2^6 = 8 * sizeof( double ), where 8 is L8
	add	pB, pB, tempOffset, lsl #6      // pB is moved L8 in its k dimension
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #2           // 2 is M2
#else
	add	tempK, tempOffset, #8           // 8 is L8
#endif
	asr 	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Ldtrmm_kernel_L8_M2_40

.Ldtrmm_kernel_L8_M2_22:

/* TODO: to be unrolled, for better fops/memops interleaving */
	KERNEL2x8_SUB
	KERNEL2x8_SUB
	KERNEL2x8_SUB
	KERNEL2x8_SUB

	KERNEL2x8_SUB
	KERNEL2x8_SUB
	KERNEL2x8_SUB
	KERNEL2x8_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L8_M2_22

.Ldtrmm_kernel_L8_M2_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L8_M2_100

.Ldtrmm_kernel_L8_M2_42:

	KERNEL2x8_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L8_M2_42

.Ldtrmm_kernel_L8_M2_100:

	SAVE2x8

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #2                // 2 is M2
#else
	sub	tempK, tempK, #8                // 8 is L8
#endif
	add	pA, pA, tempK, lsl #4		// #4 = 2^4 = 2 * sizeof(double), where 2 is M2
	add	pB, pB, tempK, lsl #6		// #6 = 2^6 = 8 * sizeof(double), where 8 is L8
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #2      // 2 is M2
#endif

.Ldtrmm_kernel_L8_M2_END:


.Ldtrmm_kernel_L8_M1_BEGIN:

	tst	counterI, #1			// counterI = counterI % 2
	ble	.Ldtrmm_kernel_L8_END           // if zero, no need to run L8_M1

.Ldtrmm_kernel_L8_M1_20:

	INIT1x8

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
						// #3 = 2^3 = 1 * sizeof( double ), where 1 is M1
	add	pA, pA, tempOffset, lsl #3      // pA is moved M1 in its k dimension
						// #6 = 2^6 = 8 * sizeof( double ), where 8 is L8
	add	pB, pB, tempOffset, lsl #6      // pB is moved L8 in its k dimension
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #1           // 1 is M1
#else
	add	tempK, tempOffset, #8           // 8 is L8
#endif

	asr 	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Ldtrmm_kernel_L8_M1_40

.Ldtrmm_kernel_L8_M1_22:

/* TODO: to be unrolled, for better fops/memops interleaving */
	KERNEL1x8_SUB
	KERNEL1x8_SUB
	KERNEL1x8_SUB
	KERNEL1x8_SUB

	KERNEL1x8_SUB
	KERNEL1x8_SUB
	KERNEL1x8_SUB
	KERNEL1x8_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L8_M1_22


.Ldtrmm_kernel_L8_M1_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L8_M1_100

.Ldtrmm_kernel_L8_M1_42:

	KERNEL1x8_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L8_M1_42

.Ldtrmm_kernel_L8_M1_100:

	SAVE1x8

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #1                // 1 is M1
#else
	sub	tempK, tempK, #8                // 8 is L8
#endif
	add	pA, pA, tempK, lsl #3           // #3 = 2^3 = 1 * sizeof(double), where 1 is M1
	add	pB, pB, tempK, lsl #6           // #6 = 2^6 = 8 * sizeof(double), where 8 is L8
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #1      // 1 is M1
#endif

.Ldtrmm_kernel_L8_END:
						// #6 =2^6 = 8 * sizeof(double), where 8 is L8
	add	origPB, origPB, origK, lsl #6	// B = B + K * 8 * 8

#if !defined(LEFT)
	add	tempOffset, tempOffset, #8      // 8 is L8
#endif

	subs	counterJ, counterJ , #1		// j--
	bgt	.Ldtrmm_kernel_L8_BEGIN

/******************************************************************************/

.Ldtrmm_kernel_L4_BEGIN:   // less than 8 left in N direction

	mov	counterJ, origN
	tst	counterJ, #7                    // 7 is (L8-1), wether L4/L2/L1 is needed.
	ble	.Ldtrmm_kernel_L999             // all zero, done

	tst	counterJ , #4                   // Is L4 needed?
	ble	.Ldtrmm_kernel_L2_BEGIN         // bit in "#4" is zero, jump to L2

	mov	pCRow0, pC			// pCRow0 = pC
                                                // pCRow1 will be calculated in SAVE12x4
	add	pC, pC, LDC, lsl #2             // pC += 4 * LDC, where 4 is L4
                                                // pC moves to the next slice in N direction

#if defined(LEFT)
	mov	tempOffset, offset
#endif
	mov	pA, origPA			// pA = A

.Ldtrmm_kernel_L4_M12_BEGIN:

	mov	counterI, origM
	subs	counterI, counterI, #12		// M12
	blt	.Ldtrmm_kernel_L4_M8_BEGIN	// if counterI < 12

.Ldtrmm_kernel_L4_M12_20:

	INIT12x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	add	pA, pA, tempOffset, lsl #6	// 12 * sizeof(double) = (8+4)*8 = 2^6 + 2^5, where 12 is M12
	add	pA, pA, tempOffset, lsl #5	// pA moves 'tempOffset' pos. in its k dimension M12
	add	pB, pB, tempOffset, lsl #5	// #5 = 2^5 = 4 * sizeof( double ), where 4 is L4
						// pB moves 'tempOffset' pos. in its k dimension, L4
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #12		// 12 is M12
#else
	add	tempK, tempOffset, #4		// 4 is L4
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8, unroll at 8
	cmp	counterL,#0
	ble	.Ldtrmm_kernel_L4_M12_40
	.align 5

.Ldtrmm_kernel_L4_M12_22:

/* TODO: optimization by interleaving fops/memops instructions */
	KERNEL12x4_SUB
	KERNEL12x4_SUB
	KERNEL12x4_SUB
	KERNEL12x4_SUB

	KERNEL12x4_SUB
	KERNEL12x4_SUB
	KERNEL12x4_SUB
	KERNEL12x4_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L4_M12_22

.Ldtrmm_kernel_L4_M12_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L4_M12_100       // remainder is zero, finish

.Ldtrmm_kernel_L4_M12_42:

	KERNEL12x4_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L4_M12_42

.Ldtrmm_kernel_L4_M12_100:

	SAVE12x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #12               // 12 is M12
#else
	sub	tempK, tempK, #4                // 4 is L4
#endif
	add	pA, pA, tempK, lsl #6		// 12 * sizeof(double) = (8+4)*8 = 2^6 + 2^5, where 12 is M12
	add	pA, pA, tempK, lsl #5		// pA moves 'tempOffset' pos. in its k dimension M12
	add	pB, pB, tempK, lsl #5           // #5 = 2^5 = 4 * sizeof(double), where 4 is L4
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #12     // M12
#endif

.Ldtrmm_kernel_L4_M12_END:
	subs	counterI, counterI, #12		// M12, 12 rows in each cycle
	bge	.Ldtrmm_kernel_L4_M12_20	// greater than or equal


.Ldtrmm_kernel_L4_M8_BEGIN:
	adds	counterI, counterI, #12		// counterI = origM % 12
	beq	.Ldtrmm_kernel_L4_END           // if zero, skip L4_M8/M4/M2/M1

	tst	counterI, #8			// 8 is M8, test bit
	ble	.Ldtrmm_kernel_L4_M4_BEGIN      // if zero, skip L4_M8, go test L4_M4

.Ldtrmm_kernel_L4_M8_20:

	INIT8x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	add	pB, pB, tempOffset, lsl #5      // #5 = 2^5 = 4 * sizeof( double ), where 4 is L4
						// pB is moved 'tempOffset' pos. in its k dimension L4
	add	pA, pA, tempOffset, lsl #6      // #6 = 2^6 = 8 * sizeof( double ), where 8 is M8
						// pA is moved 'tempOffset' pos. in its k dimension M4
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #8           // 8 is M8
#else
	add	tempK, tempOffset, #4           // 4 is L4
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL,#0
	ble	.Ldtrmm_kernel_L4_M8_40
	.align 5

.Ldtrmm_kernel_L4_M8_22:
	KERNEL8x4_SUB
	KERNEL8x4_SUB
	KERNEL8x4_SUB
	KERNEL8x4_SUB

	KERNEL8x4_SUB
	KERNEL8x4_SUB
	KERNEL8x4_SUB
	KERNEL8x4_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L4_M8_22

.Ldtrmm_kernel_L4_M8_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L4_M8_100

.Ldtrmm_kernel_L4_M8_42:

	KERNEL8x4_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L4_M8_42

.Ldtrmm_kernel_L4_M8_100:

	SAVE8x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #8                // 8 is M8
#else
	sub	tempK, tempK, #4                // L4
#endif
				
	add	pA, pA, tempK, lsl #6		// #6 = 2^6 = 8 * sizeof(double), where 8 is M8
	add	pB, pB, tempK, lsl #5           // #5 = 2^5 = 4 * sizeof(double), where 4 is L4
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #8      // 8 is M8
#endif

.Ldtrmm_kernel_L4_M8_END:


.Ldtrmm_kernel_L4_M4_BEGIN:
						// counterI = origM % 12
	tst	counterI , #7                   // 7 is (M8-1)
	ble	.Ldtrmm_kernel_L4_END           // if zero, skip L4_M4/M2/M1

	tst	counterI, #4			// 4 is M4
	                                        // test bit (#4)
	ble	.Ldtrmm_kernel_L4_M2_BEGIN      // If zero, skip M4; jump to M2

.Ldtrmm_kernel_L4_M4_20:

	INIT4x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	add	pA, pA, tempOffset, lsl #5           // #5 = 2^5 = 4 * sizeof(double), where 4 is M4/L4
	add	pB, pB, tempOffset, lsl #5           // #5 = 2^5 = 4 * sizeof(double), where 4 is M4/L4
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #4           // 4 is M4
#else
	add	tempK, tempOffset, #4           // 4 is L4
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL,#0
	ble	.Ldtrmm_kernel_L4_M4_40
	.align 5

.Ldtrmm_kernel_L4_M4_22:
/* TODO: to be unrolled for better fops/memops interleaving */
	KERNEL4x4_SUB
	KERNEL4x4_SUB
	KERNEL4x4_SUB
	KERNEL4x4_SUB

	KERNEL4x4_SUB
	KERNEL4x4_SUB
	KERNEL4x4_SUB
	KERNEL4x4_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L4_M4_22


.Ldtrmm_kernel_L4_M4_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L4_M4_100

.Ldtrmm_kernel_L4_M4_42:

	KERNEL4x4_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L4_M4_42

.Ldtrmm_kernel_L4_M4_100:

	SAVE4x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #4                // 4 is M4
#else
	sub	tempK, tempK, #4                // 4 is L4
#endif
	add	pA, pA, tempK, lsl #5           // #5 = 2^5 = 4 * sizeof(double), where 4 is M4/L4
	add	pB, pB, tempK, lsl #5           // #5 = 2^5 = 4 * sizeof(double), where 4 is M4/L4
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #4      // 4 is M4
#endif

.Ldtrmm_kernel_L4_M4_END:


.Ldtrmm_kernel_L4_M2_BEGIN:
						// counterI = origM % 12
	tst	counterI , #3                   // 3 is (M4-1)
	ble	.Ldtrmm_kernel_L4_END           // if zero, skip L4_M2/M1

	tst	counterI, #2			// 2 is M2; counterI = counterI / 2
	ble	.Ldtrmm_kernel_L4_M1_BEGIN      // test bit. If zero, skip L2_M2 and jump to L2_M1

.Ldtrmm_kernel_L4_M2_20:

	INIT2x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	add	pB, pB, tempOffset, lsl #5	// #5 = 2^5 = 4 * sizeof(double), where 4 is L4
	add	pA, pA, tempOffset, lsl #4	// #4 = 2^4 = 2 * sizeof(double), where 2 is M2
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #2           // 2 is M2
#else
	add	tempK, tempOffset, #4           // 2 is L4
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8
        cmp	counterL,#0
	ble	.Ldtrmm_kernel_L4_M2_40

.Ldtrmm_kernel_L4_M2_22:
/* TODO: to be unrolled, for better fops/memops interleaving */
	KERNEL2x4_SUB
	KERNEL2x4_SUB
	KERNEL2x4_SUB
	KERNEL2x4_SUB

	KERNEL2x4_SUB
	KERNEL2x4_SUB
	KERNEL2x4_SUB
	KERNEL2x4_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L4_M2_22

.Ldtrmm_kernel_L4_M2_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L4_M2_100

.Ldtrmm_kernel_L4_M2_42:

	KERNEL2x4_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L4_M2_42

.Ldtrmm_kernel_L4_M2_100:

	SAVE2x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #2                // 2 is M2
#else
	sub	tempK, tempK, #4                // 4 is L4
#endif
	add	pB, pB, tempK, lsl #5		// #5 = 2^5 = 4 * sizeof(double), where 4 is L4
	add	pA, pA, tempK, lsl #4		// #4 = 2^4 = 2 * sizeof(double), where 2 is M2
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #2      // 2 is M2
#endif

.Ldtrmm_kernel_L4_M2_END:


.Ldtrmm_kernel_L4_M1_BEGIN:
						// counterI = origM % 12
	tst	counterI, #1
	ble	.Ldtrmm_kernel_L2_END           // if zero, done L4.

.Ldtrmm_kernel_L4_M1_20:

	INIT1x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, temp
	add	pA, pA, tempOffset, lsl #3	// #3 = 2^3 = 1 * sizeof( double ), where 1 is M1
						// pA is moved M2 in its k dimension
	add	pB, pB, tempOffset, lsl #5	// #5 = 2^5 = 4 * sizeof( double ), where 4 is L4
						// pB is moved L4 in its k dimension
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #1           // 1 is M1
#else
	add	tempK, tempOffset, #4           // L4
#endif

	asr 	counterL , tempK, #3		// counterL = counterL / 8
        cmp     counterL, #0
	ble	.Ldtrmm_kernel_L4_M1_40

.Ldtrmm_kernel_L4_M1_22:
/* TODO: to be unrolled, for better fops/memops interleaving */
	KERNEL1x4_SUB
	KERNEL1x4_SUB
	KERNEL1x4_SUB
	KERNEL1x4_SUB

	KERNEL1x4_SUB
	KERNEL1x4_SUB
	KERNEL1x4_SUB
	KERNEL1x4_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L2_M1_22

.Ldtrmm_kernel_L4_M1_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L4_M1_100

.Ldtrmm_kernel_L4_M1_42:

	KERNEL1x4_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L4_M1_42

.Ldtrmm_kernel_L4_M1_100:

	SAVE1x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #1                // 1 is M1
#else
	sub	tempK, tempK, #4                // 4 is L4
#endif
	add	pA, pA, tempK, lsl #3           // #3 = 2^3 = 1 * sizeof(double), where 1 is M1
	add	pB, pB, tempK, lsl #5           // #5 = 2^5 = 4 * sizeof(double), where 4 is L4
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #1      // #1 is M1
#endif

.Ldtrmm_kernel_L4_END:
#if !defined(LEFT)
	add	tempOffset, tempOffset, #4      // 4 is L4
#endif
	add	origPB, origPB, origK, lsl #5	// B = B + K * 4 * 8, 4 is L4, 8 is sizeof(double)

/******************************************************************************/

.Ldtrmm_kernel_L2_BEGIN:   // less than 4 left in N direction

	mov	counterJ , origN		// TODO: this instruction is redundant.
						//    counterJ = origN since L4_BEGIN
	tst	counterJ , #3                   // 3 is (L4-1), wether L2/L1 is needed.
	ble	.Ldtrmm_kernel_L999             // all zero, done

	tst	counterJ , #2                   // 2 is L2. Is L2 needed?
	ble	.Ldtrmm_kernel_L1_BEGIN         // bit "#2" is zero, jump to L1

	mov	pCRow0, pC			// pCRow0 = pC
                                                // pCRow1 will be calculated in SAVE16x2
	add	pC,pC,LDC, lsl #1               // pC += 2 * LDC, where 2 is L2
                                                // pC moves to the next slice in N direction

#if defined(LEFT)
	mov	tempOffset, offset
#endif
	mov	pA, origPA			// pA = A

.Ldtrmm_kernel_L2_M12_BEGIN:

	mov	counterI, origM
	subs	counterI, counterI, #12		// M12
	blt	.Ldtrmm_kernel_L2_M8_BEGIN	// if counterI < 12
	
.Ldtrmm_kernel_L2_M12_20:

	INIT12x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	add	pA, pA, tempOffset, lsl #6	// 12 * sizeof(double) = (8+4)*8 = 2^6 + 2^5, where 12 is M12
	add	pA, pA, tempOffset, lsl #5	// pA moves 'tempOffset' pos. in its k dimension M12
	add	pB, pB, tempOffset, lsl #4	// #4 = 2^4 = 2 * sizeof( double ), where 2 is L2
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #12  // M12
#else
	add	tempK, tempOffset, #2   // 2 is L2
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8, unroll at 8
	cmp	counterL,#0
	ble	.Ldtrmm_kernel_L2_M12_40
	.align 5

.Ldtrmm_kernel_L2_M12_22:

/* TODO: optimization by interleaving fops/memops instructions */
	KERNEL12x2_SUB
	KERNEL12x2_SUB
	KERNEL12x2_SUB
	KERNEL12x2_SUB

	KERNEL12x2_SUB
	KERNEL12x2_SUB
	KERNEL12x2_SUB
	KERNEL12x2_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L2_M12_22


.Ldtrmm_kernel_L2_M12_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L2_M12_100       // remainder is zero, finish

.Ldtrmm_kernel_L2_M12_42:

	KERNEL12x2_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L2_M12_42

.Ldtrmm_kernel_L2_M12_100:

	SAVE12x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #12               // M12
#else
	sub	tempK, tempK, #2                // L2
#endif
	add	pA, pA, tempK, lsl #6		// 12 * sizeof(double) = (8+4)*8 = 2^6 + 2^5, where 12 is M12
	add	pA, pA, tempK, lsl #5		// pA moves 'tempOffset' pos. in its k dimension M12
	add	pB, pB, tempK, lsl #4           // #4 = 2^4 = 2 * sizeof(double), where 2 is L2
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #12     // M12
#endif

.Ldtrmm_kernel_L2_M12_END:
	subs	counterI, counterI, #12		// M12, 12 rows in each cycle
	bge	.Ldtrmm_kernel_L2_M12_20	// greater than or equal

.Ldtrmm_kernel_L2_M8_BEGIN:

	adds	counterI, counterI, #12		// counterI = origM % 12
	beq	.Ldtrmm_kernel_L2_END           // if zero, skip L2_M8/M4/M2/M1

	tst	counterI, #8			// 8 is M8, test bit
	ble	.Ldtrmm_kernel_L2_M4_BEGIN      // if zero, skip L2_M8, go test L2_M4

.Ldtrmm_kernel_L2_M8_20:

	INIT8x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	/* todo: lsl and add can be combined into one 'add' instruction */
	lsl	temp, tempOffset, #4            // #4 = 2^4 = 2 * sizeof( double ), where 2 is L2
	add	pB, pB, temp                    // pB is moved 'tempOffset' pos. in its k dimension L4
	lsl	temp, tempOffset, #6            // #6 = 2^6 = 8 * sizeof( double ), where 8 is M8
	add	pA, pA, temp                    // pA is moved 'tempOffset' pos. in its k dimension M4
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #8           // 8 is M8
#else
	add	tempK, tempOffset, #2           // 2 is L2
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL,#0
	ble	.Ldtrmm_kernel_L2_M8_40
	.align 5

.Ldtrmm_kernel_L2_M8_22:
	KERNEL8x2_SUB
	KERNEL8x2_SUB
	KERNEL8x2_SUB
	KERNEL8x2_SUB

	KERNEL8x2_SUB
	KERNEL8x2_SUB
	KERNEL8x2_SUB
	KERNEL8x2_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L2_M8_22


.Ldtrmm_kernel_L2_M8_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L2_M8_100

.Ldtrmm_kernel_L2_M8_42:

	KERNEL8x2_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L2_M8_42

.Ldtrmm_kernel_L2_M8_100:

	SAVE8x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #8                // 8 is M8
#else
	sub	tempK, tempK, #2                // 2 is L2
#endif
	/* todo: combined lsl and add into one */
	lsl	temp, tempK, #6                 // #6 = 2^6 = 8 * sizeof(double), where 8 is M8
	add	pA, pA, temp
	lsl	temp, tempK, #4                 // #4 = 2^4 = 2 * sizeof(double), where 2 is L2
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #8      // 8 is M8
#endif

.Ldtrmm_kernel_L2_M8_END:


.Ldtrmm_kernel_L2_M4_BEGIN:

						// counterI = origM % 12
	tst	counterI , #7                   // 7 is (M8-1)
	ble	.Ldtrmm_kernel_L2_END           // if zero, skip L2_M4/M2/M1

	tst	counterI, #4			// 4 is M4
	                                        // test bit (#4)
	ble	.Ldtrmm_kernel_L2_M2_BEGIN      // If zero, skip M4; jump to M2

.Ldtrmm_kernel_L2_M4_20:

	INIT4x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	/* todo: combined lsl and add into one */
	lsl	temp, tempOffset, #4            // #4=2^4=2*sizeof(double), 2 is L2
	add	pB, pB, temp
	lsl	temp, tempOffset, #5            // #5=2^5=4*sizeof(double), 4 is M4
	add	pA, pA, temp
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #4           // 4 is M4
#else
	add	tempK, tempOffset, #2           // 2 is L2
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL,#0
	ble	.Ldtrmm_kernel_L2_M4_40
	.align 5

.Ldtrmm_kernel_L2_M4_22:
/* TODO: to be unrolled for better fops/memops interleaving */
	KERNEL4x2_SUB
	KERNEL4x2_SUB
	KERNEL4x2_SUB
	KERNEL4x2_SUB

	KERNEL4x2_SUB
	KERNEL4x2_SUB
	KERNEL4x2_SUB
	KERNEL4x2_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L2_M4_22


.Ldtrmm_kernel_L2_M4_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L2_M4_100

.Ldtrmm_kernel_L2_M4_42:

	KERNEL4x2_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L2_M4_42

.Ldtrmm_kernel_L2_M4_100:

	SAVE4x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #4                // 4 is M4
#else
	sub	tempK, tempK, #2                // 2 is L2
#endif
	/* todo: combined lsl and add into one */
	lsl	temp, tempK, #5                 // #5 = 2^5 = 4 * sizeof(double), where 4 is M4
	add	pA, pA, temp
	lsl	temp, tempK, #4                 // #4 = 2^4 = 2 * sizeof(double), where 2 is L2
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #4      // 4 is M4
#endif

.Ldtrmm_kernel_L2_M4_END:


.Ldtrmm_kernel_L2_M2_BEGIN:
						// counterI = origM % 12
	tst	counterI , #3                   // 3 is (M4-1)
	ble	.Ldtrmm_kernel_L2_END           // if zero, skip L2_M2/M1

	tst	counterI, #2			// 2 is M2; counterI = counterI / 2
	ble	.Ldtrmm_kernel_L2_M1_BEGIN      // test bit. If zero, skip L2_M2 and jump to L2_M1

.Ldtrmm_kernel_L2_M2_20:

	INIT2x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	/* todo: combined lsl and add into one */
	lsl	temp, tempOffset, #4            // #4 = 2^4 = 2 * sizeof(double), where 2 is M2, and L2
	add	pB, pB, temp
	add	pA, pA, temp
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #2           // 2 is M2
#else
	add	tempK, tempOffset, #2           // 2 is L2
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8
        cmp	counterL,#0
	ble	.Ldtrmm_kernel_L2_M2_40

.Ldtrmm_kernel_L2_M2_22:
/* TODO: to be unrolled, for better fops/memops interleaving */
	KERNEL2x2_SUB
	KERNEL2x2_SUB
	KERNEL2x2_SUB
	KERNEL2x2_SUB

	KERNEL2x2_SUB
	KERNEL2x2_SUB
	KERNEL2x2_SUB
	KERNEL2x2_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L2_M2_22

.Ldtrmm_kernel_L2_M2_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L2_M2_100

.Ldtrmm_kernel_L2_M2_42:

	KERNEL2x2_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L2_M2_42

.Ldtrmm_kernel_L2_M2_100:

	SAVE2x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #2                // 2 is M2
#else
	sub	tempK, tempK, #2                // 2 is L2
#endif
	/* todo: combined lsl and add into one */
	lsl	temp, tempK, #4                 // #4 = 2^4 = 2 * sizeof(double), where 2 is M2/L2
	add	pA, pA, temp
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #2      // 2 is M2
#endif

.Ldtrmm_kernel_L2_M2_END:


.Ldtrmm_kernel_L2_M1_BEGIN:
						// counterI = origM % 12
	tst	counterI, #1			// counterI = counterI % 2
	ble	.Ldtrmm_kernel_L2_END           // if zero, done L2.

.Ldtrmm_kernel_L2_M1_20:

	INIT1x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	/* todo: combined lsl and add into one */
	lsl	temp, tempOffset, #3            // #3 = 2^3 = 1 * sizeof( double ), where 1 is M1
	add	pA, pA, temp                    // pA is moved M2 in its k dimension
	lsl	temp, tempOffset, #4            // #4 = 2^4 = 2 * sizeof( double ), where 2 is L2
	add	pB, pB, temp                    // pB is moved L4 in its k dimension
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #1           // 1 is M1
#else
	add	tempK, tempOffset, #2           // 2 is L2
#endif

	asr 	counterL , tempK, #3		// counterL = counterL / 8
        cmp     counterL, #0
	ble	.Ldtrmm_kernel_L2_M1_40

.Ldtrmm_kernel_L2_M1_22:
/* TODO: to be unrolled, for better fops/memops interleaving */
	KERNEL1x2_SUB
	KERNEL1x2_SUB
	KERNEL1x2_SUB
	KERNEL1x2_SUB

	KERNEL1x2_SUB
	KERNEL1x2_SUB
	KERNEL1x2_SUB
	KERNEL1x2_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L2_M1_22


.Ldtrmm_kernel_L2_M1_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L2_M1_100

.Ldtrmm_kernel_L2_M1_42:

	KERNEL1x2_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L2_M1_42

.Ldtrmm_kernel_L2_M1_100:

	SAVE1x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #1                // 1 is M1
#else
	sub	tempK, tempK, #2                // 2 is L2
#endif
	/* todo: combined lsl and add into one */
	lsl	temp, tempK, #3                 // #3 = 2^3 = 1 * sizeof(double), where 1 is M1
	add	pA, pA, temp
	lsl	temp, tempK, #4                 // #4 = 2^4 = 2 * sizeof(double), where 2 is L2
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #1      // #1 is M1
#endif

.Ldtrmm_kernel_L2_END:
#if !defined(LEFT)
	add	tempOffset, tempOffset, #2      // 2 is L2
#endif
	add	origPB, origPB, origK, lsl #4	// B = B + K * 2 * 8, 2 is L2, 8 is sizeof(double)

/******************************************************************************/

.Ldtrmm_kernel_L1_BEGIN:

	mov	counterJ , origN		// TODO: this instruction is redundant.
						//    counterJ = origN since L4_BEGIN
	tst	counterJ , #1                   // 1 is (L2-1)
	ble	.Ldtrmm_kernel_L999 // done     // if zero, we don't need L1

	mov	pCRow0, pC			// pCRow0 = C
	add	pC , pC , LDC			// Update pC to point to next

#if defined(LEFT)
	mov	tempOffset, offset
#endif
	mov	pA, origPA			// pA = A

.Ldtrmm_kernel_L1_M12_BEGIN:

	mov	counterI, origM
	subs	counterI, counterI, #12		// M12
	blt	.Ldtrmm_kernel_L1_M8_BEGIN	// if counterI < 12

.Ldtrmm_kernel_L1_M12_20:

	INIT12x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	add	pA, pA, tempOffset, lsl #6	// 12 * sizeof(double) = (8+4)*8 = 2^6 + 2^5, where 12 is M12
	add	pA, pA, tempOffset, lsl #5	// pA moves 'tempOffset' pos. in its k dimension M12
	add	pB, pB, tempOffset, lsl #3	// #3 = 2^3 = 1 * sizeof( double ), where 1 is L1
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #12          // 12 is M12
#else
	add	tempK, tempOffset, #1           // 1 is L1
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Ldtrmm_kernel_L1_M12_40
	.align 5

.Ldtrmm_kernel_L1_M12_22:
	KERNEL12x1_SUB
	KERNEL12x1_SUB
	KERNEL12x1_SUB
	KERNEL12x1_SUB

	KERNEL12x1_SUB
	KERNEL12x1_SUB
	KERNEL12x1_SUB
	KERNEL12x1_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L1_M12_22

.Ldtrmm_kernel_L1_M12_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L1_M12_100       // no remainders, skip

.Ldtrmm_kernel_L1_M12_42:

	KERNEL12x1_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L1_M12_42

.Ldtrmm_kernel_L1_M12_100:

	SAVE12x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #12               // 12 is M12
#else
	sub	tempK, tempK, #1                // 1 is L1
#endif
	add	pA, pA, tempK, lsl #6		// 12 * sizeof(double) = (8+4)*8 = 2^6 + 2^5, where 12 is M12
	add	pA, pA, tempK, lsl #5		// pA moves 'tempOffset' pos. in its k dimension M12
	add	pB, pB, tempK, lsl #3           // #3 = 2^3 = 1 * sizeof(double), where 1 is L1
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #12     // 12 is M12
#endif

.Ldtrmm_kernel_L1_M12_END:

	subs	counterI, counterI, #1
	bgt	.Ldtrmm_kernel_L1_M12_20


.Ldtrmm_kernel_L1_M8_BEGIN:
	adds	counterI, counterI, #12		// counterI = origM % 12
	beq	.Ldtrmm_kernel_L1_END           // if zero, skip L1_M8/M4/M2/M1

	tst	counterI, #8                    // 8 is from M8, to test wether L1_M8 need to run
	ble	.Ldtrmm_kernel_L1_M4_BEGIN      // if zero, skip L1_M4, go test L1_M4

.Ldtrmm_kernel_L1_M8_20:

	INIT8x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	/* todo: combined lsl and add into one */
	lsl	temp, tempOffset, #3            // #3 = 2^3 = 1 * sizeof( double ), where 1 is L1
	add	pB, pB, temp                    // pB is moved 'tempOffset' pos. in its k dimension L4
	lsl	temp, tempOffset, #6            // #6 = 2^6 = 8 * sizeof( double ), where 8 is M8
	add	pA, pA, temp                    // pA is moved 'tempOffset' pos. in its k dimension M4
#endif
#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #8           // 8 is M8
#else
	add	tempK, tempOffset, #1           // 1 is L1
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Ldtrmm_kernel_L1_M8_40
	.align 5

.Ldtrmm_kernel_L1_M8_22:
	KERNEL8x1_SUB
	KERNEL8x1_SUB
	KERNEL8x1_SUB
	KERNEL8x1_SUB

	KERNEL8x1_SUB
	KERNEL8x1_SUB
	KERNEL8x1_SUB
	KERNEL8x1_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L1_M8_22


.Ldtrmm_kernel_L1_M8_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L1_M8_100

.Ldtrmm_kernel_L1_M8_42:

	KERNEL8x1_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L1_M8_42

.Ldtrmm_kernel_L1_M8_100:

	SAVE8x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #8                // 8 is M8
#else
	sub	tempK, tempK, #1                // 1 is L1
#endif
	/* todo: combined lsl and add into one */
	lsl	temp, tempK, #6                 // #6 = 2^6 = 8 * sizeof(double), where 8 is M8
	add	pA, pA, temp
	lsl	temp, tempK, #3                 // #3 = 2^3 = 1 * sizeof(double), where 1 is L1
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #8      // 8 is M8
#endif

.Ldtrmm_kernel_L1_M8_END:


.Ldtrmm_kernel_L1_M4_BEGIN:
						// counterI = origM % 12
	tst	counterI , #7                   // 7 is (M8-1)
	ble	.Ldtrmm_kernel_L1_END           // if zero, skip L1_M4/M2/M1

	tst	counterI, #4			// 4 is M4, test bit (#4)
	ble	.Ldtrmm_kernel_L1_M2_BEGIN      // if zero, skip M4; jump to M2

.Ldtrmm_kernel_L1_M4_20:

	INIT4x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	/* todo: combined lsl and add into one */
	lsl	temp, tempOffset, #3            // #3=2^3=1*sizeof(double), 1 is L1
	add	pB, pB, temp
	lsl	temp, tempOffset, #5            // #5=2^5=4*sizeof(double), 4 is M4
	add	pA, pA, temp
#endif
#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #4           // 4 is M4
#else
	add	tempK, tempOffset, #1           // 1 is L1
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Ldtrmm_kernel_L1_M4_40
	.align 5

.Ldtrmm_kernel_L1_M4_22:
	KERNEL4x1_SUB
	KERNEL4x1_SUB
	KERNEL4x1_SUB
	KERNEL4x1_SUB

	KERNEL4x1_SUB
	KERNEL4x1_SUB
	KERNEL4x1_SUB
	KERNEL4x1_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L1_M4_22


.Ldtrmm_kernel_L1_M4_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L1_M4_100

.Ldtrmm_kernel_L1_M4_42:

	KERNEL4x1_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L1_M4_42

.Ldtrmm_kernel_L1_M4_100:

	SAVE4x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #4                // 4 is M4
#else
	sub	tempK, tempK, #1                // 1 is L1
#endif
	/* todo: combined lsl and add into one */
	lsl	temp, tempK, #5                 // #5 = 2^5 = 4 * sizeof(double), where 4 is M4
	add	pA, pA, temp
	lsl	temp, tempK, #3                 // #3 = 2^3 = 1 * sizeof(double), where 1 is L1
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #4      // 4 is M4
#endif

.Ldtrmm_kernel_L1_M4_END:

.Ldtrmm_kernel_L1_M2_BEGIN:
						// counterI = origM % 12
	tst	counterI , #3                   // 3 is (M4-1)
	ble	.Ldtrmm_kernel_L1_END           // if zero, no need to run L1_M2/M1

	tst	counterI, #2			// 2 is M2, test bit (#2)
	ble	.Ldtrmm_kernel_L1_M1_BEGIN      // if zero, skip L1_M2, go test L1_M1

.Ldtrmm_kernel_L1_M2_20:

	INIT2x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	/* todo: lsl and add can be combined into one 'add' instruction */
	lsl	temp, tempOffset, #3            // #3 = 2^3 = 1 * sizeof(double), where 1 is L1
	add	pB, pB, temp
	lsl	temp, tempOffset, #4            // #4 = 2^4 = 2 * sizeof( double ), where 2 is M2
	add	pA, pA, temp
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #2           // 2 is M2
#else
	add	tempK, tempOffset, #1           // 1 is L1
#endif

	asr 	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Ldtrmm_kernel_L1_M2_40

.Ldtrmm_kernel_L1_M2_22:

	KERNEL2x1_SUB
	KERNEL2x1_SUB
	KERNEL2x1_SUB
	KERNEL2x1_SUB

	KERNEL2x1_SUB
	KERNEL2x1_SUB
	KERNEL2x1_SUB
	KERNEL2x1_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L1_M2_22


.Ldtrmm_kernel_L1_M2_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L1_M2_100

.Ldtrmm_kernel_L1_M2_42:

	KERNEL2x1_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L1_M2_42

.Ldtrmm_kernel_L1_M2_100:

	SAVE2x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #2                // 2 is M2
#else
	sub	tempK, tempK, #1                // 1 is L1
#endif
	/* todo: lsl and add can be combined into one 'add' instruction */
	lsl	temp, tempK, #4                 // #4 = 2^4 = 2 * sizeof(double), where 2 is M2
	add	pA, pA, temp
	lsl	temp, tempK, #3                 // #3 = 2^3 = 1 * sizeof(double), where 1 is L1
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #2      // 2 is M2
#endif

.Ldtrmm_kernel_L1_M2_END:


.Ldtrmm_kernel_L1_M1_BEGIN:
						// counterI = origM % 12
	tst	counterI, #1
	ble	.Ldtrmm_kernel_L1_END           // if zero, no need to run L1_M1

.Ldtrmm_kernel_L1_M1_20:

	INIT1x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	/* todo: lsl and add can be combined into one 'add' instruction */
	lsl	temp, tempOffset, #3            // #3 = 2^3 = 1 * sizeof( double ), where 1 is M1/L1
	add	pB, pB, temp
	add	pA, pA, temp
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #1           // 1 is M1
#else
	add	tempK, tempOffset, #1           // 1 is L1
#endif

	asr 	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Ldtrmm_kernel_L1_M1_40

.Ldtrmm_kernel_L1_M1_22:
	KERNEL1x1_SUB
	KERNEL1x1_SUB
	KERNEL1x1_SUB
	KERNEL1x1_SUB

	KERNEL1x1_SUB
	KERNEL1x1_SUB
	KERNEL1x1_SUB
	KERNEL1x1_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L1_M1_22


.Ldtrmm_kernel_L1_M1_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L1_M1_100

.Ldtrmm_kernel_L1_M1_42:

	KERNEL1x1_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L1_M1_42

.Ldtrmm_kernel_L1_M1_100:

	SAVE1x1


.Ldtrmm_kernel_L1_END:


.Ldtrmm_kernel_L999:
	mov	x0, #0				// set return value
	ldp	d8, d9, [sp, #(0 * 16)]
	ldp	d10, d11, [sp, #(1 * 16)]
	ldp	d12, d13, [sp, #(2 * 16)]
	ldp	d14, d15, [sp, #(3 * 16)]
	ldp	d16, d17, [sp, #(4 * 16)]
	ldp	x18, x19, [sp, #(5 * 16)]
	ldp	x20, x21, [sp, #(6 * 16)]
	ldp	x22, x23, [sp, #(7 * 16)]
	ldp	x24, x25, [sp, #(8 * 16)]
	ldp	x26, x27, [sp, #(9 * 16)]
	ldr	x28, [sp, #(10 * 16)]
	add	sp, sp, #(11*16)
	ret

	EPILOGUE

