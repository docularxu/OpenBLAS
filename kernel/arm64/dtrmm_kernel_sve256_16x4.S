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

/*                   X0          X1          X2          d0        X3        x4       x5           x6            x7*/
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
#define temp		x18   // must save
#define tempOffset	x19   // must save
#define tempK		x20   // must save
// x21 must save
// x22 must save
// x23 must save
// x24 must save
// x25 must save
// x26 must save
// x27 must save
// x28 must save
// x29 frame
// x30 link
// x31 sp

/* TODO: adjust these memory prefetch sizes */
#define A_PRE_SIZE	2560
#define B_PRE_SIZE	448
#define C_PRE_SIZE	128

// z0, z1, z2, z3: elements from one column of *a
// z4, z5        : elements from one row of *b
// z16 ~ z31     : elements for all 64 doubles of *c
// z12           : alpha, duplicated to each element of z12

// z6, z7, z8, z9: replacement vectors used in loop unroll, for elements of *a for next iteration.
// z10, z11      : replacement vectors used in loop unroll, for elements of *b for next iteration.
// z13 ~ z15     : not used, free.

// p0            : governing predicate register, all true

//v08 must save
//v09 must save
//v10 must save
//v11 must save
//v12 must save
//v13 must save
//v14 must save
//v15 must save
//v16 must save
//v17 must save

/*******************************************************************************
* Macro definitions
*******************************************************************************/

/******************************************************************************
 * L4
 ******************************************************************************/

/* set z16 ~ z31 to zero */
.macro INIT16x4
        mov     z16.d, #0
        mov     z17.d, #0
        mov     z18.d, #0
        mov     z19.d, #0

/* TODO: we can use 'dup z20.d, z16.d[0]' to initilize the other vectors.
 *       Question is: will that be more efficient than using 'mov z20.d, #0'?
 */
        mov     z20.d, #0
        mov     z21.d, #0
        mov     z22.d, #0
        mov     z23.d, #0

        mov     z24.d, #0
        mov     z25.d, #0
        mov     z26.d, #0
        mov     z27.d, #0

        mov     z28.d, #0
        mov     z29.d, #0
        mov     z30.d, #0
        mov     z31.d, #0
.endm

/* There are only three possible sequences of iterations: KERNEL16x4_I/M2/M1/E.
 *   Namely they are:
 *     Sequence 1: _I,  _M2, _M1, _M2
 *     Sequence 2: _M1, _M2, _M1, _M2
 *     Sequence 3: _M1, _M2, _M1, _E
 *
 * Seq.1 is the start;
 * Seq.2 follows, which can loop from 0 to many times;
 * Seq.3 is the end.
 */

/* Design logic:
 * _I is the first in the sequence.
 * _I is always followed by _M2.
 *
 * Input:
 *   pA: address of next element of *a
 *   pB: address of next element of *b
 * Vector registers:
 *   z16 ~ z31: c  (not initialized, don't care)
 *   z0 ~ z3  : a
 *   z4, z5   : b
 * Output:
 *   z16 ~ z31            : updated c
 *   z6, z7?, z8?, z9?    : next column of a
 *   z10, z11             : next row    of b
 * Note: z7, z8, z9 are to be loaded in next iteration.
 *
 * Note:
 *   fmul implies z16 ~ z31 can be any initial value, no need to zero out.
 *        That's different with using fmla.
 */
.macro KERNEL16x4_I
	/* load a */
	ldr  z0, [pA]                        // load a(  0:3, l )
	/* load b */
	ld1rqd  {z4.d}, p0/z, [pB]           // load b( l,0:1 )

	fmul  z16.d, z0.d, z4.d[0]           //

	ldr  z1, [pA, #1, MUL VL]            // load a(  4:7, l )
	fmul  z20.d, z0.d, z4.d[1]           //
	fmul  z17.d, z1.d, z4.d[0]           //

	ld1rqd  {z5.d}, p0/z, [pB, #16]      // load b( l,2:3 )
	fmul  z21.d, z1.d, z4.d[1]           //
	fmul  z24.d, z0.d, z5.d[0]           //
	fmul  z28.d, z0.d, z5.d[1]           //

	ldr  z2, [pA, #2, MUL VL]            // load a(  8:11,l )
	fmul  z25.d, z1.d, z5.d[0]           //
	fmul  z29.d, z1.d, z5.d[1]           //

	ldr  z3, [pA, #3, MUL VL]            // load a( 12:15,l )
	fmul  z18.d, z2.d, z4.d[0]           //
	fmul  z22.d, z2.d, z4.d[1]           //

	/* load next a */
	ldr  z6, [pA, #4, MUL VL]            // load a(  0:3, l+1 )
	fmul  z26.d, z2.d, z5.d[0]           //
	fmul  z30.d, z2.d, z5.d[1]           //

	/* load next b */
	ld1rqd  {z10.d}, p0/z, [pB, #32]      // load b( l+1,0:1 )
	fmul  z19.d, z3.d, z4.d[0]           //
	fmul  z23.d, z3.d, z4.d[1]           //

	ld1rqd  {z11.d}, p0/z, [pB, #48]      // load b( l+1,2:3 )
	fmul  z27.d, z3.d, z5.d[0]           //
	fmul  z31.d, z3.d, z5.d[1]           //

        /* Move the following to _M2 */
	// ldr  z7, [pA, #5, MUL VL]            // load a(  4:7, l+1 )
	// ldr  z8, [pA, #6, MUL VL]            // load a(  8:11,l+1 )
	// ldr  z9, [pA, #7, MUL VL]            // load a( 12:15,l+1 )

.endm

/*
 * Vector registers:
 *   z16 ~ z31         : c
 * Input
 *   z6, z7?, z8?, z9? : a
 *   z10, z11          : b
 * Note: z7, z8, z9 are to be loaded in this iteration.
 *
 * Output:
 *   z16 ~ z31         : updated c
 *   z0, z1, z2, z3?   : next column of a
 *   z4, z5            : next row    of b
 * Note: z3 will be load in the next iteration.
 */
.macro KERNEL16x4_M2
	ldr  z7, [pA, #5, MUL VL]             // load a(  4:7, l+1 )
	fmla  z16.d, z6.d, z10.d[0]           //
	fmla  z20.d, z6.d, z10.d[1]           //

	ldr  z8, [pA, #6, MUL VL]             // load a(  8:11,l+1 )
	fmla  z24.d, z6.d, z11.d[0]           //
	fmla  z28.d, z6.d, z11.d[1]           //

	ldr  z9, [pA, #7, MUL VL]             // load a( 12:15,l+1 )
	fmla  z17.d, z7.d, z10.d[0]           //
	fmla  z21.d, z7.d, z10.d[1]           //

	/* load next a */
	ldr  z0, [pA, #8,  MUL VL]            // load a(  0:3, l+2 )
	fmla  z25.d, z7.d, z11.d[0]           //
	fmla  z29.d, z7.d, z11.d[1]           //

	/* load next b */
	ld1rqd  {z4.d},  p0/z, [pB, #64]      // load b( l+2,0:1 )
	fmla  z18.d, z8.d, z10.d[0]           //
	fmla  z22.d, z8.d, z10.d[1]           //

	ld1rqd  {z5.d},  p0/z, [pB, #80]      // load b( l+2,2:3 )
	fmla  z26.d, z8.d, z11.d[0]           //
	fmla  z30.d, z8.d, z11.d[1]           //

	ldr  z1, [pA, #9,  MUL VL]            // load a(  4:7, l+2 )
	fmla  z19.d, z9.d, z10.d[0]           //
	fmla  z23.d, z9.d, z10.d[1]           //

	ldr  z2, [pA, #10, MUL VL]            // load a(  8:11,l+2 )
	fmla  z27.d, z9.d, z11.d[0]           //
	fmla  z31.d, z9.d, z11.d[1]           //
.endm

/* Design logic: 
 * _M1 follows _M2
 * _M1 is followed by either _M2 or _E. 
 *  So, _M1 needs to:
 *   1. ldr z3.
 *   2. do the accumulation. ('l+2')
 *   3. for loading next round of data, consider M2 is following M1.
 *      - so don't load z7, z8, z9 in M1
 */
/* Vector registers:
 *   z16 ~ z31            : c
 * Input
 *   z0, z1, z2, z3?      : a
 *   z4, z5               : b
 * Note: z3 is to be loaded in this iteration.
 *
 * Output:
 *   z16 ~ z31            : updated c
 *   z6, z7?, z8?, z9?    : next column of a
 *   z10, z11             : next row    of b
 *   pA: adjusted, so it can loop with _M2/_E
 *   pB: adjusted, so it can loop with _M2/_E
 */
 .macro KERNEL16x4_M1

	fmla  z16.d, z0.d, z4.d[0]           //
	fmla  z20.d, z0.d, z4.d[1]           //
	fmla  z24.d, z0.d, z5.d[0]           //
	fmla  z28.d, z0.d, z5.d[1]           //

	ldr  z3, [pA, #11, MUL VL]            // load a( 12:15,l+2 )

	fmla  z17.d, z1.d, z4.d[0]           //
	fmla  z21.d, z1.d, z4.d[1]           //

	/* load a */
	ldr  z6, [pA, #12, MUL VL]           // load a(  0:3, l+3 )
	fmla  z25.d, z1.d, z5.d[0]           //
	fmla  z29.d, z1.d, z5.d[1]           //

        /* adjust pA */
        add  pA, pA, #256                    // 256 = 8 * (MUL VL)
	                                     // next element of a has address: [pA, #5, MUL VL]
	fmla  z18.d, z2.d, z4.d[0]           //
	fmla  z22.d, z2.d, z4.d[1]           //

	/* load b */
	ld1rqd  {z10.d}, p0/z, [pB, #96]     // load b( l+3,0:1 )
	fmla  z26.d, z2.d, z5.d[0]           //
	fmla  z30.d, z2.d, z5.d[1]           //
 
	ld1rqd  {z11.d}, p0/z, [pB, #112]    // load b( l+3,2:3 )
	fmla  z19.d, z3.d, z4.d[0]           //
	fmla  z23.d, z3.d, z4.d[1]           //

        /* adjust pB */
        add  pB, pB, #64
	                                     // next element of b has address: [pB, #64]
	fmla  z27.d, z3.d, z5.d[0]           //
	fmla  z31.d, z3.d, z5.d[1]           //

        /* TO be loaded in next iteration, which could be either a _M2 or a _E */
	// ldr  z7, [pA, #13, MUL VL]           // load a(  4:7, l+3 )
	// ldr  z8, [pA, #14, MUL VL]           // load a(  8:11,l+3 )
	// ldr  z9, [pA, #15, MUL VL]           // load a( 12:15,l+3 )
 .endm
 
/* Design logic:
 * _E always follows _M1. And at the end of _E, pA and pB should point to
 *   next elements.
 * So it needs to:
 *   1. ldr z7, z8, z9, in the beginning. then,
 *   2. do the accumulation. ('l+3')
 *   3. adjust pA and pB to next to-be-calculated elements.
 *
 * Input:
 *   pA: address of loading from *a
 *   pB: address of loading from *b
 * Vector registers:
 *   z16 ~ z31         : c
 * Input
 *   z6, z7?, z8?, z9? : a
 *   z10, z11          : b
 * Note: z7, z8, z9 are to be loaded in this iteration.
 *
 * Output:
 *   z16 ~ z31         : updated c
 *   pA: pointing to next elements of *a
 *   pB: pointing to next elements of *b
 */ 
.macro KERNEL16x4_E
	ldr  z7, [pA, #5, MUL VL]             // load a(  4:7, l+3 )
	fmla  z16.d, z6.d, z10.d[0]           //
	fmla  z20.d, z6.d, z10.d[1]           //
	fmla  z24.d, z6.d, z11.d[0]           //
	fmla  z28.d, z6.d, z11.d[1]           //

	ldr  z8, [pA, #6, MUL VL]             // load a(  8:11,l+3 )
	fmla  z17.d, z7.d, z10.d[0]           //
	fmla  z21.d, z7.d, z10.d[1]           //
	fmla  z25.d, z7.d, z11.d[0]           //
	fmla  z29.d, z7.d, z11.d[1]           //

	ldr  z9, [pA, #7, MUL VL]             // load a( 12:15,l+3 )
	fmla  z18.d, z8.d, z10.d[0]           //
	fmla  z22.d, z8.d, z10.d[1]           //
	fmla  z26.d, z8.d, z11.d[0]           //
	fmla  z30.d, z8.d, z11.d[1]           //

        /* adjust pA */
        add pA, pA, #256                      // 256 = 8 * (MUL VL)
	fmla  z19.d, z9.d, z10.d[0]           //
	fmla  z23.d, z9.d, z10.d[1]           //

        /* adjust pB */
        add pB, pB, #64
	fmla  z27.d, z9.d, z11.d[0]           //
	fmla  z31.d, z9.d, z11.d[1]           //
 .endm

/* Input:
 *   pA: address of next element of *a
 *   pB: address of next element of *b
 * Vector registers:
 *   z16 ~ z31: c  (not initialized, don't care)
 *   z0 ~ z3  : a
 *   z4, z5   : b
 */
 .macro KERNEL16x4_SUB
	/* load a */
	ldr  z0, [pA]                        // load a(  0:3, l )
	/* load b */
	ld1rqd  {z4.d}, p0/z, [pB]           // load b( l,0:1 )
	fmla  z16.d, z0.d, z4.d[0]           //

	ld1rqd  {z5.d}, p0/z, [pB, #16]      // load b( l,2:3 )
	fmla  z20.d, z0.d, z4.d[1]           //
	fmla  z24.d, z0.d, z5.d[0]           //

	ldr  z1, [pA, #1, MUL VL]            // load a(  4:7, l )
	fmla  z28.d, z0.d, z5.d[1]           //
	fmla  z17.d, z1.d, z4.d[0]           //
	fmla  z21.d, z1.d, z4.d[1]           //

	ldr  z2, [pA, #2, MUL VL]            // load a(  8:11,l )
	fmla  z25.d, z1.d, z5.d[0]           //
	fmla  z29.d, z1.d, z5.d[1]           //

	ldr  z3, [pA, #3, MUL VL]            // load a( 12:15,l )
	fmla  z18.d, z2.d, z4.d[0]           //
	fmla  z22.d, z2.d, z4.d[1]           //
	fmla  z26.d, z2.d, z5.d[0]           //
	fmla  z30.d, z2.d, z5.d[1]           //

        /* adjust pA */
        add pA, pA, #128                     // 128 = 4 * (MUL VL) = 16 * sizeof(double)
	fmla  z19.d, z3.d, z4.d[0]           //
	fmla  z23.d, z3.d, z4.d[1]           //

        /* adjust pB */
        add pB, pB, #32                      // 32 = 4 * sizeof(double)
	fmla  z27.d, z3.d, z5.d[0]           //
	fmla  z31.d, z3.d, z5.d[1]           //

/* TODO: adjust prfm */
//	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
//	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
//	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
.endm

/* SAVE16x4
 * Input:
 *   z12        : alpha, and duplicated into each elements of z12
 *   z16 ~ z31  : c
 *   pCRow0, pCRow1, pCRow2, pCRow3:  pointer to next position in *c each Column
 * Progress:
 *   z0 ~ z7    : for storing alpha-multiplied result
 *   prfm L2    : prefetch before 'str'.
 *                Instead of ‘PLD', use 'PST' for store
 *
 */
.macro SAVE16x4
	prfm	PSTL2KEEP, [pCRow0, #C_PRE_SIZE]
	fmul	z0.d, z16.d, z12.d              // scale by alpha
	fmul	z1.d, z17.d, z12.d
	str 	z0, [pCRow0]                    // store column 0
	str 	z1, [pCRow0, #1, MUL VL]

	fmul	z2.d, z18.d, z12.d
	fmul	z3.d, z19.d, z12.d
	str 	z2, [pCRow0, #2, MUL VL]
	str 	z3, [pCRow0, #3, MUL VL]

	add	pCRow0, pCRow0, #128

	prfm	PSTL2KEEP, [pCRow1, #C_PRE_SIZE]
	fmul	z4.d, z20.d, z12.d              // scale by alpha
	fmul	z5.d, z21.d, z12.d
	str 	z4, [pCRow1]                    // store column 1
	str 	z5, [pCRow1, #1, MUL VL]

	fmul	z6.d, z22.d, z12.d
	fmul	z7.d, z23.d, z12.d
	str 	z6, [pCRow1, #2, MUL VL]
	str 	z7, [pCRow1, #3, MUL VL]

	add	pCRow1, pCRow1, #128

	prfm	PSTL2KEEP, [pCRow2, #C_PRE_SIZE]
	fmul	z0.d, z24.d, z12.d              // scale by alpha
	fmul	z1.d, z25.d, z12.d
	str 	z0, [pCRow2]                    // store column 2
	str 	z1, [pCRow2, #1, MUL VL]

	fmul	z2.d, z26.d, z12.d
	fmul	z3.d, z27.d, z12.d
	str 	z2, [pCRow2, #2, MUL VL]
	str 	z3, [pCRow2, #3, MUL VL]

	add	pCRow2, pCRow2, #128

	prfm	PSTL2KEEP, [pCRow3, #C_PRE_SIZE]
	fmul	z4.d, z28.d, z12.d              // scale by alpha
	fmul	z5.d, z29.d, z12.d
	str 	z4, [pCRow3]                    // store column 3
	str 	z5, [pCRow3, #1, MUL VL]

	fmul	z6.d, z30.d, z12.d
	fmul	z7.d, z31.d, z12.d
	str 	z6, [pCRow3, #2, MUL VL]
	str 	z7, [pCRow3, #3, MUL VL]

	add	pCRow3, pCRow3, #128
.endm

/******************************************************************************/

/* set
 *   z16, z17   : c( 0:7,0 )
 *   z20, z21   : c( 0:7,1 )
 *   z24, z25   : c( 0:7,2 )
 *   z28, z29   : c( 0:7,3 )
 * to zero
 */
.macro INIT8x4
        mov     z16.d, #0
        mov     z17.d, #0

/* TODO: we can use 'dup z20.d, z16.d[0]' to initilize the other vectors.
 *       Question is: will that be more efficient than using 'mov z20.d, #0'?
 */
        mov     z20.d, #0
        mov     z21.d, #0

        mov     z24.d, #0
        mov     z25.d, #0

        mov     z28.d, #0
        mov     z29.d, #0
.endm

/* Input:
 *   pA: address of next element of *a
 *   pB: address of next element of *b
 * Vector registers:
 *   z16, z17   : c( 0:7,0 )
 *   z20, z21   : c( 0:7,1 )
 *   z24, z25   : c( 0:7,2 )
 *   z28, z29   : c( 0:7,3 )
 *
 *   z0, z1     : a
 *   z4, z5     : b
 */
 .macro KERNEL8x4_SUB
	/* load a */
	ldr  z0, [pA]                        // load a(  0:3, l )
	/* load b */
	ld1rqd  {z4.d}, p0/z, [pB]           // load b( l,0:1 )
	fmla  z16.d, z0.d, z4.d[0]           //

	ld1rqd  {z5.d}, p0/z, [pB, #16]      // load b( l,2:3 )
	fmla  z20.d, z0.d, z4.d[1]           //
	fmla  z24.d, z0.d, z5.d[0]           //

	ldr  z1, [pA, #1, MUL VL]            // load a(  4:7, l )
	fmla  z28.d, z0.d, z5.d[1]           //
	fmla  z17.d, z1.d, z4.d[0]           //

        /* adjust pA */
        add pA, pA, #64                      // 64 = 2 * (MUL VL) = 8 * sizeof(double)
	fmla  z21.d, z1.d, z4.d[1]           //
	fmla  z25.d, z1.d, z5.d[0]           //

        /* adjust pB */
        add pB, pB, #32                      // 32 = 4 * sizeof(double)
	fmla  z29.d, z1.d, z5.d[1]           //
.endm

/* SAVE8x4
 * Input:
 *   z12        : alpha, and duplicated into each elements of z12
 *   z16, z17   : c( 0:7,0 )
 *   z20, z21   : c( 0:7,1 )
 *   z24, z25   : c( 0:7,2 )
 *   z28, z29   : c( 0:7,3 )
 *   pCRow0     : pointer to next position in *c Column 0
 * Progress:
 *   z0, z1, z4, z5
 *              : for storing alpha-multiplied result
 *   pCRow1, pCRow2, pCRow3
 *              : pointer to next position in *c Column 1/2/3
 * Output:
 *   pCRow0     : is adjusted to point to next
 */
.macro SAVE8x4
	fmul	z0.d, z16.d, z12.d              // scale by alpha
	fmul	z1.d, z17.d, z12.d
	str 	z0, [pCRow0]                    // store column 0
	str 	z1, [pCRow0, #1, MUL VL]

	add	pCRow1, pCRow0, LDC

	fmul	z4.d, z20.d, z12.d              // scale by alpha
	fmul	z5.d, z21.d, z12.d
	str 	z4, [pCRow1]                    // store column 1
	str 	z5, [pCRow1, #1, MUL VL]

	add	pCRow2, pCRow1, LDC

	fmul	z0.d, z24.d, z12.d              // scale by alpha
	fmul	z1.d, z25.d, z12.d
	str 	z0, [pCRow2]                    // store column 2
	str 	z1, [pCRow2, #1, MUL VL]

	add	pCRow3, pCRow2, LDC

	fmul	z4.d, z28.d, z12.d              // scale by alpha
	fmul	z5.d, z29.d, z12.d
	str 	z4, [pCRow3]                    // store column 3
	str 	z5, [pCRow3, #1, MUL VL]

	add	pCRow0, pCRow0, #64             // pCRow0 moves down
.endm

/******************************************************************************/

/* set
 *   z16   : c( 0:3,0 )
 *   z20   : c( 0:3,1 )
 *   z24   : c( 0:3,2 )
 *   z28   : c( 0:3,3 )
 * to zero
 */
.macro INIT4x4
        mov     z16.d, #0

/* TODO: we can use 'dup z20.d, z16.d[0]' to initilize the other vectors.
 *       Question is: will that be more efficient than using 'mov z20.d, #0'?
 */
        mov     z20.d, #0

        mov     z24.d, #0

        mov     z28.d, #0
.endm

/* Input:
 *   pA: address of next element of *a
 *   pB: address of next element of *b
 * Vector registers:
 *   z16   : c( 0:3,0 )
 *   z20   : c( 0:3,1 )
 *   z24   : c( 0:3,2 )
 *   z28   : c( 0:3,3 )
 *
 *   z0    : a
 *   z4, z5: b
 */
 .macro KERNEL4x4_SUB
	/* load a */
	ldr  z0, [pA]                        // load a(  0:3, l )
	/* load b */
	ld1rqd  {z4.d}, p0/z, [pB]           // load b( l,0:1 )

	fmla  z16.d, z0.d, z4.d[0]           //
	fmla  z20.d, z0.d, z4.d[1]           //
        /* adjust pA */
        add pA, pA, #32                      // 32 = 4 * sizeof(double)

	ld1rqd  {z5.d}, p0/z, [pB, #16]      // load b( l,2:3 )
	fmla  z24.d, z0.d, z5.d[0]           //

        /* adjust pB */
        add pB, pB, #32                      // 32 = 4 * sizeof(double)

	fmla  z28.d, z0.d, z5.d[1]           //
.endm

/* SAVE4x4
 * Input:
 *   z12   : alpha, and duplicated into each elements of z12
 *   z16   : c( 0:3,0 )
 *   z20   : c( 0:3,1 )
 *   z24   : c( 0:3,2 )
 *   z28   : c( 0:3,3 )
 *   pCRow0     : pointer to next position in *c Column 0
 * Progress:
 *   z0, z4     : for storing alpha-multiplied result
 *   pCRow1, pCRow2, pCRow3
 *              : pointer to next position in *c Column 1/2/3
 * Output:
 *   pCRow0     : is adjusted to point to next
 */
.macro SAVE4x4
	fmul	z0.d, z16.d, z12.d              // scale by alpha
	str 	z0, [pCRow0]                    // store column 0

	add	pCRow1, pCRow0, LDC

	fmul	z4.d, z20.d, z12.d              // scale by alpha
	str 	z4, [pCRow1]                    // store column 1

	add	pCRow2, pCRow1, LDC

	fmul	z0.d, z24.d, z12.d              // scale by alpha
	str 	z0, [pCRow2]                    // store column 2

	add	pCRow3, pCRow2, LDC

	fmul	z4.d, z28.d, z12.d              // scale by alpha
	str 	z4, [pCRow3]                    // store column 3

	add	pCRow0, pCRow0, #32            // pCRow0 moves down
.endm

/******************************************************************************/

/* set
 *   z16   : c( 0,0:3 )
 *   z20   : c( 1,0:3 )
 * to zero
 */
.macro INIT2x4
        mov     z16.d, #0
        mov     z20.d, #0
.endm

/* Input:
 *   pA: address of next element of *a
 *   pB: address of next element of *b
 *   p0    : predicate register, all true
 * Vector registers:
 *   z16   : c( 0,0:3 )
 *   z20   : c( 1,0:3 )
 *
 *   z0    : a( 0:1,l )  // col. l
 *   z4    : b( l,0:3 )  // row  l
 */
 .macro KERNEL2x4_SUB
	/* load a */
	ld1rqd  {z0.d}, p0/z, [pA]           // load a( 0,0:1 )
	/* load b */
	ldr  z4, [pB]                        // load b( l,0:3 )

	fmla  z16.d, z4.d, z0.d[0]           //
	fmla  z20.d, z4.d, z0.d[1]           //

        /* adjust pA */
        add pA, pA, #16                      // 16 = 2 * sizeof(double)
        /* adjust pB */
        add pB, pB, #32                      // 32 = 4 * sizeof(double)
.endm

/* SAVE2x4
 * Input:
 *   z12   : alpha, and duplicated into each elements of z12
 *   z16   : c( 0,0:3 )
 *   z20   : c( 1,0:3 )
 *   p0    : predicate register, all true
 *   pCRow0     : pointer to next position in *c Column 0
 * Progress:
 *   z0, z4     : for storing alpha-multiplied result
 *   z1         : index, scatter store addr.
 * Output:
 *   pCRow0     : is adjusted to point to next
 */
.macro SAVE2x4
	index	z1.d, pCRow0, LDC               // pointers of [pCRow0, pCRow1, pCRow2, pCRow3]

	fmul	z0.d, z16.d, z12.d              // scale by alpha
	fmul	z4.d, z20.d, z12.d              // scale by alpha

        /* scatter store */
	st1d	z0.d, p0, [z1.d]                // store c( 0,0:3 )
	st1d	z4.d, p0, [z1.d, #8]            // store c( 1,0:3 )

	add	pCRow0, pCRow0, #16             // pCRow0 moved down
.endm

/******************************************************************************/

/* set
 *   z16   : c( 0,0:3 )
 * to zero
 */
.macro INIT1x4
        mov     z16.d, #0
.endm

/* Input:
 *   pA: address of next element of *a
 *   pB: address of next element of *b
 *   p0    : predicate register, all true
 * Vector registers:
 *   z16     : c( 0,0:3 )
 *
 *   z0.d[0] : a( 0,l )                         // col. l
 *   z4      : b( l,0:3 )                       // row  l
 */
 .macro KERNEL1x4_SUB
	/* load a */
	ld1rd  {z0.d}, p0/z, [pA]                // load a( 0,0 )
	/* load b */
	ldr    z4, [pB]                        // load b( l,0:3 )

	fmla   z16.d, z4.d, z0.d[0]            //

        /* adjust pA */
        add    pA, pA, #8                      //  8 = 1 * sizeof(double)
        /* adjust pB */
        add    pB, pB, #32                     // 32 = 4 * sizeof(double)
.endm

/* SAVE1x4
 * Input:
 *   z12   : alpha, and duplicated into each elements of z12
 *   z16   : c( 0,0:3 )
 *   p0    : predicate register, all true
 *   pCRow0     : pointer to next position in *c Column 0
 * Progress:
 *   z0         : for storing alpha-multiplied result
 *   z1         : index, scatter store addr.
 * Output:
 *   pCRow0     : is adjusted to point to next
 */
.macro SAVE1x4
	index	z1.d, pCRow0, LDC               // pointers of [pCRow0, pCRow1, pCRow2, pCRow3]

	fmul	z0.d, z16.d, z12.d              // scale by alpha

        /* scatter store */
	st1d	z0.d, p0, [z1.d]                // store c( 0,0:3 )

	add	pCRow0, pCRow0, #8              // pCRow0 moved down
.endm

/******************************************************************************
 * L2
 ******************************************************************************/

/* set z16 ~ z23 to zero */
.macro INIT16x2
        mov     z16.d, #0
        mov     z17.d, #0
        mov     z18.d, #0
        mov     z19.d, #0

        mov     z20.d, #0
        mov     z21.d, #0
        mov     z22.d, #0
        mov     z23.d, #0
.endm

/* Input:
 *   pA: address of next element of *a
 *   pB: address of next element of *b
 * Vector registers:
 *   z16, z17, z18, z19: c( 0:15,0 )
 *   z20, z21, z22, z23: c( 0:15,1 )
 *   z0 ~ z3  : a
 *   z4       : b
 */
.macro KERNEL16x2_SUB
	/* load a */
	ldr  z0, [pA]                        // load a(  0:3, l )
	/* load b */
	ld1rqd  {z4.d}, p0/z, [pB]           // load b( l,0:1 )
        /* adjust pB */
        add pB, pB, #16                      // 16 = 2 * sizeof(double)

	fmla  z16.d, z0.d, z4.d[0]           //
	fmla  z20.d, z0.d, z4.d[1]           //

	ldr  z1, [pA, #1, MUL VL]            // load a(  4:7, l )
	fmla  z17.d, z1.d, z4.d[0]           //

	ldr  z2, [pA, #2, MUL VL]            // load a(  8:11,l )
	fmla  z21.d, z1.d, z4.d[1]           //
	fmla  z18.d, z2.d, z4.d[0]           //

	ldr  z3, [pA, #3, MUL VL]            // load a( 12:15,l )
	fmla  z22.d, z2.d, z4.d[1]           //

        /* adjust pA */
        add pA, pA, #128                     // 128 = 4 * (MUL VL) = 16 * sizeof(double)
	fmla  z19.d, z3.d, z4.d[0]           //
	fmla  z23.d, z3.d, z4.d[1]           //
.endm

/* SAVE16x2
 * Input:
 *   z12               : alpha, and duplicated into each elements of z12
 *   z16, z17, z18, z19: c( 0:15,0 )
 *   z20, z21, z22, z23: c( 0:15,1 )
 *   pCRow0     : pointer to next position in *c Column 0
 * Progress:
 *   z0 ~ z7    : for storing alpha-multiplied result
 */
.macro SAVE16x2
	fmul	z0.d, z16.d, z12.d              // scale by alpha
	fmul	z1.d, z17.d, z12.d
	str 	z0, [pCRow0]                    // store column 0
	str 	z1, [pCRow0, #1, MUL VL]

	add	pCRow1, pCRow0, LDC

	fmul	z2.d, z18.d, z12.d
	fmul	z3.d, z19.d, z12.d
	str 	z2, [pCRow0, #2, MUL VL]
	str 	z3, [pCRow0, #3, MUL VL]

	add	pCRow0, pCRow0, #128            // 128=4*(MUL VL)=16*sizeof(double)

	fmul	z4.d, z20.d, z12.d              // scale by alpha
	fmul	z5.d, z21.d, z12.d
	str 	z4, [pCRow1]                    // store column 1
	str 	z5, [pCRow1, #1, MUL VL]

	fmul	z6.d, z22.d, z12.d
	fmul	z7.d, z23.d, z12.d
	str 	z6, [pCRow1, #2, MUL VL]
	str 	z7, [pCRow1, #3, MUL VL]
.endm

/******************************************************************************/

/* set
 *   z16, z17   : c( 0:7,0 )
 *   z20, z21   : c( 0:7,1 )
 * to zero
 */
.macro INIT8x2
        mov     z16.d, #0
        mov     z17.d, #0

/* TODO: we can use 'dup z20.d, z16.d[0]' to initilize the other vectors.
 *       Question is: will that be more efficient than using 'mov z20.d, #0'?
 */
        mov     z20.d, #0
        mov     z21.d, #0
.endm

/* Input:
 *   pA: address of next element of *a
 *   pB: address of next element of *b
 * Vector registers:
 *   z16, z17   : c( 0:7,0 )
 *   z20, z21   : c( 0:7,1 )
 *
 *   z0, z1     : a
 *   z4         : b
 */
 .macro KERNEL8x2_SUB
	/* load a */
	ldr   z0, [pA]                       // load a(  0:3, l )
	/* load b */
	ld1rqd  {z4.d}, p0/z, [pB]           // load b( l,0:1 )
        /* adjust pB */
        add   pB, pB, #16                    // 16 = 2 * sizeof(double)

	fmla  z16.d, z0.d, z4.d[0]           //
	fmla  z20.d, z0.d, z4.d[1]           //

	ldr   z1, [pA, #1, MUL VL]           // load a(  4:7, l )
	fmla  z17.d, z1.d, z4.d[0]           //
        /* adjust pA */
        add   pA, pA, #64                    // 64 = 2 * (MUL VL) = 8 * sizeof(double)

	fmla  z21.d, z1.d, z4.d[1]           //
.endm

/* SAVE8x2
 * Input:
 *   z12        : alpha, and duplicated into each elements of z12
 *   z16, z17   : c( 0:7,0 )
 *   z20, z21   : c( 0:7,1 )
 *   pCRow0     : pointer to next position in *c Column 0
 * Progress:
 *   z0, z1, z4, z5
 *              : for storing alpha-multiplied result
 *   pCRow1     : pointer to next position in *c Column 1
 * Output:
 *   pCRow0     : is adjusted to point to next
 */
.macro SAVE8x2
	fmul	z0.d, z16.d, z12.d              // scale by alpha
	fmul	z1.d, z17.d, z12.d
	str 	z0, [pCRow0]                    // store column 0
	str 	z1, [pCRow0, #1, MUL VL]

	add	pCRow1, pCRow0, LDC

	fmul	z4.d, z20.d, z12.d              // scale by alpha
	fmul	z5.d, z21.d, z12.d
	str 	z4, [pCRow1]                    // store column 1
	str 	z5, [pCRow1, #1, MUL VL]

	add	pCRow0, pCRow0, #64             // pCRow0 moves down
.endm

/******************************************************************************/

/* set
 *   z16   : c( 0:3,0 )
 *   z20   : c( 0:3,1 )
 * to zero
 */
.macro INIT4x2
        mov     z16.d, #0
        mov     z20.d, #0
.endm

/* Input:
 *   pA: address of next element of *a
 *   pB: address of next element of *b
 * Vector registers:
 *   z16   : c( 0:3,0 )
 *   z20   : c( 0:3,1 )
 *   z0    : a
 *   z4    : b
 */
 .macro KERNEL4x2_SUB
	/* load a */
	ldr  z0, [pA]                        // load a(  0:3, l )
	/* load b */
	ld1rqd  {z4.d}, p0/z, [pB]           // load b( l,0:1 )

	fmla  z16.d, z0.d, z4.d[0]           //
	fmla  z20.d, z0.d, z4.d[1]           //

        /* adjust pA */
        add pA, pA, #32                      // 32 = 4 * sizeof(double)
        /* adjust pB */
        add pB, pB, #16                      // 32 = 2 * sizeof(double)
.endm

/* SAVE4x2
 * Input:
 *   z12   : alpha, and duplicated into each elements of z12
 *   z16   : c( 0:3,0 )
 *   z20   : c( 0:3,1 )
 *   pCRow0     : pointer to next position in *c Column 0
 * Progress:
 *   z0, z4     : for storing alpha-multiplied result
 *   pCRow1
 *              : pointer to next position in *c Column 1/2/3
 * Output:
 *   pCRow0     : is adjusted to point to next
 */
.macro SAVE4x2
	fmul	z0.d, z16.d, z12.d              // scale by alpha
	str 	z0, [pCRow0]                    // store column 0

	add	pCRow1, pCRow0, LDC

	fmul	z4.d, z20.d, z12.d              // scale by alpha
	str 	z4, [pCRow1]                    // store column 1

	add	pCRow0, pCRow0, #32            // pCRow0 moves down
.endm

/******************************************************************************/
/* Note: macros 2x2 and 1x2 references code from arm64 SIMD&FP.
 *       Due to small data size, there is no added value to use SVE instructions.
 */

.macro INIT2x2
	dup	v16.2d, xzr
	dup	v20.2d, xzr
.endm

/* Vector registers:
 *   v16.2d, v20.2d: c
 *   v0.2d         : a
 *   v4.2d         : b
 */
.macro KERNEL2x2_SUB
	ld1	{v4.2d}, [pB]
	add	pB, pB, #16

	ld1	{v0.2d}, [pA]
	add	pA, pA, #16

	fmla	v16.2d, v0.2d, v4.d[0]
	fmla	v20.2d, v0.2d, v4.d[1]
.endm

/* Vector registers:
 *   v16.2d, v20.2d: c
 *   z12           : duplicated alpha
 *   v0.2d,  v4.2d : elements multiplied by alpha
 */
.macro SAVE2x2
	fmul	v0.2d, v16.2d, v12.2d
	st1	{v0.2d}, [pCRow0]

	add	pCRow1 , pCRow0, LDC

	fmul	v4.2d, v20.2d, v12.2d
	st1	{v4.2d}, [pCRow1]

	add	pCRow0, pCRow0, #16
.endm

/******************************************************************************/
/* Note: macros 2x2 and 1x2 references code from arm64 SIMD&FP.
 *       Due to small data size, there is no added value to use SVE instructions.
 */

.macro INIT1x2
	dup	v16.2d, xzr
.endm

/* Vector registers:
 *   v16.2d        : c
 *   d0            : a
 *   v4.2d         : b
 */
.macro KERNEL1x2_SUB
	ld1	{v4.2d} , [pB]
	add	pB , pB, #16

	ldr	d0 , [pA]
	add	pA, pA, #8

	fmla	v16.2d, v4.2d, v0.d[0]
.endm

/* Vector registers:
 *   v16.2d        : c
 *   z12           : duplicated alpha
 *   v0.2d         : elements multiplied by alpha
 */
.macro SAVE1x2
	add	pCRow1 , pCRow0, LDC

	fmul	v0.2d, v16.2d, v12.2d
	st1	{v0.d}[0], [pCRow0]
	st1	{v0.d}[1], [pCRow1]

	add	pCRow0, pCRow0, #8
.endm

/******************************************************************************
 * L1
 ******************************************************************************/

/* set z16 ~ z19 to zero */
.macro INIT16x1
        mov     z16.d, #0
        mov     z17.d, #0
        mov     z18.d, #0
        mov     z19.d, #0
.endm

/* Input:
 *   pA: address of next element of *a
 *   pB: address of next element of *b
 * Vector registers:
 *   z16, z17, z18, z19: c( 0:15,0 )
 *   z0 ~ z3  : a
 *   z4       : b
 */
.macro KERNEL16x1_SUB
	ldr  	z0, [pA]                        // load a(  0:3, l )
	ld1rd 	{z4.d}, p0/z, [pB]              // load b( l,0 )
	add	pB , pB, #8
	fmla 	z16.d, p0/m, z0.d, z4.d         //

	ldr  	z1, [pA, #1, MUL VL]            // load a(  4:7, l )
	fmla 	z17.d, p0/m, z1.d, z4.d         //

	ldr  	z2, [pA, #2, MUL VL]            // load a(  8:11,l )
	ldr  	z3, [pA, #3, MUL VL]            // load a( 12:15,l )
	add	pA , pA, #128                   // 16 * sizeof( double )

	fmla 	z18.d, p0/m, z2.d, z4.d         //
	fmla 	z19.d, p0/m, z3.d, z4.d         //
.endm

/* SAVE16x1
 * Input:
 *   z12               : alpha, and duplicated into each elements of z12
 *   z16, z17, z18, z19: c( 0:15,0 )
 *   pCRow0     : pointer to next position in *c Column 0
 * Progress:
 *   z0 ~ z3    : for storing alpha-multiplied result
 */
.macro SAVE16x1
	fmul	z0.d, z16.d, z12.d              // scale by alpha
	fmul	z1.d, z17.d, z12.d
	str 	z0, [pCRow0]                    // store column 0
	str 	z1, [pCRow0, #1, MUL VL]

	fmul	z2.d, z18.d, z12.d
	fmul	z3.d, z19.d, z12.d
	str 	z2, [pCRow0, #2, MUL VL]
	str 	z3, [pCRow0, #3, MUL VL]

	add	pCRow0, pCRow0, #128            // 16 * sizeof( double )
.endm

/******************************************************************************/

/* set
 *   z16, z17   : c( 0:7,0 )
 * to zero
 */
.macro INIT8x1
        mov     z16.d, #0
        mov     z17.d, #0
.endm

/* Input:
 *   pA: address of next element of *a
 *   pB: address of next element of *b
 * Vector registers:
 *   z16, z17   : c( 0:7,0 )
 *
 *   z0, z1     : a
 *   z4         : b
 */
 .macro KERNEL8x1_SUB
	/* load a */
	ldr   z0, [pA]                       // load a(  0:3, l )
	/* load b */
	ld1rd 	{z4.d}, p0/z, [pB]           // load b( l,0 )
	ldr   z1, [pA, #1, MUL VL]           // load a(  4:7, l )

	fmla  z16.d, p0/m, z0.d, z4.d        //
        /* adjust pB */
        add   pB, pB, #8                     // 8 = 1 * sizeof(double)

	fmla  z17.d, p0/m, z1.d, z4.d        //
        /* adjust pA */
        add   pA, pA, #64                    // 64 = 2 * (MUL VL) = 8 * sizeof(double)
.endm

/* SAVE8x1
 * Input:
 *   z12        : alpha, and duplicated into each elements of z12
 *   z16, z17   : c( 0:7,0 )
 *   pCRow0     : pointer to next position in *c Column 0
 * Progress:
 *   z0, z1     : for storing alpha-multiplied result
 */
.macro SAVE8x1
	fmul	z0.d, z16.d, z12.d              // scale by alpha
	fmul	z1.d, z17.d, z12.d
	str 	z0, [pCRow0]                    // store column 0
	str 	z1, [pCRow0, #1, MUL VL]

	add	pCRow0, pCRow0, #64             // 8 * sizeof( double )
.endm

/******************************************************************************/

/* set
 *   z16        : c( 0:3,0 )
 * to zero
 */
.macro INIT4x1
        mov     z16.d, #0
.endm

/* Input:
 *   pA: address of next element of *a
 *   pB: address of next element of *b
 * Vector registers:
 *   z16        : c( 0:3,0 )
 *
 *   z0         : a
 *   z4         : b
 */
 .macro KERNEL4x1_SUB
	/* load a */
	ldr   z0, [pA]                       // load a(  0:3, l )
	/* load b */
	ld1rd 	{z4.d}, p0/z, [pB]           // load b( l,0 )

	fmla  z16.d, p0/m, z0.d, z4.d        //
        /* adjust pA */
        add   pA, pA, #32                    // 32 = 4 * sizeof(double)
        /* adjust pB */
        add   pB, pB, #8                     // 8 = 1 * sizeof(double)

.endm

/* SAVE4x1
 * Input:
 *   z12        : alpha, and duplicated into each elements of z12
 *   z16        : c( 0:3,0 )
 *   pCRow0     : pointer to next position in *c Column 0
 * Progress:
 *   z1         : for storing alpha-multiplied result
 */
.macro SAVE4x1
	fmul	z1.d, z16.d, z12.d              // scale by alpha
	str 	z1, [pCRow0]                    // store column 0

	add	pCRow0, pCRow0, #32             // 4 * sizeof( double )
.endm

/******************************************************************************/
/* Note: macros 2x1 and 1x1 references code from arm64 SIMD&FP.
 *       Due to small data size, there is no added value to use SVE instructions.
 */

.macro INIT2x1
	dup	v16.2d, xzr
.endm

/* Vector registers:
 *   v16.2d        : c
 *   v0.2d         : a
 *   d8            : b
 */
.macro KERNEL2x1_SUB
	ldr	d8, [pB]
	add	pB , pB, #8

	ld1	{v0.2d}, [pA]
	add	pA , pA, #16

	fmla	v16.2d, v0.2d, v8.d[0]
.endm

/* Vector registers:
 *   v16.2d        : c
 *   z12           : duplicated alpha
 *   v8.2d         : elements multiplied by alpha
 */
.macro SAVE2x1
	fmul	v8.2d, v16.2d, v12.2d
	st1	{v8.2d}, [pCRow0]

	add	pCRow0, pCRow0, #16             // 16 = 2 * sizeof(double)
.endm

/******************************************************************************/

.macro INIT1x1
	fmov	d16, xzr
.endm

/* Vector registers:
 *   d16        : c
 *   d0         : a
 *   d8         : b
 */
.macro KERNEL1x1_SUB
	ldr	d8, [pB]
	add	pB , pB, #8

	ldr	d0, [pA]
	add	pA , pA, #8

	fmadd 	d16, d0, d8, d16  
.endm

/* Vector registers:
 *   d16        : c
 *   d12        : duplicated alpha
 *   d8         : elements multiplied by alpha
 */
.macro SAVE1x1
	fmul	d8, d16, d12
	str 	d8, [pCRow0]

	add	pCRow0, pCRow0, #8              // 8 = 1 * sizeof(double)
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

	dup	z12.d, z0.d[0]                    // z12 is duplicated alpha
        ptrue   p0.d, all

	lsl	LDC, LDC, #3			// ldc = ldc * 8

#if !defined(LEFT)
	neg	tempOffset, offset
#endif
	mov	pB, origPB

	mov	counterJ, origN
	asr 	counterJ, counterJ, #2		// J = J / 4
	cmp 	counterJ, #0
	ble	.Ldtrmm_kernel_L2_BEGIN

/******************************************************************************/

.Ldtrmm_kernel_L4_BEGIN:
	mov	pCRow0, pC
	add	pCRow1, pCRow0, LDC
	add	pCRow2, pCRow1, LDC
	add	pCRow3, pCRow2, LDC

	add	pC, pCRow3, LDC


#if defined(LEFT)
	mov	tempOffset, offset
#endif
	mov	pA, origPA			// pA = start of A array

.Ldtrmm_kernel_L4_M16_BEGIN:

	mov	counterI, origM
	asr 	counterI, counterI, #4		// counterI = counterI / 16
	cmp 	counterI, #0
	ble	.Ldtrmm_kernel_L4_M8_BEGIN

	.align 5
.Ldtrmm_kernel_L4_M16_20:

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #7    // #7 = 2^7 = 16 * sizeof(double), where 16 is M16
	add	pA, pA, temp            // pA moves 'tempOffset' pos. in its k dimension M8
	lsl	temp, tempOffset, #5    // #5 = 2^5 = 4 * sizeof( double ), where 4 is L4
	add	pB, pB, temp            // pB moves 'tempOffset' pos. in its k dimension, L4
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset        // tempK is the number of non-zero rows/cols.
#elif defined(LEFT)                             //  defined(LEFT) &&  defined(TRANSA)
	add	tempK, tempOffset, #16          // 16 is M16
#else                                           // !defined(LEFT) && !defined(TRANSA)
	add	tempK, tempOffset, #4           // 4 is L4
#endif

	asr 	counterL , tempK, #3		// L = K / 8, unroll at 8
	cmp	counterL , #2			// is there at least 4 to do?
	blt	.Ldtrmm_kernel_L4_M16_32

	KERNEL16x4_I				// do one in the K
	KERNEL16x4_M2				// do another in the K
	KERNEL16x4_M1
	KERNEL16x4_M2
	KERNEL16x4_M1
	KERNEL16x4_M2
	KERNEL16x4_M1
	KERNEL16x4_M2

	subs	counterL, counterL, #2		// subtract 2
	ble	.Ldtrmm_kernel_L4_M16_22a

	.align 5
.Ldtrmm_kernel_L4_M16_22:

	KERNEL16x4_M1
	KERNEL16x4_M2
	KERNEL16x4_M1
	KERNEL16x4_M2
	KERNEL16x4_M1
	KERNEL16x4_M2
	KERNEL16x4_M1
	KERNEL16x4_M2

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L4_M16_22

	.align 5
.Ldtrmm_kernel_L4_M16_22a:

	KERNEL16x4_M1
	KERNEL16x4_M2
	KERNEL16x4_M1
	KERNEL16x4_M2
	KERNEL16x4_M1
	KERNEL16x4_M2
	KERNEL16x4_M1
	KERNEL16x4_E

	b	 .Ldtrmm_kernel_L4_M16_44

	.align 5
.Ldtrmm_kernel_L4_M16_32:

	tst	counterL, #1
	ble	.Ldtrmm_kernel_L4_M16_40

	KERNEL16x4_I
	KERNEL16x4_M2
	KERNEL16x4_M1
	KERNEL16x4_M2
	KERNEL16x4_M1
	KERNEL16x4_M2
	KERNEL16x4_M1
	KERNEL16x4_E

	b	.Ldtrmm_kernel_L4_M16_44

.Ldtrmm_kernel_L4_M16_40:

	INIT16x4

.Ldtrmm_kernel_L4_M16_44:

	ands	counterL , tempK, #7
	ble	.Ldtrmm_kernel_L4_M16_100

	.align 5
.Ldtrmm_kernel_L4_M16_46:

	KERNEL16x4_SUB

	subs	counterL, counterL, #1
	bne	.Ldtrmm_kernel_L4_M16_46

.Ldtrmm_kernel_L4_M16_100:

	SAVE16x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)                               //  defined(LEFT) && defined(TRANSA)
	sub	tempK, tempK, #16               // 16 is M16
#else                                           // !defined(LEFT) && !defined(TRANSA)
	sub	tempK, tempK, #4                // 4 is L4
#endif
	lsl	temp, tempK, #7                 // #7 = 2^7 = 16 * sizeof(double), where 16 is M16
	add	pA, pA, temp
	lsl	temp, tempK, #5                 // #5 = 2^5 = 4 * sizeof(double), where 4 is L4
	add	pB, pB, temp
#endif

#if defined(LEFT)
	add	tempOffset, tempOffset, #16     // 16 is M16
#endif
/* TODO: adjust prfm
 * 	prfm	PLDL1KEEP, [pA]
 * 	prfm	PLDL1KEEP, [pA, #64]
 * 	prfm	PLDL1KEEP, [origPB]
 */
.Ldtrmm_kernel_L4_M16_END:
	subs	counterI, counterI, #1
	bne	.Ldtrmm_kernel_L4_M16_20

.Ldtrmm_kernel_L4_M8_BEGIN:

	mov	counterI, origM
	tst	counterI , #15                  // 15 is from (M16-1), to test whehter (origM % 16) == 0
	ble	.Ldtrmm_kernel_L4_END           // if zero, no need to run L4_M8/M4/M2/M1

	tst	counterI, #8                    // 8 is from M8, to test wether L4_M8 need to run
	ble	.Ldtrmm_kernel_L4_M4_BEGIN      // if zero, skip L4_M8, go test next bit (#4) in L4_M4

.Ldtrmm_kernel_L4_M8_20:

	INIT8x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #5            // #5 = 2^5 = 4 * sizeof( double ), where 4 is L4
	add	pB, pB, temp                    // pB is moved 'tempOffset' pos. in its k dimension L4
	lsl	temp, tempOffset, #6            // #6 = 2^6 = 8 * sizeof( double ), where 8 is M8
	add	pA, pA, temp                    // pA is moved 'tempOffset' pos. in its k dimension M4
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #8           // 8 is M8
#else
	add	tempK, tempOffset, #4           // 4 is L4
#endif

	asr 	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Ldtrmm_kernel_L4_M8_40

.Ldtrmm_kernel_L4_M8_22:

/* TODO: it worths to write below kernels,
 *          unroll them, and interleave fops & memops
 */
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
	sub	tempK, tempK, #4                // 4 is L4
#endif
	lsl	temp, tempK, #6                 // #6 = 2^6 = 8 * sizeof(double), where 8 is M8
	add	pA, pA, temp
	lsl	temp, tempK, #5                 // #5 = 2^5 = 4 * sizeof(double), where 4 is L4
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #8      // 8 is M8
#endif

.Ldtrmm_kernel_L4_M8_END:


.Ldtrmm_kernel_L4_M4_BEGIN:

	mov	counterI, origM
	tst	counterI , #7                   // 7 is (M8-1)
	ble	.Ldtrmm_kernel_L4_END           // if zero, no need to run L4_M4/M2/M1

	tst	counterI, #4	        	// 4 is M4, to test whether L4_M4 need to run
	ble	.Ldtrmm_kernel_L4_M2_BEGIN      // if zero, skip L4_M4, go test next bit (#2) in L4_M2

.Ldtrmm_kernel_L4_M4_20:

	INIT4x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #5            // #5 = 2^5 = 4 * sizeof( double ), where 4 is L4
	                                        // #5 = 2^5 = 4 * sizeof( double ), where 4 is M4
	add	pB, pB, temp                    // pB is moved 'tempOffset' pos. in its k dimension L4
	add	pA, pA, temp                    // pA is moved 'tempOffset' pos. in its k dimension M4
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #4           // 4 is M4
#else
	add	tempK, tempOffset, #4           // 4 is L4
#endif
	asr 	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Ldtrmm_kernel_L4_M4_40

.Ldtrmm_kernel_L4_M4_22:

/* TODO: to be unrolled, for better fops/memops interleaving */
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
	lsl	temp, tempK, #5                 // #5 = 2^5 = 4 * sizeof(double), where 4 is L4
	                                        // #5 = 2^5 = 4 * sizeof(double), where 4 is M4
	add	pA, pA, temp
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #4      // 4 is M4
#endif

.Ldtrmm_kernel_L4_M4_END:

.Ldtrmm_kernel_L4_M2_BEGIN:

	mov	counterI, origM
	tst	counterI , #3                   // 3 is (M4-1)
	ble	.Ldtrmm_kernel_L4_END           // if zero, no need to run L4_M2/M1

	tst	counterI, #2			// 2 is M2
	ble	.Ldtrmm_kernel_L4_M1_BEGIN      // if zero, skip L4_M2, go test next bit (#1) for L4_M1

.Ldtrmm_kernel_L4_M2_20:

	INIT2x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #4            // #4 = 2^4 = 2 * sizeof( double ), where 2 is M2
	add	pA, pA, temp                    // pA is moved M2 in its k dimension
	lsl	temp, tempOffset, #5            // #5 = 2^5 = 4 * sizeof( double ), where 4 is L4
	add	pB, pB, temp                    // pB is moved L4 in its k dimension
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #2           // 2 is M2
#else
	add	tempK, tempOffset, #4           // 4 is L4
#endif
	asr 	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
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
	lsl	temp, tempK, #4                 // #4 = 2^4 = 2 * sizeof(double), where 2 is M2
	add	pA, pA, temp
	lsl	temp, tempK, #5                 // #5 = 2^5 = 4 * sizeof(double), where 4 is L4
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #2      // 2 is M2
#endif

.Ldtrmm_kernel_L4_M2_END:


.Ldtrmm_kernel_L4_M1_BEGIN:

	tst	counterI, #1			// counterI = counterI % 2
	ble	.Ldtrmm_kernel_L4_END           // if zero, no need to run L2_M1

.Ldtrmm_kernel_L4_M1_20:

	INIT1x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #3            // #3 = 2^3 = 1 * sizeof( double ), where 1 is M1
	add	pA, pA, temp                    // pA is moved M2 in its k dimension
	lsl	temp, tempOffset, #5            // #5 = 2^5 = 4 * sizeof( double ), where 4 is L4
	add	pB, pB, temp                    // pB is moved L4 in its k dimension
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #1           // 1 is M1
#else
	add	tempK, tempOffset, #4           // 4 is L4
#endif

	asr 	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
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
	bgt	.Ldtrmm_kernel_L4_M1_22


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
	lsl	temp, tempK, #3                 // #3 = 2^3 = 1 * sizeof(double), where 1 is M1
	add	pA, pA, temp
	lsl	temp, tempK, #5                 // #5 = 2^5 = 4 * sizeof(double), where 4 is L4
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #1      // 1 is M1
#endif

.Ldtrmm_kernel_L4_END:

	lsl	temp, origK, #5                 // #5 =2^5 = 4 * sizeof(double), where 4 is L4
	add	origPB, origPB, temp		// B = B + K * 4 * 8

#if !defined(LEFT)
	add	tempOffset, tempOffset, #4      // 4 is L4
#endif

	subs	counterJ, counterJ , #1		// j--
	bgt	.Ldtrmm_kernel_L4_BEGIN


/******************************************************************************/

.Ldtrmm_kernel_L2_BEGIN:   // less than 4 left in N direction

	mov	counterJ , origN
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

.Ldtrmm_kernel_L2_M16_BEGIN:

	mov	counterI, origM
	asr 	counterI, counterI, #4		// counterI = counterI / 16
	cmp	counterI, #0
	ble	.Ldtrmm_kernel_L2_M8_BEGIN      // if zero, jump to M8

.Ldtrmm_kernel_L2_M16_20:

	INIT16x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #7    // #7 = 2^7 = 16 * sizeof(double), where 16 is M16
	add	pA, pA, temp            // pA moves 'tempOffset' pos. in its k dimension M8
	lsl	temp, tempOffset, #4    // #4 = 2^4 = 2 * sizeof( double ), where 2 is L2
	add	pB, pB, temp            // pB moves 'tempOffset' pos. in its k dimension, L4
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #16  // 16 is M16
#else
	add	tempK, tempOffset, #2   // 2 is L2
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8, unroll at 8
	cmp	counterL,#0
	ble	.Ldtrmm_kernel_L2_M16_40
	.align 5

.Ldtrmm_kernel_L2_M16_22:

/* TODO: optimization by interleaving fops/memops instructions */
	KERNEL16x2_SUB
	KERNEL16x2_SUB
	KERNEL16x2_SUB
	KERNEL16x2_SUB

	KERNEL16x2_SUB
	KERNEL16x2_SUB
	KERNEL16x2_SUB
	KERNEL16x2_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L2_M16_22


.Ldtrmm_kernel_L2_M16_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L2_M16_100       // remainder is zero, finish

.Ldtrmm_kernel_L2_M16_42:

	KERNEL16x2_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L2_M16_42

.Ldtrmm_kernel_L2_M16_100:

	SAVE16x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #16               // 16 is M16
#else
	sub	tempK, tempK, #2                // 2 is L2
#endif
	lsl	temp, tempK, #7                 // #7 = 2^7 = 16 * sizeof(double), where 16 is M16
	add	pA, pA, temp
	lsl	temp, tempK, #4                 // #4 = 2^4 = 2 * sizeof(double), where 2 is L2
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #16     // 16 is M16
#endif

.Ldtrmm_kernel_L2_M16_END:
	subs	counterI, counterI, #1
	bgt	.Ldtrmm_kernel_L2_M16_20


.Ldtrmm_kernel_L2_M8_BEGIN:

	mov	counterI, origM
	tst	counterI , #15                  // 15 is from (M16-1)
	ble	.Ldtrmm_kernel_L2_END           // if zero, skip L2_M8/M4/M2/M1

	tst	counterI, #8			// 8 is M8, test bit
	ble	.Ldtrmm_kernel_L2_M4_BEGIN      // if zero, skip L2_M8, go test L2_M4

.Ldtrmm_kernel_L2_M8_20:

	INIT8x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
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

	mov	counterI, origM
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

	mov	counterI, origM
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
	lsl	temp, tempK, #4                 // #4 = 2^4 = 2 * sizeof(double), where 2 is M2/L2
	add	pA, pA, temp
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #2      // 2 is M2
#endif

.Ldtrmm_kernel_L2_M2_END:


.Ldtrmm_kernel_L2_M1_BEGIN:

	tst	counterI, #1			// counterI = counterI % 2
	ble	.Ldtrmm_kernel_L2_END           // if zero, done L2.

.Ldtrmm_kernel_L2_M1_20:

	INIT1x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
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
	add	origPB, origPB, origK, lsl #4	// B = B + K * 2 * 8

/******************************************************************************/

.Ldtrmm_kernel_L1_BEGIN:

	mov	counterJ , origN
	tst	counterJ , #1                   // 1 is (L2-1)
	ble	.Ldtrmm_kernel_L999 // done     // if zero, we don't need L1

	mov	pCRow0, pC			// pCRow0 = C
	add	pC , pC , LDC			// Update pC to point to next

#if defined(LEFT)
	mov	tempOffset, offset
#endif
	mov	pA, origPA			// pA = A

.Ldtrmm_kernel_L1_M16_BEGIN:

	mov	counterI, origM
	asr	counterI, counterI, #4		// counterI = counterI / 16
	cmp	counterI, #0
	ble	.Ldtrmm_kernel_L1_M8_BEGIN      // if zero, skip M16 and jump to L1_M8

.Ldtrmm_kernel_L1_M16_20:

	INIT16x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #7    // #7 = 2^7 = 16 * sizeof(double), where 16 is M16
	add	pA, pA, temp            // pA moves 'tempOffset' pos. in its k dimension M8
	lsl	temp, tempOffset, #3    // #3 = 2^3 = 1 * sizeof( double ), where 1 is L1
	add	pB, pB, temp            // pB moves 'tempOffset' pos. in its k dimension, L4
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #16          // 16 is M16
#else
	add	tempK, tempOffset, #1           // 1 is L1
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Ldtrmm_kernel_L1_M16_40
	.align 5

.Ldtrmm_kernel_L1_M16_22:
	KERNEL16x1_SUB
	KERNEL16x1_SUB
	KERNEL16x1_SUB
	KERNEL16x1_SUB

	KERNEL16x1_SUB
	KERNEL16x1_SUB
	KERNEL16x1_SUB
	KERNEL16x1_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L1_M16_22


.Ldtrmm_kernel_L1_M16_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Ldtrmm_kernel_L1_M16_100       // no remainders, skip

.Ldtrmm_kernel_L1_M16_42:

	KERNEL16x1_SUB

	subs	counterL, counterL, #1
	bgt	.Ldtrmm_kernel_L1_M16_42

.Ldtrmm_kernel_L1_M16_100:

	SAVE16x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #16               // 16 is M16
#else
	sub	tempK, tempK, #1                // 1 is L1
#endif
	lsl	temp, tempK, #7                 // #7 = 2^7 = 16 * sizeof(double), where 16 is M16
	add	pA, pA, temp
	lsl	temp, tempK, #3                 // #3 = 2^3 = 1 * sizeof(double), where 1 is L1
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #16     // 16 is M16
#endif

.Ldtrmm_kernel_L1_M16_END:

	subs	counterI, counterI, #1
	bgt	.Ldtrmm_kernel_L1_M16_20


.Ldtrmm_kernel_L1_M8_BEGIN:

	mov	counterI, origM
	tst	counterI , #15                  // 15 is (M16-1)
	ble	.Ldtrmm_kernel_L1_END           // if zero, skip L1_M8/M4/M2/M1

	tst	counterI, #8                    // 8 is from M8, to test wether L1_M8 need to run
	ble	.Ldtrmm_kernel_L1_M4_BEGIN      // if zero, skip L1_M4, go test L1_M4

.Ldtrmm_kernel_L1_M8_20:

	INIT8x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
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

	mov	counterI, origM
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

	mov	counterI, origM
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

	tst	counterI, #1			// counterI = counterI % 2
	ble	.Ldtrmm_kernel_L1_END           // if zero, no need to run L1_M1

.Ldtrmm_kernel_L1_M1_20:

	INIT1x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
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

