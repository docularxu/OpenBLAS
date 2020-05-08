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
	/* TODO: split below ldr and ld1rqd into the fmla's. */
	/* load a */
	ldr  z0, [pA]                        // load a(  0:3, l )
	ldr  z1, [pA, #1, MUL VL]            // load a(  4:7, l )
	ldr  z2, [pA, #2, MUL VL]            // load a(  8:11,l )

	/* load b */
	ld1rqd  {z4.d}, p0/z, [pB]           // load b( l,0:1 )
	ld1rqd  {z5.d}, p0/z, [pB, #16]      // load b( l,2:3 )
	ld1rqd  {z6.d}, p0/z, [pB, #32]      // load b( l,4:5 )
	ld1rqd  {z7.d}, p0/z, [pB, #48]      // load b( l,6:7 )

	/* accumulate */
	/* Iteration 1 */
	/* Input:
	 *        z0, z1, z2    : a( 0:11,l )
	 *        z4, z5, z6, z7: b( l,0:7  )
	 *        z3            : free
	 * Output:
	 *        z3, z1, z0      : a( 0:11,l )
	 *        z4, z5, z6, z7? : b( l,0:7  )
	 *        z2              : free
	 */
	fmla  z8.d,  z0.d, z4.d[0]           //
	fmla  z9.d,  z1.d, z4.d[0]           //
	fmla  z10.d, z2.d, z4.d[0]           //
	ldr  z3, [pA, #3, MUL VL]            // load a(  0:3, l )

	fmla  z11.d, z0.d, z4.d[1]           //
	fmla  z12.d, z1.d, z4.d[1]           //
	fmla  z13.d, z2.d, z4.d[1]           //
	fmla  z14.d, z0.d, z5.d[0]           //
	ld1rqd  {z4.d}, p0/z, [pB, #64]      // load b( l+1,0:1 )

	fmla  z17.d, z0.d, z5.d[1]           //
	fmla  z20.d, z0.d, z6.d[0]           //
	fmla  z23.d, z0.d, z6.d[1]           //
	fmla  z26.d, z0.d, z7.d[0]           //
	fmla  z29.d, z0.d, z7.d[1]           //
	fmla  z16.d, z2.d, z5.d[0]           //
	ldr  z0, [pA, #5, MUL VL]            // load a(  8:11,l )

	fmla  z15.d, z1.d, z5.d[0]           //
	fmla  z18.d, z1.d, z5.d[1]           //
	fmla  z19.d, z2.d, z5.d[1]           //
	fmla  z21.d, z1.d, z6.d[0]           //
	ld1rqd  {z5.d}, p0/z, [pB, #80]      // load b( l+1,2:3 )

	fmla  z22.d, z2.d, z6.d[0]           //
	fmla  z24.d, z1.d, z6.d[1]           //
	fmla  z25.d, z2.d, z6.d[1]           //
	fmla  z27.d, z1.d, z7.d[0]           //
	ld1rqd  {z6.d}, p0/z, [pB, #96]      // load b( l+1,4:5 )

	fmla  z28.d, z2.d, z7.d[0]           //
	fmla  z30.d, z1.d, z7.d[1]           //
	fmla  z31.d, z2.d, z7.d[1]           //
	ldr  z1, [pA, #4, MUL VL]            // load a(  4:7, l )

	/* TODO: move to next iteration */
	// ld1rqd  {z7.d}, p0/z, [pB, #112]     // load b( l+1,6:7 )

	/* Iteration 2 */
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

	fmla  z11.d, z3.d, z4.d[1]           //
	fmla  z12.d, z1.d, z4.d[1]           //
	fmla  z13.d, z0.d, z4.d[1]           //
	fmla  z14.d, z3.d, z5.d[0]           //
	ld1rqd  {z4.d}, p0/z, [pB, #128]     // load b( l+2,0:1 )

	fmla  z15.d, z1.d, z5.d[0]           //
	fmla  z16.d, z0.d, z5.d[0]           //
	fmla  z17.d, z3.d, z5.d[1]           //
	ld1rqd  {z7.d}, p0/z, [pB, #112]     // load b( l+1,6:7 )

	fmla  z18.d, z1.d, z5.d[1]           //
	fmla  z19.d, z0.d, z5.d[1]           //
	fmla  z20.d, z3.d, z6.d[0]           //
	ld1rqd  {z5.d}, p0/z, [pB, #144]     // load b( l+2,2:3 )

	fmla  z21.d, z1.d, z6.d[0]           //
	fmla  z22.d, z0.d, z6.d[0]           //
	fmla  z23.d, z3.d, z6.d[1]           //
	fmla  z24.d, z1.d, z6.d[1]           //
	fmla  z25.d, z0.d, z6.d[1]           //

	ld1rqd  {z6.d}, p0/z, [pB, #160]     // load b( l+2,4:5 )

	fmla  z28.d, z0.d, z7.d[0]           //
	fmla  z31.d, z0.d, z7.d[1]           //
	fmla  z26.d, z3.d, z7.d[0]           //
	ldr  z0, [pA, #6, MUL VL]            // load a(  0:3, l+2 )

	fmla  z30.d, z1.d, z7.d[1]           //
	fmla  z27.d, z1.d, z7.d[0]           //
	fmla  z29.d, z3.d, z7.d[1]           //
	ldr  z1, [pA, #7, MUL VL]            // load a(  4:7, l+2 )

	/* todo: move to next iteration */
	// ld1rqd  {z7.d}, p0/z, [pB, #176]     // load b( l+2,6:7 )

	/* Iteration 3 */
.endm

/* Input:
 *   pA: address of loading from *a
 *   pB: address of loading from *b
 *   p0.d: all true predicate register
 * Vector registers:
 *   z8 ~ z31  : c
 *   z0 ~ z3   : a
 *   z4, z5, z6: b
 *   z7        : free, used in loop unrolling
 * Output:
 *   z8 ~ z31          : updated c
 *   pA: pointing to next elements of *a
 *   pB: pointing to next elements of *b
 */
.macro KERNEL16x6_I
	/* load a */
	ldr  z0, [pA]                        // load a(  0:3, l )
	ldr  z1, [pA, #1, MUL VL]            // load a(  4:7, l )
	ldr  z2, [pA, #2, MUL VL]            // load a(  8:11,l )
	ldr  z3, [pA, #3, MUL VL]            // load a( 12:15,l )

	/* load b */
	ld1rqd  {z4.d}, p0/z, [pB]           // load b( l,0:1 )
	ld1rqd  {z5.d}, p0/z, [pB, #16]      // load b( l,2:3 )
	ld1rqd  {z6.d}, p0/z, [pB, #32]      // load b( l,4:5 )

	/* accumulate */
	/* Iteration 1 */
	/* Input:
	 *        z0, z1, z2, z3: a( 0:15,l )
	 *        z4, z5, z6    : b( l,0:5  )
	 *        z7            : free
	 * Output:
	 *        z7, z1, z2, z0 : a( 0:15,l+1 )
	 *        z4, z5, z6?    : b( l+1,0:5  )
	 *        z3             : free
	 *        z6             : to be loaded in next iteration
	 */
	fmla  z8.d,  z0.d, z4.d[0]           //
	fmla  z9.d,  z1.d, z4.d[0]           //
	fmla  z10.d, z2.d, z4.d[0]           //
	ldr   z7, [pA, #4, MUL VL]           // load a(  0:3, l+1 )

	fmla  z11.d, z3.d, z4.d[0]           //
	fmla  z12.d, z0.d, z4.d[1]           //
	fmla  z13.d, z1.d, z4.d[1]           //
	fmla  z14.d, z2.d, z4.d[1]           //
	fmla  z15.d, z3.d, z4.d[1]           //
	fmla  z17.d, z1.d, z5.d[0]           //
	ld1rqd  {z4.d}, p0/z, [pB, #48]      // load b( l+1,0:1 )

	fmla  z21.d, z1.d, z5.d[1]           //
	fmla  z25.d, z1.d, z6.d[0]           //
	fmla  z29.d, z1.d, z6.d[1]           //
	fmla  z16.d, z0.d, z5.d[0]           //
	ldr   z1, [pA, #5, MUL VL]           // load a(  4:7, l+1 )

	fmla  z20.d, z0.d, z5.d[1]           //
	fmla  z18.d, z2.d, z5.d[0]           //
	fmla  z22.d, z2.d, z5.d[1]           //
	fmla  z30.d, z2.d, z6.d[1]           //
	fmla  z24.d, z0.d, z6.d[0]           //
	ldr   z2, [pA, #6, MUL VL]           // load a(  8:11, l+1 )

	fmla  z23.d, z3.d, z5.d[1]           //
	fmla  z19.d, z3.d, z5.d[0]           //
	fmla  z26.d, z2.d, z6.d[0]           //
	ld1rqd  {z5.d}, p0/z, [pB, #64]      // load b( l+1,2:3 )

	fmla  z28.d, z0.d, z6.d[1]           //
	fmla  z27.d, z3.d, z6.d[0]           //
	fmla  z31.d, z3.d, z6.d[1]           //
	ldr   z0, [pA, #7, MUL VL]           // load a( 12:15, l+1 )
	/* z6 not loaded in this iteration */

	/* Iteration 2 */
	/* Input:
	 *        z7, z1, z2, z0 : a( 0:15,l )
	 *        z4, z5, z6?    : b( l,0:5  )
	 *        z3             : free
	 *        z6             : to be loaded in this iteration
	 * Output:
	 *        z0, z1, z2, z3 : a( 0:15,l+1 )
	 *        z4, z5, z6?    : b( l+1,0:5  )
	 *        z7             : free
	 *        z6             : to be loaded in next iteration
	 */
	fmla  z8.d,  z7.d, z4.d[0]           //
	fmla  z9.d,  z1.d, z4.d[0]           //
	fmla  z10.d, z2.d, z4.d[0]           //
	ldr   z3, [pA, #11, MUL VL]          // load a( 12:15, l+2 )

	fmla  z11.d, z0.d, z4.d[0]           //
	fmla  z12.d, z7.d, z4.d[1]           //
	fmla  z13.d, z1.d, z4.d[1]           //
	fmla  z14.d, z2.d, z4.d[1]           //
	fmla  z15.d, z0.d, z4.d[1]           //
	fmla  z16.d, z7.d, z5.d[0]           //
	ld1rqd  {z4.d}, p0/z, [pB, #96]      // load b( l+1,0:1 )

	fmla  z17.d, z1.d, z5.d[0]           //
	fmla  z18.d, z2.d, z5.d[0]           //
	fmla  z19.d, z0.d, z5.d[0]           //
	fmla  z20.d, z7.d, z5.d[1]           //
	ld1rqd  {z6.d}, p0/z, [pB, #80]      // load b( l+1,4:5 )

	fmla  z21.d, z1.d, z5.d[1]           //
	fmla  z22.d, z2.d, z5.d[1]           //
	fmla  z23.d, z0.d, z5.d[1]           //
	fmla  z24.d, z7.d, z6.d[0]           //
	ld1rqd  {z5.d}, p0/z, [pB, #112]     // load b( l,2:3 )

	fmla  z25.d, z1.d, z6.d[0]           //
	fmla  z29.d, z1.d, z6.d[1]           //
	fmla  z30.d, z2.d, z6.d[1]           //
	ldr   z1, [pA, #9, MUL VL]           // load a(  4:7, l+2 )


	fmla  z26.d, z2.d, z6.d[0]           //
	fmla  z27.d, z0.d, z6.d[0]           //
	ldr   z2, [pA, #10, MUL VL]          // load a(  8:11, l+2 )

	fmla  z31.d, z0.d, z6.d[1]           //
	fmla  z28.d, z7.d, z6.d[1]           //
	ldr   z0, [pA, #8, MUL VL]           // load a(  0:3, l+2 )

	/* at the end of Iteration 2,
	 *   z0, z1, z2, z3: loaded for *a next column.
	 *   z4, z5        : loaded for *b next row
	 *   z6            : not loaded yet.
	 * Need to do: loading z6 for *b next row
	 */
	/* move to next interation */
	// ld1rqd  {z6.d}, p0/z, [pB, #128]     // load b( l+1,4:5 )


	// prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]

	/* Iteration 3 */
.endm
