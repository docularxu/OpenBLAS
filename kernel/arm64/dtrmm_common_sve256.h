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

/* Description:
 *   This is common macros for use by dtrmm on SVE256 hardware.
 *   Before include this file, the following assumptions/defines must be met.
 *
 * Assumptions/defines:
 *   p0         : all true predicate
 *   z12        : duplicated alpha
 *   pA         : address of next element of *a
 *   pB         : address of next element of *b
 *   pCRow0     : pointer to next position in *c Column 0
 *   pCRow1, pCRow2, pCRow3
 *              : defined but not initialized.
 *              : to be used as pointers to next position in *c Column 1/2/3
 */

/******************************************************************************/
/* L4                                                                         */
/******************************************************************************/

/* L4_M8 */
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
	ld1d  {z0.d}, p0/z, [pA]             // load a(  0:3, l )
	/* load b */
	ld1rqd  {z4.d}, p0/z, [pB]           // load b( l,0:1 )
	fmla  z16.d, z0.d, z4.d[0]           //

	ld1rqd  {z5.d}, p0/z, [pB, #16]      // load b( l,2:3 )
	fmla  z20.d, z0.d, z4.d[1]           //
	fmla  z24.d, z0.d, z5.d[0]           //

	ld1d  {z1.d}, p0/z, [pA, #1, MUL VL] // load a(  4:7, l )
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
	st1d	{z0.d}, p0, [pCRow0]            // store column 0
	st1d	{z1.d}, p0, [pCRow0, #1, MUL VL]

	add	pCRow1, pCRow0, LDC

	fmul	z4.d, z20.d, z12.d              // scale by alpha
	fmul	z5.d, z21.d, z12.d
	st1d	{z4.d}, p0, [pCRow1]            // store column 1
	st1d	{z5.d}, p0, [pCRow1, #1, MUL VL]

	add	pCRow2, pCRow1, LDC

	fmul	z0.d, z24.d, z12.d              // scale by alpha
	fmul	z1.d, z25.d, z12.d
	st1d	{z0.d}, p0, [pCRow2]            // store column 2
	st1d	{z1.d}, p0, [pCRow2, #1, MUL VL]

	add	pCRow3, pCRow2, LDC

	fmul	z4.d, z28.d, z12.d              // scale by alpha
	fmul	z5.d, z29.d, z12.d
	st1d	{z4.d}, p0, [pCRow3]            // store column 3
	st1d 	{z5.d}, p0, [pCRow3, #1, MUL VL]

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
	ld1d  {z0.d}, p0/z, [pA]             // load a(  0:3, l )
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
	st1d	{z0.d}, p0, [pCRow0]            // store column 0

	add	pCRow1, pCRow0, LDC

	fmul	z4.d, z20.d, z12.d              // scale by alpha
	st1d	{z4.d}, p0, [pCRow1]            // store column 1

	add	pCRow2, pCRow1, LDC

	fmul	z0.d, z24.d, z12.d              // scale by alpha
	st1d	{z0.d}, p0, [pCRow2]            // store column 2

	add	pCRow3, pCRow2, LDC

	fmul	z4.d, z28.d, z12.d              // scale by alpha
	st1d	{z4.d}, p0, [pCRow3]            // store column 3

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
	ld1d  {z4.d}, p0/z, [pB]             // load b( l,0:3 )

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
 *   z0.d[0] : a( 0,l )                        // col. l
 *   z4      : b( l,0:3 )                      // row  l
 */
 .macro KERNEL1x4_SUB
	/* load a */
	ld1rd  {z0.d}, p0/z, [pA]              // load a( 0,0 )
	/* load b */
	ld1d   {z4.d}, p0/z, [pB]              // load b( l,0:3 )

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

	/* todo: this instruction is unnecessary
	 *       unlike previous SAVEs, SAVE1x4 is the last in L4,
	 *       pCRow0 will be reset to pC
	 */
	add	pCRow0, pCRow0, #8              // pCRow0 moved down
.endm

/******************************************************************************/
/* L2                                                                         */ 
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
	ld1d  {z0.d}, p0/z, [pA]             // load a(  0:3, l )
	/* load b */
	ld1rqd  {z4.d}, p0/z, [pB]           // load b( l,0:1 )
        /* adjust pB */
        add   pB, pB, #16                    // 16 = 2 * sizeof(double)

	fmla  z16.d, z0.d, z4.d[0]           //
	fmla  z20.d, z0.d, z4.d[1]           //

	ld1d  {z1.d}, p0/z, [pA, #1, MUL VL] // load a(  4:7, l )
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
	st1d	{z0.d}, p0, [pCRow0]            // store column 0
	st1d	{z1.d}, p0, [pCRow0, #1, MUL VL]

	add	pCRow1, pCRow0, LDC

	fmul	z4.d, z20.d, z12.d              // scale by alpha
	fmul	z5.d, z21.d, z12.d
	st1d	{z4.d}, p0, [pCRow1]            // store column 1
	st1d	{z5.d}, p0, [pCRow1, #1, MUL VL]

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
	ld1d  {z0.d}, p0/z, [pA]             // load a(  0:3, l )
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
	st1d	{z0.d}, p0, [pCRow0]            // store column 0

	add	pCRow1, pCRow0, LDC

	fmul	z4.d, z20.d, z12.d              // scale by alpha
	st1d	{z4.d}, p0, [pCRow1]            // store column 1

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

	ld1	{v0.d}[0], [pA]
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
	ld1d  {z0.d}, p0/z, [pA]             // load a(  0:3, l )
	/* load b */
	ld1rd 	{z4.d}, p0/z, [pB]           // load b( l,0 )
	ld1d  {z1.d}, p0/z, [pA, #1, MUL VL] // load a(  4:7, l )

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
	st1d	{z0.d}, p0, [pCRow0]            // store column 0
	st1d	{z1.d}, p0, [pCRow0, #1, MUL VL]

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
	ld1d  {z0.d}, p0/z, [pA]             // load a(  0:3, l )
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
	st1d	{z1.d}, p0, [pCRow0]            // store column 0

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
	ld1	{v8.d}[0], [pB]
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
	ld1	{v8.d}[0], [pB]
	add	pB , pB, #8

	ld1	{v0.d}[0], [pA]
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
	st1	{v8.d}[0], [pCRow0]

	add	pCRow0, pCRow0, #8              // 8 = 1 * sizeof(double)
.endm


