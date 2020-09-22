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

#define	N	x0	/* vector length */
#define	X	x1	/* X vector address */
#define	INC_X	x2	/* X stride */
#define	Y	x3	/* Y vector address */
#define	INC_Y	x4	/* Y stride */
#define I	x5	/* loop variable */

/*******************************************************************************
* Macro definitions
*******************************************************************************/

#if !defined(DOUBLE)
#if !defined(DSDOT)
#define REG0	wzr
#define DOTF	s0
#else // DSDOT
#define REG0	xzr
#define DOTF	d0
#endif
#define DOTI	s1
#define TMPX	s2
#define LD1VX	{v2.s}[0]
#define TMPY	s3
#define LD1VY	{v3.s}[0]
#define TMPVY	v3.s[0]
#define SZ	4
#else
#define REG0	xzr
#define DOTF	d0
#define DOTI	d1
#define TMPX	d2
#define LD1VX	{v2.d}[0]
#define TMPY	d3
#define LD1VY	{v3.d}[0]
#define TMPVY	v3.d[0]
#define SZ	8
#endif

/******************************************************************************/

.macro KERNEL_F1
	ldr	TMPX, [X], #SZ
	ldr	TMPY, [Y], #SZ
#if !defined(DSDOT)
	fmadd	DOTF, TMPX, TMPY, DOTF
#else // DSDOT
	fcvt	d3, TMPY
	fcvt	d2, TMPX
	fmul	d2, d2, d3
	fadd	DOTF, DOTF, d2
#endif
.endm

/* SVE-256 version
 * Pre-condition:
 *  - p0: all true
 */
.macro KERNEL_SVE256_F4
#if !defined(DOUBLE) && defined(DSDOT)
	// load 4 singles
	ld1w	{z2.d}, p0/z, [X]
	ld1w	{z4.d}, p0/z, [Y]
	// convert them to doubles
	fcvt	z2.d, p0/m, z2.s
	fcvt	z4.d, p0/m, z4.s
	// multiply and add up
	fmla	z0.d, p0/m, z2.d, z4.d
	// address update
	add	X, X, #16
	add	Y, Y, #16
#elif !defined(DOUBLE) && !defined(DSDOT)
	// TODO: update to SVE code. However, at minimum, SVE256 can load 8 singles.
	//       Need to update to 8 singles. Before that, I keep NEON version.
	// load 4 singles
	ld1	{v2.4s}, [X], #16
	ld1	{v3.4s}, [Y], #16
	// multiply and add up
	fmla	v0.4s, v2.4s, v3.4s
#else //DOUBLE
	ld1d	{z2.d}, p0/z, [X]		// load 4 doubles
	ld1d	{z4.d}, p0/z, [Y]		// load 4 doubles

	fmla	z0.d, p0/m, z2.d, z4.d	// multiply and add up
	add	X, X, #32		// address updating
	add	Y, Y, #32		// address updating
#endif
	PRFM	PLDL1KEEP, [X, #1024]
	PRFM	PLDL1KEEP, [Y, #1024]
.endm

.macro KERNEL_SVE256_F4_FINALIZE
#if !defined(DOUBLE) && !defined(DSDOT)
	ext	v1.16b, v0.16b, v0.16b, #8
	fadd	v0.2s, v0.2s, v1.2s
	faddp	DOTF, v0.2s
#else //DOUBLE || DSDOT
	faddv	DOTF, p0, z0.d		// double precsion: floating-point add recursive reduction to scalar.
#endif
.endm

.macro INIT_S
#if !defined(DOUBLE)
	lsl	INC_X, INC_X, #2
	lsl	INC_Y, INC_Y, #2
#else
	lsl	INC_X, INC_X, #3
	lsl	INC_Y, INC_Y, #3
#endif
.endm

.macro KERNEL_S1
	ld1	LD1VX, [X], INC_X
	ld1	LD1VY, [Y], INC_Y
#if !defined(DSDOT)
	fmadd	DOTF, TMPX, TMPY, DOTF
#else // DSDOT
	fcvt	d3, TMPY
	fcvt	d2, TMPX
	fmul	d2, d2, d3
	fadd	DOTF, DOTF, d2
#endif
.endm

/*******************************************************************************
* End of macro definitions
*******************************************************************************/

	PROLOGUE

	fmov	DOTF, REG0
#if defined(DOUBLE)
	fmov	d6, DOTF
#endif

	cmp	N, xzr
	ble	.Ldot_kernel_L999

	cmp	INC_X, #1
	bne	.Ldot_kernel_S_BEGIN
	cmp	INC_Y, #1
	bne	.Ldot_kernel_S_BEGIN

.Ldot_kernel_F_BEGIN:

	// initialize p0 to all true
	// initialize z0 to all #0.0
#if !defined(DOUBLE)
	ptrue	p0.s, all
	dup 	z0.s, #0
#else
	ptrue	p0.d, all
	dup	z0.d, #0
#endif

	asr	I, N, #2
	cmp	I, xzr
	beq	.Ldot_kernel_F1

.Ldot_kernel_F4:

	KERNEL_SVE256_F4

	subs	I, I, #1
	bne	.Ldot_kernel_F4

	KERNEL_SVE256_F4_FINALIZE

.Ldot_kernel_F1:

	ands	I, N, #3
	ble	.Ldot_kernel_L999

.Ldot_kernel_F10:

	KERNEL_F1

	subs    I, I, #1
        bne     .Ldot_kernel_F10

	ret

.Ldot_kernel_S_BEGIN:

	INIT_S

	asr	I, N, #2
	cmp	I, xzr
	ble	.Ldot_kernel_S1

.Ldot_kernel_S4:

	KERNEL_S1
	KERNEL_S1
	KERNEL_S1
	KERNEL_S1

	subs	I, I, #1
	bne	.Ldot_kernel_S4

.Ldot_kernel_S1:

	ands	I, N, #3
	ble	.Ldot_kernel_L999

.Ldot_kernel_S10:

	KERNEL_S1

	subs    I, I, #1
        bne     .Ldot_kernel_S10

.Ldot_kernel_L999:

	ret

	EPILOGUE
