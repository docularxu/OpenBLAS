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

/*                   X0          X1          X2          s0        X3        x4       x5           x6               x7 */
/*int CNAME(BLASLONG bm,BLASLONG bn,BLASLONG bk,FLOAT alpha,FLOAT* ba,FLOAT* bb,FLOAT* C,BLASLONG ldc, BLASLONG offset) */

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
#define alpha		w17
#define pCRow4		x18
#define pCRow5		x19
#define pCRow6		x20
#define pCRow7		x21
#define pCRow8		x22
#define pCRow9		x23
#define pCRow10		x24
#define pCRow11		x25
#define temp		x26
#define tempOffset	x27
#define tempK		x28
// 29 frame
// 30 link
// 31 sp

/* z31 is used as alpha in all macros except 16x12 */
#define alpha0		z31.s
#define alphaV0		z31.s[0]

/* z3 is used as alpha in macro 16x12 */
#define alpha0_t	z3.s
#define alphaV0_t	z3.s[0]

/* predicate registers */
// p0: all true
// p1: VL1
// p2: VL2
// p4: VL4

#define A_PRE_SIZE	2560
#define B_PRE_SIZE	224
#define C_PRE_SIZE	160

//v00 ALPHA -> pA0_00, pA0_01, pA0_02, pA0_03
//v01 pA0_04, pA0_05, pA0_06, pA0_07
//v02 pA0_08, pA0_09, pA0_10, pA0_11
//v03 pA0_12, pA0_13, pA0_14, pA0_15
//v04 pA1_00, pA1_01, pA1_02, pA1_03
//v05 pA1_04, pA1_05, pA1_06, pA1_07
//v06 pA1_08, pA1_09, pA1_10, pA1_11
//v07 pA1_12, pA1_13, pA1_14, pA1_15
//v08 must save pB00
//v09 must save pB01
//v10 must save pB02
//v11 must save pB03
//v12 must save pB10
//v13 must save pB11
//v14 must save pB12
//v15 must save pB13
//v16 must save C00, C01, C02, C03
//v17 must save C04, C05, C06, C07
//v18 C08, C09, C10, C11
//v19 C12, C13, C14, C15
//v20 C16, C17, C18, C19
//v21 C20, C21, C22, C23
//v22 C24, C25, C26, C27
//v23 C28, C29, C30, C31
//v24 C32, C33, C34, C35
//v25 C36, C37, C38, C39
//v26 C40, C41, C42, C43
//v27 C44, C45, C46, C47
//v28 C48, C49, C50, C51
//v29 C52, C53, C54, C55
//v30 C56, C57, C58, C59
//v31 C60, C61, C62, C63

/*******************************************************************************
* Macro definitions
*******************************************************************************/

/* Register Allocation:
 *  z0,z1	: A
 *  z2,z3	: A (alternative, in _I, _M1)
 * 
 *  z4,z5,z6	: B
 *  z4,z5,z7	: B (alternative, in _M2, _E)
 *
 *  z8,z9...z31	: C
 *
 *  p0		: all enabled
 */

.macro INIT16x12
        mov     z8.s,  #0
        mov     z9.s,  #0
        mov     z10.s, #0
        mov     z11.s, #0

        mov     z12.s, #0
        mov     z13.s, #0
        mov     z14.s, #0
        mov     z15.s, #0

        mov     z16.s, #0
        mov     z17.s, #0
        mov     z18.s, #0
        mov     z19.s, #0

/* TODO: we can use .sup z20.s, z16.s[0]' to initilize the other vectors.
 *       Question is: will that be more efficient than using 'mov z20.s, #0'?
 */
        mov     z20.s, #0
        mov     z21.s, #0
        mov     z22.s, #0
        mov     z23.s, #0

        mov     z24.s, #0
        mov     z25.s, #0
        mov     z26.s, #0
        mov     z27.s, #0

        mov     z28.s, #0
        mov     z29.s, #0
        mov     z30.s, #0
        mov     z31.s, #0
.endm

.macro KERNEL16x12_SUB
        ld1w    {z0.s}, p0/z, [pA]
        ld1rqw  {z4.s}, p0/z, [pB]

        fmla    z8.s,  z0.s,  z4.s[0]
        fmla    z10.s,  z0.s,  z4.s[1]
        fmla    z12.s,  z0.s,  z4.s[2]
        fmla    z14.s,  z0.s,  z4.s[3]

        ld1w    {z1.s}, p0/z, [pA, #1, MUL VL]
        fmla    z9.s,  z1.s,  z4.s[0]
        fmla    z11.s,  z1.s,  z4.s[1]
        fmla    z13.s,  z1.s,  z4.s[2]
        fmla    z15.s,  z1.s,  z4.s[3]

        ld1rqw  {z5.s}, p0/z, [pB, #16]
        fmla    z16.s,  z0.s,  z5.s[0]
        fmla    z17.s,  z1.s,  z5.s[0]
        fmla    z18.s,  z0.s,  z5.s[1]
        fmla    z19.s,  z1.s,  z5.s[1]
        incb    pA, all, MUL #2

        ld1rqw  {z6.s}, p0/z, [pB, #32]
        fmla    z20.s,  z0.s,  z5.s[2]
        fmla    z21.s,  z1.s,  z5.s[2]
        fmla    z22.s,  z0.s,  z5.s[3]
        fmla    z23.s,  z1.s,  z5.s[3]
        add     pB, pB, #48

        fmla    z24.s,  z0.s,  z6.s[0]
        fmla    z25.s,  z1.s,  z6.s[0]
        fmla    z26.s,  z0.s,  z6.s[1]
        fmla    z27.s,  z1.s,  z6.s[1]

        fmla    z28.s,  z0.s,  z6.s[2]
        fmla    z29.s,  z1.s,  z6.s[2]
        fmla    z30.s,  z0.s,  z6.s[3]
        fmla    z31.s,  z1.s,  z6.s[3]
.endm

/* macros KERNEL16x12_I/_M2/_M1/_E are copied from
 *  sgemm_kernel_sve256_16x12.S
 */
.macro KERNEL16x12_I
        ld1w    {z0.s}, p0/z, [pA]
        ld1w    {z1.s}, p0/z, [pA, #1, MUL VL]

        ld1rqw  {z4.s}, p0/z, [pB]
        ld1rqw  {z5.s}, p0/z, [pB, #16]
        ld1rqw  {z6.s}, p0/z, [pB, #32]

        fmul    z8.s,  z0.s,  z4.s[0]
        fmul    z10.s,  z0.s,  z4.s[1]
        fmul    z12.s,  z0.s,  z4.s[2]
        fmul    z14.s,  z0.s,  z4.s[3]

        ld1w    {z2.s}, p0/z, [pA, #2, MUL VL]

        fmul    z9.s,  z1.s,  z4.s[0]
        fmul    z11.s,  z1.s,  z4.s[1]
        fmul    z13.s,  z1.s,  z4.s[2]
        fmul    z15.s,  z1.s,  z4.s[3]

        ld1rqw  {z4.s}, p0/z, [pB, #48]

        fmul    z16.s,  z0.s,  z5.s[0]
        fmul    z18.s,  z0.s,  z5.s[1]
        fmul    z17.s,  z1.s,  z5.s[0]
        fmul    z19.s,  z1.s,  z5.s[1]

        ld1w    {z3.s}, p0/z, [pA, #3, MUL VL]

        fmul    z20.s,  z0.s,  z5.s[2]
        fmul    z22.s,  z0.s,  z5.s[3]
        fmul    z21.s,  z1.s,  z5.s[2]
        fmul    z23.s,  z1.s,  z5.s[3]

        ld1rqw  {z5.s}, p0/z, [pB, #64]
        fmul    z24.s,  z0.s,  z6.s[0]
        fmul    z26.s,  z0.s,  z6.s[1]
        fmul    z25.s,  z1.s,  z6.s[0]
        fmul    z27.s,  z1.s,  z6.s[1]

        ld1rqw  {z7.s}, p0/z, [pB, #80]

        fmul    z28.s,  z0.s,  z6.s[2]
        fmul    z30.s,  z0.s,  z6.s[3]
        fmul    z29.s,  z1.s,  z6.s[2]
        fmul    z31.s,  z1.s,  z6.s[3]

        incb    pA, all, MUL #4
        add     pB, pB, #96
.endm

.macro KERNEL16x12_M1
        ld1w    {z2.s}, p0/z, [pA]
        fmla    z8.s, z0.s,  z4.s[0]
        fmla    z9.s, z1.s,  z4.s[0]
        fmla    z10.s, z0.s,  z4.s[1]
        fmla    z11.s, z1.s,  z4.s[1]

        ld1w    {z3.s}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        fmla    z12.s,  z0.s,  z4.s[2]
        fmla    z13.s,  z1.s,  z4.s[2]
        fmla    z14.s,  z0.s,  z4.s[3]
        fmla    z15.s,  z1.s,  z4.s[3]

        ld1rqw  {z4.s}, p0/z, [pB]

        fmla    z16.s,  z0.s,  z5.s[0]
        fmla    z17.s,  z1.s,  z5.s[0]
        fmla    z18.s,  z0.s,  z5.s[1]
        fmla    z19.s,  z1.s,  z5.s[1]
        fmla    z20.s,  z0.s,  z5.s[2]
        fmla    z21.s,  z1.s,  z5.s[2]
        fmla    z22.s,  z0.s,  z5.s[3]
        fmla    z23.s,  z1.s,  z5.s[3]

        ld1rqw  {z5.s}, p0/z, [pB, #16]

        fmla    z24.s,  z0.s,  z6.s[0]
        fmla    z25.s,  z1.s,  z6.s[0]
        fmla    z26.s,  z0.s,  z6.s[1]
        fmla    z27.s,  z1.s,  z6.s[1]

        ld1rqw  {z7.s}, p0/z, [pB, #32]
        add     pB, pB, #48

        fmla    z28.s,  z0.s,  z6.s[2]
        fmla    z29.s,  z1.s,  z6.s[2]
        fmla    z30.s,  z0.s,  z6.s[3]
        fmla    z31.s,  z1.s,  z6.s[3]
.endm

.macro KERNEL16x12_M2
        ld1w    {z0.s}, p0/z, [pA]
        fmla    z8.s,  z2.s,  z4.s[0]
        fmla    z9.s,  z3.s,  z4.s[0]
        fmla    z10.s,  z2.s,  z4.s[1]
        fmla    z11.s,  z3.s,  z4.s[1]

        ld1w    {z1.s}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        fmla    z12.s,  z2.s,  z4.s[2]
        fmla    z13.s,  z3.s,  z4.s[2]
        fmla    z14.s,  z2.s,  z4.s[3]
        fmla    z15.s,  z3.s,  z4.s[3]

        ld1rqw  {z4.s}, p0/z, [pB]

        fmla    z16.s,  z2.s,  z5.s[0]
        fmla    z17.s,  z3.s,  z5.s[0]
        fmla    z18.s,  z2.s,  z5.s[1]
        fmla    z19.s,  z3.s,  z5.s[1]
        fmla    z20.s,  z2.s,  z5.s[2]
        fmla    z21.s,  z3.s,  z5.s[2]
        fmla    z22.s,  z2.s,  z5.s[3]
        fmla    z23.s,  z3.s,  z5.s[3]

        ld1rqw  {z5.s}, p0/z, [pB, #16]

        fmla    z24.s,  z2.s,  z7.s[0]
        fmla    z25.s,  z3.s,  z7.s[0]
        fmla    z26.s,  z2.s,  z7.s[1]
        fmla    z27.s,  z3.s,  z7.s[1]

        ld1rqw  {z6.s}, p0/z, [pB, #32]
        add     pB, pB, #48

        fmla    z28.s,  z2.s,  z7.s[2]
        fmla    z29.s,  z3.s,  z7.s[2]
        fmla    z30.s,  z2.s,  z7.s[3]
        fmla    z31.s,  z3.s,  z7.s[3]
.endm

.macro KERNEL16x12_E
        fmla    z8.s,  z2.s,  z4.s[0]
        fmla    z9.s,  z3.s,  z4.s[0]
        fmla    z10.s,  z2.s,  z4.s[1]
        fmla    z11.s,  z3.s,  z4.s[1]

        fmla    z12.s, z2.s, z4.s[2]
        fmla    z13.s, z3.s, z4.s[2]
        fmla    z14.s, z2.s, z4.s[3]
        fmla    z15.s, z3.s, z4.s[3]

        fmla    z16.s, z2.s, z5.s[0]
        fmla    z17.s, z3.s, z5.s[0]
        fmla    z18.s, z2.s, z5.s[1]
        fmla    z19.s, z3.s, z5.s[1]

        fmla    z20.s, z2.s, z5.s[2]
        fmla    z21.s, z3.s, z5.s[2]
        fmla    z22.s, z2.s, z5.s[3]
        fmla    z23.s, z3.s, z5.s[3]

        fmla    z24.s, z2.s, z7.s[0]
        fmla    z25.s, z3.s, z7.s[0]
        fmla    z26.s, z2.s, z7.s[1]
        fmla    z27.s, z3.s, z7.s[1]

        fmla    z28.s, z2.s, z7.s[2]
        fmla    z29.s, z3.s, z7.s[2]
        fmla    z30.s, z2.s, z7.s[3]
        fmla    z31.s, z3.s, z7.s[3]
.endm

/* Register allocation:
 *  z0,z1,z4,z5,z6,z7	: temporary to hold fmla results
 *  z3 (alpha0_t)	: alpha
 *  z8,z9,........z31	: C
 */
.macro SAVE16x12
	dup	alpha0_t, alpha			// z3 (alpha0_t) as alpha in L12_M16
						// z31 (alpha0) as alpha elsewhere

/* TODO: comparing 
 *      'fmul z0.s, z8.s, z3.s'
 *      'fmul z0.s, z8.s, z3.s[0]'
 *   which one is more efficient?
 */
        fmul    z0.s, z8.s, alphaV0_t		// scale by alpha
        fmul    z1.s, z9.s, alphaV0_t
        st1w    {z0.s}, p0, [pCRow0]		// store column 0
        st1w    {z1.s}, p0, [pCRow0, #1, MUL VL]
        incb    pCRow0, all, MUL #2
/* TODO: other alternative of 'incb' are:
 *        'addvl pCRow0, pCRow0, #2'
 *        'add   pCRow0, pCRow0, #128', but this is not VL agnostic.
 */

        fmul    z4.s, z10.s, alphaV0_t
        fmul    z5.s, z11.s, alphaV0_t
        st1w    {z4.s}, p0, [pCRow1]
        st1w    {z5.s}, p0, [pCRow1, #1, MUL VL]
        incb    pCRow1, all, MUL #2

        fmul    z6.s, z12.s, alphaV0_t
        fmul    z7.s, z13.s, alphaV0_t
        st1w    {z6.s}, p0, [pCRow2]
        st1w    {z7.s}, p0, [pCRow2, #1, MUL VL]
        incb    pCRow2, all, MUL #2

        fmul    z0.s, z14.s, alphaV0_t
        fmul    z1.s, z15.s, alphaV0_t
        st1w    {z0.s}, p0, [pCRow3]
        st1w    {z1.s}, p0, [pCRow3, #1, MUL VL]
        incb    pCRow3, all, MUL #2

        fmul    z4.s, z16.s, alphaV0_t
        fmul    z5.s, z17.s, alphaV0_t
        st1w    {z4.s}, p0, [pCRow4]
        st1w    {z5.s}, p0, [pCRow4, #1, MUL VL]
        incb    pCRow4, all, MUL #2

        fmul    z6.s, z18.s, alphaV0_t
        fmul    z7.s, z19.s, alphaV0_t
        st1w    {z6.s}, p0, [pCRow5]
        st1w    {z7.s}, p0, [pCRow5, #1, MUL VL]
        incb    pCRow5, all, MUL #2

        fmul    z0.s, z20.s, alphaV0_t
        fmul    z1.s, z21.s, alphaV0_t
        st1w    {z0.s}, p0, [pCRow6]
        st1w    {z1.s}, p0, [pCRow6, #1, MUL VL]
        incb    pCRow6, all, MUL #2

        fmul    z4.s, z22.s, alphaV0_t
        fmul    z5.s, z23.s, alphaV0_t
        st1w    {z4.s}, p0, [pCRow7]
        st1w    {z5.s}, p0, [pCRow7, #1, MUL VL]
        incb    pCRow7, all, MUL #2

        fmul    z6.s, z24.s, alphaV0_t
        fmul    z7.s, z25.s, alphaV0_t
        st1w    {z6.s}, p0, [pCRow8]
        st1w    {z7.s}, p0, [pCRow8, #1, MUL VL]
        incb    pCRow8, all, MUL #2

        fmul    z0.s, z26.s, alphaV0_t
        fmul    z1.s, z27.s, alphaV0_t
        st1w    {z0.s}, p0, [pCRow9]
        st1w    {z1.s}, p0, [pCRow9, #1, MUL VL]
        incb    pCRow9, all, MUL #2

        fmul    z4.s, z28.s, alphaV0_t
        fmul    z5.s, z29.s, alphaV0_t
        st1w    {z4.s}, p0, [pCRow10]
        st1w    {z5.s}, p0, [pCRow10, #1, MUL VL]
        incb    pCRow10, all, MUL #2

        fmul    z6.s, z30.s, alphaV0_t
        fmul    z7.s, z31.s, alphaV0_t
        st1w    {z6.s}, p0, [pCRow11]
        st1w    {z7.s}, p0, [pCRow11, #1, MUL VL]
        incb    pCRow11, all, MUL #2
.endm

/******************************************************************************/
.macro INIT8x12
        mov     z8.s,  #0
        mov     z9.s,  #0
        mov     z10.s, #0
        mov     z11.s, #0

        mov     z12.s, #0
        mov     z13.s, #0
        mov     z14.s, #0
        mov     z15.s, #0

        mov     z16.s, #0
        mov     z17.s, #0
        mov     z18.s, #0
        mov     z19.s, #0
.endm

/* macros KERNEL8x12_I/_M2/_M1/_E/_SUB are copied from
 *  sgemm_kernel_sve256_16x12.S
 */

/* register allocation:
 *   z0			: A
 *   z1			: A (alternative, _M2, _E)
 *
 *   z2,z3,z4		: B
 *   z5,z6,z7		: B (alternative, _M2, _E)
 *
 *   z8,z9......z19	: C
 */
.macro KERNEL8x12_I
        ld1w    {z0.s}, p0/z, [pA]

        ld1rqw  {z2.s}, p0/z, [pB]
        ld1rqw  {z3.s}, p0/z, [pB, #16]
        ld1rqw  {z4.s}, p0/z, [pB, #32]

        fmul    z8.s,  z0.s,  z2.s[0]
        fmul    z9.s,  z0.s,  z2.s[1]

        ld1w    {z1.s}, p0/z, [pA, #1, MUL VL]

        fmul    z10.s,  z0.s,  z2.s[2]
        fmul    z11.s,  z0.s,  z2.s[3]

        ld1rqw  {z5.s}, p0/z, [pB, #48]
        fmul    z12.s,  z0.s,  z3.s[0]
        fmul    z13.s,  z0.s,  z3.s[1]

        ld1rqw  {z6.s}, p0/z, [pB, #64]
        fmul    z14.s,  z0.s,  z3.s[2]
        fmul    z15.s,  z0.s,  z3.s[3]

        ld1rqw  {z7.s}, p0/z, [pB, #80]
        fmul    z16.s,  z0.s,  z4.s[0]
        fmul    z17.s,  z0.s,  z4.s[1]
        fmul    z18.s,  z0.s,  z4.s[2]
        fmul    z19.s,  z0.s,  z4.s[3]

        incb    pA, all, MUL #2
        add     pB, pB, #96
.endm

.macro KERNEL8x12_M1
        ld1w    {z1.s}, p0/z, [pA]
        add     pA, pA, #32

        fmla    z8.s,  z0.s,  z2.s[0]
        fmla    z9.s,  z0.s,  z2.s[1]
        fmla    z10.s,  z0.s,  z2.s[2]
        fmla    z11.s,  z0.s,  z2.s[3]

        ld1rqw  {z5.s}, p0/z, [pB]

        fmla    z12.s,  z0.s,  z3.s[0]
        fmla    z13.s,  z0.s,  z3.s[1]
        fmla    z14.s,  z0.s,  z3.s[2]
        fmla    z15.s,  z0.s,  z3.s[3]

        ld1rqw  {z6.s}, p0/z, [pB, #16]

        fmla    z16.s,  z0.s,  z4.s[0]
        fmla    z17.s,  z0.s,  z4.s[1]
        fmla    z18.s,  z0.s,  z4.s[2]
        fmla    z19.s,  z0.s,  z4.s[3]

        ld1rqw  {z7.s}, p0/z, [pB, #32]
        add     pB, pB, #48
.endm

.macro KERNEL8x12_M2
        ld1w    {z0.s}, p0/z, [pA]
        add     pA, pA, #32

        fmla    z8.s,  z1.s,  z5.s[0]
        fmla    z9.s,  z1.s,  z5.s[1]
        fmla    z10.s,  z1.s,  z5.s[2]
        fmla    z11.s,  z1.s,  z5.s[3]

        ld1rqw  {z2.s}, p0/z, [pB]

        fmla    z12.s,  z1.s,  z6.s[0]
        fmla    z13.s,  z1.s,  z6.s[1]
        fmla    z14.s,  z1.s,  z6.s[2]
        fmla    z15.s,  z1.s,  z6.s[3]

        ld1rqw  {z3.s}, p0/z, [pB, #16]

        fmla    z16.s,  z1.s,  z7.s[0]
        fmla    z17.s,  z1.s,  z7.s[1]
        fmla    z18.s,  z1.s,  z7.s[2]
        fmla    z19.s,  z1.s,  z7.s[3]

        ld1rqw  {z4.s}, p0/z, [pB, #32]
        add     pB, pB, #48
.endm

.macro KERNEL8x12_E
        fmla    z8.s,  z1.s,  z5.s[0]
        fmla    z9.s,  z1.s,  z5.s[1]
        fmla    z10.s,  z1.s,  z5.s[2]
        fmla    z11.s,  z1.s,  z5.s[3]

        fmla    z12.s,  z1.s,  z6.s[0]
        fmla    z13.s,  z1.s,  z6.s[1]
        fmla    z14.s,  z1.s,  z6.s[2]
        fmla    z15.s,  z1.s,  z6.s[3]

        fmla    z16.s,  z1.s,  z7.s[0]
        fmla    z17.s,  z1.s,  z7.s[1]
        fmla    z18.s,  z1.s,  z7.s[2]
        fmla    z19.s,  z1.s,  z7.s[3]
.endm

.macro KERNEL8x12_SUB
        ld1w    {z0.s}, p0/z, [pA]
        add     pA, pA, #32

        ld1rqw  {z2.s}, p0/z, [pB]
        ld1rqw  {z3.s}, p0/z, [pB, #16]
        ld1rqw  {z4.s}, p0/z, [pB, #32]
        add     pB, pB, #48

        fmla    z8.s,  z0.s,  z2.s[0]
        fmla    z9.s,  z0.s,  z2.s[1]
        fmla    z10.s,  z0.s,  z2.s[2]
        fmla    z11.s,  z0.s,  z2.s[3]

        fmla    z12.s,  z0.s,  z3.s[0]
        fmla    z13.s,  z0.s,  z3.s[1]
        fmla    z14.s,  z0.s,  z3.s[2]
        fmla    z15.s,  z0.s,  z3.s[3]

        fmla    z16.s,  z0.s,  z4.s[0]
        fmla    z17.s,  z0.s,  z4.s[1]
        fmla    z18.s,  z0.s,  z4.s[2]
        fmla    z19.s,  z0.s,  z4.s[3]
.endm

.macro SAVE8x12
        fmul    z0.s, z8.s, alpha0
        st1w    {z0.s}, p0, [pCRow0]
        add     pCRow0, pCRow0, #32

        fmul    z1.s, z9.s, alpha0
        st1w    {z1.s}, p0, [pCRow1]
        add     pCRow1, pCRow1, #32

        fmul    z2.s, z10.s, alpha0
        st1w    {z2.s}, p0, [pCRow2]
        add     pCRow2, pCRow2, #32

        fmul    z4.s, z11.s, alpha0
        st1w    {z4.s}, p0, [pCRow3]
        add     pCRow3, pCRow3, #32

        fmul    z5.s, z12.s, alpha0
        st1w    {z5.s}, p0, [pCRow4]
        add     pCRow4, pCRow4, #32

        fmul    z6.s, z13.s, alpha0
        st1w    {z6.s}, p0, [pCRow5]
        add     pCRow5, pCRow5, #32

        fmul    z7.s, z14.s, alpha0
        st1w    {z7.s}, p0, [pCRow6]
        add     pCRow6, pCRow6, #32

        fmul    z0.s, z15.s, alpha0
        st1w    {z0.s}, p0, [pCRow7]
        add     pCRow7, pCRow7, #32
        
        fmul    z1.s, z16.s, alpha0
        st1w    {z1.s}, p0, [pCRow8]
        add     pCRow8, pCRow8, #32

        fmul    z2.s, z17.s, alpha0
        st1w    {z2.s}, p0, [pCRow9]
        add     pCRow9, pCRow9, #32

        fmul    z4.s, z18.s, alpha0
        st1w    {z4.s}, p0, [pCRow10]
        add     pCRow10, pCRow10, #32

        fmul    z5.s, z19.s, alpha0
        st1w    {z5.s}, p0, [pCRow11]
        add     pCRow11, pCRow11, #32
.endm

/******************************************************************************/
.macro INIT4x12
        mov     z8.s,  #0
        mov     z9.s,  #0
        mov     z10.s, #0
        mov     z11.s, #0

        mov     z12.s, #0
        mov     z13.s, #0
        mov     z14.s, #0
        mov     z15.s, #0

        mov     z16.s, #0
        mov     z17.s, #0
        mov     z18.s, #0
        mov     z19.s, #0
.endm

/* macros KERNEL4x12_I/_M2/_M1/_E/_SUB are copied from
 *  sgemm_kernel_sve256_16x12.S with modifications:
 *    - added predicated limit into fmul and fmla 
 */

/* register allocation:
 *   z0			: A
 *   z1			: A (alternative, _M2, _E)
 *
 *   z2,z3,z4		: B
 *   z5,z6,z7		: B (alternative, _M2, _E)
 *
 *   z8,z9......z19	: C
 *
 *   p4			: VL4
 */

/* TODO: the main problem of this implementation is:
 *    - It doesn't make full use of a 256-bit vector length. In each fmul/fmla
 *      only 4-singles are calculated, which on its full width SVE-256 can calculate
 *      8.
 * TODO: need to fix it. We are wasting the SVE calculation power.
 */
.macro KERNEL4x12_I
        ld1rqw  {z0.s}, p4/z, [pA]

        ld1rqw  {z2.s}, p4/z, [pB]
        ld1rqw  {z3.s}, p4/z, [pB, #16]
        ld1rqw  {z4.s}, p4/z, [pB, #32]

        fmul    z8.s,   z0.s,  z2.s[0]
        fmul    z9.s,   z0.s,  z2.s[1]

        ld1rqw  {z1.s}, p4/z, [pA, #16]

        fmul    z10.s,   z0.s,  z2.s[2]
        fmul    z11.s,   z0.s,  z2.s[3]

        ld1rqw  {z5.s}, p4/z, [pB, #48]
        fmul    z12.s,   z0.s,  z3.s[0]
        fmul    z13.s,   z0.s,  z3.s[1]

        ld1rqw  {z6.s}, p4/z, [pB, #64]
        fmul    z14.s,   z0.s,  z3.s[2]
        fmul    z15.s,   z0.s,  z3.s[3]

        ld1rqw  {z7.s}, p4/z, [pB, #80]
        fmul    z16.s,   z0.s,  z4.s[0]
        fmul    z17.s,   z0.s,  z4.s[1]
        fmul    z18.s,   z0.s,  z4.s[2]
        fmul    z19.s,   z0.s,  z4.s[3]

        add     pA, pA, #32
        add     pB, pB, #96
.endm

.macro KERNEL4x12_M1
        ld1rqw  {z1.s}, p4/z, [pA]
        fmla    z8.s,   z0.s,  z2.s[0]
        fmla    z9.s,   z0.s,  z2.s[1]
        fmla    z10.s,   z0.s,  z2.s[2]
        fmla    z11.s,   z0.s,  z2.s[3]

        ld1rqw  {z5.s}, p4/z, [pB]

        fmla    z12.s,   z0.s,  z3.s[0]
        fmla    z13.s,   z0.s,  z3.s[1]
        fmla    z14.s,   z0.s,  z3.s[2]
        fmla    z15.s,   z0.s,  z3.s[3]

        ld1rqw  {z6.s}, p4/z, [pB, #16]
        fmla    z16.s,   z0.s,  z4.s[0]
        fmla    z17.s,   z0.s,  z4.s[1]
        fmla    z18.s,   z0.s,  z4.s[2]
        fmla    z19.s,   z0.s,  z4.s[3]

        ld1rqw  {z7.s}, p4/z, [pB, #32]

        add     pA, pA, #16
        add     pB, pB, #48
.endm

.macro KERNEL4x12_M2
        ld1rqw {z0.s}, p4/z, [pA]

        fmla    z8.s,   z1.s,  z5.s[0]
        fmla    z9.s,   z1.s,  z5.s[1]
        fmla    z10.s,   z1.s,  z5.s[2]
        fmla    z11.s,   z1.s,  z5.s[3]

        ld1rqw  {z2.s}, p4/z, [pB]

        fmla    z12.s,   z1.s,  z6.s[0]
        fmla    z13.s,   z1.s,  z6.s[1]
        fmla    z14.s,   z1.s,  z6.s[2]
        fmla    z15.s,   z1.s,  z6.s[3]

        ld1rqw  {z3.s}, p4/z, [pB, #16]

        fmla    z16.s,   z1.s,  z7.s[0]
        fmla    z17.s,   z1.s,  z7.s[1]
        fmla    z18.s,   z1.s,  z7.s[2]
        fmla    z19.s,   z1.s,  z7.s[3]

        ld1rqw  {z4.s}, p4/z, [pB, #32]

        add     pA, pA, #16
        add     pB, pB, #48
.endm

.macro KERNEL4x12_E
        fmla    z8.s,   z1.s,  z5.s[0]
        fmla    z9.s,   z1.s,  z5.s[1]
        fmla    z10.s,   z1.s,  z5.s[2]
        fmla    z11.s,   z1.s,  z5.s[3]

        fmla    z12.s,   z1.s,  z6.s[0]
        fmla    z13.s,   z1.s,  z6.s[1]
        fmla    z14.s,   z1.s,  z6.s[2]
        fmla    z15.s,   z1.s,  z6.s[3]

        fmla    z16.s,   z1.s,  z7.s[0]
        fmla    z17.s,   z1.s,  z7.s[1]
        fmla    z18.s,   z1.s,  z7.s[2]
        fmla    z19.s,   z1.s,  z7.s[3]
.endm

.macro KERNEL4x12_SUB
        ld1rqw  {z0.s}, p4/z, [pA]
        ld1rqw  {z2.s}, p4/z, [pB]
        ld1rqw  {z3.s}, p4/z, [pB, #16]
        ld1rqw  {z4.s}, p4/z, [pB, #32]
       
        fmla    z8.s,   z0.s,  z2.s[0]
        fmla    z9.s,   z0.s,  z2.s[1]
        fmla    z10.s,   z0.s,  z2.s[2]
        fmla    z11.s,   z0.s,  z2.s[3]
    
        fmla    z12.s,   z0.s,  z3.s[0]
        fmla    z13.s,   z0.s,  z3.s[1]
        fmla    z14.s,   z0.s,  z3.s[2]
        fmla    z15.s,   z0.s,  z3.s[3]
    
        fmla    z16.s,   z0.s,  z4.s[0]
        fmla    z17.s,   z0.s,  z4.s[1]
        fmla    z18.s,   z0.s,  z4.s[2]
        fmla    z19.s,   z0.s,  z4.s[3]
    
        add     pA, pA, #16
        add     pB, pB, #48
.endm

.macro SAVE4x12
        fmul    z8.s, p4/m, z8.s, alpha0
        st1w    {z8.s}, p4, [pCRow0]
        add     pCRow0, pCRow0, #16

        fmul    z9.s, p4/m, z9.s, alpha0
        st1w    {z9.s}, p4, [pCRow1]
        add     pCRow1, pCRow1, #16

        fmul    z10.s, p4/m, z10.s, alpha0
        st1w    {z10.s}, p4, [pCRow2]
        add     pCRow2, pCRow2, #16

        fmul    z11.s, p4/m, z11.s, alpha0
        st1w    {z11.s}, p4, [pCRow3]
        add     pCRow3, pCRow3, #16

        fmul    z12.s, p4/m, z12.s, alpha0
        st1w    {z12.s}, p4, [pCRow4]
        add     pCRow4, pCRow4, #16

        fmul    z13.s, p4/m, z13.s, alpha0
        st1w    {z13.s}, p4, [pCRow5]
        add     pCRow5, pCRow5, #16

        fmul    z14.s, p4/m, z14.s, alpha0
        st1w    {z14.s}, p4, [pCRow6]
        add     pCRow6, pCRow6, #16

        fmul    z15.s, p4/m, z15.s, alpha0
        st1w    {z15.s}, p4, [pCRow7]
        add     pCRow7, pCRow7, #16

        fmul    z16.s, p4/m, z16.s, alpha0
        st1w    {z16.s}, p4, [pCRow8]
        add     pCRow8, pCRow8, #16

        fmul    z17.s, p4/m, z17.s, alpha0
        st1w    {z17.s}, p4, [pCRow9]
        add     pCRow9, pCRow9, #16

        fmul    z18.s, p4/m, z18.s, alpha0
        st1w    {z18.s}, p4, [pCRow10]
        add     pCRow10, pCRow10, #16

        fmul    z19.s, p4/m, z19.s, alpha0
        st1w    {z19.s}, p4, [pCRow11]
        add     pCRow11, pCRow11, #16
.endm

/******************************************************************************/
.macro INIT2x12
        mov     z9.s,  #0
        mov     z10.s, #0
        mov     z11.s, #0
        mov     z12.s, #0

        mov     z13.s, #0
        mov     z14.s, #0
        mov     z15.s, #0
        mov     z16.s, #0

        mov     z17.s, #0
        mov     z18.s, #0
        mov     z19.s, #0
        mov     z20.s, #0
.endm

/* macros KERNEL2x12_SUB and SAVE2x12 are copied from
 *  sgemm_kernel_sve256_16x12.S with modifications:
 *    - added predicated limit into fmul and fmla 
 */

/* register allocation:
 *   z0			: A
 *   z2,z3,z4		: B
 *   z9,z10......z20	: C
 *
 *   p4			: VL4
 *   p2			: VL2
 */

/* TODO: the main problem of this implementation is:
 *    - It doesn't make full use of a 256-bit vector length. In each fmul/fmla
 *      only 2-singles are calculated, which on its full width SVE-256 can calculate 8.
 */
.macro KERNEL2x12_SUB
        ld1rqw  {z0.s}, p2/z, [pA]

        ld1rqw  {z2.s}, p4/z, [pB]
        ld1rqw  {z3.s}, p4/z, [pB, #16]
        ld1rqw  {z4.s}, p4/z, [pB, #32]

        fmla    z9.s,   z0.s,  z2.s[0]
        fmla    z10.s,   z0.s,  z2.s[1]
        fmla    z11.s,   z0.s,  z2.s[2]
        fmla    z12.s,   z0.s,  z2.s[3]

        fmla    z13.s,   z0.s,  z3.s[0]
        fmla    z14.s,   z0.s,  z3.s[1]
        fmla    z15.s,   z0.s,  z3.s[2]
        fmla    z16.s,   z0.s,  z3.s[3]

        fmla    z17.s,   z0.s,  z4.s[0]
        fmla    z18.s,   z0.s,  z4.s[1]
        fmla    z19.s,   z0.s,  z4.s[2]
        fmla    z20.s,   z0.s,  z4.s[3]

        add     pA, pA, #8
        add     pB, pB, #48
.endm


.macro SAVE2x12
        fmul    z9.s, p2/m, z9.s, alpha0
        st1w    {z9.s}, p2, [pCRow0]
        add     pCRow0, pCRow0, #8

        fmul    z10.s, p2/m, z10.s, alpha0
        st1w    {z10.s}, p2, [pCRow1]
        add     pCRow1, pCRow1, #8
        
        fmul    z11.s, p2/m, z11.s, alpha0
        st1w    {z11.s}, p2, [pCRow2]
        add     pCRow2, pCRow2, #8

        fmul    z12.s, p2/m, z12.s, alpha0
        st1w    {z12.s}, p2, [pCRow3]
        add     pCRow3, pCRow3, #8

        fmul    z13.s, p2/m, z13.s, alpha0
        st1w    {z13.s}, p2, [pCRow4]
        add     pCRow4, pCRow4, #8

        fmul    z14.s, p2/m, z14.s, alpha0
        st1w    {z14.s}, p2, [pCRow5]
        add     pCRow5, pCRow5, #8

        fmul    z15.s, p2/m, z15.s, alpha0
        st1w    {z15.s}, p2, [pCRow6]
        add     pCRow6, pCRow6, #8

        fmul    z16.s, p2/m, z16.s, alpha0
        st1w    {z16.s}, p2, [pCRow7]
        add     pCRow7, pCRow7, #8

        fmul    z17.s, p2/m, z17.s, alpha0
        st1w    {z17.s}, p2, [pCRow8]
        add     pCRow8, pCRow8, #8

        fmul    z18.s, p2/m, z18.s, alpha0
        st1w    {z18.s}, p2, [pCRow9]
        add     pCRow9, pCRow9, #8
        
        fmul    z19.s, p2/m, z19.s, alpha0
        st1w    {z19.s}, p2, [pCRow10]
        add     pCRow10, pCRow10, #8

        fmul    z20.s, p2/m, z20.s, alpha0
        st1w    {z20.s}, p2, [pCRow11]
        add     pCRow11, pCRow11, #8
.endm

/******************************************************************************/
.macro INIT1x12
        mov     z8.s,  #0
        mov     z9.s,  #0
        mov     z10.s, #0
        mov     z11.s, #0

        mov     z12.s, #0
        mov     z13.s, #0
        mov     z14.s, #0
        mov     z15.s, #0

        mov     z16.s, #0
        mov     z17.s, #0
        mov     z18.s, #0
        mov     z19.s, #0
.endm

/* macros KERNEL1x12_SUB and SAVE1x12 are copied from
 *  sgemm_kernel_sve256_16x12.S with modifications:
 *    - added predicated limit into fmul and fmla 
 */

/* register allocation:
 *   z0			: A
 *   z2,z3,z4		: B
 *   z8,z9......z19	: C
 *
 *   p4			: VL4
 *   p1			: VL1
 */

/* TODO: the main problem of this implementation is:
 *    - It doesn't make full use of a 256-bit vector length. In each fmul/fmla
 *      only 1-singles are calculated, which on its full width SVE-256 can calculate 8.
 */
.macro KERNEL1x12_SUB
        ld1rw   {z0.s}, p1/z, [pA]

        ld1rqw  {z2.s}, p4/z, [pB]
        ld1rqw  {z3.s}, p4/z, [pB, #16]
        ld1rqw  {z4.s}, p4/z, [pB, #32]

        fmla    z8.s,   z0.s,  z2.s[0]
        fmla    z9.s,   z0.s,  z2.s[1]
        fmla    z10.s,   z0.s,  z2.s[2]
        fmla    z11.s,   z0.s,  z2.s[3]

        fmla    z12.s,   z0.s,  z3.s[0]
        fmla    z13.s,   z0.s,  z3.s[1]
        fmla    z14.s,   z0.s,  z3.s[2]
        fmla    z15.s,   z0.s,  z3.s[3]

        fmla    z16.s,   z0.s,  z4.s[0]
        fmla    z17.s,   z0.s,  z4.s[1]
        fmla    z18.s,   z0.s,  z4.s[2]
        fmla    z19.s,   z0.s,  z4.s[3]

        add     pA, pA, #4
        add     pB, pB, #48
.endm

.macro SAVE1x12
        fmul    z8.s, p1/m, z8.s, alpha0
        st1w    {z8.s}, p1, [pCRow0]
        add     pCRow0, pCRow0, #4

        fmul    z9.s, p1/m, z9.s, alpha0
        st1w    {z9.s}, p1, [pCRow1]
        add    pCRow1, pCRow1, #4

        fmul    z10.s, p1/m, z10.s, alpha0
        st1w    {z10.s}, p1, [pCRow2]
        add     pCRow2, pCRow2, #4

        fmul    z11.s, p1/m, z11.s, alpha0
        st1w    {z11.s}, p1, [pCRow3]
        add     pCRow3, pCRow3, #4

        fmul    z12.s, p1/m, z12.s, alpha0
        st1w    {z12.s}, p1, [pCRow4]
        add     pCRow4, pCRow4, #4

        fmul    z13.s, p1/m, z13.s, alpha0
        st1w    {z13.s}, p1, [pCRow5]
        add     pCRow5, pCRow5, #4

        fmul    z14.s, p1/m, z14.s, alpha0
        st1w    {z14.s}, p1, [pCRow6]
        add     pCRow6, pCRow6, #4

        fmul    z15.s, p1/m, z15.s, alpha0
        st1w    {z15.s}, p1, [pCRow7]
        add     pCRow7, pCRow7, #4

        fmul    z16.s, p1/m, z16.s, alpha0
        st1w    {z16.s}, p1, [pCRow8]
        add     pCRow8, pCRow8, #4

        fmul    z17.s, p1/m, z17.s, alpha0
        st1w    {z17.s}, p1, [pCRow9]
        add     pCRow9, pCRow9, #4

        fmul    z18.s, p1/m, z18.s, alpha0
        st1w    {z18.s}, p1, [pCRow10]
        add     pCRow10, pCRow10, #4

        fmul    z19.s, p1/m, z19.s, alpha0
        st1w    {z19.s}, p1, [pCRow11]
        add     pCRow11, pCRow11, #4
.endm

/******************************************************************************/

.macro INIT16x8
        mov     z8.s,  #0
        mov     z9.s,  #0
        mov     z10.s, #0
        mov     z11.s, #0

        mov     z12.s, #0
        mov     z13.s, #0
        mov     z14.s, #0
        mov     z15.s, #0

        mov     z16.s, #0
        mov     z17.s, #0
        mov     z18.s, #0
        mov     z19.s, #0

        mov     z20.s, #0
        mov     z21.s, #0
        mov     z22.s, #0
        mov     z23.s, #0
.endm

/* macros KERNEL16x8_SUB and SAVE16x8 are copied from
 *  sgemm_kernel_sve256_16x12.S
 */

/* register allocation:
 *   z0,z1		: A
 *   z4,z5		: B
 *   z8,z9......z23	: C
 *
 *   p0			: all true
 */

.macro KERNEL16x8_SUB
        ld1w    {z0.s}, p0/z, [pA]
        ld1w    {z1.s}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        ld1rqw  {z4.s}, p0/z, [pB]

        fmla    z8.s, z0.s, z4.s[0]
        fmla    z10.s, z0.s, z4.s[1]
        fmla    z12.s, z0.s, z4.s[2]
        fmla    z14.s, z0.s, z4.s[3]

        ld1rqw  {z5.s}, p0/z, [pB, #16]
        add     pB, pB, #32

        fmla    z9.s, z1.s, z4.s[0]
        fmla    z11.s, z1.s, z4.s[1]
        fmla    z13.s, z1.s, z4.s[2]
        fmla    z15.s, z1.s, z4.s[3]

        fmla    z16.s, z0.s, z5.s[0]
        fmla    z18.s, z0.s, z5.s[1]
        fmla    z20.s, z0.s, z5.s[2]
        fmla    z22.s, z0.s, z5.s[3]
        fmla    z17.s, z1.s, z5.s[0]
        fmla    z19.s, z1.s, z5.s[1]
        fmla    z21.s, z1.s, z5.s[2]
        fmla    z23.s, z1.s, z5.s[3]
.endm

.macro SAVE16x8
        fmul    z0.s, z8.s, alpha0
        fmul    z1.s, z9.s, alpha0
        st1w    {z0.s}, p0, [pCRow0]
        st1w    {z1.s}, p0, [pCRow0, #1, MUL VL]
        incb    pCRow0, all, MUL #2

        fmul    z4.s, z10.s, alpha0
        fmul    z5.s, z11.s, alpha0
        st1w    {z4.s}, p0, [pCRow1]
        st1w    {z5.s}, p0, [pCRow1, #1, MUL VL]
        incb    pCRow1, all, MUL #2

        fmul    z0.s, z12.s, alpha0
        fmul    z1.s, z13.s, alpha0
        st1w    {z0.s}, p0, [pCRow2]
        st1w    {z1.s}, p0, [pCRow2, #1, MUL VL]
        incb    pCRow2, all, MUL #2

        fmul    z4.s, z14.s, alpha0
        fmul    z5.s, z15.s, alpha0
        st1w    {z4.s}, p0, [pCRow3]
        st1w    {z5.s}, p0, [pCRow3, #1, MUL VL]
        incb    pCRow3, all, MUL #2

        fmul    z0.s, z16.s, alpha0
        fmul    z1.s, z17.s, alpha0
        st1w    {z0.s}, p0, [pCRow4]
        st1w    {z1.s}, p0, [pCRow4, #1, MUL VL]
        incb    pCRow4, all, MUL #2

        fmul    z4.s, z18.s, alpha0
        fmul    z5.s, z19.s, alpha0
        st1w    {z4.s}, p0, [pCRow5]
        st1w    {z5.s}, p0, [pCRow5, #1, MUL VL]
        incb    pCRow5, all, MUL #2

        fmul    z0.s, z20.s, alpha0
        fmul    z1.s, z21.s, alpha0
        st1w    {z0.s}, p0, [pCRow6]
        st1w    {z1.s}, p0, [pCRow6, #1, MUL VL]
        incb    pCRow6, all, MUL #2

        fmul    z4.s, z22.s, alpha0
        fmul    z5.s, z23.s, alpha0
        st1w    {z4.s}, p0, [pCRow7]
        st1w    {z5.s}, p0, [pCRow7, #1, MUL VL]
        incb    pCRow7, all, MUL #2
.endm

/******************************************************************************/

.macro INIT8x8
        mov     z8.s,  #0
        mov     z9.s,  #0
        mov     z10.s, #0
        mov     z11.s, #0

        mov     z12.s, #0
        mov     z13.s, #0
        mov     z14.s, #0
        mov     z15.s, #0
.endm

/* macros KERNEL8x8_SUB and SAVE8x8 are copied from
 *  sgemm_kernel_sve256_16x12.S
 */

/* register allocation:
 *   z0			: A
 *   z2,z3		: B
 *   z8,z9......z15	: C
 *
 *   p0			: all true
 */

.macro KERNEL8x8_SUB
        ld1w    {z0.s}, p0/z, [pA]
        add     pA, pA, #32
        ld1rqw  {z2.s}, p0/z, [pB]
        ld1rqw  {z3.s}, p0/z, [pB, #16]
        add     pB, pB, #32

        fmla    z8.s, z0.s, z2.s[0]
        fmla    z9.s, z0.s, z2.s[1]
        fmla    z10.s, z0.s, z2.s[2]
        fmla    z11.s, z0.s, z2.s[3]
        fmla    z12.s, z0.s, z3.s[0]
        fmla    z13.s, z0.s, z3.s[1]
        fmla    z14.s, z0.s, z3.s[2]
        fmla    z15.s, z0.s, z3.s[3]
.endm

.macro SAVE8x8
        fmul    z0.s, z8.s, alpha0
        st1w    {z0.s}, p0, [pCRow0]
        add     pCRow0, pCRow0, #32

        fmul    z1.s, z9.s, alpha0
        st1w    {z1.s}, p0, [pCRow1]
        add     pCRow1, pCRow1, #32

        fmul    z2.s, z10.s, alpha0
        st1w    {z2.s}, p0, [pCRow2]
        add     pCRow2, pCRow2, #32
        
        fmul    z4.s, z11.s, alpha0
        st1w    {z4.s}, p0, [pCRow3]
        add     pCRow3, pCRow3, #32

        fmul    z0.s, z12.s, alpha0
        st1w    {z0.s}, p0, [pCRow4]
        add     pCRow4, pCRow4, #32

        fmul    z1.s, z13.s, alpha0
        st1w    {z1.s}, p0, [pCRow5]
        add     pCRow5, pCRow5, #32

        fmul    z2.s, z14.s, alpha0
        st1w    {z2.s}, p0, [pCRow6]
        add     pCRow6, pCRow6, #32

        fmul    z4.s, z15.s, alpha0
        st1w    {z4.s}, p0, [pCRow7]
        add     pCRow7, pCRow7, #32
.endm

/******************************************************************************/

.macro INIT4x8
        mov     z8.s,  #0
        mov     z9.s,  #0
        mov     z10.s, #0
        mov     z11.s, #0

        mov     z12.s, #0
        mov     z13.s, #0
        mov     z14.s, #0
        mov     z15.s, #0
.endm

/* macros KERNEL4x8_SUB and SAVE4x8 are copied from
 *  sgemm_kernel_sve256_16x12.S
 */

/* register allocation:
 *   z0			: A
 *   z2,z3		: B
 *   z8,z9......z15	: C
 *
 *   p4			: VL4
 */
.macro KERNEL4x8_SUB
        ld1rqw  {z0.s}, p4/z, [pA]
        ld1rqw  {z2.s}, p4/z, [pB]
        ld1rqw  {z3.s}, p4/z, [pB, #16]

        fmla    z8.s, z0.s, z2.s[0]
        fmla    z9.s, z0.s, z2.s[1]
        fmla    z10.s, z0.s, z2.s[2]
        fmla    z11.s, z0.s, z2.s[3]
        fmla    z12.s, z0.s, z3.s[0]
        fmla    z13.s, z0.s, z3.s[1]
        fmla    z14.s, z0.s, z3.s[2]
        fmla    z15.s, z0.s, z3.s[3]

        add pA, pA, #16
        add pB, pB, #32
.endm

.macro SAVE4x8
        fmul    z8.s, p4/m, z8.s, alpha0
        st1w    {z8.s}, p4, [pCRow0]
        add     pCRow0, pCRow0, #16
        
        fmul    z9.s, p4/m, z9.s, alpha0
        st1w    {z9.s}, p4, [pCRow1]
        add     pCRow1, pCRow1, #16

        fmul    z10.s, p4/m, z10.s, alpha0
        st1w    {z10.s}, p4, [pCRow2]
        add     pCRow2, pCRow2, #16

        fmul    z11.s, p4/m, z11.s, alpha0
        st1w    {z11.s}, p4, [pCRow3]
        add     pCRow3, pCRow3, #16

        fmul    z12.s, p4/m, z12.s, alpha0
        st1w    {z12.s}, p4, [pCRow4]
        add     pCRow4, pCRow4, #16

        fmul    z13.s, p4/m, z13.s, alpha0
        st1w    {z13.s}, p4, [pCRow5]
        add     pCRow5, pCRow5, #16

        fmul    z14.s, p4/m, z14.s, alpha0
        st1w    {z14.s}, p4, [pCRow6]
        add     pCRow6, pCRow6, #16

        fmul    z15.s, p4/m, z15.s, alpha0
        st1w    {z15.s}, p4, [pCRow7]
        add     pCRow7, pCRow7, #16
.endm

/******************************************************************************/

.macro INIT2x8
        mov     z9.s,  #0
        mov     z10.s, #0
        mov     z11.s, #0
        mov     z12.s, #0

        mov     z13.s, #0
        mov     z14.s, #0
        mov     z15.s, #0
        mov     z16.s, #0
.endm

/* macros KERNEL2x8_SUB and SAVE2x8 are copied from
 *  sgemm_kernel_sve256_16x12.S
 */

/* register allocation:
 *   z0			: A
 *   z2,z3		: B
 *   z9,z10......z16	: C
 *
 *   p2			: VL2
 *   p4			: VL4
 */
.macro KERNEL2x8_SUB
        ld1rqw  z0.s, p2/z, [pA]
        ld1rqw  z2.s, p4/z, [pB]
        ld1rqw  z3.s, p4/z, [pB, #16]

        fmla    z9.s, z0.s, z2.s[0]
        fmla    z10.s, z0.s, z2.s[1]
        fmla    z11.s, z0.s, z2.s[2]
        fmla    z12.s, z0.s, z2.s[3]
        fmla    z13.s, z0.s, z3.s[0]
        fmla    z14.s, z0.s, z3.s[1]
        fmla    z15.s, z0.s, z3.s[2]
        fmla    z16.s, z0.s, z3.s[3]

        add     pA, pA, #8
        add     pB, pB, #32
.endm

.macro SAVE2x8
        fmul    z9.s, p2/m, z9.s, alpha0
        st1w    {z9.s}, p2, [pCRow0]
        add     pCRow0, pCRow0, #8

        fmul    z10.s, p2/m, z10.s, alpha0
        st1w    {z10.s}, p2, [pCRow1]
        add     pCRow1, pCRow1, #8

        fmul    z11.s, p2/m, z11.s, alpha0
        st1w    {z11.s}, p2, [pCRow2]
        add     pCRow2, pCRow2, #8

        fmul    z12.s, p2/m, z12.s, alpha0
        st1w    {z12.s}, p2, [pCRow3]
        add     pCRow3, pCRow3, #8

        fmul    z13.s, p2/m, z13.s, alpha0
        st1w    {z13.s}, p2, [pCRow4]
        add     pCRow4, pCRow4, #8

        fmul    z14.s, p2/m, z14.s, alpha0
        st1w    {z14.s}, p2, [pCRow5]
        add    pCRow5, pCRow5, #8

        fmul    z15.s, p2/m, z15.s, alpha0
        st1w    {z15.s}, p2, [pCRow6]
        add     pCRow6, pCRow6, #8

        fmul    z16.s, p2/m, z16.s, alpha0
        st1w    {z16.s}, p2, [pCRow7]
        add     pCRow7, pCRow7, #8
.endm

/******************************************************************************/

.macro INIT1x8
        mov     z8.s,  #0
        mov     z9.s,  #0
        mov     z10.s, #0
        mov     z11.s, #0

        mov     z12.s, #0
        mov     z13.s, #0
        mov     z14.s, #0
        mov     z15.s, #0
.endm

/* macros KERNEL1x8_SUB and SAVE1x8 are copied from
 *  sgemm_kernel_sve256_16x12.S
 */

/* register allocation:
 *   z0			: A
 *   z2,z3		: B
 *   z8,z9......z15	: C
 *
 *   p1			: VL1
 *   p4			: VL4
 */
.macro KERNEL1x8_SUB
        ld1rw   {z0.s}, p1/z, [pA]

        ld1rqw  {z2.s}, p4/z, [pB]
        ld1rqw  {z3.s}, p4/z, [pB, #16]

        fmla    z8.s, z0.s, z2.s[0]
        fmla    z9.s, z0.s, z2.s[1]
        fmla    z10.s, z0.s, z2.s[2]
        fmla    z11.s, z0.s, z2.s[3]
        fmla    z12.s, z0.s, z3.s[0]
        fmla    z13.s, z0.s, z3.s[1]
        fmla    z14.s, z0.s, z3.s[2]
        fmla    z15.s, z0.s, z3.s[3]

        add     pA, pA, #4
        add     pB, pB, #32
.endm

.macro SAVE1x8
        fmul    z8.s, p1/m, z8.s, alpha0
        st1w    {z8.s}, p1, [pCRow0]
        add     pCRow0, pCRow0, #4

        fmul    z9.s, p1/m, z9.s, alpha0
        st1w    {z9.s}, p1, [pCRow1]
        add     pCRow1, pCRow1, #4

        fmul    z10.s, p1/m, z10.s, alpha0
        st1w    {z10.s}, p1, [pCRow2]
        add     pCRow2, pCRow2, #4

        fmul    z11.s, p1/m, z11.s, alpha0
        st1w    {z11.s}, p1, [pCRow3]
        add     pCRow3, pCRow3, #4

        fmul    z12.s, p1/m, z12.s, alpha0
        st1w    {z12.s}, p1, [pCRow4]
        add     pCRow4, pCRow4, #4

        fmul    z13.s, p1/m, z13.s, alpha0
        st1w    {z13.s}, p1, [pCRow5]
        add     pCRow5, pCRow5, #4

        fmul    z14.s, p1/m, z14.s, alpha0
        st1w    {z14.s}, p1, [pCRow6]
        add     pCRow6, pCRow6, #4

        fmul    z15.s, p1/m, z15.s, alpha0
        st1w    {z15.s}, p1, [pCRow7]
        add     pCRow7, pCRow7, #4
.endm

/******************************************************************************/

.macro INIT16x4
        mov     z8.s,  #0
        mov     z9.s,  #0
        mov     z10.s, #0
        mov     z11.s, #0

        mov     z12.s, #0
        mov     z13.s, #0
        mov     z14.s, #0
        mov     z15.s, #0
.endm

/* macros KERNEL16x4_SUB and SAVE16x4 are copied from
 *  sgemm_kernel_sve256_16x12.S
 */

/* register allocation:
 *   z0,z1		: A
 *   z4			: B
 *   z8,z9......z15	: C
 *
 *   p0			: all true
 */

.macro KERNEL16x4_SUB
        ld1w    {z0.s}, p0/z, [pA]
        ld1w    {z1.s}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        ld1rqw  {z4.s}, p0/z, [pB]
        add     pB, pB, #16

        fmla    z8.s, z0.s, z4.s[0]
        fmla    z10.s, z0.s, z4.s[1]
        fmla    z12.s, z0.s, z4.s[2]
        fmla    z14.s, z0.s, z4.s[3]
        fmla    z9.s, z1.s, z4.s[0]
        fmla    z11.s, z1.s, z4.s[1]
        fmla    z13.s, z1.s, z4.s[2]
        fmla    z15.s, z1.s, z4.s[3]
.endm

.macro SAVE16x4
        fmul    z0.s, z8.s, alpha0
        fmul    z1.s, z9.s, alpha0
        st1w    {z0.s}, p0, [pCRow0]
        st1w    {z1.s}, p0, [pCRow0, #1, MUL VL]
        incb    pCRow0, all, MUL #2

        fmul    z4.s, z10.s, alpha0
        fmul    z5.s, z11.s, alpha0
        st1w    {z4.s}, p0, [pCRow1]
        st1w    {z5.s}, p0, [pCRow1, #1, MUL VL]
        incb    pCRow1, all, MUL #2

        fmul    z0.s, z12.s, alpha0
        fmul    z1.s, z13.s, alpha0
        st1w    {z0.s}, p0, [pCRow2]
        st1w    {z1.s}, p0, [pCRow2, #1, MUL VL]
        incb    pCRow2, all, MUL #2

        fmul    z4.s, z14.s, alpha0
        fmul    z5.s, z15.s, alpha0
        st1w    {z4.s}, p0, [pCRow3]
        st1w    {z5.s}, p0, [pCRow3, #1, MUL VL]
        incb    pCRow3, all, MUL #2
.endm

/******************************************************************************/

.macro INIT8x4
        mov     z8.s,  #0
        mov     z9.s,  #0
        mov     z10.s, #0
        mov     z11.s, #0
.endm

/* macros KERNEL8x4_SUB and SAVE8x4 are copied from
 *  sgemm_kernel_sve256_16x12.S
 */

/* register allocation:
 *   z0			: A
 *   z4			: B
 *   z8,z9......z11	: C
 *
 *   p0			: all true
 */
.macro KERNEL8x4_SUB
        ld1w    {z0.s}, p0/z, [pA]
        ld1rqw  {z4.s}, p0/z, [pB]

        fmla    z8.s, z0.s, z4.s[0]
        fmla    z9.s, z0.s, z4.s[1]
        fmla    z10.s, z0.s, z4.s[2]
        fmla    z11.s, z0.s, z4.s[3]

        add     pA, pA, #32
        add     pB, pB, #16
.endm

.macro SAVE8x4
        fmul    z0.s, z8.s, alpha0
        st1w    {z0.s}, p0, [pCRow0]
        add     pCRow0, pCRow0, #32

        fmul    z4.s, z9.s, alpha0
        st1w    {z4.s}, p0, [pCRow1]
        add     pCRow1, pCRow1, #32

        fmul    z0.s, z10.s, alpha0
        st1w    {z0.s}, p0, [pCRow2]
        add     pCRow2, pCRow2, #32

        fmul    z4.s, z11.s, alpha0
        st1w    {z4.s}, p0, [pCRow3]
        add     pCRow3, pCRow3, #32
.endm

/******************************************************************************/

.macro INIT4x4
        mov     z8.s,  #0
        mov     z9.s,  #0
        mov     z10.s, #0
        mov     z11.s, #0
.endm

/* macros KERNEL8x4_SUB and SAVE8x4 are copied from
 *  sgemm_kernel_sve256_16x12.S
 */

/* register allocation:
 *   z0			: A
 *   z4			: B
 *   z8,z9......z11	: C
 *
 *   p4			: VL4
 */
.macro KERNEL4x4_SUB
        ld1rqw  {z0.s}, p4/z, [pA]
        ld1rqw  {z4.s}, p4/z, [pB]

        fmla    z8.s, z0.s, z4.s[0]
        fmla    z9.s, z0.s, z4.s[1]
        fmla    z10.s, z0.s, z4.s[2]
        fmla    z11.s, z0.s, z4.s[3]

        add     pA, pA, #16
        add     pB, pB, #16
.endm

.macro SAVE4x4
        fmul    z8.s, p4/m, z8.s, alpha0
        st1w    {z8.s}, p4, [pCRow0]
        add     pCRow0, pCRow0, #16
        
        fmul    z9.s, p4/m, z9.s, alpha0
        st1w    {z9.s}, p4, [pCRow1]
        add     pCRow1, pCRow1, #16

        fmul    z10.s, p4/m, z10.s, alpha0
        st1w    {z10.s}, p4, [pCRow2]
        add     pCRow2, pCRow2, #16

        fmul    z11.s, p4/m, z11.s, alpha0
        st1w    {z11.s}, p4, [pCRow3]
        add     pCRow3, pCRow3, #16
.endm

/******************************************************************************/

.macro INIT2x4
        mov     z9.s,  #0
        mov     z10.s, #0
        mov     z11.s, #0
        mov     z12.s, #0
.endm

/* macros KERNEL2x4_SUB and SAVE2x4 are copied from
 *  sgemm_kernel_sve256_16x12.S
 */

/* register allocation:
 *   z0			: A
 *   z4			: B
 *   z9,z10......z12	: C
 *
 *   p4			: VL4
 *   p2			: VL2
 */

.macro KERNEL2x4_SUB
        ld1rqw  {z0.s}, p2/z, [pA]
        ld1rqw  {z4.s}, p4/z, [pB]

        fmla    z9.s, z0.s, z4.s[0]
        fmla    z10.s, z0.s, z4.s[1]
        fmla    z11.s, z0.s, z4.s[2]
        fmla    z12.s, z0.s, z4.s[3]

        add     pA, pA, #8
        add     pB, pB, #16
.endm

.macro SAVE2x4
        fmul    z9.s, p2/m, z9.s, alpha0
        st1w    {z9.s}, p2, [pCRow0]
        add     pCRow0, pCRow0, #8

        fmul    z10.s, p2/m, z10.s, alpha0
        st1w    {z10.s}, p2, [pCRow1]
        add     pCRow1, pCRow1, #8

        fmul    z11.s, p2/m, z11.s, alpha0
        st1w    {z11.s}, p2, [pCRow2]
        add     pCRow2, pCRow2, #8

        fmul    z12.s, p2/m, z12.s, alpha0
        st1w    {z12.s}, p2, [pCRow3]
        add     pCRow3, pCRow3, #8
.endm

/******************************************************************************/

.macro INIT1x4
        mov     z8.s,  #0
        mov     z9.s,  #0
        mov     z10.s, #0
        mov     z11.s, #0
.endm

/* macros KERNEL1x4_SUB and SAVE1x4 are copied from
 *  sgemm_kernel_sve256_16x12.S
 */

/* register allocation:
 *   z0			: A
 *   z4			: B
 *   z8,z9......z11	: C
 *
 *   p4			: VL4
 *   p1			: VL1
 */

.macro KERNEL1x4_SUB
        ld1rw   {z0.s}, p1/z, [pA]
        ld1rqw  {z4.s}, p4/z, [pB]

        fmla    z8.s, z0.s, z4.s[0]
        fmla    z9.s, z0.s, z4.s[1]
        fmla    z10.s, z0.s, z4.s[2]
        fmla    z11.s, z0.s, z4.s[3]

        add     pA, pA, #4
        add     pB, pB, #16
.endm

.macro SAVE1x4
        fmul    z8.s, p1/m, z8.s, alpha0
        st1w    {z8.s}, p1, [pCRow0]
        add     pCRow0, pCRow0, #4

        fmul    z9.s, p1/m, z9.s, alpha0
        st1w    {z9.s}, p1, [pCRow1]
        add     pCRow1, pCRow1, #4

        fmul    z10.s, p1/m, z10.s, alpha0
        st1w    {z10.s}, p1, [pCRow2]
        add     pCRow2, pCRow2, #4

        fmul    z11.s, p1/m, z11.s, alpha0
        st1w    {z11.s}, p1, [pCRow3]
        add     pCRow3, pCRow3, #4
.endm

/******************************************************************************/

.macro INIT16x2
        mov     z8.s,  #0
        mov     z9.s,  #0
        mov     z10.s, #0
        mov     z11.s, #0
.endm

/* macros KERNEL16x2_SUB and SAVE16x2 are copied from
 *  sgemm_kernel_sve256_16x12.S
 *  , with modifications of adding p2/z into [pB] ld1rqw.
 */

/* register allocation:
 *   z0,z1		: A
 *   z4			: B
 *   z8,z9......z11	: C
 *
 *   p0			: all true
 *   p2			: VL2
 */

.macro KERNEL16x2_SUB
        ld1w    {z0.s}, p0/z, [pA]
        ld1w    {z1.s}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        ld1rqw  {z4.s}, p2/z, [pB]
        add     pB, pB, #8

        fmla    z8.s, z0.s, z4.s[0]
        fmla    z10.s, z0.s, z4.s[1]
        fmla    z9.s, z1.s, z4.s[0]
        fmla    z11.s, z1.s, z4.s[1]
.endm

.macro SAVE16x2
        fmul    z0.s, z8.s, alpha0
        fmul    z1.s, z9.s, alpha0
        st1w    {z0.s}, p0, [pCRow0]
        st1w    {z1.s}, p0, [pCRow0, #1, MUL VL]
        incb    pCRow0, all, MUL #2

        fmul    z4.s, z10.s, alpha0
        fmul    z5.s, z11.s, alpha0
        st1w    {z4.s}, p0, [pCRow1]
        st1w    {z5.s}, p0, [pCRow1, #1, MUL VL]
        incb    pCRow1, all, MUL #2
.endm

/******************************************************************************/

.macro INIT8x2
        mov     z8.s,  #0
        mov     z9.s,  #0
.endm

/* macros KERNEL8x2_SUB and SAVE8x2 are copied from
 *  sgemm_kernel_sve256_16x12.S
 *  , with modifications of adding p2/z into [pB] ld1rqw.
 */

/* register allocation:
 *   z0			: A
 *   z4			: B
 *   z8,z9		: C
 *
 *   p0			: all true
 *   p2			: VL2
 */

.macro KERNEL8x2_SUB
        ld1w    {z0.s}, p0/z, [pA]
        ld1rqw  {z4.s}, p2/z, [pB]

        fmla    z8.s, z0.s, z4.s[0]
        fmla    z9.s, z0.s, z4.s[1]

        add     pA, pA, #32
        add     pB, pB, #8
.endm

.macro SAVE8x2
        fmul    z0.s, z8.s, alpha0
        st1w    {z0.s}, p0, [pCRow0]
        add     pCRow0, pCRow0, #32

        fmul    z4.s, z9.s, alpha0
        st1w    {z4.s}, p0, [pCRow1]
        add     pCRow1, pCRow1, #32
.endm

/******************************************************************************/

.macro INIT4x2
        mov     z8.s,  #0
        mov     z9.s,  #0
.endm

/* macros KERNEL4x2_SUB and SAVE4x2 are copied from
 *  sgemm_kernel_sve256_16x12.S
 *  , with modifications of adding p2/z into [pB] ld1rqw.
 */

/* register allocation:
 *   z0			: A
 *   z4			: B
 *   z8,z9		: C
 *
 *   p4			: VL4
 *   p2			: VL2
 */

.macro KERNEL4x2_SUB
        ld1rqw  {z0.s}, p4/z, [pA]
        ld1rqw  {z4.s}, p2/z, [pB]

        fmla    z8.s, z0.s, z4.s[0]
        fmla    z9.s, z0.s, z4.s[1]

        add     pA, pA, #16
        add     pB, pB, #8
.endm

.macro SAVE4x2
        fmul    z8.s, p4/m, z8.s, alpha0
        st1w    {z8.s}, p4, [pCRow0]
        add     pCRow0, pCRow0, #16

        fmul    z9.s, p4/m, z9.s, alpha0
        st1w    {z9.s}, p4, [pCRow1]
        add     pCRow1, pCRow1, #16
.endm

/******************************************************************************/

.macro INIT2x2
        mov     z8.s,  #0
        mov     z9.s,  #0
.endm

/* macros KERNEL2x2_SUB and SAVE2x2 are copied from
 *  sgemm_kernel_sve256_16x12.S
 *  , with modifications of adding p2/z into [pB] ld1rqw.
 */

/* register allocation:
 *   z0			: A
 *   z4			: B
 *   z8,z9		: C
 *
 *   p2			: VL2
 */

.macro KERNEL2x2_SUB
        ld1rqw  {z0.s}, p2/z, [pA]
        ld1rqw  {z4.s}, p2/z, [pB]

        fmla    z8.s, z0.s, z4.s[0]
        fmla    z9.s, z0.s, z4.s[1]

        add     pA, pA, #8
        add     pB, pB, #8
.endm

.macro SAVE2x2
        fmul    z8.s, p2/m, z8.s, alpha0
        st1w    {z8.s}, p2, [pCRow0]
        add     pCRow0, pCRow0, #8

        fmul    z9.s, p2/m, z9.s, alpha0
        st1w    {z9.s}, p2, [pCRow1]
        add     pCRow1, pCRow1, #8
.endm

/******************************************************************************/

.macro INIT1x2
        mov     z8.s,  #0
        mov     z9.s,  #0
.endm

/* macros KERNEL1x2_SUB and SAVE1x2 are copied from
 *  sgemm_kernel_sve256_16x12.S
 */

/* register allocation:
 *   z0			: A
 *   z4			: B
 *   z8,z9		: C
 *
 *   p1			: VL1
 *   p2			: VL2
 */

.macro KERNEL1x2_SUB
        ld1rw   {z0.s}, p1/z, [pA]
        ld1rqw  {z4.s}, p2/z, [pB]

        fmla    z8.s, z0.s, z4.s[0]
        fmla    z9.s, z0.s, z4.s[1]

        add     pA, pA, #4
        add     pB, pB, #8
.endm

.macro SAVE1x2
        fmul    z8.s, p1/m, z8.s, alpha0
        st1w    {z8.s}, p1, [pCRow0]
        add     pCRow0, pCRow0, #4

        fmul    z9.s, p1/m, z9.s, alpha0
        st1w    {z9.s}, p1, [pCRow1]
        add     pCRow1, pCRow1, #4
.endm

/******************************************************************************/

.macro INIT16x1
        mov     z8.s,  #0
        mov     z9.s,  #0
.endm

/* macros KERNEL16x1_SUB and SAVE16x1 are copied from
 *  sgemm_kernel_sve256_16x12.S
 */

/* register allocation:
 *   z0,z1		: A
 *   z4			: B
 *   z8,z9		: C
 *
 *   p0			: all true
 *   p1			: VL1
 */

.macro KERNEL16x1_SUB
        ld1w    {z0.s}, p0/z, [pA]
        ld1w    {z1.s}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        ld1rqw  {z4.s}, p1/z, [pB]
        add     pB, pB, #4

        fmla    z8.s, z0.s, z4.s[0]
        fmla    z9.s, z1.s, z4.s[0]
.endm

.macro SAVE16x1
        fmul    z0.s, z8.s, alpha0
        fmul    z1.s, z9.s, alpha0
        st1w    {z0.s}, p0, [pCRow0]
        st1w    {z1.s}, p0, [pCRow0, #1, MUL VL]
        incb    pCRow0, all, MUL #2
.endm

/******************************************************************************/

.macro INIT8x1
        mov     z8.s,  #0
.endm

/* macros KERNEL8x1_SUB and SAVE8x1 are copied from
 *  sgemm_kernel_sve256_16x12.S
 */

/* register allocation:
 *   z0			: A
 *   z4			: B
 *   z8			: C
 *
 *   p0			: all true
 *   p1			: VL1
 */
.macro KERNEL8x1_SUB
        ld1w    {z0.s}, p0/z, [pA]
        ld1rqw  {z4.s}, p1/z, [pB]

        fmla    z8.s, z0.s, z4.s[0]

        add     pA, pA, #32
        add     pB, pB, #4
.endm

.macro SAVE8x1
        fmul    z0.s, z8.s, alpha0
        st1w    {z0.s}, p0, [pCRow0]
        add     pCRow0, pCRow0, #32
.endm

/******************************************************************************/

.macro INIT4x1
        mov     z8.s,  #0
.endm

/* macros KERNEL4x1_SUB and SAVE4x1 are copied from
 *  sgemm_kernel_sve256_16x12.S
 */

/* register allocation:
 *   z0			: A
 *   z4			: B
 *   z8			: C
 *
 *   p4			: VL4
 *   p1			: VL1
 */

.macro KERNEL4x1_SUB
        ld1rqw {z0.s}, p4/z, [pA]
        ld1rw  {z4.s}, p1/z, [pB]

        fmla    z8.s, z0.s, z4.s[0]

        add     pA, pA, #16
        add     pB, pB, #4
.endm

.macro SAVE4x1
        fmul    z8.s, p4/m, z8.s, alpha0
        st1w    {z8.s}, p4, [pCRow0]
        add     pCRow0, pCRow0, #16
.endm

/******************************************************************************/

.macro INIT2x1
        mov     z8.s,  #0
.endm

/* macros KERNEL2x1_SUB and SAVE2x1 are copied from
 *  sgemm_kernel_sve256_16x12.S
 */

/* register allocation:
 *   z0			: A
 *   z4			: B
 *   z8			: C
 *
 *   p2			: VL2
 *   p1			: VL1
 */
.macro KERNEL2x1_SUB
        ld1rqw  {z0.s}, p2/z, [pA]
        ld1rw   {z4.s}, p1/z, [pB]

        fmla    z8.s, z0.s, z4.s[0]

        add     pA, pA, #8
        add     pB, pB, #4
.endm

.macro SAVE2x1
        fmul    z8.s, p2/m, z8.s, alpha0
        st1w    {z8.s}, p2, [pCRow0]
        add     pCRow0, pCRow0, #8
.endm

/******************************************************************************/

.macro INIT1x1
        mov     z16.s,  #0
.endm

/* macros KERNEL1x1_SUB and SAVE1x1 are copied from
 *  sgemm_kernel_sve256_16x12.S
 */

/* register allocation:
 *   z0			: A
 *   z8			: B
 *   z16		: C
 *
 *   p1			: VL1
 */
.macro KERNEL1x1_SUB
        ld1rw   z0.s, p1/z, [pA]
        ld1rw   z8.s, p1/z, [pB]

        add     pA , pA, #4
        add     pB , pB, #4

        fmla    z16.s, p1/m, z0.s, z8.s
.endm

.macro SAVE1x1
        fmul    z16.s, p1/m, z16.s, alpha0
        st1w    {z16.s}, p1, [pCRow0]

        add     pCRow0, pCRow0, #4
.endm

/*******************************************************************************
* End of macro definitions
*******************************************************************************/

	PROLOGUE

.Lstrmm_kernel_begin:

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

	fmov	alpha, s0			// load in general-purpose register
	ptrue	p0.s, all			// p0 is all enabled
	ptrue   p4.s, VL4			// p4 is VL4
	ptrue   p2.s, VL2			// p2 is VL2
	ptrue   p1.s, VL1			// p1 is VL1

	lsl	LDC, LDC, #2			// ldc = ldc * 4
						// 4 bytes per single float

#if !defined(LEFT)
	neg	tempOffset, offset
#endif
	mov	pB, origPB

	mov	counterJ, origN
	cmp	counterJ, #11			// 11 = (12-1)
	ble	.Lstrmm_kernel_L8_BEGIN		// If less than 12 column, skip L12

/******************************************************************************/

.Lstrmm_kernel_L12_BEGIN:
        mov    pCRow0, pC
        add    pCRow1, pCRow0, LDC
        add    pCRow2, pCRow1, LDC
        add    pCRow3, pCRow2, LDC
        add    pCRow4, pCRow3, LDC
        add    pCRow5, pCRow4, LDC
        add    pCRow6, pCRow5, LDC
        add    pCRow7, pCRow6, LDC
        add    pCRow8, pCRow7, LDC
        add    pCRow9, pCRow8, LDC
        add    pCRow10, pCRow9, LDC
        add    pCRow11, pCRow10, LDC

	add	pC, pCRow11, LDC

#if defined(LEFT)
	mov	tempOffset, offset
#endif
	mov	pA, origPA			// pA = start of A array

.Lstrmm_kernel_L12_M16_BEGIN:

	mov	counterI, origM
	asr 	counterI, counterI, #4		// counterI = counterI / 16
	cmp 	counterI, #0
	ble	.Lstrmm_kernel_L12_M8_BEGIN	// skip L12_M16

	.align 5
.Lstrmm_kernel_L12_M16_20:

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #6	// 2^6 = 16 * sizeof(single), where 16 is M16
	add	pA, pA, temp
					// 12 * sizeof(single), where 12 is L12
					// 12 = (16-4), with that to change mutiplication of 12 to
					//      substraction of (times of 16) and (times of 4)
					// 2^6 = 16 * sizeof(single)
					// 2^4 = 4  * sizeof(single)
	add	pB, pB, temp
	sub	pB, pB, tempOffset, lsl #4
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #16	// M16
#else
	add	tempK, tempOffset, #12	// L12
#endif

	asr	counterL , tempK, #3		// L = K / 8, unroll at 8
	cmp	counterL , #2			// is there at least 4 to do?
	blt	.Lstrmm_kernel_L12_M16_32

	KERNEL16x12_I
	KERNEL16x12_M2
	KERNEL16x12_M1
	KERNEL16x12_M2

	KERNEL16x12_M1
	KERNEL16x12_M2
	KERNEL16x12_M1
	KERNEL16x12_M2

	subs	counterL, counterL, #2
	ble	.Lstrmm_kernel_L12_M16_22a

	.align 5
.Lstrmm_kernel_L12_M16_22:

	KERNEL16x12_M1
	KERNEL16x12_M2
	KERNEL16x12_M1
	KERNEL16x12_M2

	KERNEL16x12_M1
	KERNEL16x12_M2
	KERNEL16x12_M1
	KERNEL16x12_M2

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L12_M16_22

	.align 5
.Lstrmm_kernel_L12_M16_22a:

	KERNEL16x12_M1
	KERNEL16x12_M2
	KERNEL16x12_M1
	KERNEL16x12_M2

	KERNEL16x12_M1
	KERNEL16x12_M2
	KERNEL16x12_M1
	KERNEL16x12_E

	b	 .Lstrmm_kernel_L12_M16_44

	.align 5
.Lstrmm_kernel_L12_M16_32:

	tst	counterL, #1
	ble	.Lstrmm_kernel_L12_M16_40

	KERNEL16x12_I
	KERNEL16x12_M2
	KERNEL16x12_M1
	KERNEL16x12_M2

	KERNEL16x12_M1
	KERNEL16x12_M2
	KERNEL16x12_M1
	KERNEL16x12_E

	b	.Lstrmm_kernel_L12_M16_44

.Lstrmm_kernel_L12_M16_40:

	INIT16x12

.Lstrmm_kernel_L12_M16_44:

	ands	counterL , tempK, #7
	ble	.Lstrmm_kernel_L12_M16_100

	.align 5
.Lstrmm_kernel_L12_M16_46:

	KERNEL16x12_SUB
	subs	counterL, counterL, #1
	bne	.Lstrmm_kernel_L12_M16_46

.Lstrmm_kernel_L12_M16_100:

	SAVE16x12

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #16	// M16
#else
	sub	tempK, tempK, #12	// L12
#endif
	lsl	temp, tempK, #6		// 2^6 = 16 * sizeof(single), 16 is M16
	add	pA, pA, temp
					// 12 * sizeof(single), where 12 is L12
					// 12 = (16-4), with that to change mutiplication of 12 to
					//      substraction of (times of 16) and (times of 4)
					// 2^6 = 16 * sizeof(single)
					// 2^4 = 4  * sizeof(single)
	add	pB, pB, temp
	sub	pB, pB, tempK, lsl #4
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #16	// 16 is M16
#endif
/* TODO: adjust prfm
	prfm	PLDL1KEEP, [pA]
	prfm	PLDL1KEEP, [pA, #64]
	prfm	PLDL1KEEP, [origPB]
 */

.Lstrmm_kernel_L12_M16_END:
	subs	counterI, counterI, #1
	bne	.Lstrmm_kernel_L12_M16_20

//------------------------------------------------------------------------------

.Lstrmm_kernel_L12_M8_BEGIN:
        dup     alpha0, alpha			// z31 (alpha0) is not used except in L12_M16
						// also Note that we use z3 (alpha0_t) as alpha in L12_M16

	mov	counterI, origM
	tst	counterI , #15			// 15 is from (M16-1)
	ble	.Lstrmm_kernel_L12_END

	tst	counterI, #8			// 8 is from M8
	ble	.Lstrmm_kernel_L12_M4_BEGIN

.Lstrmm_kernel_L12_M8_20:

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
						// 12*sizeof(single) =
						//     (2^6 - 2^4)
						// 12 is L12
	add	pB, pB, tempOffset, lsl #6	// pB for L
	sub	pB, pB, tempOffset, lsl #4
						// 2^5 = 8*sizeof(single), 8 is M8
	add	pA, pA, tempOffset, lsl #5	// pA for M
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #8		// 8 is M8
#else
	add	tempK, tempOffset, #12		// 12 is L12
#endif

	asr 	counterL , tempK, #1		// L = K / 2
	cmp	counterL , #2			// is there at least 4 to do?
	blt	.Lstrmm_kernel_L12_M8_32

	KERNEL8x12_I				// do one in the K
	KERNEL8x12_M2				// do another in the K

	subs	counterL, counterL, #2
	ble	.Lstrmm_kernel_L12_M8_22a
	.align 5

.Lstrmm_kernel_L12_M8_22:

	KERNEL8x12_M1
	KERNEL8x12_M2

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L12_M8_22

.Lstrmm_kernel_L12_M8_22a:

	KERNEL8x12_M1
	KERNEL8x12_E

	b	 .Lstrmm_kernel_L12_M8_44

.Lstrmm_kernel_L12_M8_32:

	tst	counterL, #1
	ble	.Lstrmm_kernel_L12_M8_40

	KERNEL8x12_I
	KERNEL8x12_E

	b	.Lstrmm_kernel_L12_M8_44

.Lstrmm_kernel_L12_M8_40:

	INIT8x12

.Lstrmm_kernel_L12_M8_44:

	ands	counterL , tempK, #1
	ble	.Lstrmm_kernel_L12_M8_100

.Lstrmm_kernel_L12_M8_46:

	KERNEL8x12_SUB

.Lstrmm_kernel_L12_M8_100:

	SAVE8x12

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #8		// 8 is M8
#else
	sub	tempK, tempK, #12		// 12 is L12
#endif
	add	pA, pA, tempK, lsl #5		// pA moves M8. 2^5=8*sizeof(single)
						// pB moves L12.
						// 12*sizeof(single) = 2^6 - 2^4
	add	pB, pB, tempK, lsl #6
	sub	pB, pB, tempK, lsl #4
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #8	// 8 is M8
#endif

.Lstrmm_kernel_L12_M8_END:

//------------------------------------------------------------------------------

.Lstrmm_kernel_L12_M4_BEGIN:

	mov	counterI, origM
	tst	counterI , #7			// 7 is (M8 -1)
	ble	.Lstrmm_kernel_L12_END

	tst	counterI, #4			// M4
	ble	.Lstrmm_kernel_L12_M2_BEGIN

.Lstrmm_kernel_L12_M4_20:

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
						// 12*sizeof(single) =
						//     (2^6 - 2^4)
						// 12 is L12
	add	pB, pB, tempOffset, lsl #6	// pB for L
	sub	pB, pB, tempOffset, lsl #4
						// pA for M
	add	pA, pA, tempOffset, lsl #4	// 2^4=4*sizeof(single), 4 is M4
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #4		// M4
#else
	add	tempK, tempOffset, #12		// L12
#endif
	asr 	counterL , tempK, #1		// L = K / 2
	cmp	counterL , #2			// is there at least 4 to do?
	blt	.Lstrmm_kernel_L12_M4_32

	KERNEL4x12_I				// do one in the K
	KERNEL4x12_M2				// do another in the K

	subs	counterL, counterL, #2
	ble	.Lstrmm_kernel_L12_M4_22a
	.align 5

.Lstrmm_kernel_L12_M4_22:

	KERNEL4x12_M1
	KERNEL4x12_M2

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L12_M4_22

.Lstrmm_kernel_L12_M4_22a:

	KERNEL4x12_M1
	KERNEL4x12_E

	b	 .Lstrmm_kernel_L12_M4_44

.Lstrmm_kernel_L12_M4_32:

	tst	counterL, #1
	ble	.Lstrmm_kernel_L12_M4_40

	KERNEL4x12_I
	KERNEL4x12_E

	b	.Lstrmm_kernel_L12_M4_44

.Lstrmm_kernel_L12_M4_40:

	INIT4x12

.Lstrmm_kernel_L12_M4_44:

	ands	counterL , tempK, #1
	ble	.Lstrmm_kernel_L12_M4_100

.Lstrmm_kernel_L12_M4_46:

	KERNEL4x12_SUB

.Lstrmm_kernel_L12_M4_100:

	SAVE4x12

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #4		// 4 is M4
#else
	sub	tempK, tempK, #12		// 12 is L12
#endif
	add	pA, pA, tempK, lsl #4		// pA moves M4. 2^4=4*sizeof(single)
						// pB moves L12.
						// 12*sizeof(single) = 2^6 - 2^4
	add	pB, pB, tempK, lsl #6
	sub	pB, pB, tempK, lsl #4
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #4	// 4 is M4
#endif
.Lstrmm_kernel_L12_M4_END:

//------------------------------------------------------------------------------

.Lstrmm_kernel_L12_M2_BEGIN:

	mov	counterI, origM
	tst	counterI , #3			// 3 is (M4 - 1)
	ble	.Lstrmm_kernel_L12_END

	tst	counterI, #2			// M2
	ble	.Lstrmm_kernel_L12_M1_BEGIN

.Lstrmm_kernel_L12_M2_20:

	INIT2x12

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	add	pA, pA, tempOffset, lsl #3	// pA for M, 2*sizeof(single) = 2^3
						// 12*sizeof(single) =
						//     (2^6 - 2^4)
						// 12 is L12
	add	pB, pB, tempOffset, lsl #6	// pB for L
	sub	pB, pB, tempOffset, lsl #4
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #2		// 2 is M2
#else
	add	tempK, tempOffset, #12		// 12 is L12
#endif
	asr 	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Lstrmm_kernel_L12_M2_40

.Lstrmm_kernel_L12_M2_22:

	KERNEL2x12_SUB
	KERNEL2x12_SUB
	KERNEL2x12_SUB
	KERNEL2x12_SUB

	KERNEL2x12_SUB
	KERNEL2x12_SUB
	KERNEL2x12_SUB
	KERNEL2x12_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L12_M2_22


.Lstrmm_kernel_L12_M2_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L12_M2_100

.Lstrmm_kernel_L12_M2_42:

	KERNEL2x12_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L12_M2_42

.Lstrmm_kernel_L12_M2_100:

	SAVE2x12

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #2		// M2
#else
	sub	tempK, tempK, #12		// L12
#endif
	add	pA, pA, tempK, lsl #3		// pA moves M2. 2^3=2*sizeof(single)
						// pB moves L12.
						// 12*sizeof(single) = 2^6 - 2^4
	add	pB, pB, tempK, lsl #6
	sub	pB, pB, tempK, lsl #4
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #2	// M2
#endif
.Lstrmm_kernel_L12_M2_END:


.Lstrmm_kernel_L12_M1_BEGIN:

	tst	counterI, #1			// counterI = counterI % 2
	ble	.Lstrmm_kernel_L12_END

.Lstrmm_kernel_L12_M1_20:

	INIT1x12

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
						// 12*sizeof(single) =
						//     (2^6 - 2^4)
						// 12 is L12
	add	pB, pB, tempOffset, lsl #6	// pB for L
	sub	pB, pB, tempOffset, lsl #4
						// pA for M
	add	pA, pA, tempOffset, lsl #2	// 2^2=1*sizeof(single), 1 is M1
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #1		// M1
#else
	add	tempK, tempOffset, #12		// L12
#endif
	asr 	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Lstrmm_kernel_L12_M1_40

.Lstrmm_kernel_L12_M1_22:
	KERNEL1x12_SUB
	KERNEL1x12_SUB
	KERNEL1x12_SUB
	KERNEL1x12_SUB

	KERNEL1x12_SUB
	KERNEL1x12_SUB
	KERNEL1x12_SUB
	KERNEL1x12_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L12_M1_22


.Lstrmm_kernel_L12_M1_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L12_M1_100

.Lstrmm_kernel_L12_M1_42:

	KERNEL1x12_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L12_M1_42

.Lstrmm_kernel_L12_M1_100:

	SAVE1x12

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #1		// M1
#else
	sub	tempK, tempK, #12		// L12
#endif
	add	pA, pA, tempK, lsl #2		// pA moves M1. 2^2=1*sizeof(single)
						// pB moves L12.
						// 12*sizeof(single) = 2^6 - 2^4
	add	pB, pB, tempK, lsl #6
	sub	pB, pB, tempK, lsl #4
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #1	// M1
#endif
.Lstrmm_kernel_L12_END:
						// B = B + K * 12 * 4, where 12 is L12, 4 is sizeof(single)
						// B = B + K*16*4 - K*4*4
	add	origPB, origPB, origK, lsl #6
	sub	origPB, origPB, origK, lsl #4
#if !defined(LEFT)
	add	tempOffset, tempOffset, #12	// L12
#endif

	sub	counterJ, counterJ , #12
	cmp	counterJ, #11
	bgt	.Lstrmm_kernel_L12_BEGIN

/******************************************************************************/

.Lstrmm_kernel_L8_BEGIN:   // less than 12 left in N direction
        dup     alpha0, alpha			// z31 (alpha0) is not used except in L12_M16
						// also Note that we use z3 (alpha0_t) as alpha in L12_M16


	cmp	counterJ , #0
	ble	.Lstrmm_kernel_L999

	cmp	counterJ , #8			// L8
	blt	.Lstrmm_kernel_L4_BEGIN

	mov	pCRow0, pC			// pCRow0 = pC
        add    pCRow1, pCRow0, LDC
        add    pCRow2, pCRow1, LDC
        add    pCRow3, pCRow2, LDC
        add    pCRow4, pCRow3, LDC
        add    pCRow5, pCRow4, LDC
        add    pCRow6, pCRow5, LDC
        add    pCRow7, pCRow6, LDC

        add    pC, pCRow7, LDC

#if defined(LEFT)
	mov	tempOffset, offset
#endif
	mov	pA, origPA			// pA = A

.Lstrmm_kernel_L8_M16_BEGIN:

	mov	counterI, origM
	asr	counterI, counterI, #4		// counterI = counterI / 16
	cmp	counterI,#0
	ble	.Lstrmm_kernel_L8_M8_BEGIN

.Lstrmm_kernel_L8_M16_20:

	INIT16x8

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
/* TODO: these lsl instructions and 'temp' register are not necessary
 *     - replaece them with: add pA, pA, tempOffset, lsl #6
 *     - And find all other similar cases in this file.
 */
	lsl	temp, tempOffset, #6		// 16*sizeof(single) = 2^6
	add	pA, pA, temp			// pA for M
	lsl	temp, tempOffset, #5		// 8 *sizeof(single) = 2^5
	add	pB, pB, temp			// pB for L
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #16		// M16
#else
	add	tempK, tempOffset, #8		// L8
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL,#0
	ble	.Lstrmm_kernel_L8_M16_40
	.align 5

.Lstrmm_kernel_L8_M16_22:
/* TODO optimize these batches */
	KERNEL16x8_SUB
	KERNEL16x8_SUB
	KERNEL16x8_SUB
	KERNEL16x8_SUB

	KERNEL16x8_SUB
	KERNEL16x8_SUB
	KERNEL16x8_SUB
	KERNEL16x8_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L8_M16_22

.Lstrmm_kernel_L8_M16_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L8_M16_100

.Lstrmm_kernel_L8_M16_42:

	KERNEL16x8_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L8_M16_42

.Lstrmm_kernel_L8_M16_100:

	SAVE16x8

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #16		// M16
#else
	sub	tempK, tempK, #8		// L8
#endif
						// 16*sizeof(single) = 2^6
	add	pA, pA, tempK, lsl #6		// PA for M
						// 8 *sizeof(single) = 2^5
	add	pB, pB, tempK, lsl #5		// pB for L
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #16	// M16
#endif

.Lstrmm_kernel_L8_M16_END:

	subs	counterI, counterI, #1
	bgt	.Lstrmm_kernel_L8_M16_20

//------------------------------------------------------------------------------

.Lstrmm_kernel_L8_M8_BEGIN:
	mov	counterI, origM
	tst	counterI , #15			// M16 - 1
	ble	.Lstrmm_kernel_L8_END

	tst	counterI, #8			// M8
	ble	.Lstrmm_kernel_L8_M4_BEGIN

.Lstrmm_kernel_L8_M8_20:

	INIT8x8

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #5		// 8*sizeof(single) = 2^5, M8 and L8
	add	pA, pA, temp			// pA for M
	add	pB, pB, temp			// pB for L
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #8		// M8
#else
	add	tempK, tempOffset, #8		// L8
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL,#0
	ble	.Lstrmm_kernel_L8_M8_40
	.align 5

.Lstrmm_kernel_L8_M8_22:
	KERNEL8x8_SUB
	KERNEL8x8_SUB
	KERNEL8x8_SUB
	KERNEL8x8_SUB

	KERNEL8x8_SUB
	KERNEL8x8_SUB
	KERNEL8x8_SUB
	KERNEL8x8_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L8_M8_22


.Lstrmm_kernel_L8_M8_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L8_M8_100

.Lstrmm_kernel_L8_M8_42:

	KERNEL8x8_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L8_M8_42

.Lstrmm_kernel_L8_M8_100:

	SAVE8x8

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #8		// M8
#else
	sub	tempK, tempK, #8		// L8
#endif
	lsl	temp, tempK, #5			// 8*sizeof(single) = 2^5, M8/L8
	add	pA, pA, temp
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #8	// M8
#endif

.Lstrmm_kernel_L8_M8_END:

//------------------------------------------------------------------------------

.Lstrmm_kernel_L8_M4_BEGIN:
	mov	counterI, origM
	tst	counterI , #7			// M8-1
	ble	.Lstrmm_kernel_L8_END

	tst	counterI, #4			// M4
	ble	.Lstrmm_kernel_L8_M2_BEGIN

.Lstrmm_kernel_L8_M4_20:

	INIT4x8

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #5		// 8*sizeof(single) = 2^5, L8
	add	pB, pB, temp			// pB for L
	lsl	temp, tempOffset, #4		// 4*sizeof(single) = 2^4, M4
	add	pA, pA, temp			// pA for M
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #4		// M4
#else
	add	tempK, tempOffset, #8		// L8
#endif
	asr	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL,#0
	ble	.Lstrmm_kernel_L8_M4_40
	.align 5

.Lstrmm_kernel_L8_M4_22:
	KERNEL4x8_SUB
	KERNEL4x8_SUB
	KERNEL4x8_SUB
	KERNEL4x8_SUB

	KERNEL4x8_SUB
	KERNEL4x8_SUB
	KERNEL4x8_SUB
	KERNEL4x8_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L8_M4_22


.Lstrmm_kernel_L8_M4_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L8_M4_100

.Lstrmm_kernel_L8_M4_42:

	KERNEL4x8_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L8_M4_42

.Lstrmm_kernel_L8_M4_100:

	SAVE4x8

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #4		// M4
#else
	sub	tempK, tempK, #8		// L8
#endif
	lsl	temp, tempK, #4			// 4*sizeof(single)=2^4, M4
	add	pA, pA, temp			// pA for M
	lsl	temp, tempK, #5			// 8*sizeof(single)=2^5, L8
	add	pB, pB, temp			// pB for L
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #4	// M4
#endif
.Lstrmm_kernel_L8_M4_END:

//------------------------------------------------------------------------------

.Lstrmm_kernel_L8_M2_BEGIN:

	mov	counterI, origM
	tst	counterI , #3			// M4 - 1
	ble	.Lstrmm_kernel_L8_END

	tst	counterI, #2			// M2
	ble	.Lstrmm_kernel_L8_M1_BEGIN

.Lstrmm_kernel_L8_M2_20:

	INIT2x8

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #3		// 2*sizeof(single)=2^3, M2
	add	pA, pA, temp			// pA for M
	lsl	temp, tempOffset, #5		// 8*sizeof(single)=2^5, L8
	add	pB, pB, temp			// pB for L
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #2		// M2
#else
	add	tempK, tempOffset, #8		// L8
#endif
	asr	counterL , tempK, #3		// counterL = counterL / 8
        cmp	counterL,#0
	ble	.Lstrmm_kernel_L8_M2_40

.Lstrmm_kernel_L8_M2_22:

	KERNEL2x8_SUB
	KERNEL2x8_SUB
	KERNEL2x8_SUB
	KERNEL2x8_SUB

	KERNEL2x8_SUB
	KERNEL2x8_SUB
	KERNEL2x8_SUB
	KERNEL2x8_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L8_M2_22


.Lstrmm_kernel_L8_M2_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L8_M2_100

.Lstrmm_kernel_L8_M2_42:

	KERNEL2x8_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L8_M2_42

.Lstrmm_kernel_L8_M2_100:

	SAVE2x8
#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #2		// M2
#else
	sub	tempK, tempK, #8		// L8
#endif
	lsl	temp, tempK, #3			// 2*sizeof(single) = 2^3
	add	pA, pA, temp			// pA for M
	lsl	temp, tempK, #5			// 8*sizeof(single) = 2^5
	add	pB, pB, temp			// pB for L
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #2	// M2
#endif

.Lstrmm_kernel_L8_M2_END:


.Lstrmm_kernel_L8_M1_BEGIN:

	tst	counterI, #1			// M1
	ble	.Lstrmm_kernel_L8_END

.Lstrmm_kernel_L8_M1_20:

	INIT1x8

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #2		// 1*sizeof(single)=2^2, M1
	add	pA, pA, temp			// pA for M
	lsl	temp, tempOffset, #5		// 8*sizeof(single)=2^5, L8
	add	pB, pB, temp
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #1		// M1
#else
	add	tempK, tempOffset, #8		// L8
#endif
	asr 	counterL , tempK, #3		// counterL = counterL / 8
        cmp     counterL, #0
	ble	.Lstrmm_kernel_L8_M1_40

.Lstrmm_kernel_L8_M1_22:
	KERNEL1x8_SUB
	KERNEL1x8_SUB
	KERNEL1x8_SUB
	KERNEL1x8_SUB

	KERNEL1x8_SUB
	KERNEL1x8_SUB
	KERNEL1x8_SUB
	KERNEL1x8_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L8_M1_22


.Lstrmm_kernel_L8_M1_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L8_M1_100

.Lstrmm_kernel_L8_M1_42:

	KERNEL1x8_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L8_M1_42

.Lstrmm_kernel_L8_M1_100:

	SAVE1x8

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #1		// M1
#else
	sub	tempK, tempK, #8		// L8
#endif
	lsl	temp, tempK, #2			// 1*sizeof(single)=2^2
	add	pA, pA, temp			// pA for M
	lsl	temp, tempK, #5			// 8*sizeof(single)=2^5
	add	pB, pB, temp			// pB for L
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #1	// M1
#endif
.Lstrmm_kernel_L8_END:
#if !defined(LEFT)
	add	tempOffset, tempOffset, #8	// L8
#endif
	add	origPB, origPB, origK, lsl #5	// B = B + K * 8 * 4, where 8 is L8, 4 is sizeof(single)

	sub	counterJ, counterJ, #8		// L8
/******************************************************************************/

.Lstrmm_kernel_L4_BEGIN:   // less than 8 left in N direction

	cmp	counterJ , #4			// L4
	blt	.Lstrmm_kernel_L2_BEGIN

	mov	pCRow0, pC			// pCRow0 = pC
        add    pCRow1, pCRow0, LDC
        add    pCRow2, pCRow1, LDC
        add    pCRow3, pCRow2, LDC

        add    pC, pCRow3, LDC

#if defined(LEFT)
	mov	tempOffset, offset
#endif
	mov	pA, origPA			// pA = A

.Lstrmm_kernel_L4_M16_BEGIN:

	mov	counterI, origM
	asr	counterI, counterI, #4		// counterI = counterI / 16
	cmp	counterI,#0
	ble	.Lstrmm_kernel_L4_M8_BEGIN

.Lstrmm_kernel_L4_M16_20:

	INIT16x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
/* TODO: these lsl instructions and 'temp' register are not necessary
 *     - replaece them with: add pA, pA, tempOffset, lsl #6
 *     - And find all other similar cases in this file.
 */
	lsl	temp, tempOffset, #6		// 16*sizeof(single) = 2^6
	add	pA, pA, temp			// pA for M
	lsl	temp, tempOffset, #4		// 4 *sizeof(single) = 2^4
	add	pB, pB, temp			// pB for L
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #16		// M16
#else
	add	tempK, tempOffset, #4		// L4
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL,#0
	ble	.Lstrmm_kernel_L4_M16_40
	.align 5

.Lstrmm_kernel_L4_M16_22:
/* TODO optimize these batches */
	KERNEL16x4_SUB
	KERNEL16x4_SUB
	KERNEL16x4_SUB
	KERNEL16x4_SUB

	KERNEL16x4_SUB
	KERNEL16x4_SUB
	KERNEL16x4_SUB
	KERNEL16x4_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L4_M16_22

.Lstrmm_kernel_L4_M16_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L4_M16_100

.Lstrmm_kernel_L4_M16_42:

	KERNEL16x4_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L4_M16_42

.Lstrmm_kernel_L4_M16_100:

	SAVE16x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #16		// M16
#else
	sub	tempK, tempK, #4		// L4
#endif
						// 16*sizeof(single) = 2^6
	add	pA, pA, tempK, lsl #6		// PA for M
						// 4 *sizeof(single) = 2^4
	add	pB, pB, tempK, lsl #4		// pB for L
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #16	// M16
#endif

.Lstrmm_kernel_L4_M16_END:

	subs	counterI, counterI, #1
	bgt	.Lstrmm_kernel_L4_M16_20

//------------------------------------------------------------------------------

.Lstrmm_kernel_L4_M8_BEGIN:
	mov	counterI, origM
	tst	counterI , #15			// M16 - 1
	ble	.Lstrmm_kernel_L4_END

	tst	counterI, #8			// M8
	ble	.Lstrmm_kernel_L4_M4_BEGIN

.Lstrmm_kernel_L4_M8_20:

	INIT8x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #5		// 8*sizeof(single) = 2^5, M8
	add	pA, pA, temp			// pA for M
	lsl	temp, tempOffset, #4		// 4*sizeof(single) = 2^4, L4
	add	pB, pB, temp			// pB for L
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #8		// M8
#else
	add	tempK, tempOffset, #4		// L4
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL,#0
	ble	.Lstrmm_kernel_L4_M8_40
	.align 5

.Lstrmm_kernel_L4_M8_22:
	KERNEL8x4_SUB
	KERNEL8x4_SUB
	KERNEL8x4_SUB
	KERNEL8x4_SUB

	KERNEL8x4_SUB
	KERNEL8x4_SUB
	KERNEL8x4_SUB
	KERNEL8x4_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L4_M8_22

.Lstrmm_kernel_L4_M8_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L4_M8_100

.Lstrmm_kernel_L4_M8_42:

	KERNEL8x4_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L4_M8_42

.Lstrmm_kernel_L4_M8_100:

	SAVE8x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #8		// M8
#else
	sub	tempK, tempK, #4		// L4
#endif
	lsl	temp, tempK, #5			// 8*sizeof(single) = 2^5, M8
	add	pA, pA, temp			// pA for M
	lsl	temp, tempK, #4			// 4*sizeof(single) = 2^4, L4
	add	pB, pB, temp			// pB for L
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #8	// M8
#endif

.Lstrmm_kernel_L4_M8_END:

//------------------------------------------------------------------------------

.Lstrmm_kernel_L4_M4_BEGIN:
	mov	counterI, origM			// todo: this instruction is not necessary.
						//  - counterI is initialized since L4_BEGIN.
	tst	counterI , #7			// M8-1
	ble	.Lstrmm_kernel_L4_END

	tst	counterI, #4			// M4
	ble	.Lstrmm_kernel_L4_M2_BEGIN

.Lstrmm_kernel_L4_M4_20:

	INIT4x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #4		// 4*sizeof(single) = 2^4, M4
	add	pA, pA, temp			// pA for M
	lsl	temp, tempOffset, #4		// 4*sizeof(single) = 2^4, L4
	add	pB, pB, temp			// pB for L
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #4		// M4
#else
	add	tempK, tempOffset, #4		// L4
#endif
	asr	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL,#0
	ble	.Lstrmm_kernel_L4_M4_40
	.align 5

.Lstrmm_kernel_L4_M4_22:
	KERNEL4x4_SUB
	KERNEL4x4_SUB
	KERNEL4x4_SUB
	KERNEL4x4_SUB

	KERNEL4x4_SUB
	KERNEL4x4_SUB
	KERNEL4x4_SUB
	KERNEL4x4_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L4_M4_22


.Lstrmm_kernel_L4_M4_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L4_M4_100

.Lstrmm_kernel_L4_M4_42:

	KERNEL4x4_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L4_M4_42

.Lstrmm_kernel_L4_M4_100:

	SAVE4x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #4		// M4
#else
	sub	tempK, tempK, #4		// L4
#endif
	lsl	temp, tempK, #4			// 4*sizeof(single)=2^4, M4/L4
	add	pA, pA, temp			// pA for M
	add	pB, pB, temp			// pB for L
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #4	// M4
#endif
.Lstrmm_kernel_L4_M4_END:

//------------------------------------------------------------------------------

.Lstrmm_kernel_L4_M2_BEGIN:

	mov	counterI, origM
	tst	counterI , #3			// M4 - 1
	ble	.Lstrmm_kernel_L4_END

	tst	counterI, #2			// M2
	ble	.Lstrmm_kernel_L4_M1_BEGIN

.Lstrmm_kernel_L4_M2_20:

	INIT2x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #3		// 2*sizeof(single)=2^3, M2
	add	pA, pA, temp			// pA for M
	lsl	temp, tempOffset, #4		// 4*sizeof(single)=2^4, L4
	add	pB, pB, temp			// pB for L
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #2		// M2
#else
	add	tempK, tempOffset, #4		// L4
#endif
	asr	counterL , tempK, #3		// counterL = counterL / 8
        cmp	counterL,#0
	ble	.Lstrmm_kernel_L4_M2_40

.Lstrmm_kernel_L4_M2_22:

	KERNEL2x4_SUB
	KERNEL2x4_SUB
	KERNEL2x4_SUB
	KERNEL2x4_SUB

	KERNEL2x4_SUB
	KERNEL2x4_SUB
	KERNEL2x4_SUB
	KERNEL2x4_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L4_M2_22


.Lstrmm_kernel_L4_M2_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L4_M2_100

.Lstrmm_kernel_L4_M2_42:

	KERNEL2x4_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L4_M2_42

.Lstrmm_kernel_L4_M2_100:

	SAVE2x4
#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #2		// M2
#else
	sub	tempK, tempK, #4		// L4
#endif
	lsl	temp, tempK, #3			// 2*sizeof(single) = 2^3, M2
	add	pA, pA, temp			// pA for M
	lsl	temp, tempK, #4			// 4*sizeof(single) = 2^4, L4
	add	pB, pB, temp			// pB for L
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #2	// M2
#endif

.Lstrmm_kernel_L4_M2_END:


.Lstrmm_kernel_L4_M1_BEGIN:

	tst	counterI, #1			// M1
	ble	.Lstrmm_kernel_L4_END

.Lstrmm_kernel_L4_M1_20:

	INIT1x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #2		// 1*sizeof(single)=2^2, M1
	add	pA, pA, temp			// pA for M
	lsl	temp, tempOffset, #4		// 4*sizeof(single)=2^4, L4
	add	pB, pB, temp
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #1		// M1
#else
	add	tempK, tempOffset, #4		// L4
#endif
	asr 	counterL , tempK, #3		// counterL = counterL / 8
        cmp     counterL, #0
	ble	.Lstrmm_kernel_L4_M1_40

.Lstrmm_kernel_L4_M1_22:
	KERNEL1x4_SUB
	KERNEL1x4_SUB
	KERNEL1x4_SUB
	KERNEL1x4_SUB

	KERNEL1x4_SUB
	KERNEL1x4_SUB
	KERNEL1x4_SUB
	KERNEL1x4_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L4_M1_22


.Lstrmm_kernel_L4_M1_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L4_M1_100

.Lstrmm_kernel_L4_M1_42:

	KERNEL1x4_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L4_M1_42

.Lstrmm_kernel_L4_M1_100:

	SAVE1x4

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #1		// M1
#else
	sub	tempK, tempK, #4		// L4
#endif
	lsl	temp, tempK, #2			// 1*sizeof(single)=2^2, M1
	add	pA, pA, temp			// pA for M
	lsl	temp, tempK, #4			// 4*sizeof(single)=2^4, L4
	add	pB, pB, temp			// pB for L
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #1	// M1
#endif
.Lstrmm_kernel_L4_END:
#if !defined(LEFT)
	add	tempOffset, tempOffset, #4	// L4
#endif
	add	origPB, origPB, origK, lsl #4	// B = B + K * 4 * 4, where 4 is L4, 4 is sizeof(single)

/******************************************************************************/

.Lstrmm_kernel_L2_BEGIN:   // less than 4 left in N direction

	// 12 can be divided by 4. So, any N before they can reach
	//   this step must be not '00' in their binary lowest two bits.
	//
	// That's why in L8/L4_BEGIN counterJ must not be reloaded as 'origN'.
	// That's why in L2_BEGIN there is no harm to reload counterJ as 'origN'

	mov	counterJ , origN
	tst	counterJ , #3
	ble	.Lstrmm_kernel_L999

	tst	counterJ , #2
	ble	.Lstrmm_kernel_L1_BEGIN

	mov	pCRow0, pC			// pCRow0 = pC
        add    pCRow1, pCRow0, LDC

	add	pC,pC,LDC, lsl #1

#if defined(LEFT)
	mov	tempOffset, offset
#endif
	mov	pA, origPA			// pA = A

.Lstrmm_kernel_L2_M16_BEGIN:

	mov	counterI, origM
	asr 	counterI, counterI, #4		// counterI = counterI / 16
	cmp	counterI,#0
	ble	.Lstrmm_kernel_L2_M8_BEGIN

.Lstrmm_kernel_L2_M16_20:

	INIT16x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #6
	add	pA, pA, temp
	lsl	temp, tempOffset, #3
	add	pB, pB, temp
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #16
#else
	add	tempK, tempOffset, #2
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL,#0
	ble	.Lstrmm_kernel_L2_M16_40
	.align 5

.Lstrmm_kernel_L2_M16_22:
	KERNEL16x2_SUB
	KERNEL16x2_SUB
	KERNEL16x2_SUB
	KERNEL16x2_SUB

	KERNEL16x2_SUB
	KERNEL16x2_SUB
	KERNEL16x2_SUB
	KERNEL16x2_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L2_M16_22


.Lstrmm_kernel_L2_M16_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L2_M16_100

.Lstrmm_kernel_L2_M16_42:

	KERNEL16x2_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L2_M16_42

.Lstrmm_kernel_L2_M16_100:

	SAVE16x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #16
#else
	sub	tempK, tempK, #2
#endif
	lsl	temp, tempK, #6
	add	pA, pA, temp
	lsl	temp, tempK, #3
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #16
#endif

.Lstrmm_kernel_L2_M16_END:

	subs	counterI, counterI, #1
	bgt	.Lstrmm_kernel_L2_M16_20

//------------------------------------------------------------------------------

.Lstrmm_kernel_L2_M8_BEGIN:
	mov	counterI, origM			// Todo: from here on, in L2_M4_BEGIN,
						//   L2_M2/M1_BEGIN, no need to reload counterI as origM.
	tst	counterI , #15
	ble	.Lstrmm_kernel_L2_END

	tst	counterI, #8
	ble	.Lstrmm_kernel_L2_M4_BEGIN

.Lstrmm_kernel_L2_M8_20:

	INIT8x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #5
	add	pA, pA, temp
	lsl	temp, tempOffset, #3
	add	pB, pB, temp
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #8
#else
	add	tempK, tempOffset, #2
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL,#0
	ble	.Lstrmm_kernel_L2_M8_40
	.align 5

.Lstrmm_kernel_L2_M8_22:
	KERNEL8x2_SUB
	KERNEL8x2_SUB
	KERNEL8x2_SUB
	KERNEL8x2_SUB

	KERNEL8x2_SUB
	KERNEL8x2_SUB
	KERNEL8x2_SUB
	KERNEL8x2_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L2_M8_22


.Lstrmm_kernel_L2_M8_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L2_M8_100

.Lstrmm_kernel_L2_M8_42:

	KERNEL8x2_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L2_M8_42

.Lstrmm_kernel_L2_M8_100:

	SAVE8x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #8
#else
	sub	tempK, tempK, #2
#endif
	lsl	temp, tempK, #5
	add	pA, pA, temp
	lsl	temp, tempK, #3
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #8
#endif

.Lstrmm_kernel_L2_M8_END:

//------------------------------------------------------------------------------

.Lstrmm_kernel_L2_M4_BEGIN:
	mov	counterI, origM
	tst	counterI , #7
	ble	.Lstrmm_kernel_L2_END

	tst	counterI, #4
	ble	.Lstrmm_kernel_L2_M2_BEGIN

.Lstrmm_kernel_L2_M4_20:

	INIT4x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #3
	add	pB, pB, temp
	lsl	temp, tempOffset, #4
	add	pA, pA, temp
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #4
#else
	add	tempK, tempOffset, #2
#endif
	asr	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL,#0
	ble	.Lstrmm_kernel_L2_M4_40
	.align 5

.Lstrmm_kernel_L2_M4_22:
	KERNEL4x2_SUB
	KERNEL4x2_SUB
	KERNEL4x2_SUB
	KERNEL4x2_SUB

	KERNEL4x2_SUB
	KERNEL4x2_SUB
	KERNEL4x2_SUB
	KERNEL4x2_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L2_M4_22


.Lstrmm_kernel_L2_M4_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L2_M4_100

.Lstrmm_kernel_L2_M4_42:

	KERNEL4x2_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L2_M4_42

.Lstrmm_kernel_L2_M4_100:

	SAVE4x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #4
#else
	sub	tempK, tempK, #2
#endif
	lsl	temp, tempK, #4
	add	pA, pA, temp
	lsl	temp, tempK, #3
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #4
#endif
.Lstrmm_kernel_L2_M4_END:

//------------------------------------------------------------------------------


.Lstrmm_kernel_L2_M2_BEGIN:

	mov	counterI, origM
	tst	counterI , #3
	ble	.Lstrmm_kernel_L2_END

	tst	counterI, #2			// counterI = counterI / 2
	ble	.Lstrmm_kernel_L2_M1_BEGIN

.Lstrmm_kernel_L2_M2_20:

	INIT2x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #3
	add	pB, pB, temp
	lsl	temp, tempOffset, #3
	add	pA, pA, temp
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #2
#else
	add	tempK, tempOffset, #2
#endif
	asr	counterL , tempK, #3		// counterL = counterL / 8
        cmp	counterL,#0
	ble	.Lstrmm_kernel_L2_M2_40

.Lstrmm_kernel_L2_M2_22:

	KERNEL2x2_SUB
	KERNEL2x2_SUB
	KERNEL2x2_SUB
	KERNEL2x2_SUB

	KERNEL2x2_SUB
	KERNEL2x2_SUB
	KERNEL2x2_SUB
	KERNEL2x2_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L2_M2_22


.Lstrmm_kernel_L2_M2_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L2_M2_100

.Lstrmm_kernel_L2_M2_42:

	KERNEL2x2_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L2_M2_42

.Lstrmm_kernel_L2_M2_100:

	SAVE2x2
#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #2
#else
	sub	tempK, tempK, #2
#endif
	lsl	temp, tempK, #3
	add	pA, pA, temp
	lsl	temp, tempK, #3
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #2
#endif

.Lstrmm_kernel_L2_M2_END:


.Lstrmm_kernel_L2_M1_BEGIN:

	tst	counterI, #1			// counterI = counterI % 2
	ble	.Lstrmm_kernel_L2_END

.Lstrmm_kernel_L2_M1_20:

	INIT1x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #3
	add	pB, pB, temp
	lsl	temp, tempOffset, #2
	add	pA, pA, temp
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #1
#else
	add	tempK, tempOffset, #2
#endif
	asr 	counterL , tempK, #3		// counterL = counterL / 8
        cmp     counterL, #0
	ble	.Lstrmm_kernel_L2_M1_40

.Lstrmm_kernel_L2_M1_22:
	KERNEL1x2_SUB
	KERNEL1x2_SUB
	KERNEL1x2_SUB
	KERNEL1x2_SUB

	KERNEL1x2_SUB
	KERNEL1x2_SUB
	KERNEL1x2_SUB
	KERNEL1x2_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L2_M1_22


.Lstrmm_kernel_L2_M1_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L2_M1_100

.Lstrmm_kernel_L2_M1_42:

	KERNEL1x2_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L2_M1_42

.Lstrmm_kernel_L2_M1_100:

	SAVE1x2

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #1
#else
	sub	tempK, tempK, #2
#endif
	lsl	temp, tempK, #2
	add	pA, pA, temp
	lsl	temp, tempK, #3
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #1
#endif
.Lstrmm_kernel_L2_END:
#if !defined(LEFT)
	add	tempOffset, tempOffset, #2
#endif
	add	origPB, origPB, origK, lsl #3	// B = B + K * 2 * 4

/******************************************************************************/

.Lstrmm_kernel_L1_BEGIN:

	mov	counterJ , origN
	tst	counterJ , #1
	ble	.Lstrmm_kernel_L999 // done


	mov	pCRow0, pC			// pCRow0 = C
	add	pC , pC , LDC			// Update pC to point to next

#if defined(LEFT)
	mov	tempOffset, offset
#endif
	mov	pA, origPA			// pA = A

.Lstrmm_kernel_L1_M16_BEGIN:

	mov	counterI, origM
	asr 	counterI, counterI, #4		// counterI = counterI / 16
	cmp	counterI, #0
	ble	.Lstrmm_kernel_L1_M8_BEGIN

.Lstrmm_kernel_L1_M16_20:

	INIT16x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #6
	add	pA, pA, temp
	lsl	temp, tempOffset, #2
	add	pB, pB, temp
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #16
#else
	add	tempK, tempOffset, #1
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Lstrmm_kernel_L1_M16_40
	.align 5

.Lstrmm_kernel_L1_M16_22:
	KERNEL16x1_SUB
	KERNEL16x1_SUB
	KERNEL16x1_SUB
	KERNEL16x1_SUB

	KERNEL16x1_SUB
	KERNEL16x1_SUB
	KERNEL16x1_SUB
	KERNEL16x1_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L1_M16_22


.Lstrmm_kernel_L1_M16_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L1_M16_100

.Lstrmm_kernel_L1_M16_42:

	KERNEL16x1_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L1_M16_42

.Lstrmm_kernel_L1_M16_100:

	SAVE16x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #16
#else
	sub	tempK, tempK, #1
#endif
	lsl	temp, tempK, #6
	add	pA, pA, temp
	lsl	temp, tempK, #2
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #16
#endif

.Lstrmm_kernel_L1_M16_END:

	subs	counterI, counterI, #1
	bgt	.Lstrmm_kernel_L1_M16_20

//------------------------------------------------------------------------------

.Lstrmm_kernel_L1_M8_BEGIN:

	mov	counterI, origM
	tst	counterI , #15
	ble	.Lstrmm_kernel_L1_END

	tst	counterI, #8
	ble	.Lstrmm_kernel_L1_M4_BEGIN

.Lstrmm_kernel_L1_M8_20:

	INIT8x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #5
	add	pA, pA, temp
	lsl	temp, tempOffset, #2
	add	pB, pB, temp
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #8
#else
	add	tempK, tempOffset, #1
#endif

	asr	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Lstrmm_kernel_L1_M8_40
	.align 5

.Lstrmm_kernel_L1_M8_22:
	KERNEL8x1_SUB
	KERNEL8x1_SUB
	KERNEL8x1_SUB
	KERNEL8x1_SUB

	KERNEL8x1_SUB
	KERNEL8x1_SUB
	KERNEL8x1_SUB
	KERNEL8x1_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L1_M8_22


.Lstrmm_kernel_L1_M8_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L1_M8_100

.Lstrmm_kernel_L1_M8_42:

	KERNEL8x1_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L1_M8_42

.Lstrmm_kernel_L1_M8_100:

	SAVE8x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #8
#else
	sub	tempK, tempK, #1
#endif
	lsl	temp, tempK, #5
	add	pA, pA, temp
	lsl	temp, tempK, #2
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #8
#endif

.Lstrmm_kernel_L1_M8_END:

//------------------------------------------------------------------------------

.Lstrmm_kernel_L1_M4_BEGIN:
	mov	counterI, origM
	tst	counterI , #7
	ble	.Lstrmm_kernel_L1_END

	tst	counterI, #4
	ble	.Lstrmm_kernel_L1_M2_BEGIN

.Lstrmm_kernel_L1_M4_20:

	INIT4x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #2
	add	pB, pB, temp
	lsl	temp, tempOffset, #4
	add	pA, pA, temp
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #4
#else
	add	tempK, tempOffset, #1
#endif
	asr	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Lstrmm_kernel_L1_M4_40
	.align 5

.Lstrmm_kernel_L1_M4_22:
	KERNEL4x1_SUB
	KERNEL4x1_SUB
	KERNEL4x1_SUB
	KERNEL4x1_SUB

	KERNEL4x1_SUB
	KERNEL4x1_SUB
	KERNEL4x1_SUB
	KERNEL4x1_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L1_M4_22


.Lstrmm_kernel_L1_M4_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L1_M4_100

.Lstrmm_kernel_L1_M4_42:

	KERNEL4x1_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L1_M4_42

.Lstrmm_kernel_L1_M4_100:

	SAVE4x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #4
#else
	sub	tempK, tempK, #1
#endif
	lsl	temp, tempK, #4
	add	pA, pA, temp
	lsl	temp, tempK, #2
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #4
#endif
.Lstrmm_kernel_L1_M4_END:

//------------------------------------------------------------------------------

.Lstrmm_kernel_L1_M2_BEGIN:

	mov	counterI, origM
	tst	counterI , #3
	ble	.Lstrmm_kernel_L1_END

	tst	counterI, #2			// counterI = counterI / 2
	ble	.Lstrmm_kernel_L1_M1_BEGIN

.Lstrmm_kernel_L1_M2_20:

	INIT2x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #2
	add	pB, pB, temp
	lsl	temp, tempOffset, #3
	add	pA, pA, temp
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #2
#else
	add	tempK, tempOffset, #1
#endif
	asr 	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Lstrmm_kernel_L1_M2_40

.Lstrmm_kernel_L1_M2_22:

	KERNEL2x1_SUB
	KERNEL2x1_SUB
	KERNEL2x1_SUB
	KERNEL2x1_SUB

	KERNEL2x1_SUB
	KERNEL2x1_SUB
	KERNEL2x1_SUB
	KERNEL2x1_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L1_M2_22


.Lstrmm_kernel_L1_M2_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L1_M2_100

.Lstrmm_kernel_L1_M2_42:

	KERNEL2x1_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L1_M2_42

.Lstrmm_kernel_L1_M2_100:

	SAVE2x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
#if defined(LEFT)
	sub	tempK, tempK, #2
#else
	sub	tempK, tempK, #1
#endif
	lsl	temp, tempK, #3
	add	pA, pA, temp
	lsl	temp, tempK, #2
	add	pB, pB, temp
#endif
#if defined(LEFT)
	add	tempOffset, tempOffset, #2
#endif
.Lstrmm_kernel_L1_M2_END:


.Lstrmm_kernel_L1_M1_BEGIN:

	tst	counterI, #1			// counterI = counterI % 2
	ble	.Lstrmm_kernel_L1_END

.Lstrmm_kernel_L1_M1_20:

	INIT1x1

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
#else
	mov	pB, origPB
	lsl	temp, tempOffset, #2
	add	pB, pB, temp
	lsl	temp, tempOffset, #2
	add	pA, pA, temp
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
#elif defined(LEFT)
	add	tempK, tempOffset, #1
#else
	add	tempK, tempOffset, #1
#endif
	asr 	counterL , tempK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Lstrmm_kernel_L1_M1_40

.Lstrmm_kernel_L1_M1_22:
	KERNEL1x1_SUB
	KERNEL1x1_SUB
	KERNEL1x1_SUB
	KERNEL1x1_SUB

	KERNEL1x1_SUB
	KERNEL1x1_SUB
	KERNEL1x1_SUB
	KERNEL1x1_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L1_M1_22


.Lstrmm_kernel_L1_M1_40:

	ands	counterL , tempK, #7		// counterL = counterL % 8
	ble	.Lstrmm_kernel_L1_M1_100

.Lstrmm_kernel_L1_M1_42:

	KERNEL1x1_SUB

	subs	counterL, counterL, #1
	bgt	.Lstrmm_kernel_L1_M1_42

.Lstrmm_kernel_L1_M1_100:

	SAVE1x1

.Lstrmm_kernel_L1_END:

.Lstrmm_kernel_L999:
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

