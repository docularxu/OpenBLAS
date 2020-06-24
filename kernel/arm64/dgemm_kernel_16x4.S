/*******************************************************************************
Copyright (c) 2015, The OpenBLAS Project
Copyright (c) 2020, Hisilicon Limited
All rights reserved.

Based on the original dgemm_kernel_8x4.S and modified by SVE instructions
Created by Jia Yuan <yuanjia11@huawei.com>
           Anjun Wu <wuanjun@huawei.com>

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

/*                   X0          X1          X2          s0         X3        x4       x5           x6 */
/*int CNAME(BLASLONG bm,BLASLONG bn,BLASLONG bk,FLOAT alpha0,FLOAT* ba,FLOAT* bb,FLOAT* C,BLASLONG ldc )*/

#define origM           x0
#define origN           x1
#define origK           x2
#define origPA          x3
#define origPB          x4
#define pC              x5
#define LDC             x6
#define temp            x7
#define counterL        x8
#define counterI        x9
#define counterJ        x10
#define pB              x11
#define pCRow0          x12
#define pCRow1          x13
#define pCRow2          x14
#define pCRow3          x15
#define pA              x16
#define alpha           x17

#define alpha0          z10.d
#define alphaV0         z10.d[0]

#define A_PRE_SIZE      1024
#define B_PRE_SIZE      1024
#define C_PRE_SIZE      1024

// 00 origM
// 01 origN
// 02 origK
// 03 origPA
// 04 origPB
// 05 pC
// 06 origLDC -> LDC
// 07 temp
// 08 counterL
// 09 counterI
// 10 counterJ
// 11 pB
// 12 pCRow0
// 13 pCRow1
// 14 pCRow2
// 15 pCRow3
// 16 pA
// 17
// 18 must save
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

.macro INIT16x4
        fmov    d16, xzr
        fmov    d17, xzr
        fmov    d18, d16
        fmov    d19, xzr
        fmov    d20, xzr
        fmov    d21, d16
        fmov    d22, d17
        fmov    d23, d18
        fmov    d24, xzr
        fmov    d25, d16
        fmov    d26, d17
        fmov    d27, d18
        fmov    d28, xzr
        fmov    d29, d16
        fmov    d30, d17
        fmov    d31, d18
.endm

.macro KERNEL16x4_I
        ptrue   p0.d, all
        ld1d    {z0.d}, p0/z, [pA]
        ld1rqd  z8.d, p0/z, [pB]

        fmul    z16.d, z0.d, z8.d[0]
        fmul    z20.d, z0.d, z8.d[1]

        ld1d    {z1.d}, p0/z, [pA, #1, MUL VL]
        ld1rqd  z9.d, p0/z, [pB, #16]

        fmul    z17.d, z1.d, z8.d[0]
        fmul    z21.d, z1.d, z8.d[1]
        fmul    z24.d, z0.d, z9.d[0]
        fmul    z28.d, z0.d, z9.d[1]

        ld1d    {z2.d}, p0/z, [pA, #2, MUL VL]
        ld1d    {z3.d}, p0/z, [pA, #3, MUL VL]

        fmul    z25.d, z1.d, z9.d[0]
        fmul    z29.d, z1.d, z9.d[1]

        ld1d    {z4.d}, p0/z, [pA, #4, MUL VL]
        ld1d    {z5.d}, p0/z, [pA, #5, MUL VL]

        fmul    z18.d, z2.d, z8.d[0]
        fmul    z22.d, z2.d, z8.d[1]

        ld1d    {z6.d}, p0/z, [pA, #6, MUL VL]
        ld1d    {z7.d}, p0/z, [pA, #7, MUL VL]

        fmul    z27.d, z3.d, z9.d[0]
        fmul    z31.d, z3.d, z9.d[1]
        fmul    z26.d, z2.d, z9.d[0]
        fmul    z30.d, z2.d, z9.d[1]

        ld1rqd  z10.d, p0/z, [pB, #32]
        ld1rqd  z11.d, p0/z, [pB, #48]

        fmul    z19.d, z3.d, z8.d[0]
        fmul    z23.d, z3.d, z8.d[1]

        incb    pA, all, MUL #8
        add     pB, pB, #64
.endm

.macro KERNEL16x4_M1
        ptrue   p0.d, all

        fmla    z16.d, z0.d, z8.d[0]
        fmla    z20.d, z0.d, z8.d[1]
        fmla    z24.d, z0.d, z9.d[0]
        fmla    z28.d, z0.d, z9.d[1]

        ld1d    {z4.d}, p0/z, [pA]
        ld1d    {z5.d}, p0/z, [pA, #1, MUL VL]

        fmla    z17.d, z1.d, z8.d[0]
        fmla    z25.d, z1.d, z9.d[0]
        fmla    z21.d, z1.d, z8.d[1]
        fmla    z29.d, z1.d, z9.d[1]

        ld1d    {z6.d}, p0/z, [pA, #2, MUL VL]
        ld1d    {z7.d}, p0/z, [pA, #3, MUL VL]

        fmla    z18.d, z2.d, z8.d[0]
        fmla    z22.d, z2.d, z8.d[1]
        fmla    z26.d, z2.d, z9.d[0]
        fmla    z30.d, z2.d, z9.d[1]

        ld1rqd  z10.d, p0/z, [pB]
        ld1rqd  z11.d, p0/z, [pB, #16]

        fmla    z19.d, z3.d, z8.d[0]
        fmla    z23.d, z3.d, z8.d[1]
        fmla    z27.d, z3.d, z9.d[0]
        fmla    z31.d, z3.d, z9.d[1]

        incb    pA, all, MUL #4
        add     pB, pB, #32
.endm

.macro KERNEL16x4_M2
        ptrue   p0.d, all

        fmla    z16.d, z4.d, z10.d[0]
        fmla    z20.d, z4.d, z10.d[1]
        fmla    z24.d, z4.d, z11.d[0]
        fmla    z28.d, z4.d, z11.d[1]

        ld1d    {z0.d}, p0/z, [pA]
        ld1d    {z1.d}, p0/z, [pA, #1, MUL VL]

        fmla    z17.d, z5.d, z10.d[0]
        fmla    z25.d, z5.d, z11.d[0]
        fmla    z21.d, z5.d, z10.d[1]
        fmla    z29.d, z5.d, z11.d[1]

        ld1d    {z2.d}, p0/z, [pA, #2, MUL VL]
        ld1d    {z3.d}, p0/z, [pA, #3, MUL VL]

        fmla    z18.d, z6.d, z10.d[0]
        fmla    z22.d, z6.d, z10.d[1]
        fmla    z26.d, z6.d, z11.d[0]
        fmla    z30.d, z6.d, z11.d[1]

        ld1rqd  z8.d, p0/z, [pB]
        ld1rqd  z9.d, p0/z, [pB, #16]

        fmla    z19.d, z7.d, z10.d[0]
        fmla    z23.d, z7.d, z10.d[1]
        fmla    z27.d, z7.d, z11.d[0]
        fmla    z31.d, z7.d, z11.d[1]

        incb    pA, all, MUL #4
        add     pB, pB, #32
.endm

.macro KERNEL16x4_E
        fmla    z16.d, z4.d, z10.d[0]
        fmla    z20.d, z4.d, z10.d[1]
        fmla    z24.d, z4.d, z11.d[0]
        fmla    z28.d, z4.d, z11.d[1]

        fmla    z17.d, z5.d, z10.d[0]
        fmla    z25.d, z5.d, z11.d[0]
        fmla    z21.d, z5.d, z10.d[1]
        fmla    z29.d, z5.d, z11.d[1]

        fmla    z18.d, z6.d, z10.d[0]
        fmla    z22.d, z6.d, z10.d[1]
        fmla    z26.d, z6.d, z11.d[0]
        fmla    z30.d, z6.d, z11.d[1]

        fmla    z19.d, z7.d, z10.d[0]
        fmla    z23.d, z7.d, z10.d[1]
        fmla    z27.d, z7.d, z11.d[0]
        fmla    z31.d, z7.d, z11.d[1]
.endm

.macro KERNEL16x4_SUB
        ptrue   p0.d, all
        ld1d    {z0.d}, p0/z, [pA]
        ld1rqd  z8.d, p0/z, [pB]

        fmla    z16.d, z0.d, z8.d[0]
        fmla    z20.d, z0.d, z8.d[1]

        ld1d    {z1.d}, p0/z, [pA, #1, MUL VL]
        ld1rqd  z9.d, p0/z, [pB, #16]

        fmla    z17.d, z1.d, z8.d[0]
        fmla    z21.d, z1.d, z8.d[1]

        ld1d    {z2.d}, p0/z, [pA, #2, MUL VL]

        fmla    z24.d, z0.d, z9.d[0]
        fmla    z28.d, z0.d, z9.d[1]
        fmla    z25.d, z1.d, z9.d[0]
        fmla    z29.d, z1.d, z9.d[1]

        ld1d    {z3.d}, p0/z, [pA, #3, MUL VL]

        fmla    z18.d, z2.d, z8.d[0]
        fmla    z22.d, z2.d, z8.d[1]
        fmla    z26.d, z2.d, z9.d[0]
        fmla    z30.d, z2.d, z9.d[1]

        fmla    z19.d, z3.d, z8.d[0]
        fmla    z23.d, z3.d, z8.d[1]
        fmla    z27.d, z3.d, z9.d[0]
        fmla    z31.d, z3.d, z9.d[1]

        incb    pA, all, MUL #4
        add     pB, pB, #32
.endm

.macro SAVE16x4
        ptrue   p0.d, all
        dup     alpha0, alpha

        ld1d    {z0.d}, p0/z, [pCRow0]
        ld1d    {z1.d}, p0/z, [pCRow0, #1, MUL VL]

        fmla    z0.d, z16.d, alphaV0
        fmla    z1.d, z17.d, alphaV0

        st1d    {z0.d}, p0, [pCRow0]
        st1d    {z1.d}, p0, [pCRow0, #1, MUL VL]

        ld1d    {z2.d}, p0/z, [pCRow0, #2, MUL VL]
        ld1d    {z3.d}, p0/z, [pCRow0, #3, MUL VL]

        fmla    z2.d, z18.d, alphaV0
        fmla    z3.d, z19.d, alphaV0

        st1d    {z2.d}, p0, [pCRow0, #2, MUL VL]
        st1d    {z3.d}, p0, [pCRow0, #3, MUL VL]

        incb    pCRow0, all, MUL #4

        ld1d    {z4.d}, p0/z, [pCRow1]
        ld1d    {z5.d}, p0/z, [pCRow1 , #1, MUL VL]

        fmla    z4.d, z20.d, alphaV0
        fmla    z5.d, z21.d, alphaV0

        st1d    {z4.d}, p0, [pCRow1]
        st1d    {z5.d}, p0, [pCRow1, #1, MUL VL]

        ld1d    {z6.d}, p0/z, [pCRow1, #2, MUL VL]
        ld1d    {z7.d}, p0/z, [pCRow1, #3, MUL VL]

        fmla    z6.d, z22.d, alphaV0
        fmla    z7.d, z23.d, alphaV0

        st1d    {z6.d}, p0, [pCRow1, #2, MUL VL]
        st1d    {z7.d}, p0, [pCRow1, #3, MUL VL]

        incb    pCRow1, all, MUL #4

        ld1d    {z0.d}, p0/z, [pCRow2]
        ld1d    {z1.d}, p0/z, [pCRow2, #1, MUL VL]

        fmla    z0.d, z24.d, alphaV0
        fmla    z1.d, z25.d, alphaV0

        st1d    {z0.d}, p0, [pCRow2]
        st1d    {z1.d}, p0, [pCRow2, #1, MUL VL]

        ld1d    {z2.d}, p0/z, [pCRow2, #2, MUL VL]
        ld1d    {z3.d}, p0/z, [pCRow2, #3, MUL VL]

        fmla    z2.d, z26.d, alphaV0
        fmla    z3.d, z27.d, alphaV0

        st1d    {z2.d}, p0, [pCRow2, #2, MUL VL]
        st1d    {z3.d}, p0, [pCRow2, #3, MUL VL]

        incb    pCRow2, all, MUL #4

        ld1d    {z4.d}, p0/z, [pCRow3]
        ld1d    {z5.d}, p0/z, [pCRow3, #1, MUL VL]

        fmla    z4.d, z28.d, alphaV0
        fmla    z5.d, z29.d, alphaV0

        st1d    {z4.d}, p0, [pCRow3]
        st1d    {z5.d}, p0, [pCRow3, #1, MUL VL]

        ld1d    {z6.d}, p0/z, [pCRow3, #2, MUL VL]
        ld1d    {z7.d}, p0/z, [pCRow3, #3, MUL VL]

        fmla    z6.d, z30.d, alphaV0
        fmla    z7.d, z31.d, alphaV0

        st1d    {z6.d}, p0, [pCRow3, #2, MUL VL]
        st1d    {z7.d}, p0, [pCRow3, #3, MUL VL]

        incb    pCRow3, all, MUL #4
.endm

/******************************************************************************/

.macro INIT8x4
        fmov    d16, xzr
        fmov    d17, d16
        fmov    d20, d17
        fmov    d21, d16
        fmov    d24, d17
        fmov    d25, d16
        fmov    d28, d17
        fmov    d29, d16
.endm

.macro KERNEL8x4_SUB
        ptrue   p0.d, all
        ld1rqd  {z8.d}, p0/z, [pB]
        ld1rqd  {z9.d}, p0/z, [pB, #16]
        add     pB, pB, #32

        ld1d    {z0.d}, p0/z, [pA]
        ld1d    {z1.d}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        fmla    z16.d, z0.d, z8.d[0]
        fmla    z29.d, z1.d, z9.d[1]

        fmla    z20.d, z0.d, z8.d[1]
        fmla    z25.d, z1.d, z9.d[0]

        fmla    z24.d, z0.d, z9.d[0]
        fmla    z21.d, z1.d, z8.d[1]

        fmla    z28.d, z0.d, z9.d[1]
        fmla    z17.d, z1.d, z8.d[0]
.endm

.macro SAVE8x4
        ptrue  p0.d, all
        dup    alpha0, alpha

        ld1d    {z8.d}, p0/z, [pCRow0]
        ld1d    {z9.d}, p0/z, [pCRow0, #1, MUL VL]

        fmla    z8.d, z16.d, alphaV0
        fmla    z9.d, z17.d, alphaV0

        st1d    {z8.d}, p0, [pCRow0]
        st1d    {z9.d}, p0, [pCRow0, #1, MUL VL]

        incb    pCRow0, all, MUL #2

        ld1d    {z12.d}, p0/z, [pCRow1]
        ld1d    {z13.d}, p0/z, [pCRow1, #1, MUL VL]
        fmla    z12.d, z20.d, alphaV0
        fmla    z13.d, z21.d, alphaV0

        st1d    {z12.d}, p0, [pCRow1]
        st1d    {z13.d}, p0, [pCRow1, #1, MUL VL]

        incb    pCRow1, all, MUL #2

        ld1d    {z8.d}, p0/z, [pCRow2]
        ld1d    {z9.d}, p0/z, [pCRow2, #1, MUL VL]
        fmla    z8.d, z24.d, alphaV0
        fmla    z9.d, z25.d, alphaV0

        st1d    {z8.d}, p0, [pCRow2]
        st1d    {z9.d}, p0, [pCRow2, #1, MUL VL]

        incb    pCRow2, all, MUL #2

        ld1d    {z12.d}, p0/z, [pCRow3]
        ld1d    {z13.d}, p0/z, [pCRow3, #1, MUL VL]
        fmla    z12.d, z28.d, alphaV0
        fmla    z13.d, z29.d, alphaV0

        st1d    {z12.d}, p0, [pCRow3]
        st1d    {z13.d}, p0, [pCRow3, #1, MUL VL]

        incb    pCRow3, all, MUL #2
.endm

/******************************************************************************/

.macro INIT4x4
        fmov    d16, xzr
        fmov    d20, d16
        fmov    d24, d20
        fmov    d28, d16
.endm

.macro KERNEL4x4_SUB
        ptrue    p0.d, all

        ld1rqd  {z8.d}, p0/z, [pB]
        ld1rqd  {z9.d}, p0/z, [pB, #16]
        add     pB, pB, #32

        ld1d    {z0.d}, p0/z, [pA]
        add     pA, pA, #32

        fmla    z16.d, z0.d, z8.d[0]
        fmla    z20.d, z0.d, z8.d[1]
        fmla    z24.d, z0.d, z9.d[0]
        fmla    z28.d, z0.d, z9.d[1]
.endm

.macro SAVE4x4
        ptrue   p0.d, all
        dup     alpha0, alpha

        ld1d    {z8.d}, p0/z, [pCRow0]
        ld1d    {z12.d}, p0/z, [pCRow1]

        fmla    z8.d, z16.d, alphaV0
        fmla    z12.d, z20.d, alphaV0

        st1d    {z8.d}, p0, [pCRow0]
        st1d    {z12.d}, p0, [pCRow1]

        add     pCRow0, pCRow0, #32
        add     pCRow1, pCRow1, #32

        ld1d    {z8.d}, p0/z, [pCRow2]
        ld1d    {z12.d}, p0/z, [pCRow3]

        fmla    z8.d, z24.d, alphaV0
        fmla    z12.d, z28.d, alphaV0

        st1d    {z8.d}, p0, [pCRow2]
        st1d    {z12.d}, p0, [pCRow3]

        add     pCRow2, pCRow2, #32
        add     pCRow3, pCRow3, #32
.endm


/******************************************************************************/

.macro INIT2x4
        fmov    d16, xzr
        fmov    d20, d16
        fmov    d24, d20
        fmov    d28, d16
.endm

.macro KERNEL2x4_SUB
        ptrue   p0.d, VL2
        ld1rqd  z0.d, p0/z, [pA]
        add     pA, pA, #16

        ld1rqd  {z8.d}, p0/z, [pB]
        ld1rqd  {z9.d}, p0/z, [pB, #16]
        add     pB, pB, #32

        fmla    z16.d, z0.d, z8.d[0]
        fmla    z20.d, z0.d, z8.d[1]
        fmla    z24.d, z0.d, z9.d[0]
        fmla    z28.d, z0.d, z9.d[1]
.endm

.macro SAVE2x4
        ptrue   p0.d, VL2
        dup     alpha0, alpha

        ld1rqd  z8.d, p0/z, [pCRow0]
        ld1rqd  z12.d, p0/z, [pCRow1]

        fmla    z8.d,  z16.d, alphaV0
        fmla    z12.d,  z20.d, alphaV0

        st1d    {z8.d}, p0, [pCRow0]
        st1d    {z12.d}, p0, [pCRow1]

        add     pCRow0, pCRow0, #16
        add     pCRow1, pCRow1, #16

        ld1rqd  z8.d, p0/z, [pCRow2]
        ld1rqd  z12.d, p0/z, [pCRow3]

        fmla    z8.d, z24.d, alphaV0
        fmla    z12.d, z28.d, alphaV0

        st1d    {z8.d}, p0, [pCRow2]
        st1d    {z12.d}, p0, [pCRow3]

        add     pCRow2, pCRow2, #16
        add     pCRow3, pCRow3, #16
.endm

/******************************************************************************/

.macro INIT1x4
        fmov    d16, xzr
        fmov    d17, d16
        fmov    d18, d17
        fmov    d19, d18
.endm

.macro KERNEL1x4_SUB
       ptrue   p0.d, VL1
       ld1rd   z0.d, p0/z, [pA]
       add     pA, pA, #8

       ld1rd   {z8.d}, p0/z, [pB]
       ld1rd   {z9.d}, p0/z, [pB, #8]
       ld1rd   {z10.d}, p0/z, [pB, #16]
       ld1rd   {z11.d}, p0/z, [pB, #24]
       add     pB, pB, #32

       fmla    z16.d, z8.d, z0.d[0]
       fmla    z17.d, z9.d, z0.d[0]
       fmla    z18.d, z10.d, z0.d[0]
       fmla    z19.d, z11.d, z0.d[0]
.endm

.macro SAVE1x4
        ptrue   p0.d, VL1
        dup     alpha0, alpha

        ld1rd   z8.d, p0/z, [pCRow0]
        ld1rd   z9.d, p0/z, [pCRow1]

        fmla    z8.d,  z16.d, alphaV0
        fmla    z9.d,  z17.d, alphaV0

        st1d    {z8.d}, p0, [pCRow0]
        st1d    {z9.d}, p0, [pCRow1]

        add     pCRow0, pCRow0, #8
        add     pCRow1, pCRow1, #8

        ld1rd   z12.d, p0/z, [pCRow2]
        ld1rd   z13.d, p0/z, [pCRow3]

        fmla    z12.d, z18.d, alphaV0
        fmla    z13.d, z19.d, alphaV0

        st1d    {z12.d}, p0, [pCRow2]
        st1d    {z13.d}, p0, [pCRow3]

        add     pCRow2, pCRow2, #8
        add     pCRow3, pCRow3, #8
.endm

/******************************************************************************/

.macro INIT16x2
        fmov    d16, xzr
        fmov    d17, xzr
        fmov    d18, d16
        fmov    d19, d17
        fmov    d20, xzr
        fmov    d21, d16
        fmov    d22, d17
        fmov    d23, d18
.endm

.macro KERNEL16x2_SUB
        ptrue   p0.d, all
        ld1d    {z0.d}, p0/z, [pA]
        ld1d    {z1.d}, p0/z, [pA, #1, MUL VL]
        ld1d    {z2.d}, p0/z, [pA, #2, MUL VL]
        ld1d    {z3.d}, p0/z, [pA, #3, MUL VL]

        ld1rqd  {z8.d}, p0/z, [pB]
        add     pB, pB, #16
        incb    pA, all, MUL #4

        fmla    z16.d, z0.d, z8.d[0]
        fmla    z17.d, z1.d, z8.d[0]
        fmla    z18.d, z2.d, z8.d[0]
        fmla    z19.d, z3.d, z8.d[0]

        fmla    z20.d, z0.d, z8.d[1]
        fmla    z21.d, z1.d, z8.d[1]
        fmla    z22.d, z2.d, z8.d[1]
        fmla    z23.d, z3.d, z8.d[1]
.endm

.macro SAVE16x2
        ptrue   p0.d, all
        dup     alpha0, alpha

        ld1d    {z0.d}, p0/z, [pCRow0]
        ld1d    {z1.d}, p0/z, [pCRow0, #1, MUL VL]
        fmla    z0.d, z16.d, alphaV0
        fmla    z1.d, z17.d, alphaV0

        st1d    {z0.d}, p0, [pCRow0]
        st1d    {z1.d}, p0, [pCRow0, #1, MUL VL]

        ld1d    {z2.d}, p0/z, [pCRow0, #2, MUL VL]
        ld1d    {z3.d}, p0/z, [pCRow0, #3, MUL VL]
        fmla    z2.d, z18.d, alphaV0
        fmla    z3.d, z19.d, alphaV0

        st1d    {z2.d}, p0, [pCRow0, #2, MUL VL]
        st1d    {z3.d}, p0, [pCRow0, #3, MUL VL]

        incb    pCRow0, all, MUL #4

        ld1d    {z4.d}, p0/z, [pCRow1]
        ld1d    {z5.d}, p0/z, [pCRow1, #1, MUL VL]
        fmla    z4.d, z20.d, alphaV0
        fmla    z5.d, z21.d, alphaV0

        st1d    {z4.d}, p0, [pCRow1]
        st1d    {z5.d}, p0, [pCRow1, #1, MUL VL]

        ld1d    {z6.d}, p0/z, [pCRow1, #2, MUL VL]
        ld1d    {z7.d}, p0/z, [pCRow1, #3, MUL VL]
        fmla    z6.d, z22.d, alphaV0
        fmla    z7.d, z23.d, alphaV0

        st1d    {z6.d}, p0, [pCRow1, #2, MUL VL]
        st1d    {z7.d}, p0, [pCRow1, #3, MUL VL]
        incb    pCRow1, all, MUL #4
.endm

/******************************************************************************/

.macro INIT8x2
        fmov    d16, xzr
        fmov    d17, xzr
        fmov    d20, xzr
        fmov    d21, d16
.endm

.macro KERNEL8x2_SUB
        ptrue   p0.d, all
        ld1d    {z0.d}, p0/z, [pA]
        ld1d    {z1.d}, p0/z, [pA, #1, MUL VL]
        ld1rqd  z8.d, p0/z, [pB]
        incb pA, all, MUL #2
        add     pB, pB, #16

        fmla    z16.d, z0.d, z8.d[0]
        fmla    z17.d, z1.d, z8.d[0]
        fmla    z20.d, z0.d, z8.d[1]
        fmla    z21.d, z1.d, z8.d[1]
.endm

.macro SAVE8x2
        ptrue   p0.d, all
        dup     alpha0, alpha

        ld1d    {z0.d}, p0/z, [pCRow0]
        ld1d    {z1.d}, p0/z, [pCRow0, #1, MUL VL]
        fmla    z0.d, z16.d, alphaV0
        fmla    z1.d, z17.d, alphaV0

        st1d    {z0.d}, p0, [pCRow0]
        st1d    {z1.d}, p0, [pCRow0, #1, MUL VL]

        incb    pCRow0, all, MUL #2

        ld1d    {z2.d}, p0/z, [pCRow1]
        ld1d    {z3.d}, p0/z, [pCRow1, #1, MUL VL]
        fmla    z2.d, z20.d, alphaV0
        fmla    z3.d, z21.d, alphaV0

        st1d    {z2.d}, p0, [pCRow1]
        st1d    {z3.d}, p0, [pCRow1, #1, MUL VL]

        incb    pCRow1, all, MUL #2
.endm

/******************************************************************************/

.macro INIT4x2
        fmov    d16, xzr
        fmov    d17, d16
.endm

.macro KERNEL4x2_SUB
        ptrue   p0.d, all
        ld1rqd  z8.d, p0/z, [pB]
        ld1d    z0.d, p0/z, [pA]
        add     pB, pB, #16
        add     pA, pA, #32

        fmla    z16.d, z0.d, z8.d[0]
        fmla    z17.d, z0.d, z8.d[1]
.endm

.macro SAVE4x2
        ptrue   p0.d, all
        dup     alpha0, alpha

        ld1d    {z8.d}, p0/z, [pCRow0]
        ld1d    {z12.d}, p0/z, [pCRow1]

        fmla    z8.d, z16.d, alphaV0
        fmla    z12.d, z17.d, alphaV0

        st1d    {z8.d}, p0, [pCRow0]
        st1d    {z12.d}, p0, [pCRow1]

        add     pCRow0, pCRow0, #32
        add     pCRow1, pCRow1, #32
.endm
/******************************************************************************/

.macro INIT2x2
        fmov    d16, xzr
        fmov    d20, d16
.endm

.macro KERNEL2x2_SUB
        ptrue   p0.d, VL2
        ld1rqd  z8.d, p0/z, [pB]
        add     pB, pB, #16

        ld1rqd  z0.d, p0/z, [pA]
        add     pA, pA, #16

        fmla    z16.d, z0.d, z8.d[0]
        fmla    z20.d, z0.d, z8.d[1]
.endm

.macro SAVE2x2
        ptrue   p0.d, VL2
        dup     alpha0, alpha

        ld1rqd  z8.d, p0/z, [pCRow0]
        ld1rqd  z12.d, p0/z, [pCRow1]

        fmla    z8.d, z16.d, alphaV0
        fmla    z12.d, z20.d, alphaV0

        st1d    {z8.d}, p0, [pCRow0]
        st1d    {z12.d}, p0, [pCRow1]

        add     pCRow0, pCRow0, #16
        add     pCRow1, pCRow1, #16
.endm

/******************************************************************************/

.macro INIT1x2
        fmov    d16, xzr
        fmov    d17, xzr
.endm

.macro KERNEL1x2_SUB
        ptrue   p0.d, VL1
        ld1rd   z8.d, p0/z, [pB]
        ld1rd   z9.d, p0/z, [pB, #8]
        add     pB , pB, #16

        ld1rd   z0.d, p0/z, [pA]
        add     pA, pA, #8

        fmla    z16.d, z8.d, z0.d[0]
        fmla    z17.d, z9.d, z0.d[0]
.endm

.macro SAVE1x2
        ptrue   p0.d, VL1
        dup     alpha0, alpha

        ld1rd   z8.d, p0/z, [pCRow0]
        ld1rd   z9.d, p0/z, [pCRow1]

        fmla    z8.d, z16.d, alphaV0
        fmla    z9.d, z17.d, alphaV0

        st1d    {z8.d}, p0, [pCRow0]
        st1d    {z9.d}, p0, [pCRow1]

        add     pCRow0, pCRow0, #8
        add     pCRow1, pCRow1, #8
.endm

/******************************************************************************/

.macro INIT16x1
        fmov    d16, xzr
        fmov    d17, xzr
        fmov    d18, d16
        fmov    d19, d17
.endm

.macro KERNEL16x1_SUB
        ptrue   p0.d, all
        ld1d    {z0.d}, p0/z, [pA]
        ld1d    {z1.d}, p0/z, [pA, #1, MUL VL]
        ld1d    {z2.d}, p0/z, [pA, #2, MUL VL]
        ld1d    {z3.d}, p0/z, [pA, #3, MUL VL]

        ld1rd   z8.d, p0/z, [pB]

        incb    pA, all, MUL #4
        add     pB , pB, #8

        fmla    z16.d, z0.d, z8.d[0]
        fmla    z17.d, z1.d, z8.d[0]
        fmla    z18.d, z2.d, z8.d[0]
        fmla    z19.d, z3.d, z8.d[0]
.endm

.macro SAVE16x1
        ptrue   p0.d, all
        dup     alpha0, alpha

        ld1d    {z0.d}, p0/z, [pCRow0]
        ld1d    {z1.d}, p0/z, [pCRow0, #1, MUL VL]
        ld1d    {z2.d}, p0/z, [pCRow0, #2, MUL VL]
        ld1d    {z3.d}, p0/z, [pCRow0, #3, MUL VL]

        fmla    z0.d, z16.d, alphaV0
        fmla    z1.d, z17.d, alphaV0
        fmla    z2.d, z18.d, alphaV0
        fmla    z3.d, z19.d, alphaV0

        st1d    {z0.d}, p0, [pCRow0]
        st1d    {z1.d}, p0, [pCRow0, #1, MUL VL]
        st1d    {z2.d}, p0, [pCRow0, #2, MUL VL]
        st1d    {z3.d}, p0, [pCRow0, #3, MUL VL]

        incb    pCRow0, all, MUL #4
.endm

/******************************************************************************/

.macro INIT8x1
        fmov    d16, xzr
        fmov    d17, xzr
.endm

.macro KERNEL8x1_SUB
        ptrue   p0.d, all
        ld1d    {z0.d}, p0/z, [pA]
        ld1d    {z1.d}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        ld1rd   z8.d, p0/z, [pB]
        add     pB, pB, #8

        fmla    z16.d, z0.d, z8.d[0]
        fmla    z17.d, z1.d, z8.d[0]
.endm

.macro SAVE8x1
        ptrue   p0.d, all
        dup     alpha0, alpha

        ld1d    {z0.d}, p0/z, [pCRow0]
        ld1d    {z1.d}, p0/z, [pCRow0, #1, MUL VL]

        fmla    z0.d, z16.d, alphaV0
        fmla    z1.d, z17.d, alphaV0

        st1d    {z0.d}, p0, [pCRow0]
        st1d    {z1.d}, p0, [pCRow0, #1, MUL VL]

        incb    pCRow0, all, MUL #2
.endm


/******************************************************************************/

.macro INIT4x1
        fmov    d16, xzr
.endm

.macro KERNEL4x1_SUB
        ptrue   p0.d, all
        ld1d    {z0.d}, p0/z, [pA]
        ld1rd   z8.d, p0/z, [pB]

        add     pA, pA, #32
        add     pB, pB, #8

        fmla    z16.d, z0.d, z8.d[0]
.endm

.macro SAVE4x1
        ptrue   p0.d, all
        dup     alpha0, alpha

        ld1d    {z8.d}, p0/z, [pCRow0]
        fmla    z8.d, z16.d, alphaV0

        st1d    {z8.d}, p0, [pCRow0]
        add     pCRow0, pCRow0, #32
.endm


/******************************************************************************/

.macro INIT2x1
        fmov    d16, xzr
.endm

.macro KERNEL2x1_SUB
        ptrue   p0.d, VL2
        ld1rqd  z0.d, p0/z, [pA]
        ld1rd   z8.d, p0/z, [pB]

        add     pA , pA, #16
        add     pB , pB, #8

        fmla    z16.d, z0.d, z8.d[0]
.endm

.macro SAVE2x1
        ptrue   p0.d, VL2
        dup     alpha0, alpha

        ld1rqd  z8.d, p0/z, [pCRow0]
        fmla    z8.d,  z16.d, alphaV0

        st1d    {z8.d}, p0, [pCRow0]
        add     pCRow0, pCRow0, #16
.endm

/******************************************************************************/

.macro INIT1x1
        fmov    d16, xzr
.endm

.macro KERNEL1x1_SUB
        ptrue   p0.d, VL1
        ld1rd   z0.d, p0/z, [pA]
        ld1rd   z8.d, p0/z, [pB]

        add     pA, pA, #8
        add     pB, pB, #8

        fmla    z16.d, p0/m, z0.d, z8.d
.endm

.macro SAVE1x1
        ptrue   p0.d, VL1
        dup     alpha0, alpha

        ld1rd   z8.d, p0/z, [pCRow0]
        fmla    z8.d, p0/m, alpha0, z16.d
        st1d    {z8.d}, p0, [pCRow0]

        add     pCRow0, pCRow0, #8
.endm

/*******************************************************************************
* End of macro definitions
*******************************************************************************/

PROLOGUE

        .align 5
        add    sp, sp, #-(11 * 16)
        stp    d8, d9, [sp, #(0 * 16)]
        stp    d10, d11, [sp, #(1 * 16)]
        stp    d12, d13, [sp, #(2 * 16)]
        stp    d14, d15, [sp, #(3 * 16)]
        stp    d16, d17, [sp, #(4 * 16)]
        stp    x18, x19, [sp, #(5 * 16)]
        stp    x20, x21, [sp, #(6 * 16)]
        stp    x22, x23, [sp, #(7 * 16)]
        stp    x24, x25, [sp, #(8 * 16)]
        stp    x26, x27, [sp, #(9 * 16)]
        str    x28, [sp, #(10 * 16)]

        prfm   PLDL1KEEP, [origPB]
        prfm   PLDL1KEEP, [origPA]

        fmov   alpha, d0

        lsl    LDC, LDC, #3            // ldc = ldc * 8

        mov    pB, origPB

        mov    counterJ, origN
        asr    counterJ, counterJ, #2        // J = J / 4
        cmp    counterJ, #0
        ble    .Ldgemm_kernel_L2_BEGIN

/******************************************************************************/

    .align 5
.Ldgemm_kernel_L4_BEGIN:

        mov    pCRow0, pC
        add    pCRow1, pCRow0, LDC
        add    pCRow2, pCRow1, LDC
        add    pCRow3, pCRow2, LDC

        add    pC, pCRow3, LDC

        mov    pA, origPA            // pA = start of A array

.Ldgemm_kernel_L4_M16_BEGIN:

        mov    counterI, origM
        asr    counterI, counterI, #4        // counterI = counterI / 16
        cmp    counterI, #0
        ble    .Ldgemm_kernel_L4_M8_BEGIN

        .align 5
.Ldgemm_kernel_L4_M16_20:

        mov    pB, origPB

        asr    counterL , origK, #3     // L = K / 8
        cmp    counterL , #2            // is there at least 4 to do?
        blt    .Ldgemm_kernel_L4_M16_32

        KERNEL16x4_I
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_M2

        subs    counterL, counterL, #2        // subtract 2
        ble    .Ldgemm_kernel_L4_M16_22a

        .align 5
.Ldgemm_kernel_L4_M16_22:

        KERNEL16x4_M1
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_M2

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L4_M16_22

        .align 5
.Ldgemm_kernel_L4_M16_22a:

        KERNEL16x4_M1
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_E

        b     .Ldgemm_kernel_L4_M16_44

        .align 5
.Ldgemm_kernel_L4_M16_32:

        tst    counterL, #1
        ble    .Ldgemm_kernel_L4_M16_40

        KERNEL16x4_I
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_E

        b    .Ldgemm_kernel_L4_M16_44

.Ldgemm_kernel_L4_M16_40:

        INIT16x4

.Ldgemm_kernel_L4_M16_44:

        ands    counterL , origK, #7
        ble    .Ldgemm_kernel_L4_M16_100

        .align 5
.Ldgemm_kernel_L4_M16_46:

        KERNEL16x4_SUB

        subs    counterL, counterL, #1
        bne    .Ldgemm_kernel_L4_M16_46

.Ldgemm_kernel_L4_M16_100:

        prfm    PLDL1KEEP, [pA]
        prfm    PLDL1KEEP, [pA, #64]
        prfm    PLDL1KEEP, [origPB]

        SAVE16x4

.Ldgemm_kernel_L4_M16_END:

        subs    counterI, counterI, #1
        bne    .Ldgemm_kernel_L4_M16_20

.Ldgemm_kernel_L4_M8_BEGIN:

        mov    counterI, origM
        tst    counterI, #15
        ble    .Ldgemm_kernel_L4_END

        tst    counterI, #8
        ble    .Ldgemm_kernel_L4_M4_BEGIN

.Ldgemm_kernel_L4_M8_20:

        INIT8x4

        mov    pB, origPB

        asr    counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL , #0
        ble    .Ldgemm_kernel_L4_M8_40

        .align 5
.Ldgemm_kernel_L4_M8_22:

        KERNEL8x4_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL8x4_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        KERNEL8x4_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL8x4_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]

        KERNEL8x4_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL8x4_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        KERNEL8x4_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL8x4_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L4_M8_22

.Ldgemm_kernel_L4_M8_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Ldgemm_kernel_L4_M8_100

.Ldgemm_kernel_L4_M8_42:

        KERNEL8x4_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L4_M8_42

.Ldgemm_kernel_L4_M8_100:

        SAVE8x4

.Ldgemm_kernel_L4_M8_END:

.Ldgemm_kernel_L4_M4_BEGIN:

        mov    counterI, origM
        tst    counterI , #7
        ble    .Ldgemm_kernel_L4_END

        tst    counterI, #4            // counterI = counterI / 2
        ble    .Ldgemm_kernel_L4_M2_BEGIN

.Ldgemm_kernel_L4_M4_20:

        INIT4x4

        mov    pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL , #0
        ble    .Ldgemm_kernel_L4_M4_40

        .align 5
.Ldgemm_kernel_L4_M4_22:

        KERNEL4x4_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL4x4_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        KERNEL4x4_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL4x4_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]

        KERNEL4x4_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL4x4_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        KERNEL4x4_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL4x4_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L4_M4_22

.Ldgemm_kernel_L4_M4_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Ldgemm_kernel_L4_M4_100

.Ldgemm_kernel_L4_M4_42:

        KERNEL4x4_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L4_M4_42

.Ldgemm_kernel_L4_M4_100:

        SAVE4x4

.Ldgemm_kernel_L4_M4_END:

.Ldgemm_kernel_L4_M2_BEGIN:

        mov    counterI, origM
        tst    counterI , #3
        ble    .Ldgemm_kernel_L4_END

        tst    counterI, #2            // counterI = counterI / 2
        ble    .Ldgemm_kernel_L4_M1_BEGIN

.Ldgemm_kernel_L4_M2_20:

        INIT2x4

        mov    pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL , #0
        ble    .Ldgemm_kernel_L4_M2_40

        .align 5
.Ldgemm_kernel_L4_M2_22:

        KERNEL2x4_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL2x4_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        KERNEL2x4_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL2x4_SUB

        KERNEL2x4_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL2x4_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        KERNEL2x4_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL2x4_SUB

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L4_M2_22

.Ldgemm_kernel_L4_M2_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Ldgemm_kernel_L4_M2_100

        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE+64]

.Ldgemm_kernel_L4_M2_42:

        KERNEL2x4_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L4_M2_42

.Ldgemm_kernel_L4_M2_100:

        SAVE2x4

.Ldgemm_kernel_L4_M2_END:


.Ldgemm_kernel_L4_M1_BEGIN:

        tst    counterI, #1            // counterI = counterI % 2
        ble    .Ldgemm_kernel_L4_END

.Ldgemm_kernel_L4_M1_20:

        INIT1x4

        mov    pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL , #0
        ble    .Ldgemm_kernel_L4_M1_40

        .align 5
.Ldgemm_kernel_L4_M1_22:

        KERNEL1x4_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL1x4_SUB
        KERNEL1x4_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL1x4_SUB

        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]

        KERNEL1x4_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL1x4_SUB
        KERNEL1x4_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL1x4_SUB

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L4_M1_22

.Ldgemm_kernel_L4_M1_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Ldgemm_kernel_L4_M1_100

        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]

.Ldgemm_kernel_L4_M1_42:

        KERNEL1x4_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L4_M1_42

.Ldgemm_kernel_L4_M1_100:

        SAVE1x4

.Ldgemm_kernel_L4_END:

        lsl    temp, origK, #5
        add    origPB, origPB, temp        // B = B + K * 4 * 8

        subs    counterJ, counterJ , #1        // j--
        bgt    .Ldgemm_kernel_L4_BEGIN


/******************************************************************************/

.Ldgemm_kernel_L2_BEGIN:   // less than 2 left in N direction

        mov    counterJ , origN
        tst    counterJ , #3
        ble    .Ldgemm_kernel_L999   // error, N was less than 4?

        tst    counterJ , #2
        ble    .Ldgemm_kernel_L1_BEGIN

        mov    pCRow0, pC
        add    pCRow1, pCRow0, LDC

        add    pC, pCRow1, LDC

        mov    pA, origPA            // pA = A

.Ldgemm_kernel_L2_M16_BEGIN: //todo

        mov    counterI, origM
        asr     counterI, counterI, #4        // counterI = counterI / 16
        cmp    counterI, #0
        ble    .Ldgemm_kernel_L2_M8_BEGIN

        .align 5
.Ldgemm_kernel_L2_M16_20:

        INIT16x2

        mov    pB, origPB

        asr    counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL,#0
        ble    .Ldgemm_kernel_L2_M16_40

        .align 5
.Ldgemm_kernel_L2_M16_22:

        KERNEL16x2_SUB
        KERNEL16x2_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL16x2_SUB
        KERNEL16x2_SUB

        KERNEL16x2_SUB
        KERNEL16x2_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL16x2_SUB
        KERNEL16x2_SUB

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L2_M16_22

.Ldgemm_kernel_L2_M16_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Ldgemm_kernel_L2_M16_100

        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE+64]

.Ldgemm_kernel_L2_M16_42:

        KERNEL16x2_SUB

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L2_M16_42

.Ldgemm_kernel_L2_M16_100:

        SAVE16x2

.Ldgemm_kernel_L2_M16_END:

        subs    counterI, counterI, #1
        bgt    .Ldgemm_kernel_L2_M16_20

.Ldgemm_kernel_L2_M8_BEGIN:  //todo

        mov    counterI, origM
        tst    counterI , #15
        ble    .Ldgemm_kernel_L2_END

        tst    counterI, #8            // counterI = counterI / 2
        ble    .Ldgemm_kernel_L2_M4_BEGIN

        .align 5
.Ldgemm_kernel_L2_M8_20:

        INIT8x2

        mov    pB, origPB

        asr    counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL,#0
        ble    .Ldgemm_kernel_L2_M8_40

        .align 5
.Ldgemm_kernel_L2_M8_22:

        KERNEL8x2_SUB
        KERNEL8x2_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL8x2_SUB
        KERNEL8x2_SUB

        KERNEL8x2_SUB
        KERNEL8x2_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL8x2_SUB
        KERNEL8x2_SUB

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L2_M8_22

.Ldgemm_kernel_L2_M8_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Ldgemm_kernel_L2_M8_100

        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE+64]

.Ldgemm_kernel_L2_M8_42:

        KERNEL8x2_SUB

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L2_M8_42

.Ldgemm_kernel_L2_M8_100:

        SAVE8x2

.Ldgemm_kernel_L2_M8_END:

.Ldgemm_kernel_L2_M4_BEGIN:

        mov    counterI, origM
        tst    counterI , #7
        ble    .Ldgemm_kernel_L2_END

        tst    counterI, #4            // counterI = counterI / 2
        ble    .Ldgemm_kernel_L2_M2_BEGIN

.Ldgemm_kernel_L2_M4_20:

        INIT4x2

        mov    pB, origPB

        asr    counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL,#0
        ble    .Ldgemm_kernel_L2_M4_40

        .align 5
.Ldgemm_kernel_L2_M4_22:

        KERNEL4x2_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        KERNEL4x2_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL4x2_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        KERNEL4x2_SUB

        KERNEL4x2_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        KERNEL4x2_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL4x2_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        KERNEL4x2_SUB

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L2_M4_22

.Ldgemm_kernel_L2_M4_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Ldgemm_kernel_L2_M4_100

        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE+64]

.Ldgemm_kernel_L2_M4_42:

        KERNEL4x2_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L2_M4_42

.Ldgemm_kernel_L2_M4_100:

        SAVE4x2

.Ldgemm_kernel_L2_M4_END:

.Ldgemm_kernel_L2_M2_BEGIN:

        mov    counterI, origM
        tst    counterI , #3
        ble    .Ldgemm_kernel_L2_END

        tst    counterI, #2            // counterI = counterI / 2
        ble    .Ldgemm_kernel_L2_M1_BEGIN

.Ldgemm_kernel_L2_M2_20:

        INIT2x2

        mov    pB, origPB

        asr    counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL,#0
        ble    .Ldgemm_kernel_L2_M2_40

.Ldgemm_kernel_L2_M2_22:

        KERNEL2x2_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL2x2_SUB
        KERNEL2x2_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        KERNEL2x2_SUB

        KERNEL2x2_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL2x2_SUB
        KERNEL2x2_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        KERNEL2x2_SUB

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L2_M2_22

        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE+64]
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE+64]

.Ldgemm_kernel_L2_M2_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Ldgemm_kernel_L2_M2_100

.Ldgemm_kernel_L2_M2_42:

        KERNEL2x2_SUB

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L2_M2_42

.Ldgemm_kernel_L2_M2_100:

        SAVE2x2

.Ldgemm_kernel_L2_M2_END:

.Ldgemm_kernel_L2_M1_BEGIN:

        tst    counterI, #1            // counterI = counterI % 2
        ble    .Ldgemm_kernel_L2_END

.Ldgemm_kernel_L2_M1_20:

        INIT1x2

        mov    pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp     counterL, #0
        ble    .Ldgemm_kernel_L2_M1_40

.Ldgemm_kernel_L2_M1_22:

        KERNEL1x2_SUB
        KERNEL1x2_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL1x2_SUB
        KERNEL1x2_SUB

        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]

        KERNEL1x2_SUB
        KERNEL1x2_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL1x2_SUB
        KERNEL1x2_SUB

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L2_M1_22

        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE+64]

.Ldgemm_kernel_L2_M1_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Ldgemm_kernel_L2_M1_100

.Ldgemm_kernel_L2_M1_42:

        KERNEL1x2_SUB

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L2_M1_42

.Ldgemm_kernel_L2_M1_100:

        SAVE1x2

.Ldgemm_kernel_L2_END:
        add    origPB, origPB, origK, lsl #4    // B = B + K * 2 * 8

/******************************************************************************/

.Ldgemm_kernel_L1_BEGIN:

        mov    counterJ , origN
        tst    counterJ , #1
        ble    .Ldgemm_kernel_L999 // done

        mov    pCRow0, pC            // pCRow0 = C
        add    pC , pC , LDC            // Update pC to point to next

        mov    pA, origPA            // pA = A

.Ldgemm_kernel_L1_M16_BEGIN:

        mov    counterI, origM
        asr     counterI, counterI, #4        // counterI = counterI /16
        cmp    counterI, #0
        ble    .Ldgemm_kernel_L1_M8_BEGIN

        .align 5
.Ldgemm_kernel_L1_M16_20:

        INIT16x1

        mov    pB, origPB
        asr    counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL , #0
        ble    .Ldgemm_kernel_L1_M16_40

        .align 5
.Ldgemm_kernel_L1_M16_22:
        KERNEL16x1_SUB
        KERNEL16x1_SUB
        KERNEL16x1_SUB
        KERNEL16x1_SUB

        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]

        KERNEL16x1_SUB
        KERNEL16x1_SUB
        KERNEL16x1_SUB
        KERNEL16x1_SUB

        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L1_M16_22

.Ldgemm_kernel_L1_M16_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Ldgemm_kernel_L1_M16_100

        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
.Ldgemm_kernel_L1_M16_42:

        KERNEL16x1_SUB

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L1_M16_42

.Ldgemm_kernel_L1_M16_100:

        SAVE16x1

.Ldgemm_kernel_L1_M16_END:

        subs    counterI, counterI, #1
        bgt    .Ldgemm_kernel_L1_M16_20


.Ldgemm_kernel_L1_M8_BEGIN:
        mov    counterI, origM
        tst    counterI , #15
        ble    .Ldgemm_kernel_L1_END

        tst    counterI, #8            // counterI = counterI / 2
        ble    .Ldgemm_kernel_L1_M4_BEGIN

        .align 5
.Ldgemm_kernel_L1_M8_20:

        INIT8x1

        mov    pB, origPB
        asr    counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL , #0
        ble    .Ldgemm_kernel_L1_M8_40

        .align 5
.Ldgemm_kernel_L1_M8_22:

        KERNEL8x1_SUB
        KERNEL8x1_SUB
        KERNEL8x1_SUB
        KERNEL8x1_SUB

        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]

        KERNEL8x1_SUB
        KERNEL8x1_SUB
        KERNEL8x1_SUB
        KERNEL8x1_SUB

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L1_M8_22

.Ldgemm_kernel_L1_M8_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Ldgemm_kernel_L1_M8_100

        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]

.Ldgemm_kernel_L1_M8_42:

        KERNEL8x1_SUB

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L1_M8_42

.Ldgemm_kernel_L1_M8_100:

        SAVE8x1

.Ldgemm_kernel_L1_M8_END:

.Ldgemm_kernel_L1_M4_BEGIN:

        mov    counterI, origM
        tst    counterI , #7
        ble    .Ldgemm_kernel_L1_END

        tst    counterI, #4            // counterI = counterI / 2
        ble    .Ldgemm_kernel_L1_M2_BEGIN

.Ldgemm_kernel_L1_M4_20:

        INIT4x1

        mov    pB, origPB
        asr    counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL , #0
        ble    .Ldgemm_kernel_L1_M4_40

        .align 5
.Ldgemm_kernel_L1_M4_22:

        KERNEL4x1_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        KERNEL4x1_SUB
        KERNEL4x1_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        KERNEL4x1_SUB

        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]

        KERNEL4x1_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        KERNEL4x1_SUB
        KERNEL4x1_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        KERNEL4x1_SUB

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L1_M4_22

.Ldgemm_kernel_L1_M4_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Ldgemm_kernel_L1_M4_100

        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
.Ldgemm_kernel_L1_M4_42:

        KERNEL4x1_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L1_M4_42

.Ldgemm_kernel_L1_M4_100:

        SAVE4x1

.Ldgemm_kernel_L1_M4_END:

.Ldgemm_kernel_L1_M2_BEGIN:

        mov    counterI, origM
        tst    counterI , #3
        ble    .Ldgemm_kernel_L1_END

        tst    counterI, #2            // counterI = counterI / 2
        ble    .Ldgemm_kernel_L1_M1_BEGIN

.Ldgemm_kernel_L1_M2_20:

        INIT2x1

        mov    pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp     counterL , #0
        ble     .Ldgemm_kernel_L1_M2_40

.Ldgemm_kernel_L1_M2_22:

        KERNEL2x1_SUB
        KERNEL2x1_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        KERNEL2x1_SUB
        KERNEL2x1_SUB

        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]

        KERNEL2x1_SUB
        KERNEL2x1_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        KERNEL2x1_SUB
        KERNEL2x1_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L1_M2_22

        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE+64]
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]

.Ldgemm_kernel_L1_M2_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble     .Ldgemm_kernel_L1_M2_100

.Ldgemm_kernel_L1_M2_42:

        KERNEL2x1_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L1_M2_42

.Ldgemm_kernel_L1_M2_100:

        SAVE2x1

.Ldgemm_kernel_L1_M2_END:


.Ldgemm_kernel_L1_M1_BEGIN:

        tst    counterI, #1            // counterI = counterI % 2
        ble    .Ldgemm_kernel_L1_END

.Ldgemm_kernel_L1_M1_20:

        INIT1x1

        mov    pB, origPB

        asr    counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL , #0
        ble    .Ldgemm_kernel_L1_M1_40

.Ldgemm_kernel_L1_M1_22:

        KERNEL1x1_SUB
        KERNEL1x1_SUB
        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        KERNEL1x1_SUB
        KERNEL1x1_SUB

        KERNEL1x1_SUB
        KERNEL1x1_SUB
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]
        KERNEL1x1_SUB
        KERNEL1x1_SUB

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L1_M1_22


.Ldgemm_kernel_L1_M1_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Ldgemm_kernel_L1_M1_100

        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        prfm    PLDL1KEEP, [pB, #B_PRE_SIZE]

.Ldgemm_kernel_L1_M1_42:

        KERNEL1x1_SUB

        subs    counterL, counterL, #1
        bgt    .Ldgemm_kernel_L1_M1_42

.Ldgemm_kernel_L1_M1_100:

        SAVE1x1


.Ldgemm_kernel_L1_END:


.Ldgemm_kernel_L999:

        mov    x0, #0                // set return value
        ldp    d8, d9, [sp, #(0 * 16)]
        ldp    d10, d11, [sp, #(1 * 16)]
        ldp    d12, d13, [sp, #(2 * 16)]
        ldp    d14, d15, [sp, #(3 * 16)]
        ldp    d16, d17, [sp, #(4 * 16)]
        ldp    x18, x19, [sp, #(5 * 16)]
        ldp    x20, x21, [sp, #(6 * 16)]
        ldp    x22, x23, [sp, #(7 * 16)]
        ldp    x24, x25, [sp, #(8 * 16)]
        ldp    x26, x27, [sp, #(9 * 16)]
        ldr    x28, [sp, #(10 * 16)]
        add    sp, sp, #(11*16)
        ret

    EPILOGUE



