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
#define pCRow4          x18
#define pCRow5          x19
#define pCRow6          x20
#define pCRow7          x21

#define alpha0          z3.d
#define alphaV0         z3.d[0]

#define A_PRE_SIZE      256
#define B_PRE_SIZE      256
#define C_PRE_SIZE      256


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
.macro INIT12x8
        fmov    d8,  xzr
        fmov    d9,  xzr
        fmov    d10, d8
        fmov    d11, d9
        fmov    d12, d8
        fmov    d13, d9
        fmov    d14, d8
        fmov    d15, d9
        fmov    d16, xzr
        fmov    d17, xzr
        fmov    d18, d16
        fmov    d19, d17
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

.macro KERNEL12x8_I
        ld1d    {z0.d}, p0/z, [pA]
        ld1d    {z1.d}, p0/z, [pA, #1, MUL VL]
        ld1d    {z2.d}, p0/z, [pA, #2, MUL VL]

        ld1rqd  {z3.d}, p2/z, [pB]
        ld1rqd  {z4.d}, p2/z, [pB, #16]
        ld1rqd  {z5.d}, p2/z, [pB, #32]
        ld1rqd  {z6.d}, p2/z, [pB, #48]

        fmul    z8.d,  z0.d, z3.d[0]
        fmul    z11.d, z0.d, z3.d[1]
        fmul    z14.d, z0.d, z4.d[0]
        fmul    z17.d, z0.d, z4.d[1]
        fmul    z20.d, z0.d, z5.d[0]
        fmul    z23.d, z0.d, z5.d[1]
        fmul    z26.d, z0.d, z6.d[0]
        fmul    z29.d, z0.d, z6.d[1]

        ld1d    {z0.d}, p0/z, [pA, #3, MUL VL]

        fmul    z9.d,  z1.d, z3.d[0]
        fmul    z12.d, z1.d, z3.d[1]
        fmul    z10.d, z2.d, z3.d[0]
        fmul    z13.d, z2.d, z3.d[1]

        ld1rqd  {z3.d}, p2/z, [pB, #64]

        fmul    z15.d, z1.d, z4.d[0]
        fmul    z18.d, z1.d, z4.d[1]
        fmul    z16.d, z2.d, z4.d[0]
        fmul    z19.d, z2.d, z4.d[1]

        ld1rqd  {z4.d}, p2/z, [pB, #80]

        fmul    z21.d, z1.d, z5.d[0]
        fmul    z24.d, z1.d, z5.d[1]
        fmul    z27.d, z1.d, z6.d[0]
        fmul    z30.d, z1.d, z6.d[1]

        ld1d    {z1.d}, p0/z, [pA, #4, MUL VL]

        fmul    z22.d, z2.d, z5.d[0]
        fmul    z25.d, z2.d, z5.d[1]
        fmul    z28.d, z2.d, z6.d[0]
        fmul    z31.d, z2.d, z6.d[1]

        ld1d    {z2.d}, p0/z, [pA, #5, MUL VL]
        incb    pA, all, MUL #6
        ld1rqd  {z5.d}, p2/z, [pB, #96]
        ld1rqd  {z6.d}, p2/z, [pB, #112]
        incb    pB, all, MUL #4
.endm

.macro KERNEL12x8_M1
        fmla    z8.d,  z0.d, z3.d[0]
        fmla    z11.d, z0.d, z3.d[1]
        fmla    z14.d, z0.d, z4.d[0]
        fmla    z17.d, z0.d, z4.d[1]
        fmla    z20.d, z0.d, z5.d[0]
        fmla    z23.d, z0.d, z5.d[1]
        fmla    z26.d, z0.d, z6.d[0]
        fmla    z29.d, z0.d, z6.d[1]

        ld1d    {z0.d}, p0/z, [pA]

        fmla    z9.d,  z1.d, z3.d[0]
        fmla    z12.d, z1.d, z3.d[1]
        fmla    z10.d, z2.d, z3.d[0]
        fmla    z13.d, z2.d, z3.d[1]

        ld1rqd  {z3.d}, p2/z, [pB]

        fmla    z15.d, z1.d, z4.d[0]
        fmla    z18.d, z1.d, z4.d[1]
        fmla    z16.d, z2.d, z4.d[0]
        fmla    z19.d, z2.d, z4.d[1]

        ld1rqd  {z4.d}, p2/z, [pB, #16]

        fmla    z21.d, z1.d, z5.d[0]
        fmla    z24.d, z1.d, z5.d[1]
        fmla    z22.d, z2.d, z5.d[0]
        fmla    z25.d, z2.d, z5.d[1]

        ld1rqd  {z5.d}, p2/z, [pB, #32]

        fmla    z27.d, z1.d, z6.d[0]
        fmla    z30.d, z1.d, z6.d[1]
        fmla    z28.d, z2.d, z6.d[0]
        fmla    z31.d, z2.d, z6.d[1]

        ld1rqd  {z6.d}, p2/z, [pB, #48]
        incb    pB, all, MUL #2

        prfm    PLDL1KEEP, [pA, #A_PRE_SIZE]
        ld1d    {z1.d}, p0/z, [pA, #1, MUL VL]
        ld1d    {z2.d}, p0/z, [pA, #2, MUL VL]
        incb    pA, all, MUL #3
.endm

.macro KERNEL12x8_E
        fmla    z8.d,  z0.d, z3.d[0]
        fmla    z11.d, z0.d, z3.d[1]
        fmla    z14.d, z0.d, z4.d[0]
        fmla    z17.d, z0.d, z4.d[1]
        fmla    z20.d, z0.d, z5.d[0]
        fmla    z23.d, z0.d, z5.d[1]
        fmla    z26.d, z0.d, z6.d[0]
        fmla    z29.d, z0.d, z6.d[1]

        fmla    z9.d,  z1.d, z3.d[0]
        fmla    z12.d, z1.d, z3.d[1]
        fmla    z10.d, z2.d, z3.d[0]
        fmla    z13.d, z2.d, z3.d[1]
        fmla    z15.d, z1.d, z4.d[0]
        fmla    z18.d, z1.d, z4.d[1]
        fmla    z16.d, z2.d, z4.d[0]
        fmla    z19.d, z2.d, z4.d[1]

        fmla    z21.d, z1.d, z5.d[0]
        fmla    z24.d, z1.d, z5.d[1]
        fmla    z22.d, z2.d, z5.d[0]
        fmla    z25.d, z2.d, z5.d[1]
        fmla    z27.d, z1.d, z6.d[0]
        fmla    z30.d, z1.d, z6.d[1]
        fmla    z28.d, z2.d, z6.d[0]
        fmla    z31.d, z2.d, z6.d[1]
.endm

.macro KERNEL12x8_SUB
        ld1d    {z0.d}, p0/z, [pA]
        ld1d    {z1.d}, p0/z, [pA, #1, MUL VL]
        ld1d    {z2.d}, p0/z, [pA, #2, MUL VL]
        incb    pA, all, MUL #3

        ld1rqd  {z3.d}, p2/z, [pB]
        ld1rqd  {z4.d}, p2/z, [pB, #16]

        fmla    z8.d,  z0.d, z3.d[0]
        fmla    z11.d, z0.d, z3.d[1]
        fmla    z9.d,  z1.d, z3.d[0]
        fmla    z12.d, z1.d, z3.d[1]
        fmla    z10.d, z2.d, z3.d[0]
        fmla    z13.d, z2.d, z3.d[1]

        ld1rqd  {z5.d}, p2/z, [pB, #32]

        fmla    z14.d, z0.d, z4.d[0]
        fmla    z17.d, z0.d, z4.d[1]
        fmla    z15.d, z1.d, z4.d[0]
        fmla    z18.d, z1.d, z4.d[1]
        fmla    z16.d, z2.d, z4.d[0]
        fmla    z19.d, z2.d, z4.d[1]

        ld1rqd  {z6.d}, p2/z, [pB, #48]
        incb    pB, all, MUL #2

        fmla    z22.d, z2.d, z5.d[0]
        fmla    z25.d, z2.d, z5.d[1]
        fmla    z20.d, z0.d, z5.d[0]
        fmla    z23.d, z0.d, z5.d[1]
        fmla    z21.d, z1.d, z5.d[0]
        fmla    z24.d, z1.d, z5.d[1]

        fmla    z26.d, z0.d, z6.d[0]
        fmla    z29.d, z0.d, z6.d[1]
        fmla    z27.d, z1.d, z6.d[0]
        fmla    z30.d, z1.d, z6.d[1]
        fmla    z28.d, z2.d, z6.d[0]
        fmla    z31.d, z2.d, z6.d[1]
.endm

.macro SAVE12x8
        dup     alpha0, alpha

        prfm    PLDL1KEEP, [pCRow0, #C_PRE_SIZE]
        ld1d    {z0.d}, p0/z, [pCRow0]
        fmla    z0.d, z8.d, alphaV0
        st1d    {z0.d}, p0, [pCRow0]

        ld1d    {z1.d}, p0/z, [pCRow0, #1, MUL VL]
        fmla    z1.d, z9.d, alphaV0
        st1d    {z1.d}, p0, [pCRow0, #1, MUL VL]

        ld1d    {z2.d}, p0/z, [pCRow0, #2, MUL VL]
        fmla    z2.d, z10.d, alphaV0
        st1d    {z2.d}, p0, [pCRow0, #2, MUL VL]
        incb    pCRow0, all, MUL #3

        prfm    PLDL1KEEP, [pCRow1, #C_PRE_SIZE]
        ld1d    {z4.d}, p0/z, [pCRow1]
        fmla    z4.d, z11.d, alphaV0
        st1d    {z4.d}, p0, [pCRow1]

        ld1d    {z5.d}, p0/z, [pCRow1, #1, MUL VL]
        fmla    z5.d, z12.d, alphaV0
        st1d    {z5.d}, p0, [pCRow1, #1, MUL VL]

        ld1d    {z6.d}, p0/z, [pCRow1, #2, MUL VL]
        fmla    z6.d, z13.d, alphaV0
        st1d    {z6.d}, p0, [pCRow1, #2, MUL VL]
        incb    pCRow1, all, MUL #3

        prfm    PLDL1KEEP, [pCRow2, #C_PRE_SIZE]
        ld1d    {z0.d}, p0/z, [pCRow2]
        fmla    z0.d, z14.d, alphaV0
        st1d    {z0.d}, p0, [pCRow2]

        ld1d    {z1.d}, p0/z, [pCRow2, #1, MUL VL]
        fmla    z1.d, z15.d, alphaV0
        st1d    {z1.d}, p0, [pCRow2, #1, MUL VL]

        ld1d    {z2.d}, p0/z, [pCRow2, #2, MUL VL]
        fmla    z2.d, z16.d, alphaV0
        st1d    {z2.d}, p0, [pCRow2, #2, MUL VL]
        incb    pCRow2, all, MUL #3

        prfm    PLDL1KEEP, [pCRow3, #C_PRE_SIZE]
        ld1d    {z4.d}, p0/z, [pCRow3]
        fmla    z4.d, z17.d, alphaV0
        st1d    {z4.d}, p0, [pCRow3]

        ld1d    {z5.d}, p0/z, [pCRow3, #1, MUL VL]
        fmla    z5.d, z18.d, alphaV0
        st1d    {z5.d}, p0, [pCRow3, #1, MUL VL]

        ld1d    {z6.d}, p0/z, [pCRow3, #2, MUL VL]
        fmla    z6.d, z19.d, alphaV0
        st1d    {z6.d}, p0, [pCRow3, #2, MUL VL]
        incb    pCRow3, all, MUL #3

        prfm    PLDL1KEEP, [pCRow4, #C_PRE_SIZE]
        ld1d    {z0.d}, p0/z, [pCRow4]
        fmla    z0.d, z20.d, alphaV0
        st1d    {z0.d}, p0, [pCRow4]

        ld1d    {z1.d}, p0/z, [pCRow4, #1, MUL VL]
        fmla    z1.d, z21.d, alphaV0
        st1d    {z1.d}, p0, [pCRow4, #1, MUL VL]

        ld1d    {z2.d}, p0/z, [pCRow4, #2, MUL VL]
        fmla    z2.d, z22.d, alphaV0
        st1d    {z2.d}, p0, [pCRow4, #2, MUL VL]
        incb    pCRow4, all, MUL #3

        prfm    PLDL1KEEP, [pCRow5, #C_PRE_SIZE]
        ld1d    {z4.d}, p0/z, [pCRow5]
        fmla    z4.d, z23.d, alphaV0
        st1d    {z4.d}, p0, [pCRow5]

        ld1d    {z5.d}, p0/z, [pCRow5, #1, MUL VL]
        fmla    z5.d, z24.d, alphaV0
        st1d    {z5.d}, p0, [pCRow5, #1, MUL VL]

        ld1d    {z6.d}, p0/z, [pCRow5, #2, MUL VL]
        fmla    z6.d, z25.d, alphaV0
        st1d    {z6.d}, p0, [pCRow5, #2, MUL VL]
        incb    pCRow5, all, MUL #3

        prfm    PLDL1KEEP, [pCRow6, #C_PRE_SIZE]
        ld1d    {z0.d}, p0/z, [pCRow6]
        fmla    z0.d, z26.d, alphaV0
        st1d    {z0.d}, p0, [pCRow6]

        ld1d    {z1.d}, p0/z, [pCRow6, #1, MUL VL]
        fmla    z1.d, z27.d, alphaV0
        st1d    {z1.d}, p0, [pCRow6, #1, MUL VL]

        ld1d    {z2.d}, p0/z, [pCRow6, #2, MUL VL]
        fmla    z2.d, z28.d, alphaV0
        st1d    {z2.d}, p0, [pCRow6, #2, MUL VL]
        incb    pCRow6, all, MUL #3

        // load miss > 100 (11)
        prfm    PLDL1KEEP, [pCRow7]
        prfm    PLDL1KEEP, [pCRow7, #C_PRE_SIZE+128]
        ld1d    {z4.d}, p0/z, [pCRow7]
        fmla    z4.d, z29.d, alphaV0
        st1d    {z4.d}, p0, [pCRow7]

        ld1d    {z5.d}, p0/z, [pCRow7, #1, MUL VL]
        fmla    z5.d, z30.d, alphaV0
        st1d    {z5.d}, p0, [pCRow7, #1, MUL VL]

        ld1d    {z6.d}, p0/z, [pCRow7, #2, MUL VL]
        fmla    z6.d, z31.d, alphaV0
        st1d    {z6.d}, p0, [pCRow7, #2, MUL VL]
        incb    pCRow7, all, MUL #3
.endm

/******************************************************************************/
.macro INIT8x8
        fmov    d16,  xzr
        fmov    d17,  xzr
        fmov    d18, d16
        fmov    d19, d16
        fmov    d20, d16
        fmov    d21, d16
        fmov    d22, xzr
        fmov    d23, d18
        fmov    d24, xzr
        fmov    d25, d19
        fmov    d26, d20
        fmov    d27, xzr
        fmov    d28, d21
        fmov    d29, d23
        fmov    d30, d25
        fmov    d31, d26
.endm

.macro KERNEL8x8_I
        ld1d    {z0.d}, p0/z, [pA]
        ld1d    {z1.d}, p0/z, [pA, #1, MUL VL]

        ld1rqd  {z2.d}, p2/z, [pB]
        ld1rqd  {z3.d}, p2/z, [pB, #16]
        ld1rqd  {z4.d}, p2/z, [pB, #32]
        ld1rqd  {z5.d}, p2/z, [pB, #48]

        fmul    z16.d, z0.d, z2.d[0]
        fmul    z18.d, z0.d, z2.d[1]
        fmul    z20.d, z0.d, z3.d[0]
        fmul    z22.d, z0.d, z3.d[1]

        ld1d    {z6.d}, p0/z, [pA, #2, MUL VL]

        fmul    z24.d, z0.d, z4.d[0]
        fmul    z26.d, z0.d, z4.d[1]

        ld1d    {z7.d}, p0/z, [pA, #3, MUL VL]
        incb    pA, all, MUL #4

        fmul    z28.d, z0.d, z5.d[0]
        fmul    z30.d, z0.d, z5.d[1]

        ld1rqd  {z8.d}, p2/z, [pB, #64]

        fmul    z17.d, z1.d, z2.d[0]
        fmul    z19.d, z1.d, z2.d[1]

        ld1rqd  {z9.d}, p2/z, [pB, #80]

        fmul    z21.d, z1.d, z3.d[0]
        fmul    z23.d, z1.d, z3.d[1]

        ld1rqd  {z10.d}, p2/z, [pB, #96]

        fmul    z25.d, z1.d, z4.d[0]
        fmul    z27.d, z1.d, z4.d[1]

        ld1rqd  {z11.d}, p2/z, [pB, #112]
        incb    pB, all, MUL #4

        fmul    z29.d, z1.d, z5.d[0]
        fmul    z31.d, z1.d, z5.d[1]
.endm

.macro KERNEL8x8_M1
        fmla    z16.d, z0.d, z2.d[0]
        fmla    z18.d, z0.d, z2.d[1]
        fmla    z20.d, z0.d, z3.d[0]
        fmla    z22.d, z0.d, z3.d[1]

        ld1d    {z6.d}, p0/z, [pA]

        fmla    z24.d, z0.d, z4.d[0]
        fmla    z26.d, z0.d, z4.d[1]

        ld1d    {z7.d}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        fmla    z28.d, z0.d, z5.d[0]
        fmla    z30.d, z0.d, z5.d[1]

        ld1rqd  {z8.d}, p2/z, [pB]

        fmla    z17.d, z1.d, z2.d[0]
        fmla    z19.d, z1.d, z2.d[1]

        ld1rqd  {z9.d}, p2/z, [pB, #16]

        fmla    z21.d, z1.d, z3.d[0]
        fmla    z23.d, z1.d, z3.d[1]

        ld1rqd  {z10.d}, p2/z, [pB, #32]

        fmla    z25.d, z1.d, z4.d[0]
        fmla    z27.d, z1.d, z4.d[1]

        ld1rqd  {z11.d}, p2/z, [pB, #48]
        incb    pB, all, MUL #2

        fmla    z29.d, z1.d, z5.d[0]
        fmla    z31.d, z1.d, z5.d[1]
.endm

.macro KERNEL8x8_M2
        fmla    z16.d, z6.d, z8.d[0]
        fmla    z18.d, z6.d, z8.d[1]
        fmla    z20.d, z6.d, z9.d[0]
        fmla    z22.d, z6.d, z9.d[1]

        ld1d    {z0.d}, p0/z, [pA]

        fmla    z24.d, z6.d, z10.d[0]
        fmla    z26.d, z6.d, z10.d[1]

        ld1d    {z1.d}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        fmla    z28.d, z6.d, z11.d[0]
        fmla    z30.d, z6.d, z11.d[1]

        ld1rqd  {z2.d}, p2/z, [pB]

        fmla    z17.d, z7.d, z8.d[0]
        fmla    z19.d, z7.d, z8.d[1]

        ld1rqd  {z3.d}, p2/z, [pB, #16]

        fmla    z21.d, z7.d, z9.d[0]
        fmla    z23.d, z7.d, z9.d[1]

        ld1rqd  {z4.d}, p2/z, [pB, #32]

        fmla    z25.d, z7.d, z10.d[0]
        fmla    z27.d, z7.d, z10.d[1]

        ld1rqd  {z5.d}, p2/z, [pB, #48]
        incb    pB, all, MUL #2

        fmla    z29.d, z7.d, z11.d[0]
        fmla    z31.d, z7.d, z11.d[1]
.endm

.macro KERNEL8x8_E
        fmla    z16.d, z6.d, z8.d[0]
        fmla    z18.d, z6.d, z8.d[1]
        fmla    z20.d, z6.d, z9.d[0]
        fmla    z22.d, z6.d, z9.d[1]
        fmla    z24.d, z6.d, z10.d[0]
        fmla    z26.d, z6.d, z10.d[1]
        fmla    z28.d, z6.d, z11.d[0]
        fmla    z30.d, z6.d, z11.d[1]

        fmla    z17.d, z7.d, z8.d[0]
        fmla    z19.d, z7.d, z8.d[1]
        fmla    z21.d, z7.d, z9.d[0]
        fmla    z23.d, z7.d, z9.d[1]
        fmla    z25.d, z7.d, z10.d[0]
        fmla    z27.d, z7.d, z10.d[1]
        fmla    z29.d, z7.d, z11.d[0]
        fmla    z31.d, z7.d, z11.d[1]
.endm

.macro KERNEL8x8_SUB
        ld1d    {z6.d}, p0/z, [pA]
        ld1rqd  {z8.d}, p2/z, [pB]
        ld1rqd  {z9.d}, p2/z, [pB, #16]

        fmla    z16.d, z6.d, z8.d[0]
        fmla    z18.d, z6.d, z8.d[1]

        ld1d    {z7.d}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        fmla    z20.d, z6.d, z9.d[0]
        fmla    z22.d, z6.d, z9.d[1]

        ld1rqd  {z10.d}, p2/z, [pB, #32]

        fmla    z17.d, z7.d, z8.d[0]
        fmla    z19.d, z7.d, z8.d[1]
        fmla    z21.d, z7.d, z9.d[0]
        fmla    z23.d, z7.d, z9.d[1]

        ld1rqd  {z11.d}, p2/z, [pB, #48]
        incb    pB, all, MUL #2

        fmla    z24.d, z6.d, z10.d[0]
        fmla    z26.d, z6.d, z10.d[1]
        fmla    z25.d, z7.d, z10.d[0]
        fmla    z27.d, z7.d, z10.d[1]

        fmla    z28.d, z6.d, z11.d[0]
        fmla    z30.d, z6.d, z11.d[1]
        fmla    z29.d, z7.d, z11.d[0]
        fmla    z31.d, z7.d, z11.d[1]
.endm

.macro SAVE8x8
        dup     alpha0, alpha

        ld1d    {z0.d}, p0/z, [pCRow0]
        ld1d    {z1.d}, p0/z, [pCRow0, #1, MUL VL]
        fmla    z0.d, z16.d, alphaV0
        fmla    z1.d, z17.d, alphaV0
        st1d    {z0.d}, p0, [pCRow0]
        st1d    {z1.d}, p0, [pCRow0, #1, MUL VL]
        incb    pCRow0, all, MUL #2

        ld1d    {z4.d}, p0/z, [pCRow1]
        ld1d    {z5.d}, p0/z, [pCRow1, #1, MUL VL]
        fmla    z4.d, z18.d, alphaV0
        fmla    z5.d, z19.d, alphaV0
        st1d    {z4.d}, p0, [pCRow1]
        st1d    {z5.d}, p0, [pCRow1, #1, MUL VL]
        incb    pCRow1, all, MUL #2

        ld1d    {z6.d}, p0/z, [pCRow2]
        ld1d    {z7.d}, p0/z, [pCRow2, #1, MUL VL]
        fmla    z6.d, z20.d, alphaV0
        fmla    z7.d, z21.d, alphaV0
        st1d    {z6.d}, p0, [pCRow2]
        st1d    {z7.d}, p0, [pCRow2, #1, MUL VL]
        incb    pCRow2, all, MUL #2

        ld1d    {z8.d}, p0/z, [pCRow3]
        ld1d    {z9.d}, p0/z, [pCRow3, #1, MUL VL]
        fmla    z8.d, z22.d, alphaV0
        fmla    z9.d, z23.d, alphaV0
        st1d    {z8.d}, p0, [pCRow3]
        st1d    {z9.d}, p0, [pCRow3, #1, MUL VL]
        incb    pCRow3, all, MUL #2

        ld1d    {z0.d}, p0/z, [pCRow4]
        ld1d    {z1.d}, p0/z, [pCRow4, #1, MUL VL]
        fmla    z0.d, z24.d, alphaV0
        fmla    z1.d, z25.d, alphaV0
        st1d    {z0.d}, p0, [pCRow4]
        st1d    {z1.d}, p0, [pCRow4, #1, MUL VL]
        incb    pCRow4, all, MUL #2

        ld1d    {z4.d}, p0/z, [pCRow5]
        ld1d    {z5.d}, p0/z, [pCRow5, #1, MUL VL]
        fmla    z4.d, z26.d, alphaV0
        fmla    z5.d, z27.d, alphaV0
        st1d    {z4.d}, p0, [pCRow5]
        st1d    {z5.d}, p0, [pCRow5, #1, MUL VL]
        incb    pCRow5, all, MUL #2

        ld1d    {z6.d}, p0/z, [pCRow6]
        ld1d    {z7.d}, p0/z, [pCRow6, #1, MUL VL]
        fmla    z6.d, z28.d, alphaV0
        fmla    z7.d, z29.d, alphaV0
        st1d    {z6.d}, p0, [pCRow6]
        st1d    {z7.d}, p0, [pCRow6, #1, MUL VL]
        incb    pCRow6, all, MUL #2

        ld1d    {z8.d}, p0/z, [pCRow7]
        ld1d    {z9.d}, p0/z, [pCRow7, #1, MUL VL]
        fmla    z8.d, z30.d, alphaV0
        fmla    z9.d, z31.d, alphaV0
        st1d    {z8.d}, p0, [pCRow7]
        st1d    {z9.d}, p0, [pCRow7, #1, MUL VL]
        incb    pCRow7, all, MUL #2
.endm

/******************************************************************************/
.macro INIT4x8
        fmov    d10, xzr
        fmov    d11, d10
        fmov    d12, d10
        fmov    d13, xzr
        fmov    d14, xzr
        fmov    d15, d11
        fmov    d16, d12
        fmov    d17, d15
.endm

.macro KERNEL4x8_I
        ld1d    {z0.d}, p0/z, [pA]

        ld1rqd  {z1.d}, p2/z, [pB]
        ld1rqd  {z2.d}, p2/z, [pB, #16]
        ld1rqd  {z3.d}, p2/z, [pB, #32]
        ld1rqd  {z4.d}, p2/z, [pB, #48]

        fmul    z10.d, z0.d, z1.d[0]
        fmul    z11.d, z0.d, z1.d[1]

        ld1d    {z5.d}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        fmul    z12.d, z0.d, z2.d[0]
        fmul    z13.d, z0.d, z2.d[1]

        ld1rqd  {z6.d}, p2/z, [pB, #64]

        fmul    z14.d, z0.d, z3.d[0]
        fmul    z15.d, z0.d, z3.d[1]

        ld1rqd  {z7.d}, p2/z, [pB, #80]

        fmul    z16.d, z0.d, z4.d[0]
        fmul    z17.d, z0.d, z4.d[1]

        ld1rqd  {z8.d}, p2/z, [pB, #96]
        ld1rqd  {z9.d}, p2/z, [pB, #112]
        add     pB, pB, #128
.endm

.macro KERNEL4x8_M1
        ld1d    {z5.d}, p0/z, [pA]
        incb    pA, all, MUL #1

        fmla    z10.d, z0.d, z1.d[0]
        fmla    z11.d, z0.d, z1.d[1]

        ld1rqd  {z6.d}, p2/z, [pB]

        fmla    z12.d, z0.d, z2.d[0]
        fmla    z13.d, z0.d, z2.d[1]

        ld1rqd  {z7.d}, p2/z, [pB, #16]

        fmla    z14.d, z0.d, z3.d[0]
        fmla    z15.d, z0.d, z3.d[1]

        ld1rqd  {z8.d}, p2/z, [pB, #32]

        fmla    z16.d, z0.d, z4.d[0]
        fmla    z17.d, z0.d, z4.d[1]

        ld1rqd  {z9.d}, p2/z, [pB, #48]
        add     pB, pB, #64
.endm

.macro KERNEL4x8_M2
        ld1d    {z0.d}, p0/z, [pA]
        incb    pA, all, MUL #1

        fmla    z10.d, z5.d, z6.d[0]
        fmla    z11.d, z5.d, z6.d[1]

        ld1rqd  {z1.d}, p2/z, [pB]

        fmla    z12.d, z5.d, z7.d[0]
        fmla    z13.d, z5.d, z7.d[1]

        ld1rqd  {z2.d}, p2/z, [pB, #16]

        fmla    z14.d, z5.d, z8.d[0]
        fmla    z15.d, z5.d, z8.d[1]

        ld1rqd  {z3.d}, p2/z, [pB, #32]

        fmla    z16.d, z5.d, z9.d[0]
        fmla    z17.d, z5.d, z9.d[1]

        ld1rqd  {z4.d}, p2/z, [pB, #48]
        add     pB, pB, #64
.endm

.macro KERNEL4x8_E
        fmla    z10.d, z5.d, z6.d[0]
        fmla    z11.d, z5.d, z6.d[1]
        fmla    z12.d, z5.d, z7.d[0]
        fmla    z13.d, z5.d, z7.d[1]
        fmla    z14.d, z5.d, z8.d[0]
        fmla    z15.d, z5.d, z8.d[1]
        fmla    z16.d, z5.d, z9.d[0]
        fmla    z17.d, z5.d, z9.d[1]
.endm

.macro KERNEL4x8_SUB
        ld1d    {z5.d}, p0/z, [pA]
        incb    pA, all, MUL #1

        ld1rqd  {z6.d}, p2/z, [pB]
        ld1rqd  {z7.d}, p2/z, [pB, #16]
        ld1rqd  {z8.d}, p2/z, [pB, #32]
        ld1rqd  {z9.d}, p2/z, [pB, #48]
        add     pB, pB, #64

        fmla    z10.d, z5.d, z6.d[0]
        fmla    z11.d, z5.d, z6.d[1]
        fmla    z12.d, z5.d, z7.d[0]
        fmla    z13.d, z5.d, z7.d[1]
        fmla    z14.d, z5.d, z8.d[0]
        fmla    z15.d, z5.d, z8.d[1]
        fmla    z16.d, z5.d, z9.d[0]
        fmla    z17.d, z5.d, z9.d[1]
.endm

.macro SAVE4x8
        dup     alpha0, alpha

        ld1d    {z0.d}, p0/z, [pCRow0]
        fmla    z0.d, z10.d, alphaV0
        st1d    {z0.d}, p0, [pCRow0]
        incb    pCRow0, all, MUL #1

        ld1d    {z1.d}, p0/z, [pCRow1]
        fmla    z1.d, z11.d, alphaV0
        st1d    {z1.d}, p0, [pCRow1]
        incb    pCRow1, all, MUL #1

        ld1d    {z2.d}, p0/z, [pCRow2]
        fmla    z2.d, z12.d, alphaV0
        st1d    {z2.d}, p0, [pCRow2]
        incb    pCRow2, all, MUL #1

        ld1d    {z4.d}, p0/z, [pCRow3]
        fmla    z4.d, z13.d, alphaV0
        st1d    {z4.d}, p0, [pCRow3]
        incb    pCRow3, all, MUL #1

        ld1d    {z5.d}, p0/z, [pCRow4]
        fmla    z5.d, z14.d, alphaV0
        st1d    {z5.d}, p0, [pCRow4]
        incb    pCRow4, all, MUL #1

        ld1d    {z6.d}, p0/z, [pCRow5]
        fmla    z6.d, z15.d, alphaV0
        st1d    {z6.d}, p0, [pCRow5]
        incb    pCRow5, all, MUL #1

        ld1d    {z7.d}, p0/z, [pCRow6]
        fmla    z7.d, z16.d, alphaV0
        st1d    {z7.d}, p0, [pCRow6]
        incb    pCRow6, all, MUL #1

        ld1d    {z8.d}, p0/z, [pCRow7]
        fmla    z8.d, z17.d, alphaV0
        st1d    {z8.d}, p0, [pCRow7]
        incb    pCRow7, all, MUL #1
.endm

/******************************************************************************/
.macro INIT2x8
        fmov    d10, xzr
        fmov    d11, d10
        fmov    d12, d10
        fmov    d13, xzr
        fmov    d14, xzr
        fmov    d15, d11
        fmov    d16, d12
        fmov    d17, d15
.endm

.macro KERNEL2x8_SUB
        ld1rqd  {z5.d}, p2/z, [pA]
        add     pA, pA, #16

        ld1rqd  {z6.d}, p2/z, [pB]
        ld1rqd  {z7.d}, p2/z, [pB, #16]
        ld1rqd  {z8.d}, p2/z, [pB, #32]
        ld1rqd  {z9.d}, p2/z, [pB, #48]
        add     pB, pB, #64

        fmla    z10.d, z5.d, z6.d[0]
        fmla    z11.d, z5.d, z6.d[1]
        fmla    z12.d, z5.d, z7.d[0]
        fmla    z13.d, z5.d, z7.d[1]
        fmla    z14.d, z5.d, z8.d[0]
        fmla    z15.d, z5.d, z8.d[1]
        fmla    z16.d, z5.d, z9.d[0]
        fmla    z17.d, z5.d, z9.d[1]
.endm

.macro SAVE2x8
        dup     alpha0, alpha

        ld1rqd  {z0.d}, p2/z, [pCRow0]
        fmla    z0.d, z10.d, alphaV0
        st1d    {z0.d}, p2, [pCRow0]
        add     pCRow0, pCRow0, #16

        ld1rqd  {z1.d}, p2/z, [pCRow1]
        fmla    z1.d, z11.d, alphaV0
        st1d    {z1.d}, p2, [pCRow1]
        add     pCRow1, pCRow1, #16

        ld1rqd  {z2.d}, p2/z, [pCRow2]
        fmla    z2.d, z12.d, alphaV0
        st1d    {z2.d}, p2, [pCRow2]
        add     pCRow2, pCRow2, #16

        ld1rqd  {z4.d}, p2/z, [pCRow3]
        fmla    z4.d, z13.d, alphaV0
        st1d    {z4.d}, p2, [pCRow3]
        add     pCRow3, pCRow3, #16

        ld1rqd  {z5.d}, p2/z, [pCRow4]
        fmla    z5.d, z14.d, alphaV0
        st1d    {z5.d}, p2, [pCRow4]
        add     pCRow4, pCRow4, #16

        ld1rqd  {z6.d}, p2/z, [pCRow5]
        fmla    z6.d, z15.d, alphaV0
        st1d    {z6.d}, p2, [pCRow5]
        add     pCRow5, pCRow5, #16

        ld1rqd  {z7.d}, p2/z, [pCRow6]
        fmla    z7.d, z16.d, alphaV0
        st1d    {z7.d}, p2, [pCRow6]
        add     pCRow6, pCRow6, #16

        ld1rqd  {z8.d}, p2/z, [pCRow7]
        fmla    z8.d, z17.d, alphaV0
        st1d    {z8.d}, p2, [pCRow7]
        add     pCRow7, pCRow7, #16
.endm

/******************************************************************************/
.macro INIT1x8
        fmov    d10, xzr
        fmov    d11, d10
        fmov    d12, d10
        fmov    d13, xzr
        fmov    d14, xzr
        fmov    d15, d11
        fmov    d16, d12
        fmov    d17, d15
.endm

.macro KERNEL1x8_SUB
        ld1rd   {z18.d}, p1/z, [pA]
        add     pA, pA, #8

        ld1rd   {z1.d}, p1/z, [pB]
        ld1rd   {z2.d}, p1/z, [pB, #8]
        ld1rd   {z3.d}, p1/z, [pB, #16]
        ld1rd   {z4.d}, p1/z, [pB, #24]
        ld1rd   {z5.d}, p1/z, [pB, #32]
        ld1rd   {z6.d}, p1/z, [pB, #40]
        ld1rd   {z7.d}, p1/z, [pB, #48]
        ld1rd   {z8.d}, p1/z, [pB, #56]
        add     pB, pB, #64

        fmla    z10.d, z18.d, z1.d[0]
        fmla    z11.d, z18.d, z2.d[0]
        fmla    z12.d, z18.d, z3.d[0]
        fmla    z13.d, z18.d, z4.d[0]
        fmla    z14.d, z18.d, z5.d[0]
        fmla    z15.d, z18.d, z6.d[0]
        fmla    z16.d, z18.d, z7.d[0]
        fmla    z17.d, z18.d, z8.d[0]
.endm

.macro SAVE1x8
        dup     alpha0, alpha

        ld1rd   {z0.d}, p1/z, [pCRow0]
        fmla    z0.d, z10.d, alphaV0
        st1d    {z0.d}, p1, [pCRow0]
        add     pCRow0, pCRow0, #8

        ld1rd   {z1.d}, p1/z, [pCRow1]
        fmla    z1.d, z11.d, alphaV0
        st1d    {z1.d}, p1, [pCRow1]
        add     pCRow1, pCRow1, #8

        ld1rd   {z2.d}, p1/z, [pCRow2]
        fmla    z2.d, z12.d, alphaV0
        st1d    {z2.d}, p1, [pCRow2]
        add     pCRow2, pCRow2, #8

        ld1rd   {z4.d}, p1/z, [pCRow3]
        fmla    z4.d, z13.d, alphaV0
        st1d    {z4.d}, p1, [pCRow3]
        add     pCRow3, pCRow3, #8

        ld1rd   {z5.d}, p1/z, [pCRow4]
        fmla    z5.d, z14.d, alphaV0
        st1d    {z5.d}, p1, [pCRow4]
        add     pCRow4, pCRow4, #8

        ld1rd   {z6.d}, p1/z, [pCRow5]
        fmla    z6.d, z15.d, alphaV0
        st1d    {z6.d}, p1, [pCRow5]
        add     pCRow5, pCRow5, #8

        ld1rd   {z7.d}, p1/z, [pCRow6]
        fmla    z7.d, z16.d, alphaV0
        st1d    {z7.d}, p1, [pCRow6]
        add     pCRow6, pCRow6, #8

        ld1rd   {z8.d}, p1/z, [pCRow7]
        fmla    z8.d, z17.d, alphaV0
        st1d    {z8.d}, p1, [pCRow7]
        add     pCRow7, pCRow7, #8
.endm

/******************************************************************************/
.macro INIT12x4
        fmov    d10, xzr
        fmov    d11, xzr
        fmov    d12, d10
        fmov    d13, d11
        fmov    d14, d10
        fmov    d15, d11
        fmov    d16, xzr
        fmov    d17, xzr
        fmov    d18, d16
        fmov    d19, d17
        fmov    d20, xzr
        fmov    d21, d16
.endm

.macro KERNEL12x4_I
        ld1d    {z0.d}, p0/z, [pA]
        ld1d    {z1.d}, p0/z, [pA, #1, MUL VL]
        ld1d    {z2.d}, p0/z, [pA, #2, MUL VL]
        ld1rqd  {z3.d}, p2/z, [pB]
        ld1rqd  {z4.d}, p2/z, [pB, #16]

        fmul    z10.d, z0.d, z3.d[0]
        fmul    z13.d, z0.d, z3.d[1]

        ld1d    {z5.d}, p0/z, [pA, #3, MUL VL]

        fmul    z16.d, z0.d, z4.d[0]
        fmul    z19.d, z0.d, z4.d[1]

        ld1d    {z6.d}, p0/z, [pA, #4, MUL VL]

        fmul    z11.d, z1.d, z3.d[0]
        fmul    z14.d, z1.d, z3.d[1]

        ld1d    {z7.d}, p0/z, [pA, #5, MUL VL]
        incb    pA, all, MUL #6

        fmul    z17.d, z1.d, z4.d[0]
        fmul    z20.d, z1.d, z4.d[1]

        ld1rqd  {z8.d}, p2/z, [pB, #32]

        fmul    z12.d, z2.d, z3.d[0]
        fmul    z15.d, z2.d, z3.d[1]

        ld1rqd  {z9.d}, p2/z, [pB, #48]
        incb    pB, all, MUL #2

        fmul    z18.d, z2.d, z4.d[0]
        fmul    z21.d, z2.d, z4.d[1]
.endm

.macro KERNEL12x4_M1
        ld1d    {z5.d}, p0/z, [pA]

        fmla    z10.d, z0.d, z3.d[0]
        fmla    z13.d, z0.d, z3.d[1]

        ld1d    {z6.d}, p0/z, [pA, #1, MUL VL]

        fmla    z16.d, z0.d, z4.d[0]
        fmla    z19.d, z0.d, z4.d[1]

        ld1rqd  {z8.d}, p2/z, [pB]

        fmla    z11.d, z1.d, z3.d[0]
        fmla    z14.d, z1.d, z3.d[1]

        ld1rqd  {z9.d}, p2/z, [pB, #16]
        add     pB, pB, #32

        fmla    z17.d, z1.d, z4.d[0]
        fmla    z20.d, z1.d, z4.d[1]

        ld1d    {z7.d}, p0/z, [pA, #2, MUL VL]
        incb    pA, all, MUL #3

        fmla    z12.d, z2.d, z3.d[0]
        fmla    z15.d, z2.d, z3.d[1]
        fmla    z18.d, z2.d, z4.d[0]
        fmla    z21.d, z2.d, z4.d[1]
.endm

.macro KERNEL12x4_M2
        ld1d    {z0.d}, p0/z, [pA]

        fmla    z10.d, z5.d, z8.d[0]
        fmla    z13.d, z5.d, z8.d[1]

        ld1d    {z1.d}, p0/z, [pA, #1, MUL VL]

        fmla    z16.d, z5.d, z9.d[0]
        fmla    z19.d, z5.d, z9.d[1]

        ld1rqd  {z3.d}, p2/z, [pB]

        fmla    z11.d, z6.d, z8.d[0]
        fmla    z14.d, z6.d, z8.d[1]

        ld1rqd  {z4.d}, p2/z, [pB, #16]
        add     pB, pB, #32

        fmla    z17.d, z6.d, z9.d[0]
        fmla    z20.d, z6.d, z9.d[1]

        ld1d    {z2.d}, p0/z, [pA, #2, MUL VL]
        incb    pA, all, MUL #3

        fmla    z12.d, z7.d, z8.d[0]
        fmla    z15.d, z7.d, z8.d[1]
        fmla    z18.d, z7.d, z9.d[0]
        fmla    z21.d, z7.d, z9.d[1]
.endm

.macro KERNEL12x4_E
        fmla    z10.d, z5.d, z8.d[0]
        fmla    z13.d, z5.d, z8.d[1]
        fmla    z16.d, z5.d, z9.d[0]
        fmla    z19.d, z5.d, z9.d[1]
        fmla    z11.d, z6.d, z8.d[0]
        fmla    z14.d, z6.d, z8.d[1]
        fmla    z17.d, z6.d, z9.d[0]
        fmla    z20.d, z6.d, z9.d[1]
        fmla    z12.d, z7.d, z8.d[0]
        fmla    z15.d, z7.d, z8.d[1]
        fmla    z18.d, z7.d, z9.d[0]
        fmla    z21.d, z7.d, z9.d[1]
.endm

.macro KERNEL12x4_SUB
        ld1d    {z5.d}, p0/z, [pA]
        ld1d    {z6.d}, p0/z, [pA, #1, MUL VL]

        ld1rqd  {z8.d}, p2/z, [pB]

        fmla    z10.d, z5.d, z8.d[0]
        fmla    z13.d, z5.d, z8.d[1]

        ld1d    {z7.d}, p0/z, [pA, #2, MUL VL]
        incb    pA, all, MUL #3

        fmla    z11.d, z6.d, z8.d[0]
        fmla    z14.d, z6.d, z8.d[1]

        ld1rqd  {z9.d}, p2/z, [pB, #16]
        add     pB, pB, #32

        fmla    z16.d, z5.d, z9.d[0]
        fmla    z19.d, z5.d, z9.d[1]
        fmla    z12.d, z7.d, z8.d[0]
        fmla    z15.d, z7.d, z8.d[1]
        fmla    z17.d, z6.d, z9.d[0]
        fmla    z20.d, z6.d, z9.d[1]
        fmla    z18.d, z7.d, z9.d[0]
        fmla    z21.d, z7.d, z9.d[1]
.endm

.macro SAVE12x4
        dup     alpha0, alpha

        ld1d    {z0.d}, p0/z, [pCRow0]
        ld1d    {z1.d}, p0/z, [pCRow0, #1, MUL VL]
        ld1d    {z2.d}, p0/z, [pCRow0, #2, MUL VL]
        fmla    z0.d, z10.d, alphaV0
        fmla    z1.d, z11.d, alphaV0
        fmla    z2.d, z12.d, alphaV0
        st1d    {z0.d}, p0, [pCRow0]
        st1d    {z1.d}, p0, [pCRow0, #1, MUL VL]
        st1d    {z2.d}, p0, [pCRow0, #2, MUL VL]
        incb    pCRow0, all, MUL #3

        ld1d    {z4.d}, p0/z, [pCRow1]
        ld1d    {z5.d}, p0/z, [pCRow1, #1, MUL VL]
        ld1d    {z6.d}, p0/z, [pCRow1, #2, MUL VL]
        fmla    z4.d, z13.d, alphaV0
        fmla    z5.d, z14.d, alphaV0
        fmla    z6.d, z15.d, alphaV0
        st1d    {z4.d}, p0, [pCRow1]
        st1d    {z5.d}, p0, [pCRow1, #1, MUL VL]
        st1d    {z6.d}, p0, [pCRow1, #2, MUL VL]
        incb    pCRow1, all, MUL #3

        ld1d    {z0.d}, p0/z, [pCRow2]
        ld1d    {z1.d}, p0/z, [pCRow2, #1, MUL VL]
        ld1d    {z2.d}, p0/z, [pCRow2, #2, MUL VL]
        fmla    z0.d, z16.d, alphaV0
        fmla    z1.d, z17.d, alphaV0
        fmla    z2.d, z18.d, alphaV0
        st1d    {z0.d}, p0, [pCRow2]
        st1d    {z1.d}, p0, [pCRow2, #1, MUL VL]
        st1d    {z2.d}, p0, [pCRow2, #2, MUL VL]
        incb    pCRow2, all, MUL #3

        ld1d    {z4.d}, p0/z, [pCRow3]
        ld1d    {z5.d}, p0/z, [pCRow3, #1, MUL VL]
        ld1d    {z6.d}, p0/z, [pCRow3, #2, MUL VL]
        fmla    z4.d, z19.d, alphaV0
        fmla    z5.d, z20.d, alphaV0
        fmla    z6.d, z21.d, alphaV0
        st1d    {z4.d}, p0, [pCRow3]
        st1d    {z5.d}, p0, [pCRow3, #1, MUL VL]
        st1d    {z6.d}, p0, [pCRow3, #2, MUL VL]
        incb    pCRow3, all, MUL #3
.endm

/******************************************************************************/
.macro INIT12x2
        fmov    d10,  xzr
        fmov    d11,  xzr
        fmov    d12, d10
        fmov    d13, d11
        fmov    d14, d10
        fmov    d15, d11
.endm

.macro KERNEL12x2_SUB
        ld1d    {z3.d}, p0/z, [pA]
        ld1d    {z4.d}, p0/z, [pA, #1, MUL VL]
        ld1rqd  {z7.d}, p2/z, [pB]
        add     pB, pB, #16

        fmla    z10.d, z3.d, z7.d[0]
        fmla    z13.d, z3.d, z7.d[1]

        ld1d    {z5.d}, p0/z, [pA, #2, MUL VL]
        incb    pA, all, MUL #3

        fmla    z11.d, z4.d, z7.d[0]
        fmla    z14.d, z4.d, z7.d[1]

        fmla    z12.d, z5.d, z7.d[0]
        fmla    z15.d, z5.d, z7.d[1]
.endm

.macro SAVE12x2
        dup     alpha0, alpha

        ld1d    {z0.d}, p0/z, [pCRow0]
        ld1d    {z1.d}, p0/z, [pCRow0, #1, MUL VL]
        ld1d    {z2.d}, p0/z, [pCRow0, #2, MUL VL]
        fmla    z0.d, z10.d, alphaV0
        fmla    z1.d, z11.d, alphaV0
        fmla    z2.d, z12.d, alphaV0
        st1d    {z0.d}, p0, [pCRow0]
        st1d    {z1.d}, p0, [pCRow0, #1, MUL VL]
        st1d    {z2.d}, p0, [pCRow0, #2, MUL VL]
        incb    pCRow0, all, MUL #3

        ld1d    {z4.d}, p0/z, [pCRow1]
        ld1d    {z5.d}, p0/z, [pCRow1, #1, MUL VL]
        ld1d    {z6.d}, p0/z, [pCRow1, #2, MUL VL]
        fmla    z4.d, z13.d, alphaV0
        fmla    z5.d, z14.d, alphaV0
        fmla    z6.d, z15.d, alphaV0
        st1d    {z4.d}, p0, [pCRow1]
        st1d    {z5.d}, p0, [pCRow1, #1, MUL VL]
        st1d    {z6.d}, p0, [pCRow1, #2, MUL VL]
        incb    pCRow1, all, MUL #3
.endm

/******************************************************************************/
.macro INIT12x1
        fmov    d10, xzr
        fmov    d11, xzr
        fmov    d12, d10
.endm

.macro KERNEL12x1_SUB
        ld1d    z3.d, p0/z, [pA]
        ld1d    z4.d, p0/z, [pA, #1, MUL VL]
        ld1d    z5.d, p0/z, [pA, #2, MUL VL]
        incb    pA, all, MUL #3

        ld1rd   {z7.d}, p0/z, [pB]
        add     pB, pB, #8

        fmla    z10.d, z3.d, z7.d[0]
        fmla    z11.d, z4.d, z7.d[0]
        fmla    z12.d, z5.d, z7.d[0]
.endm

.macro SAVE12x1
        dup     alpha0, alpha

        ld1d    {z0.d}, p0/z, [pCRow0]
        ld1d    {z1.d}, p0/z, [pCRow0, #1, MUL VL]
        ld1d    {z2.d}, p0/z, [pCRow0, #2, MUL VL]
        fmla    z0.d, z10.d, alphaV0
        fmla    z1.d, z11.d, alphaV0
        fmla    z2.d, z12.d, alphaV0
        st1d    {z0.d}, p0, [pCRow0]
        st1d    {z1.d}, p0, [pCRow0, #1, MUL VL]
        st1d    {z2.d}, p0, [pCRow0, #2, MUL VL]
        incb    pCRow0, all, MUL #3
.endm

/******************************************************************************/
.macro INIT8x4
        fmov    d16, xzr
        fmov    d17, xzr
        fmov    d18, xzr
        fmov    d19, d16
        fmov    d20, xzr
        fmov    d21, d16
        fmov    d22, xzr
        fmov    d23, d16
.endm

.macro KERNEL8x4_I
        ld1rd   {z2.d}, p0/z, [pB]
        ld1rd   {z3.d}, p0/z, [pB, #8]
        ld1rd   {z4.d}, p0/z, [pB, #16]
        ld1rd   {z5.d}, p0/z, [pB, #24]
        ld1d    {z0.d}, p0/z, [pA]
        ld1d    {z1.d}, p0/z, [pA, #1, MUL VL]

        fmul    z16.d, z0.d, z2.d[0]
        fmul    z18.d, z0.d, z3.d[0]

        ld1d    {z6.d}, p0/z, [pA, #2, MUL VL]
        ld1d    {z7.d}, p0/z, [pA, #3, MUL VL]
        incb    pA, all, MUL #4

        fmul    z20.d, z0.d, z4.d[0]
        fmul    z22.d, z0.d, z5.d[0]

        ld1rd   {z10.d}, p0/z, [pB, #32]
        ld1rd   {z11.d}, p0/z, [pB, #40]

        fmul    z17.d, z1.d, z2.d[0]
        fmul    z19.d, z1.d, z3.d[0]

        ld1rd   {z12.d}, p0/z, [pB, #48]
        ld1rd   {z13.d}, p0/z, [pB, #56]
        add     pB, pB, #64

        fmul    z21.d, z1.d, z4.d[0]
        fmul    z23.d, z1.d, z5.d[0]
.endm

.macro KERNEL8x4_M1
        ld1d    {z6.d}, p0/z, [pA]
        ld1d    {z7.d}, p0/z, [pA, #1, MUL VL]
        ld1rd   {z10.d}, p0/z, [pB]
        incb    pA, all, MUL #2

        fmla    z16.d, z0.d, z2.d[0]
        fmla    z17.d, z1.d, z2.d[0]

        ld1rd   {z11.d}, p0/z, [pB, #8]

        fmla    z18.d, z0.d, z3.d[0]
        fmla    z19.d, z1.d, z3.d[0]

        ld1rd   {z12.d}, p0/z, [pB, #16]

        fmla    z20.d, z0.d, z4.d[0]
        fmla    z21.d, z1.d, z4.d[0]

        ld1rd   {z13.d}, p0/z, [pB, #24]
        add     pB, pB, #32

        fmla    z22.d, z0.d, z5.d[0]
        fmla    z23.d, z1.d, z5.d[0]
.endm

.macro KERNEL8x4_M2
        ld1d    {z0.d}, p0/z, [pA]
        ld1d    {z1.d}, p0/z, [pA, #1, MUL VL]
        ld1rd   {z2.d}, p0/z, [pB]
        incb    pA, all, MUL #2

        fmla    z16.d, z6.d, z10.d[0]
        fmla    z17.d, z7.d, z10.d[0]

        ld1rd   {z3.d}, p0/z, [pB, #8]

        fmla    z18.d, z6.d, z11.d[0]
        fmla    z19.d, z7.d, z11.d[0]

        ld1rd   {z4.d}, p0/z, [pB, #16]

        fmla    z20.d, z6.d, z12.d[0]
        fmla    z21.d, z7.d, z12.d[0]

        ld1rd   {z5.d}, p0/z, [pB, #24]
        add     pB, pB, #32

        fmla    z22.d, z6.d, z13.d[0]
        fmla    z23.d, z7.d, z13.d[0]
.endm

.macro KERNEL8x4_E
        fmla    z16.d, z6.d, z10.d[0]
        fmla    z17.d, z7.d, z10.d[0]
        fmla    z18.d, z6.d, z11.d[0]
        fmla    z19.d, z7.d, z11.d[0]
        fmla    z20.d, z6.d, z12.d[0]
        fmla    z21.d, z7.d, z12.d[0]
        fmla    z22.d, z6.d, z13.d[0]
        fmla    z23.d, z7.d, z13.d[0]
.endm

.macro KERNEL8x4_SUB
        ld1d    {z6.d}, p0/z, [pA]
        ld1d    {z7.d}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        ld1rd   {z10.d}, p0/z, [pB]
        ld1rd   {z11.d}, p0/z, [pB, #8]

        fmla    z16.d, z6.d, z10.d[0]
        fmla    z17.d, z7.d, z10.d[0]

        ld1rd   {z12.d}, p0/z, [pB, #16]

        fmla    z18.d, z6.d, z11.d[0]
        fmla    z19.d, z7.d, z11.d[0]

        ld1rd   {z13.d}, p0/z, [pB, #24]
        add     pB, pB, #32

        fmla    z20.d, z6.d, z12.d[0]
        fmla    z21.d, z7.d, z12.d[0]
        fmla    z22.d, z6.d, z13.d[0]
        fmla    z23.d, z7.d, z13.d[0]
.endm

.macro SAVE8x4
        dup     alpha0, alpha

        ld1d    {z0.d}, p0/z, [pCRow0]
        ld1d    {z1.d}, p0/z, [pCRow0, #1, MUL VL]
        fmla    z0.d, z16.d, alphaV0
        fmla    z1.d, z17.d, alphaV0
        st1d    {z0.d}, p0, [pCRow0]
        st1d    {z1.d}, p0, [pCRow0, #1, MUL VL]
        incb    pCRow0, all, MUL #2

        ld1d    {z4.d}, p0/z, [pCRow1]
        ld1d    {z5.d}, p0/z, [pCRow1, #1, MUL VL]
        fmla    z4.d, z18.d, alphaV0
        fmla    z5.d, z19.d, alphaV0
        st1d    {z4.d}, p0, [pCRow1]
        st1d    {z5.d}, p0, [pCRow1, #1, MUL VL]
        incb    pCRow1, all, MUL #2

        ld1d    {z6.d}, p0/z, [pCRow2]
        ld1d    {z7.d}, p0/z, [pCRow2, #1, MUL VL]
        fmla    z6.d, z20.d, alphaV0
        fmla    z7.d, z21.d, alphaV0
        st1d    {z6.d}, p0, [pCRow2]
        st1d    {z7.d}, p0, [pCRow2, #1, MUL VL]
        incb    pCRow2, all, MUL #2

        ld1d    {z8.d}, p0/z, [pCRow3]
        ld1d    {z9.d}, p0/z, [pCRow3, #1, MUL VL]
        fmla    z8.d, z22.d, alphaV0
        fmla    z9.d, z23.d, alphaV0
        st1d    {z8.d}, p0, [pCRow3]
        st1d    {z9.d}, p0, [pCRow3, #1, MUL VL]
        incb    pCRow3, all, MUL #2
.endm

/******************************************************************************/
.macro INIT4x4
        fmov    d16, xzr
        fmov    d20, xzr
        fmov    d24, xzr
        fmov    d28, xzr
.endm

.macro KERNEL4x4_I
        ld1d    {z0.d}, p0/z, [pA]
        ld1d    {z1.d}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        ld1rqd  {z2.d}, p2/z, [pB]
        ld1rqd  {z3.d}, p2/z, [pB, #16]

        fmul    z16.d, z0.d, z2.d[0]
        fmul    z20.d, z0.d, z2.d[1]

        ld1rqd  {z4.d}, p2/z, [pB, #32]
        ld1rqd  {z5.d}, p2/z, [pB, #48]
        incb    pB, all, MUL #2

        fmul    z24.d, z0.d, z3.d[0]
        fmul    z28.d, z0.d, z3.d[1]
.endm

.macro KERNEL4x4_M1
        ld1d    {z1.d}, p0/z, [pA]
        ld1rqd  {z4.d}, p2/z, [pB]
        ld1rqd  {z5.d}, p2/z, [pB, #16]
        incb    pA, all, MUL #1
        incb    pB, all, MUL #1

        fmla    z16.d, z0.d, z2.d[0]
        fmla    z20.d, z0.d, z2.d[1]
        fmla    z24.d, z0.d, z3.d[0]
        fmla    z28.d, z0.d, z3.d[1]
.endm

.macro KERNEL4x4_M2
        ld1d    {z0.d}, p0/z, [pA]
        ld1rqd  {z2.d}, p2/z, [pB]
        ld1rqd  {z3.d}, p2/z, [pB, #16]
        incb    pA, all, MUL #1
        incb    pB, all, MUL #1

        fmla    z16.d, z1.d, z4.d[0]
        fmla    z20.d, z1.d, z4.d[1]
        fmla    z24.d, z1.d, z5.d[0]
        fmla    z28.d, z1.d, z5.d[1]
.endm

.macro KERNEL4x4_E
        fmla    z16.d, z1.d, z4.d[0]
        fmla    z20.d, z1.d, z4.d[1]
        fmla    z24.d, z1.d, z5.d[0]
        fmla    z28.d, z1.d, z5.d[1]
.endm

.macro KERNEL4x4_SUB
        ld1d    {z1.d}, p0/z, [pA]
        ld1rqd  {z4.d}, p2/z, [pB]
        ld1rqd  {z5.d}, p2/z, [pB, #16]
        incb    pA, all, MUL #1
        incb    pB, all, MUL #1

        fmla    z16.d, z1.d, z4.d[0]
        fmla    z20.d, z1.d, z4.d[1]
        fmla    z24.d, z1.d, z5.d[0]
        fmla    z28.d, z1.d, z5.d[1]
.endm

.macro SAVE4x4
        dup     alpha0, alpha

        ld1d    {z0.d}, p0/z, [pCRow0]
        fmla    z0.d, z16.d, alphaV0
        st1d    {z0.d}, p0, [pCRow0]
        add     pCRow0, pCRow0, #32

        ld1d    {z1.d}, p0/z, [pCRow1]
        fmla    z1.d, z20.d, alphaV0
        st1d    {z1.d}, p0, [pCRow1]
        add     pCRow1, pCRow1, #32

        ld1d    {z2.d}, p0/z, [pCRow2]
        fmla    z2.d, z24.d, alphaV0
        st1d    {z2.d}, p0, [pCRow2]
        add     pCRow2, pCRow2, #32

        ld1d    {z4.d}, p0/z, [pCRow3]
        fmla    z4.d, z28.d, alphaV0
        st1d    {z4.d}, p0, [pCRow3]
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
        ld1rqd    {z0.d}, p2/z, [pA]
        ld1rqd  {z1.d}, p2/z, [pB]
        ld1rqd  {z2.d}, p2/z, [pB, #16]
        add     pA, pA, #16
        add     pB, pB, #32

        fmla    z16.d, z0.d, z1.d[0]
        fmla    z20.d, z0.d, z1.d[1]
        fmla    z24.d, z0.d, z2.d[0]
        fmla    z28.d, z0.d, z2.d[1]
.endm

.macro SAVE2x4
        dup     alpha0, alpha

        ld1d    {z0.d}, p2/z, [pCRow0]
        fmla    z0.d, z16.d, alphaV0
        st1d    {z0.d}, p2, [pCRow0]
        add     pCRow0, pCRow0, #16

        ld1d    {z1.d}, p2/z, [pCRow1]
        fmla    z1.d, z20.d, alphaV0
        st1d    {z1.d}, p2, [pCRow1]
        add     pCRow1, pCRow1, #16

        ld1d    {z2.d}, p2/z, [pCRow2]
        fmla    z2.d, z24.d, alphaV0
        st1d    {z2.d}, p2, [pCRow2]
        add     pCRow2, pCRow2, #16

        ld1d    {z4.d}, p2/z, [pCRow3]
        fmla    z4.d, z28.d, alphaV0
        st1d    {z4.d}, p2, [pCRow3]
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
       ld1rd   z0.d, p1/z, [pA]
       add     pA, pA, #8

       ld1rd   {z8.d}, p1/z, [pB]
       ld1rd   {z9.d}, p1/z, [pB, #8]
       ld1rd   {z10.d}, p1/z, [pB, #16]
       ld1rd   {z11.d}, p1/z, [pB, #24]
       add     pB, pB, #32

       fmla    z16.d,  z8.d, z0.d[0]
       fmla    z17.d,  z9.d, z0.d[0]
       fmla    z18.d,  z10.d, z0.d[0]
       fmla    z19.d,  z11.d, z0.d[0]
.endm

.macro SAVE1x4
        dup     alpha0, alpha

        ld1rd   z8.d, p1/z, [pCRow0]
        fmla    z8.d,  z16.d, alphaV0
        st1d    {z8.d}, p1, [pCRow0]
        add     pCRow0, pCRow0, #8

        ld1rd   z9.d, p1/z, [pCRow1]
        fmla    z9.d,  z17.d, alphaV0
        st1d    {z9.d}, p1, [pCRow1]
        add     pCRow1, pCRow1, #8


        ld1rd   z12.d, p1/z, [pCRow2]
        fmla    z12.d, z18.d, alphaV0
        st1d    {z12.d}, p1, [pCRow2]
        add     pCRow2, pCRow2, #8

        ld1rd   z13.d, p1/z, [pCRow3]
        fmla    z13.d, z19.d, alphaV0
        st1d    {z13.d}, p1, [pCRow3]
        add     pCRow3, pCRow3, #8
.endm

/******************************************************************************/
.macro INIT8x2
        fmov    d16, xzr
        fmov    d17, d16
        fmov    d20, d17
        fmov    d21, d16
.endm

.macro KERNEL8x2_SUB
        ld1d    {z0.d}, p0/z, [pA]
        ld1d    {z1.d}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        ld1rqd  {z4.d}, p2/z, [pB]
        add     pB, pB, #16

        fmla    z16.d, z0.d, z4.d[0]
        fmla    z17.d, z1.d, z4.d[0]
        fmla    z20.d, z0.d, z4.d[1]
        fmla    z21.d, z1.d, z4.d[1]
.endm

.macro SAVE8x2
        dup     alpha0, alpha

        ld1d    z0.d, p0/z, [pCRow0]
        fmla    z0.d,  z16.d, alphaV0
        st1d    {z0.d}, p0, [pCRow0]

        ld1d    z1.d, p0/z, [pCRow0, #1, MUL VL]
        fmla    z1.d,  z17.d, alphaV0
        st1d    {z1.d}, p0, [pCRow0, #1, MUL VL]

        add     pCRow0, pCRow0, #64

        ld1d    z4.d, p0/z, [pCRow1]
        fmla    z4.d,  z20.d, alphaV0
        st1d    {z4.d}, p0, [pCRow1]

        ld1d    z5.d, p0/z, [pCRow1, #1, MUL VL]
        fmla    z5.d,  z21.d, alphaV0
        st1d    {z5.d}, p0, [pCRow1, #1, MUL VL]
        add     pCRow1, pCRow1, #64
.endm

/******************************************************************************/
.macro INIT4x2
        fmov    d16, xzr
        fmov    d17, d16
.endm

.macro KERNEL4x2_SUB
        ld1d    {z0.d}, p0/z, [pA]
        incb    pA, all, MUL #1

        ld1rqd  {z4.d}, p2/z, [pB]
        add     pB, pB, #16

        fmla    z16.d, z0.d, z4.d[0]
        fmla    z17.d, z0.d, z4.d[1]
.endm

.macro SAVE4x2
        dup     alpha0, alpha

        ld1d    z8.d, p0/z, [pCRow0]
        fmla    z8.d, z16.d, alphaV0
        st1d    {z8.d}, p0, [pCRow0]
        add     pCRow0, pCRow0, #32

        ld1d    z9.d, p0/z, [pCRow1]
        fmla    z9.d, z17.d, alphaV0
        st1d    {z9.d}, p0, [pCRow1]
        add     pCRow1, pCRow1, #32
.endm

/******************************************************************************/
.macro INIT2x2
        fmov    d16, xzr
        fmov    d20, d16
.endm

.macro KERNEL2x2_SUB
        ld1rqd  {z0.d}, p2/z, [pA]
        add     pA, pA, #16

        ld1rqd  {z4.d}, p2/z, [pB]
        add     pB, pB, #16

        fmla    z16.d, z0.d, z4.d[0]
        fmla    z20.d, z0.d, z4.d[1]
.endm

.macro SAVE2x2
        dup     alpha0, alpha

        ld1rqd  z8.d, p2/z, [pCRow0]
        fmla    z8.d, z16.d, alphaV0
        st1d    {z8.d}, p2, [pCRow0]
        add     pCRow0, pCRow0, #16

        ld1rqd  z12.d, p2/z, [pCRow1]
        fmla    z12.d, z20.d, alphaV0
        st1d    {z12.d}, p2, [pCRow1]
        add     pCRow1, pCRow1, #16
.endm

/******************************************************************************/
.macro INIT1x2
        fmov    d16, xzr
        fmov    d17, xzr
.endm

.macro KERNEL1x2_SUB
        ld1rd   z8.d, p1/z, [pB]
        ld1rd   z9.d, p1/z, [pB, #8]
        add     pB , pB, #16

        ld1rd   z0.d, p1/z, [pA]
        add     pA, pA, #8

        fmla    z16.d, z0.d, z8.d[0]
        fmla    z17.d, z0.d, z9.d[0]
.endm

.macro SAVE1x2
        dup     alpha0,   alpha

        ld1rd   z8.d, p1/z, [pCRow0]
        ld1rd   z9.d, p1/z, [pCRow1]

        fmla    z8.d,  z16.d, alphaV0
        fmla    z9.d,  z17.d, alphaV0

        st1d    {z8.d}, p1, [pCRow0]
        st1d    {z9.d}, p1, [pCRow1]

        add     pCRow0, pCRow0, #8
        add     pCRow1, pCRow1, #8
.endm

/******************************************************************************/
.macro INIT8x1
        fmov    d16, xzr
        fmov    d17, xzr
.endm

.macro KERNEL8x1_SUB
        ld1d    z0.d, p0/z, [pA]
        ld1d    z1.d, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        ld1rd   z3.d, p0/z, [pB]
        add     pB , pB, #8

        fmla    z16.d, z0.d, z3.d[0]
        fmla    z17.d, z1.d, z3.d[0]
.endm

.macro SAVE8x1
        dup     alpha0, alpha

        ld1d    z8.d, p0/z, [pCRow0]
        ld1d    z9.d, p0/z, [pCRow0, #1, MUL VL]
        fmla    z8.d,  z16.d, alphaV0
        fmla    z9.d,  z17.d, alphaV0
        st1d    {z8.d}, p0, [pCRow0]
        st1d    {z9.d}, p0, [pCRow0, #1, MUL VL]

        add     pCRow0, pCRow0, #64
.endm

/******************************************************************************/
.macro INIT4x1
        fmov    d16, xzr
.endm

.macro KERNEL4x1_SUB
        ld1d    z0.d, p0/z, [pA]
        incb    pA, all, MUL #1

        ld1rd   z1.d, p0/z, [pB]
        add     pB , pB, #8

        fmla    z16.d, z0.d, z1.d[0]
.endm

.macro SAVE4x1
        dup     alpha0, alpha

        ld1d    z8.d, p0/z, [pCRow0]
        fmla    z8.d,  z16.d, alphaV0
        st1d    {z8.d}, p0, [pCRow0]

        add     pCRow0, pCRow0, #32
.endm

/******************************************************************************/
.macro INIT2x1
        fmov    d16, xzr
.endm

.macro KERNEL2x1_SUB
        ld1rqd  z0.d, p2/z, [pA]
        ld1rd   z8.d, p1/z, [pB]

        add     pA , pA, #16
        add     pB , pB, #8

        fmla    z16.d,  z0.d, z8.d[0]
.endm

.macro SAVE2x1
        dup     alpha0, alpha

        ld1rqd  z8.d, p2/z, [pCRow0]
        fmla    z8.d, z16.d, alphaV0

        st1d    {z8.d}, p2, [pCRow0]
        add     pCRow0, pCRow0, #16
.endm

/******************************************************************************/
.macro INIT1x1
        fmov    d16, xzr
.endm

.macro KERNEL1x1_SUB
        ld1rd   z0.d, p1/z, [pA]
        add     pA , pA, #8

        ld1rd   z8.d, p1/z, [pB]
        add     pB , pB, #8

        fmla    z16.d, p1/m, z0.d, z8.d
.endm

.macro SAVE1x1
        dup     alpha0, alpha

        ld1rd   z8.d, p1/z, [pCRow0]
        fmla    z8.d, p1/m, alpha0, z16.d
        st1d    {z8.d}, p1, [pCRow0]

        add     pCRow0, pCRow0, #8
.endm

/*******************************************************************************
* End of macro definitions
*******************************************************************************/

        PROLOGUE

.Ldgemm_kernel_begin:

        .align 4
        add     sp, sp, #-(11 * 16)
        stp     d8, d9, [sp, #(0 * 16)]
        stp     d10, d11, [sp, #(1 * 16)]
        stp     d12, d13, [sp, #(2 * 16)]
        stp     d14, d15, [sp, #(3 * 16)]
        stp     d16, d17, [sp, #(4 * 16)]
        stp     x18, x19, [sp, #(5 * 16)]
        stp     x20, x21, [sp, #(6 * 16)]
        stp     x22, x23, [sp, #(7 * 16)]
        stp     x24, x25, [sp, #(8 * 16)]
        stp     x26, x27, [sp, #(9 * 16)]
        str     x28, [sp, #(10 * 16)]

        prfm    PLDL1KEEP, [origPB]
        prfm    PLDL1KEEP, [origPA]

        fmov    alpha, d0

        lsl     LDC, LDC, #3            // ldc = ldc * 8

        mov     pB, origPB
        
        ptrue   p0.d, all
        ptrue   p1.d, VL1
        ptrue   p2.d, VL2

        mov     counterJ, origN
        asr     counterJ, counterJ, #3      // J = J / 8
        cmp     counterJ, #0
        ble     .Ldgemm_kernel_L4_BEGIN

/******************************************************************************/

.Ldgemm_kernel_L8_BEGIN:
        mov     pCRow0, pC
        add     pCRow1, pCRow0, LDC
        add     pCRow2, pCRow1, LDC
        add     pCRow3, pCRow2, LDC
        add     pCRow4, pCRow3, LDC
        add     pCRow5, pCRow4, LDC
        add     pCRow6, pCRow5, LDC
        add     pCRow7, pCRow6, LDC

        add     pC, pCRow7, LDC

        mov     pA, origPA          // pA = start of A array

.Ldgemm_kernel_L8_M12_BEGIN:
        

        mov     counterI, origM
        cmp     counterI, #11
        ble     .Ldgemm_kernel_L8_M8_BEGIN

        .align 4
.Ldgemm_kernel_L8_M12_20:

        mov     pB, origPB

        asr     counterL, origK, #3
        cmp     counterL, #2
        blt     .Ldgemm_kernel_L8_M12_32

        KERNEL12x8_I
        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_M1

        subs    counterL, counterL, #2
        ble     .Ldgemm_kernel_L8_M12_22a

        .align 4
.Ldgemm_kernel_L8_M12_22:

        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_M1

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L8_M12_22

        .align 4
.Ldgemm_kernel_L8_M12_22a:

        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_E

        b       .Ldgemm_kernel_L8_M12_44

        .align 4
.Ldgemm_kernel_L8_M12_32:

        tst     counterL, #1
        ble     .Ldgemm_kernel_L8_M12_40

        KERNEL12x8_I
        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_M1
        KERNEL12x8_E

        b       .Ldgemm_kernel_L8_M12_44

.Ldgemm_kernel_L8_M12_40:

        INIT12x8

.Ldgemm_kernel_L8_M12_44:

        ands    counterL , origK, #7
        ble     .Ldgemm_kernel_L8_M12_100

        .align 4
.Ldgemm_kernel_L8_M12_46:

        KERNEL12x8_SUB
        subs    counterL, counterL, #1
        bne     .Ldgemm_kernel_L8_M12_46

.Ldgemm_kernel_L8_M12_100:
        prfm    PLDL1KEEP, [pA]
        prfm    PLDL1KEEP, [pA, #64]
        prfm    PLDL1KEEP, [origPB]

        SAVE12x8

.Ldgemm_kernel_L8_M12_END:
        subs    counterI, counterI, #12
        cmp     counterI, #11
        bgt     .Ldgemm_kernel_L8_M12_20

//------------------------------------------------------------------------------

.Ldgemm_kernel_L8_M8_BEGIN:

        cmp     counterI, #1
        blt     .Ldgemm_kernel_L8_END

        cmp     counterI, #8
        blt     .Ldgemm_kernel_L8_M4_BEGIN

.Ldgemm_kernel_L8_M8_20:

        mov     pB, origPB

        asr     counterL , origK, #1        // L = K / 2
        cmp     counterL , #2            // is there at least 4 to do?
        blt     .Ldgemm_kernel_L8_M8_32

        KERNEL8x8_I                // do one in the K
        KERNEL8x8_M2                // do another in the K

        subs    counterL, counterL, #2
        ble     .Ldgemm_kernel_L8_M8_22a
        .align 4

.Ldgemm_kernel_L8_M8_22:

        KERNEL8x8_M1
        KERNEL8x8_M2

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L8_M8_22

.Ldgemm_kernel_L8_M8_22a:

        KERNEL8x8_M1
        KERNEL8x8_E

        b       .Ldgemm_kernel_L8_M8_44

.Ldgemm_kernel_L8_M8_32:

        tst     counterL, #1
        ble     .Ldgemm_kernel_L8_M8_40

        KERNEL8x8_I
        KERNEL8x8_E

        b       .Ldgemm_kernel_L8_M8_44

.Ldgemm_kernel_L8_M8_40:

        INIT8x8

.Ldgemm_kernel_L8_M8_44:

        ands    counterL , origK, #1
        ble     .Ldgemm_kernel_L8_M8_100

.Ldgemm_kernel_L8_M8_46:

        KERNEL8x8_SUB

.Ldgemm_kernel_L8_M8_100:

        SAVE8x8

.Ldgemm_kernel_L8_M8_END:

        subs    counterI, counterI, #8

//------------------------------------------------------------------------------

.Ldgemm_kernel_L8_M4_BEGIN:

        cmp     counterI ,#1
        blt     .Ldgemm_kernel_L8_END

        cmp     counterI, #4
        blt     .Ldgemm_kernel_L8_M2_BEGIN

.Ldgemm_kernel_L8_M4_20:

        mov     pB, origPB

        asr     counterL , origK, #1        // L = K / 2
        cmp     counterL , #2            // is there at least 4 to do?
        blt     .Ldgemm_kernel_L8_M4_32

        KERNEL4x8_I                // do one in the K
        KERNEL4x8_M2                // do another in the K

        subs    counterL, counterL, #2
        ble     .Ldgemm_kernel_L8_M4_22a
        .align 4

.Ldgemm_kernel_L8_M4_22:

        KERNEL4x8_M1
        KERNEL4x8_M2

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L8_M4_22

.Ldgemm_kernel_L8_M4_22a:

        KERNEL4x8_M1
        KERNEL4x8_E

        b       .Ldgemm_kernel_L8_M4_44

.Ldgemm_kernel_L8_M4_32:

        tst     counterL, #1
        ble     .Ldgemm_kernel_L8_M4_40

        KERNEL4x8_I
        KERNEL4x8_E

        b       .Ldgemm_kernel_L8_M4_44

.Ldgemm_kernel_L8_M4_40:

        INIT4x8

.Ldgemm_kernel_L8_M4_44:

        ands    counterL , origK, #1
        ble     .Ldgemm_kernel_L8_M4_100

.Ldgemm_kernel_L8_M4_46:

        KERNEL4x8_SUB

.Ldgemm_kernel_L8_M4_100:

        SAVE4x8

.Ldgemm_kernel_L8_M4_END:

        subs    counterI, counterI, #4

//------------------------------------------------------------------------------

.Ldgemm_kernel_L8_M2_BEGIN:

        cmp     counterI, #1
        blt     .Ldgemm_kernel_L8_END
        cmp     counterI, #2            // counterI = counterI / 2
        blt     .Ldgemm_kernel_L8_M1_BEGIN

.Ldgemm_kernel_L8_M2_20:

        INIT2x8

        mov     pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp     counterL , #0
        ble     .Ldgemm_kernel_L8_M2_40

.Ldgemm_kernel_L8_M2_22:

        KERNEL2x8_SUB
        KERNEL2x8_SUB
        KERNEL2x8_SUB
        KERNEL2x8_SUB

        KERNEL2x8_SUB
        KERNEL2x8_SUB
        KERNEL2x8_SUB
        KERNEL2x8_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L8_M2_22


.Ldgemm_kernel_L8_M2_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble     .Ldgemm_kernel_L8_M2_100

.Ldgemm_kernel_L8_M2_42:

        KERNEL2x8_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L8_M2_42

.Ldgemm_kernel_L8_M2_100:

        SAVE2x8

.Ldgemm_kernel_L8_M2_END:

        subs    counterI, counterI, #2


.Ldgemm_kernel_L8_M1_BEGIN:

        cmp     counterI, #1            // counterI = counterI % 2
        blt     .Ldgemm_kernel_L8_END

.Ldgemm_kernel_L8_M1_20:

        INIT1x8

        mov     pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp     counterL , #0
        ble     .Ldgemm_kernel_L8_M1_40

.Ldgemm_kernel_L8_M1_22:

        KERNEL1x8_SUB
        KERNEL1x8_SUB
        KERNEL1x8_SUB
        KERNEL1x8_SUB

        KERNEL1x8_SUB
        KERNEL1x8_SUB
        KERNEL1x8_SUB
        KERNEL1x8_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L8_M1_22

.Ldgemm_kernel_L8_M1_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble     .Ldgemm_kernel_L8_M1_100

.Ldgemm_kernel_L8_M1_42:

        KERNEL1x8_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L8_M1_42

.Ldgemm_kernel_L8_M1_100:

        SAVE1x8

.Ldgemm_kernel_L8_END:

        add     origPB, origPB, origK, lsl #6    // B = B + K * 8 * 8

        subs    counterJ, counterJ , #1        // j--
        bgt     .Ldgemm_kernel_L8_BEGIN


/******************************************************************************/
/******************************************************************************/

.Ldgemm_kernel_L4_BEGIN:

        mov     counterJ , origN
        tst     counterJ , #7
        ble     .Ldgemm_kernel_L999

        tst     counterJ , #4
        ble     .Ldgemm_kernel_L2_BEGIN

        mov     pCRow0, pC
        add     pCRow1, pCRow0, LDC
        add     pCRow2, pCRow1, LDC
        add     pCRow3, pCRow2, LDC

        add     pC, pCRow3, LDC

        mov     pA, origPA            // pA = start of A array

.Ldgemm_kernel_L4_M12_BEGIN:

        mov     counterI, origM
        cmp     counterI, #11
        ble     .Ldgemm_kernel_L4_M8_BEGIN

        .align 4
.Ldgemm_kernel_L4_M12_20:

        mov     pB, origPB

        asr     counterL , origK, #3
        cmp     counterL , #2
        blt     .Ldgemm_kernel_L4_M12_32

        KERNEL12x4_I
        KERNEL12x4_M2
        KERNEL12x4_M1
        KERNEL12x4_M2
        KERNEL12x4_M1
        KERNEL12x4_M2
        KERNEL12x4_M1
        KERNEL12x4_M2

        subs    counterL, counterL, #2
        ble     .Ldgemm_kernel_L4_M12_22a

        .align 4
.Ldgemm_kernel_L4_M12_22:

        KERNEL12x4_M1
        KERNEL12x4_M2
        KERNEL12x4_M1
        KERNEL12x4_M2
        KERNEL12x4_M1
        KERNEL12x4_M2
        KERNEL12x4_M1
        KERNEL12x4_M2

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L4_M12_22

        .align 4
.Ldgemm_kernel_L4_M12_22a:

        KERNEL12x4_M1
        KERNEL12x4_M2
        KERNEL12x4_M1
        KERNEL12x4_M2
        KERNEL12x4_M1
        KERNEL12x4_M2
        KERNEL12x4_M1
        KERNEL12x4_E

        b    .Ldgemm_kernel_L4_M12_44

        .align 4
.Ldgemm_kernel_L4_M12_32:

        tst     counterL, #1
        ble     .Ldgemm_kernel_L4_M12_40

        KERNEL12x4_I
        KERNEL12x4_M2
        KERNEL12x4_M1
        KERNEL12x4_M2
        KERNEL12x4_M1
        KERNEL12x4_M2
        KERNEL12x4_M1
        KERNEL12x4_E

        b       .Ldgemm_kernel_L4_M12_44

.Ldgemm_kernel_L4_M12_40:

        INIT12x4

.Ldgemm_kernel_L4_M12_44:

        ands    counterL , origK, #7
        ble     .Ldgemm_kernel_L4_M12_100

        .align 4
.Ldgemm_kernel_L4_M12_46:

        KERNEL12x4_SUB
        subs    counterL, counterL, #1
        bne     .Ldgemm_kernel_L4_M12_46

.Ldgemm_kernel_L4_M12_100:

        prfm    PLDL1KEEP, [pA]
        prfm    PLDL1KEEP, [pA, #64]
        prfm    PLDL1KEEP, [origPB]

        SAVE12x4

.Ldgemm_kernel_L4_M12_END:

        subs    counterI, counterI, #12
        cmp     counterI, #11
        bgt     .Ldgemm_kernel_L4_M12_20

//------------------------------------------------------------------------------

.Ldgemm_kernel_L4_M8_BEGIN:

        cmp     counterI, #1
        blt     .Ldgemm_kernel_L4_END

        cmp     counterI, #8
        blt     .Ldgemm_kernel_L4_M4_BEGIN

.Ldgemm_kernel_L4_M8_20:

        mov     pB, origPB

        asr     counterL , origK, #1        // L = K / 2
        cmp     counterL , #2            // is there at least 4 to do?
        blt     .Ldgemm_kernel_L4_M8_32

        KERNEL8x4_I                // do one in the K
        KERNEL8x4_M2                // do another in the K

        subs    counterL, counterL, #2
        ble     .Ldgemm_kernel_L4_M8_22a
        .align 4

.Ldgemm_kernel_L4_M8_22:

        KERNEL8x4_M1
        KERNEL8x4_M2

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L4_M8_22

.Ldgemm_kernel_L4_M8_22a:

        KERNEL8x4_M1
        KERNEL8x4_E

        b       .Ldgemm_kernel_L4_M8_44

.Ldgemm_kernel_L4_M8_32:

        tst     counterL, #1
        ble     .Ldgemm_kernel_L4_M8_40

        KERNEL8x4_I
        KERNEL8x4_E

        b       .Ldgemm_kernel_L4_M8_44

.Ldgemm_kernel_L4_M8_40:

        INIT8x4

.Ldgemm_kernel_L4_M8_44:

        ands    counterL , origK, #1
        ble     .Ldgemm_kernel_L4_M8_100

.Ldgemm_kernel_L4_M8_46:

        KERNEL8x4_SUB

.Ldgemm_kernel_L4_M8_100:

        SAVE8x4

.Ldgemm_kernel_L4_M8_END:

        subs    counterI, counterI, #8

//------------------------------------------------------------------------------

.Ldgemm_kernel_L4_M4_BEGIN:

        cmp     counterI , #1
        blt     .Ldgemm_kernel_L4_END

        cmp     counterI, #4
        blt     .Ldgemm_kernel_L4_M2_BEGIN

.Ldgemm_kernel_L4_M4_20:

        mov     pB, origPB

        asr     counterL , origK, #1        // L = K / 2
        cmp     counterL , #2            // is there at least 4 to do?
        blt     .Ldgemm_kernel_L4_M4_32

        KERNEL4x4_I                // do one in the K
        KERNEL4x4_M2                // do another in the K

        subs    counterL, counterL, #2
        ble     .Ldgemm_kernel_L4_M4_22a
        .align 4

.Ldgemm_kernel_L4_M4_22:

        KERNEL4x4_M1
        KERNEL4x4_M2

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L4_M4_22

.Ldgemm_kernel_L4_M4_22a:

        KERNEL4x4_M1
        KERNEL4x4_E

        b       .Ldgemm_kernel_L4_M4_44

.Ldgemm_kernel_L4_M4_32:

        tst     counterL, #1
        ble     .Ldgemm_kernel_L4_M4_40

        KERNEL4x4_I
        KERNEL4x4_E

        b       .Ldgemm_kernel_L4_M4_44

.Ldgemm_kernel_L4_M4_40:

        INIT4x4

.Ldgemm_kernel_L4_M4_44:

        ands    counterL , origK, #1
        ble     .Ldgemm_kernel_L4_M4_100

.Ldgemm_kernel_L4_M4_46:

        KERNEL4x4_SUB

.Ldgemm_kernel_L4_M4_100:

        SAVE4x4

.Ldgemm_kernel_L4_M4_END:

        subs    counterI, counterI, #4

//------------------------------------------------------------------------------

.Ldgemm_kernel_L4_M2_BEGIN:

        cmp     counterI , #1
        blt     .Ldgemm_kernel_L4_END

        cmp     counterI, #2            // counterI = counterI / 2
        blt     .Ldgemm_kernel_L4_M1_BEGIN

.Ldgemm_kernel_L4_M2_20:

        INIT2x4

        mov     pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp     counterL , #0
        ble     .Ldgemm_kernel_L4_M2_40

.Ldgemm_kernel_L4_M2_22:

        KERNEL2x4_SUB
        KERNEL2x4_SUB
        KERNEL2x4_SUB
        KERNEL2x4_SUB

        KERNEL2x4_SUB
        KERNEL2x4_SUB
        KERNEL2x4_SUB
        KERNEL2x4_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L4_M2_22

.Ldgemm_kernel_L4_M2_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble     .Ldgemm_kernel_L4_M2_100

.Ldgemm_kernel_L4_M2_42:

        KERNEL2x4_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L4_M2_42

.Ldgemm_kernel_L4_M2_100:

        SAVE2x4

.Ldgemm_kernel_L4_M2_END:

        subs    counterI, counterI, #2

.Ldgemm_kernel_L4_M1_BEGIN:

        cmp     counterI, #1            // counterI = counterI % 2
        blt     .Ldgemm_kernel_L4_END

.Ldgemm_kernel_L4_M1_20:

        INIT1x4

        mov     pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp     counterL , #0
        ble     .Ldgemm_kernel_L4_M1_40

.Ldgemm_kernel_L4_M1_22:

        KERNEL1x4_SUB
        KERNEL1x4_SUB
        KERNEL1x4_SUB
        KERNEL1x4_SUB

        KERNEL1x4_SUB
        KERNEL1x4_SUB
        KERNEL1x4_SUB
        KERNEL1x4_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L4_M1_22

.Ldgemm_kernel_L4_M1_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble     .Ldgemm_kernel_L4_M1_100

.Ldgemm_kernel_L4_M1_42:

        KERNEL1x4_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L4_M1_42

.Ldgemm_kernel_L4_M1_100:

        SAVE1x4

.Ldgemm_kernel_L4_END:

        add     origPB, origPB, origK, lsl #5    // B = B + K * 8 * 4


/******************************************************************************/

.Ldgemm_kernel_L2_BEGIN:   // less than 2 left in N direction

        mov     counterJ , origN
        tst     counterJ , #3
        ble     .Ldgemm_kernel_L999

        tst     counterJ , #2
        ble     .Ldgemm_kernel_L1_BEGIN

        mov     pCRow0, pC            // pCRow0 = pC
        add     pCRow1, pCRow0, LDC
        add     pC, pCRow1, LDC

        mov     pA, origPA            // pA = A

.Ldgemm_kernel_L2_M12_BEGIN:

        mov     counterI, origM

        cmp     counterI, #12
        blt     .Ldgemm_kernel_L2_M8_BEGIN

.Ldgemm_kernel_L2_M12_20:

        INIT12x2

        mov     pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp     counterL,#0
        ble     .Ldgemm_kernel_L2_M12_40
        .align 4

.Ldgemm_kernel_L2_M12_22:

        KERNEL12x2_SUB
        KERNEL12x2_SUB
        KERNEL12x2_SUB
        KERNEL12x2_SUB

        KERNEL12x2_SUB
        KERNEL12x2_SUB
        KERNEL12x2_SUB
        KERNEL12x2_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L2_M12_22

.Ldgemm_kernel_L2_M12_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble     .Ldgemm_kernel_L2_M12_100

.Ldgemm_kernel_L2_M12_42:

        KERNEL12x2_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L2_M12_42

.Ldgemm_kernel_L2_M12_100:

        SAVE12x2

.Ldgemm_kernel_L2_M12_END:

        subs    counterI, counterI, #12

        cmp     counterI, #11
        bgt     .Ldgemm_kernel_L2_M12_20

//------------------------------------------------------------------------------

.Ldgemm_kernel_L2_M8_BEGIN:

        cmp     counterI, #0
        ble     .Ldgemm_kernel_L2_END

        cmp     counterI, #7
        ble     .Ldgemm_kernel_L2_M4_BEGIN

.Ldgemm_kernel_L2_M8_20:

        INIT8x2

        mov     pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp     counterL,#0
        ble     .Ldgemm_kernel_L2_M8_40
        .align 4

.Ldgemm_kernel_L2_M8_22:

        KERNEL8x2_SUB
        KERNEL8x2_SUB
        KERNEL8x2_SUB
        KERNEL8x2_SUB

        KERNEL8x2_SUB
        KERNEL8x2_SUB
        KERNEL8x2_SUB
        KERNEL8x2_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L2_M8_22

.Ldgemm_kernel_L2_M8_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble     .Ldgemm_kernel_L2_M8_100

.Ldgemm_kernel_L2_M8_42:

        KERNEL8x2_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L2_M8_42

.Ldgemm_kernel_L2_M8_100:

        SAVE8x2

.Ldgemm_kernel_L2_M8_END:

        subs    counterI, counterI, #8

//------------------------------------------------------------------------------

.Ldgemm_kernel_L2_M4_BEGIN:

        cmp     counterI, #1
        blt     .Ldgemm_kernel_L2_END

        cmp     counterI, #4
        blt     .Ldgemm_kernel_L2_M2_BEGIN

.Ldgemm_kernel_L2_M4_20:

        INIT4x2

        mov     pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp     counterL,#0
        ble     .Ldgemm_kernel_L2_M4_40
        .align 4

.Ldgemm_kernel_L2_M4_22:

        KERNEL4x2_SUB
        KERNEL4x2_SUB
        KERNEL4x2_SUB
        KERNEL4x2_SUB

        KERNEL4x2_SUB
        KERNEL4x2_SUB
        KERNEL4x2_SUB
        KERNEL4x2_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L2_M4_22

.Ldgemm_kernel_L2_M4_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble     .Ldgemm_kernel_L2_M4_100

.Ldgemm_kernel_L2_M4_42:

        KERNEL4x2_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L2_M4_42

.Ldgemm_kernel_L2_M4_100:

        SAVE4x2

.Ldgemm_kernel_L2_M4_END:

        subs    counterI, counterI, #4

//------------------------------------------------------------------------------

.Ldgemm_kernel_L2_M2_BEGIN:

        cmp     counterI , #1
        blt     .Ldgemm_kernel_L2_END

        cmp     counterI, #2            // counterI = counterI / 2
        blt     .Ldgemm_kernel_L2_M1_BEGIN

.Ldgemm_kernel_L2_M2_20:

        INIT2x2

        mov     pB, origPB

        asr     counterL, origK, #3        // counterL = counterL / 8
        cmp     counterL, #0
        ble     .Ldgemm_kernel_L2_M2_40

.Ldgemm_kernel_L2_M2_22:

        KERNEL2x2_SUB
        KERNEL2x2_SUB
        KERNEL2x2_SUB
        KERNEL2x2_SUB

        KERNEL2x2_SUB
        KERNEL2x2_SUB
        KERNEL2x2_SUB
        KERNEL2x2_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L2_M2_22

.Ldgemm_kernel_L2_M2_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble     .Ldgemm_kernel_L2_M2_100

.Ldgemm_kernel_L2_M2_42:

        KERNEL2x2_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L2_M2_42

.Ldgemm_kernel_L2_M2_100:

        SAVE2x2

.Ldgemm_kernel_L2_M2_END:

        subs    counterI, counterI, #2

.Ldgemm_kernel_L2_M1_BEGIN:

        cmp     counterI, #1            // counterI = counterI % 2
        blt     .Ldgemm_kernel_L2_END

.Ldgemm_kernel_L2_M1_20:

        INIT1x2

        mov     pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp     counterL, #0
        ble     .Ldgemm_kernel_L2_M1_40

.Ldgemm_kernel_L2_M1_22:

        KERNEL1x2_SUB
        KERNEL1x2_SUB
        KERNEL1x2_SUB
        KERNEL1x2_SUB

        KERNEL1x2_SUB
        KERNEL1x2_SUB
        KERNEL1x2_SUB
        KERNEL1x2_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L2_M1_22

.Ldgemm_kernel_L2_M1_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble     .Ldgemm_kernel_L2_M1_100

.Ldgemm_kernel_L2_M1_42:

        KERNEL1x2_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L2_M1_42

.Ldgemm_kernel_L2_M1_100:

        SAVE1x2

.Ldgemm_kernel_L2_END:

        add     origPB, origPB, origK, lsl #4    // B = B + K * 2 * 8

/******************************************************************************/

.Ldgemm_kernel_L1_BEGIN:

        mov     counterJ , origN
        tst     counterJ , #1
        ble     .Ldgemm_kernel_L999 // done


        mov     pCRow0, pC            // pCRow0 = C
        add     pC , pC , LDC            // Update pC to point to next

        mov     pA, origPA            // pA = A

.Ldgemm_kernel_L1_M12_BEGIN:

        mov     counterI, origM

        cmp     counterI, #12
        blt     .Ldgemm_kernel_L1_M8_BEGIN

.Ldgemm_kernel_L1_M12_20:

        INIT12x1

        mov     pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp     counterL , #0
        ble     .Ldgemm_kernel_L1_M12_40
        .align 4

.Ldgemm_kernel_L1_M12_22:

        KERNEL12x1_SUB
        KERNEL12x1_SUB
        KERNEL12x1_SUB
        KERNEL12x1_SUB

        KERNEL12x1_SUB
        KERNEL12x1_SUB
        KERNEL12x1_SUB
        KERNEL12x1_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L1_M12_22


.Ldgemm_kernel_L1_M12_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble     .Ldgemm_kernel_L1_M12_100

.Ldgemm_kernel_L1_M12_42:

        KERNEL12x1_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L1_M12_42

.Ldgemm_kernel_L1_M12_100:

        SAVE12x1

.Ldgemm_kernel_L1_M12_END:

        subs    counterI, counterI, #12

        cmp     counterI, #11
        bgt     .Ldgemm_kernel_L1_M12_20

//------------------------------------------------------------------------------

.Ldgemm_kernel_L1_M8_BEGIN:

        cmp     counterI, #1
        blt     .Ldgemm_kernel_L1_END

        cmp     counterI, #8
        blt     .Ldgemm_kernel_L1_M4_BEGIN

.Ldgemm_kernel_L1_M8_20:

        INIT8x1

        mov     pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp     counterL , #0
        ble     .Ldgemm_kernel_L1_M8_40
        .align 4

.Ldgemm_kernel_L1_M8_22:

        KERNEL8x1_SUB
        KERNEL8x1_SUB
        KERNEL8x1_SUB
        KERNEL8x1_SUB

        KERNEL8x1_SUB
        KERNEL8x1_SUB
        KERNEL8x1_SUB
        KERNEL8x1_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L1_M8_22

.Ldgemm_kernel_L1_M8_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble     .Ldgemm_kernel_L1_M8_100

.Ldgemm_kernel_L1_M8_42:

        KERNEL8x1_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L1_M8_42

.Ldgemm_kernel_L1_M8_100:

        SAVE8x1

.Ldgemm_kernel_L1_M8_END:

        subs    counterI, counterI, #8

//------------------------------------------------------------------------------

.Ldgemm_kernel_L1_M4_BEGIN:

        cmp     counterI, #1
        blt     .Ldgemm_kernel_L1_END

        cmp     counterI, #4
        blt     .Ldgemm_kernel_L1_M2_BEGIN

.Ldgemm_kernel_L1_M4_20:

        INIT4x1

        mov     pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp     counterL , #0
        ble     .Ldgemm_kernel_L1_M4_40
        .align 4

.Ldgemm_kernel_L1_M4_22:

        KERNEL4x1_SUB
        KERNEL4x1_SUB
        KERNEL4x1_SUB
        KERNEL4x1_SUB

        KERNEL4x1_SUB
        KERNEL4x1_SUB
        KERNEL4x1_SUB
        KERNEL4x1_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L1_M4_22

.Ldgemm_kernel_L1_M4_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble     .Ldgemm_kernel_L1_M4_100

.Ldgemm_kernel_L1_M4_42:

        KERNEL4x1_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L1_M4_42

.Ldgemm_kernel_L1_M4_100:

        SAVE4x1

.Ldgemm_kernel_L1_M4_END:

        subs    counterI, counterI, #4

//------------------------------------------------------------------------------

.Ldgemm_kernel_L1_M2_BEGIN:

        cmp     counterI, #1
        blt     .Ldgemm_kernel_L1_END

        cmp     counterI, #2            // counterI = counterI / 2
        blt     .Ldgemm_kernel_L1_M1_BEGIN

.Ldgemm_kernel_L1_M2_20:

        INIT2x1

        mov     pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp     counterL , #0
        ble     .Ldgemm_kernel_L1_M2_40

.Ldgemm_kernel_L1_M2_22:

        KERNEL2x1_SUB
        KERNEL2x1_SUB
        KERNEL2x1_SUB
        KERNEL2x1_SUB

        KERNEL2x1_SUB
        KERNEL2x1_SUB
        KERNEL2x1_SUB
        KERNEL2x1_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L1_M2_22

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

        subs    counterI, counterI, #2

.Ldgemm_kernel_L1_M1_BEGIN:

        cmp     counterI, #1            // counterI = counterI % 2
        blt     .Ldgemm_kernel_L1_END

.Ldgemm_kernel_L1_M1_20:

        INIT1x1

        mov     pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp     counterL , #0
        ble     .Ldgemm_kernel_L1_M1_40

.Ldgemm_kernel_L1_M1_22:

        KERNEL1x1_SUB
        KERNEL1x1_SUB
        KERNEL1x1_SUB
        KERNEL1x1_SUB

        KERNEL1x1_SUB
        KERNEL1x1_SUB
        KERNEL1x1_SUB
        KERNEL1x1_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L1_M1_22

.Ldgemm_kernel_L1_M1_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble     .Ldgemm_kernel_L1_M1_100

.Ldgemm_kernel_L1_M1_42:

        KERNEL1x1_SUB

        subs    counterL, counterL, #1
        bgt     .Ldgemm_kernel_L1_M1_42

.Ldgemm_kernel_L1_M1_100:

        SAVE1x1

.Ldgemm_kernel_L1_END:

.Ldgemm_kernel_L999:

        mov     x0, #0                // set return value
        ldp     d8, d9, [sp, #(0 * 16)]
        ldp     d10, d11, [sp, #(1 * 16)]
        ldp     d12, d13, [sp, #(2 * 16)]
        ldp     d14, d15, [sp, #(3 * 16)]
        ldp     d16, d17, [sp, #(4 * 16)]
        ldp     x18, x19, [sp, #(5 * 16)]
        ldp     x20, x21, [sp, #(6 * 16)]
        ldp     x22, x23, [sp, #(7 * 16)]
        ldp     x24, x25, [sp, #(8 * 16)]
        ldp     x26, x27, [sp, #(9 * 16)]
        ldr     x28, [sp, #(10 * 16)]
        add     sp, sp, #(11*16)
        ret

        EPILOGUE

