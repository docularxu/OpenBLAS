/*******************************************************************************
Copyright (c) 2015, The OpenBLAS Project
Copyright (c) 2020, Hisilicon Limited
All rights reserved.

Based on the original sgemm_kernel_12x8.S and modified by SVE instructions
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

/*                   X0          X1          X2          s0        X3        x4       x5           x6  */
/*int CNAME(BLASLONG bm,BLASLONG bn,BLASLONG bk,FLOAT alpha,FLOAT* ba,FLOAT* bb,FLOAT* C,BLASLONG ldc) */

#define origM         x0
#define origN         x1
#define origK         x2
#define origPA        x3
#define origPB        x4
#define pC            x5
#define LDC           x6
#define temp          x7
#define counterL      x8
#define counterI      x9
#define counterJ      x10
#define pB            x11
#define pCRow0        x12
#define pCRow1        x13
#define pCRow2        x14
#define pCRow3        x15
#define pA            x16
#define alpha         w17
#define pCRow4        x18
#define pCRow5        x19
#define pCRow6        x20
#define pCRow7        x21
#define pCRow8        x22
#define pCRow9        x23
#define pCRow10       x24
#define pCRow11       x25
#define temp1         x26
#define temp2         x27

#define alpha0        z3.s
#define alphaV0       z3.s[0]

#define A_PRE_SIZE    896
#define B_PRE_SIZE    416
#define C_PRE_SIZE    160

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
.macro INIT16x12
        fmov    s8,  wzr
        fmov    s9,  wzr
        fmov    s10, s8
        fmov    s11, s9
        fmov    s12, s8
        fmov    s13, s9
        fmov    s14, s8
        fmov    s15, s9
        fmov    s16, wzr
        fmov    s17, wzr
        fmov    s18, s16
        fmov    s19, s17
        fmov    s20, wzr
        fmov    s21, s16
        fmov    s22, s17
        fmov    s23, s18
        fmov    s24, wzr
        fmov    s25, s16
        fmov    s26, s17
        fmov    s27, s18
        fmov    s28, wzr
        fmov    s29, s16
        fmov    s30, s17
        fmov    s31, s18
.endm

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

.macro KERNEL16x12_SUB
        ld1w    {z0.s}, p0/z, [pA]
        ld1w    {z1.s}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        ld1rqw  {z4.s}, p0/z, [pB]

        fmla    z8.s,  z0.s,  z4.s[0]
        fmla    z9.s,  z1.s,  z4.s[0]
        fmla    z10.s,  z0.s,  z4.s[1]
        fmla    z11.s,  z1.s,  z4.s[1]

        ld1rqw  {z5.s}, p0/z, [pB, #16]

        fmla    z12.s,  z0.s,  z4.s[2]
        fmla    z13.s,  z1.s,  z4.s[2]
        fmla    z14.s,  z0.s,  z4.s[3]
        fmla    z15.s,  z1.s,  z4.s[3]

        ld1rqw  {z6.s}, p0/z, [pB, #32]
        add     pB, pB, #48

        fmla    z16.s,  z0.s,  z5.s[0]
        fmla    z17.s,  z1.s,  z5.s[0]
        fmla    z18.s,  z0.s,  z5.s[1]
        fmla    z19.s,  z1.s,  z5.s[1]

        fmla    z20.s,  z0.s,  z5.s[2]
        fmla    z21.s,  z1.s,  z5.s[2]
        fmla    z22.s,  z0.s,  z5.s[3]
        fmla    z23.s,  z1.s,  z5.s[3]

        fmla    z24.s,  z0.s,  z6.s[0]
        fmla    z25.s,  z1.s,  z6.s[0]
        fmla    z26.s,  z0.s,  z6.s[1]
        fmla    z27.s,  z1.s,  z6.s[1]

        fmla    z28.s,  z0.s,  z6.s[2]
        fmla    z29.s,  z1.s,  z6.s[2]
        fmla    z30.s,  z0.s,  z6.s[3]
        fmla    z31.s,  z1.s,  z6.s[3]
.endm

.macro SAVE16x12
        dup     alpha0, alpha

        ld1w    {z0.s}, p0/z, [pCRow0]
        ld1w    {z1.s}, p0/z, [pCRow0, #1, MUL VL]
        fmla    z0.s, z8.s, alphaV0
        fmla    z1.s, z9.s, alphaV0
        st1w    {z0.s}, p0, [pCRow0]
        st1w    {z1.s}, p0, [pCRow0, #1, MUL VL]
        incb    pCRow0, all, MUL #2

        ld1w    {z4.s}, p0/z, [pCRow1]
        ld1w    {z5.s}, p0/z, [pCRow1, #1, MUL VL]
        fmla    z4.s, z10.s, alphaV0
        fmla    z5.s, z11.s, alphaV0
        st1w    {z4.s}, p0, [pCRow1]
        st1w    {z5.s}, p0, [pCRow1, #1, MUL VL]
        incb    pCRow1, all, MUL #2

        ld1w    {z6.s}, p0/z, [pCRow2]
        ld1w    {z7.s}, p0/z, [pCRow2, #1, MUL VL]
        fmla    z6.s, z12.s, alphaV0
        fmla    z7.s, z13.s, alphaV0
        st1w    {z6.s}, p0, [pCRow2]
        st1w    {z7.s}, p0, [pCRow2, #1, MUL VL]
        incb    pCRow2, all, MUL #2

        ld1w    {z0.s}, p0/z, [pCRow3]
        ld1w    {z1.s}, p0/z, [pCRow3, #1, MUL VL]
        fmla    z0.s, z14.s, alphaV0
        fmla    z1.s, z15.s, alphaV0
        st1w    {z0.s}, p0, [pCRow3]
        st1w    {z1.s}, p0, [pCRow3, #1, MUL VL]
        incb    pCRow3, all, MUL #2

        ld1w    {z4.s}, p0/z, [pCRow4]
        ld1w    {z5.s}, p0/z, [pCRow4, #1, MUL VL]
        fmla    z4.s, z16.s, alphaV0
        fmla    z5.s, z17.s, alphaV0
        st1w    {z4.s}, p0, [pCRow4]
        st1w    {z5.s}, p0, [pCRow4, #1, MUL VL]
        incb    pCRow4, all, MUL #2

        ld1w    {z6.s}, p0/z, [pCRow5]
        ld1w    {z7.s}, p0/z, [pCRow5, #1, MUL VL]
        fmla    z6.s, z18.s, alphaV0
        fmla    z7.s, z19.s, alphaV0
        st1w    {z6.s}, p0, [pCRow5]
        st1w    {z7.s}, p0, [pCRow5, #1, MUL VL]
        incb    pCRow5, all, MUL #2

        ld1w    {z0.s}, p0/z, [pCRow6]
        ld1w    {z1.s}, p0/z, [pCRow6, #1, MUL VL]
        fmla    z0.s, z20.s, alphaV0
        fmla    z1.s, z21.s, alphaV0
        st1w    {z0.s}, p0, [pCRow6]
        st1w    {z1.s}, p0, [pCRow6, #1, MUL VL]
        incb    pCRow6, all, MUL #2

        ld1w    {z4.s}, p0/z, [pCRow7]
        ld1w    {z5.s}, p0/z, [pCRow7, #1, MUL VL]
        fmla    z4.s, z22.s, alphaV0
        fmla    z5.s, z23.s, alphaV0
        st1w    {z4.s}, p0, [pCRow7]
        st1w    {z5.s}, p0, [pCRow7, #1, MUL VL]
        incb    pCRow7, all, MUL #2

        ld1w    {z6.s}, p0/z, [pCRow8]
        ld1w    {z7.s}, p0/z, [pCRow8, #1, MUL VL]
        fmla    z6.s, z24.s, alphaV0
        fmla    z7.s, z25.s, alphaV0
        st1w    {z6.s}, p0, [pCRow8]
        st1w    {z7.s}, p0, [pCRow8, #1, MUL VL]
        incb    pCRow8, all, MUL #2

        ld1w    {z0.s}, p0/z, [pCRow9]
        ld1w    {z1.s}, p0/z, [pCRow9, #1, MUL VL]
        fmla    z0.s, z26.s, alphaV0
        fmla    z1.s, z27.s, alphaV0
        st1w    {z0.s}, p0, [pCRow9]
        st1w    {z1.s}, p0, [pCRow9, #1, MUL VL]
        incb    pCRow9, all, MUL #2

        ld1w    {z4.s}, p0/z, [pCRow10]
        ld1w    {z5.s}, p0/z, [pCRow10, #1, MUL VL]
        fmla    z4.s, z28.s, alphaV0
        fmla    z5.s, z29.s, alphaV0
        st1w    {z4.s}, p0, [pCRow10]
        st1w    {z5.s}, p0, [pCRow10, #1, MUL VL]
        incb    pCRow10, all, MUL #2

        ld1w    {z6.s}, p0/z, [pCRow11]
        ld1w    {z7.s}, p0/z, [pCRow11, #1, MUL VL]
        fmla    z6.s, z30.s, alphaV0
        fmla    z7.s, z31.s, alphaV0
        st1w    {z6.s}, p0, [pCRow11]
        st1w    {z7.s}, p0, [pCRow11, #1, MUL VL]
        incb    pCRow11, all, MUL #2
.endm

/******************************************************************************/
.macro INIT16x8
        fmov    s8,  wzr
        fmov    s9,  wzr
        fmov    s10, s9
        fmov    s11, s8
        fmov    s12, s8
        fmov    s13, s9
        fmov    s14, wzr
        fmov    s15, s11
        fmov    s16, wzr
        fmov    s17, s11
        fmov    s18, s12
        fmov    s19, wzr
        fmov    s20, s12
        fmov    s21, s14
        fmov    s22, s15
        fmov    s23, s17
.endm

.macro KERNEL16x8_I
        ld1w    {z0.s}, p0/z, [pA]
        ld1rqw  {z4.s}, p0/z, [pB]
        ld1rqw  {z5.s}, p0/z, [pB, #16]

        fmul    z8.s, z0.s, z4.s[0]
        fmul    z10.s, z0.s, z4.s[1]
        fmul    z12.s, z0.s, z4.s[2]
        fmul    z14.s, z0.s, z4.s[3]

        ld1w    {z1.s}, p0/z, [pA, #1, MUL VL]

        fmul    z16.s, z0.s, z5.s[0]
        fmul    z18.s, z0.s, z5.s[1]
        fmul    z20.s, z0.s, z5.s[2]
        fmul    z22.s, z0.s, z5.s[3]

        ld1w    {z6.s}, p0/z, [pA, #2, MUL VL]

        fmul    z9.s, z1.s, z4.s[0]
        fmul    z11.s, z1.s, z4.s[1]
        fmul    z13.s, z1.s, z4.s[2]
        fmul    z15.s, z1.s, z4.s[3]

        ld1w    {z7.s}, p0/z, [pA, #3, MUL VL]
        incb    pA, all, MUL #4

        fmul    z17.s, z1.s, z5.s[0]
        fmul    z19.s, z1.s, z5.s[1]
        fmul    z21.s, z1.s, z5.s[2]
        fmul    z23.s, z1.s, z5.s[3]

        ld1rqw  {z2.s}, p0/z, [pB,#32]
        ld1rqw  {z3.s}, p0/z, [pB, #48]
        add     pB, pB, #64
.endm

.macro KERNEL16x8_M1
        fmla    z8.s, z0.s,  z4.s[0]
        fmla    z10.s, z0.s, z4.s[1]
        fmla    z12.s, z0.s, z4.s[2]
        fmla    z14.s, z0.s, z4.s[3]

        ld1w    {z6.s}, p0/z, [pA]

        fmla    z16.s, z0.s, z5.s[0]
        fmla    z18.s, z0.s, z5.s[1]
        fmla    z20.s, z0.s, z5.s[2]
        fmla    z22.s, z0.s, z5.s[3]

        ld1rqw  {z2.s}, p0/z, [pB]

        fmla    z9.s, z1.s, z4.s[0]
        fmla    z11.s, z1.s, z4.s[1]
        fmla    z13.s, z1.s, z4.s[2]
        fmla    z15.s, z1.s, z4.s[3]

        ld1rqw  {z3.s}, p0/z, [pB, #16]
        add     pB, pB, #32

        fmla    z17.s, z1.s, z5.s[0]
        fmla    z19.s, z1.s, z5.s[1]
        fmla    z21.s, z1.s, z5.s[2]
        fmla    z23.s, z1.s, z5.s[3]

        ld1w    {z7.s}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2
.endm

.macro KERNEL16x8_M2
        fmla    z8.s, z6.s, z2.s[0]
        fmla    z10.s, z6.s, z2.s[1]
        fmla    z12.s, z6.s, z2.s[2]
        fmla    z14.s, z6.s, z2.s[3]

        ld1w    {z0.s}, p0/z, [pA]

        fmla    z16.s, z6.s, z3.s[0]
        fmla    z18.s, z6.s, z3.s[1]
        fmla    z20.s, z6.s, z3.s[2]
        fmla    z22.s, z6.s, z3.s[3]

        ld1rqw  {z4.s}, p0/z, [pB]

        fmla    z9.s, z7.s, z2.s[0]
        fmla    z11.s, z7.s, z2.s[1]
        fmla    z13.s, z7.s, z2.s[2]
        fmla    z15.s, z7.s, z2.s[3]

        ld1rqw  {z5.s}, p0/z, [pB, #16]
        add     pB, pB, #32

        fmla    z17.s, z7.s, z3.s[0]
        fmla    z19.s, z7.s, z3.s[1]
        fmla    z21.s, z7.s, z3.s[2]
        fmla    z23.s, z7.s, z3.s[3]

        ld1w    {z1.s}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2
.endm

.macro KERNEL16x8_E
        fmla    z8.s, z6.s, z2.s[0]
        fmla    z10.s, z6.s, z2.s[1]
        fmla    z12.s, z6.s, z2.s[2]
        fmla    z14.s, z6.s, z2.s[3]
        fmla    z16.s, z6.s, z3.s[0]
        fmla    z18.s, z6.s, z3.s[1]
        fmla    z20.s, z6.s, z3.s[2]
        fmla    z22.s, z6.s, z3.s[3]
        fmla    z9.s, z7.s, z2.s[0]
        fmla    z11.s, z7.s, z2.s[1]
        fmla    z13.s, z7.s, z2.s[2]
        fmla    z15.s, z7.s, z2.s[3]
        fmla    z17.s, z7.s, z3.s[0]
        fmla    z19.s, z7.s, z3.s[1]
        fmla    z21.s, z7.s, z3.s[2]
        fmla    z23.s, z7.s, z3.s[3]
.endm


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
        dup     alpha0,   alpha
        
        ld1w    {z0.s}, p0/z, [pCRow0]
        ld1w    {z1.s}, p0/z, [pCRow0, #1, MUL VL]
        fmla    z0.s, z8.s, alphaV0
        fmla    z1.s, z9.s, alphaV0
        st1w    {z0.s}, p0, [pCRow0]
        st1w    {z1.s}, p0, [pCRow0, #1, MUL VL]
        incb    pCRow0, all, MUL #2

        ld1w    {z4.s}, p0/z, [pCRow1]
        ld1w    {z5.s}, p0/z, [pCRow1, #1, MUL VL]
        fmla    z4.s, z10.s, alphaV0
        fmla    z5.s, z11.s, alphaV0
        st1w    {z4.s}, p0, [pCRow1]
        st1w    {z5.s}, p0, [pCRow1, #1, MUL VL]
        incb    pCRow1, all, MUL #2

        ld1w    {z0.s}, p0/z, [pCRow2]
        ld1w    {z1.s}, p0/z, [pCRow2, #1, MUL VL]
        fmla    z0.s, z12.s, alphaV0
        fmla    z1.s, z13.s, alphaV0
        st1w    {z0.s}, p0, [pCRow2]
        st1w    {z1.s}, p0, [pCRow2, #1, MUL VL]
        incb    pCRow2, all, MUL #2

        ld1w    {z4.s}, p0/z, [pCRow3]
        ld1w    {z5.s}, p0/z, [pCRow3, #1, MUL VL]
        fmla    z4.s, z14.s, alphaV0
        fmla    z5.s, z15.s, alphaV0
        st1w    {z4.s}, p0, [pCRow3]
        st1w    {z5.s}, p0, [pCRow3, #1, MUL VL]
        incb    pCRow3, all, MUL #2

        ld1w    {z0.s}, p0/z, [pCRow4]
        ld1w    {z1.s}, p0/z, [pCRow4, #1, MUL VL]
        fmla    z0.s, z16.s, alphaV0
        fmla    z1.s, z17.s, alphaV0
        st1w    {z0.s}, p0, [pCRow4]
        st1w    {z1.s}, p0, [pCRow4, #1, MUL VL]
        incb    pCRow4, all, MUL #2

        ld1w    {z4.s}, p0/z, [pCRow5]
        ld1w    {z5.s}, p0/z, [pCRow5, #1, MUL VL]
        fmla    z4.s, z18.s, alphaV0
        fmla    z5.s, z19.s, alphaV0
        st1w    {z4.s}, p0, [pCRow5]
        st1w    {z5.s}, p0, [pCRow5, #1, MUL VL]
        incb    pCRow5, all, MUL #2

        ld1w    {z0.s}, p0/z, [pCRow6]
        ld1w    {z1.s}, p0/z, [pCRow6, #1, MUL VL]
        fmla    z0.s, z20.s, alphaV0
        fmla    z1.s, z21.s, alphaV0
        st1w    {z0.s}, p0, [pCRow6]
        st1w    {z1.s}, p0, [pCRow6, #1, MUL VL]
        incb    pCRow6, all, MUL #2

        ld1w    {z4.s}, p0/z, [pCRow7]
        ld1w    {z5.s}, p0/z, [pCRow7, #1, MUL VL]
        fmla    z4.s, z22.s, alphaV0
        fmla    z5.s, z23.s, alphaV0
        st1w    {z4.s}, p0, [pCRow7]
        st1w    {z5.s}, p0, [pCRow7, #1, MUL VL]
        incb    pCRow7, all, MUL #2
.endm

/******************************************************************************/
.macro INIT16x4
        fmov    s8,  wzr
        fmov    s9,  wzr
        fmov    s10, s9
        fmov    s11, s8
        fmov    s12, s8
        fmov    s13, s9
        fmov    s14, wzr
        fmov    s15, s11
.endm

.macro KERNEL16x4_I
        ld1w    {z0.s}, p0/z, [pA]
        ld1w    {z1.s}, p0/z, [pA, #1, MUL VL]

        ld1rqw  {z4.s}, p0/z, [pB]

        fmul    z8.s, z0.s, z4.s[0]
        fmul    z10.s, z0.s, z4.s[1]

        ld1w    {z6.s}, p0/z, [pA, #2, MUL VL]

        fmul    z12.s, z0.s, z4.s[2]
        fmul    z14.s, z0.s, z4.s[3]

        ld1w    {z7.s}, p0/z, [pA, #3, MUL VL]
        incb    pA, all, MUL #4

        fmul    z9.s, z1.s, z4.s[0]
        fmul    z11.s, z1.s, z4.s[1]
        fmul    z13.s, z1.s, z4.s[2]
        fmul    z15.s, z1.s, z4.s[3]

        ld1rqw  {z2.s}, p0/z, [pB, #16]
        add     pB, pB, #32
.endm

.macro KERNEL16x4_M1
        fmla    z8.s, z0.s,  z4.s[0]
        fmla    z10.s, z0.s, z4.s[1]

        ld1w    {z6.s}, p0/z, [pA]

        fmla    z12.s, z0.s, z4.s[2]
        fmla    z14.s, z0.s, z4.s[3]

        ld1rqw  {z2.s}, p0/z, [pB]
        add     pB, pB, #16

        fmla    z9.s, z1.s, z4.s[0]
        fmla    z11.s, z1.s, z4.s[1]

        ld1w    {z7.s}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        fmla    z13.s, z1.s, z4.s[2]
        fmla    z15.s, z1.s, z4.s[3]
.endm

.macro KERNEL16x4_M2
        fmla    z8.s, z6.s, z2.s[0]
        fmla    z10.s, z6.s, z2.s[1]

        ld1w    {z0.s}, p0/z, [pA]

        fmla    z12.s, z6.s, z2.s[2]
        fmla    z14.s, z6.s, z2.s[3]

        ld1rqw  {z4.s}, p0/z, [pB]
        add     pB, pB, #16

        fmla    z9.s, z7.s, z2.s[0]
        fmla    z11.s, z7.s, z2.s[1]

        ld1w    {z1.s}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        fmla    z13.s, z7.s, z2.s[2]
        fmla    z15.s, z7.s, z2.s[3]
.endm

.macro KERNEL16x4_E
        fmla    z8.s, z6.s, z2.s[0]
        fmla    z10.s, z6.s, z2.s[1]
        fmla    z12.s, z6.s, z2.s[2]
        fmla    z14.s, z6.s, z2.s[3]
        fmla    z9.s, z7.s, z2.s[0]
        fmla    z11.s, z7.s, z2.s[1]
        fmla    z13.s, z7.s, z2.s[2]
        fmla    z15.s, z7.s, z2.s[3]
.endm


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
        dup     alpha0, alpha

        ld1w    {z0.s}, p0/z, [pCRow0]
        ld1w    {z1.s}, p0/z, [pCRow0, #1, MUL VL]
        fmla    z0.s, z8.s, alphaV0
        fmla    z1.s, z9.s, alphaV0
        st1w    {z0.s}, p0, [pCRow0]
        st1w    {z1.s}, p0, [pCRow0, #1, MUL VL]
        incb    pCRow0, all, MUL #2

        ld1w    {z4.s}, p0/z, [pCRow1]
        ld1w    {z5.s}, p0/z, [pCRow1, #1, MUL VL]
        fmla    z4.s, z10.s, alphaV0
        fmla    z5.s, z11.s, alphaV0
        st1w    {z4.s}, p0, [pCRow1]
        st1w    {z5.s}, p0, [pCRow1, #1, MUL VL]
        incb    pCRow1, all, MUL #2

        ld1w    {z0.s}, p0/z, [pCRow2]
        ld1w    {z1.s}, p0/z, [pCRow2, #1, MUL VL]
        fmla    z0.s, z12.s, alphaV0
        fmla    z1.s, z13.s, alphaV0
        st1w    {z0.s}, p0, [pCRow2]
        st1w    {z1.s}, p0, [pCRow2, #1, MUL VL]
        incb    pCRow2, all, MUL #2

        ld1w    {z4.s}, p0/z, [pCRow3]
        ld1w    {z5.s}, p0/z, [pCRow3, #1, MUL VL]
        fmla    z4.s, z14.s, alphaV0
        fmla    z5.s, z15.s, alphaV0
        st1w    {z4.s}, p0, [pCRow3]
        st1w    {z5.s}, p0, [pCRow3, #1, MUL VL]
        incb    pCRow3, all, MUL #2
.endm

/******************************************************************************/
.macro INIT16x2
        fmov    s8,  wzr
        fmov    s9,  wzr
        fmov    s10, s9
        fmov    s11, s8
.endm

.macro KERNEL16x2_I
        ld1w    {z0.s}, p0/z, [pA]
        ld1w    {z1.s}, p0/z, [pA, #1, MUL VL]
        ld1rqw  {z4.s}, p0/z, [pB]

        fmul    z8.s, z0.s, z4.s[0]
        fmul    z10.s, z0.s, z4.s[1]
        ld1w    {z6.s}, p0/z, [pA, #2, MUL VL]

        fmul    z9.s, z1.s, z4.s[0]
        fmul    z11.s, z1.s, z4.s[1]

        ld1w    {z7.s}, p0/z, [pA, #3, MUL VL]
        incb    pA, all, MUL #4
        ld1rqw  {z2.s}, p0/z, [pB, #8]
        add     pB, pB, #16
.endm

.macro KERNEL16x2_M1
        ld1w    {z6.s}, p0/z, [pA]
        fmla    z8.s, z0.s,  z4.s[0]
        fmla    z10.s, z0.s, z4.s[1]

        ld1rqw  {z2.s}, p0/z, [pB]
        add     pB, pB, #8 

        fmla    z9.s, z1.s, z4.s[0]
        fmla    z11.s, z1.s, z4.s[1]

        ld1w    {z7.s}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2
.endm

.macro KERNEL16x2_M2
        ld1w    {z0.s}, p0/z, [pA]

        fmla    z8.s, z6.s, z2.s[0]
        fmla    z10.s, z6.s, z2.s[1]

        ld1rqw  {z4.s}, p0/z, [pB]
        add     pB, pB, #8

        fmla    z9.s, z7.s, z2.s[0]
        fmla    z11.s, z7.s, z2.s[1]

        ld1w    {z1.s}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2
.endm

.macro KERNEL16x2_E
        fmla    z8.s, z6.s, z2.s[0]
        fmla    z10.s, z6.s, z2.s[1]
        fmla    z9.s, z7.s, z2.s[0]
        fmla    z11.s, z7.s, z2.s[1]
.endm

.macro KERNEL16x2_SUB
        ld1w    {z0.s}, p0/z, [pA]
        ld1w    {z1.s}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        ld1rqw  {z4.s}, p0/z, [pB]
        add     pB, pB, #8

        fmla    z8.s, z0.s, z4.s[0]
        fmla    z10.s, z0.s, z4.s[1]
        fmla    z9.s, z1.s, z4.s[0]
        fmla    z11.s, z1.s, z4.s[1]
.endm

.macro SAVE16x2
        dup     alpha0, alpha

        ld1w    {z0.s}, p0/z, [pCRow0]
        ld1w    {z1.s}, p0/z, [pCRow0, #1, MUL VL]
        fmla    z0.s, z8.s, alphaV0
        fmla    z1.s, z9.s, alphaV0
        st1w    {z0.s}, p0, [pCRow0]
        st1w    {z1.s}, p0, [pCRow0, #1, MUL VL]
        incb    pCRow0, all, MUL #2

        ld1w    {z4.s}, p0/z, [pCRow1]
        ld1w    {z5.s}, p0/z, [pCRow1, #1, MUL VL]
        fmla    z4.s, z10.s, alphaV0
        fmla    z5.s, z11.s, alphaV0
        st1w    {z4.s}, p0, [pCRow1]
        st1w    {z5.s}, p0, [pCRow1, #1, MUL VL]
        incb    pCRow1, all, MUL #2
.endm

/******************************************************************************/
.macro INIT16x1
        fmov    s8,  wzr
        fmov    s9,  wzr
.endm

.macro KERNEL16x1_SUB
        ld1w    {z0.s}, p0/z, [pA]
        ld1w    {z1.s}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        ld1rw  {z4.s}, p0/z, [pB]
        add     pB, pB, #4

        fmla    z8.s, z0.s, z4.s[0]
        fmla    z9.s, z1.s, z4.s[0]
.endm

.macro SAVE16x1
        dup     alpha0, alpha

        ld1w    {z0.s}, p0/z, [pCRow0]
        ld1w    {z1.s}, p0/z, [pCRow0, #1, MUL VL]
        fmla    z0.s, z8.s, alphaV0
        fmla    z1.s, z9.s, alphaV0
        st1w    {z0.s}, p0, [pCRow0]
        st1w    {z1.s}, p0, [pCRow0, #1, MUL VL]
        incb    pCRow0, all, MUL #2
.endm


/******************************************************************************/
.macro INIT8x12
        fmov    s8, wzr
        fmov    s9, wzr
        fmov    s10, s8
        fmov    s11, s9
        fmov    s12, s8
        fmov    s13, s9
        fmov    s14, s8
        fmov    s15, s9
        fmov    s16, wzr
        fmov    s17, wzr
        fmov    s18, s16
        fmov    s19, s17
.endm

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
        dup     alpha0,   alpha
        
        ld1w    {z0.s}, p0/z, [pCRow0]
        fmla    z0.s, z8.s, alphaV0
        st1w    {z0.s}, p0, [pCRow0]
        add     pCRow0, pCRow0, #32

        ld1w    {z1.s}, p0/z, [pCRow1]
        fmla    z1.s, z9.s, alphaV0
        st1w    {z1.s}, p0, [pCRow1]
        add     pCRow1, pCRow1, #32

        ld1w    {z2.s}, p0/z, [pCRow2]
        fmla    z2.s, z10.s, alphaV0
        st1w    {z2.s}, p0, [pCRow2]
        add     pCRow2, pCRow2, #32

        ld1w    {z4.s}, p0/z, [pCRow3]
        fmla    z4.s, z11.s, alphaV0
        st1w    {z4.s}, p0, [pCRow3]
        add     pCRow3, pCRow3, #32

        ld1w    {z5.s}, p0/z, [pCRow4]
        fmla    z5.s, z12.s, alphaV0
        st1w    {z5.s}, p0, [pCRow4]
        add     pCRow4, pCRow4, #32

        ld1w    {z6.s}, p0/z, [pCRow5]
        fmla    z6.s, z13.s, alphaV0
        st1w    {z6.s}, p0, [pCRow5]
        add     pCRow5, pCRow5, #32

        ld1w    {z7.s}, p0/z, [pCRow6]
        fmla    z7.s, z14.s, alphaV0
        st1w    {z7.s}, p0, [pCRow6]
        add     pCRow6, pCRow6, #32

        ld1w    {z0.s}, p0/z, [pCRow7]
        fmla    z0.s, z15.s, alphaV0
        st1w    {z0.s}, p0, [pCRow7]
        add     pCRow7, pCRow7, #32
        
        ld1w    {z1.s}, p0/z, [pCRow8]
        fmla    z1.s, z16.s, alphaV0
        st1w    {z1.s}, p0, [pCRow8]
        add     pCRow8, pCRow8, #32

        ld1w    {z2.s}, p0/z, [pCRow9]
        fmla    z2.s, z17.s, alphaV0
        st1w    {z2.s}, p0, [pCRow9]
        add     pCRow9, pCRow9, #32

        ld1w    {z4.s}, p0/z, [pCRow10]
        fmla    z4.s, z18.s, alphaV0
        st1w    {z4.s}, p0, [pCRow10]
        add     pCRow10, pCRow10, #32

        ld1w    {z5.s}, p0/z, [pCRow11]
        fmla    z5.s, z19.s, alphaV0
        st1w    {z5.s}, p0, [pCRow11]
        add     pCRow11, pCRow11, #32
.endm

/******************************************************************************/
.macro INIT8x8
        fmov    s8,  wzr
        fmov    s9, s8
        fmov    s10, s8
        fmov    s11, wzr
        fmov    s12, wzr
        fmov    s13, s9
        fmov    s14, s10
        fmov    s15, s11
.endm

.macro KERNEL8x8_I
        ld1w    {z0.s}, p0/z, [pA]

        ld1rqw  {z2.s}, p0/z, [pB]
        ld1rqw  {z3.s}, p0/z, [pB, #16]

        fmul    z8.s, z0.s, z2.s[0]
        fmul    z9.s, z0.s, z2.s[1]

        ld1w    {z1.s}, p0/z, [pA, #1, MUL VL]
        incb    pA, all, MUL #2

        fmul    z10.s, z0.s, z2.s[2]
        fmul    z11.s, z0.s, z2.s[3]

        ld1rqw  {z5.s}, p0/z, [pB, #32]

        fmul    z12.s, z0.s, z3.s[0]
        fmul    z13.s, z0.s, z3.s[1]

        ld1rqw  {z6.s}, p0/z, [pB, #48]
        add     pB, pB, #64

        fmul    z14.s, z0.s, z3.s[2]
        fmul    z15.s, z0.s, z3.s[3]
.endm


.macro KERNEL8x8_M1
        fmla    z8.s, z0.s, z2.s[0]
        fmla    z9.s, z0.s, z2.s[1]

        ld1w    {z1.s}, p0/z, [pA]
        add     pA, pA, #32

        fmla    z10.s, z0.s, z2.s[2]
        fmla    z11.s, z0.s, z2.s[3]

        ld1rqw  {z5.s}, p0/z, [pB]

        fmla    z12.s, z0.s, z3.s[0]
        fmla    z13.s, z0.s, z3.s[1]

        ld1rqw  {z6.s}, p0/z, [pB, #16]
        add     pB, pB, #32

        fmla    z14.s, z0.s, z3.s[2]
        fmla    z15.s, z0.s, z3.s[3]
.endm

.macro KERNEL8x8_M2
        fmla    z8.s, z1.s, z5.s[0]
        fmla    z9.s, z1.s, z5.s[1]

        ld1w    {z0.s}, p0/z, [pA]
        add     pA, pA, #32

        fmla    z10.s, z1.s, z5.s[2]
        fmla    z11.s, z1.s, z5.s[3]

        ld1rqw  {z2.s}, p0/z, [pB]

        fmla    z12.s, z1.s, z6.s[0]
        fmla    z13.s, z1.s, z6.s[1]

        ld1rqw  {z3.s}, p0/z, [pB, #16]
        add     pB, pB, #32

        fmla    z14.s, z1.s, z6.s[2]
        fmla    z15.s, z1.s, z6.s[3]
.endm


.macro KERNEL8x8_E
        fmla    z8.s, z1.s, z5.s[0]
        fmla    z9.s, z1.s, z5.s[1]
        fmla    z10.s, z1.s, z5.s[2]
        fmla    z11.s, z1.s, z5.s[3]
        fmla    z12.s, z1.s, z6.s[0]
        fmla    z13.s, z1.s, z6.s[1]
        fmla    z14.s, z1.s, z6.s[2]
        fmla    z15.s, z1.s, z6.s[3]
.endm


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
        dup     alpha0,   alpha

        ld1w    {z0.s}, p0/z, [pCRow0]
        fmla    z0.s, z8.s, alphaV0
        st1w    {z0.s}, p0, [pCRow0]
        add     pCRow0, pCRow0, #32

        ld1w    {z1.s}, p0/z, [pCRow1]
        fmla    z1.s, z9.s, alphaV0
        st1w    {z1.s}, p0, [pCRow1]
        add     pCRow1, pCRow1, #32

        ld1w    {z2.s}, p0/z, [pCRow2]
        fmla    z2.s, z10.s, alphaV0
        st1w    {z2.s}, p0, [pCRow2]
        add     pCRow2, pCRow2, #32
        
        ld1w    {z4.s}, p0/z, [pCRow3]
        fmla    z4.s, z11.s, alphaV0
        st1w    {z4.s}, p0, [pCRow3]
        add     pCRow3, pCRow3, #32

        ld1w    {z0.s}, p0/z, [pCRow4]
        fmla    z0.s, z12.s, alphaV0
        st1w    {z0.s}, p0, [pCRow4]
        add     pCRow4, pCRow4, #32

        ld1w    {z1.s}, p0/z, [pCRow5]
        fmla    z1.s, z13.s, alphaV0
        st1w    {z1.s}, p0, [pCRow5]
        add     pCRow5, pCRow5, #32

        ld1w    {z2.s}, p0/z, [pCRow6]
        fmla    z2.s, z14.s, alphaV0
        st1w    {z2.s}, p0, [pCRow6]
        add     pCRow6, pCRow6, #32

        ld1w    {z4.s}, p0/z, [pCRow7]
        fmla    z4.s, z15.s, alphaV0
        st1w    {z4.s}, p0, [pCRow7]
        add     pCRow7, pCRow7, #32
.endm

/******************************************************************************/
.macro INIT8x4
        fmov    s8,  wzr
        fmov    s9,  wzr
        fmov    s10, s9
        fmov    s11, s8
    
.endm

.macro KERNEL8x4_I
        ld1w    {z0.s}, p0/z, [pA]
        ld1rqw  {z4.s}, p0/z, [pB]

        fmul    z8.s, z0.s, z4.s[0]
        fmul    z9.s, z0.s, z4.s[1]
        fmul    z10.s, z0.s, z4.s[2]
        fmul    z11.s, z0.s, z4.s[3]

        ld1w    {z1.s}, p0/z, [pA, #1, MUL VL]
        ld1rqw  {z5.s}, p0/z, [pB, #16]

        incb    pA, all, MUL #2
        add     pB, pB, #32
.endm

.macro KERNEL8x4_M1
        fmla    z8.s, z0.s, z4.s[0]
        fmla    z9.s, z0.s, z4.s[1]
        fmla    z10.s, z0.s, z4.s[2]
        fmla    z11.s, z0.s, z4.s[3]

        ld1w    {z1.s}, p0/z, [pA]
        ld1rqw  {z5.s}, p0/z, [pB]

        add     pA, pA, #32
        add     pB, pB, #16
.endm

.macro KERNEL8x4_M2
        fmla    z8.s, z1.s, z5.s[0]
        fmla    z9.s, z1.s, z5.s[1]
        fmla    z10.s, z1.s, z5.s[2]
        fmla    z11.s, z1.s, z5.s[3]

        ld1w    {z0.s}, p0/z, [pA]
        ld1rqw  {z4.s}, p0/z, [pB]

        add     pA, pA, #32
        add     pB, pB, #16
.endm

.macro KERNEL8x4_E
        fmla    z8.s, z1.s, z5.s[0]
        fmla    z9.s, z1.s, z5.s[1]
        fmla    z10.s, z1.s, z5.s[2]
        fmla    z11.s, z1.s, z5.s[3]
.endm

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
        dup     alpha0,   alpha

        ld1w    {z0.s}, p0/z, [pCRow0]
        fmla    z0.s, z8.s, alphaV0
        st1w    {z0.s}, p0, [pCRow0]
        add     pCRow0, pCRow0, #32

        ld1w    {z4.s}, p0/z, [pCRow1]
        fmla    z4.s, z9.s, alphaV0
        st1w    {z4.s}, p0, [pCRow1]
        add     pCRow1, pCRow1, #32

        ld1w    {z0.s}, p0/z, [pCRow2]
        fmla    z0.s, z10.s, alphaV0
        st1w    {z0.s}, p0, [pCRow2]
        add     pCRow2, pCRow2, #32

        ld1w    {z4.s}, p0/z, [pCRow3]
        fmla    z4.s, z11.s, alphaV0
        st1w    {z4.s}, p0, [pCRow3]
        add     pCRow3, pCRow3, #32
.endm

/******************************************************************************/
.macro INIT8x2
        fmov    s8,  wzr
        fmov    s9,  wzr
.endm

.macro KERNEL8x2_SUB
        ld1w    {z0.s}, p0/z, [pA]
        ld1rqw  {z4.s}, p0/z, [pB]

        fmla    z8.s, z0.s, z4.s[0]
        fmla    z9.s, z0.s, z4.s[1]

        add     pA, pA, #32
        add     pB, pB, #8
.endm

.macro SAVE8x2
        dup     alpha0,   alpha

        ld1w    {z0.s}, p0/z, [pCRow0]
        fmla    z0.s, z8.s, alphaV0
        st1w    {z0.s}, p0, [pCRow0]
        add     pCRow0, pCRow0, #32

        ld1w    {z4.s}, p0/z, [pCRow1]
        fmla    z4.s, z9.s, alphaV0
        st1w    {z4.s}, p0, [pCRow1]
        add     pCRow1, pCRow1, #32
.endm
    
/******************************************************************************/
.macro INIT8x1
        fmov    s8,  wzr
.endm

.macro KERNEL8x1_SUB
        ld1w    {z0.s}, p0/z, [pA]
        ld1rw   {z4.s}, p0/z, [pB]

        fmla    z8.s, z0.s, z4.s[0]

        add     pA, pA, #32
        add     pB, pB, #4
.endm

.macro SAVE8x1
        dup     alpha0,   alpha

        ld1w    {z0.s}, p0/z, [pCRow0]
        fmla    z0.s, z8.s, alphaV0
        st1w    {z0.s}, p0, [pCRow0]
        add     pCRow0, pCRow0, #32
.endm

/******************************************************************************/
.macro INIT4x12
        fmov    s8,  wzr
        fmov    s9,  wzr
        fmov    s10, s8
        fmov    s11, s9
        fmov    s12, s8
        fmov    s13, s9
        fmov    s14, s8
        fmov    s15, s9
        fmov    s16, wzr
        fmov    s17, wzr
        fmov    s18, s16
        fmov    s19, s17
.endm

.macro KERNEL4x12_I
        ld1rqw  {z0.s}, p0/z, [pA]

        ld1rqw  {z2.s}, p0/z, [pB]
        ld1rqw  {z3.s}, p0/z, [pB, #16]
        ld1rqw  {z4.s}, p0/z, [pB, #32]

        fmul    z8.s,  z0.s,  z2.s[0]
        fmul    z9.s,  z0.s,  z2.s[1]

        ld1rqw  {z1.s}, p0/z, [pA, #16]

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

        add     pA, pA, #32
        add     pB, pB, #96
.endm

.macro KERNEL4x12_M1
        ld1rqw  {z1.s}, p0/z, [pA]
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

        add     pA, pA, #16
        add     pB, pB, #48
.endm

.macro KERNEL4x12_M2
        ld1rqw {z0.s}, p0/z, [pA]

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

        add     pA, pA, #16
        add     pB, pB, #48
.endm

.macro KERNEL4x12_E
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

.macro KERNEL4x12_SUB
        ld1rqw  {z0.s}, p0/z, [pA]
        ld1rqw  {z2.s}, p0/z, [pB]
        ld1rqw  {z3.s}, p0/z, [pB, #16]
        ld1rqw  {z4.s}, p0/z, [pB, #32]
       
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
    
        add     pA, pA, #16
        add     pB, pB, #48
.endm

.macro SAVE4x12
        dup     alpha0,   alpha

        ld1rqw    {z0.s}, p0/z, [pCRow0]
        ld1rqw    {z4.s}, p0/z, [pCRow1]

        fmla    z0.s, z8.s, alphaV0
        st1w    {z0.s}, p0, [pCRow0]
        add     pCRow0, pCRow0, #16

        fmla    z4.s, z9.s, alphaV0
        st1w    {z4.s}, p0, [pCRow1]
        add     pCRow1, pCRow1, #16

        ld1rqw    {z0.s}, p0/z, [pCRow2]
        ld1rqw    {z4.s}, p0/z, [pCRow3]

        fmla    z0.s, z10.s, alphaV0
        st1w    {z0.s}, p0, [pCRow2]
        add     pCRow2, pCRow2, #16

        fmla    z4.s, z11.s, alphaV0
        st1w    {z4.s}, p0, [pCRow3]
        add     pCRow3, pCRow3, #16

        ld1rqw    {z0.s}, p0/z, [pCRow4]
        ld1rqw    {z4.s}, p0/z, [pCRow5]

        fmla    z0.s, z12.s, alphaV0
        st1w    {z0.s}, p0, [pCRow4]
        add     pCRow4, pCRow4, #16

        fmla    z4.s, z13.s, alphaV0
        st1w    {z4.s}, p0, [pCRow5]
        add     pCRow5, pCRow5, #16

        ld1rqw    {z0.s}, p0/z, [pCRow6]
        ld1rqw    {z4.s}, p0/z, [pCRow7]
   
        fmla    z0.s, z14.s, alphaV0
        st1w    {z0.s}, p0, [pCRow6]
        add     pCRow6, pCRow6, #16

        fmla    z4.s, z15.s, alphaV0
        st1w    {z4.s}, p0, [pCRow7]
        add     pCRow7, pCRow7, #16

        ld1rqw    {z0.s}, p0/z, [pCRow8]
        ld1rqw    {z4.s}, p0/z, [pCRow9]

        fmla    z0.s, z16.s, alphaV0
        st1w    {z0.s}, p0, [pCRow8]
        add     pCRow8, pCRow8, #16

        fmla    z4.s, z17.s, alphaV0
        st1w    {z4.s}, p0, [pCRow9]
        add     pCRow9, pCRow9, #16

        ld1rqw    {z0.s}, p0/z, [pCRow10]
        ld1rqw    {z4.s}, p0/z, [pCRow11]

        fmla    z0.s, z18.s, alphaV0
        st1w    {z0.s}, p0, [pCRow10]
        add     pCRow10, pCRow10, #16

        fmla    z4.s, z19.s, alphaV0
        st1w    {z4.s}, p0, [pCRow11]
        add     pCRow11, pCRow11, #16
.endm

/******************************************************************************/
.macro INIT4x8
        fmov    s8,  wzr
        fmov    s9, s8
        fmov    s10, s8
        fmov    s11, wzr
        fmov    s12, wzr
        fmov    s13, s9
        fmov    s14, s10
        fmov    s15, s13
.endm

.macro KERNEL4x8_I
        ld1rqw  {z0.s}, p0/z, [pA]

        ld1rqw  {z2.s}, p0/z, [pB]
        ld1rqw  {z3.s}, p0/z, [pB, #16]

        fmul    z8.s, z0.s, z2.s[0]
        fmul    z9.s, z0.s, z2.s[1]
        fmul    z10.s, z0.s, z2.s[2]
        fmul    z11.s, z0.s, z2.s[3]

        ld1rqw  {z1.s}, p0/z, [pA, #16]

        fmul    z12.s, z0.s, z3.s[0]
        fmul    z13.s, z0.s, z3.s[1]
        fmul    z14.s, z0.s, z3.s[2]
        fmul    z15.s, z0.s, z3.s[3]

        ld1rqw  {z5.s}, p0/z, [pB, #32]
        ld1rqw  {z6.s}, p0/z, [pB, #48]

        add     pA, pA, #32 
        add     pB, pB, #64 
.endm


.macro KERNEL4x8_M1
        fmla    z8.s, z0.s, z2.s[0]
        fmla    z9.s, z0.s, z2.s[1]
        fmla    z10.s, z0.s, z2.s[2]
        fmla    z11.s, z0.s, z2.s[3]
        fmla    z12.s, z0.s, z3.s[0]
        fmla    z13.s, z0.s, z3.s[1]

        ld1rqw  {z1.s}, p0/z, [pA]

        fmla    z14.s, z0.s, z3.s[2]
        fmla    z15.s, z0.s, z3.s[3]

        ld1rqw  {z5.s}, p0/z, [pB]
        ld1rqw  {z6.s}, p0/z, [pB, #16]

        add     pA, pA, #16
        add     pB, pB, #32
.endm

.macro KERNEL4x8_M2
        fmla    z8.s, z1.s, z5.s[0]
        fmla    z9.s, z1.s, z5.s[1]
        fmla    z10.s, z1.s, z5.s[2]
        fmla    z11.s, z1.s, z5.s[3]
        fmla    z12.s, z1.s, z6.s[0]
        fmla    z13.s, z1.s, z6.s[1]

        ld1rqw  {z2.s}, p0/z, [pB]
        ld1rqw  {z3.s}, p0/z, [pB, #16]

        fmla    z14.s, z1.s, z6.s[2]
        fmla    z15.s, z1.s, z6.s[3]

        ld1rqw  {z0.s}, p0/z, [pA]
        add     pA, pA, #16
        add     pB, pB, #32
.endm

.macro KERNEL4x8_E
        fmla    z8.s, z1.s, z5.s[0]
        fmla    z9.s, z1.s, z5.s[1]
        fmla    z10.s, z1.s, z5.s[2]
        fmla    z11.s, z1.s, z5.s[3]
        fmla    z12.s, z1.s, z6.s[0]
        fmla    z13.s, z1.s, z6.s[1]
        fmla    z14.s, z1.s, z6.s[2]
        fmla    z15.s, z1.s, z6.s[3]
.endm

.macro KERNEL4x8_SUB
        ld1rqw  {z0.s}, p0/z, [pA]
        ld1rqw  {z2.s}, p0/z, [pB]
        ld1rqw  {z3.s}, p0/z, [pB, #16]

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
        dup     alpha0,   alpha

        ld1rqw  {z0.s}, p0/z, [pCRow0]
        fmla    z0.s, z8.s, alphaV0
        st1w    {z0.s}, p0, [pCRow0]
        add     pCRow0, pCRow0, #16
        
        ld1rqw  {z4.s}, p0/z, [pCRow1]
        fmla    z4.s, z9.s, alphaV0
        st1w    {z4.s}, p0, [pCRow1]
        add     pCRow1, pCRow1, #16

        ld1rqw  {z0.s}, p0/z, [pCRow2]
        fmla    z0.s, z10.s, alphaV0
        st1w    {z0.s}, p0, [pCRow2]
        add     pCRow2, pCRow2, #16

        ld1rqw  {z4.s}, p0/z, [pCRow3]
        fmla    z4.s, z11.s, alphaV0
        st1w    {z4.s}, p0, [pCRow3]
        add     pCRow3, pCRow3, #16

        ld1rqw  {z0.s}, p0/z, [pCRow4]
        fmla    z0.s, z12.s, alphaV0
        st1w    {z0.s}, p0, [pCRow4]
        add     pCRow4, pCRow4, #16

        ld1rqw  {z4.s}, p0/z, [pCRow5]
        fmla    z4.s, z13.s, alphaV0
        st1w    {z4.s}, p0, [pCRow5]
        add     pCRow5, pCRow5, #16

        ld1rqw  {z0.s}, p0/z, [pCRow6]
        fmla    z0.s, z14.s, alphaV0
        st1w    {z0.s}, p0, [pCRow6]
        add     pCRow6, pCRow6, #16

        ld1rqw  {z4.s}, p0/z, [pCRow7]
        fmla    z4.s, z15.s, alphaV0
        st1w    {z4.s}, p0, [pCRow7]
        add     pCRow7, pCRow7, #16
.endm

/******************************************************************************/
.macro INIT4x4
        fmov    s8,  wzr
        fmov    s9,  wzr
        fmov    s10, s9
        fmov    s11, s8
.endm

.macro KERNEL4x4_I
        ld1rqw  {z0.s}, p0/z, [pA]
        ld1rqw  {z4.s}, p0/z, [pB]

        fmul    z8.s, z0.s, z4.s[0]
        fmul    z9.s, z0.s, z4.s[1]
        fmul    z10.s, z0.s, z4.s[2]
        fmul    z11.s, z0.s, z4.s[3]

        ld1rqw  {z1.s}, p0/z, [pA, #16]
        ld1rqw  {z5.s}, p0/z, [pB,#16]

        add     pA, pA, #32
        add     pB, pB, #32
.endm

.macro KERNEL4x4_M1
        fmla    z8.s, z0.s, z4.s[0]
        fmla    z9.s, z0.s, z4.s[1]
        fmla    z10.s, z0.s, z4.s[2]
        fmla    z11.s, z0.s, z4.s[3]

        ld1rqw  {z1.s}, p0/z, [pA]
        ld1rqw  {z5.s}, p0/z, [pB]

        add     pA, pA, #16
        add     pB, pB, #16
.endm

.macro KERNEL4x4_M2
        fmla    z8.s, z1.s, z5.s[0]
        fmla    z9.s, z1.s, z5.s[1]
        fmla    z10.s, z1.s, z5.s[2]
        fmla    z11.s, z1.s, z5.s[3]

        ld1rqw  {z4.s}, p0/z, [pB]
        ld1rqw  {z0.s}, p0/z, [pA]

        add     pA, pA, #16
        add     pB, pB, #16
.endm

.macro KERNEL4x4_E
        fmla    z8.s, z1.s, z5.s[0]
        fmla    z9.s, z1.s, z5.s[1]
        fmla    z10.s, z1.s, z5.s[2]
        fmla    z11.s, z1.s, z5.s[3]
.endm

.macro KERNEL4x4_SUB
        ld1rqw  {z0.s}, p0/z, [pA]
        ld1rqw  {z4.s}, p0/z, [pB]

        fmla    z8.s, z0.s, z4.s[0]
        fmla    z9.s, z0.s, z4.s[1]
        fmla    z10.s, z0.s, z4.s[2]
        fmla    z11.s, z0.s, z4.s[3]

        add     pA, pA, #16
        add     pB, pB, #16
.endm

.macro SAVE4x4
        dup     alpha0,   alpha

        ld1rqw  {z0.s}, p0/z, [pCRow0]
        fmla    z0.s, z8.s, alphaV0
        st1w    {z0.s}, p0, [pCRow0]
        add     pCRow0, pCRow0, #16
        
        ld1rqw    {z4.s}, p0/z, [pCRow1]
        fmla    z4.s, z9.s, alphaV0
        st1w    {z4.s}, p0, [pCRow1]
        add     pCRow1, pCRow1, #16

        ld1rqw  {z0.s}, p0/z, [pCRow2]
        fmla    z0.s, z10.s, alphaV0
        st1w    {z0.s}, p0, [pCRow2]
        add     pCRow2, pCRow2, #16

        ld1rqw  {z4.s}, p0/z, [pCRow3]
        fmla    z4.s, z11.s, alphaV0
        st1w    {z4.s}, p0, [pCRow3]
        add     pCRow3, pCRow3, #16
.endm

/******************************************************************************/
.macro INIT4x2
        fmov    s8,  wzr
        fmov    s9,  wzr
.endm

.macro KERNEL4x2_SUB
        ld1rqw  {z0.s}, p0/z, [pA]
        ld1rqw  {z4.s}, p0/z, [pB]

        fmla    z8.s, z0.s, z4.s[0]
        fmla    z9.s, z0.s, z4.s[1]

        add     pA, pA, #16
        add     pB, pB, #8
.endm

.macro SAVE4x2
        dup     alpha0,   alpha

        ld1rqw  {z0.s}, p0/z, [pCRow0]
        fmla    z0.s, z8.s, alphaV0
        st1w    {z0.s}, p0, [pCRow0]
        add     pCRow0, pCRow0, #16

        ld1rqw  {z4.s}, p0/z, [pCRow1]
        fmla    z4.s, z9.s, alphaV0
        st1w    {z4.s}, p0, [pCRow1]
        add     pCRow1, pCRow1, #16
.endm
    
/******************************************************************************/
.macro INIT4x1
        fmov    s8,  wzr
.endm

.macro KERNEL4x1_SUB
        ld1rqw {z0.s}, p0/z, [pA]
        ld1rw  {z4.s}, p0/z, [pB]

        fmla    z8.s, z0.s, z4.s[0]

        add     pA, pA, #16
        add     pB, pB, #4
.endm

.macro SAVE4x1
        dup     alpha0,   alpha

        ld1rqw  {z0.s}, p0/z, [pCRow0]
        fmla    z0.s, z8.s, alphaV0
        st1w    {z0.s}, p0, [pCRow0]
        add     pCRow0, pCRow0, #16
.endm

/******************************************************************************/
.macro INIT2x12
        fmov    s9,  wzr
        fmov    s10,  wzr
        fmov    s11, s9
        fmov    s12, s10
        fmov    s13, s9
        fmov    s14, s10
        fmov    s15, s9
        fmov    s16, s10
        fmov    s17, wzr
        fmov    s18, wzr
        fmov    s19, s16
        fmov    s20, s17
.endm

.macro KERNEL2x12_SUB
        ld1rqw  {z0.s}, p1/z, [pA]

        ld1rqw  {z2.s}, p0/z, [pB]
        ld1rqw  {z3.s}, p0/z, [pB, #16]
        ld1rqw  {z4.s}, p0/z, [pB, #32]

        fmla    z9.s,  z0.s,  z2.s[0]
        fmla    z10.s,  z0.s,  z2.s[1]
        fmla    z11.s,  z0.s,  z2.s[2]
        fmla    z12.s,  z0.s,  z2.s[3]

        fmla    z13.s,  z0.s,  z3.s[0]
        fmla    z14.s,  z0.s,  z3.s[1]
        fmla    z15.s,  z0.s,  z3.s[2]
        fmla    z16.s,  z0.s,  z3.s[3]

        fmla    z17.s,  z0.s,  z4.s[0]
        fmla    z18.s,  z0.s,  z4.s[1]
        fmla    z19.s,  z0.s,  z4.s[2]
        fmla    z20.s,  z0.s,  z4.s[3]

        add     pA, pA, #8
        add     pB, pB, #48
.endm


.macro SAVE2x12
        dup     alpha0,   alpha

        ld1rqw  {z0.s}, p1/z, [pCRow0]
        ld1rqw  {z4.s}, p1/z, [pCRow1]

        fmla    z0.s, z9.s, alphaV0
        st1w    {z0.s}, p1, [pCRow0]
        add     pCRow0, pCRow0, #8

        fmla    z4.s, z10.s, alphaV0
        st1w    {z4.s}, p1, [pCRow1]
        add     pCRow1, pCRow1, #8
        
        ld1rqw  {z0.s}, p1/z, [pCRow2]
        ld1rqw  {z4.s}, p1/z, [pCRow3]

        fmla    z0.s, z11.s, alphaV0
        st1w    {z0.s}, p1, [pCRow2]
        add     pCRow2, pCRow2, #8

        fmla    z4.s, z12.s, alphaV0
        st1w    {z4.s}, p1, [pCRow3]
        add     pCRow3, pCRow3, #8

        ld1rqw  {z0.s}, p1/z, [pCRow4]
        ld1rqw  {z4.s}, p1/z, [pCRow5]

        fmla    z0.s, z13.s, alphaV0
        st1w    {z0.s}, p1, [pCRow4]
        add     pCRow4, pCRow4, #8

        fmla    z4.s, z14.s, alphaV0
        st1w    {z4.s}, p1, [pCRow5]
        add     pCRow5, pCRow5, #8

        ld1rqw  {z0.s}, p1/z, [pCRow6]
        ld1rqw  {z4.s}, p1/z, [pCRow7]

        fmla    z0.s, z15.s, alphaV0
        st1w    {z0.s}, p1, [pCRow6]
        add     pCRow6, pCRow6, #8

        fmla    z4.s, z16.s, alphaV0
        st1w    {z4.s}, p1, [pCRow7]
        add     pCRow7, pCRow7, #8

        ld1rqw  {z0.s}, p1/z, [pCRow8]
        ld1rqw  {z4.s}, p1/z, [pCRow9]

        fmla    z0.s, z17.s, alphaV0
        st1w    {z0.s}, p1, [pCRow8]
        add     pCRow8, pCRow8, #8

        fmla    z4.s, z18.s, alphaV0
        st1w    {z4.s}, p1, [pCRow9]
        add     pCRow9, pCRow9, #8
        
        ld1rqw  {z0.s}, p1/z, [pCRow10]
        ld1rqw  {z4.s}, p1/z, [pCRow11]

        fmla    z0.s, z19.s, alphaV0
        st1w    {z0.s}, p1, [pCRow10]
        add     pCRow10, pCRow10, #8

        fmla    z4.s, z20.s, alphaV0
        st1w    {z4.s}, p1, [pCRow11]
        add     pCRow11, pCRow11, #8
.endm

/******************************************************************************/
.macro INIT2x8
        fmov    s9,  wzr
        fmov    s10, wzr
        fmov    s11, s9
        fmov    s12, wzr
        fmov    s13, wzr
        fmov    s14, s11
        fmov    s15, s13
        fmov    s16, s10
.endm


.macro KERNEL2x8_SUB
        ld1rqw  z0.s, p1/z, [pA]
        ld1rqw  z2.s, p0/z, [pB]
        ld1rqw  z3.s, p0/z, [pB, #16]

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
        dup     alpha0,   alpha

        ld1rqw  {z0.s}, p1/z, [pCRow0]
        fmla    z0.s, z9.s, alphaV0
        st1w    {z0.s}, p1, [pCRow0]
        add     pCRow0, pCRow0, #8

        ld1rqw  {z4.s}, p1/z, [pCRow1]
        fmla    z4.s, z10.s, alphaV0
        st1w    {z4.s}, p1, [pCRow1]
        add     pCRow1, pCRow1, #8

        ld1rqw  {z0.s}, p1/z, [pCRow2]
        fmla    z0.s, z11.s, alphaV0
        st1w    {z0.s}, p1, [pCRow2]
        add     pCRow2, pCRow2, #8

        ld1rqw  {z4.s}, p1/z, [pCRow3]
        fmla    z4.s, z12.s, alphaV0
        st1w    {z4.s}, p1, [pCRow3]
        add     pCRow3, pCRow3, #8

        ld1rqw  {z0.s}, p1/z, [pCRow4]
        fmla    z0.s, z13.s, alphaV0
        st1w    {z0.s}, p1, [pCRow4]
        add     pCRow4, pCRow4, #8

        ld1rqw    {z4.s}, p1/z, [pCRow5]
        fmla    z4.s, z14.s, alphaV0
        st1w    {z4.s}, p1, [pCRow5]
        add    pCRow5, pCRow5, #8

        ld1rqw  {z0.s}, p1/z, [pCRow6]
        fmla    z0.s, z15.s, alphaV0
        st1w    {z0.s}, p1, [pCRow6]
        add     pCRow6, pCRow6, #8

        ld1rqw  {z4.s}, p1/z, [pCRow7]
        fmla    z4.s, z16.s, alphaV0
        st1w    {z4.s}, p1, [pCRow7]
        add     pCRow7, pCRow7, #8
.endm

/******************************************************************************/
.macro INIT2x4
        fmov    s9,  wzr
        fmov    s10,  wzr
        fmov    s11, s9
        fmov    s12, s10
.endm

.macro KERNEL2x4_SUB
        ld1rqw  {z0.s}, p1/z, [pA]
        ld1rqw  {z4.s}, p0/z, [pB]

        fmla    z9.s, z0.s, z4.s[0]
        fmla    z10.s, z0.s, z4.s[1]
        fmla    z11.s, z0.s, z4.s[2]
        fmla    z12.s, z0.s, z4.s[3]

        add     pA, pA, #8
        add     pB, pB, #16
.endm

.macro SAVE2x4
        dup     alpha0,   alpha

        ld1rqw  {z0.s}, p1/z, [pCRow0]
        fmla    z0.s, z9.s, alphaV0
        st1w    {z0.s}, p1, [pCRow0]
        add     pCRow0, pCRow0, #8

        ld1rqw  {z4.s}, p1/z, [pCRow1]
        fmla    z4.s, z10.s, alphaV0
        st1w    {z4.s}, p1, [pCRow1]
        add     pCRow1, pCRow1, #8

        ld1rqw  {z0.s}, p1/z, [pCRow2]
        fmla    z0.s, z11.s, alphaV0
        st1w    {z0.s}, p1, [pCRow2]
        add     pCRow2, pCRow2, #8

        ld1rqw  {z4.s}, p1/z, [pCRow3]
        fmla    z4.s, z12.s, alphaV0
        st1w    {z4.s}, p1, [pCRow3]
        add     pCRow3, pCRow3, #8
.endm

/******************************************************************************/
.macro INIT2x2
        fmov    s8,  wzr
        fmov    s9,  wzr
.endm

.macro KERNEL2x2_SUB
        ld1rqw  {z0.s}, p1/z, [pA]
        ld1rqw  {z4.s}, p1/z, [pB]

        fmla    z8.s, z0.s, z4.s[0]
        fmla    z9.s, z0.s, z4.s[1]

        add     pA, pA, #8
        add     pB, pB, #8
.endm

.macro SAVE2x2
        dup     alpha0,   alpha

        ld1rqw  {z0.s}, p1/z, [pCRow0]
        fmla    z0.s, z8.s, alphaV0
        st1w    {z0.s}, p1, [pCRow0]
        add     pCRow0, pCRow0, #8
        
        ld1rqw  {z4.s}, p1/z, [pCRow1]
        fmla    z4.s, z9.s, alphaV0
        st1w    {z4.s}, p1, [pCRow1]
        add     pCRow1, pCRow1, #8
.endm
    
/******************************************************************************/
.macro INIT2x1
        fmov    s8,  wzr
.endm

.macro KERNEL2x1_SUB
        ld1rqw  {z0.s}, p1/z, [pA]
        ld1rw   {z4.s}, p1/z, [pB]

        fmla    z8.s, z0.s, z4.s[0]

        add     pA, pA, #8
        add     pB, pB, #4
.endm

.macro SAVE2x1
        dup     alpha0,   alpha

        ld1rqw  {z0.s}, p1/z, [pCRow0]
        fmla    z0.s, z8.s, alphaV0
        st1w    {z0.s}, p1, [pCRow0]
        add     pCRow0, pCRow0, #8
.endm

/******************************************************************************/
.macro INIT1x12
        fmov    s8,  wzr
        fmov    s9,  wzr
        fmov    s10, s8
        fmov    s11, s9
        fmov    s12, s8
        fmov    s13, s9
        fmov    s14, s8
        fmov    s15, s9
        fmov    s16, wzr
        fmov    s17, wzr
        fmov    s18, s16
        fmov    s19, s17
.endm

.macro KERNEL1x12_SUB
        ld1rw   {z0.s}, p1/z, [pA]

        ld1rqw  {z2.s}, p0/z, [pB]
        ld1rqw  {z3.s}, p0/z, [pB, #16]
        ld1rqw  {z4.s}, p0/z, [pB, #32]

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

        add     pA, pA, #4
        add     pB, pB, #48
.endm


.macro SAVE1x12
        dup     alpha0,   alpha

        ld1rw    {z0.s}, p1/z, [pCRow0]
        ld1rw    {z4.s}, p1/z, [pCRow1]
 
        fmla    z0.s, z8.s, alphaV0
        st1w    {z0.s}, p1, [pCRow0]
        add     pCRow0, pCRow0, #4

        fmla    z4.s, z9.s, alphaV0
        st1w    {z4.s}, p1, [pCRow1]
        add    pCRow1, pCRow1, #4

        ld1rw    {z0.s}, p1/z, [pCRow2]
        ld1rw    {z4.s}, p1/z, [pCRow3]

        fmla    z0.s, z10.s, alphaV0
        st1w    {z0.s}, p1, [pCRow2]
        add     pCRow2, pCRow2, #4

        fmla    z4.s, z11.s, alphaV0
        st1w    {z4.s}, p1, [pCRow3]
        add     pCRow3, pCRow3, #4

        ld1rw    {z0.s}, p1/z, [pCRow4]
        ld1rw    {z4.s}, p1/z, [pCRow5]

        fmla    z0.s, z12.s, alphaV0
        st1w    {z0.s}, p1, [pCRow4]
        add     pCRow4, pCRow4, #4

        fmla    z4.s, z13.s, alphaV0
        st1w    {z4.s}, p1, [pCRow5]
        add     pCRow5, pCRow5, #4

        ld1rw    {z0.s}, p1/z, [pCRow6]
        ld1rw    {z4.s}, p1/z, [pCRow7]

        fmla    z0.s, z14.s, alphaV0
        st1w    {z0.s}, p1, [pCRow6]
        add     pCRow6, pCRow6, #4

        fmla    z4.s, z15.s, alphaV0
        st1w    {z4.s}, p1, [pCRow7]
        add     pCRow7, pCRow7, #4

        ld1rw    {z0.s}, p1/z, [pCRow8]
        ld1rw    {z4.s}, p1/z, [pCRow9]

        fmla    z0.s, z16.s, alphaV0
        st1w    {z0.s}, p1, [pCRow8]
        add     pCRow8, pCRow8, #4

        fmla    z4.s, z17.s, alphaV0
        st1w    {z4.s}, p1, [pCRow9]
        add     pCRow9, pCRow9, #4

        ld1rw    {z0.s}, p1/z, [pCRow10]
        ld1rw    {z4.s}, p1/z, [pCRow11]

        fmla    z0.s, z18.s, alphaV0
        st1w    {z0.s}, p1, [pCRow10]
        add     pCRow10, pCRow10, #4

        fmla    z4.s, z19.s, alphaV0
        st1w    {z4.s}, p1, [pCRow11]
        add     pCRow11, pCRow11, #4
.endm

/******************************************************************************/
.macro INIT1x8
        fmov    s8,  wzr
        fmov    s9, s8
        fmov    s10, s8
        fmov    s11, wzr
        fmov    s12, wzr
        fmov    s13, s9
        fmov    s14, s11
        fmov    s15, s12
.endm

.macro KERNEL1x8_SUB
        ld1rw   {z0.s}, p1/z, [pA]

        ld1rqw  {z2.s}, p0/z, [pB]
        ld1rqw  {z3.s}, p0/z, [pB, #16]

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
        dup     alpha0,   alpha
        
        ld1rw   {z0.s}, p1/z, [pCRow0]
        fmla    z0.s, z8.s, alphaV0
        st1w    {z0.s}, p1, [pCRow0]
        add     pCRow0, pCRow0, #4

        ld1rw   {z4.s}, p1/z, [pCRow1]
        fmla    z4.s, z9.s, alphaV0
        st1w    {z4.s}, p1, [pCRow1]
        add     pCRow1, pCRow1, #4

        ld1rw   {z0.s}, p1/z, [pCRow2]
        fmla    z0.s, z10.s, alphaV0
        st1w    {z0.s}, p1, [pCRow2]
        add     pCRow2, pCRow2, #4

        ld1rw   {z4.s}, p1/z, [pCRow3]
        fmla    z4.s, z11.s, alphaV0
        st1w    {z4.s}, p1, [pCRow3]
        add     pCRow3, pCRow3, #4

        ld1rw   {z0.s}, p1/z, [pCRow4]
        fmla    z0.s, z12.s, alphaV0
        st1w    {z0.s}, p1, [pCRow4]
        add     pCRow4, pCRow4, #4

        ld1rw   {z4.s}, p1/z, [pCRow5]
        fmla    z4.s, z13.s, alphaV0
        st1w    {z4.s}, p1, [pCRow5]
        add     pCRow5, pCRow5, #4

        ld1rw   {z0.s}, p1/z, [pCRow6]
        fmla    z0.s, z14.s, alphaV0
        st1w    {z0.s}, p1, [pCRow6]
        add     pCRow6, pCRow6, #4

        ld1rw   {z4.s}, p1/z, [pCRow7]
        fmla    z4.s, z15.s, alphaV0
        st1w    {z4.s}, p1, [pCRow7]
        add     pCRow7, pCRow7, #4
.endm

/******************************************************************************/
.macro INIT1x4
        fmov    s8,  wzr
        fmov    s9,  wzr
        fmov    s10, s9
        fmov    s11, s8
.endm

.macro KERNEL1x4_SUB
        ld1rw   {z0.s}, p1/z, [pA]
        ld1rqw  {z4.s}, p0/z, [pB]

        fmla    z8.s, z0.s, z4.s[0]
        fmla    z9.s, z0.s, z4.s[1]
        fmla    z10.s, z0.s, z4.s[2]
        fmla    z11.s, z0.s, z4.s[3]

        add     pA, pA, #4
        add     pB, pB, #16
.endm

.macro SAVE1x4
        dup     alpha0,   alpha

        ld1rw   {z0.s}, p1/z, [pCRow0]
        fmla    z0.s, z8.s, alphaV0
        st1w    {z0.s}, p1, [pCRow0]
        add     pCRow0, pCRow0, #4

        ld1rw   {z4.s}, p1/z, [pCRow1]
        fmla    z4.s, z9.s, alphaV0
        st1w    {z4.s}, p1, [pCRow1]
        add     pCRow1, pCRow1, #4

        ld1rw   {z0.s}, p1/z, [pCRow2]
        fmla    z0.s, z10.s, alphaV0
        st1w    {z0.s}, p1, [pCRow2]
        add     pCRow2, pCRow2, #4

        ld1rw   {z4.s}, p1/z, [pCRow3]
        fmla    z4.s, z11.s, alphaV0
        st1w    {z4.s}, p1, [pCRow3]
        add     pCRow3, pCRow3, #4
.endm

/******************************************************************************/
.macro INIT1x2
        fmov    s8,  wzr
        fmov    s9,  wzr
.endm

.macro KERNEL1x2_SUB
        ld1rw   {z0.s}, p1/z, [pA]
        ld1rqw  {z4.s}, p0/z, [pB]

        fmla    z8.s, z0.s, z4.s[0]
        fmla    z9.s, z0.s, z4.s[1]

        add     pA, pA, #4
        add     pB, pB, #8
.endm

.macro SAVE1x2
        dup     alpha0,   alpha

        ld1rqw  {z0.s}, p1/z, [pCRow0]
        fmla    z0.s, z8.s, alphaV0
        st1w    {z0.s}, p1, [pCRow0]
        add     pCRow0, pCRow0, #4

        ld1rqw  {z4.s}, p1/z, [pCRow1]
        fmla    z4.s, z9.s, alphaV0
        st1w    {z4.s}, p1, [pCRow1]
        add     pCRow1, pCRow1, #4
.endm
    
/******************************************************************************/
.macro INIT1x1
        fmov    s16, wzr
.endm

.macro KERNEL1x1_SUB
        ld1rw   z0.s, p0/z, [pA]
        ld1rw   z8.s, p0/z, [pB]

        add     pA , pA, #4
        add     pB , pB, #4

        fmla    z16.s, p0/m, z0.s, z8.s
.endm

.macro SAVE1x1
        dup     alpha0, alpha

        ld1rw   z1.s, p0/z, [pCRow0]
        fmla    z1.s, p0/m, alpha0, z16.s
        st1w    {z1.s}, p0, [pCRow0]

        add     pCRow0, pCRow0, #4
.endm


/*******************************************************************************
* End of macro definitions
*******************************************************************************/

        PROLOGUE

.Lsgemm_kernel_begin:

        .align 4
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

        prfm    PLDL1KEEP, [origPB]
        prfm    PLDL1KEEP, [origPA]

        fmov    alpha, s0

        lsl    LDC, LDC, #2            // ldc = ldc * 4

        mov    pB, origPB

        mov    counterJ, origN
        cmp     counterJ, #11
        ble    .Lsgemm_kernel_L8_BEGIN

/******************************************************************************/
.Lsgemm_kernel_L12_BEGIN:
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

        add    pC, pCRow11, LDC

        mov    pA, origPA            // pA = start of A array

.Lsgemm_kernel_L12_M16_BEGIN:

        mov    counterI, origM
        cmp    counterI, #15
        ble    .Lsgemm_kernel_L12_M8_BEGIN

        .align 4
.Lsgemm_kernel_L12_M16_20:
        ptrue  p0.s, all

        mov    pB, origPB

        asr    counterL , origK, #3
        cmp    counterL , #2
        blt    .Lsgemm_kernel_L12_M16_32

        KERNEL16x12_I
        KERNEL16x12_M2
        KERNEL16x12_M1
        KERNEL16x12_M2
        KERNEL16x12_M1
        KERNEL16x12_M2
        KERNEL16x12_M1
        KERNEL16x12_M2

        subs    counterL, counterL, #2
        ble    .Lsgemm_kernel_L12_M16_22a

        .align 4
.Lsgemm_kernel_L12_M16_22:

        KERNEL16x12_M1
        KERNEL16x12_M2
        KERNEL16x12_M1
        KERNEL16x12_M2
        KERNEL16x12_M1
        KERNEL16x12_M2
        KERNEL16x12_M1
        KERNEL16x12_M2

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L12_M16_22

        .align 4
.Lsgemm_kernel_L12_M16_22a:

        KERNEL16x12_M1
        KERNEL16x12_M2
        KERNEL16x12_M1
        KERNEL16x12_M2
        KERNEL16x12_M1
        KERNEL16x12_M2
        KERNEL16x12_M1
        KERNEL16x12_E

        b     .Lsgemm_kernel_L12_M16_44

        .align 4
.Lsgemm_kernel_L12_M16_32:

        tst    counterL, #1
        ble    .Lsgemm_kernel_L12_M16_40

        KERNEL16x12_I
        KERNEL16x12_M2
        KERNEL16x12_M1
        KERNEL16x12_M2
        KERNEL16x12_M1
        KERNEL16x12_M2
        KERNEL16x12_M1
        KERNEL16x12_E

        b    .Lsgemm_kernel_L12_M16_44

.Lsgemm_kernel_L12_M16_40:

        INIT16x12

.Lsgemm_kernel_L12_M16_44:

        ands    counterL , origK, #7
        ble    .Lsgemm_kernel_L12_M16_100

        .align 4
.Lsgemm_kernel_L12_M16_46:

        KERNEL16x12_SUB
        subs    counterL, counterL, #1
        bne    .Lsgemm_kernel_L12_M16_46

.Lsgemm_kernel_L12_M16_100:
        prfm    PLDL1KEEP, [pA]
        prfm    PLDL1KEEP, [pA, #64]
        prfm    PLDL1KEEP, [origPB]

        SAVE16x12

.Lsgemm_kernel_L12_M16_END:
        subs    counterI, counterI, #16
        cmp    counterI, #15
        bgt    .Lsgemm_kernel_L12_M16_20

//------------------------------------------------------------------------------

.Lsgemm_kernel_L12_M8_BEGIN:

        cmp     counterI, #1
        blt    .Lsgemm_kernel_L12_END

        cmp    counterI, #8
        blt    .Lsgemm_kernel_L12_M4_BEGIN

.Lsgemm_kernel_L12_M8_20:
        ptrue   p0.s, all

        mov    pB, origPB

        asr    counterL , origK, #1        // L = K / 2
        cmp    counterL , #2            // is there at least 4 to do?
        blt    .Lsgemm_kernel_L12_M8_32

        KERNEL8x12_I                // do one in the K
        KERNEL8x12_M2                // do another in the K

        subs    counterL, counterL, #2
        ble    .Lsgemm_kernel_L12_M8_22a
        .align 4

.Lsgemm_kernel_L12_M8_22:

        KERNEL8x12_M1
        KERNEL8x12_M2

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L12_M8_22

.Lsgemm_kernel_L12_M8_22a:

        KERNEL8x12_M1
        KERNEL8x12_E

        b     .Lsgemm_kernel_L12_M8_44

.Lsgemm_kernel_L12_M8_32:

        tst    counterL, #1
        ble    .Lsgemm_kernel_L12_M8_40

        KERNEL8x12_I
        KERNEL8x12_E

        b    .Lsgemm_kernel_L12_M8_44

.Lsgemm_kernel_L12_M8_40:

        INIT8x12

.Lsgemm_kernel_L12_M8_44:

        ands    counterL , origK, #1
        ble    .Lsgemm_kernel_L12_M8_100

.Lsgemm_kernel_L12_M8_46:

        KERNEL8x12_SUB

.Lsgemm_kernel_L12_M8_100:

        SAVE8x12

.Lsgemm_kernel_L12_M8_END:

        subs    counterI, counterI, #8

//------------------------------------------------------------------------------

.Lsgemm_kernel_L12_M4_BEGIN:

        cmp    counterI ,#1
        blt    .Lsgemm_kernel_L12_END

        cmp    counterI, #4
        blt    .Lsgemm_kernel_L12_M2_BEGIN

.Lsgemm_kernel_L12_M4_20:
        ptrue  p0.s, VL4

        mov    pB, origPB

        asr     counterL , origK, #1        // L = K / 2
        cmp    counterL , #2            // is there at least 4 to do?
        blt    .Lsgemm_kernel_L12_M4_32

        KERNEL4x12_I                // do one in the K
        KERNEL4x12_M2                // do another in the K

        subs    counterL, counterL, #2
        ble    .Lsgemm_kernel_L12_M4_22a
        .align 4

.Lsgemm_kernel_L12_M4_22:

        KERNEL4x12_M1
        KERNEL4x12_M2

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L12_M4_22

.Lsgemm_kernel_L12_M4_22a:

        KERNEL4x12_M1
        KERNEL4x12_E

        b     .Lsgemm_kernel_L12_M4_44

.Lsgemm_kernel_L12_M4_32:

        tst    counterL, #1
        ble    .Lsgemm_kernel_L12_M4_40

        KERNEL4x12_I
        KERNEL4x12_E

        b    .Lsgemm_kernel_L12_M4_44

.Lsgemm_kernel_L12_M4_40:

        INIT4x12

.Lsgemm_kernel_L12_M4_44:

        ands    counterL , origK, #1
        ble    .Lsgemm_kernel_L12_M4_100

.Lsgemm_kernel_L12_M4_46:

        KERNEL4x12_SUB

.Lsgemm_kernel_L12_M4_100:

        SAVE4x12

.Lsgemm_kernel_L12_M4_END:

        subs    counterI, counterI, #4

//------------------------------------------------------------------------------

.Lsgemm_kernel_L12_M2_BEGIN:

        cmp    counterI, #1
        blt    .Lsgemm_kernel_L12_END
        cmp    counterI, #2            // counterI = counterI / 2
        blt    .Lsgemm_kernel_L12_M1_BEGIN

.Lsgemm_kernel_L12_M2_20:
        ptrue   p0.s, VL4
        ptrue   p1.s, VL2

        INIT2x12

        mov    pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL , #0
        ble    .Lsgemm_kernel_L12_M2_40

.Lsgemm_kernel_L12_M2_22:

        KERNEL2x12_SUB
        KERNEL2x12_SUB
        KERNEL2x12_SUB
        KERNEL2x12_SUB

        KERNEL2x12_SUB
        KERNEL2x12_SUB
        KERNEL2x12_SUB
        KERNEL2x12_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L12_M2_22


.Lsgemm_kernel_L12_M2_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Lsgemm_kernel_L12_M2_100

.Lsgemm_kernel_L12_M2_42:

        KERNEL2x12_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L12_M2_42

.Lsgemm_kernel_L12_M2_100:

        SAVE2x12

.Lsgemm_kernel_L12_M2_END:

        subs    counterI, counterI, #2


.Lsgemm_kernel_L12_M1_BEGIN:

        cmp    counterI, #1            // counterI = counterI % 2
        blt    .Lsgemm_kernel_L12_END

.Lsgemm_kernel_L12_M1_20:
        ptrue   p0.s, VL4
        ptrue   p1.s, VL1

        INIT1x12

        mov    pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL , #0
        ble    .Lsgemm_kernel_L12_M1_40

.Lsgemm_kernel_L12_M1_22:
        KERNEL1x12_SUB
        KERNEL1x12_SUB
        KERNEL1x12_SUB
        KERNEL1x12_SUB

        KERNEL1x12_SUB
        KERNEL1x12_SUB
        KERNEL1x12_SUB
        KERNEL1x12_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L12_M1_22

.Lsgemm_kernel_L12_M1_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Lsgemm_kernel_L12_M1_100

.Lsgemm_kernel_L12_M1_42:

        KERNEL1x12_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L12_M1_42

.Lsgemm_kernel_L12_M1_100:

        SAVE1x12

.Lsgemm_kernel_L12_END:

        lsl     temp1, origK, #5
        lsl     temp2, origK, #4
        add     temp, temp1, temp2
        add     origPB, origPB, temp   // B = B + K * 12 * 4

        subs     counterJ, counterJ , #12
        cmp     counterJ, #11
        bgt     .Lsgemm_kernel_L12_BEGIN


/******************************************************************************/
/******************************************************************************/

.Lsgemm_kernel_L8_BEGIN:
        cmp     counterJ, #0
        ble     .Lsgemm_kernel_L999

        cmp     counterJ, #8
        blt     .Lsgemm_kernel_L4_BEGIN

        mov    pCRow0, pC
        add    pCRow1, pCRow0, LDC
        add    pCRow2, pCRow1, LDC
        add    pCRow3, pCRow2, LDC
        add    pCRow4, pCRow3, LDC
        add    pCRow5, pCRow4, LDC
        add    pCRow6, pCRow5, LDC
        add    pCRow7, pCRow6, LDC

        add    pC, pCRow7, LDC

        mov    pA, origPA            // pA = start of A array

.Lsgemm_kernel_L8_M16_BEGIN:

        mov    counterI, origM
        cmp     counterI, #15
        ble    .Lsgemm_kernel_L8_M8_BEGIN

        .align 4
.Lsgemm_kernel_L8_M16_20:
        ptrue   p0.s, all

        mov    pB, origPB

        asr     counterL , origK, #3
        cmp    counterL , #2
        blt    .Lsgemm_kernel_L8_M16_32

        KERNEL16x8_I
        KERNEL16x8_M2
        KERNEL16x8_M1
        KERNEL16x8_M2
        KERNEL16x8_M1
        KERNEL16x8_M2
        KERNEL16x8_M1
        KERNEL16x8_M2

        subs    counterL, counterL, #2
        ble    .Lsgemm_kernel_L8_M16_22a

        .align 4
.Lsgemm_kernel_L8_M16_22:

        KERNEL16x8_M1
        KERNEL16x8_M2
        KERNEL16x8_M1
        KERNEL16x8_M2
        KERNEL16x8_M1
        KERNEL16x8_M2
        KERNEL16x8_M1
        KERNEL16x8_M2

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L8_M16_22

        .align 4
.Lsgemm_kernel_L8_M16_22a:

        KERNEL16x8_M1
        KERNEL16x8_M2
        KERNEL16x8_M1
        KERNEL16x8_M2
        KERNEL16x8_M1
        KERNEL16x8_M2
        KERNEL16x8_M1
        KERNEL16x8_E

        b     .Lsgemm_kernel_L8_M16_44

        .align 4
.Lsgemm_kernel_L8_M16_32:

        tst    counterL, #1
        ble    .Lsgemm_kernel_L8_M16_40

        KERNEL16x8_I
        KERNEL16x8_M2
        KERNEL16x8_M1
        KERNEL16x8_M2
        KERNEL16x8_M1
        KERNEL16x8_M2
        KERNEL16x8_M1
        KERNEL16x8_E

        b    .Lsgemm_kernel_L8_M16_44

.Lsgemm_kernel_L8_M16_40:

        INIT16x8

.Lsgemm_kernel_L8_M16_44:

        ands    counterL , origK, #7
        ble    .Lsgemm_kernel_L8_M16_100

        .align 4
.Lsgemm_kernel_L8_M16_46:

        KERNEL16x8_SUB
        subs    counterL, counterL, #1
        bne    .Lsgemm_kernel_L8_M16_46

.Lsgemm_kernel_L8_M16_100:
        prfm    PLDL1KEEP, [pA]
        prfm    PLDL1KEEP, [pA, #64]
        prfm    PLDL1KEEP, [origPB]

        SAVE16x8

.Lsgemm_kernel_L8_M16_END:
        subs    counterI, counterI, #16
        cmp    counterI, #15
        bgt    .Lsgemm_kernel_L8_M16_20

//------------------------------------------------------------------------------

.Lsgemm_kernel_L8_M8_BEGIN:

        cmp counterI, #1
        blt    .Lsgemm_kernel_L8_END

        cmp    counterI, #8
        blt    .Lsgemm_kernel_L8_M4_BEGIN

.Lsgemm_kernel_L8_M8_20:
        ptrue   p0.s, all

        mov    pB, origPB

        asr     counterL , origK, #1        // L = K / 2
        cmp    counterL , #2            // is there at least 4 to do?
        blt    .Lsgemm_kernel_L8_M8_32

        KERNEL8x8_I                // do one in the K
        KERNEL8x8_M2                // do another in the K

        subs    counterL, counterL, #2
        ble    .Lsgemm_kernel_L8_M8_22a
        .align 4

.Lsgemm_kernel_L8_M8_22:

        KERNEL8x8_M1
        KERNEL8x8_M2

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L8_M8_22

.Lsgemm_kernel_L8_M8_22a:

        KERNEL8x8_M1
        KERNEL8x8_E

        b     .Lsgemm_kernel_L8_M8_44

.Lsgemm_kernel_L8_M8_32:

        tst    counterL, #1
        ble    .Lsgemm_kernel_L8_M8_40

        KERNEL8x8_I
        KERNEL8x8_E

        b    .Lsgemm_kernel_L8_M8_44

.Lsgemm_kernel_L8_M8_40:

        INIT8x8

.Lsgemm_kernel_L8_M8_44:

        ands    counterL , origK, #1
        ble    .Lsgemm_kernel_L8_M8_100

.Lsgemm_kernel_L8_M8_46:

        KERNEL8x8_SUB

.Lsgemm_kernel_L8_M8_100:

        SAVE8x8

.Lsgemm_kernel_L8_M8_END:

        subs    counterI, counterI, #8

//------------------------------------------------------------------------------

.Lsgemm_kernel_L8_M4_BEGIN:

        cmp    counterI ,#1
        blt    .Lsgemm_kernel_L8_END

        cmp    counterI, #4
        blt    .Lsgemm_kernel_L8_M2_BEGIN

.Lsgemm_kernel_L8_M4_20:
        ptrue   p0.s, VL4

        mov    pB, origPB

        asr     counterL , origK, #1        // L = K / 2
        cmp    counterL , #2            // is there at least 4 to do?
        blt    .Lsgemm_kernel_L8_M4_32

        KERNEL4x8_I                // do one in the K
        KERNEL4x8_M2                // do another in the K

        subs    counterL, counterL, #2
        ble    .Lsgemm_kernel_L8_M4_22a
        .align 4

.Lsgemm_kernel_L8_M4_22:

        KERNEL4x8_M1
        KERNEL4x8_M2

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L8_M4_22

.Lsgemm_kernel_L8_M4_22a:

        KERNEL4x8_M1
        KERNEL4x8_E

        b     .Lsgemm_kernel_L8_M4_44

.Lsgemm_kernel_L8_M4_32:

        tst    counterL, #1
        ble    .Lsgemm_kernel_L8_M4_40

        KERNEL4x8_I
        KERNEL4x8_E

        b    .Lsgemm_kernel_L8_M4_44

.Lsgemm_kernel_L8_M4_40:

        INIT4x8

.Lsgemm_kernel_L8_M4_44:

        ands    counterL , origK, #1
        ble    .Lsgemm_kernel_L8_M4_100

.Lsgemm_kernel_L8_M4_46:

        KERNEL4x8_SUB

.Lsgemm_kernel_L8_M4_100:

        SAVE4x8

.Lsgemm_kernel_L8_M4_END:

        subs    counterI, counterI, #4

//------------------------------------------------------------------------------

.Lsgemm_kernel_L8_M2_BEGIN:

        cmp    counterI, #1
        blt    .Lsgemm_kernel_L8_END
        cmp    counterI, #2            // counterI = counterI / 2
        blt    .Lsgemm_kernel_L8_M1_BEGIN

.Lsgemm_kernel_L8_M2_20:
        ptrue   p0.s, VL4
        ptrue   p1.s, VL2

        INIT2x8

        mov    pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL , #0
        ble    .Lsgemm_kernel_L8_M2_40

.Lsgemm_kernel_L8_M2_22:

        KERNEL2x8_SUB
        KERNEL2x8_SUB
        KERNEL2x8_SUB
        KERNEL2x8_SUB

        KERNEL2x8_SUB
        KERNEL2x8_SUB
        KERNEL2x8_SUB
        KERNEL2x8_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L8_M2_22

.Lsgemm_kernel_L8_M2_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Lsgemm_kernel_L8_M2_100

.Lsgemm_kernel_L8_M2_42:

        KERNEL2x8_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L8_M2_42

.Lsgemm_kernel_L8_M2_100:

        SAVE2x8

.Lsgemm_kernel_L8_M2_END:

        subs    counterI, counterI, #2

.Lsgemm_kernel_L8_M1_BEGIN:

        cmp    counterI, #1            // counterI = counterI % 2
        blt    .Lsgemm_kernel_L8_END

.Lsgemm_kernel_L8_M1_20:
        ptrue   p0.s, VL4
        ptrue   p1.s, VL1

        INIT1x8

        mov    pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL , #0
        ble    .Lsgemm_kernel_L8_M1_40

.Lsgemm_kernel_L8_M1_22:
        KERNEL1x8_SUB
        KERNEL1x8_SUB
        KERNEL1x8_SUB
        KERNEL1x8_SUB

        KERNEL1x8_SUB
        KERNEL1x8_SUB
        KERNEL1x8_SUB
        KERNEL1x8_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L8_M1_22

.Lsgemm_kernel_L8_M1_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Lsgemm_kernel_L8_M1_100

.Lsgemm_kernel_L8_M1_42:

        KERNEL1x8_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L8_M1_42

.Lsgemm_kernel_L8_M1_100:

        SAVE1x8

.Lsgemm_kernel_L8_END:
        add    origPB, origPB, origK, lsl #5    // B = B + K * 4 * 8

        subs    counterJ, counterJ , #8
        bgt    .Lsgemm_kernel_L8_BEGIN


/******************************************************************************/
/******************************************************************************/

.Lsgemm_kernel_L4_BEGIN:

        cmp     counterJ, #0
        ble     .Lsgemm_kernel_L999

        cmp     counterJ, #4
        blt     .Lsgemm_kernel_L2_BEGIN

        mov    pCRow0, pC
        add    pCRow1, pCRow0, LDC
        add    pCRow2, pCRow1, LDC
        add    pCRow3, pCRow2, LDC

        add    pC, pCRow3, LDC

        mov    pA, origPA            // pA = start of A array

.Lsgemm_kernel_L4_M16_BEGIN:

        mov    counterI, origM
        cmp    counterI, #15
        ble    .Lsgemm_kernel_L4_M8_BEGIN

        .align 4
.Lsgemm_kernel_L4_M16_20:
        ptrue  p0.s, all

        mov    pB, origPB

        asr    counterL , origK, #3
        cmp    counterL , #2
        blt    .Lsgemm_kernel_L4_M16_32

        KERNEL16x4_I
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_M2

        subs    counterL, counterL, #2
        ble    .Lsgemm_kernel_L4_M16_22a

        .align 4
.Lsgemm_kernel_L4_M16_22:

        KERNEL16x4_M1
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_M2

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L4_M16_22

        .align 4
.Lsgemm_kernel_L4_M16_22a:

        KERNEL16x4_M1
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_E

        b     .Lsgemm_kernel_L4_M16_44

        .align 4
.Lsgemm_kernel_L4_M16_32:

        tst    counterL, #1
        ble    .Lsgemm_kernel_L4_M16_40

        KERNEL16x4_I
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_M2
        KERNEL16x4_M1
        KERNEL16x4_E

        b    .Lsgemm_kernel_L4_M16_44

.Lsgemm_kernel_L4_M16_40:

        INIT16x4

.Lsgemm_kernel_L4_M16_44:

        ands    counterL , origK, #7
        ble    .Lsgemm_kernel_L4_M16_100

        .align 4
.Lsgemm_kernel_L4_M16_46:

        KERNEL16x4_SUB
        subs    counterL, counterL, #1
        bne    .Lsgemm_kernel_L4_M16_46

.Lsgemm_kernel_L4_M16_100:
        prfm    PLDL1KEEP, [pA]
        prfm    PLDL1KEEP, [pA, #64]
        prfm    PLDL1KEEP, [origPB]

        SAVE16x4

.Lsgemm_kernel_L4_M16_END:
        subs    counterI, counterI, #16
        cmp     counterI, #15
        bgt    .Lsgemm_kernel_L4_M16_20

//------------------------------------------------------------------------------

.Lsgemm_kernel_L4_M8_BEGIN:

        cmp    counterI, #1
        blt    .Lsgemm_kernel_L4_END

        cmp    counterI, #8
        blt    .Lsgemm_kernel_L4_M4_BEGIN

.Lsgemm_kernel_L4_M8_20:
        ptrue   p0.s, all

        mov    pB, origPB

        asr     counterL , origK, #1        // L = K / 2
        cmp    counterL , #2            // is there at least 4 to do?
        blt    .Lsgemm_kernel_L4_M8_32

        KERNEL8x4_I                // do one in the K
        KERNEL8x4_M2                // do another in the K

        subs    counterL, counterL, #2
        ble    .Lsgemm_kernel_L4_M8_22a
        .align 4

.Lsgemm_kernel_L4_M8_22:

        KERNEL8x4_M1
        KERNEL8x4_M2

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L4_M8_22

.Lsgemm_kernel_L4_M8_22a:

        KERNEL8x4_M1
        KERNEL8x4_E

        b     .Lsgemm_kernel_L4_M8_44

.Lsgemm_kernel_L4_M8_32:

        tst    counterL, #1
        ble    .Lsgemm_kernel_L4_M8_40

        KERNEL8x4_I
        KERNEL8x4_E

        b    .Lsgemm_kernel_L4_M8_44

.Lsgemm_kernel_L4_M8_40:

        INIT8x4

.Lsgemm_kernel_L4_M8_44:

        ands    counterL , origK, #1
        ble    .Lsgemm_kernel_L4_M8_100

.Lsgemm_kernel_L4_M8_46:

        KERNEL8x4_SUB

.Lsgemm_kernel_L4_M8_100:

        SAVE8x4

.Lsgemm_kernel_L4_M8_END:

        subs    counterI, counterI, #8

//------------------------------------------------------------------------------

.Lsgemm_kernel_L4_M4_BEGIN:

        cmp    counterI , #1
        blt    .Lsgemm_kernel_L4_END

        cmp    counterI, #4
        blt    .Lsgemm_kernel_L4_M2_BEGIN

.Lsgemm_kernel_L4_M4_20:
        ptrue  p0.s, VL4

        mov    pB, origPB

        asr    counterL , origK, #1        // L = K / 2
        cmp    counterL , #2            // is there at least 4 to do?
        blt    .Lsgemm_kernel_L4_M4_32

        KERNEL4x4_I                // do one in the K
        KERNEL4x4_M2                // do another in the K

        subs    counterL, counterL, #2
        ble    .Lsgemm_kernel_L4_M4_22a
        .align 4

.Lsgemm_kernel_L4_M4_22:

        KERNEL4x4_M1
        KERNEL4x4_M2

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L4_M4_22

.Lsgemm_kernel_L4_M4_22a:

        KERNEL4x4_M1
        KERNEL4x4_E

        b     .Lsgemm_kernel_L4_M4_44

.Lsgemm_kernel_L4_M4_32:

        tst    counterL, #1
        ble    .Lsgemm_kernel_L4_M4_40

        KERNEL4x4_I
        KERNEL4x4_E

        b    .Lsgemm_kernel_L4_M4_44

.Lsgemm_kernel_L4_M4_40:

        INIT4x4

.Lsgemm_kernel_L4_M4_44:

        ands    counterL , origK, #1
        ble    .Lsgemm_kernel_L4_M4_100

.Lsgemm_kernel_L4_M4_46:

        KERNEL4x4_SUB

.Lsgemm_kernel_L4_M4_100:

        SAVE4x4

.Lsgemm_kernel_L4_M4_END:

        subs    counterI, counterI, #4

//------------------------------------------------------------------------------

.Lsgemm_kernel_L4_M2_BEGIN:

        cmp    counterI , #1
        blt    .Lsgemm_kernel_L4_END

        cmp    counterI, #2            // counterI = counterI / 2
        blt    .Lsgemm_kernel_L4_M1_BEGIN

.Lsgemm_kernel_L4_M2_20:
        ptrue   p0.s, VL4
        ptrue   p1.s, VL2

        INIT2x4

        mov    pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL , #0
        ble    .Lsgemm_kernel_L4_M2_40

.Lsgemm_kernel_L4_M2_22:

        KERNEL2x4_SUB
        KERNEL2x4_SUB
        KERNEL2x4_SUB
        KERNEL2x4_SUB

        KERNEL2x4_SUB
        KERNEL2x4_SUB
        KERNEL2x4_SUB
        KERNEL2x4_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L4_M2_22

.Lsgemm_kernel_L4_M2_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Lsgemm_kernel_L4_M2_100

.Lsgemm_kernel_L4_M2_42:

        KERNEL2x4_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L4_M2_42

.Lsgemm_kernel_L4_M2_100:

        SAVE2x4

.Lsgemm_kernel_L4_M2_END:

        subs    counterI, counterI, #2


.Lsgemm_kernel_L4_M1_BEGIN:

         cmp    counterI, #1            // counterI = counterI % 2
         blt    .Lsgemm_kernel_L4_END

.Lsgemm_kernel_L4_M1_20:
        ptrue   p0.s, VL4
        ptrue   p1.s, VL1

        INIT1x4

        mov    pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL , #0
        ble    .Lsgemm_kernel_L4_M1_40

.Lsgemm_kernel_L4_M1_22:
        KERNEL1x4_SUB
        KERNEL1x4_SUB
        KERNEL1x4_SUB
        KERNEL1x4_SUB

        KERNEL1x4_SUB
        KERNEL1x4_SUB
        KERNEL1x4_SUB
        KERNEL1x4_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L4_M1_22

.Lsgemm_kernel_L4_M1_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Lsgemm_kernel_L4_M1_100

.Lsgemm_kernel_L4_M1_42:

        KERNEL1x4_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L4_M1_42

.Lsgemm_kernel_L4_M1_100:

        SAVE1x4

.Lsgemm_kernel_L4_END:
        add    origPB, origPB, origK, lsl #4    // B = B + K * 4 * 4
        sub     counterJ, counterJ, #4

/******************************************************************************/

.Lsgemm_kernel_L2_BEGIN:

        mov    counterJ , origN
        tst    counterJ , #3
        ble    .Lsgemm_kernel_L999

        tst    counterJ , #2
        ble    .Lsgemm_kernel_L1_BEGIN

        mov    pCRow0, pC            // pCRow0 = pC
        add    pCRow1, pCRow0, LDC
        add    pC, pC, LDC, lsl #1

        mov    pA, origPA            // pA = A

.Lsgemm_kernel_L2_M16_BEGIN:

        mov    counterI, origM

        cmp    counterI, #16
        blt    .Lsgemm_kernel_L2_M8_BEGIN

.Lsgemm_kernel_L2_M16_20:
        ptrue   p0.s, all

        INIT16x2

        mov    pB, origPB

        asr    counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL,#0
        ble    .Lsgemm_kernel_L2_M16_40
        .align 4

.Lsgemm_kernel_L2_M16_22:
        KERNEL16x2_SUB
        KERNEL16x2_SUB
        KERNEL16x2_SUB
        KERNEL16x2_SUB

        KERNEL16x2_SUB
        KERNEL16x2_SUB
        KERNEL16x2_SUB
        KERNEL16x2_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L2_M16_22

.Lsgemm_kernel_L2_M16_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Lsgemm_kernel_L2_M16_100

.Lsgemm_kernel_L2_M16_42:

        KERNEL16x2_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L2_M16_42

.Lsgemm_kernel_L2_M16_100:

        SAVE16x2

.Lsgemm_kernel_L2_M16_END:

        subs    counterI, counterI, #16

        cmp    counterI, #15
        bgt    .Lsgemm_kernel_L2_M16_20

//------------------------------------------------------------------------------

.Lsgemm_kernel_L2_M8_BEGIN:

        cmp     counterI, #0
        ble    .Lsgemm_kernel_L2_END

        cmp     counterI, #7
        ble    .Lsgemm_kernel_L2_M4_BEGIN

.Lsgemm_kernel_L2_M8_20:
        ptrue   p0.s, all

        INIT8x2

        mov     pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp     counterL,#0
        ble     .Lsgemm_kernel_L2_M8_40
        .align 4

.Lsgemm_kernel_L2_M8_22:
        KERNEL8x2_SUB
        KERNEL8x2_SUB
        KERNEL8x2_SUB
        KERNEL8x2_SUB

        KERNEL8x2_SUB
        KERNEL8x2_SUB
        KERNEL8x2_SUB
        KERNEL8x2_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L2_M8_22


.Lsgemm_kernel_L2_M8_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Lsgemm_kernel_L2_M8_100

.Lsgemm_kernel_L2_M8_42:

        KERNEL8x2_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L2_M8_42

.Lsgemm_kernel_L2_M8_100:

        SAVE8x2

.Lsgemm_kernel_L2_M8_END:

        subs    counterI, counterI, #8

//------------------------------------------------------------------------------

.Lsgemm_kernel_L2_M4_BEGIN:

        cmp    counterI, #1
        blt    .Lsgemm_kernel_L2_END

        cmp    counterI, #4
        blt    .Lsgemm_kernel_L2_M2_BEGIN

.Lsgemm_kernel_L2_M4_20:
        ptrue   p0.s, VL4

        INIT4x2

        mov    pB, origPB

        asr    counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL,#0
        ble    .Lsgemm_kernel_L2_M4_40
        .align 4

.Lsgemm_kernel_L2_M4_22:
        KERNEL4x2_SUB
        KERNEL4x2_SUB
        KERNEL4x2_SUB
        KERNEL4x2_SUB

        KERNEL4x2_SUB
        KERNEL4x2_SUB
        KERNEL4x2_SUB
        KERNEL4x2_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L2_M4_22

.Lsgemm_kernel_L2_M4_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Lsgemm_kernel_L2_M4_100

.Lsgemm_kernel_L2_M4_42:

        KERNEL4x2_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L2_M4_42

.Lsgemm_kernel_L2_M4_100:

        SAVE4x2

.Lsgemm_kernel_L2_M4_END:

        subs    counterI, counterI, #4

//------------------------------------------------------------------------------

.Lsgemm_kernel_L2_M2_BEGIN:

        cmp    counterI , #1
        blt    .Lsgemm_kernel_L2_END

        cmp    counterI, #2            // counterI = counterI / 2
        blt    .Lsgemm_kernel_L2_M1_BEGIN

.Lsgemm_kernel_L2_M2_20:
        ptrue   p1.s, VL2

        INIT2x2

        mov    pB, origPB

        asr    counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL,#0
        ble    .Lsgemm_kernel_L2_M2_40

.Lsgemm_kernel_L2_M2_22:

        KERNEL2x2_SUB
        KERNEL2x2_SUB
        KERNEL2x2_SUB
        KERNEL2x2_SUB

        KERNEL2x2_SUB
        KERNEL2x2_SUB
        KERNEL2x2_SUB
        KERNEL2x2_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L2_M2_22

.Lsgemm_kernel_L2_M2_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Lsgemm_kernel_L2_M2_100

.Lsgemm_kernel_L2_M2_42:

        KERNEL2x2_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L2_M2_42

.Lsgemm_kernel_L2_M2_100:

        SAVE2x2

.Lsgemm_kernel_L2_M2_END:

        subs    counterI, counterI, #2

.Lsgemm_kernel_L2_M1_BEGIN:

        cmp    counterI, #1            // counterI = counterI % 2
        blt    .Lsgemm_kernel_L2_END

.Lsgemm_kernel_L2_M1_20:
        ptrue   p0.s, VL2
        ptrue   p1.s, VL1

        INIT1x2

        mov    pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp     counterL, #0
        ble    .Lsgemm_kernel_L2_M1_40

.Lsgemm_kernel_L2_M1_22:
        KERNEL1x2_SUB
        KERNEL1x2_SUB
        KERNEL1x2_SUB
        KERNEL1x2_SUB

        KERNEL1x2_SUB
        KERNEL1x2_SUB
        KERNEL1x2_SUB
        KERNEL1x2_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L2_M1_22

.Lsgemm_kernel_L2_M1_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Lsgemm_kernel_L2_M1_100

.Lsgemm_kernel_L2_M1_42:

        KERNEL1x2_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L2_M1_42

.Lsgemm_kernel_L2_M1_100:

        SAVE1x2

.Lsgemm_kernel_L2_END:
        add    origPB, origPB, origK, lsl #3    // B = B + K * 2 * 4

/******************************************************************************/

.Lsgemm_kernel_L1_BEGIN:

        mov    counterJ , origN
        tst    counterJ , #1
        ble    .Lsgemm_kernel_L999 // done


        mov    pCRow0, pC            // pCRow0 = C
        add    pC , pC , LDC            // Update pC to point to next

        mov    pA, origPA            // pA = A

.Lsgemm_kernel_L1_M16_BEGIN:

        mov    counterI, origM

        cmp    counterI, #16
        blt    .Lsgemm_kernel_L1_M8_BEGIN

.Lsgemm_kernel_L1_M16_20:
        ptrue   p0.s, all

        INIT16x1

        mov    pB, origPB

        asr    counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL , #0
        ble    .Lsgemm_kernel_L1_M16_40
        .align 4

.Lsgemm_kernel_L1_M16_22:
        KERNEL16x1_SUB
        KERNEL16x1_SUB
        KERNEL16x1_SUB
        KERNEL16x1_SUB

        KERNEL16x1_SUB
        KERNEL16x1_SUB
        KERNEL16x1_SUB
        KERNEL16x1_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L1_M16_22

.Lsgemm_kernel_L1_M16_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Lsgemm_kernel_L1_M16_100

.Lsgemm_kernel_L1_M16_42:

        KERNEL16x1_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L1_M16_42

.Lsgemm_kernel_L1_M16_100:

        SAVE16x1

.Lsgemm_kernel_L1_M16_END:

        subs    counterI, counterI, #16

        cmp    counterI, #15
        bgt    .Lsgemm_kernel_L1_M16_20

//------------------------------------------------------------------------------

.Lsgemm_kernel_L1_M8_BEGIN:

        cmp    counterI, #1
        blt    .Lsgemm_kernel_L1_END

        cmp    counterI, #8
        blt    .Lsgemm_kernel_L1_M4_BEGIN

.Lsgemm_kernel_L1_M8_20:
        ptrue   p0.s, all

        INIT8x1

        mov    pB, origPB

        asr    counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL , #0
        ble    .Lsgemm_kernel_L1_M8_40
        .align 4

.Lsgemm_kernel_L1_M8_22:
        KERNEL8x1_SUB
        KERNEL8x1_SUB
        KERNEL8x1_SUB
        KERNEL8x1_SUB

        KERNEL8x1_SUB
        KERNEL8x1_SUB
        KERNEL8x1_SUB
        KERNEL8x1_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L1_M8_22

.Lsgemm_kernel_L1_M8_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Lsgemm_kernel_L1_M8_100

.Lsgemm_kernel_L1_M8_42:

        KERNEL8x1_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L1_M8_42

.Lsgemm_kernel_L1_M8_100:

        SAVE8x1

.Lsgemm_kernel_L1_M8_END:

        subs    counterI, counterI, #8

//------------------------------------------------------------------------------

.Lsgemm_kernel_L1_M4_BEGIN:

        cmp    counterI, #1
        blt    .Lsgemm_kernel_L1_END

        cmp    counterI, #4
        blt    .Lsgemm_kernel_L1_M2_BEGIN

.Lsgemm_kernel_L1_M4_20:
        ptrue   p0.s, VL4

        INIT4x1

        mov    pB, origPB

        asr    counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL , #0
        ble    .Lsgemm_kernel_L1_M4_40
        .align 4

.Lsgemm_kernel_L1_M4_22:
        KERNEL4x1_SUB
        KERNEL4x1_SUB
        KERNEL4x1_SUB
        KERNEL4x1_SUB

        KERNEL4x1_SUB
        KERNEL4x1_SUB
        KERNEL4x1_SUB
        KERNEL4x1_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L1_M4_22


.Lsgemm_kernel_L1_M4_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Lsgemm_kernel_L1_M4_100

.Lsgemm_kernel_L1_M4_42:

        KERNEL4x1_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L1_M4_42

.Lsgemm_kernel_L1_M4_100:

        SAVE4x1

.Lsgemm_kernel_L1_M4_END:

        subs    counterI, counterI, #4

//------------------------------------------------------------------------------

.Lsgemm_kernel_L1_M2_BEGIN:

        cmp    counterI, #1
        blt    .Lsgemm_kernel_L1_END

        cmp    counterI, #2            // counterI = counterI / 2
        blt    .Lsgemm_kernel_L1_M1_BEGIN

.Lsgemm_kernel_L1_M2_20:
        ptrue   p1.s, VL2

        INIT2x1

        mov    pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL , #0
        ble    .Lsgemm_kernel_L1_M2_40

.Lsgemm_kernel_L1_M2_22:

        KERNEL2x1_SUB
        KERNEL2x1_SUB
        KERNEL2x1_SUB
        KERNEL2x1_SUB

        KERNEL2x1_SUB
        KERNEL2x1_SUB
        KERNEL2x1_SUB
        KERNEL2x1_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L1_M2_22

.Lsgemm_kernel_L1_M2_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Lsgemm_kernel_L1_M2_100

.Lsgemm_kernel_L1_M2_42:

        KERNEL2x1_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L1_M2_42

.Lsgemm_kernel_L1_M2_100:

        SAVE2x1

.Lsgemm_kernel_L1_M2_END:

        subs    counterI, counterI, #2

.Lsgemm_kernel_L1_M1_BEGIN:

        cmp    counterI, #1            // counterI = counterI % 2
        blt    .Lsgemm_kernel_L1_END

.Lsgemm_kernel_L1_M1_20:
        ptrue   p0.s, VL1

        INIT1x1

        mov    pB, origPB

        asr     counterL , origK, #3        // counterL = counterL / 8
        cmp    counterL , #0
        ble    .Lsgemm_kernel_L1_M1_40

.Lsgemm_kernel_L1_M1_22:
        KERNEL1x1_SUB
        KERNEL1x1_SUB
        KERNEL1x1_SUB
        KERNEL1x1_SUB

        KERNEL1x1_SUB
        KERNEL1x1_SUB
        KERNEL1x1_SUB
        KERNEL1x1_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L1_M1_22

.Lsgemm_kernel_L1_M1_40:

        ands    counterL , origK, #7        // counterL = counterL % 8
        ble    .Lsgemm_kernel_L1_M1_100

.Lsgemm_kernel_L1_M1_42:

        KERNEL1x1_SUB

        subs    counterL, counterL, #1
        bgt    .Lsgemm_kernel_L1_M1_42

.Lsgemm_kernel_L1_M1_100:

        SAVE1x1

.Lsgemm_kernel_L1_END:

.Lsgemm_kernel_L999:
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

