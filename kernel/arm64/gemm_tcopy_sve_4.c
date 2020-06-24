/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* Copyright (c) 2020, Hisilicon Limited                              */
/* All rights reserved.                                              */
/*                                                                   */
/* Based on gemm_tcopy_4.c and modified by ACLE instructions         */
/* Created by Jia Yuan <yuanjia11@huawei.com>                         */
/*            Anjun Wu <wuanjun@huawei.com>                          */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/*   1. Redistributions of source code must retain the above         */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer.                                                  */
/*                                                                   */
/*   2. Redistributions in binary form must reproduce the above      */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer in the documentation and/or other materials       */
/*      provided with the distribution.                              */
/*                                                                   */
/*    THIS  SOFTWARE IS PROVIDED  BY THE  UNIVERSITY OF  TEXAS AT    */
/*    AUSTIN  ``AS IS''  AND ANY  EXPRESS OR  IMPLIED WARRANTIES,    */
/*    INCLUDING, BUT  NOT LIMITED  TO, THE IMPLIED  WARRANTIES OF    */
/*    MERCHANTABILITY  AND FITNESS FOR  A PARTICULAR  PURPOSE ARE    */
/*    DISCLAIMED.  IN  NO EVENT SHALL THE UNIVERSITY  OF TEXAS AT    */
/*    AUSTIN OR CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT,    */
/*    INCIDENTAL,  SPECIAL, EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES    */
/*    (INCLUDING, BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE    */
/*    GOODS  OR  SERVICES; LOSS  OF  USE,  DATA,  OR PROFITS;  OR    */
/*    BUSINESS INTERRUPTION) HOWEVER CAUSED  AND ON ANY THEORY OF    */
/*    LIABILITY, WHETHER  IN CONTRACT, STRICT  LIABILITY, OR TORT    */
/*    (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY WAY OUT    */
/*    OF  THE  USE OF  THIS  SOFTWARE,  EVEN  IF ADVISED  OF  THE    */
/*    POSSIBILITY OF SUCH DAMAGE.                                    */
/*                                                                   */
/* The views and conclusions contained in the software and           */
/* documentation are those of the authors and should not be          */
/* interpreted as representing official policies, either expressed   */
/* or implied, of The University of Texas at Austin.                 */
/*********************************************************************/

#include <stdio.h>
#include "common.h"

#ifdef __ARM_FEATURE_SVE
#include <arm_sve.h>
#endif /* __ARM_FEATURE_SVE */

/* ARM SVE Intrinsic */
int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, FLOAT *b) {

    BLASLONG i, j;

    FLOAT *a_offset, *a_offset1, *a_offset2, *a_offset3, *a_offset4;
    FLOAT *b_offset, *b_offset1, *b_offset2, *b_offset3;
    FLOAT ctemp1, ctemp2, ctemp3, ctemp4;

    svbool_t p64_all = svptrue_b64();
    svbool_t p64_2 = svwhilelt_b64(0, 2);

    svfloat64_t svtemp11, svtemp12, svtemp13, svtemp14;

    a_offset = a;
    b_offset = b;

    b_offset2 = b + m * (n & ~3);
    b_offset3 = b + m * (n & ~1);

    j = (m >> 2); // j = m / 4
    if (j > 0) {
        do {
            a_offset1 = a_offset;
            a_offset2 = a_offset1 + lda;
            a_offset3 = a_offset2 + lda;
            a_offset4 = a_offset3 + lda;
            a_offset += 4 * lda;

            b_offset1 = b_offset;
            b_offset += 16;

            i = (n >> 2); // i = n / 4
            if (i > 0) {
                do { // 4 x 4
                    svtemp11 = svld1(p64_all, a_offset1);
                    svtemp12 = svld1(p64_all, a_offset2);
                    svtemp13 = svld1(p64_all, a_offset3);
                    svtemp14 = svld1(p64_all, a_offset4);

                    svst1(p64_all, b_offset1, svtemp11);
                    svst1_vnum(p64_all, b_offset1, 1, svtemp12);
                    svst1_vnum(p64_all, b_offset1, 2, svtemp13);
                    svst1_vnum(p64_all, b_offset1, 3, svtemp14);

                    a_offset1 += 4;
                    a_offset2 += 4;
                    a_offset3 += 4;
                    a_offset4 += 4;

                    b_offset1 += m * 4;
                    i--;
                } while (i > 0);
            }

            if (n & 2) { // 4 x 2
                svtemp11 = svld1rq(p64_2, a_offset1);
                svtemp12 = svld1rq(p64_2, a_offset2);
                svtemp13 = svld1rq(p64_2, a_offset3);
                svtemp14 = svld1rq(p64_2, a_offset4);

                svst1(p64_2, b_offset2, svtemp11);
                svst1(p64_2, b_offset2 + 2, svtemp12);
                svst1(p64_2, b_offset2 + 4, svtemp13);
                svst1(p64_2, b_offset2 + 6, svtemp14);

                a_offset1 += 2;
                a_offset2 += 2;
                a_offset3 += 2;
                a_offset4 += 2;
                b_offset2 += 8;
            }

            if (n & 1) { // 4 x 1
                ctemp1 = *(a_offset1 + 0);
                ctemp2 = *(a_offset2 + 0);
                ctemp3 = *(a_offset3 + 0);
                ctemp4 = *(a_offset4 + 0);

                *(b_offset3 + 0) = ctemp1;
                *(b_offset3 + 1) = ctemp2;
                *(b_offset3 + 2) = ctemp3;
                *(b_offset3 + 3) = ctemp4;

                b_offset3 += 4;
            }

            j--;
        } while (j > 0);
    }

    if (m & 2) {
        a_offset1 = a_offset;
        a_offset2 = a_offset1 + lda;
        a_offset += 2 * lda;

        b_offset1 = b_offset;
        b_offset += 8;

        i = (n >> 2);
        if (i > 0) { // 2 x 4
            do {
                svtemp11 = svld1(p64_all, a_offset1);
                svtemp12 = svld1(p64_all, a_offset2);
                svst1(p64_all, b_offset1, svtemp11);
                svst1_vnum(p64_all, b_offset1, 1, svtemp12);

                a_offset1 += 4;
                a_offset2 += 4;
                b_offset1 += m * 4;
                i--;
            } while (i > 0);
        }

        if (n & 2) { // 2 x 2
            svtemp11 = svld1rq(p64_2, a_offset1);
            svtemp12 = svld1rq(p64_2, a_offset2);
            svst1(p64_2, b_offset2, svtemp11);
            svst1(p64_2, b_offset2 + 2, svtemp12);
            a_offset1 += 2;
            a_offset2 += 2;
            b_offset2 += 4;
        }

        if (n & 1) { // 2 x 1
            ctemp1 = *(a_offset1 + 0);
            ctemp2 = *(a_offset2 + 0);

            *(b_offset3 + 0) = ctemp1;
            *(b_offset3 + 1) = ctemp2;
            b_offset3 += 2;
        }
    }

    if (m & 1) {
        a_offset1 = a_offset;
        b_offset1 = b_offset;

        i = (n >> 2);
        if (i > 0) { // 1 x 4
            do {
                svtemp11 = svld1(p64_all, a_offset1);
                svst1(p64_all, b_offset1, svtemp11);
                a_offset1 += 4;
                b_offset1 += 4 * m;
                i--;
            } while (i > 0);
        }

        if (n & 2) { // 1 x 2

            svtemp11 = svld1rq(p64_2, a_offset1);
            svst1(p64_2, b_offset2, svtemp11);
            a_offset1 += 2;
        }

        if (n & 1) { // 1 x 1
            ctemp1 = *(a_offset1 + 0);
            *(b_offset3 + 0) = ctemp1;
        }
    }

    return 0;
}
