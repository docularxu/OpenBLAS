/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* Copyright (c) 2020, Hisilicon Limited                              */
/* All rights reserved.                                              */
/*                                                                   */
/* Based on gemm_ncopy_4.c and modified by ACLE instructions         */
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
    FLOAT *a_offset, *aoffset1, *aoffset2, *aoffset3, *aoffset4;
    FLOAT *b_offset;
    FLOAT ctemp1, ctemp3, ctemp4, ctemp5, ctemp9, ctemp13;

    svbool_t p64_all = svptrue_b64();
    svbool_t p64_2 = svwhilelt_b64(0, 2);

    svfloat64_t svtemp01, svtemp02, svtemp03, svtemp04;

    svfloat64_t svins01, svins02, svins03, svins04;
    svfloat64_t svp01, svp02, svp03, svp04;

    a_offset = a;
    b_offset = b;

    j = (n >> 2);
    if (j > 0) {
        do {
            aoffset1 = a_offset;
            aoffset2 = aoffset1 + lda;
            aoffset3 = aoffset2 + lda;
            aoffset4 = aoffset3 + lda;
            a_offset += 4 * lda;

            i = (m >> 2);
            if (i > 0) {
                do {
                    svtemp01 = svld1(p64_all, aoffset1);
                    svtemp02 = svld1(p64_all, aoffset2);
                    svtemp03 = svld1(p64_all, aoffset3);
                    svtemp04 = svld1(p64_all, aoffset4);

                    svins01 = svuzp1(svtemp01, svtemp02);
                    svins02 = svuzp2(svtemp01, svtemp02);
                    svins03 = svuzp1(svtemp03, svtemp04);
                    svins04 = svuzp2(svtemp03, svtemp04);

                    svp01 = svuzp1(svins01, svins03);
                    svp02 = svuzp2(svins01, svins03);
                    svp03 = svuzp1(svins02, svins04);
                    svp04 = svuzp2(svins02, svins04);

                    svst1(p64_all, b_offset, svp01);
                    svst1_vnum(p64_all, b_offset, 1, svp03);
                    svst1_vnum(p64_all, b_offset, 2, svp02);
                    svst1_vnum(p64_all, b_offset, 3, svp04);

                    aoffset1 += 4;
                    aoffset2 += 4;
                    aoffset3 += 4;
                    aoffset4 += 4;

                    b_offset += 16;
                    i--;
                } while (i > 0);
            }

            if (m & 2) {
                svtemp01 = svld1rq(p64_2, aoffset1);
                svtemp02 = svld1rq(p64_2, aoffset2);
                svtemp03 = svld1rq(p64_2, aoffset3);
                svtemp04 = svld1rq(p64_2, aoffset4);

                svins01 = svzip1(svtemp01, svtemp03);
                svins02 = svzip1(svtemp02, svtemp04);

                svp01 = svzip1(svins01, svins02);
                svp02 = svzip2(svins01, svins02);

                svst1(p64_all, b_offset, svp01);
                svst1_vnum(p64_all, b_offset, 1, svp02);

                aoffset1 += 2;
                aoffset2 += 2;
                aoffset3 += 2;
                aoffset4 += 2;

                b_offset += 8;
            }

            if (m & 1) {
                ctemp1 = *(aoffset1 + 0);
                ctemp5 = *(aoffset2 + 0);
                ctemp9 = *(aoffset3 + 0);
                ctemp13 = *(aoffset4 + 0);

                *(b_offset + 0) = ctemp1;
                *(b_offset + 1) = ctemp5;
                *(b_offset + 2) = ctemp9;
                *(b_offset + 3) = ctemp13;

                b_offset += 4;
            }
            j--;
        } while (j > 0);
    } /* end of if(j > 0) */

    if (n & 2) {
        aoffset1 = a_offset;
        aoffset2 = aoffset1 + lda;
        a_offset += 2 * lda;

        i = (m >> 2);
        if (i > 0) {
            do {
                svtemp01 = svld1(p64_all, aoffset1);
                svtemp02 = svld1(p64_all, aoffset2);

                svins01 = svzip1(svtemp01, svtemp02);
                svins02 = svzip2(svtemp01, svtemp02);

                svst1(p64_all, b_offset, svins01);
                svst1_vnum(p64_all, b_offset, 1, svins02);

                aoffset1 += 4;
                aoffset2 += 4;
                b_offset += 8;
                i--;
            } while (i > 0);
        }

        if (m & 2) {
            svtemp01 = svld1rq(p64_2, aoffset1);
            svtemp02 = svld1rq(p64_2, aoffset2);

            svins01 = svzip1(svtemp01, svtemp02);

            svst1(p64_all, b_offset, svins01);

            aoffset1 += 2;
            aoffset2 += 2;

            b_offset += 4;
        }

        if (m & 1) {
            ctemp1 = *(aoffset1 + 0);
            ctemp5 = *(aoffset2 + 0);

            *(b_offset + 0) = ctemp1;
            *(b_offset + 1) = ctemp5;

            b_offset += 2;

        }
    } /* end of if(j > 0) */

    if (n & 1) {
        aoffset1 = a_offset;

        i = (m >> 2);
        if (i > 0) {
            do {
                svtemp01 = svld1(p64_all, aoffset1);
                svst1(p64_all, b_offset, svtemp01);

                aoffset1 += 4;
                b_offset += 4;
                i--;
            } while (i > 0);
        }

        i = (m & 3);
        if (i > 0) {
            do {
                ctemp1 = *(aoffset1 + 0);
                *(b_offset + 0) = ctemp1;
                aoffset1++;
                b_offset += 1;
                i--;
            } while (i > 0);
        }
    } /* end of if(j > 0) */

    return 0;
}
