/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* Copyright (c) 2020, Hisilicon Limited                              */
/* All rights reserved.                                              */
/*                                                                   */
/* Based on gemm_tcopy_16.c and modified by ACLE instructions         */
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

    FLOAT *aoffset;
    FLOAT *aoffset1, *aoffset2;
    FLOAT *boffset;

    FLOAT ctemp01, ctemp02;

    svbool_t p64_all = svptrue_b64();
    svbool_t p64_2 = svwhilelt_b64(0, 2);

    svfloat64x2_t svtemp01, svtemp02, svtemp03, svtemp04;
    svfloat64_t svtemp11, svtemp12;

    aoffset = a;
    boffset = b;

    j = (n >> 4);
    if (j > 0) {
        do {
            aoffset1 = aoffset;
            aoffset2 = aoffset + lda;
            aoffset += 16;

            i = (m >> 1);
            if (i > 0) {
                do {
                    svtemp01 = svld2(p64_all, aoffset1);
                    svtemp02 = svld2_vnum(p64_all, aoffset1, 2);

                    svtemp03 = svld2(p64_all, aoffset2);
                    svtemp04 = svld2_vnum(p64_all, aoffset2, 2);

                    svst2(p64_all, boffset, svtemp01);
                    svst2_vnum(p64_all, boffset, 2, svtemp02);
                    svst2_vnum(p64_all, boffset, 4, svtemp03);
                    svst2_vnum(p64_all, boffset, 6, svtemp04);

                    aoffset1 += 2 * lda;
                    aoffset2 += 2 * lda;
                    boffset += 32;

                    i--;
                } while (i > 0);
            }

            if (m & 1) {
                svtemp01 = svld2(p64_all, aoffset1);
                svtemp02 = svld2_vnum(p64_all, aoffset1, 2);

                svst2(p64_all, boffset, svtemp01);
                svst2_vnum(p64_all, boffset, 2, svtemp02);

                boffset += 16;
            }

            j--;
        } while (j > 0);
    } /* end of if(j > 0) */

    if (n & 8) {
        aoffset1 = aoffset;
        aoffset2 = aoffset + lda;
        aoffset += 8;

        i = (m >> 1);
        if (i > 0) {
            do {
                svtemp01 = svld2(p64_all, aoffset1);
                svtemp02 = svld2(p64_all, aoffset2);

                svst2(p64_all, boffset, svtemp01);
                svst2_vnum(p64_all, boffset, 2, svtemp02);

                aoffset1 += 2 * lda;
                aoffset2 += 2 * lda;
                boffset += 16;

                i--;
            } while (i > 0);
        }

        if (m & 1) {
            svtemp01 = svld2(p64_all, aoffset1);

            svst2(p64_all, boffset, svtemp01);

            boffset += 8;
        }
    }

    if (n & 4) {
        aoffset1 = aoffset;
        aoffset2 = aoffset + lda;
        aoffset += 4;

        i = (m >> 1);
        if (i > 0) {
            do {
                svtemp11 = svld1(p64_all, aoffset1);
                svtemp12 = svld1(p64_all, aoffset2);

                svst1(p64_all, boffset, svtemp11);
                svst1_vnum(p64_all, boffset, 1, svtemp12);
                aoffset1 += 2 * lda;
                aoffset2 += 2 * lda;
                boffset += 8;

                i--;
            } while (i > 0);
        }

        if (m & 1) {
            svtemp11 = svld1(p64_all, aoffset1);

            svst1(p64_all, boffset, svtemp11);

            boffset += 4;
        }
    }

    if (n & 2) {
        aoffset1 = aoffset;
        aoffset2 = aoffset + lda;
        aoffset += 2;

        i = (m >> 1);
        if (i > 0) {
            do {
                svtemp11 = svld1rq(p64_2, aoffset1);
                svtemp12 = svld1rq(p64_2, aoffset2);

                svst1(p64_2, boffset, svtemp11);
                svst1(p64_2, boffset + 2, svtemp12);

                aoffset1 += 2 * lda;
                aoffset2 += 2 * lda;
                boffset += 4;

                i--;
            } while (i > 0);
        }

        if (m & 1) {
            svtemp11 = svld1rq(p64_2, aoffset1);

            svst1(p64_2, boffset, svtemp11);
            boffset += 2;
        }
    }

    if (n & 1) {
        aoffset1 = aoffset;
        aoffset2 = aoffset + lda;

        i = (m >> 1);
        if (i > 0) {
            do {
                ctemp01 = *(aoffset1 + 0);
                ctemp02 = *(aoffset2 + 0);

                *(boffset + 0) = ctemp01;
                *(boffset + 1) = ctemp02;

                aoffset1 += 2 * lda;
                aoffset2 += 2 * lda;
                boffset += 2;

                i--;
            } while (i > 0);
        }

        if (m & 1) {
            ctemp01 = *(aoffset1 + 0);
            *(boffset + 0) = ctemp01;
            // boffset   += 1;
        }
    }

    return 0;
}
