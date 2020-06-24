/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* Copyright (c) 2020, Hisilicon Limited                              */
/* All rights reserved.                                              */
/* Based on gemm_ncopy_16.c and modified by ACLE instructions        */
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
    FLOAT *aoffset1, *aoffset2, *aoffset3, *aoffset4;
    FLOAT *aoffset5, *aoffset6, *aoffset7, *aoffset8;
    FLOAT *aoffset9, *aoffset10, *aoffset11, *aoffset12;
    FLOAT *aoffset13, *aoffset14, *aoffset15, *aoffset16;

    FLOAT *boffset;
    FLOAT ctemp01, ctemp02, ctemp03, ctemp05, ctemp07, ctemp09, ctemp11, ctemp13, ctemp15, ctemp17, ctemp19, ctemp21, ctemp23, ctemp25, ctemp27, ctemp29, ctemp31;

    svbool_t p64_all = svptrue_b64();
    svbool_t p64_2 = svwhilelt_b64(0, 2);

    svfloat64_t svtemp01, svtemp02, svtemp03, svtemp04, svtemp05, svtemp06, svtemp07, svtemp08, svtemp09, svtemp10, svtemp11, svtemp12, svtemp13, svtemp14, svtemp15, svtemp16;

    svfloat64_t svins01, svins02, svins03, svins04, svins05, svins06, svins07, svins08, svins09, svins10, svins11, svins12, svins13, svins14, svins15, svins16;
    svfloat64_t svp01, svp02, svp03, svp04, svp05, svp06, svp07, svp08, svp09, svp10, svp11, svp12, svp13, svp14, svp15, svp16;

    aoffset = a;
    boffset = b;

    j = (n >> 4);
    if (j > 0) {
        do {
            aoffset1 = aoffset;
            aoffset2 = aoffset1 + lda;
            aoffset3 = aoffset2 + lda;
            aoffset4 = aoffset3 + lda;
            aoffset5 = aoffset4 + lda;
            aoffset6 = aoffset5 + lda;
            aoffset7 = aoffset6 + lda;
            aoffset8 = aoffset7 + lda;
            aoffset9 = aoffset8 + lda;
            aoffset10 = aoffset9 + lda;
            aoffset11 = aoffset10 + lda;
            aoffset12 = aoffset11 + lda;
            aoffset13 = aoffset12 + lda;
            aoffset14 = aoffset13 + lda;
            aoffset15 = aoffset14 + lda;
            aoffset16 = aoffset15 + lda;
            aoffset += 16 * lda;

            i = (m >> 2);
            if (i > 0) {
                do {
                    svtemp01 = svld1(p64_all, aoffset1);
                    svtemp02 = svld1(p64_all, aoffset2);
                    svtemp03 = svld1(p64_all, aoffset3);
                    svtemp04 = svld1(p64_all, aoffset4);
                    svtemp05 = svld1(p64_all, aoffset5);
                    svtemp06 = svld1(p64_all, aoffset6);
                    svtemp07 = svld1(p64_all, aoffset7);
                    svtemp08 = svld1(p64_all, aoffset8);
                    svtemp09 = svld1(p64_all, aoffset9);
                    svtemp10 = svld1(p64_all, aoffset10);
                    svtemp11 = svld1(p64_all, aoffset11);
                    svtemp12 = svld1(p64_all, aoffset12);
                    svtemp13 = svld1(p64_all, aoffset13);
                    svtemp14 = svld1(p64_all, aoffset14);
                    svtemp15 = svld1(p64_all, aoffset15);
                    svtemp16 = svld1(p64_all, aoffset16);

                    svins01 = svuzp1(svtemp01, svtemp02);
                    svins02 = svuzp2(svtemp01, svtemp02);
                    svins03 = svuzp1(svtemp03, svtemp04);
                    svins04 = svuzp2(svtemp03, svtemp04);

                    svp01 = svuzp1(svins01, svins03);
                    svp02 = svuzp2(svins01, svins03);
                    svp03 = svuzp1(svins02, svins04);
                    svp04 = svuzp2(svins02, svins04);

                    svins05 = svuzp1(svtemp05, svtemp06);
                    svins06 = svuzp2(svtemp05, svtemp06);
                    svins07 = svuzp1(svtemp07, svtemp08);
                    svins08 = svuzp2(svtemp07, svtemp08);

                    svp05 = svuzp1(svins05, svins07);
                    svp06 = svuzp2(svins05, svins07);
                    svp07 = svuzp1(svins06, svins08);
                    svp08 = svuzp2(svins06, svins08);

                    svins09 = svuzp1(svtemp09, svtemp10);
                    svins10 = svuzp2(svtemp09, svtemp10);
                    svins11 = svuzp1(svtemp11, svtemp12);
                    svins12 = svuzp2(svtemp11, svtemp12);

                    svp09 = svuzp1(svins09, svins11);
                    svp10 = svuzp2(svins09, svins11);
                    svp11 = svuzp1(svins10, svins12);
                    svp12 = svuzp2(svins10, svins12);

                    svins13 = svuzp1(svtemp13, svtemp14);
                    svins14 = svuzp2(svtemp13, svtemp14);
                    svins15 = svuzp1(svtemp15, svtemp16);
                    svins16 = svuzp2(svtemp15, svtemp16);

                    svp13 = svuzp1(svins13, svins15);
                    svp14 = svuzp2(svins13, svins15);
                    svp15 = svuzp1(svins14, svins16);
                    svp16 = svuzp2(svins14, svins16);

                    svst1(p64_all, boffset, svp01);
                    svst1_vnum(p64_all, boffset, 1, svp05);
                    svst1_vnum(p64_all, boffset, 2, svp09);
                    svst1_vnum(p64_all, boffset, 3, svp13);
                    svst1_vnum(p64_all, boffset, 4, svp03);
                    svst1_vnum(p64_all, boffset, 5, svp07);
                    svst1_vnum(p64_all, boffset, 6, svp11);
                    svst1_vnum(p64_all, boffset, 7, svp15);
                    svst1_vnum(p64_all, boffset, 8, svp02);
                    svst1_vnum(p64_all, boffset, 9, svp06);
                    svst1_vnum(p64_all, boffset, 10, svp10);
                    svst1_vnum(p64_all, boffset, 11, svp14);
                    svst1_vnum(p64_all, boffset, 12, svp04);
                    svst1_vnum(p64_all, boffset, 13, svp08);
                    svst1_vnum(p64_all, boffset, 14, svp12);
                    svst1_vnum(p64_all, boffset, 15, svp16);

                    aoffset1 += 4;
                    aoffset2 += 4;
                    aoffset3 += 4;
                    aoffset4 += 4;
                    aoffset5 += 4;
                    aoffset6 += 4;
                    aoffset7 += 4;
                    aoffset8 += 4;

                    aoffset9 += 4;
                    aoffset10 += 4;
                    aoffset11 += 4;
                    aoffset12 += 4;
                    aoffset13 += 4;
                    aoffset14 += 4;
                    aoffset15 += 4;
                    aoffset16 += 4;
                    boffset += 64;

                    i--;
                } while (i > 0);
            }

            if (m & 2) {
                svtemp01 = svld1rq(p64_2, aoffset1);
                svtemp02 = svld1rq(p64_2, aoffset2);
                svtemp03 = svld1rq(p64_2, aoffset3);
                svtemp04 = svld1rq(p64_2, aoffset4);
                svtemp05 = svld1rq(p64_2, aoffset5);
                svtemp06 = svld1rq(p64_2, aoffset6);
                svtemp07 = svld1rq(p64_2, aoffset7);
                svtemp08 = svld1rq(p64_2, aoffset8);
                svtemp09 = svld1rq(p64_2, aoffset9);
                svtemp10 = svld1rq(p64_2, aoffset10);
                svtemp11 = svld1rq(p64_2, aoffset11);
                svtemp12 = svld1rq(p64_2, aoffset12);
                svtemp13 = svld1rq(p64_2, aoffset13);
                svtemp14 = svld1rq(p64_2, aoffset14);
                svtemp15 = svld1rq(p64_2, aoffset15);
                svtemp16 = svld1rq(p64_2, aoffset16);

                svins01 = svzip1(svtemp01, svtemp03);
                svins02 = svzip1(svtemp02, svtemp04);

                svp01 = svzip1(svins01, svins02);
                svp02 = svzip2(svins01, svins02);

                svins03 = svzip1(svtemp05, svtemp07);
                svins04 = svzip1(svtemp06, svtemp08);

                svp03 = svzip1(svins03, svins04);
                svp04 = svzip2(svins03, svins04);

                svins05 = svzip1(svtemp09, svtemp11);
                svins06 = svzip1(svtemp10, svtemp12);

                svp05 = svzip1(svins05, svins06);
                svp06 = svzip2(svins05, svins06);

                svins07 = svzip1(svtemp13, svtemp15);
                svins08 = svzip1(svtemp14, svtemp16);

                svp07 = svzip1(svins07, svins08);
                svp08 = svzip2(svins07, svins08);

                svst1(p64_all, boffset, svp01);
                svst1_vnum(p64_all, boffset, 1, svp03);
                svst1_vnum(p64_all, boffset, 2, svp05);
                svst1_vnum(p64_all, boffset, 3, svp07);
                svst1_vnum(p64_all, boffset, 4, svp02);
                svst1_vnum(p64_all, boffset, 5, svp04);
                svst1_vnum(p64_all, boffset, 6, svp06);
                svst1_vnum(p64_all, boffset, 7, svp08);

                aoffset1 += 2;
                aoffset2 += 2;
                aoffset3 += 2;
                aoffset4 += 2;
                aoffset5 += 2;
                aoffset6 += 2;
                aoffset7 += 2;
                aoffset8 += 2;

                aoffset9 += 2;
                aoffset10 += 2;
                aoffset11 += 2;
                aoffset12 += 2;
                aoffset13 += 2;
                aoffset14 += 2;
                aoffset15 += 2;
                aoffset16 += 2;
                boffset += 32;
            }

            if (m & 1) {
                ctemp01 = *(aoffset1 + 0);
                ctemp03 = *(aoffset2 + 0);
                ctemp05 = *(aoffset3 + 0);
                ctemp07 = *(aoffset4 + 0);
                ctemp09 = *(aoffset5 + 0);
                ctemp11 = *(aoffset6 + 0);
                ctemp13 = *(aoffset7 + 0);
                ctemp15 = *(aoffset8 + 0);

                ctemp17 = *(aoffset9 + 0);
                ctemp19 = *(aoffset10 + 0);
                ctemp21 = *(aoffset11 + 0);
                ctemp23 = *(aoffset12 + 0);
                ctemp25 = *(aoffset13 + 0);
                ctemp27 = *(aoffset14 + 0);
                ctemp29 = *(aoffset15 + 0);
                ctemp31 = *(aoffset16 + 0);

                *(boffset + 0) = ctemp01;
                *(boffset + 1) = ctemp03;
                *(boffset + 2) = ctemp05;
                *(boffset + 3) = ctemp07;
                *(boffset + 4) = ctemp09;
                *(boffset + 5) = ctemp11;
                *(boffset + 6) = ctemp13;
                *(boffset + 7) = ctemp15;

                *(boffset + 8) = ctemp17;
                *(boffset + 9) = ctemp19;
                *(boffset + 10) = ctemp21;
                *(boffset + 11) = ctemp23;
                *(boffset + 12) = ctemp25;
                *(boffset + 13) = ctemp27;
                *(boffset + 14) = ctemp29;
                *(boffset + 15) = ctemp31;

                boffset += 16;
            }
            j--;
        } while (j > 0);
    } /* end of if(j > 0) */

    if (n & 8) {
        aoffset1 = aoffset;
        aoffset2 = aoffset1 + lda;
        aoffset3 = aoffset2 + lda;
        aoffset4 = aoffset3 + lda;
        aoffset5 = aoffset4 + lda;
        aoffset6 = aoffset5 + lda;
        aoffset7 = aoffset6 + lda;
        aoffset8 = aoffset7 + lda;
        aoffset += 8 * lda;

        i = (m >> 2);
        if (i > 0) {
            do {
                svtemp01 = svld1(p64_all, aoffset1);
                svtemp02 = svld1(p64_all, aoffset2);
                svtemp03 = svld1(p64_all, aoffset3);
                svtemp04 = svld1(p64_all, aoffset4);
                svtemp05 = svld1(p64_all, aoffset5);
                svtemp06 = svld1(p64_all, aoffset6);
                svtemp07 = svld1(p64_all, aoffset7);
                svtemp08 = svld1(p64_all, aoffset8);

                svins01 = svuzp1(svtemp01, svtemp02);
                svins02 = svuzp2(svtemp01, svtemp02);
                svins03 = svuzp1(svtemp03, svtemp04);
                svins04 = svuzp2(svtemp03, svtemp04);

                svp01 = svuzp1(svins01, svins03);
                svp02 = svuzp2(svins01, svins03);
                svp03 = svuzp1(svins02, svins04);
                svp04 = svuzp2(svins02, svins04);

                svins05 = svuzp1(svtemp05, svtemp06);
                svins06 = svuzp2(svtemp05, svtemp06);
                svins07 = svuzp1(svtemp07, svtemp08);
                svins08 = svuzp2(svtemp07, svtemp08);

                svp05 = svuzp1(svins05, svins07);
                svp06 = svuzp2(svins05, svins07);
                svp07 = svuzp1(svins06, svins08);
                svp08 = svuzp2(svins06, svins08);

                svst1(p64_all, boffset, svp01);
                svst1_vnum(p64_all, boffset, 1, svp05);
                svst1_vnum(p64_all, boffset, 2, svp03);
                svst1_vnum(p64_all, boffset, 3, svp07);
                svst1_vnum(p64_all, boffset, 4, svp02);
                svst1_vnum(p64_all, boffset, 5, svp06);
                svst1_vnum(p64_all, boffset, 6, svp04);
                svst1_vnum(p64_all, boffset, 7, svp08);

                aoffset1 += 4;
                aoffset2 += 4;
                aoffset3 += 4;
                aoffset4 += 4;
                aoffset5 += 4;
                aoffset6 += 4;
                aoffset7 += 4;
                aoffset8 += 4;

                boffset += 32;
                i--;
            } while (i > 0);
        }

        if (m & 2) {
            svtemp01 = svld1rq(p64_2, aoffset1);
            svtemp02 = svld1rq(p64_2, aoffset2);
            svtemp03 = svld1rq(p64_2, aoffset3);
            svtemp04 = svld1rq(p64_2, aoffset4);
            svtemp05 = svld1rq(p64_2, aoffset5);
            svtemp06 = svld1rq(p64_2, aoffset6);
            svtemp07 = svld1rq(p64_2, aoffset7);
            svtemp08 = svld1rq(p64_2, aoffset8);

            svins01 = svzip1(svtemp01, svtemp03);
            svins02 = svzip1(svtemp02, svtemp04);

            svp01 = svzip1(svins01, svins02);
            svp02 = svzip2(svins01, svins02);

            svins03 = svzip1(svtemp05, svtemp07);
            svins04 = svzip1(svtemp06, svtemp08);

            svp03 = svzip1(svins03, svins04);
            svp04 = svzip2(svins03, svins04);

            svst1(p64_all, boffset, svp01);
            svst1_vnum(p64_all, boffset, 1, svp03);
            svst1_vnum(p64_all, boffset, 2, svp02);
            svst1_vnum(p64_all, boffset, 3, svp04);

            aoffset1 += 2;
            aoffset2 += 2;
            aoffset3 += 2;
            aoffset4 += 2;
            aoffset5 += 2;
            aoffset6 += 2;
            aoffset7 += 2;
            aoffset8 += 2;

            boffset += 16;
        }

        if (m & 1) {
            ctemp01 = *(aoffset1 + 0);
            ctemp03 = *(aoffset2 + 0);
            ctemp05 = *(aoffset3 + 0);
            ctemp07 = *(aoffset4 + 0);
            ctemp09 = *(aoffset5 + 0);
            ctemp11 = *(aoffset6 + 0);
            ctemp13 = *(aoffset7 + 0);
            ctemp15 = *(aoffset8 + 0);

            *(boffset + 0) = ctemp01;
            *(boffset + 1) = ctemp03;
            *(boffset + 2) = ctemp05;
            *(boffset + 3) = ctemp07;
            *(boffset + 4) = ctemp09;
            *(boffset + 5) = ctemp11;
            *(boffset + 6) = ctemp13;
            *(boffset + 7) = ctemp15;

            boffset += 8;
        }
    }

    if (n & 4) {
        aoffset1 = aoffset;
        aoffset2 = aoffset1 + lda;
        aoffset3 = aoffset2 + lda;
        aoffset4 = aoffset3 + lda;
        aoffset += 4 * lda;

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

                svst1(p64_all, boffset, svp01);
                svst1_vnum(p64_all, boffset, 1, svp03);
                svst1_vnum(p64_all, boffset, 2, svp02);
                svst1_vnum(p64_all, boffset, 3, svp04);

                aoffset1 += 4;
                aoffset2 += 4;
                aoffset3 += 4;
                aoffset4 += 4;

                boffset += 16;

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

            svst1(p64_all, boffset, svp01);
            svst1_vnum(p64_all, boffset, 1, svp02);

            aoffset1 += 2;
            aoffset2 += 2;
            aoffset3 += 2;
            aoffset4 += 2;

            boffset += 8;
        }

        if (m & 1) {
            ctemp01 = *(aoffset1 + 0);
            ctemp03 = *(aoffset2 + 0);
            ctemp05 = *(aoffset3 + 0);
            ctemp07 = *(aoffset4 + 0);

            *(boffset + 0) = ctemp01;
            *(boffset + 1) = ctemp03;
            *(boffset + 2) = ctemp05;
            *(boffset + 3) = ctemp07;
            boffset += 4;
        }
    }

    if (n & 2) {

        aoffset1 = aoffset;
        aoffset2 = aoffset1 + lda;
        aoffset += 2 * lda;

        i = (m >> 2);
        if (i > 0) {
            do {
                svtemp01 = svld1(p64_all, aoffset1);
                svtemp02 = svld1(p64_all, aoffset2);

                svins01 = svzip1(svtemp01, svtemp02);
                svins02 = svzip2(svtemp01, svtemp02);

                svst1(p64_all, boffset, svins01);
                svst1_vnum(p64_all, boffset, 1, svins02);

                aoffset1 += 4;
                aoffset2 += 4;
                boffset += 8;
                i--;
            } while (i > 0);
        }

        if (m & 2) {
            svtemp01 = svld1rq(p64_2, aoffset1);
            svtemp02 = svld1rq(p64_2, aoffset2);

            svins01 = svzip1(svtemp01, svtemp02);

            svst1(p64_all, boffset, svins01);

            aoffset1 += 2;
            aoffset2 += 2;

            boffset += 4;
        }

        if (m & 1) {

            ctemp01 = *(aoffset1 + 0);
            ctemp03 = *(aoffset2 + 0);

            *(boffset + 0) = ctemp01;
            *(boffset + 1) = ctemp03;
            boffset += 2;
        }
    }

    if (n & 1) {
        aoffset1 = aoffset;

        i = (m >> 2);
        if (i > 0) {
            do {
                svtemp01 = svld1(p64_all, aoffset1);
                svst1(p64_all, boffset, svtemp01);

                aoffset1 += 4;
                boffset += 4;
                i--;
            } while (i > 0);
        }

        if (m & 2) {
            ctemp01 = *(aoffset1 + 0);
            ctemp02 = *(aoffset1 + 1);

            *(boffset + 0) = ctemp01;
            *(boffset + 1) = ctemp02;

            aoffset1 += 2;
            boffset += 2;
        }

        if (m & 1) {
            ctemp01 = *(aoffset1 + 0);

            *(boffset + 0) = ctemp01;
            // boffset   += 1;
        }
    }

    return 0;
}
