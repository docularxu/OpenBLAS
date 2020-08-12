/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* Copyright (c) 2020, Hisilicon Limited.                            */
/* All rights reserved.                                              */
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

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, BLASLONG posX, BLASLONG posY, FLOAT *b){

  BLASLONG i, js, ii, mi;
  BLASLONG X;

  FLOAT *a01, *a02, *a03 ,*a04, *a05, *a06, *a07, *a08;
  FLOAT *a09, *a10, *a11, *a12;

  printf("Enter %s()\n", __func__);
  printf("m=%ld, n=%ld, lda=%ld, posX=%ld, posY=%ld\n", m, n, lda, posX, posY);
  js = (n / 12);

  if (js > 0){
    do {
      X = posX;

      if (posX <= posY) {
	a01 = a + posY + (posX +  0) * lda;
	a02 = a + posY + (posX +  1) * lda;
	a03 = a + posY + (posX +  2) * lda;
	a04 = a + posY + (posX +  3) * lda;
	a05 = a + posY + (posX +  4) * lda;
	a06 = a + posY + (posX +  5) * lda;
	a07 = a + posY + (posX +  6) * lda;
	a08 = a + posY + (posX +  7) * lda;
	a09 = a + posY + (posX +  8) * lda;
	a10 = a + posY + (posX +  9) * lda;
	a11 = a + posY + (posX + 10) * lda;
	a12 = a + posY + (posX + 11) * lda;
	
      } else {
	a01 = a + posX + (posY +  0) * lda;
	a02 = a + posX + (posY +  1) * lda;
	a03 = a + posX + (posY +  2) * lda;
	a04 = a + posX + (posY +  3) * lda;
	a05 = a + posX + (posY +  4) * lda;
	a06 = a + posX + (posY +  5) * lda;
	a07 = a + posX + (posY +  6) * lda;
	a08 = a + posX + (posY +  7) * lda;
	a09 = a + posX + (posY +  8) * lda;
	a10 = a + posX + (posY +  9) * lda;
	a11 = a + posX + (posY + 10) * lda;
	a12 = a + posX + (posY + 11) * lda;
	
      }

      i = m / 12;
      if (i > 0) {
	do {
	  if (X > posY) {
	    a01 += 12;
	    a02 += 12;
	    a03 += 12;
	    a04 += 12;
	    a05 += 12;
	    a06 += 12;
	    a07 += 12;
	    a08 += 12;
	    a09 += 12;
	    a10 += 12;
	    a11 += 12;
	    a12 += 12;
	    
	    b += 144;
	  } else
	    if (X < posY) {

	      for (ii = 0; ii < 12; ii++){

		b[  0] = *(a01 +  0);
		b[  1] = *(a01 +  1);
		b[  2] = *(a01 +  2);
		b[  3] = *(a01 +  3);
		b[  4] = *(a01 +  4);
		b[  5] = *(a01 +  5);
		b[  6] = *(a01 +  6);
		b[  7] = *(a01 +  7);

		b[  8] = *(a01 +  8);
		b[  9] = *(a01 +  9);
		b[ 10] = *(a01 + 10);
		b[ 11] = *(a01 + 11);
		

		a01 += lda;
		b += 12;
	      }

	      a02 += 12 * lda;
	      a03 += 12 * lda;
	      a04 += 12 * lda;
	      a05 += 12 * lda;
	      a06 += 12 * lda;
	      a07 += 12 * lda;
	      a08 += 12 * lda;
	      a09 += 12 * lda;
	      a10 += 12 * lda;
	      a11 += 12 * lda;
	      a12 += 12 * lda;
	      

	    } else {
#ifdef UNIT
	      b[  0] = ONE;
#else
	      b[  0] = *(a01 +  0);
#endif
	      b[  1] = *(a01 +  1);
	      b[  2] = *(a01 +  2);
	      b[  3] = *(a01 +  3);
	      b[  4] = *(a01 +  4);
	      b[  5] = *(a01 +  5);
	      b[  6] = *(a01 +  6);
	      b[  7] = *(a01 +  7);
	      b[  8] = *(a01 +  8);
	      b[  9] = *(a01 +  9);
	      b[ 10] = *(a01 + 10);
	      b[ 11] = *(a01 + 11);
	      

	      b[ 12] = ZERO;
#ifdef UNIT
	      b[ 13] = ONE;
#else
	      b[ 13] = *(a02 +  1);
#endif
		  b[ 14] = *(a02 +  2);
	      b[ 15] = *(a02 +  3);
	      b[ 16] = *(a02 +  4);
	      b[ 17] = *(a02 +  5);
	      b[ 18] = *(a02 +  6);
	      b[ 19] = *(a02 +  7);
	      b[ 20] = *(a02 +  8);
	      b[ 21] = *(a02 +  9);
	      b[ 22] = *(a02 +  10);
	      b[ 23] = *(a02 +  11);
	      

	      b[ 24] = ZERO;
	      b[ 25] = ZERO;
#ifdef UNIT
	      b[ 26] = ONE;
#else
	      b[ 26] = *(a03 +  2);
#endif
	      b[ 27] = *(a03 +  3);
	      b[ 28] = *(a03 +  4);
	      b[ 29] = *(a03 +  5);
	      b[ 30] = *(a03 +  6);
	      b[ 31] = *(a03 +  7);
	      b[ 32] = *(a03 +  8);
	      b[ 33] = *(a03 +  9);
	      b[ 34] = *(a03 +  10);
	      b[ 35] = *(a03 +  11);
	      

	      b[ 36] = ZERO;
	      b[ 37] = ZERO;
	      b[ 38] = ZERO;
#ifdef UNIT
	      b[ 39] = ONE;
#else
	      b[ 39] = *(a04 +  3);
#endif
	      b[ 40] = *(a04 +  4);
	      b[ 41] = *(a04 +  5);
	      b[ 42] = *(a04 +  6);
	      b[ 43] = *(a04 +  7);
	      b[ 44] = *(a04 +  8);
	      b[ 45] = *(a04 +  9);
	      b[ 46] = *(a04 +  10);
	      b[ 47] = *(a04 +  11);
	      
	      b[ 48] = ZERO;
	      b[ 49] = ZERO;
	      b[ 50] = ZERO;
	      b[ 51] = ZERO;
#ifdef UNIT
	      b[ 52] = ONE;
#else
	      b[ 52] = *(a05 +  4);
#endif
	      b[ 53] = *(a05 +  5);
	      b[ 54] = *(a05 +  6);
	      b[ 55] = *(a05 +  7);
	      b[ 56] = *(a05 +  8);
	      b[ 57] = *(a05 +  9);
	      b[ 58] = *(a05 +  10);
	      b[ 59] = *(a05 +  11);
	      
	      b[ 60] = ZERO;
	      b[ 61] = ZERO;
	      b[ 62] = ZERO;
	      b[ 63] = ZERO;
	      b[ 64] = ZERO;
#ifdef UNIT
	      b[ 65] = ONE;
#else
	      b[ 65] = *(a06 +  5);
#endif
	      b[ 66] = *(a06 +  6);
	      b[ 67] = *(a06 +  7);
	      b[ 68] = *(a06 +  8);
	      b[ 69] = *(a06 +  9);
	      b[ 70] = *(a06 +  10);
	      b[ 71] = *(a06 +  11);
	      
	      b[ 72] = ZERO;
	      b[ 73] = ZERO;
	      b[ 74] = ZERO;
	      b[ 75] = ZERO;
	      b[ 76] = ZERO;
	      b[ 77] = ZERO;

#ifdef UNIT
	      b[ 78] = ONE;
#else
	      b[ 78] = *(a07 +  6);
#endif
	      b[ 79] = *(a07 +  7);
	      b[ 80] = *(a07 +  8);
	      b[ 81] = *(a07 +  9);
	      b[ 82] = *(a07 +  10);
	      b[ 83] = *(a07 +  11);
	   
	      b[ 84] = ZERO;
	      b[ 85] = ZERO;
	      b[ 86] = ZERO;
	      b[ 87] = ZERO;
	      b[ 88] = ZERO;
	      b[ 89] = ZERO;
	      b[ 90] = ZERO;
#ifdef UNIT
	      b[ 91] = ONE;
#else
	      b[ 91] = *(a08 +  7);
#endif
	      b[ 92] = *(a08 +  8);
	      b[ 93] = *(a08 +  9);
	      b[ 94] = *(a08 +  10);
	      b[ 95] = *(a08 +  11);
	      
	      b[ 96] = ZERO;
	      b[ 97] = ZERO;
	      b[ 98] = ZERO;
	      b[ 99] = ZERO;
	      b[100] = ZERO;
	      b[101] = ZERO;
	      b[102] = ZERO;
	      b[103] = ZERO;
#ifdef UNIT
	      b[104] = ONE;
#else
	      b[104] = *(a09 +  8);
#endif
	      b[105] = *(a09 +  9);
	      b[106] = *(a09 + 10);
	      b[107] = *(a09 + 11);
	      
	      b[108] = ZERO;
	      b[109] = ZERO;
	      b[110] = ZERO;
	      b[111] = ZERO;
	      b[112] = ZERO;
	      b[113] = ZERO;
	      b[114] = ZERO;
	      b[115] = ZERO;
	      b[116] = ZERO;
#ifdef UNIT
	      b[117] = ONE;
#else
	      b[117] = *(a10 +  9);
#endif
	      b[118] = *(a10 + 10);
	      b[119] = *(a10 + 11);
	      
	      b[120] = ZERO;
	      b[121] = ZERO;
	      b[122] = ZERO;
	      b[123] = ZERO;
	      b[124] = ZERO;
	      b[125] = ZERO;
	      b[126] = ZERO;
	      b[127] = ZERO;
	      b[128] = ZERO;
	      b[129] = ZERO;
#ifdef UNIT
	      b[130] = ONE;
#else
	      b[130] = *(a11 + 10);
#endif
	      b[131] = *(a11 + 11);
	      
	      b[132] = ZERO;
	      b[133] = ZERO;
	      b[134] = ZERO;
	      b[135] = ZERO;
	      b[136] = ZERO;
	      b[137] = ZERO;
	      b[138] = ZERO;
	      b[139] = ZERO;
	      b[140] = ZERO;
	      b[141] = ZERO;
	      b[142] = ZERO;
#ifdef UNIT
	      b[143] = ONE;
#else
	      b[143] = *(a12 + 11);
#endif

	      a01 += 12;
	      a02 += 12;
	      a03 += 12;
	      a04 += 12;
	      a05 += 12;
	      a06 += 12;
	      a07 += 12;
	      a08 += 12;
	      a09 += 12;
	      a10 += 12;
	      a11 += 12;
	      a12 += 12;
	     
	      b += 144;
	    }

	  X += 12;
	  i --;
	} while (i > 0);
      }

	  //i = (m & 15);
	  // m & 15 equal to m % 16
	  i = m % 12;
      if (i > 0) {
	if (X > posY) {
	  /* a01 += i;
	  a02 += i;
	  a03 += i;
	  a04 += i;
	  a05 += i;
	  a06 += i;
	  a07 += i;
	  a08 += i;
	  a09 += i;
	  a10 += i;
	  a11 += i;
	  a12 += i;
	  a13 += i;
	  a14 += i;
	  a15 += i;
	  a16 += i; */
	  b += 12 * i;
	} else
	  if (X < posY) {

	    for (ii = 0; ii < i; ii++){

	      b[  0] = *(a01 +  0);
	      b[  1] = *(a01 +  1);
	      b[  2] = *(a01 +  2);
	      b[  3] = *(a01 +  3);
	      b[  4] = *(a01 +  4);
	      b[  5] = *(a01 +  5);
	      b[  6] = *(a01 +  6);
	      b[  7] = *(a01 +  7);

	      b[  8] = *(a01 +  8);
	      b[  9] = *(a01 +  9);
	      b[ 10] = *(a01 + 10);
	      b[ 11] = *(a01 + 11);
	      

	      a01 += lda;
	      a02 += lda;
	      a03 += lda;
	      a04 += lda;
	      a05 += lda;
	      a06 += lda;
	      a07 += lda;
	      a08 += lda;
	      a09 += lda;
	      a10 += lda;
	      a11 += lda;
	      a12 += lda;
	      
	      b += 12;
	    }
	  } else {
#ifdef UNIT
	    b[  0] = ONE;
#else
	    b[  0] = *(a01 +  0);
#endif
	    b[  1] = *(a01 +  1);
	    b[  2] = *(a01 +  2);
	    b[  3] = *(a01 +  3);
	    b[  4] = *(a01 +  4);
	    b[  5] = *(a01 +  5);
	    b[  6] = *(a01 +  6);
	    b[  7] = *(a01 +  7);
	    b[  8] = *(a01 +  8);
	    b[  9] = *(a01 +  9);
	    b[ 10] = *(a01 + 10);
	    b[ 11] = *(a01 + 11);
	    b += 12;

	    if (i >= 2) {
	      b[ 0] = ZERO;
#ifdef UNIT
	      b[ 1] = ONE;
#else
	      b[ 1] = *(a02 +  1);
#endif
	      b[ 2] = *(a02 +  2);
	      b[ 3] = *(a02 +  3);
	      b[ 4] = *(a02 +  4);
	      b[ 5] = *(a02 +  5);
	      b[ 6] = *(a02 +  6);
	      b[ 7] = *(a02 +  7);
	      b[ 8] = *(a02 +  8);
	      b[ 9] = *(a02 +  9);
	      b[10] = *(a02 + 10);
	      b[11] = *(a02 + 11);
	      
	      b += 12;
	    }

	    if (i >= 3) {
	      b[ 0] = ZERO;
	      b[ 1] = ZERO;
#ifdef UNIT
	      b[ 2] = ONE;
#else
	      b[ 2] = *(a03 +  2);
#endif
	      b[ 3] = *(a03 +  3);
	      b[ 4] = *(a03 +  4);
	      b[ 5] = *(a03 +  5);
	      b[ 6] = *(a03 +  6);
	      b[ 7] = *(a03 +  7);
	      b[ 8] = *(a03 +  8);
	      b[ 9] = *(a03 +  9);
	      b[10] = *(a03 + 10);
	      b[11] = *(a03 + 11);
	      
	      b += 12;
	    }

	    if (i >= 4) {
	      b[ 0] = ZERO;
	      b[ 1] = ZERO;
	      b[ 2] = ZERO;
#ifdef UNIT
	      b[ 3] = ONE;
#else
	      b[ 3] = *(a04 +  3);
#endif
	      b[ 4] = *(a04 +  4);
	      b[ 5] = *(a04 +  5);
	      b[ 6] = *(a04 +  6);
	      b[ 7] = *(a04 +  7);
	      b[ 8] = *(a04 +  8);
	      b[ 9] = *(a04 +  9);
	      b[10] = *(a04 + 10);
	      b[11] = *(a04 + 11);
	      
	      b += 12;
	    }

	    if (i >= 5) {
	      b[ 0] = ZERO;
	      b[ 1] = ZERO;
	      b[ 2] = ZERO;
	      b[ 3] = ZERO;
#ifdef UNIT
	      b[ 4] = ONE;
#else
	      b[ 4] = *(a05 +  4);
#endif
	      b[ 5] = *(a05 +  5);
	      b[ 6] = *(a05 +  6);
	      b[ 7] = *(a05 +  7);
	      b[ 8] = *(a05 +  8);
	      b[ 9] = *(a05 +  9);
	      b[10] = *(a05 + 10);
	      b[11] = *(a05 + 11);
	      
	      b += 12;
	    }

	    if (i >= 6) {
	      b[ 0] = ZERO;
	      b[ 1] = ZERO;
	      b[ 2] = ZERO;
	      b[ 3] = ZERO;
	      b[ 4] = ZERO;
#ifdef UNIT
	      b[ 5] = ONE;
#else
	      b[ 5] = *(a06 +  5);
#endif
	      b[ 6] = *(a06 +  6);
	      b[ 7] = *(a06 +  7);
	      b[ 8] = *(a06 +  8);
	      b[ 9] = *(a06 +  9);
	      b[10] = *(a06 + 10);
	      b[11] = *(a06 + 11);
	      
	      b += 12;
	    }

	    if (i >= 7) {
	      b[ 0] = ZERO;
	      b[ 1] = ZERO;
	      b[ 2] = ZERO;
	      b[ 3] = ZERO;
	      b[ 4] = ZERO;
	      b[ 5] = ZERO;
#ifdef UNIT
	      b[ 6] = ONE;
#else
	      b[ 6] = *(a07 +  6);
#endif
	      b[ 7] = *(a07 +  7);
	      b[ 8] = *(a07 +  8);
	      b[ 9] = *(a07 +  9);
	      b[10] = *(a07 + 10);
	      b[11] = *(a07 + 11);
	      
	      b += 12;
	    }

	    if (i >= 8) {
	      b[ 0] = ZERO;
	      b[ 1] = ZERO;
	      b[ 2] = ZERO;
	      b[ 3] = ZERO;
	      b[ 4] = ZERO;
	      b[ 5] = ZERO;
	      b[ 6] = ZERO;
#ifdef UNIT
	      b[ 7] = ONE;
#else
	      b[ 7] = *(a08 +  7);
#endif
	      b[ 8] = *(a08 +  8);
	      b[ 9] = *(a08 +  9);
	      b[10] = *(a08 + 10);
	      b[11] = *(a08 + 11);
	      
	      b += 12;
	    }

	    if (i >= 9) {
	      b[ 0] = ZERO;
	      b[ 1] = ZERO;
	      b[ 2] = ZERO;
	      b[ 3] = ZERO;
	      b[ 4] = ZERO;
	      b[ 5] = ZERO;
	      b[ 6] = ZERO;
	      b[ 7] = ZERO;
#ifdef UNIT
	      b[ 8] = ONE;
#else
	      b[ 8] = *(a09 +  8);
#endif
	      b[ 9] = *(a09 +  9);
	      b[10] = *(a09 + 10);
	      b[11] = *(a09 + 11);
	      
	      b += 12;
	    }

	    if (i >= 10) {
	      b[ 0] = ZERO;
	      b[ 1] = ZERO;
	      b[ 2] = ZERO;
	      b[ 3] = ZERO;
	      b[ 4] = ZERO;
	      b[ 5] = ZERO;
	      b[ 6] = ZERO;
	      b[ 7] = ZERO;
	      b[ 8] = ZERO;
#ifdef UNIT
	      b[ 9] = ONE;
#else
	      b[ 9] = *(a10 +  9);
#endif
	      b[10] = *(a10 + 10);
	      b[11] = *(a10 + 11);
	      
	      b += 12;
	    }

	    if (i >= 11) {
	      b[ 0] = ZERO;
	      b[ 1] = ZERO;
	      b[ 2] = ZERO;
	      b[ 3] = ZERO;
	      b[ 4] = ZERO;
	      b[ 5] = ZERO;
	      b[ 6] = ZERO;
	      b[ 7] = ZERO;
	      b[ 8] = ZERO;
	      b[ 9] = ZERO;
#ifdef UNIT
	      b[10] = ONE;
#else
	      b[10] = *(a11 + 10);
#endif
	      b[11] = *(a11 + 11);
	      
	      b += 12;
	    }

	  }
      }

      posY += 12;
      js --;
    } while (js > 0);
  } /* End of main loop */

  js = n - (( n / 12) * 12);
  mi = m - (( m / 12) * 12);
  if (js & 8){
    X = posX;

    if (posX <= posY) {
      a01 = a + posY + (posX +  0) * lda;
      a02 = a + posY + (posX +  1) * lda;
      a03 = a + posY + (posX +  2) * lda;
      a04 = a + posY + (posX +  3) * lda;
      a05 = a + posY + (posX +  4) * lda;
      a06 = a + posY + (posX +  5) * lda;
      a07 = a + posY + (posX +  6) * lda;
      a08 = a + posY + (posX +  7) * lda;
    } 
	else {
      a01 = a + posX + (posY +  0) * lda;
      a02 = a + posX + (posY +  1) * lda;
      a03 = a + posX + (posY +  2) * lda;
      a04 = a + posX + (posY +  3) * lda;
      a05 = a + posX + (posY +  4) * lda;
      a06 = a + posX + (posY +  5) * lda;
      a07 = a + posX + (posY +  6) * lda;
      a08 = a + posX + (posY +  7) * lda;
    }

    i = (m >> 3);
    if (i > 0) {
      do {
	if (X > posY) {
	  a01 += 8;
	  a02 += 8;
	  a03 += 8;
	  a04 += 8;
	  a05 += 8;
	  a06 += 8;
	  a07 += 8;
	  a08 += 8;
	  b += 64;
	} 
	else if (X < posY) {
		// Handle cases where 12 is not divisible by 8
		if (posY - X == 4){
			for (ii = 0; ii < 4; ii++){
			  b[  0] = *(a01 +  0);
			  b[  1] = *(a01 +  1);
			  b[  2] = *(a01 +  2);
			  b[  3] = *(a01 +  3);
			  b[  4] = *(a01 +  4);
			  b[  5] = *(a01 +  5);
			  b[  6] = *(a01 +  6);
			  b[  7] = *(a01 +  7);
			  a01 += lda;
			  b += 8;
			}

			a02 += 4 * lda;
			a03 += 4 * lda;
			a04 += 4 * lda;
			a05 += 4 * lda;
			a06 += 4 * lda;
			a07 += 4 * lda;
			a08 += 4 * lda;
		}
		else{
			for (ii = 0; ii < 8; ii++){
			  b[  0] = *(a01 +  0);
			  b[  1] = *(a01 +  1);
			  b[  2] = *(a01 +  2);
			  b[  3] = *(a01 +  3);
			  b[  4] = *(a01 +  4);
			  b[  5] = *(a01 +  5);
			  b[  6] = *(a01 +  6);
			  b[  7] = *(a01 +  7);
			  a01 += lda;
			  b += 8;
			}

			a02 += 8 * lda;
			a03 += 8 * lda;
			a04 += 8 * lda;
			a05 += 8 * lda;
			a06 += 8 * lda;
			a07 += 8 * lda;
			a08 += 8 * lda;
		}

	    
	  } 
	  else {

#ifdef UNIT
	    b[  0] = ONE;
#else
	    b[  0] = *(a01 +  0);
#endif
	    b[  1] = *(a01 +  1);
	    b[  2] = *(a01 +  2);
	    b[  3] = *(a01 +  3);
	    b[  4] = *(a01 +  4);
	    b[  5] = *(a01 +  5);
	    b[  6] = *(a01 +  6);
	    b[  7] = *(a01 +  7);

	    b[  8] = ZERO;
#ifdef UNIT
	    b[ 9] = ONE;
#else
	    b[ 9] = *(a02 +  1);
#endif
	    b[ 10] = *(a02 +  2);
	    b[ 11] = *(a02 +  3);
	    b[ 12] = *(a02 +  4);
	    b[ 13] = *(a02 +  5);
	    b[ 14] = *(a02 +  6);
	    b[ 15] = *(a02 +  7);

	    b[ 16] = ZERO;
	    b[ 17] = ZERO;
#ifdef UNIT
	    b[ 18] = ONE;
#else
	    b[ 18] = *(a03 +  2);
#endif
	    b[ 19] = *(a03 +  3);
	    b[ 20] = *(a03 +  4);
	    b[ 21] = *(a03 +  5);
	    b[ 22] = *(a03 +  6);
	    b[ 23] = *(a03 +  7);

	    b[ 24] = ZERO;
	    b[ 25] = ZERO;
	    b[ 26] = ZERO;
#ifdef UNIT
	    b[ 27] = ONE;
#else
	    b[ 27] = *(a04 +  3);
#endif
	    b[ 28] = *(a04 +  4);
	    b[ 29] = *(a04 +  5);
	    b[ 30] = *(a04 +  6);
	    b[ 31] = *(a04 +  7);

	    b[ 32] = ZERO;
	    b[ 33] = ZERO;
	    b[ 34] = ZERO;
	    b[ 35] = ZERO;
#ifdef UNIT
	    b[ 36] = ONE;
#else
	    b[ 36] = *(a05 +  4);
#endif
	    b[ 37] = *(a05 +  5);
	    b[ 38] = *(a05 +  6);
	    b[ 39] = *(a05 +  7);

	    b[ 40] = ZERO;
	    b[ 41] = ZERO;
	    b[ 42] = ZERO;
	    b[ 43] = ZERO;
	    b[ 44] = ZERO;
#ifdef UNIT
	    b[ 45] = ONE;
#else
	    b[ 45] = *(a06 +  5);
#endif
	    b[ 46] = *(a06 +  6);
	    b[ 47] = *(a06 +  7);

	    b[ 48] = ZERO;
	    b[ 49] = ZERO;
	    b[ 50] = ZERO;
	    b[ 51] = ZERO;
	    b[ 52] = ZERO;
	    b[ 53] = ZERO;
#ifdef UNIT
	    b[ 54] = ONE;
#else
	    b[ 54] = *(a07 +  6);
#endif
	    b[ 55] = *(a07 +  7);

	    b[ 56] = ZERO;
	    b[ 57] = ZERO;
	    b[ 58] = ZERO;
	    b[ 59] = ZERO;
	    b[ 60] = ZERO;
	    b[ 61] = ZERO;
	    b[ 62] = ZERO;
#ifdef UNIT
	    b[ 63] = ONE;
#else
	    b[ 63] = *(a08 +  7);
#endif

	      a01 += 8;
	      a02 += 8;
	      a03 += 8;
	      a04 += 8;
	      a05 += 8;
	      a06 += 8;
	      a07 += 8;
	      a08 += 8;
	      b += 64;

	  }

	  if (posY - X == 4){
		  X += 4;
	  }
	  else{
	    X += 8;
		i --;
	  }
	  
	} while (i > 0);
      }

    i = (mi & 7);
    if (i > 0) {
      if (X > posY) {
	  /* a01 += i;
	  a02 += i;
	  a03 += i;
	  a04 += i;
	  a05 += i;
	  a06 += i;
	  a07 += i;
	  a08 += i; */
	  b += 8 * i;
      } else
	if (X < posY) {

	  for (ii = 0; ii < i; ii++){

	    b[  0] = *(a01 +  0);
	    b[  1] = *(a01 +  1);
	    b[  2] = *(a01 +  2);
	    b[  3] = *(a01 +  3);
	    b[  4] = *(a01 +  4);
	    b[  5] = *(a01 +  5);
	    b[  6] = *(a01 +  6);
	    b[  7] = *(a01 +  7);

	    a01 += lda;
	    a02 += lda;
	    a03 += lda;
	    a04 += lda;
	    a05 += lda;
	    a06 += lda;
	    a07 += lda;
	    a08 += lda;
	    b += 8;
	  }
	} else {
#ifdef UNIT
	    b[  0] = ONE;
#else
	    b[  0] = *(a01 +  0);
#endif
	    b[  1] = *(a01 +  1);
	    b[  2] = *(a01 +  2);
	    b[  3] = *(a01 +  3);
	    b[  4] = *(a01 +  4);
	    b[  5] = *(a01 +  5);
	    b[  6] = *(a01 +  6);
	    b[  7] = *(a01 +  7);
	    b += 8;

	    if (i >= 2) {
	      b[ 0] = ZERO;
#ifdef UNIT
	      b[ 1] = ONE;
#else
	      b[ 1] = *(a02 +  1);
#endif
	      b[ 2] = *(a02 +  2);
	      b[ 3] = *(a02 +  3);
	      b[ 4] = *(a02 +  4);
	      b[ 5] = *(a02 +  5);
	      b[ 6] = *(a02 +  6);
	      b[ 7] = *(a02 +  7);
	      b += 8;
	    }

	    if (i >= 3) {
	      b[ 0] = ZERO;
	      b[ 1] = ZERO;
#ifdef UNIT
	      b[ 2] = ONE;
#else
	      b[ 2] = *(a03 +  2);
#endif
	      b[ 3] = *(a03 +  3);
	      b[ 4] = *(a03 +  4);
	      b[ 5] = *(a03 +  5);
	      b[ 6] = *(a03 +  6);
	      b[ 7] = *(a03 +  7);
	      b += 8;
	    }

	    if (i >= 4) {
	      b[ 0] = ZERO;
	      b[ 1] = ZERO;
	      b[ 2] = ZERO;
#ifdef UNIT
	      b[ 3] = ONE;
#else
	      b[ 3] = *(a04 +  3);
#endif
	      b[ 4] = *(a04 +  4);
	      b[ 5] = *(a04 +  5);
	      b[ 6] = *(a04 +  6);
	      b[ 7] = *(a04 +  7);
	      b += 8;
	    }

	    if (i >= 5) {
	      b[ 0] = ZERO;
	      b[ 1] = ZERO;
	      b[ 2] = ZERO;
	      b[ 3] = ZERO;
#ifdef UNIT
	      b[ 4] = ONE;
#else
	      b[ 4] = *(a05 +  4);
#endif
	      b[ 5] = *(a05 +  5);
	      b[ 6] = *(a05 +  6);
	      b[ 7] = *(a05 +  7);
	      b += 8;
	    }

	    if (i >= 6) {
	      b[ 0] = ZERO;
	      b[ 1] = ZERO;
	      b[ 2] = ZERO;
	      b[ 3] = ZERO;
	      b[ 4] = ZERO;
#ifdef UNIT
	      b[ 5] = ONE;
#else
	      b[ 5] = *(a06 +  5);
#endif
	      b[ 6] = *(a06 +  6);
	      b[ 7] = *(a06 +  7);
	      b += 8;
	    }

	    if (i >= 7) {
	      b[ 0] = ZERO;
	      b[ 1] = ZERO;
	      b[ 2] = ZERO;
	      b[ 3] = ZERO;
	      b[ 4] = ZERO;
	      b[ 5] = ZERO;
#ifdef UNIT
	      b[ 6] = ONE;
#else
	      b[ 6] = *(a07 +  6);
#endif
	      b[ 7] = *(a07 +  7);
	      b += 8;
	    }
	}
    }
    posY += 8;
  }

  if (js & 4){
    X = posX;

    if (posX <= posY) {
      a01 = a + posY + (posX +  0) * lda;
      a02 = a + posY + (posX +  1) * lda;
      a03 = a + posY + (posX +  2) * lda;
      a04 = a + posY + (posX +  3) * lda;
    } else {
      a01 = a + posX + (posY +  0) * lda;
      a02 = a + posX + (posY +  1) * lda;
      a03 = a + posX + (posY +  2) * lda;
      a04 = a + posX + (posY +  3) * lda;
    }

    i = (m >> 2);
    if (i > 0) {
      do {
	if (X > posY) {
	  a01 += 4;
	  a02 += 4;
	  a03 += 4;
	  a04 += 4;
	  b += 16;
	} else
	  if (X < posY) {

	    for (ii = 0; ii < 4; ii++){
	      b[  0] = *(a01 +  0);
	      b[  1] = *(a01 +  1);
	      b[  2] = *(a01 +  2);
	      b[  3] = *(a01 +  3);
	      a01 += lda;
	      b += 4;
	    }

	    a02 += 4 * lda;
	    a03 += 4 * lda;
	    a04 += 4 * lda;
	  } else {

#ifdef UNIT
	    b[  0] = ONE;
#else
	    b[  0] = *(a01 +  0);
#endif
	    b[  1] = *(a01 +  1);
	    b[  2] = *(a01 +  2);
	    b[  3] = *(a01 +  3);

	    b[  4] = ZERO;
#ifdef UNIT
	    b[  5] = ONE;
#else
	    b[  5] = *(a02 +  1);
#endif
	    b[  6] = *(a02 +  2);
	    b[  7] = *(a02 +  3);

	    b[  8] = ZERO;
	    b[  9] = ZERO;
#ifdef UNIT
	    b[ 10] = ONE;
#else
	    b[ 10] = *(a03 +  2);
#endif
	    b[ 11] = *(a03 +  3);

	    b[ 12] = ZERO;
	    b[ 13] = ZERO;
	    b[ 14] = ZERO;
#ifdef UNIT
	    b[ 15] = ONE;
#else
	    b[ 15] = *(a04 +  3);
#endif

	      a01 += 4;
	      a02 += 4;
	      a03 += 4;
	      a04 += 4;
	      b += 16;
	  }

	X += 4;
	i --;
      } while (i > 0);
    }

    i = (m & 3);
    if (i > 0) {
      if (X > posY) {
	  /* a01 += i;
	  a02 += i;
	  a03 += i;
	  a04 += i; */
	  b += 4 * i;
      } else
	if (X < posY) {

	  for (ii = 0; ii < i; ii++){

	    b[  0] = *(a01 +  0);
	    b[  1] = *(a01 +  1);
	    b[  2] = *(a01 +  2);
	    b[  3] = *(a01 +  3);

	    a01 += lda;
	    a02 += lda;
	    a03 += lda;
	    a04 += lda;
	    b += 4;
	  }
	} else {

#ifdef UNIT
	    b[  0] = ONE;
#else
	    b[  0] = *(a01 +  0);
#endif
	    b[  1] = *(a01 +  1);
	    b[  2] = *(a01 +  2);
	    b[  3] = *(a01 +  3);
	    b += 4;

	    if (i >= 2) {
	      b[ 0] = ZERO;
#ifdef UNIT
	      b[ 1] = ONE;
#else
	      b[ 1] = *(a02 +  1);
#endif
	      b[ 2] = *(a02 +  2);
	      b[ 3] = *(a02 +  3);
	      b += 4;
	    }

	    if (i >= 3) {
	      b[ 0] = ZERO;
	      b[ 1] = ZERO;
#ifdef UNIT
	      b[ 2] = ONE;
#else
	      b[ 2] = *(a03 +  2);
#endif
	      b[ 3] = *(a03 +  3);
	      b += 4;
	    }
	}
    }
    posY += 4;
  }

  if (js & 2){
    X = posX;

    if (posX <= posY) {
      a01 = a + posY + (posX +  0) * lda;
      a02 = a + posY + (posX +  1) * lda;
    } else {
      a01 = a + posX + (posY +  0) * lda;
      a02 = a + posX + (posY +  1) * lda;
    }

    i = (m >> 1);
    if (i > 0) {
      do {
	if (X > posY) {
	  a01 += 2;
	  a02 += 2;
	  b += 4;
	} else
	  if (X < posY) {
	    b[  0] = *(a01 +  0);
	    b[  1] = *(a01 +  1);
	    b[  2] = *(a02 +  0);
	    b[  3] = *(a02 +  1);
	    a01 += 2 * lda;
	    a02 += 2 * lda;
	    b += 4;
	  } else {
#ifdef UNIT
	    b[  0] = ONE;
#else
	    b[  0] = *(a01 +  0);
#endif
	    b[  1] = *(a01 +  1);

	    b[  2] = ZERO;
#ifdef UNIT
	    b[  3] = ONE;
#else
	    b[  3] = *(a02 +  1);
#endif

	      a01 += 2;
	      a02 += 2;
	      b += 4;
	  }

	X += 2;
	i --;
      } while (i > 0);
    }

    if (m & 1) {
      if (X > posY) {
	a01 ++;
	a02 ++;
	b += 2;
      } else
	if (X < posY) {
	  b[  0] = *(a01 +  0);
	  b[  1] = *(a01 +  1);
	  /* a01 += lda;
	  a02 += lda; */
	  b += 2;
	  }
    } else {
#ifdef UNIT
	    b[  0] = ONE;
#else
	    b[  0] = *(a01 +  0);
#endif
	    b[  1] = *(a01 +  1);
      b += 2;
    }
    posY += 2;
  }

  if (js & 1){
    X = posX;

    if (posX <= posY) {
      a01 = a + posY + (posX +  0) * lda;
    } else {
      a01 = a + posX + (posY +  0) * lda;
    }

    i = m;
    if (i > 0) {
      do {

	if (X > posY) {
	  b ++;
	  a01 ++;
	} else
	  if (X < posY) {
	    b[  0] = *(a01 + 0);
	    a01 += lda;
	    b ++;
	  } else {
#ifdef UNIT
	    b[  0] = ONE;
#else
	    b[  0] = *(a01 +  0);
#endif
	    a01 ++;
	    b ++;
	  }
	X += 1;
	i --;
      } while (i > 0);
    }
    // posY += 1;
  }

  return 0;
}
