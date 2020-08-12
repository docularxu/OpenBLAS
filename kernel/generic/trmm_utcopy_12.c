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
	
      } else {
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
	
      }

      //i = (m >> 4);
	  i = m / 12;
      if (i > 0) {
	do {
	  if (X < posY) {
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
	    if (X > posY) {

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
	      b[  1] = ZERO;
	      b[  2] = ZERO;
	      b[  3] = ZERO;
	      b[  4] = ZERO;
	      b[  5] = ZERO;
	      b[  6] = ZERO;
	      b[  7] = ZERO;
	      b[  8] = ZERO;
	      b[  9] = ZERO;
	      b[ 10] = ZERO;
	      b[ 11] = ZERO;
	      

	      b[ 12] = *(a02 +  0);
#ifdef UNIT
	      b[ 13] = ONE;
#else
	      b[ 13] = *(a02 +  1);
#endif
	      b[ 14] = ZERO;
	      b[ 15] = ZERO;
	      b[ 16] = ZERO;
	      b[ 17] = ZERO;
	      b[ 18] = ZERO;
	      b[ 19] = ZERO;
	      b[ 20] = ZERO;
	      b[ 21] = ZERO;
	      b[ 22] = ZERO;
	      b[ 23] = ZERO;
	      

	      b[ 24] = *(a03 +  0);
	      b[ 25] = *(a03 +  1);
#ifdef UNIT
	      b[ 26] = ONE;
#else
	      b[ 26] = *(a03 +  2);
#endif
	      b[ 27] = ZERO;
	      b[ 28] = ZERO;
	      b[ 29] = ZERO;
	      b[ 30] = ZERO;
	      b[ 31] = ZERO;
	      b[ 32] = ZERO;
	      b[ 33] = ZERO;
	      b[ 34] = ZERO;
	      b[ 35] = ZERO;
	      
		  
		  b[ 36] = *(a04 +  0);
	      b[ 37] = *(a04 +  1);
	      b[ 38] = *(a04 +  2);
	      
#ifdef UNIT
	      b[ 39] = ONE;
#else
	      b[ 39] = *(a04 +  3);
#endif
	      b[ 40] = ZERO;
	      b[ 41] = ZERO;
	      b[ 42] = ZERO;
	      b[ 43] = ZERO;
	      b[ 44] = ZERO;
	      b[ 45] = ZERO;
	      b[ 46] = ZERO;
	      b[ 47] = ZERO;
	      

	      b[ 48] = *(a05 +  0);
	      b[ 49] = *(a05 +  1);
	      b[ 50] = *(a05 +  2);
	      b[ 51] = *(a05 +  3);
#ifdef UNIT
	      b[ 52] = ONE;
#else
	      b[ 52] = *(a05 +  4);
#endif
	      b[ 53] = ZERO;
	      b[ 54] = ZERO;
	      b[ 55] = ZERO;
	      b[ 56] = ZERO;
	      b[ 57] = ZERO;
	      b[ 58] = ZERO;
	      b[ 59] = ZERO;

	      b[ 60] = *(a06 +  0);
	      b[ 61] = *(a06 +  1);
	      b[ 62] = *(a06 +  2);
	      b[ 63] = *(a06 +  3);
	      b[ 64] = *(a06 +  4);
#ifdef UNIT
	      b[ 65] = ONE;
#else
	      b[ 65] = *(a06 +  5);
#endif
	      b[ 66] = ZERO;
	      b[ 67] = ZERO;
	      b[ 68] = ZERO;
	      b[ 69] = ZERO;
	      b[ 70] = ZERO;
	      b[ 71] = ZERO;

	      b[ 72] = *(a07 +  0);
	      b[ 73] = *(a07 +  1);
	      b[ 74] = *(a07 +  2);
	      b[ 75] = *(a07 +  3);
	      b[ 76] = *(a07 +  4);
	      b[ 77] = *(a07 +  5);
#ifdef UNIT
	      b[ 78] = ONE;
#else
	      b[ 78] = *(a07 +  6);
#endif
	      b[ 79] = ZERO;
	      b[ 80] = ZERO;
	      b[ 81] = ZERO;
	      b[ 82] = ZERO;
	      b[ 83] = ZERO;

	      b[ 84] = *(a08 +  0);
	      b[ 85] = *(a08 +  1);
	      b[ 86] = *(a08 +  2);
	      b[ 87] = *(a08 +  3);
	      b[ 88] = *(a08 +  4);
	      b[ 89] = *(a08 +  5);
	      b[ 90] = *(a08 +  6);
#ifdef UNIT
	      b[ 91] = ONE;
#else
	      b[ 91] = *(a08 +  7);
#endif
	      b[ 92] = ZERO;
	      b[ 93] = ZERO;
	      b[ 94] = ZERO;
	      b[ 95] = ZERO;

	      b[ 96] = *(a09 +  0);
	      b[ 97] = *(a09 +  1);
	      b[ 98] = *(a09 +  2);
	      b[ 99] = *(a09 +  3);
	      b[100] = *(a09 +  4);
	      b[101] = *(a09 +  5);
	      b[102] = *(a09 +  6);
	      b[103] = *(a09 +  7);
#ifdef UNIT
	      b[104] = ONE;
#else
	      b[104] = *(a09 +  8);
#endif
	      b[105] = ZERO;
	      b[106] = ZERO;
	      b[107] = ZERO;


	      b[108] = *(a10 +  0);
	      b[109] = *(a10 +  1);
	      b[110] = *(a10 +  2);
	      b[111] = *(a10 +  3);
	      b[112] = *(a10 +  4);
	      b[113] = *(a10 +  5);
	      b[114] = *(a10 +  6);
	      b[115] = *(a10 +  7);
	      b[116] = *(a10 +  8);
#ifdef UNIT
	      b[117] = ONE;
#else
	      b[117] = *(a10 +  9);
#endif
	      b[118] = ZERO;
	      b[119] = ZERO;

	      b[120] = *(a11 +  0);
	      b[121] = *(a11 +  1);
	      b[122] = *(a11 +  2);
	      b[123] = *(a11 +  3);
	      b[124] = *(a11 +  4);
	      b[125] = *(a11 +  5);
	      b[126] = *(a11 +  6);
	      b[127] = *(a11 +  7);
	      b[128] = *(a11 +  8);
	      b[129] = *(a11 +  9);
#ifdef UNIT
	      b[130] = ONE;
#else
	      b[130] = *(a11 +10);
#endif
	      b[131] = ZERO;

	      b[132] = *(a12 +  0);
	      b[133] = *(a12 +  1);
	      b[134] = *(a12 +  2);
	      b[135] = *(a12 +  3);
	      b[136] = *(a12 +  4);
	      b[137] = *(a12 +  5);
	      b[138] = *(a12 +  6);
	      b[139] = *(a12 +  7);
	      b[140] = *(a12 +  8);
	      b[141] = *(a12 +  9);
	      b[142] = *(a12 + 10);
#ifdef UNIT
	      b[143] = ONE;
#else
	      b[143] = *(a12 + 11);
#endif

	      a01 += 12 * lda;
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
	if (X < posY) {		// 
	  a01 += i;
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
	  
	  b += 12 * i;
	} 
	else if (X > posY) {

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

	  } 
	  else {
#ifdef UNIT
	      b[ 0] = ONE;
#else
	      b[ 0] = *(a01 +  0);
#endif
	      b[ 1] = ZERO;
	      b[ 2] = ZERO;
	      b[ 3] = ZERO;
	      b[ 4] = ZERO;
	      b[ 5] = ZERO;
	      b[ 6] = ZERO;
	      b[ 7] = ZERO;
	      b[ 8] = ZERO;
	      b[ 9] = ZERO;
	      b[10] = ZERO;
	      b[11] = ZERO;
	      

	      if (i >= 2) {
		b[ 0] = *(a02 +  0);
#ifdef UNIT
		b[ 1] = ONE;
#else
		b[ 1] = *(a02 +  1);
#endif
		b[ 2] = ZERO;
		b[ 3] = ZERO;
		b[ 4] = ZERO;
		b[ 5] = ZERO;
		b[ 6] = ZERO;
		b[ 7] = ZERO;
		b[ 8] = ZERO;
		b[ 9] = ZERO;
		b[10] = ZERO;
		b[11] = ZERO;
		
		b += 12;
	      }

	      if (i >= 3) {
		b[ 0] = *(a03 +  0);
		b[ 1] = *(a03 +  1);
#ifdef UNIT
		b[ 2] = ONE;
#else
		b[ 2] = *(a03 +  2);
#endif
		b[ 3] = ZERO;
		b[ 4] = ZERO;
		b[ 5] = ZERO;
		b[ 6] = ZERO;
		b[ 7] = ZERO;
		b[ 8] = ZERO;
		b[ 9] = ZERO;
		b[10] = ZERO;
		b[11] = ZERO;
		
		b += 12;
	      }

	      if (i >= 4) {
		b[ 0] = *(a04 +  0);
		b[ 1] = *(a04 +  1);
		b[ 2] = *(a04 +  2);
#ifdef UNIT
		b[ 3] = ONE;
#else
		b[ 3] = *(a04 +  3);
#endif
		b[ 4] = ZERO;
		b[ 5] = ZERO;
		b[ 6] = ZERO;
		b[ 7] = ZERO;
		b[ 8] = ZERO;
		b[ 9] = ZERO;
		b[10] = ZERO;
		b[11] = ZERO;
		b[12] = ZERO;
		b += 12;
	      }


	      if (i >= 5) {
		b[ 0] = *(a05 +  0);
		b[ 1] = *(a05 +  1);
		b[ 2] = *(a05 +  2);
		b[ 3] = *(a05 +  3);
#ifdef UNIT
		b[ 4] = ONE;
#else
		b[ 4] = *(a05 +  4);
#endif
		b[ 5] = ZERO;
		b[ 6] = ZERO;
		b[ 7] = ZERO;
		b[ 8] = ZERO;
		b[ 9] = ZERO;
		b[10] = ZERO;
		b[11] = ZERO;
	
		b += 12;
	      }

	      if (i >= 6) {
		b[ 0] = *(a06 +  0);
		b[ 1] = *(a06 +  1);
		b[ 2] = *(a06 +  2);
		b[ 3] = *(a06 +  3);
		b[ 4] = *(a06 +  4);
#ifdef UNIT
		b[ 5] = ONE;
#else
		b[ 5] = *(a06 +  5);
#endif
		b[ 6] = ZERO;
		b[ 7] = ZERO;
		b[ 8] = ZERO;
		b[ 9] = ZERO;
		b[10] = ZERO;
		b[11] = ZERO;
		b += 12;
	      }

	      if (i >= 7) {
		b[ 0] = *(a07 +  0);
		b[ 1] = *(a07 +  1);
		b[ 2] = *(a07 +  2);
		b[ 3] = *(a07 +  3);
		b[ 4] = *(a07 +  4);
		b[ 5] = *(a07 +  5);
#ifdef UNIT
		b[ 6] = ONE;
#else
		b[ 6] = *(a07 +  6);
#endif
		b[ 7] = ZERO;
		b[ 8] = ZERO;
		b[ 9] = ZERO;
		b[10] = ZERO;
		b[11] = ZERO;
		
		b += 12;
	      }

	      if (i >= 8) {
		b[ 0] = *(a08 +  0);
		b[ 1] = *(a08 +  1);
		b[ 2] = *(a08 +  2);
		b[ 3] = *(a08 +  3);
		b[ 4] = *(a08 +  4);
		b[ 5] = *(a08 +  5);
		b[ 6] = *(a08 +  6);
#ifdef UNIT
		b[ 7] = ONE;
#else
		b[ 7] = *(a08 +  7);
#endif
		b[ 8] = ZERO;
		b[ 9] = ZERO;
		b[10] = ZERO;
		b[11] = ZERO;
		b[12] = ZERO;
		b += 12;
	      }

	      if (i >= 9) {
		b[ 0] = *(a09 +  0);
		b[ 1] = *(a09 +  1);
		b[ 2] = *(a09 +  2);
		b[ 3] = *(a09 +  3);
		b[ 4] = *(a09 +  4);
		b[ 5] = *(a09 +  5);
		b[ 6] = *(a09 +  6);
		b[ 7] = *(a09 +  7);
#ifdef UNIT
		b[ 8] = ONE;
#else
		b[ 8] = *(a09 +  8);
#endif
		b[ 9] = ZERO;
		b[10] = ZERO;
		b[11] = ZERO;
		
		b += 12;
	      }

	      if (i >= 10) {
		b[ 0] = *(a10 +  0);
		b[ 1] = *(a10 +  1);
		b[ 2] = *(a10 +  2);
		b[ 3] = *(a10 +  3);
		b[ 4] = *(a10 +  4);
		b[ 5] = *(a10 +  5);
		b[ 6] = *(a10 +  6);
		b[ 7] = *(a10 +  7);
		b[ 8] = *(a10 +  8);
#ifdef UNIT
		b[ 9] = ONE;
#else
		b[ 9] = *(a10 +  9);
#endif
		b[10] = ZERO;
		b[11] = ZERO;
		
		b += 12;
	      }

	      if (i >= 11) {
		b[ 0] = *(a11 +  0);
		b[ 1] = *(a11 +  1);
		b[ 2] = *(a11 +  2);
		b[ 3] = *(a11 +  3);
		b[ 4] = *(a11 +  4);
		b[ 5] = *(a11 +  5);
		b[ 6] = *(a11 +  6);
		b[ 7] = *(a11 +  7);
		b[ 8] = *(a11 +  8);
		b[ 9] = *(a11 +  9);
#ifdef UNIT
		b[10] = ONE;
#else
		b[10] = *(a11 + 10);
#endif
		b[11] = ZERO;
		
		b += 12;
	      }
  	  }
      }
      posY += 12;
      js --;
    } while (js > 0);
  } /* End of main loop */

  js = n - (( n / 12) * 12);
  // equal to js = n % 12;  
  mi = m - (( m / 12) * 12);
  if (js & 8){
    X = posX;

    if (posX <= posY) {
      a01 = a + posX + (posY +  0) * lda;
      a02 = a + posX + (posY +  1) * lda;
      a03 = a + posX + (posY +  2) * lda;
      a04 = a + posX + (posY +  3) * lda;
      a05 = a + posX + (posY +  4) * lda;
      a06 = a + posX + (posY +  5) * lda;
      a07 = a + posX + (posY +  6) * lda;
      a08 = a + posX + (posY +  7) * lda;
    } else {
      a01 = a + posY + (posX +  0) * lda;
      a02 = a + posY + (posX +  1) * lda;
      a03 = a + posY + (posX +  2) * lda;
      a04 = a + posY + (posX +  3) * lda;
      a05 = a + posY + (posX +  4) * lda;
      a06 = a + posY + (posX +  5) * lda;
      a07 = a + posY + (posX +  6) * lda;
      a08 = a + posY + (posX +  7) * lda;
    }

	
    i = (m >> 3);
    if (i > 0) {
      do {
	if (X < posY) {
		// Handle cases where 12 is not divisible by 8
		if (posY - X == 4){
		  a01 += 4;
		  a02 += 4;
		  a03 += 4;
		  a04 += 4;
		  a05 += 4;
		  a06 += 4;
		  a07 += 4;
		  a08 += 4;
		  b += 32;
		}
		else{
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
	} 
	else if (X > posY) {

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
	  else {
#ifdef UNIT
	      b[  0] = ONE;
#else
	      b[  0] = *(a01 +  0);
#endif
	      b[  1] = ZERO;
	      b[  2] = ZERO;
	      b[  3] = ZERO;
	      b[  4] = ZERO;
	      b[  5] = ZERO;
	      b[  6] = ZERO;
	      b[  7] = ZERO;

	      b[  8] = *(a02 +  0);
#ifdef UNIT
	      b[  9] = ONE;
#else
	      b[  9] = *(a02 +  1);
#endif
	      b[ 10] = ZERO;
	      b[ 11] = ZERO;
	      b[ 12] = ZERO;
	      b[ 13] = ZERO;
	      b[ 14] = ZERO;
	      b[ 15] = ZERO;

	      b[ 16] = *(a03 +  0);
	      b[ 17] = *(a03 +  1);
#ifdef UNIT
	      b[ 18] = ONE;
#else
	      b[ 18] = *(a03 +  2);
#endif
	      b[ 19] = ZERO;
	      b[ 20] = ZERO;
	      b[ 21] = ZERO;
	      b[ 22] = ZERO;
	      b[ 23] = ZERO;

	      b[ 24] = *(a04 +  0);
	      b[ 25] = *(a04 +  1);
	      b[ 26] = *(a04 +  2);
#ifdef UNIT
	      b[ 27] = ONE;
#else
	      b[ 27] = *(a04 +  3);
#endif
	      b[ 28] = ZERO;
	      b[ 29] = ZERO;
	      b[ 30] = ZERO;
	      b[ 31] = ZERO;

	      b[ 32] = *(a05 +  0);
	      b[ 33] = *(a05 +  1);
	      b[ 34] = *(a05 +  2);
	      b[ 35] = *(a05 +  3);
#ifdef UNIT
	      b[ 36] = ONE;
#else
	      b[ 36] = *(a05 +  4);
#endif
	      b[ 37] = ZERO;
	      b[ 38] = ZERO;
	      b[ 39] = ZERO;

	      b[ 40] = *(a06 +  0);
	      b[ 41] = *(a06 +  1);
	      b[ 42] = *(a06 +  2);
	      b[ 43] = *(a06 +  3);
	      b[ 44] = *(a06 +  4);
#ifdef UNIT
	      b[ 45] = ONE;
#else
	      b[ 45] = *(a06 +  5);
#endif
	      b[ 46] = ZERO;
	      b[ 47] = ZERO;

	      b[ 48] = *(a07 +  0);
	      b[ 49] = *(a07 +  1);
	      b[ 50] = *(a07 +  2);
	      b[ 51] = *(a07 +  3);
	      b[ 52] = *(a07 +  4);
	      b[ 53] = *(a07 +  5);
#ifdef UNIT
	      b[ 54] = ONE;
#else
	      b[ 54] = *(a07 +  6);
#endif
	      b[ 55] = ZERO;

	      b[ 56] = *(a08 +  0);
	      b[ 57] = *(a08 +  1);
	      b[ 58] = *(a08 +  2);
	      b[ 59] = *(a08 +  3);
	      b[ 60] = *(a08 +  4);
	      b[ 61] = *(a08 +  5);
	      b[ 62] = *(a08 +  6);
#ifdef UNIT
	      b[ 63] = ONE;
#else
	      b[ 63] = *(a08 +  7);
#endif

	      a01 += 8 * lda;
	      a02 += 8 * lda;
	      a03 += 8 * lda;
	      a04 += 8 * lda;
	      a05 += 8 * lda;
	      a06 += 8 * lda;
	      a07 += 8 * lda;
	      a08 += 8 * lda;

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
	if (X < posY) {
	  a01 += i;
	  a02 += i;
	  a03 += i;
	  a04 += i;
	  a05 += i;
	  a06 += i;
	  a07 += i;
	  a08 += i;
	  b += 8 * i;
	} else
	  if (X > posY) {
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
	      b += 8;
	    }

	    a02 += i * lda;
	    a03 += i * lda;
	    a04 += i * lda;
	    a05 += i * lda;
	    a06 += i * lda;
	    a07 += i * lda;
	    a08 += i * lda;
	  } else {
#ifdef UNIT
	    b[ 0] = ONE;
#else
	    b[ 0] = *(a01 +  0);
#endif
	    b[ 1] = ZERO;
	    b[ 2] = ZERO;
	    b[ 3] = ZERO;
	    b[ 4] = ZERO;
	    b[ 5] = ZERO;
	    b[ 6] = ZERO;
	    b[ 7] = ZERO;
	    b += 8;

	    if (i >= 2) {
	      b[ 0] = *(a02 +  0);
#ifdef UNIT
	      b[ 1] = ONE;
#else
	      b[ 1] = *(a02 +  1);
#endif
	      b[ 2] = ZERO;
	      b[ 3] = ZERO;
	      b[ 4] = ZERO;
	      b[ 5] = ZERO;
	      b[ 6] = ZERO;
	      b[ 7] = ZERO;
	      b += 8;
	    }

	    if (i >= 3) {
	      b[ 0] = *(a03 +  0);
	      b[ 1] = *(a03 +  1);
#ifdef UNIT
	      b[ 2] = ONE;
#else
	      b[ 2] = *(a03 +  2);
#endif
	      b[ 3] = ZERO;
	      b[ 4] = ZERO;
	      b[ 5] = ZERO;
	      b[ 6] = ZERO;
	      b[ 7] = ZERO;
	      b += 8;
	    }

	    if (i >= 4) {
	      b[ 0] = *(a04 +  0);
	      b[ 1] = *(a04 +  1);
	      b[ 2] = *(a04 +  2);
#ifdef UNIT
	      b[ 3] = ONE;
#else
	      b[ 3] = *(a04 +  3);
#endif
	      b[ 4] = ZERO;
	      b[ 5] = ZERO;
	      b[ 6] = ZERO;
	      b[ 7] = ZERO;
	      b += 8;
	    }

	    if (i >= 5) {
	      b[ 0] = *(a05 +  0);
	      b[ 1] = *(a05 +  1);
	      b[ 2] = *(a05 +  2);
	      b[ 3] = *(a05 +  3);
#ifdef UNIT
	      b[ 4] = ONE;
#else
	      b[ 4] = *(a05 +  4);
#endif
	      b[ 5] = ZERO;
	      b[ 6] = ZERO;
	      b[ 7] = ZERO;
	      b += 8;
	    }

	    if (i >= 6) {
	      b[ 0] = *(a06 +  0);
	      b[ 1] = *(a06 +  1);
	      b[ 2] = *(a06 +  2);
	      b[ 3] = *(a06 +  3);
	      b[ 4] = *(a06 +  4);
#ifdef UNIT
	      b[ 5] = ONE;
#else
	      b[ 5] = *(a06 +  5);
#endif
	      b[ 6] = ZERO;
	      b[ 7] = ZERO;
	      b += 8;
	    }

	    if (i >= 7) {
	      b[ 0] = *(a07 +  0);
	      b[ 1] = *(a07 +  1);
	      b[ 2] = *(a07 +  2);
	      b[ 3] = *(a07 +  3);
	      b[ 4] = *(a07 +  4);
	      b[ 5] = *(a07 +  5);
#ifdef UNIT
	      b[ 6] = ONE;
#else
	      b[ 6] = *(a07 +  6);
#endif
	      b[ 7] = ZERO;
	      b += 8;
	    }
	  }
      }
      posY += 8;
  }

  if (js & 4){
    X = posX;

    if (posX <= posY) {
      a01 = a + posX + (posY +  0) * lda;
      a02 = a + posX + (posY +  1) * lda;
      a03 = a + posX + (posY +  2) * lda;
      a04 = a + posX + (posY +  3) * lda;
    } else {
      a01 = a + posY + (posX +  0) * lda;
      a02 = a + posY + (posX +  1) * lda;
      a03 = a + posY + (posX +  2) * lda;
      a04 = a + posY + (posX +  3) * lda;
    }

    i = (m >> 2);
    if (i > 0) {
      do {
	if (X < posY) {
	  a01 += 4;
	  a02 += 4;
	  a03 += 4;
	  a04 += 4;
	  b += 16;
	} else
	  if (X > posY) {

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
	    b[  1] = ZERO;
	    b[  2] = ZERO;
	    b[  3] = ZERO;

	    b[  4] = *(a02 +  0);
#ifdef UNIT
	    b[  5] = ONE;
#else
	    b[  5] = *(a02 +  1);
#endif
	    b[  6] = ZERO;
	    b[  7] = ZERO;

	    b[  8] = *(a03 +  0);
	    b[  9] = *(a03 +  1);
#ifdef UNIT
	    b[ 10] = ONE;
#else
	    b[ 10] = *(a03 +  2);
#endif
	    b[ 11] = ZERO;

	    b[ 12] = *(a04 +  0);
	    b[ 13] = *(a04 +  1);
	    b[ 14] = *(a04 +  2);
#ifdef UNIT
	    b[ 15] = ONE;
#else
	    b[ 15] = *(a04 +  3);
#endif

	    a01 += 4 * lda;
	    a02 += 4 * lda;
	    a03 += 4 * lda;
	    a04 += 4 * lda;
	    b += 16;
	  }

	  X += 4;
	  i --;
	} while (i > 0);
      }

      i = (m & 3);
      if (i > 0) {
	if (X < posY) {
	  a01 += i;
	  a02 += i;
	  a03 += i;
	  a04 += i;
	  b += 4 * i;
	} else
	  if (X > posY) {
	    for (ii = 0; ii < i; ii++){

	      b[  0] = *(a01 +  0);
	      b[  1] = *(a01 +  1);
	      b[  2] = *(a01 +  2);
	      b[  3] = *(a01 +  3);
	      a01 += lda;
	      b += 4;
	    }
	    a02 += lda;
	    a03 += lda;
	    a04 += lda;
	  } else {

#ifdef UNIT
	    b[ 0] = ONE;
#else
	    b[ 0] = *(a01 +  0);
#endif
	    b[ 1] = ZERO;
	    b[ 2] = ZERO;
	    b[ 3] = ZERO;
	    b += 4;

	    if (i >= 2) {
	      b[ 0] = *(a02 +  0);
#ifdef UNIT
	      b[ 1] = ONE;
#else
	      b[ 1] = *(a02 +  1);
#endif
	      b[ 2] = ZERO;
	      b[ 3] = ZERO;
	      b += 4;
	    }

	    if (i >= 3) {
	      b[ 0] = *(a03 +  0);
	      b[ 1] = *(a03 +  1);
#ifdef UNIT
	      b[ 2] = ONE;
#else
	      b[ 2] = *(a03 +  2);
#endif
	      b[ 3] = ZERO;
	      b += 4;
	    }
	  }
      }
      posY += 4;
  }

  if (js & 2){
    X = posX;

    if (posX <= posY) {
      a01 = a + posX + (posY +  0) * lda;
      a02 = a + posX + (posY +  1) * lda;
    } else {
      a01 = a + posY + (posX +  0) * lda;
      a02 = a + posY + (posX +  1) * lda;
    }

    i = (m >> 1);
    if (i > 0) {
      do {
	if (X < posY) {
	  a01 += 2;
	  a02 += 2;
	  b += 4;
	} else
	  if (X > posY) {
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
	    b[  1] = ZERO;

	    b[  2] = *(a02 +  0);
#ifdef UNIT
	    b[  3] = ONE;
#else
	    b[  3] = *(a02 +  1);
#endif

	    a01 += 2 * lda;
	    a02 += 2 * lda;
	    b += 4;
	  }

	  X += 2;
	  i --;
	} while (i > 0);
      }

    if (m & 1) {
      if (X < posY) {
	a01 ++;
	a02 ++;
	b += 2;
      } else
	if (X > posY) {
	  b[  0] = *(a01 +  0);
	  b[  1] = *(a01 +  1);
	  a01 += lda;
	  b += 2;
	} else {
#ifdef UNIT
	  b[  0] = ONE;
	  b[  1] = *(a01 +  1);
#else
	  b[  0] = *(a01 +  0);
	  b[  1] = *(a01 +  1);
#endif
	  b += 2;
	}
    }
    posY += 2;
  }


  if (js & 1){
    X = posX;

    if (posX <= posY) {
      a01 = a + posX + (posY +  0) * lda;
    } else {
      a01 = a + posY + (posX +  0) * lda;
    }

    i = m;
    if (i > 0) {
      do {
	if (X < posY) {
	  a01 += 1;
	  b ++;
	} else
	  if (X > posY) {
	    b[  0] = *(a01 +  0);
	    a01 += lda;
	    b ++;
	  } else {
#ifdef UNIT
	    b[ 0] = ONE;
#else
	    b[ 0] = *(a01 +  0);
#endif
	    a01 += lda;
	    b ++;
	  }

	  X += 1;
	  i --;
      } while (i > 0);
    }
    posY += 1;
  }

  return 0;
}
