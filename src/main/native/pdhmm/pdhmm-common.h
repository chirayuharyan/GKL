/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2023-2024 Intel Corporation
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#ifndef AVX_PDHMM_COMMON
#define AVX_PDHMM_COMMON
#define TRANS_PROB_ARRAY_LENGTH 6
#define MAX_QUAL 254

#define OFF 1

#define ROW_UNROLL 4
#define ALIGN_SIZE 64

#define CAT(X, Y) X##Y
#define CONCAT(X, Y) CAT(X, Y)

#define PDHMM_SUCCESS 0
#define PDHMM_MEMORY_ALLOCATION_FAILED 1
#define PDHMM_INPUT_DATA_ERROR 2
#define PDHMM_FAILURE 3
#define PDHMM_MEMORY_ACCESS_ERROR 4

#include "MathUtils.h"
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <debug.h>
#include <immintrin.h>
#include <algorithm>

extern double g_matchToMatchLog10[(((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1)] __attribute__((aligned(64)));

extern double g_matchToMatchProb[(((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1)] __attribute__((aligned(64)));

extern double g_qualToErrorProbCache[(MAX_QUAL + 1)] __attribute__((aligned(64)));

extern double g_qualToProbLog10Cache[(MAX_QUAL + 1)] __attribute__((aligned(64)));

enum HMMState
{
    // The regular state
    NORMAL,

    // Indicating that we must be copying the array elements to the right
    INSIDE_DEL,

    // Indicating that we must handle the special case for merging events after
    // the del
    AFTER_DEL,
};

enum ProbIndex
{
    matchToMatch,
    indelToMatch,
    matchToInsertion,
    insertionToInsertion,
    matchToDeletion,
    deletionToDeletion,
};

double qualToErrorProb(double qual, int32_t &status);

int32_t init(double *&matchToMatchLog10, double *&matchToMatchProb,
             double *&qualToErrorProbCache,
             double *&qualToProbLog10Cache);

int32_t initializeCache();

#endif
