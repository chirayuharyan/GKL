/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2016-2024 Intel Corporation
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
/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class com_intel_gkl_smithwaterman_IntelSmithWaterman */

#ifndef _Included_com_intel_gkl_smithwaterman_IntelSmithWaterman
#define _Included_com_intel_gkl_smithwaterman_IntelSmithWaterman
#ifdef __cplusplus
extern "C" {
#endif

JNIEXPORT void JNICALL Java_com_intel_gkl_smithwaterman_IntelSmithWaterman_initNative
  (JNIEnv * env, jclass obj );
/*
 * Class:     com_intel_gkl_smithwaterman_IntelSmithWaterman
 * Method:    alignNative
 * Signature: ([B[B[BIIIII)I
 */
JNIEXPORT jint JNICALL Java_com_intel_gkl_smithwaterman_IntelSmithWaterman_alignNative
  (JNIEnv *, jclass, jbyteArray, jbyteArray, jbyteArray, jint, jint, jint, jint, jbyte);

/*
 * Class:     com_intel_gkl_smithwaterman_IntelSmithWaterman
 * Method:    doneNative
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_com_intel_gkl_smithwaterman_IntelSmithWaterman_doneNative
  (JNIEnv *, jclass);

#ifdef __cplusplus
}
#endif
#endif
