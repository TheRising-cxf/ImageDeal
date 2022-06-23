#pragma once
#include "common.h"
DLL_EXPORT_IMPORT int F_FastGaussFilter(unsigned char* srcData, int width, int height, int stride, int r);
DLL_EXPORT_IMPORT int F_SmartBlur(unsigned char* srcData, int nWidth, int nHeight, int nStride, int radius, int threshold);
DLL_EXPORT_IMPORT int F_Filter512(unsigned  char*  srcData, int  width, int  height, int  stride, unsigned  char*Map);
DLL_EXPORT_IMPORT int F_SurfaceBlur(unsigned char* srcData, int width, int height, int stride, int radius, int threshold);