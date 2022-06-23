#pragma once
#include "common.h"
#include "ColorConver.h"
#include <cmath>
#include <string.h>
DLL_EXPORT_IMPORT int F_FastGaussFilter(unsigned char* srcData, int width, int height, int stride, int r);
DLL_EXPORT_IMPORT int F_SmartBlur(unsigned char* srcData, int nWidth, int nHeight, int nStride, int radius, int threshold);