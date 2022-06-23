#pragma once
#include"common.h"
DLL_EXPORT_IMPORT int F_SkinWhite(unsigned char* srcData, int width, int height, int stride, unsigned char* lutData, int ratio);
DLL_EXPORT_IMPORT int F_SkinWhiteCurve(unsigned char* srcData, int width, int height, int stride, int belta, int ratio);
DLL_EXPORT_IMPORT int F_SkinWhitePS(unsigned char* srcData, int width, int height, int stride, int ratio);