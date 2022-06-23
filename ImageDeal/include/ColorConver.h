#pragma once
#include "common.h"
DLL_EXPORT_IMPORT void RGB2HSV(unsigned char R, unsigned char G, unsigned char B, float* h, float* s, float* v);
DLL_EXPORT_IMPORT void HSV2RGB(float h, float s, float v, unsigned char* R, unsigned char* G, unsigned char* B);
DLL_EXPORT_IMPORT void RGBToYCbCr(int R, int G, int B, int* Y, int* Cb, int* Cr);
DLL_EXPORT_IMPORT void YCbCrToRGB(int Y, int Cb, int Cr, int* R, int* G, int* B);