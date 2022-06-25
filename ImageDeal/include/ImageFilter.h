#pragma once
#include "common.h"
DLL_EXPORT_IMPORT int F_StdBilateralFilter(unsigned char* srcData, int nWidth, int nHeight, int channels, int radius, int deltaD, int deltaR);
DLL_EXPORT_IMPORT int F_FastGaussFilter(unsigned char* srcData, int width, int height, int channels, int r);
DLL_EXPORT_IMPORT int F_SmartBlur(unsigned char* srcData, int nWidth, int nHeight, int channels, int radius, int threshold);
DLL_EXPORT_IMPORT int F_Filter512(unsigned  char*  srcData, int  width, int  height, int  channels, unsigned  char*Map);
DLL_EXPORT_IMPORT int F_SurfaceBlur(unsigned char* srcData, int width, int height, int channels, int radius, int threshold);
DLL_EXPORT_IMPORT int F_GuidedFilter(unsigned char* srcData, int nWidth, int nHeight, int channels, int radius, float delta);
DLL_EXPORT_IMPORT int F_LSNFilter(unsigned char* srcData, int nWidth, int nHeight, int channels, int radius, int delta2);
DLL_EXPORT_IMPORT int F_AnisotropicFilter(unsigned char* srcData, int width, int height, int channels, int iter, float k, float lambda, int offset);
DLL_EXPORT_IMPORT int F_MeanShiftFilter(unsigned char* srcData, int nWidth, int nHeight, int channels, int radius, int threshold, int maxIter);
DLL_EXPORT_IMPORT int F_BeepsFilter(unsigned char* srcData, int nWidth, int nHeight, int channels, float delta, float delta_s);