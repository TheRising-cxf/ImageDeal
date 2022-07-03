#pragma once
#include "common.h"
DLL_EXPORT_IMPORT int F_IDW(unsigned char* srcData, int width, int height, int channels, int inputPoints[], int outputPoints[], int pointNum);
DLL_EXPORT_IMPORT int F_FBIM(unsigned char* srcData, int width, int height, int channels, int inputlinePoints[], int outputlinePoints[], int lineNum);
DLL_EXPORT_IMPORT void F_MLSImageWrapping(unsigned char* oriImg, int width, int height, int channels, int srcPoint[], int dragPoint[], int pointNum, double transRatio, int preScale, int gridSize, int method);