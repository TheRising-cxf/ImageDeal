#pragma once
#include "common.h"
DLL_EXPORT_IMPORT int F_IDW(unsigned char* srcData, int width, int height, int channels, int inputPoints[], int outputPoints[], int pointNum);
DLL_EXPORT_IMPORT int F_FBIM(unsigned char* srcData, int width, int height, int channels, int inputlinePoints[], int outputlinePoints[], int lineNum);