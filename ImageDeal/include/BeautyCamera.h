#pragma once
#include "common.h"
DLL_EXPORT_IMPORT int F_BeautyCamera(unsigned char* srcData, int width, int height, int channels, unsigned char* curveMap, int skinGrindRatio, int skinWhiteRatio, int skinColorRatio, int sharpenRatio, int skinGrindMode, int skinWhiteMode, int skinColorMode, int sharpenMode);