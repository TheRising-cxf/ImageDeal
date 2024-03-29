#pragma once
#include "common.h"
DLL_EXPORT_IMPORT int F_GrindskinByGaussFilter(unsigned char* srcData, int width, int height, int channels, int skinMode, int ratio);
DLL_EXPORT_IMPORT int F_Softskin_ChannelMethod(unsigned char* srcData, int width, int height, int channels, unsigned char* lightMap, int ratio);
DLL_EXPORT_IMPORT int F_Softskin_DetailsAddingMethod(unsigned char* srcData, int width, int height, int channels, int ratio, float K);
DLL_EXPORT_IMPORT int F_Softskin_HP(unsigned char* srcData, int width, int height, int channels, int textureRatio, int ratio);
DLL_EXPORT_IMPORT int F_Softskin_A(unsigned char* srcData, int width, int height, int channels, unsigned char* lightMap, int ratio);
DLL_EXPORT_IMPORT int F_Softskin_MixMethod(unsigned char* srcData, int width, int height, int channels, unsigned char* lightMap, int textureRatio, int ratio);