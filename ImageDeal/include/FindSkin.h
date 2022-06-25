#pragma once
#include "common.h"

//根据RGB对皮肤进行检测
DLL_EXPORT_IMPORT int F_SkinDetectionBGR(unsigned char* src, int width, int high,int channels);

//根据HSV对皮肤进行检测
DLL_EXPORT_IMPORT int F_SkinDetectionHSV(unsigned char* src, int width, int high, int channels);

//根据YCrCb对皮肤进行检测
DLL_EXPORT_IMPORT int F_SkinDetectionYCgCr(unsigned char* src, int width, int high, int channels);

DLL_EXPORT_IMPORT int F_SkinProbability(unsigned char* src, int width, int high, int channels);
