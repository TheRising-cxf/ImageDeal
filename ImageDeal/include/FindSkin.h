#pragma once
#include "common.h"

//����RGB��Ƥ�����м��
DLL_EXPORT_IMPORT int F_SkinDetectionBGR(unsigned char* src, int width, int high,int channels);

//����HSV��Ƥ�����м��
DLL_EXPORT_IMPORT int F_SkinDetectionHSV(unsigned char* src, int width, int high, int channels);

//����YCrCb��Ƥ�����м��
DLL_EXPORT_IMPORT int F_SkinDetectionYCgCr(unsigned char* src, int width, int high, int channels);

DLL_EXPORT_IMPORT int F_SkinProbability(unsigned char* src, int width, int high, int channels);
