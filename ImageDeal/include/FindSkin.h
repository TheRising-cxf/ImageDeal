#pragma once
#include "common.h"

//����RGB��Ƥ�����м��
DLL_EXPORT_IMPORT int F_SkinDetectionBGR(unsigned char* src, int width, int high,int stride);

//����HSV��Ƥ�����м��
DLL_EXPORT_IMPORT int F_SkinDetectionHSV(unsigned char* src, int width, int high, int stride);

//����YCrCb��Ƥ�����м��
DLL_EXPORT_IMPORT int F_SkinDetectionYCgCr(unsigned char* src, int width, int high, int stride);

DLL_EXPORT_IMPORT int F_SkinProbability(unsigned char* src, int width, int high, int stride);
