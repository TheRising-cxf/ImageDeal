#pragma once
#ifdef EXPORTS
#define DLL_EXPORT_IMPORT  _declspec(dllexport)
#else
#define DLL_EXPORT_IMPORT  _declspec(dllimport)
#endif

//����RGB��Ƥ�����м��
DLL_EXPORT_IMPORT int F_SkinDetectionBGR(unsigned char* src,int width,int high);

//����HSV��Ƥ�����м��
DLL_EXPORT_IMPORT int F_SkinDetectionHSV(unsigned char* src, int width, int high);

//����YCrCb��Ƥ�����м��
DLL_EXPORT_IMPORT int F_SkinDetectionYCgCr(unsigned char* src, int width, int high);

DLL_EXPORT_IMPORT int F_SkinProbability(unsigned char* src, int width, int high);
