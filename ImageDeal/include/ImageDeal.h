#pragma once
#ifdef EXPORTS
#define DLL_EXPORT_IMPORT  _declspec(dllexport)
#else
#define DLL_EXPORT_IMPORT  _declspec(dllimport)
#endif

//根据RGB对皮肤进行检测
DLL_EXPORT_IMPORT int F_SkinDetectionBGR(unsigned char* src,int width,int high);

//根据HSV对皮肤进行检测
DLL_EXPORT_IMPORT int F_SkinDetectionHSV(unsigned char* src, int width, int high);

//根据YCrCb对皮肤进行检测
DLL_EXPORT_IMPORT int F_SkinDetectionYCgCr(unsigned char* src, int width, int high);

DLL_EXPORT_IMPORT int F_SkinProbability(unsigned char* src, int width, int high);
