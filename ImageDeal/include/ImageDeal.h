#pragma once
#ifdef EXPORTS
#define DLL_EXPORT_IMPORT  _declspec(dllexport)
#else
#define DLL_EXPORT_IMPORT  _declspec(dllimport)
#endif
DLL_EXPORT_IMPORT int add(int a, int b);