#pragma once
#include "common.h"
#include "FindSkin.h"
#include "ColorConver.h"
#include "ImageFilter.h"
#include <cmath>
#include <string.h>
DLL_EXPORT_IMPORT int F_GrindskinByGaussFilter(unsigned char* srcData, int width, int height, int stride, int skinMode, int ratio);