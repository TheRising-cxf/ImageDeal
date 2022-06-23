#include "GrindSkin.h"
int F_GrindskinByGaussFilter(unsigned char* srcData, int width, int height, int stride, int skinMode, int ratio)
{
	int ret = 0;
	int length = width * height * stride;
	unsigned char* smoothData = (unsigned char*)malloc(sizeof(unsigned char) * length);
	memcpy(smoothData, srcData, sizeof(unsigned char) * length);
	unsigned char* skinProbability = (unsigned char*)malloc(sizeof(unsigned char) * length);
	memcpy(skinProbability, srcData, sizeof(unsigned char) * length);
	int smoothRadius = 8;
	int smoothThreshold = 38;
	int maskSmoothRadius = 3;
	ret = F_SmartBlur(smoothData, width, height, stride, smoothRadius, smoothThreshold);
	if (skinMode == 0)
	{
		ret = F_SkinDetectionYCgCr(skinProbability, width, height, stride);
		maskSmoothRadius = 6;
	}
	else
	{
		ret = F_SkinProbability(skinProbability, width, height, stride);
	}
	ret = F_FastGaussFilter(skinProbability, width, height, stride, maskSmoothRadius);
	unsigned char* pSrc = srcData;
	unsigned char* pMask = skinProbability;
	unsigned char* pSmooth = smoothData;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			int mask = (pMask[0] + pMask[1] + pMask[2]) / 3;
			int b = (pSrc[0] * (255 - mask) + pSmooth[0] * mask) / 255;
			int g = (pSrc[1] * (255 - mask) + pSmooth[1] * mask) / 255;
			int r = (pSrc[2] * (255 - mask) + pSmooth[2] * mask) / 255;
			pSrc[0] = CLIP3((b * ratio + pSrc[0] * (100 - ratio)) / 100, 0, 255);
			pSrc[1] = CLIP3((g * ratio + pSrc[1] * (100 - ratio)) / 100, 0, 255);
			pSrc[2] = CLIP3((r * ratio + pSrc[2] * (100 - ratio)) / 100, 0, 255);
			pSrc += stride;
			pSmooth += stride;
			pMask += stride;
		}
	}
	free(skinProbability);
	free(smoothData);
	return ret;
};