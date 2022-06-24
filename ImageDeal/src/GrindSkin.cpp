#include "GrindSkin.h"
#include "FindSkin.h"
#include "ColorConver.h"
#include "ImageFilter.h"
#include <cmath>
#include <string.h>
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
int F_Softskin_ChannelMethod(unsigned char* srcData, int width, int height, int stride, unsigned char* lightMap, int ratio)
{
	int ret = 0;
	unsigned char* greenData = (unsigned char*)malloc(sizeof(unsigned char)* width * stride * height);
	unsigned char* gaussData = (unsigned char*)malloc(sizeof(unsigned char) * width * stride * height);
	unsigned char* curveData = (unsigned char*)malloc(sizeof(unsigned char) * width * stride * height);
	unsigned char* skinData = (unsigned char*)malloc(sizeof(unsigned char) * width * stride * height);
	unsigned char* pSrc = srcData;
	unsigned char* pGreen = greenData;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			pGreen[0] = pSrc[0];
			pGreen[1] = pSrc[0];
			pGreen[2] = pSrc[0];
			pSrc += 4;
			pGreen += 4;
		}
	}
	memcpy(gaussData, greenData, sizeof(unsigned char) * width * height * stride);
	memcpy(curveData, srcData, sizeof(unsigned char) * width * height * stride);
	ret = F_Filter512(curveData, width, height, stride, lightMap);
	float hpRadius = 10.0f * width * height / (594 * 677);
	ret = F_FastGaussFilter(gaussData, width, height, stride, hpRadius);
	pSrc = srcData;
	pGreen = greenData;
	unsigned char* pCurve = curveData;
	unsigned char* pGauss = gaussData;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			int t;
			t = CLIP3(pGauss[0] - pGreen[0] + 128, 0, 255);
			t = ModeSuperposition(t, t);
			t = ModeSuperposition(t, t);
			t = t * 220 / 255;
			pGreen[0] = CLIP3((pCurve[0] * t + (255 - t) * pSrc[0]) / 255, 0, 255);
			pGreen[1] = CLIP3((pCurve[1] * t + (255 - t) * pSrc[1]) / 255, 0, 255);
			pGreen[2] = CLIP3((pCurve[2] * t + (255 - t) * pSrc[2]) / 255, 0, 255);
			pGreen += stride;
			pGauss += stride;
			pSrc += stride;
			pCurve += stride;
		}
	}
	memcpy(skinData, greenData, sizeof(unsigned char)* width * height * stride);
	int maskSmoothRadius = 3 * width * height / (594 * 677);
	ret = F_SkinProbability(skinData, width, height, stride);
	ret = F_FastGaussFilter(skinData, width, height, stride, maskSmoothRadius);
	pGauss = skinData;
	pSrc = srcData;
	pGreen = greenData;
	int k = ratio * 128 / 100;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			int mask = (pGauss[0] + pGauss[1] + pGauss[2]) / 3;
			int tb = CLIP3((pSrc[0] * (255 - mask) + pGreen[0] * mask) / 255, 0, 255);
			int tg = CLIP3((pSrc[1] * (255 - mask) + pGreen[1] * mask) / 255, 0, 255);
			int tr = CLIP3((pSrc[2] * (255 - mask) + pGreen[2] * mask) / 255, 0, 255);
			pSrc[0] = CLIP3((pSrc[0] * (128 - k) + tb * k) >> 7, 0, 255);
			pSrc[1] = CLIP3((pSrc[1] * (128 - k) + tg * k) >> 7, 0, 255);
			pSrc[2] = CLIP3((pSrc[2] * (128 - k) + tr * k) >> 7, 0, 255);
			pSrc += stride;
			pGauss += stride;
			pGreen += stride;
		}
	}
	free(gaussData);
	free(greenData);
	free(curveData);
	free(skinData);
	return ret;
};
int F_Softskin_DetailsAddingMethod(unsigned char* srcData, int width, int height, int stride, int ratio, float K)
{
	int ret = 0;
	int len = sizeof(unsigned char) * width * height * stride;
	unsigned char*coarseSmoothData = (unsigned char*)malloc(len);
	unsigned char*fineSmoothData = (unsigned char*)malloc(len);
	memcpy(coarseSmoothData, srcData, len);
	memcpy(fineSmoothData, srcData, len);
	int std_fine = 5;
	int std_coarse = 10;
	F_SmartBlur(fineSmoothData, width, height, stride, std_fine, 30);
	F_SmartBlur(coarseSmoothData, width, height, stride, std_coarse, 30);
	unsigned char* skinPDF = (unsigned char*)malloc(sizeof(unsigned char) * len);
	memcpy(skinPDF, coarseSmoothData, sizeof(unsigned char) * len);
	ret = F_SkinProbability(skinPDF, width, height, stride);
	float maskSmoothRadius = 3;
	ret = F_FastGaussFilter(skinPDF, width, height, stride, maskSmoothRadius);
	unsigned char* pSrc = srcData;
	unsigned char* pCoarse = coarseSmoothData;
	unsigned char* pFine = fineSmoothData;
	unsigned char* pSkin = skinPDF;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			int alpha = *pSkin * ratio / 100;
			int detailsFine = pSrc[0] - pFine[0];
			int detailsCoarse = pFine[0] - pCoarse[0];
			float K0 = alpha / 255.0f;
			pSrc[0] = (unsigned char)CLIP3(pCoarse[0] + (255 - alpha) * detailsCoarse / 255 + (1.0f - K0 * K) * detailsFine, 0, 255);

			detailsFine = pSrc[1] - pFine[1];
			detailsCoarse = pFine[1] - pCoarse[1];
			pSrc[1] = (unsigned char)CLIP3(pCoarse[1] + (255 - alpha) * detailsCoarse / 255 + (1.0f - K0 * K) * detailsFine, 0, 255);

			detailsFine = pSrc[2] - pFine[2];
			detailsCoarse = pFine[2] - pCoarse[2];
			pSrc[2] = (unsigned char)CLIP3(pCoarse[2] + (255 - alpha) * detailsCoarse / 255 + (1.0f - K0 * K) * detailsFine, 0, 255);

			pSrc += stride;
			pCoarse += stride;
			pFine += stride;
			pSkin += stride;
		}
	}
	free(skinPDF);
	free(coarseSmoothData);
	free(fineSmoothData);
	return ret;
};
int F_Softskin_HP(unsigned char* srcData, int width, int height, int stride, int textureRatio, int ratio)
{
	int ret = 0;
	int length = width * height * stride;
	unsigned char* smoothData = (unsigned char*)malloc(sizeof(unsigned char) * length);
	unsigned char* hpData = (unsigned char*)malloc(sizeof(unsigned char) * length);
	memcpy(smoothData, srcData, sizeof(unsigned char) * length);
	unsigned char* skinPDF = (unsigned char*)malloc(sizeof(unsigned char) * length);
	int smoothRadius = 8;
	int smoothThreshold = 38;
	int maskSmoothRadius = 3;
	ret = F_SurfaceBlur(smoothData, width, height, stride, smoothRadius, smoothThreshold);
	memcpy(skinPDF, smoothData, sizeof(unsigned char) * length);
	ret = F_SkinProbability(skinPDF, width, height, stride);
	ret = F_FastGaussFilter(skinPDF, width, height, stride, maskSmoothRadius);
	unsigned char* pSrc = srcData;
	unsigned char* pSkin = skinPDF;
	unsigned char* pSmooth = smoothData;
	unsigned char* pHP = hpData;
	int k = ratio * 128 / 100;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			pHP[0] = CLIP3(pSmooth[0] - pSrc[0] + 128, 0, 255);
			pHP[1] = CLIP3(pSmooth[1] - pSrc[1] + 128, 0, 255);
			pHP[2] = CLIP3(pSmooth[2] - pSrc[2] + 128, 0, 255);

			pHP += stride;
			pSmooth += stride;
			pSrc += stride;
		}
	}
	float hpRadius = 3.5f * textureRatio / 100.0f;
	ret = F_FastGaussFilter(hpData, width, height, stride, hpRadius);
	pSmooth = smoothData;
	pHP = hpData;
	pSrc = srcData;

	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			int hpb = pHP[0];
			int hpg = pHP[1];
			int hpr = pHP[2];
			hpb = ModeLinearLight(pSrc[0], hpb);
			hpg = ModeLinearLight(pSrc[1], hpg);
			hpr = ModeLinearLight(pSrc[2], hpr);
			int mask = (pSkin[0] + pSkin[1] + pSkin[2]) / 3;
			hpb = CLIP3((hpb * mask + pSmooth[0] * (255 - mask)) / 255, 0, 255);
			hpg = CLIP3((hpg * mask + pSmooth[1] * (255 - mask)) / 255, 0, 255);
			hpr = CLIP3((hpr * mask + pSmooth[2] * (255 - mask)) / 255, 0, 255);
			pSrc[0] = CLIP3((hpb * k + pSrc[0] * (128 - k)) >> 7, 0, 255);
			pSrc[1] = CLIP3((hpg * k + pSrc[1] * (128 - k)) >> 7, 0, 255);
			pSrc[2] = CLIP3((hpr * k + pSrc[2] * (128 - k)) >> 7, 0, 255);
			pSrc += stride;
			pHP += stride;
			pSmooth += stride;
			pSkin += stride;
		}
	}
	free(skinPDF);
	free(smoothData);
	free(hpData);
	return ret = 0;
};
int F_Softskin_A(unsigned char* srcData, int width, int height, int stride, unsigned char* lightMap, int ratio)
{
	int ret = 0;
	unsigned char* greenData = (unsigned char*)malloc(sizeof(unsigned char)* width * stride * height);
	unsigned char* gaussData = (unsigned char*)malloc(sizeof(unsigned char) * width * stride * height);
	unsigned char* curveData = (unsigned char*)malloc(sizeof(unsigned char) * width * stride * height);
	unsigned char* skinData = (unsigned char*)malloc(sizeof(unsigned char) * width * stride * height);
	unsigned char* smoothData = (unsigned char*)malloc(sizeof(unsigned char) * width * stride * height);
	unsigned char* pSrc = srcData;
	unsigned char* pGreen = greenData;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			pGreen[0] = pSrc[0];
			pGreen[1] = pSrc[0];
			pGreen[2] = pSrc[0];
			pSrc += stride;
			pGreen += stride;
		}
	}
	memcpy(gaussData, greenData, sizeof(unsigned char) * height * stride);
	memcpy(curveData, srcData, sizeof(unsigned char) * height * stride);
	ret = F_Filter512(curveData, width, height, stride, lightMap);
	float hpRadius = 6.0f * width * height / (594 * 677);
	ret = F_FastGaussFilter(gaussData, width, height, stride, hpRadius);
	pSrc = srcData;
	pGreen = greenData;
	unsigned char* pCurve = curveData;
	unsigned char* pGauss = gaussData;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			int t;
			t = CLIP3(pGauss[0] - pGreen[0] + 128, 0, 255);
			t = ModeSuperposition(t, t);
			t = ModeSuperposition(t, t);
			t = t * 200 / 255;
			pGreen[0] = CLIP3((pCurve[0] * t + (255 - t) * pSrc[0]) / 255, 0, 255);
			pGreen[1] = CLIP3((pCurve[1] * t + (255 - t) * pSrc[1]) / 255, 0, 255);
			pGreen[2] = CLIP3((pCurve[2] * t + (255 - t) * pSrc[2]) / 255, 0, 255);
			pGreen += stride;
			pGauss += stride;
			pSrc += stride;
			pCurve += stride;
		}
	}

	int k = ratio * 128 / 100;
	memcpy(smoothData, greenData, sizeof(unsigned char) * stride * height);
	int smoothRadius = 6;
	int smoothThreshold = 38;
	ret = F_SmartBlur(smoothData, width, height, stride, smoothRadius, smoothThreshold);
	memcpy(skinData, smoothData, sizeof(unsigned char) * height * stride);
	int maskSmoothRadius = 3;
	ret = F_SkinProbability(skinData, width, height, stride);
	ret = F_FastGaussFilter(skinData, width, height, stride, maskSmoothRadius);
	pGauss = skinData;
	pSrc = srcData;
	pGreen = greenData;
	unsigned char* pSmooth = smoothData;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			int mask = (pGauss[0] + pGauss[1] + pGauss[2]) / 3;
			int tb = CLIP3((pSrc[0] * (255 - mask) + pGreen[0] * mask) / 255, 0, 255);
			int tg = CLIP3((pSrc[1] * (255 - mask) + pGreen[1] * mask) / 255, 0, 255);
			int tr = CLIP3((pSrc[2] * (255 - mask) + pGreen[2] * mask) / 255, 0, 255);
			tb = CLIP3((tb * (255 - mask) + pSmooth[0] * mask) / 255, 0, 255);
			tg = CLIP3((tg * (255 - mask) + pSmooth[1] * mask) / 255, 0, 255);
			tr = CLIP3((tr * (255 - mask) + pSmooth[2] * mask) / 255, 0, 255);

			pSrc[0] = CLIP3((pSrc[0] * (128 - k) + tb * k) >> 7, 0, 255);
			pSrc[1] = CLIP3((pSrc[1] * (128 - k) + tg * k) >> 7, 0, 255);
			pSrc[2] = CLIP3((pSrc[2] * (128 - k) + tr * k) >> 7, 0, 255);
			pSrc += stride;
			pGauss += stride;
			pGreen += stride;
			pSmooth += stride;
		}
	}
	free(gaussData);
	free(greenData);
	free(curveData);
	free(skinData);
	free(smoothData);
	return ret;
}
int F_Softskin_MixMethod(unsigned char* srcData, int width, int height, int stride, unsigned char* lightMap, int textureRatio, int ratio)
{
	int ret = 0;
	int length = width * height * stride;
	unsigned char* smoothData = (unsigned char*)malloc(sizeof(unsigned char) * length);
	unsigned char* tempData = (unsigned char*)malloc(sizeof(unsigned char) * length);
	unsigned char* hpData = (unsigned char*)malloc(sizeof(unsigned char) * length);
	memcpy(smoothData, srcData, sizeof(unsigned char) * length);
	memcpy(tempData, srcData, sizeof(unsigned char) * length);
	unsigned char* skinPDF = (unsigned char*)malloc(sizeof(unsigned char) * length);
	int smoothRadius = 8;
	int smoothThreshold = 38;
	int maskSmoothRadius = 3;
	ret = F_Softskin_A(smoothData, width, height, stride, lightMap, 95);
	memcpy(skinPDF, smoothData, sizeof(unsigned char) * length);
	ret = F_SkinProbability(skinPDF, width, height, stride);
	unsigned char* pSrc = srcData;
	unsigned char* pSkin = skinPDF;
	unsigned char* pSmooth = smoothData;
	unsigned char* pHP = hpData;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			int mask = 128 - textureRatio * 128 / 100;
			pSrc[0] = CLIP3((pSrc[0] * (128 - mask) + pSmooth[0] * mask) >> 7, 0, 255);
			pSrc[1] = CLIP3((pSrc[1] * (128 - mask) + pSmooth[1] * mask) >> 7, 0, 255);
			pSrc[2] = CLIP3((pSrc[2] * (128 - mask) + pSmooth[2] * mask) >> 7, 0, 255);
			pHP[0] = CLIP3(pSmooth[0] - pSrc[0] + 128, 0, 255);
			pHP[1] = CLIP3(pSmooth[1] - pSrc[1] + 128, 0, 255);
			pHP[2] = CLIP3(pSmooth[2] - pSrc[2] + 128, 0, 255);
			pHP += stride;
			pSmooth += stride;
			pSrc += stride;
			pSkin += stride;
		}
	}
	float hpRadius = 3.1;
	ret = F_FastGaussFilter(hpData, width, height, stride, hpRadius);
	pSrc = srcData;
	pSkin = skinPDF;
	pSmooth = smoothData;
	pHP = hpData;
	int k = ratio * 128 / 100;
	unsigned char* pTemp = tempData;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			int hpb = pHP[0];
			int hpg = pHP[1];
			int hpr = pHP[2];
			hpb = ModeLinearLight(pSrc[0], hpb);
			hpg = ModeLinearLight(pSrc[1], hpg);
			hpr = ModeLinearLight(pSrc[2], hpr);
			int mask = (pSkin[0] + pSkin[1] + pSkin[2]) / 3;
			hpb = CLIP3((hpb * mask + pSmooth[0] * (255 - mask)) / 255, 0, 255);
			hpg = CLIP3((hpg * mask + pSmooth[1] * (255 - mask)) / 255, 0, 255);
			hpr = CLIP3((hpr * mask + pSmooth[2] * (255 - mask)) / 255, 0, 255);

			pSrc[0] = CLIP3((hpb * k + pTemp[0] * (128 - k)) >> 7, 0, 255);
			pSrc[1] = CLIP3((hpg * k + pTemp[1] * (128 - k)) >> 7, 0, 255);
			pSrc[2] = CLIP3((hpr * k + pTemp[2] * (128 - k)) >> 7, 0, 255);
			pSrc += stride;
			pHP += stride;
			pSmooth += stride;
			pSkin += stride;
			pTemp += stride;
		}
	}
	free(skinPDF);
	free(smoothData);
	free(hpData);
	free(tempData);
	return ret;
};