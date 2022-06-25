#pragma once
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"FindSkin.h"
#include"ImageFilter.h"
#include"ColorSkin.h"
int F_SkinWhite(unsigned char* srcData, int width, int height, int channels, unsigned char* lutData, int ratio)
{
	int ret = 0;
	int length = width * height * channels;
	unsigned char* tempData = (unsigned char*)malloc(sizeof(unsigned char) * length);
	memcpy(tempData, srcData, sizeof(unsigned char) * length);
	unsigned char* skinPDF = (unsigned char*)malloc(sizeof(unsigned char) * length);
	memcpy(skinPDF, srcData, sizeof(unsigned char) * length);
	ret = F_SkinProbability(skinPDF, width, height, channels);
	int maskSmoothRadius = 3;
	ret = F_FastGaussFilter(skinPDF, width, height, channels, maskSmoothRadius);
	ret = F_Filter512(tempData, width, height, channels, lutData);
	unsigned char* pSrc = srcData;
	unsigned char* pLut = tempData;
	unsigned char* pSkin = skinPDF;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			int r, g, b, a;
			b = CLIP3((pSrc[0] * (100 - ratio) + pLut[0] * ratio) / 100, 0, 255);
			g = CLIP3((pSrc[1] * (100 - ratio) + pLut[1] * ratio) / 100, 0, 255);
			r = CLIP3((pSrc[2] * (100 - ratio) + pLut[2] * ratio) / 100, 0, 255);
			a = (pSkin[0] + pSkin[1] + pSkin[2]) / 3;
			pSrc[0] = CLIP3((b * a + pSrc[0] * (255 - a)) / 255, 0, 255);
			pSrc[1] = CLIP3((g * a + pSrc[1] * (255 - a)) / 255, 0, 255);
			pSrc[2] = CLIP3((r * a + pSrc[2] * (255 - a)) / 255, 0, 255);
			pSrc += channels;
			pLut += channels;
			pSkin += channels;
		}
	}
	free(tempData);
	free(skinPDF);
	return ret;
};
int F_SkinWhiteCurve(unsigned char* srcData, int width, int height, int channels, int belta, int ratio)
{
	int ret = 0;
	int length = width * height * channels;
	unsigned char* skinPDF = (unsigned char*)malloc(sizeof(unsigned char) * length);
	memcpy(skinPDF, srcData, sizeof(unsigned char) * length);
	ret = F_SkinProbability(skinPDF, width, height, channels);
	int maskSmoothRadius = 3;
	ret = F_FastGaussFilter(skinPDF, width, height, channels, maskSmoothRadius);
	unsigned char* pSrc = srcData;
	unsigned char* pSkin = skinPDF;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			int r, g, b, a;
			//skin white curve
			b = CLIP3(log((float)pSrc[0] * (belta - 1) / 255.0f + 1) / log((float)belta) * 255.0f, 0, 255);
			g = CLIP3(log((float)pSrc[1] * (belta - 1) / 255.0f + 1) / log((float)belta) * 255.0f, 0, 255);
			r = CLIP3(log((float)pSrc[2] * (belta - 1) / 255.0f + 1) / log((float)belta) * 255.0f, 0, 255);
			b = CLIP3((b * ratio + pSrc[0] * (100 - ratio)) / 100, 0, 255);
			g = CLIP3((g * ratio + pSrc[1] * (100 - ratio)) / 100, 0, 255);
			r = CLIP3((r * ratio + pSrc[2] * (100 - ratio)) / 100, 0, 255);
			//skin pdf
			a = (pSkin[0] + pSkin[1] + pSkin[2]) / 3;
			pSrc[0] = CLIP3((b * a + pSrc[0] * (255 - a)) / 255, 0, 255);
			pSrc[1] = CLIP3((g * a + pSrc[1] * (255 - a)) / 255, 0, 255);
			pSrc[2] = CLIP3((r * a + pSrc[2] * (255 - a)) / 255, 0, 255);
			pSrc += channels;
			pSkin += channels;
		}
	}
	free(skinPDF);
	return ret;
}
int F_SkinWhitePS(unsigned char* srcData, int width, int height, int channels, int ratio)
{
	int ret = 0;
	int length = width * height * channels;
	unsigned char* skinPDF = (unsigned char*)malloc(sizeof(unsigned char) * length);
	memcpy(skinPDF, srcData, sizeof(unsigned char) * length);
	ret = F_SkinProbability(skinPDF, width, height, channels);
	int maskSmoothRadius = 3;
	ret = F_FastGaussFilter(skinPDF, width, height, channels, maskSmoothRadius);
	unsigned char* pSrc = srcData;

	unsigned char* pSkin = skinPDF;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			int r, g, b, a;
			//skin white using smoothlight of ps
			b = ModeSmoothLight(pSrc[0], 255);
			g = ModeSmoothLight(pSrc[1], 255);
			r = ModeSmoothLight(pSrc[2], 255);
			b = CLIP3((b * ratio + pSrc[0] * (100 - ratio)) / 100, 0, 255);
			g = CLIP3((g * ratio + pSrc[1] * (100 - ratio)) / 100, 0, 255);
			r = CLIP3((r * ratio + pSrc[2] * (100 - ratio)) / 100, 0, 255);
			//skin pdf
			a = (pSkin[0] + pSkin[1] + pSkin[2]) / 3;
			pSrc[0] = CLIP3((b * a + pSrc[0] * (255 - a)) / 255, 0, 255);
			pSrc[1] = CLIP3((g * a + pSrc[1] * (255 - a)) / 255, 0, 255);
			pSrc[2] = CLIP3((r * a + pSrc[2] * (255 - a)) / 255, 0, 255);
			pSrc += channels;
			pSkin += channels;
		}
	}
	free(skinPDF);
	return ret;
}