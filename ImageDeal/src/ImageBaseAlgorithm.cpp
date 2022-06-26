#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdio.h>
#include"common.h"

int f_Saturation(unsigned char *srcData, int width, int height, int stride, int saturation)
{
	int ret = 0;
	if (saturation == 0)
		return ret;
	unsigned char* pSrc = srcData;
	int r, g, b, rgbMin, rgbMax;
	saturation = CLIP3(saturation, -100, 100);
	int k = saturation / 100.0f * 128;
	int alpha = 0;
	int offset = stride - width * 4;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			r = pSrc[2];
			g = pSrc[1];
			b = pSrc[0];
			rgbMin = MIN2(MIN2(r, g), b);
			rgbMax = MAX2(MAX2(r, g), b);
			int delta = (rgbMax - rgbMin);
			int value = (rgbMax + rgbMin);
			if (delta == 0)
			{
				pSrc += 4;
				continue;
			}
			int L = value >> 1;
			int S = L < 128 ? (delta << 7) / value : (delta << 7) / (510 - value);
			if (k >= 0)
			{
				alpha = k + S >= 128 ? S : 128 - k;
				alpha = 128 * 128 / alpha - 128;
			}
			else
				alpha = k;
			r = r + ((r - L) * alpha >> 7);
			g = g + ((g - L) * alpha >> 7);
			b = b + ((b - L) * alpha >> 7);
			pSrc[0] = CLIP3(b, 0, 255);
			pSrc[1] = CLIP3(g, 0, 255);
			pSrc[2] = CLIP3(r, 0, 255);
			pSrc += 4;
		}
		pSrc += offset;
	}
	return ret;
};

int f_BrightContrast(unsigned char *srcData, int width, int height, int stride, int bright, int contrast)
{
	int ret = 0;
	bright = CLIP3(bright, -100, 100);
	contrast = CLIP3(contrast, -100, 100);
	//compute average light of image
	int Average = 0;
	int offset = stride - width * 4;
	unsigned char* pSrc = srcData;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			Average += (299 * pSrc[2] + 587 * pSrc[1] + 114 * pSrc[0]) / 1000;
			pSrc += 4;
		}
		pSrc += offset;
	}
	Average = Average / (width * height);
	pSrc = srcData;
	unsigned char BC_MAP[256];
	int temp = 0;
	for (int i = 0; i < 256; i++)
	{
		int temp = contrast > 0 ? CLIP3(i + bright, 0, 255) : i;
		if (contrast > 0)
		{
			temp = CLIP3(i + bright, 0, 255);
			temp = CLIP3(Average + (temp - Average) * (1.0f / (1.0f - contrast / 100.0f)), 0, 255);
		}
		else
		{
			temp = i;
			temp = CLIP3(Average + (temp - Average) * (1.0f + contrast / 100.0f), 0, 255);
			temp = CLIP3(temp + bright, 0, 255);
		}
		BC_MAP[i] = temp;
	}
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			pSrc[0] = BC_MAP[pSrc[0]];
			pSrc[1] = BC_MAP[pSrc[1]];
			pSrc[2] = BC_MAP[pSrc[2]];
			pSrc += 4;
		}
		pSrc += offset;
	}
	return ret;
};
int f_Histagram(unsigned char *srcData, int width, int height, int stride, int hist[256], int mode)
{
	int ret = 0;
	int i, j, gray, offset;
	offset = stride - width * 4;
	unsigned char* pSrc = srcData;
	switch (mode)
	{
	case 0://Gray histagram
		for (j = 0; j < height; j++)
		{
			for (i = 0; i < width; i++)
			{
				gray = (pSrc[0] + pSrc[1] + pSrc[2]) / 3;
				hist[gray]++;
				pSrc += 4;
			}
			pSrc += offset;
		}
		break;
	case 1://Red histagram
		for (j = 0; j < height; j++)
		{
			for (i = 0; i < width; i++)
			{
				hist[pSrc[2]]++;
				pSrc += 4;
			}
			pSrc += offset;
		}
		break;
	case 2://Green histagram
		for (j = 0; j < height; j++)
		{
			for (i = 0; i < width; i++)
			{
				hist[pSrc[1]]++;
				pSrc += 4;
			}
			pSrc += offset;
		}
		break;
	case 3://Blue histagram
		for (j = 0; j < height; j++)
		{
			for (i = 0; i < width; i++)
			{
				hist[pSrc[0]]++;
				pSrc += 4;
			}
			pSrc += offset;
		}
		break;
	default:
		break;
	}

	return ret;
};
int f_Threshold(unsigned char *srcData, int width, int height, int stride, int T)
{
	int ret = 0;
	int i, j, gray, offset;
	offset = stride - width * 4;
	unsigned char* pSrc = srcData;
	for (j = 0; j < height; j++)
	{
		for (i = 0; i < width; i++)
		{
			gray = (pSrc[0] + pSrc[1] + pSrc[2]) / 3;
			gray = gray < T ? 0 : 255;
			pSrc[0] = gray;
			pSrc[1] = gray;
			pSrc[2] = gray;
			pSrc += 4;
		}
		pSrc += offset;
	}
	return ret;
};
int f_Gray(unsigned char *srcData, int width, int height, int stride, int mode)
{
	int ret = 0;
	int i, j, gray, offset;
	offset = stride - width * 4;
	unsigned char* pSrc = srcData;
	switch (mode)
	{
	case 0://mean gray method
		for (j = 0; j < height; j++)
		{
			for (i = 0; i < width; i++)
			{
				gray = (pSrc[0] + pSrc[1] + pSrc[2]) / 3;
				pSrc[0] = gray;
				pSrc[1] = gray;
				pSrc[2] = gray;
				pSrc += 4;
			}
			pSrc += offset;
		}
		break;
	case 1://classic gray method
		for (j = 0; j < height; j++)
		{
			for (i = 0; i < width; i++)
			{
				gray = (299 * pSrc[2] + 587 * pSrc[1] + 114 * pSrc[0]) / 1000;
				pSrc[0] = gray;
				pSrc[1] = gray;
				pSrc[2] = gray;
				pSrc += 4;
			}
			pSrc += offset;
		}
		break;
	case 2://photoshop gray method
		for (j = 0; j < height; j++)
		{
			for (i = 0; i < width; i++)
			{
				gray = (MAX2(pSrc[0], MAX2(pSrc[1], pSrc[2])) + MIN2(pSrc[0], MIN2(pSrc[1], pSrc[2]))) / 2;
				pSrc[0] = gray;
				pSrc[1] = gray;
				pSrc[2] = gray;
				pSrc += 4;
			}
			pSrc += offset;
		}
		break;
	default:
		break;
	}
	return ret;
};