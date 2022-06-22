#include "ImageDeal.h"
#include<string.h>
#include<cmath>
#define MIN2(a, b) ((a) < (b) ? (a) : (b))
#define MAX2(a, b) ((a) > (b) ? (a) : (b))
#define CLIP3(x, a, b) MIN2(MAX2(a,x), b)
const float YCbCrYRF = 0.299F;              // RGB转YCbCr的系数(浮点类型）
const float YCbCrYGF = 0.587F;
const float YCbCrYBF = 0.114F;
const float YCbCrCbRF = -0.168736F;
const float YCbCrCbGF = -0.331264F;
const float YCbCrCbBF = 0.500000F;
const float YCbCrCrRF = 0.500000F;
const float YCbCrCrGF = -0.418688F;
const float YCbCrCrBF = -0.081312F;

const float RGBRYF = 1.00000F;            // YCbCr转RGB的系数(浮点类型）
const float RGBRCbF = 0.0000F;
const float RGBRCrF = 1.40200F;
const float RGBGYF = 1.00000F;
const float RGBGCbF = -0.34414F;
const float RGBGCrF = -0.71414F;
const float RGBBYF = 1.00000F;
const float RGBBCbF = 1.77200F;
const float RGBBCrF = 0.00000F;

const int Shift = 20;
const int HalfShiftValue = 1 << (Shift - 1);

const int YCbCrYRI = (int)(YCbCrYRF * (1 << Shift) + 0.5);         // RGB转YCbCr的系数(整数类型）
const int YCbCrYGI = (int)(YCbCrYGF * (1 << Shift) + 0.5);
const int YCbCrYBI = (int)(YCbCrYBF * (1 << Shift) + 0.5);
const int YCbCrCbRI = (int)(YCbCrCbRF * (1 << Shift) + 0.5);
const int YCbCrCbGI = (int)(YCbCrCbGF * (1 << Shift) + 0.5);
const int YCbCrCbBI = (int)(YCbCrCbBF * (1 << Shift) + 0.5);
const int YCbCrCrRI = (int)(YCbCrCrRF * (1 << Shift) + 0.5);
const int YCbCrCrGI = (int)(YCbCrCrGF * (1 << Shift) + 0.5);
const int YCbCrCrBI = (int)(YCbCrCrBF * (1 << Shift) + 0.5);

const int RGBRYI = (int)(RGBRYF * (1 << Shift) + 0.5);              // YCbCr转RGB的系数(整数类型）
const int RGBRCbI = (int)(RGBRCbF * (1 << Shift) + 0.5);
const int RGBRCrI = (int)(RGBRCrF * (1 << Shift) + 0.5);
const int RGBGYI = (int)(RGBGYF * (1 << Shift) + 0.5);
const int RGBGCbI = (int)(RGBGCbF * (1 << Shift) + 0.5);
const int RGBGCrI = (int)(RGBGCrF * (1 << Shift) + 0.5);
const int RGBBYI = (int)(RGBBYF * (1 << Shift) + 0.5);
const int RGBBCbI = (int)(RGBBCbF * (1 << Shift) + 0.5);
const int RGBBCrI = (int)(RGBBCrF * (1 << Shift) + 0.5);
inline int abs(int a) {
	if (a < 0)return - a;
	return a;
}
float GetGaussianValue(float x, float mean, float var)
{
	float t = -0.5f * (x - mean) * (x - mean) / var;
	return exp(t);
}
void RGB2HSV(unsigned char R, unsigned char G, unsigned char B, float* h, float* s, float * v)
{
	float min, max;
	float r = R / 255.0f;
	float g = G / 255.0f;
	float b = B / 255.0f;
	min = MIN2(r, MIN2(g, b));
	max = MAX2(r, MAX2(g, b));
	if (max == min)
		*h = 0;
	else if (max == r && g >= b)
		*h = 60.0f * (g - b) / (max - min);
	else if (max == r && g < b)
		*h = 60.0f * (g - b) / (max - min) + 360.0f;
	else if (max == g)
		*h = 60.0f * (b - r) / (max - min) + 120.0f;
	else if (max == b)
		*h = 60.0f * (r - g) / (max - min) + 240.0f;
	*h = CLIP3(*h, 0, 360);
	if (max == 0)
		*s = 0;
	else
		*s = (max - min) / max;
	*v = max;
};
void HSV2RGB(float h, float s, float v, unsigned char* R, unsigned char *G, unsigned char *B)
{
	float q = 0, p = 0, t = 0, f = 0, r = 0, g = 0, b = 0;
	int hN = 0;
	if (h == 360)
		h = 0;
	if (h < 0)
		h = 360 + h;
	hN = (int)((int)h / 60);
	f = h / 60.0f - hN;
	p = v * (1.0f - s);
	q = v * (1.0f - f * s);
	t = v * (1.0f - (1.0f - f) * s);
	switch (hN)
	{
	case 0:
		r = v;
		g = t;
		b = p;
		break;
	case 1:
		r = q;
		g = v;
		b = p;
		break;
	case 2:
		r = p;
		g = v;
		b = t;
		break;
	case 3:
		r = p;
		g = q;
		b = v;
		break;
	case 4:
		r = t;
		g = p;
		b = v;
		break;
	case 5:
		r = v;
		g = p;
		b = q;
		break;
	default:
		break;
	}
	*R = (unsigned char)CLIP3((r * 255.0f), 0, 255);
	*G = (unsigned char)CLIP3((g * 255.0f), 0, 255);
	*B = (unsigned char)CLIP3((b * 255.0f), 0, 255);
};

void RGBToYCbCr(int R, int G, int B, int*Y, int*Cb, int* Cr)
{
	*Y = ((YCbCrYRI * R + YCbCrYGI * G + YCbCrYBI * B + HalfShiftValue) >> Shift);
	*Cb = (128 + ((YCbCrCbRI * R + YCbCrCbGI * G + YCbCrCbBI * B + HalfShiftValue) >> Shift));
	*Cr = (128 + ((YCbCrCrRI * R + YCbCrCrGI * G + YCbCrCrBI * B + HalfShiftValue) >> Shift));
}

void YCbCrToRGB(int Y, int Cb, int Cr, int*R, int*G, int* B)
{
	Cb = Cb - 128; Cr = Cr - 128;
	*R = Y + ((RGBRCrI * Cr + HalfShiftValue) >> Shift);
	*G = Y + ((RGBGCbI * Cb + RGBGCrI * Cr + HalfShiftValue) >> Shift);
	*B = Y + ((RGBBCbI * Cb + HalfShiftValue) >> Shift);
	if (*R > 255) *R = 255; else if (*R < 0) *R = 0;
	if (*G > 255) *G = 255; else if (*G < 0) *G = 0;
	if (*B > 255) *B = 255; else if (*B < 0) *B = 0;
}

int F_SkinDetectionBGR(unsigned char* srcData, int width, int height)
{
	int ret = 0;
	unsigned char* pSrc = srcData;
	int R, G, B;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			B = pSrc[0];
			G = pSrc[1];
			R = pSrc[2];
			if (!((R > 95) && (G > 40) && (B > 20) && (R > G) && (R > B) && (MAX2(R, MAX2(G, B)) - MIN2(R, MAX2(G, B)) > 15) && (abs(R - G) > 15)))
			{
				pSrc[0] = 0;
				pSrc[1] = 0;
				pSrc[2] = 0;
			}
			pSrc += 3;
		}
	}
	return ret;
}
int F_SkinDetectionHSV(unsigned char* srcData, int width, int height)
{
	int ret = 0;
	unsigned char* pSrc = srcData;
	int R, G, B;
	float H, S, V;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			B = pSrc[0];
			G = pSrc[1];
			R = pSrc[2];
			RGB2HSV(R, G, B, &H, &S, &V);
			if (!(((S >= 0.1) && (S <= 0.68)) && ((V >= 0.13 && V <= 0.25 && H >= ((0.4 - V) / 0.014) && H <= ((V + 0.062) / 0.01)) || ((V > 0.25) && (V <= 0.38) && (H >= (0.4 - V) / 0.014) && (H <= (0.67 - V) / 0.014)) || ((V > 0.38) && (V <= 0.46) && (H >= (V - 0.34) / 0.03) && (H <= (0.67 - V) / 0.014)) ||
				((V > 0.46) && (V <= 0.6) && (H >= (V - 0.34) / 0.03) && (H <= (V - 0.31) / 0.009)) || ((V > 0.6) && (V <= 0.76) && (H >= (0.91 - V) / 0.14) && (H <= (V - 0.31) / 0.009)) || ((V > 0.76) && (V <= 0.91) && (H >= (0.91 - V) / 0.14) && (H <= (1.17 - V) / 0.0082)) || ((V > 0.91) && (V <= 1) && (H >= (V - 0.91) / 0.0041) && (H <= (1.17 - V) / 0.0082)))))
			{
				pSrc[0] = 0;
				pSrc[1] = 0;
				pSrc[2] = 0;
			}
			pSrc += 3;
		}
	}
	return ret;
}
int F_SkinDetectionYCgCr(unsigned char* srcData, int width, int height)
{
	int ret = 0;
	unsigned char* pSrc = srcData;
	int R, G, B;
	float Cr, Cg;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			B = pSrc[0];
			G = pSrc[1];
			R = pSrc[2];
			Cg = 128 - 0.318f * R + 0.4392f * G - 0.1212f * B;
			Cr = 128 + 0.4392f * R - 0.3677f * G - 0.0714f * B;
			if (!((Cg >= 85) && (Cg <= 135) && ((Cr <= (280 - Cg)) && (Cr >= (260 - Cg)))))
			{
				pSrc[0] = 0;
				pSrc[1] = 0;
				pSrc[2] = 0;
			}
			pSrc += 3;
		}
	}
	return ret;
}
float GetProbability(int R, int G, int B, float meanCb, float varCb, float meanCr, float varCr)
{
	int Y, Cb, Cr;
	RGBToYCbCr(R, G, B, &Y, &Cb, &Cr);
	float pcb = GetGaussianValue(Cb, meanCb, varCb);
	float pcr = GetGaussianValue(Cr, meanCr, varCr);
	return 2.0f * pcb * pcr;
};
void GaussMask(int r, double sigma, double gaussMask[])
{
	double PI = 3.1415926;
	double sum = 0;
	int stride = 2 * r + 1;
	for (int y = -r, h = 0; y <= r; y++, h++)
	{
		for (int x = -r, w = 0; x <= r; x++, w++)
		{
			gaussMask[w + h * stride] = (1.0 / (2.0 * PI * sigma * sigma)) * (exp(-((double)x * (double)x + (double)y * (double)y) / (2.0 * sigma * sigma)));
			sum += gaussMask[w + h * stride];
		}
	}
	for (int i = 0; i < stride * stride; i++)
	{
		gaussMask[i] = gaussMask[i] / sum;
	}
};
int F_FastGaussFilter(unsigned char* srcData, int width, int height, int stride, float r)
{
	int ret = 0;
	int radius = (int)r;
	if (r == 0)
		return ret;
	unsigned char* dstData = (unsigned char*)malloc(sizeof(unsigned char)*height*stride);
	unsigned char* tempData = (unsigned char*)malloc(sizeof(unsigned char)*height*stride);
	int totalWei = 0;
	int i, j, k;
	float sigma = r;
	unsigned char *kernel = (unsigned char *)malloc(2 * radius + 1);
	for (i = -radius; i <= radius; i++)
	{
		kernel[i + radius] = (unsigned char)(exp(-(float)i*i / (2 * sigma*sigma)) * 128);
		totalWei += kernel[i + radius];
	}
	int tempR = 0, tempG = 0, tempB = 0;
	int v = 0;
	int K = 0;
	int rem = 0;
	int t = 0;
	int offset = stride - width * 4;
	for (j = 0; j < height; j++)
	{
		for (i = 0; i < width; i++)
		{
			tempR = 0; tempG = 0; tempB = 0;
			for (k = -radius; k <= radius; k++)
			{
				rem = (abs(i + k) % width);
				t = rem * 4 + j * stride;
				K = kernel[k + radius];
				tempB += srcData[t] * K;
				tempG += srcData[t + 1] * K;
				tempR += srcData[t + 2] * K;
			}
			v = i * 4 + j * stride;
			tempData[v] = tempB / totalWei;
			tempData[v + 1] = tempG / totalWei;
			tempData[v + 2] = tempR / totalWei;
		}
	}
	for (i = 0; i < width; i++)
	{
		for (j = 0; j < height; j++)
		{
			tempR = 0; tempG = 0; tempB = 0;
			for (k = -radius; k <= radius; k++)
			{
				rem = (abs(j + k) % height);
				t = rem * stride + i * 4;
				K = kernel[k + radius];
				tempB += tempData[t] * K;
				tempG += tempData[t + 1] * K;
				tempR += tempData[t + 2] * K;
			}
			v = i * 4 + j * stride;
			dstData[v] = tempB / totalWei;
			dstData[v + 1] = tempG / totalWei;
			dstData[v + 2] = tempR / totalWei;
		}
	}
	memcpy(srcData, dstData, sizeof(unsigned char) * height * stride);
	free(dstData);
	free(tempData);
	return ret;
};
int F_SkinProbability(unsigned char* srcData, int width, int height)
{
	int ret = 0;
	float sum = 0, mean = 0, variance = 0;
	unsigned char* pSrc = srcData;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			//default setting is computed using special skin data.
			//meanCb-varCb:102-196
			//meanCr-varCr:143-196
			int gray = CLIP3(GetProbability(pSrc[2], pSrc[1], pSrc[0], 102, 196, 143, 196) * 255.0f, 0, 255);
			pSrc[0] = pSrc[1] = pSrc[2] = gray;
			pSrc += 3;
		}
	}
	return ret;
};
int SmartBlurOneChannel(unsigned char* srcData, int width, int height, int radius, int threshold)
{
	int len = sizeof(unsigned long) * width * height;
	int i, j;
	int gray = 0;
	unsigned char* tempData = (unsigned char*)malloc(sizeof(unsigned char) * height * width);
	memcpy(tempData, srcData, sizeof(unsigned char) * height * width);
	for (j = 0; j < height; j++)
	{
		for (i = 0; i < width; i++)
		{
			len = i + j * width;
			gray = tempData[len];
			int low = CLIP3(gray - threshold, 0, 255);
			int high = CLIP3(gray + threshold, 0, 255);
			int sum = 0;
			int count = 0;
			for (int n = -radius; n <= radius; n++)
			{
				for (int m = -radius; m <= radius; m++)
				{
					int x = CLIP3(i + m, 0, width - 1);
					int y = CLIP3(j + n, 0, height - 1);
					int pos = x + y * width;
					gray = tempData[pos];
					if (gray > low && gray < high)
					{
						sum += gray;
						count++;
					}
				}
			}
			gray = count == 0 ? srcData[len] : sum / count;//sum / MAX2(count, 1);						
			srcData[len] = CLIP3(gray, 0, 255);
		}
	}
	free(tempData);
	return 0;
};
int F_SmartBlur(unsigned char* srcData, int nWidth, int nHeight, int nStride, int radius, int threshold)
{
	int ret = 0;
	if (srcData == NULL)
	{
		return ret;
	}
	if (radius == 0 || threshold == 0)
		return ret;
	unsigned char* yData = (unsigned char*)malloc(sizeof(unsigned char) * nWidth * nHeight);
	unsigned char* cbData = (unsigned char*)malloc(sizeof(unsigned char) * nWidth * nHeight);
	unsigned char* crData = (unsigned char*)malloc(sizeof(unsigned char) * nWidth * nHeight);
	unsigned char* pSrc = srcData;
	int Y, CB, CR;
	unsigned char* pY = yData;
	unsigned char* pCb = cbData;
	unsigned char* pCr = crData;
	for (int j = 0; j < nHeight; j++)
	{
		for (int i = 0; i < nWidth; i++)
		{

			RGBToYCbCr(pSrc[2], pSrc[1], pSrc[0], &Y, &CB, &CR);
			*pY = Y;
			*pCb = CB;
			*pCr = CR;
			pY++;
			pCb++;
			pCr++;
			pSrc += 4;
		}
	}
	SmartBlurOneChannel(yData, nWidth, nHeight, radius, threshold);
	pSrc = srcData;
	pY = yData;
	pCb = cbData;
	pCr = crData;
	int R, G, B;
	for (int j = 0; j < nHeight; j++)
	{
		for (int i = 0; i < nWidth; i++)
		{
			YCbCrToRGB(*pY, *pCb, *pCr, &R, &G, &B);
			pSrc[0] = B;
			pSrc[1] = G;
			pSrc[2] = R;
			pY++;
			pCb++;
			pCr++;
			pSrc += 4;
		}
	}
	free(yData);
	free(cbData);
	free(crData);
	return ret;
}
int f_Softskin(unsigned char* srcData, int width, int height, int stride, int skinMode, int ratio)
{
	int ret = 0;
	int length = height * stride;
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
		ret = F_SkinDetectionYCgCr(skinProbability, width, height);
		maskSmoothRadius = 6;
	}
	else
	{
		ret = F_SkinProbability(skinProbability, width, height);
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
			pSrc += 4;
			pSmooth += 4;
			pMask += 4;
		}
	}
	free(skinProbability);
	free(smoothData);
	return ret;
};