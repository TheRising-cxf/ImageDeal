#include "FindSkin.h"
#include "ColorConver.h"
#include <cmath>
int F_SkinDetectionBGR(unsigned char* srcData, int width, int height, int channels)
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
			pSrc += channels;
		}
	}
	return ret;
}
int F_SkinDetectionHSV(unsigned char* srcData, int width, int height, int channels)
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
			pSrc += channels;
		}
	}
	return ret;
}
int F_SkinDetectionYCgCr(unsigned char* srcData, int width, int height, int channels)
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
			pSrc += channels;
		}
	}
	return ret;
}
int F_SkinProbability(unsigned char* srcData, int width, int height, int channels)
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
			pSrc += channels;
		}
	}
	return ret;
};