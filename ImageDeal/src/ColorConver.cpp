#include "ColorConver.h"
void RGB2HSV(unsigned char R, unsigned char G, unsigned char B, float* h, float* s, float* v)
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
void HSV2RGB(float h, float s, float v, unsigned char* R, unsigned char* G, unsigned char* B)
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
void RGBToYCbCr(int R, int G, int B, int* Y, int* Cb, int* Cr)
{
	*Y = ((YCbCrYRI * R + YCbCrYGI * G + YCbCrYBI * B + HalfShiftValue) >> Shift);
	*Cb = (128 + ((YCbCrCbRI * R + YCbCrCbGI * G + YCbCrCbBI * B + HalfShiftValue) >> Shift));
	*Cr = (128 + ((YCbCrCrRI * R + YCbCrCrGI * G + YCbCrCrBI * B + HalfShiftValue) >> Shift));
}
void YCbCrToRGB(int Y, int Cb, int Cr, int* R, int* G, int* B)
{
	Cb = Cb - 128; Cr = Cr - 128;
	*R = Y + ((RGBRCrI * Cr + HalfShiftValue) >> Shift);
	*G = Y + ((RGBGCbI * Cb + RGBGCrI * Cr + HalfShiftValue) >> Shift);
	*B = Y + ((RGBBCbI * Cb + HalfShiftValue) >> Shift);
	if (*R > 255) *R = 255; else if (*R < 0) *R = 0;
	if (*G > 255) *G = 255; else if (*G < 0) *G = 0;
	if (*B > 255) *B = 255; else if (*B < 0) *B = 0;
}