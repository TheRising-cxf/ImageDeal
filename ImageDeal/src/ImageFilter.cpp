#include "ImageFilter.h"
int F_FastGaussFilter(unsigned char* srcData, int width, int height, int stride, int r)
{
	int ret = 0;
	int radius = r;
	if (r == 0)
		return ret;
	unsigned char* dstData = (unsigned char*)malloc(sizeof(unsigned char)* width * height * stride);
	unsigned char* tempData = (unsigned char*)malloc(sizeof(unsigned char) * width * height * stride);
	int totalWei = 0;
	int i, j, k;
	float sigma = r;
	unsigned char* kernel = (unsigned char*)malloc(2 * radius + 1);
	for (i = -radius; i <= radius; i++)
	{
		kernel[i + radius] = (unsigned char)(exp(-(float)i * i / (2 * sigma * sigma)) * 128);
		totalWei += kernel[i + radius];
	}
	int tempR = 0, tempG = 0, tempB = 0;
	int v = 0;
	int K = 0;
	int rem = 0;
	int t = 0;
	int offset = width * stride;
	for (j = 0; j < height; j++)
	{
		for (i = 0; i < width; i++)
		{
			tempR = 0; tempG = 0; tempB = 0;
			for (k = -radius; k <= radius; k++)
			{
				rem = (abs(i + k) % width);
				t = rem * stride + j * offset;
				K = kernel[k + radius];
				tempB += srcData[t] * K;
				tempG += srcData[t + 1] * K;
				tempR += srcData[t + 2] * K;
			}
			v = i * stride + j * offset;
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
				t = rem * offset + i * stride;
				K = kernel[k + radius];
				tempB += tempData[t] * K;
				tempG += tempData[t + 1] * K;
				tempR += tempData[t + 2] * K;
			}
			v = i * stride + j * offset;
			dstData[v] = tempB / totalWei;
			dstData[v + 1] = tempG / totalWei;
			dstData[v + 2] = tempR / totalWei;
		}
	}
	memcpy(srcData, dstData, sizeof(unsigned char) * width * height * stride);
	free(dstData);
	free(tempData);
	return ret;
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
	if (nStride == 1) {
		int len = sizeof(unsigned long) * nWidth * nHeight;
		int i, j;
		int gray = 0;
		unsigned char* tempData = (unsigned char*)malloc(sizeof(unsigned char) * nWidth * nHeight);
		memcpy(tempData, srcData, sizeof(unsigned char) * nWidth * nHeight);
		for (j = 0; j < nHeight; j++)
		{
			for (i = 0; i < nWidth; i++)
			{
				len = i + j * nWidth;
				gray = tempData[len];
				int low = CLIP3(gray - threshold, 0, 255);
				int high = CLIP3(gray + threshold, 0, 255);
				int sum = 0;
				int count = 0;
				for (int n = -radius; n <= radius; n++)
				{
					for (int m = -radius; m <= radius; m++)
					{
						int x = CLIP3(i + m, 0, nWidth - 1);
						int y = CLIP3(j + n, 0, nHeight - 1);
						int pos = x + y * nWidth;
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
	}
	else if (nStride == 3) {

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
				pSrc += nStride;
			}
		}
		F_SmartBlur(yData, nWidth, nHeight, 1, radius, threshold);
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
				pSrc += nStride;
			}
		}
		free(yData);
		free(cbData);
		free(crData);
	}
	return ret;
}
