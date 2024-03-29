#include "ImageFilter.h"
#include "ColorConver.h"
#include <cmath>
#include <string.h>
int F_BeautyCamera_FilterMode(unsigned char* srcData, int width, int height, int channels, int Mode) {
	int ret = 0;
	if (Mode == 0) {
		int radius = 3;
		int deltaD = 6;
		int deltaR = 10;
		F_StdBilateralFilter(srcData, width, height, channels, radius, deltaD, deltaR);
	}
	else if (Mode == 1) {
		int hpRadius = 10.0f * width * height / (594 * 677);
		F_FastGaussFilter(srcData, width, height, channels, hpRadius);
	}
	else if (Mode == 2) {
		F_SmartBlur(srcData, width, height, channels, 5, 30);
	}
	else if (Mode == 3) {
		F_SurfaceBlur(srcData, width, height, channels, 8, 38);
	}
	else if (Mode == 4) {
		F_GuidedFilter(srcData, width, height, channels, 10, 0.02);
	}
	else if (Mode == 5) {
		F_LSNFilter(srcData, width, height, channels, 5, 2500);
	}
	else if (Mode == 6) {
		F_AnisotropicFilter(srcData, width, height, channels, 10, 10, 0.33, 1);
	}
	else if (Mode == 7) {
		F_MeanShiftFilter(srcData, width, height, channels, 5, 3, 3);
	}
	else if (Mode == 9) {
		F_BeepsFilter(srcData, width, height, channels, 10, 20);
	}
	return ret;
}
int F_StdBilateralFilter(unsigned char* srcData, int nWidth, int nHeight, int channels, int radius,int deltaD,int deltaR) {
	int ret = 0;
	if (srcData == NULL)
	{
		return ret;
	}
	if (deltaD == 0 || deltaR == 0)
		return ret;
	if (channels == 1) {
		int len = sizeof(unsigned long) * nWidth * nHeight;
		unsigned char* tempData = (unsigned char*)malloc(sizeof(unsigned char) * nWidth * nHeight);
		memcpy(tempData, srcData, sizeof(unsigned char) * nWidth * nHeight);
		int cx, cy;
		float sum_w = 0;
		float w = 0;
		double sum = 0;
		float w1 = 1.0f / (nWidth*nWidth);
		float h1 = 1.0f / (nHeight*nHeight);
		float deltad1 = (float)1.0f / (2.0f*deltaD*deltaD);
		float deltar1 = (float)1.0f / (2.0f*deltaR*deltaR);
		for (int j = 0; j < nHeight; j++)
		{
			for (int i = 0; i < nWidth; i++)
			{
				sum_w = 0;
				sum = 0;
				for (int n = -radius; n < radius; n++) {
					for (int m = -radius; m < radius; m++) {
						cx = CLIP3(i + m, 0, nWidth - 1);
						cy = CLIP3(j + n, 0, nHeight - 1);
						float md = (i - m)*(i - m)*w1 + (j - n)*(j - n)*h1;
						float mr = (tempData[cx + cy * nWidth] - tempData[i + j * nWidth])*(tempData[cx + cy * nWidth] - tempData[i + j * nWidth]) / (255.0f*255.0f);
						md = -md * deltad1;
						mr = -mr * deltar1;
						w = exp(md + mr);
						sum_w += w;
						sum += w * tempData[cx + cy * nWidth];
					}
				}
				srcData[i + j * nWidth] = (unsigned char)(sum_w == 0 ? srcData[i + j * nWidth] : CLIP3(sum / sum_w, 0, 255));
			}
		}
		free(tempData);
	}
	else if (channels == 3) {

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
				pSrc += channels;
			}
		}
		F_StdBilateralFilter(yData, nWidth, nHeight, 1, radius,deltaD, deltaR);
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
				pSrc += channels;
			}
		}
		free(yData);
		free(cbData);
		free(crData);
	}
	return ret;
}
int F_SurfaceBlur(unsigned char* srcData, int width, int height, int channels, int radius, int threshold)
{
	if (srcData == NULL || radius == 0 || threshold == 0)
	{
		return 0;
	}
	if (channels == 1) {
		float matrixItems[511];//255*2+1
		float* map = &matrixItems[255];
		float fv = threshold * 2.5f;
		for (int i = 0; i < 256; i++)
		{
			map[-i] = map[i] = MAX2(1 - i / fv, 0);
		}
		int ret = 0;
		int i, j, n, m, len;
		len = sizeof(unsigned char) * width * height;
		unsigned char* tempData = (unsigned char*)malloc(len);
		int* tmap = (int*)malloc(sizeof(int) * height);
		if (tempData == NULL || tmap == NULL)
		{
			ret = 1;
			return ret;
		}
		if (NULL == tempData || NULL == tmap)
			return 1;
		for (i = 0; i < height; i++)
		{
			tmap[i] = i * width;
		}
		memcpy(tempData, srcData, len);
		int kx, ky;
		len = (radius << 1) + 1;
		int gray = 0;
		float sum, sum_a;
		int pos, pos_a;
		unsigned char val;
		for (j = 0; j < height; j++)
		{
			for (i = 0; i < width; i++)
			{
				pos = i + tmap[j];
				val = tempData[pos];
				sum = 0;
				sum_a = 0;
				for (n = -radius; n <= radius; n++)
				{
					ky = CLIP3(j + n, 0, height - 1);
					pos_a = tmap[ky];
					for (m = -radius; m <= radius; m++)
					{
						kx = CLIP3(i + m, 0, width - 1);
						gray = tempData[kx + pos_a];
						sum_a += map[gray - val];
						sum += gray * map[gray - val];
					}
				}
				gray = sum_a == 0 ? gray : sum / sum_a;//(int)(sum / MAX2(sum_a, 0.1));
				srcData[pos] = gray;//CLIP3(gray, 0 , 255);
			}
		}
		free(tempData);
		free(tmap);
	}
	else if (channels == 3) {
		unsigned char* yData = (unsigned char*)malloc(sizeof(unsigned char) * width * height);
		unsigned char* cbData = (unsigned char*)malloc(sizeof(unsigned char) * width * height);
		unsigned char* crData = (unsigned char*)malloc(sizeof(unsigned char) * width * height);
		unsigned char* pSrc = srcData;
		int Y, CB, CR;
		unsigned char* pY = yData;
		unsigned char* pCb = cbData;
		unsigned char* pCr = crData;
		for (int j = 0; j < height; j++)
		{
			for (int i = 0; i < width; i++)
			{
				RGBToYCbCr(pSrc[2], pSrc[1], pSrc[0], &Y, &CB, &CR);
				*pY = Y;
				*pCb = CB;
				*pCr = CR;
				pY++;
				pCb++;
				pCr++;
				pSrc += channels;
			}
		}
		F_SurfaceBlur(yData, width, height, 1, radius, threshold);
		F_SurfaceBlur(cbData, width, height, 1, radius, threshold);
		F_SurfaceBlur(crData, width, height, 1, radius, threshold);
		pSrc = srcData;
		pY = yData;
		pCb = cbData;
		pCr = crData;
		int R, G, B;
		for (int j = 0; j < height; j++)
		{
			for (int i = 0; i < width; i++)
			{
				YCbCrToRGB(*pY, *pCb, *pCr, &R, &G, &B);
				pSrc[0] = B;
				pSrc[1] = G;
				pSrc[2] = R;
				pY++;
				pCb++;
				pCr++;
				pSrc += channels;
			}
		}
		free(yData);
		free(cbData);
		free(crData);
	}
	return 0;
}
//��˹�˲�
int F_FastGaussFilter(unsigned char* srcData, int width, int height, int channels, int r)
{
	int ret = 0;
	int radius = r;
	if (r == 0)
		return ret;
	unsigned char* dstData = (unsigned char*)malloc(sizeof(unsigned char)* width * height * channels);
	unsigned char* tempData = (unsigned char*)malloc(sizeof(unsigned char) * width * height * channels);
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
	int offset = width * channels;
	for (j = 0; j < height; j++)
	{
		for (i = 0; i < width; i++)
		{
			tempR = 0; tempG = 0; tempB = 0;
			for (k = -radius; k <= radius; k++)
			{
				rem = (abs(i + k) % width);
				t = rem * channels + j * offset;
				K = kernel[k + radius];
				tempB += srcData[t] * K;
				tempG += srcData[t + 1] * K;
				tempR += srcData[t + 2] * K;
			}
			v = i * channels + j * offset;
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
				t = rem * offset + i * channels;
				K = kernel[k + radius];
				tempB += tempData[t] * K;
				tempG += tempData[t + 1] * K;
				tempR += tempData[t + 2] * K;
			}
			v = i * channels + j * offset;
			dstData[v] = tempB / totalWei;
			dstData[v + 1] = tempG / totalWei;
			dstData[v + 2] = tempR / totalWei;
		}
	}
	memcpy(srcData, dstData, sizeof(unsigned char) * width * height * channels);
	free(dstData);
	free(tempData);
	return ret;
};
int F_SmartBlur(unsigned char* srcData, int nWidth, int nHeight, int channels, int radius, int threshold)
{
	int ret = 0;
	if (srcData == NULL)
	{
		return ret;
	}
	if (radius == 0 || threshold == 0)
		return ret;
	if (channels == 1) {
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
	else if (channels == 3) {

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
				pSrc += channels;
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
				pSrc += channels;
			}
		}
		free(yData);
		free(cbData);
		free(crData);
	}
	return ret;
}
int F_Filter512(unsigned  char*  srcData, int  width, int  height, int  channels, unsigned  char*Map)
{
	int  i, j, r, g, b, offset, pos, nx, ny, k;
	unsigned  char*  pSrc = srcData;
	offset = channels * width;
	for (j = 0; j < height; j++)
	{
		for (i = 0; i < width; i++)
		{
			b = pSrc[0];
			g = pSrc[1];
			r = pSrc[2];
			k = (b >> 2);
			nx = (int)(r >> 2) + ((k - ((k >> 3) << 3)) << 6);
			ny = (int)(((b >> 5) << 6) + (g >> 2));
			pos = (nx * channels) + (ny * 512 * channels);
			pSrc[0] = Map[pos];
			pSrc[1] = Map[pos + 1];
			pSrc[2] = Map[pos + 2];
			pSrc += channels;
		}
	}
	return    0;
};
int F_MeanCovMapCalculate(float* srcData, int width, int height, float* meanData, int radius)
{
	float* stdData = (float*)malloc(sizeof(float) * width * height);
	int i, j, k;
	int pos = 0;
	for (j = 0; j < height; j++)
	{
		for (i = 0; i < width; i++)
		{
			pos = i + j * width;
			stdData[pos] = srcData[pos] * srcData[pos];
		}
	}
	if (srcData == NULL || meanData == NULL)
		return 1;
	if (radius == 0)
		return 0;
	if (radius > MIN2(width, height) / 2)
		radius = (int)(MIN2(width, height) / 2 - 0.5);
	int t = 0, t1 = 0;
	int block = (radius << 1) + 1;
	int winSize = block * block;
	float sumMean = 0;
	float* pSrc = srcData;
	float* tempMean = (float*)malloc(sizeof(float)* width);
	memset(tempMean, 0, sizeof(float) * width);
	for (k = -radius; k <= radius; k++)
	{
		for (j = 0; j < width; j++)
		{
			t1 = abs(k) * width;
			tempMean[j] += pSrc[j + t1];
		}
	}
	for (i = 0; i < height; i++)
	{
		sumMean = 0;
		for (j = -radius; j <= radius; j++)
		{
			t = abs(j);
			sumMean += tempMean[t];
		}
		for (j = 0; j < width; j++)
		{
			t = j + i * width;
			meanData[t] = (float)(sumMean / winSize);
			if (j < width - 1)
			{
				t = abs(j - radius);
				t1 = (j + radius + 1) % width;
				sumMean = sumMean - tempMean[t] + tempMean[t1];
			}
		}
		if (i < height - 1)
		{
			for (k = 0; k < width; k++)
			{
				t = k + abs(i - radius) * width;
				t1 = k + (i + radius + 1) % height * width;
				tempMean[k] = tempMean[k] - pSrc[t] + pSrc[t1];
			}
		}
	}
	free(tempMean);
	return 0;
};
int F_MeanCovMapCalculate(unsigned char* srcData, int width, int height,unsigned long* meanData, unsigned long* covData, int radius)
{
	unsigned long* stdData = (unsigned long*)malloc(sizeof(unsigned long) * width * height);
	int i, j, k;
	int pos = 0;
	for (j = 0; j < height; j++)
	{
		for (i = 0; i < width; i++)
		{
			pos = i + j * width;
			stdData[pos] = srcData[pos] * srcData[pos];
		}
	}
	if (srcData == NULL || meanData == NULL)
		return 1;
	if (radius == 0)
		return 0;
	if (radius > MIN2(width, height) / 2)
		radius = (int)(MIN2(width, height) / 2 - 0.5);
	int t = 0, t1 = 0;
	int block = (radius << 1) + 1;
	int winSize = block * block;
	long sumMean = 0;
	unsigned char* pSrc = srcData;
	int* tempMean = (int*)malloc(sizeof(int)* width);
	memset(tempMean, 0, sizeof(int) * width);
	long sumCov = 0;
	int* tempCov = (int*)malloc(sizeof(int)* width);
	memset(tempCov, 0, sizeof(int) * width);
	for (k = -radius; k <= radius; k++)
	{
		for (j = 0; j < width; j++)
		{
			t1 = abs(k) * width;
			tempMean[j] += pSrc[j + t1];
			tempCov[j] += stdData[j + t1];
		}
	}
	for (i = 0; i < height; i++)
	{
		sumMean = 0, sumCov = 0;
		for (j = -radius; j <= radius; j++)
		{
			t = abs(j);
			sumMean += tempMean[t];
			sumCov += tempCov[t];
		}
		for (j = 0; j < width; j++)
		{
			t = j + i * width;
			meanData[t] = (int)(sumMean / winSize);
			covData[t] = (int)(sumCov / winSize);
			if (j < width - 1)
			{
				t = abs(j - radius);
				t1 = (j + radius + 1) % width;
				sumMean = sumMean - tempMean[t] + tempMean[t1];
				sumCov = sumCov - tempCov[t] + tempCov[t1];
			}
		}
		if (i < height - 1)
		{
			for (k = 0; k < width; k++)
			{
				t = k + abs(i - radius) * width;
				t1 = k + (i + radius + 1) % height * width;
				tempMean[k] = tempMean[k] - pSrc[t] + pSrc[t1];
				tempCov[k] = tempCov[k] - stdData[t] + stdData[t1];
			}
		}
	}
	free(tempMean);
	free(tempCov);
	return 0;
};
int F_GuidedFilter(unsigned char* srcData, int nWidth, int nHeight, int channels, int radius, float delta)
{
	if (srcData == NULL)
	{
		return 0;
	}
	if (channels == 1) {
		float *data = (float*)malloc(sizeof(float) * nWidth * nHeight);
		float *meanIP = (float*)malloc(sizeof(float) * nWidth * nHeight);
		float *corrIP = (float*)malloc(sizeof(float) * nWidth * nHeight);
		float *varIP = (float*)malloc(sizeof(float) * nWidth * nHeight);
		float *a = (float*)malloc(sizeof(float) * nWidth * nHeight);
		float *b = (float*)malloc(sizeof(float) * nWidth * nHeight);
		for (int i = 0; i < nWidth * nHeight; i++)
		{
			data[i] = (float)srcData[i] / 255.0;
		}
		//mean and cov compute
		F_MeanCovMapCalculate(data, nWidth, nHeight, meanIP, radius);
		for (int i = 0; i < nWidth * nHeight; i++)
		{
			data[i] *= data[i];
		}
		//mean and cov compute
		F_MeanCovMapCalculate(data, nWidth, nHeight, corrIP, radius);
		for (int i = 0; i < nWidth * nHeight; i++)
		{
			varIP[i] = corrIP[i] - meanIP[i] * meanIP[i];
		}
		for (int i = 0; i < nWidth * nHeight; i++)
		{
			a[i] = varIP[i] / (varIP[i] + delta);
			b[i] = meanIP[i] - a[i] * meanIP[i];
		}
		//mean and cov compute
		F_MeanCovMapCalculate(a, nWidth, nHeight, meanIP, radius);
		F_MeanCovMapCalculate(b, nWidth, nHeight, corrIP, radius);
		for (int i = 0; i < nWidth * nHeight; i++)
		{
			srcData[i] = (unsigned char)(CLIP3((meanIP[i] * srcData[i] / 255.0f + corrIP[i])*255.0f, 0, 255));
		}
		free(data);
		free(meanIP);
		free(corrIP);
		free(varIP);
		free(a);
		free(b);
	}
	else if(channels == 3){
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
				pSrc += channels;
			}
		}
		F_GuidedFilter(yData, nWidth, nHeight,1, radius, delta);
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
				pSrc += channels;
			}
		}
	}
	return 0;
}
int F_LSNFilter(unsigned char* srcData, int nWidth, int nHeight, int channels, int radius, int delta2)
{
	if (srcData == NULL)
	{
		return 0;
	}
	if (channels == 1) {
		int len = sizeof(unsigned long) * nWidth * nHeight;
		unsigned char* dstData = (unsigned char*)malloc(len);
		unsigned long* meanData = (unsigned long*)malloc(len);
		unsigned long* covData = (unsigned long*)malloc(len);
		memset(meanData, 0, len);
		memset(covData, 0, len);
		F_MeanCovMapCalculate(srcData, nWidth, nHeight, meanData, covData, radius);
		float mean = 0, cov = 0, K = 0;
		int i, j, num = 2 * radius + 1;
		num = num * num;
		int gray = 0;
		float det = 0;
		for (int i = 0; i < len; i++)
		{
			det += (float)covData[i] / len;
		}
		for (j = 0; j < nHeight; j++)
		{
			for (i = 0; i < nWidth; i++)
			{
				len = i + j * nWidth;
				mean = (float)meanData[len];
				cov = (float)covData[len];
				cov = cov - mean * mean;
				K = cov / (cov + delta2);
				gray = (int)((1.0 - K) * mean + K * (float)srcData[len]);
				dstData[len] = CLIP3(gray, 0, 255);
			}
		}
		memcpy(srcData, dstData, len);
		free(meanData);
		free(covData);
		free(dstData);
	}
	else if (channels == 3) {
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
				pSrc += channels;
			}
		}
		F_LSNFilter(yData, nWidth, nHeight,1, radius, delta2);
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
				pSrc += channels;
			}
		}
	}
	return 0;
}
int F_AnisotropicFilter(unsigned char* srcData, int width, int height, int channels, int iter, float k, float lambda, int offset)
{
	int ret = 0;
	if (iter <= 0 || k == 0 || lambda == 0 || offset == 0)
		return 1;
	int i, j, pos1, pos2, pos3, pos4, n, pos_src;
	int NI, SI, EI, WI;
	float cN, cS, cE, cW;
	unsigned char* grayData = (unsigned char*)malloc(sizeof(unsigned char) * width* channels * height);
	if (grayData == NULL)
		return 1;
	unsigned char* pSrc = srcData;
	int stride = width * channels;
	float MAP[512];
	float kk = 1.0f / (k * k);
	for (i = -255; i <= 255; i++)
	{
		MAP[i + 255] = exp(-i * i * kk) * lambda * i;
	}
	int r, g, b;
	for (n = 0; n < iter; n++)
	{
		memcpy(grayData, srcData, sizeof(unsigned char) * height * width* channels);
		pSrc = srcData;
		for (j = 0; j < height; j++)
		{
			for (i = 0; i < width; i++)
			{
				pos_src = i*channels + j * stride;
				pos1 = i * channels + CLIP3((j - offset), 0, height - 1) * stride;
				pos2 = i * channels + CLIP3((j + offset), 0, height - 1) * stride;
				pos3 = (CLIP3((i - offset), 0, width - 1) * channels) + j * stride;
				pos4 = (CLIP3((i + offset), 0, width - 1) * channels) + j * stride;
				b = grayData[pos_src];
				NI = grayData[pos1] - b;
				SI = grayData[pos2] - b;
				EI = grayData[pos3] - b;
				WI = grayData[pos4] - b;
				cN = MAP[NI + 255];// opt:exp(-NI*NI / (k * k));
				cS = MAP[SI + 255];
				cE = MAP[EI + 255];
				cW = MAP[WI + 255];
				pSrc[0] = (int)(CLIP3((b + (cN + cS + cE + cW)), 0, 255));

				pos_src = pos_src + 1;
				pos1 = pos1 + 1;
				pos2 = pos2 + 1;
				pos3 = pos3 + 1;
				pos4 = pos4 + 1;
				g = grayData[pos_src];
				NI = grayData[pos1] - g;
				SI = grayData[pos2] - g;
				EI = grayData[pos3] - g;
				WI = grayData[pos4] - g;
				cN = MAP[NI + 255];
				cS = MAP[SI + 255];
				cE = MAP[EI + 255];
				cW = MAP[WI + 255];
				pSrc[1] = (int)(CLIP3((g + (cN + cS + cE + cW)), 0, 255));

				pos_src = pos_src + 1;
				pos1 = pos1 + 1;
				pos2 = pos2 + 1;
				pos3 = pos3 + 1;
				pos4 = pos4 + 1;
				r = grayData[pos_src];
				NI = grayData[pos1] - r;
				SI = grayData[pos2] - r;
				EI = grayData[pos3] - r;
				WI = grayData[pos4] - r;
				cN = MAP[NI + 255];
				cS = MAP[SI + 255];
				cE = MAP[EI + 255];
				cW = MAP[WI + 255];
				pSrc[2] = (int)(CLIP3((r + (cN + cS + cE + cW)), 0, 255));
				pSrc += channels;
			}
		}
	}
	free(grayData);
	return ret;
};//
int F_MeanShiftFilter(unsigned char* srcData, int nWidth, int nHeight, int channels, int radius, int threshold, int maxIter)
{
	if (srcData == NULL)
	{
		return 0;
	}
	if (radius == 0 || threshold == 0)
		return 0;
	if (channels == 1) {
		int len = sizeof(unsigned char) * nWidth * nHeight;
		int i, j;
		int gray = 0, sum = 0, srcGray = 0, count = 0;
		unsigned char* tempData = (unsigned char*)malloc(sizeof(unsigned char) * len);
		memcpy(tempData, srcData, sizeof(unsigned char) * len);
		for (j = 0; j < nHeight; j++)
		{
			for (i = 0; i < nWidth; i++)
			{
				len = i + j * nWidth;
				int nIter = 0, cx = 0, cy = 0, sumx = 0, sumy = 0;
				srcGray = tempData[len];
				cx = i;
				cy = j;

				while (nIter < maxIter)
				{
					sum = 0;
					sumx = 0;
					sumy = 0;
					count = 0;
					for (int y = cy - radius; y <= cy + radius; y++)
					{
						for (int x = cx - radius; x <= cx + radius; x++)
						{
							int px = CLIP3(x, 0, nWidth - 1);
							int py = CLIP3(y, 0, nHeight - 1);
							len = px + py * nWidth;
							gray = tempData[len];
							if (abs(gray - srcGray) < threshold)
							{
								count++;
								sum += gray;
								sumx += x;
								sumy += y;
							}
						}
					}
					if (count == 0)
						break;
					srcGray = sum / count;
					cx = sumx / count;
					cy = sumy / count;
					nIter++;
				}
				srcData[i + j * nWidth] = srcGray;
			}
		}
		free(tempData);
	}
	else if (channels == 3) {
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
				pSrc += channels;
			}
		}
		F_MeanShiftFilter(yData, nWidth, nHeight,1, radius, threshold, maxIter);
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
				pSrc += channels;
			}
		}
		free(yData);
		free(cbData);
		free(crData);
	}
	return 0;
}
int BEEPSHorizontal(unsigned char* srcPtr, int width, int height, unsigned char* outData, double sigma, int c)
{
	unsigned char* F = (unsigned char*)malloc(sizeof(unsigned char) * width * height);
	int *s = (int*)malloc(sizeof(int) * width);
	int *v = (int*)malloc(sizeof(int) * width);
	int pos = 0, X = 0, Y = 0;
	int p = 0;
	memset(F, 0, width * height);
	memset(outData, 0, width * height);
	memset(s, 0, width);
	memset(v, 0, width);
	unsigned char* D = outData;
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			X = width - 1 - x;
			Y = height - 1 - y;
			if (x == 0)
			{
				pos = x + y * width;
				F[pos] = srcPtr[pos];
				s[0] = srcPtr[pos];
				p = X;
				pos = p + Y * width;
				v[p] = srcPtr[pos];
				D[pos] = srcPtr[pos];
			}
			else
			{
				p = x;
				pos = p + y * width;
				s[p] = (int)(10.0 * GetGaussianValue(srcPtr[pos], F[pos - 1], sigma));
				F[pos] = CLIP3((((100 - s[p] * c) * srcPtr[pos] + s[p] * c * F[pos - 1]) / 100), 0, 255);

				p = X;
				pos = p + Y * width;
				v[p] = (int)(10.0 * GetGaussianValue(srcPtr[pos], D[pos + 1], sigma));
				D[pos] = CLIP3((((100 - v[p] * c) * srcPtr[pos] + v[p] * c * D[pos + 1]) / 100), 0, 255);
			}

		}
	}
	for (int i = 0; i < height * width; i++)
	{
		D[i] = CLIP3(((10 * F[i] - (10 - c) * (srcPtr[i]) + 10 * D[i]) / (10 + c)), 0, 255);
	}
	free(F);
	free(s);
	free(v);
	return 0;
}
int BEEPSVertical(unsigned char* srcPtr, int width, int height, unsigned char* outData, double sigma, int c)
{
	unsigned char* F = (unsigned char*)malloc(sizeof(unsigned char) * width * height);
	unsigned char* D = outData;
	int* s = (int*)malloc(sizeof(int) * height);
	int* v = (int*)malloc(sizeof(int) * height);
	int pos = 0, X = 0, Y = 0;
	memset(s, 0, height);
	memset(v, 0, height);
	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			X = width - 1 - x;
			Y = height - 1 - y;
			if (y == 0)
			{
				pos = x + y * width;
				F[pos] = srcPtr[pos];
				s[y] = srcPtr[pos];

				pos = X + Y * width;
				D[pos] = srcPtr[pos];
				v[Y] = srcPtr[pos];
			}
			else
			{
				pos = x + y * width;
				s[y] = (int)(10.0 * GetGaussianValue(srcPtr[pos], F[pos - width], sigma));
				F[pos] = CLIP3((((100 - s[y] * c) * srcPtr[pos] + s[y] * c * F[pos - width]) / 100), 0, 255);


				pos = X + Y * width;
				v[Y] = (int)(10.0 * GetGaussianValue(srcPtr[pos], D[pos + width], sigma));
				D[pos] = CLIP3((((100 - v[Y] * c) * srcPtr[pos] + v[Y] * c * D[pos + width]) / 100), 0, 255);
			}

		}
	}
	for (int i = 0; i < height*width; i++)
	{
		D[i] = CLIP3(((10 * F[i] - (10 - c) * (srcPtr[i]) + 10 * D[i]) / (10 + c)), 0, 255);
	}
	free(F);
	free(s);
	free(v);
	return 0;
}
int F_BeepsFilter(unsigned char* srcData, int nWidth, int nHeight, int channels, float delta, float delta_s)
{
	if (srcData == NULL)
	{
		return 0;
	}
	if (delta == 0 || delta_s == 0)
		return 0;
	if (channels == 1) {
		float* GMAP = (float*)malloc(sizeof(float) * 256 * 256);
		for (int j = 0; j < 256; j++)
		{
			for (int i = 0; i < 256; i++)
			{
				GMAP[i + j * 256] = GetGaussianValue(i, j, delta);
			}
		}
		delta = delta > 50 ? 50 : delta;
		delta = delta * delta * 2.0f;
		float Lamba = 10.0f * (float)(1 - (sqrt(2.0f * delta_s * delta_s + 1) - 1) / (delta_s * delta_s));
		unsigned char* pSrc = srcData;
		unsigned char* hValue = (unsigned char*)malloc(sizeof(unsigned char) * nWidth * nHeight);
		unsigned char* vValue = (unsigned char*)malloc(sizeof(unsigned char) * nWidth * nHeight);
		unsigned char* dstValue = (unsigned char*)malloc(sizeof(unsigned char) * nWidth * nHeight);
		BEEPSHorizontal(pSrc, nWidth, nHeight, hValue, delta, Lamba);
		BEEPSVertical(hValue, nWidth, nHeight, vValue, delta, Lamba);
		BEEPSVertical(pSrc, nWidth, nHeight, hValue, delta, Lamba);
		BEEPSHorizontal(hValue, nWidth, nHeight, dstValue, delta, Lamba);
		for (int i = 0; i < nWidth * nHeight; i++)
		{
			*pSrc++ = CLIP3(((vValue[i] + dstValue[i]) / 2), 0, 255);
		}
		free(hValue);
		free(vValue);
		free(dstValue);
	}
	else if (channels == 3) {
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
				pSrc += channels;
			}
		}
		F_BeepsFilter(yData, nWidth, nHeight, 1,delta, delta_s);
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
				pSrc += channels;
			}
		}
		free(yData);
		free(cbData);
		free(crData);
	}
	return 0;
}
//Fast mean filter based histagram computation
int f_FastMeanFilter(unsigned char* srcData, int width, int height, int stride, int radius)
{
	int ret = 0;
	if (radius == 0)
		return ret;
	if (radius > MIN2(width, height) / 2)
		radius = (MIN2(width, height) / 2 - 0.5);
	unsigned char* dstData = (unsigned char*)malloc(sizeof(unsigned char) * height * stride);
	memset(dstData, 255, sizeof(unsigned char) * height * stride);
	int unit = 4, t = 0, t1 = 0;
	int i, j, k, len = width * height * unit;
	int block = (radius << 1) + 1;
	int winSize = block * block;
	long sumB = 0, sumG = 0, sumR = 0;
	unsigned char* pSrc = srcData;
	int* temp = (int*)malloc(sizeof(int)* width * unit);
	memset(temp, 0, sizeof(int) * width * unit);
	for (k = -radius; k <= radius; k++)
	{
		for (j = 0; j < width; j++)
		{
			t = j * unit;
			t1 = abs(k) * stride;
			temp[t] += pSrc[t + t1];
			temp[t + 1] += pSrc[t + 1 + t1];
			temp[t + 2] += pSrc[t + 2 + t1];
		}
	}
	for (i = 0; i < height; i++)
	{
		sumB = sumG = sumR = 0;
		for (j = -radius; j <= radius; j++)
		{
			t = abs(j) * unit;
			sumB += temp[t];
			sumG += temp[t + 1];
			sumR += temp[t + 2];
		}
		for (j = 0; j < width; j++)
		{
			t = j * unit + i * stride;
			dstData[t] = (sumB / winSize);
			dstData[t + 1] = (sumG / winSize);
			dstData[t + 2] = (sumR / winSize);
			if (j < width - 1)
			{
				t = abs(j - radius) * unit;
				t1 = (j + radius + 1) % width * unit;
				sumB = sumB - temp[t] + temp[t1];
				sumG = sumG - temp[t + 1] + temp[t1 + 1];
				sumR = sumR - temp[t + 2] + temp[t1 + 2];
			}
		}
		if (i < height - 1)
		{
			for (k = 0; k < width; k++)
			{
				t = k * unit + abs(i - radius) * stride;
				t1 = k * unit + (i + radius + 1) % height * stride;
				temp[k * unit] = temp[k * unit] - pSrc[t] + pSrc[t1];
				temp[k * unit + 1] = temp[k * unit + 1] - pSrc[t + 1] + pSrc[t1 + 1];
				temp[k * unit + 2] = temp[k * unit + 2] - pSrc[t + 2] + pSrc[t1 + 2];
			}
		}
	}
	memcpy(srcData, dstData, sizeof(unsigned char) * height * stride);
	free(dstData);
	free(temp);
	return ret;
};
//std mean filter
int f_MeanFilter(unsigned char *srcData, int width, int height, int stride, int radius)
{
	int ret = 0;
	if (radius == 0)
		return ret;
	int offset = stride - width * 4;
	unsigned char* temp = (unsigned char*)malloc(sizeof(unsigned char) * height * stride);
	memcpy(temp, srcData, sizeof(unsigned char) * height * stride);
	int M = (radius * 2 + 1) * (radius * 2 + 1);
	int sumr = 0, sumg = 0, sumb = 0;
	unsigned char* pSrc = srcData;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			sumr = sumg = sumb = 0;
			for (int n = -radius; n <= radius; n++)
			{
				for (int m = -radius; m <= radius; m++)
				{
					int ny = CLIP3(j + n, 0, height - 1);
					int nx = CLIP3(i + m, 0, width - 1);
					int pos = nx * 4 + ny * stride;
					sumb += temp[pos];
					sumg += temp[pos + 1];
					sumr += temp[pos + 2];
				}
			}
			pSrc[0] = sumb / M;
			pSrc[1] = sumg / M;
			pSrc[2] = sumr / M;
			pSrc += 4;
		}
		pSrc += offset;
	}
	free(temp);
	return ret;
};
int f_Instagram1977LutFilter(unsigned char *srcData, int width, int height, int stride, unsigned char* Map)
{
	unsigned char* pSrc = srcData;
	int offset = stride - width * 4;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			pSrc[0] = Map[pSrc[0] * 4];
			pSrc[1] = Map[pSrc[1] * 4 + 1];
			pSrc[2] = Map[pSrc[2] * 4 + 2];
			pSrc += 4;
		}
		pSrc += offset;
	}
	return 0;
};
int f_Instagram1977Filter(unsigned char* srcData, int width, int height, int stride)
{
	int ret = 0;
	unsigned char MAP_R[256] = { 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 82, 84, 85, 86, 88, 88, 89, 92, 93, 93, 95, 96, 98, 99, 100, 101, 103, 104, 106, 107, 107, 108, 110, 111, 112, 114, 115, 116, 117, 118, 119, 121, 122, 123, 125, 126, 126, 128, 130, 131, 132, 133, 134, 135, 136, 137, 138, 140, 141, 143, 144, 145, 145, 147, 148, 149, 150, 151, 152, 154, 154, 155, 156, 159, 159, 161, 162, 163, 164, 164, 165, 166, 167, 168, 169, 170, 171, 173, 173, 174, 176, 177, 178, 178, 180, 181, 182, 183, 183, 184, 185, 186, 187, 188, 189, 191, 192, 192, 193, 195, 196, 196, 198, 199, 200, 201, 202, 203, 204, 206, 206, 207, 207, 208, 210, 210, 211, 213, 214, 215, 215, 217, 218, 219, 220, 221, 222, 223, 223, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225, 225 };
	unsigned char MAP_G[256] = { 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 58, 59, 60, 61, 63, 64, 65, 66, 67, 69, 70, 71, 73, 73, 74, 75, 76, 78, 79, 80, 81, 83, 84, 85, 86, 87, 88, 89, 90, 91, 94, 95, 96, 96, 97, 98, 99, 100, 102, 104, 105, 106, 107, 108, 108, 109, 110, 112, 113, 114, 115, 116, 117, 118, 119, 120, 122, 123, 124, 125, 126, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 142, 143, 143, 144, 144, 145, 146, 147, 148, 149, 150, 152, 153, 154, 155, 155, 156, 157, 158, 159, 160, 162, 163, 164, 165, 166, 167, 167, 167, 168, 169, 170, 172, 173, 174, 175, 176, 177, 178, 179, 179, 181, 182, 183, 183, 184, 185, 186, 187, 188, 189, 190, 191, 191, 192, 193, 194, 195, 196, 197, 197, 198, 199, 200, 202, 202, 202, 203, 204, 205, 206, 207, 208, 209, 209, 210, 212, 213, 214, 214, 214, 215, 216, 217, 217, 218, 219, 221, 222, 223, 223, 224, 225, 226, 226, 226, 227, 228, 229, 229, 231, 232, 233, 234, 234, 235, 236, 237, 237, 237, 238, 238, 239, 241, 242, 243, 243, 244, 245, 245, 246, 247, 248, 248, 249, 249, 251, 251, 252, 253, 253 };
	unsigned char MAP_B[256] = { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 66, 68, 69, 70, 71, 73, 75, 76, 77, 78, 79, 80, 82, 84, 85, 86, 87, 89, 90, 92, 93, 95, 96, 97, 98, 99, 102, 103, 104, 105, 106, 107, 109, 110, 111, 113, 115, 116, 117, 118, 119, 120, 122, 123, 124, 125, 126, 127, 129, 130, 131, 133, 134, 136, 137, 138, 139, 139, 140, 142, 143, 144, 145, 146, 147, 149, 150, 151, 152, 153, 154, 156, 157, 158, 159, 160, 160, 161, 163, 164, 165, 166, 167, 169, 170, 171, 172, 173, 174, 176, 177, 178, 179, 180, 181, 183, 184, 185, 186, 187, 187, 189, 190, 191, 192, 193, 194, 196, 197, 198, 199, 200, 201, 203, 204, 205, 205, 206, 207, 208, 210, 211, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212, 212 };
	unsigned char* pSrc = srcData;
	int offset = stride - width * 4;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			pSrc[0] = MAP_B[pSrc[0]];
			pSrc[1] = MAP_G[pSrc[1]];
			pSrc[2] = MAP_R[pSrc[2]];
			pSrc += 4;
		}
		pSrc += offset;
	}
	return ret;
};
int f_ConcaveMirrorFilter(unsigned char* srcData, int width, int height, int stride, int x, int y, int k)
{
	x = CLIP3(x, 0, width - 1);
	y = CLIP3(y, 0, height - 1);
	k = MAX2(k, 0);
	int radius = 0;
	float theta = 0;
	int tX = 0;
	int tY = 0;
	int mapX = 0;
	int mapY = 0;
	int mapR = 0;
	unsigned char* pSrc = srcData;
	unsigned char* tempData = (unsigned char*)malloc(sizeof(unsigned char) * height * stride);
	memcpy(tempData, srcData, sizeof(unsigned char) * height * stride);
	int offset = stride - width * 4;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			tX = i - x;
			tY = j - y;
			radius = k;;
			float distance = (i - x) * (i - x) + (j - y) * (j - y);
			float dis = sqrt(distance);
			if (distance <= k * k && distance > 0)
			{
				mapX = floor(dis * (i - x) / k + x);
				mapY = floor(dis * (j - y) / k + y);
				pSrc[0] = tempData[mapX * 4 + mapY * stride];
				pSrc[1] = tempData[mapX * 4 + mapY * stride + 1];
				pSrc[2] = tempData[mapX * 4 + mapY * stride + 2];
			}
			pSrc += 4;
		}
		pSrc += offset;
	}
	free(tempData);
	return 0;
}
int f_ConvexMirrorFilter(unsigned char* srcData, int width, int height, int stride, int x, int y, int k)
{
	x = CLIP3(x, 0, width - 1);
	y = CLIP3(y, 0, height - 1);
	k = MAX2(k, 0);
	int radius = 0;
	float theta = 0;
	int tX = 0;
	int tY = 0;
	int mapX = 0;
	int mapY = 0;
	int mapR = 0;
	unsigned char* pSrc = srcData;
	unsigned char* tempData = (unsigned char*)malloc(sizeof(unsigned char) * height * stride);
	memcpy(tempData, srcData, sizeof(unsigned char) * height * stride);
	int offset = stride - width * 4;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			tX = i - x;
			tY = j - y;
			theta = atan2((float)tY, (float)tX);
			radius = (int)sqrt((float)(tX * tX + tY * tY));
			mapR = (int)(sqrt((float)radius * k));
			mapX = CLIP3(x + (int)(mapR * cos(theta)), 0, width - 1);
			mapY = CLIP3(y + (int)(mapR * sin(theta)), 0, height - 1);
			pSrc[0] = tempData[mapX * 4 + mapY * stride];
			pSrc[1] = tempData[mapX * 4 + mapY * stride + 1];
			pSrc[2] = tempData[mapX * 4 + mapY * stride + 2];
			pSrc += 4;
		}
		pSrc += offset;
	}
	free(tempData);
	return 0;
}
