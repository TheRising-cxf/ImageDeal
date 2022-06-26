#pragma once
#include<string.h>
#include<stdlib.h>
#include "ImageSharpen.h"
int F_LaplaceSharpen(unsigned char* srcData, int width, int height, int channels, int mode)
{
	int ret = 0;
	int stride = channels * width;
	unsigned char* tempData = (unsigned char*)malloc(sizeof(unsigned char) * height * stride);
	int pos;
	memcpy(tempData, srcData, sizeof(unsigned char) * height * stride);
	if (mode == 0)
	{
		for (int j = 1; j < height - 1; j++)
		{
			for (int i = 1; i < width - 1; i++)
			{
				pos = i * channels + j * stride;
				srcData[pos] = CLIP3(tempData[pos] + (tempData[pos] * 4 - tempData[pos - stride] - tempData[pos - channels] - tempData[pos + channels] - tempData[pos + stride]), 0, 255);
				pos++;
				srcData[pos] = CLIP3(tempData[pos] + (tempData[pos] * 4 - tempData[pos - stride] - tempData[pos - channels] - tempData[pos + channels] - tempData[pos + stride]), 0, 255);
				pos++;
				srcData[pos] = CLIP3(tempData[pos] + (tempData[pos] * 4 - tempData[pos - stride] - tempData[pos - channels] - tempData[pos + channels] - tempData[pos + stride]), 0, 255);
			}
		}
	}
	else
	{
		for (int j = 1; j < height - 1; j++)
		{
			for (int i = 1; i < width - 1; i++)
			{
				pos = i * channels + j * stride;
				srcData[pos] = CLIP3(tempData[pos] + (tempData[pos] * 8 - tempData[pos - stride] - tempData[pos - channels] - tempData[pos + channels] - tempData[pos + stride] - tempData[pos - channels - stride] - tempData[pos +channels - stride] - tempData[pos - channels + stride] - tempData[pos + channels + stride]), 0, 255);
				pos++;																								
				srcData[pos] = CLIP3(tempData[pos] + (tempData[pos] * 8 - tempData[pos - stride] - tempData[pos - channels] - tempData[pos + channels] - tempData[pos + stride] - tempData[pos - channels - stride] - tempData[pos +channels - stride] - tempData[pos - channels + stride] - tempData[pos + channels + stride]), 0, 255);
				pos++;																													
				srcData[pos] = CLIP3(tempData[pos] + (tempData[pos] * 8 - tempData[pos - stride] - tempData[pos - channels] - tempData[pos + channels] - tempData[pos + stride] - tempData[pos - channels - stride] - tempData[pos +channels - stride] - tempData[pos - channels + stride] - tempData[pos + channels + stride]), 0, 255);
			}
		}
	}
	free(tempData);
	return ret;
};