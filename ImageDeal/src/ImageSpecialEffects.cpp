#include"ImageSpecialEffects.h"
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<iostream>
#include<vector>
int F_IDW(unsigned char* srcData, int width, int height, int channels, int inputPoints[], int outputPoints[], int pointNum)
{
	int ret = 0;
	int stride = channels * width;
	unsigned char* pSrc = srcData;
	int aa, bb, cc, dd, pos, pos1, xx, yy;
	double r1, r2, r3;
	unsigned char* pSrcL1;
	unsigned char* pSrcL2;
	unsigned char* tempData = (unsigned char*)malloc(sizeof(unsigned char) * height * stride);
	memcpy(tempData, srcData, sizeof(unsigned char) * height * stride);
	unsigned char* p = tempData;
	double w = 0, x_in, y_in, sumw, v;
	double u = 1;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			x_in = 0, y_in = 0;
			//F function compute
			sumw = 0;
			for (int k = 0; k < pointNum; k++)
			{
				int t = (k << 1);
				v = 1.0 / (pow((double)((i - inputPoints[t]) * (i - inputPoints[t]) + (j - inputPoints[t + 1]) * (j - inputPoints[t + 1])), u));
				sumw += v;
				if (i == inputPoints[t] && j == inputPoints[t + 1])
					w = 1.0;
				else
					w = v;
				x_in += w * (outputPoints[t] + i - inputPoints[t]);
				y_in += w * (outputPoints[t + 1] + j - inputPoints[t + 1]);
			}
			w = 1.0 / sumw;
			x_in = x_in * w;
			y_in = y_in * w;
			//interpolation
			x_in = CLIP3(x_in, 0, width - 1);
			y_in = CLIP3(y_in, 0, height - 1);
			xx = (int)x_in;
			yy = (int)y_in;
			pSrcL1 = p + yy * stride;
			pSrcL2 = p + CLIP3((yy + 1), 0, height - 1) * stride;
			pos = (xx * channels);
			for (int n = 0; n < channels; n++) {
				aa = pSrcL1[pos + n];
				bb = pSrcL1[pos + channels + n];
				cc = pSrcL2[pos + n];
				dd = pSrcL2[pos + channels + n];
				r1 = aa + (bb - aa) * (x_in - xx);
				r2 = cc + (dd - cc) * (x_in - xx);
				r3 = r1 + (r2 - r1) * (y_in - yy);
				pSrc[n] = (unsigned char)(CLIP3(r3, 0, 255));//B
			}
			pSrc += channels;
		}
	}
	free(tempData);
	return ret;
};
float lineGetU(int x, int y, int px, int py, int qx, int qy, float len)
{
	int xpx = x - px;
	int xpy = y - py;
	int qpx = qx - px;
	int qpy = qy - py;
	return (float)(xpx * qpx + xpy * qpy) / (len * len);
}
float lineGetV(int x, int y, int px, int py, int qx, int qy, float len)
{
	int xpx = x - px;
	int xpy = y - py;
	int qpx = qx - px;
	int qpy = qy - py;
	int per_qpx = qpy;
	int per_qpy = -qpx;
	return (float)(xpx * per_qpx + xpy * per_qpy) / len;
}
void getPointX(float u, float v, int px, int py, int qx, int qy, float len_dst, float *outPointX, float* outPointY)
{
	int qpx = qx - px;
	int qpy = qy - py;
	int per_qpx = qpy;
	int per_qpy = -qpx;
	*outPointX = px + u * (qx - px) + (v * per_qpx) / len_dst;
	*outPointY = py + u * (qy - py) + (v * per_qpy) / len_dst;
}
float getWeight(float x, float y, float px, float py, float qx, float qy, float len, float a, float b, float p, float getU, float getV)
{
	float d = 0;
	float u = getU;//lineGetU(x, y, px, py, qx, qy, len);
	if (u > 1.0f)
		d = sqrt((x - qx) * (x - qx) + (y - qy) * (y - qy));
	else if (u < 0)
		d = sqrt((x - px) * (x - px) + (y - py) * (y - py));
	else
		d = fabs(getV);//abs(lineGetV(x, y, px, py,  qx, qy, len));
	return pow(pow(len, p) / (a + d), b);
}
int F_FBIM(unsigned char* srcData, int width, int height, int channels, int inputlinePoints[], int outputlinePoints[], int lineNum)
{
	int ret = 0;
	int *line_p_x = (int*)malloc(sizeof(int) * lineNum);
	int *line_p_y = (int*)malloc(sizeof(int) * lineNum);
	int *line_q_x = (int*)malloc(sizeof(int) * lineNum);
	int *line_q_y = (int*)malloc(sizeof(int) * lineNum);
	int *line_c_p_x = (int*)malloc(sizeof(int) * lineNum);
	int *line_c_p_y = (int*)malloc(sizeof(int) * lineNum);
	int *line_c_q_x = (int*)malloc(sizeof(int) * lineNum);
	int *line_c_q_y = (int*)malloc(sizeof(int) * lineNum);
	int index = 0;
	for (int i = 0; i < lineNum; i++)
	{
		line_p_x[i] = inputlinePoints[2 * index];
		line_p_y[i] = inputlinePoints[2 * index + 1];
		line_q_x[i] = inputlinePoints[2 * (index + 1)];
		line_q_y[i] = inputlinePoints[2 * (index + 1) + 1];

		line_c_p_x[i] = outputlinePoints[2 * index];
		line_c_p_y[i] = outputlinePoints[2 * index + 1];
		line_c_q_x[i] = outputlinePoints[2 * (index + 1)];
		line_c_q_y[i] = outputlinePoints[2 * (index + 1) + 1];
		index += 2;
	}
	int i, j, pos, nx, ny, k;
	int stride = channels * width;
	unsigned char*  tempData = (unsigned char*)malloc(sizeof(unsigned char) *  height * stride);
	memcpy(tempData, srcData, sizeof(unsigned char) *  height * stride);
	float leftXSum_x, leftXSum_y, leftWeightSum, nu, nv, srcWeight, left_src_x, left_src_y, srcPx, srcPy, len_src, len_dst;
	float t1, t2;
	int pos_x_1, pos_y_1;
	float *MAP_Len_src = (float*)malloc(sizeof(float) * lineNum);
	float* MAP_Len_dst = (float*)malloc(sizeof(float) * lineNum);
	for (k = 0; k < lineNum; k++)
	{
		MAP_Len_src[k] = sqrt((float)(line_c_p_x[k] - line_c_q_x[k]) * (line_c_p_x[k] - line_c_q_x[k]) + (line_c_p_y[k] - line_c_q_y[k]) * (line_c_p_y[k] - line_c_q_y[k]));
		MAP_Len_dst[k] = sqrt((float)(line_p_x[k] - line_q_x[k]) * (line_p_x[k] - line_q_x[k]) + (line_p_y[k] - line_q_y[k]) * (line_p_y[k] - line_q_y[k]));
	}
	int* MAP_S = (int*)malloc(sizeof(int) * height);
	for (k = 0; k < height; k++)
	{
		MAP_S[k] = k * stride;
	}
	float kw = 0, DX, DY;
	for (j = 0; j < height; j++)
	{
		for (i = 0; i < width; i++)
		{
			leftXSum_x = 0.0f;
			leftXSum_y = 0.0f;
			leftWeightSum = 0.0f;
			for (k = 0; k < lineNum; k++)
			{
				len_src = MAP_Len_src[k];
				len_dst = MAP_Len_dst[k];
				nu = lineGetU(i, j, line_c_p_x[k], line_c_p_y[k], line_c_q_x[k], line_c_q_y[k], len_src);
				nv = lineGetV(i, j, line_c_p_x[k], line_c_p_y[k], line_c_q_x[k], line_c_q_y[k], len_src);
				getPointX(nu, nv, line_p_x[k], line_p_y[k], line_q_x[k], line_q_y[k], len_dst, &srcPx, &srcPy);
				srcWeight = getWeight((float)i, (float)j, (float)line_c_p_x[k], (float)line_c_p_y[k], (float)line_c_q_x[k], (float)line_c_q_y[k], (float)len_src, 1.0f, 2.0f, 0.0f, nu, nv);
				leftXSum_x = leftXSum_x + srcPx * srcWeight;
				leftXSum_y = leftXSum_y + srcPy * srcWeight;
				leftWeightSum += srcWeight;
			}
			kw = (1.0f / leftWeightSum);
			left_src_x = leftXSum_x * kw;
			left_src_y = leftXSum_y * kw;
			left_src_x = CLIP3(left_src_x, 0, width - 2);
			left_src_y = CLIP3(left_src_y, 0, height - 2);
			nx = (int)left_src_x;
			ny = (int)left_src_y;
			pos = i*channels + MAP_S[j];

			DX = left_src_x - nx;
			DY = left_src_y - ny;
			pos_x_1 = nx + 1;
			pos_y_1 = ny + 1;

			t1 = (1.f - DX) * tempData[(nx * channels) + MAP_S[ny]] + DX * tempData[(pos_x_1 * channels) + MAP_S[ny]];
			t2 = (1.f - DX) * tempData[(nx * channels) + MAP_S[pos_y_1]] + DX * tempData[(pos_x_1 * channels) + MAP_S[pos_y_1]];
			t1 = (1.f - DY) * t1 + DY * t2;

			srcData[pos] = (unsigned char)CLIP3(t1, 0, 255);
			t1 = (1.f - DX) * tempData[(nx * channels) + 1 + MAP_S[ny]] + DX * tempData[(pos_x_1 * channels) + 1 + MAP_S[ny]];
			t2 = (1.f - DX) * tempData[(nx * channels) + 1 + MAP_S[pos_y_1]] + DX * tempData[(pos_x_1 * channels) + 1 + MAP_S[pos_y_1]];
			t1 = (1.f - DY) * t1 + DY * t2;
			srcData[pos + 1] = (unsigned char)CLIP3(t1, 0, 255);
			t1 = (1.f - DX) * tempData[(nx * channels) + 2 + MAP_S[ny]] + DX * tempData[(pos_x_1 * channels) + 2 + MAP_S[ny]];
			t2 = (1.f - DX) * tempData[(nx * channels) + 2 + MAP_S[pos_y_1]] + DX * tempData[(pos_x_1 * channels) + 2 + MAP_S[pos_y_1]];
			t1 = (1.f - DY) * t1 + DY * t2;
			srcData[pos + 2] = (unsigned char)CLIP3(t1, 0, 255);
		}
	}
	free(MAP_S);
	free(tempData);
	free(line_p_x);
	free(line_p_y);
	free(line_q_x);
	free(line_q_y);
	free(line_c_p_x);
	free(line_c_p_y);
	free(line_c_q_x);
	free(line_c_q_y);
	free(MAP_Len_src);
	free(MAP_Len_dst);
	return ret;
};

struct PointD
{
	double x;
	double y;
};
int f_AffineTransform(unsigned char* srcData, int width, int height, int stride, float H[6], unsigned char* dstData, int dWidth, int dHeight, int dStride)
{
	int ret = 0;
	unsigned char* pSrc = dstData;
	int tx, ty, pos;
	int offset[2];
	offset[0] = ((dWidth / 2.0) - (H[0] * (width / 2.0) + H[1] * (height / 2.0) + H[2]) + 0.5);
	offset[1] = ((dHeight / 2.0) - (H[3] * (width / 2.0) + H[4] * (height / 2.0) + H[5]) + 0.5);
	H[2] += offset[0];
	H[5] += offset[1];
	pSrc = srcData;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			tx = CLIP3(((H[0] * (float)(i)+H[1] * (float)(j)+H[2])), 0, dWidth - 1);
			ty = CLIP3(((H[3] * (float)(i)+H[4] * (float)(j)+H[5])), 0, dHeight - 1);
			pos = (tx << 2) + ty * dStride;
			dstData[pos] = pSrc[0];
			dstData[pos + 1] = pSrc[1];
			dstData[pos + 2] = pSrc[2];
			dstData[pos + 3] = 255;
			pSrc += 4;
		}
	}
	return ret;
};
int f_PerspectiveTransform(unsigned char* srcData, int width, int height, int stride, float H[6])
{
	int ret = 0;
	int tx, ty, pos;
	unsigned char* tempData = (unsigned char*)malloc(sizeof(unsigned char) * height * stride);
	memset(tempData, 0, sizeof(unsigned char) * height * stride);
	unsigned char* pSrc = tempData;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			float w = H[6] * i / width + H[7] * j / height + 1.0f;
			tx = (H[0] * i / width + H[1] * j / height + H[2] / width) / w * width;
			ty = (H[3] * i / width + H[4] * j / height + H[5] / height) / w * height;
			tx = CLIP3(tx, 0, width - 1);
			ty = CLIP3(ty, 0, height - 1);
			if (tx == width - 1 || ty == height - 1 || tx == 0 || ty == 0)
			{
				pSrc[0] = 0;
				pSrc[1] = 0;
				pSrc[2] = 0;
			}
			else
			{
				pos = (tx << 2) + ty * stride;
				pSrc[0] = srcData[pos];
				pSrc[1] = srcData[pos + 1];
				pSrc[2] = srcData[pos + 2];
			}
			pSrc[3] = 255;
			pSrc += 4;
		}
	}
	memcpy(srcData, tempData, sizeof(unsigned char) * height * stride);
	free(tempData);
	return ret;
};
static double bilinear_interp(double x, double y, double v11, double v12,
	double v21, double v22) {
	return (v11 * (1 - y) + v12 * y) * (1 - x) + (v21 * (1 - y) + v22 * y) * x;
}

static double calcArea(const std::vector<PointD> &V) {
	PointD lt, rb;
	lt.x = lt.y = 1e10;
	rb.x = rb.y = -1e10;
	for (std::vector<PointD >::const_iterator i = V.begin(); i != V.end();
		i++) {
		if (i->x < lt.x) lt.x = i->x;
		if (i->x > rb.x) rb.x = i->x;
		if (i->y < lt.y) lt.y = i->y;
		if (i->y > rb.y) rb.y = i->y;
	}
	return (rb.x - lt.x) * (rb.y - lt.y);
}
static void calcDelta_rigid(int srcW, int srcH, int tarW, int tarH, double alpha, int gridSize, int nPoint, int preScale, double *rDx, double *rDy, std::vector<PointD> &oldDotL, std::vector<PointD> &newDotL)
{
	int i, j, k;
	PointD swq, qstar, newP, tmpP;
	double sw;

	double ratio;

	if (preScale) {
		ratio = sqrt(calcArea(newDotL) / calcArea(oldDotL));
		for (i = 0; i < nPoint; i++) {
			newDotL[i].x *= 1 / ratio;
			newDotL[i].y *= 1 / ratio;
		}
	}
	double *w = new double[nPoint];

	if (nPoint < 2) {
		return;
	}
	PointD swp, pstar, curV, curVJ, Pi, PiJ, Qi;
	double miu_r;

	for (i = 0;; i += gridSize) {
		if (i >= tarW && i < tarW + gridSize - 1)
			i = tarW - 1;
		else if (i >= tarW)
			break;
		for (j = 0;; j += gridSize) {
			if (j >= tarH && j < tarH + gridSize - 1)
				j = tarH - 1;
			else if (j >= tarH)
				break;
			sw = 0;
			swp.x = swp.y = 0;
			swq.x = swq.y = 0;
			newP.x = newP.y = 0;
			curV.x = i;
			curV.y = j;
			for (k = 0; k < nPoint; k++) {
				if ((i == oldDotL[k].x) && j == oldDotL[k].y) break;
				if (alpha == 1)
					w[k] = 1 / ((i - oldDotL[k].x) * (i - oldDotL[k].x) +
					(j - oldDotL[k].y) * (j - oldDotL[k].y));
				else
					w[k] = pow((i - oldDotL[k].x) * (i - oldDotL[k].x) +
					(j - oldDotL[k].y) * (j - oldDotL[k].y),
						-alpha);
				sw = sw + w[k];
				swp.x = swp.x + w[k] * oldDotL[k].x;
				swp.y = swp.y + w[k] * oldDotL[k].y;
				swq.x = swq.x + w[k] * newDotL[k].x;
				swq.y = swq.y + w[k] * newDotL[k].y;
			}
			if (k == nPoint) {
				pstar.x = (1 / sw) * swp.x;
				pstar.y = (1 / sw) * swp.y;
				qstar.x = 1 / sw * swq.x;
				qstar.y = 1 / sw * swq.y;
				// Calc miu_r
				double s1 = 0, s2 = 0;
				for (k = 0; k < nPoint; k++) {
					if (i == oldDotL[k].x && j == oldDotL[k].y) continue;
					Pi.x = oldDotL[k].x - pstar.x;
					Pi.y = oldDotL[k].y - pstar.y;
					PiJ.x = -Pi.y, PiJ.y = Pi.x;
					Qi.x = newDotL[k].x - qstar.x;
					Qi.y = newDotL[k].y - qstar.y;
					s1 += w[k] * (Qi.x*Pi.x + Qi.y*Pi.y);
					s2 += w[k] * (Qi.x*PiJ.x + Qi.y*PiJ.y);
				}
				miu_r = sqrt(s1 * s1 + s2 * s2);
				curV.x -= pstar.x;
				curV.y -= pstar.y;

				curVJ.x = -curV.y, curVJ.y = curV.x;

				for (k = 0; k < nPoint; k++) {
					if (i == oldDotL[k].x && j == oldDotL[k].y) continue;
					Pi.x = oldDotL[k].x - pstar.x;
					Pi.y = oldDotL[k].y - pstar.y;
					PiJ.x = -Pi.y, PiJ.y = Pi.x;
					Qi.x = newDotL[k].x - qstar.x;
					Qi.y = newDotL[k].y - qstar.y;
					tmpP.x = (Pi.x*curV.x + Pi.y*curV.y)* Qi.x -
						(PiJ.x*curV.x + PiJ.y*curV.y)* Qi.y;
					tmpP.y = -(Pi.x*curVJ.x + Pi.y*curVJ.y) * Qi.x +
						(PiJ.x*curVJ.x + PiJ.y*curVJ.y) * Qi.y;
					tmpP.x *= w[k] / miu_r;
					tmpP.y *= w[k] / miu_r;
					newP.x += tmpP.x;
					newP.y += tmpP.y;
				}
				newP.x += qstar.x;
				newP.y += qstar.y;
			}
			else {
				newP = newDotL[k];
			}

			if (preScale) {
				rDx[j * tarW + i] = newP.x * ratio - i;
				rDy[j * tarW + i] = newP.y * ratio - j;
			}
			else {
				rDx[j * tarW + i] = newP.x - i;
				rDy[j * tarW + i] = newP.y - j;
			}
		}
	}
	delete[] w;

	if (preScale != 0) {
		for (i = 0; i < nPoint; i++) {
			newDotL[i].x *= ratio;
			newDotL[i].y *= ratio;
		}
	}
}
static void calcDelta_Similarity(int srcW, int srcH, int tarW, int tarH, double alpha, int gridSize, int nPoint, int preScale, double *rDx, double *rDy, std::vector<PointD> &oldDotL, std::vector<PointD> &newDotL)
{
	int i, j, k;

	PointD swq, qstar, newP, tmpP;
	double sw;

	double ratio;

	if (preScale) {
		ratio = sqrt(calcArea(newDotL) / calcArea(oldDotL));
		for (i = 0; i < nPoint; i++) {
			newDotL[i].x *= 1 / ratio;
			newDotL[i].y *= 1 / ratio;
		}
	}
	double *w = new double[nPoint];

	if (nPoint < 2) {
		return;
	}

	PointD swp, pstar, curV, curVJ, Pi, PiJ;
	double miu_s;

	for (i = 0;; i += gridSize) {
		if (i >= tarW && i < tarW + gridSize - 1)
			i = tarW - 1;
		else if (i >= tarW)
			break;
		for (j = 0;; j += gridSize) {
			if (j >= tarH && j < tarH + gridSize - 1)
				j = tarH - 1;
			else if (j >= tarH)
				break;
			sw = 0;
			swp.x = swp.y = 0;
			swq.x = swq.y = 0;
			newP.x = newP.y = 0;
			curV.x = i;
			curV.y = j;
			for (k = 0; k < nPoint; k++) {
				if ((i == oldDotL[k].x) && j == oldDotL[k].y) break;
				w[k] = 1 / ((i - oldDotL[k].x) * (i - oldDotL[k].x) +
					(j - oldDotL[k].y) * (j - oldDotL[k].y));
				sw = sw + w[k];
				swp.x = swp.x + w[k] * oldDotL[k].x;
				swp.y = swp.y + w[k] * oldDotL[k].y;
				swq.x = swq.x + w[k] * newDotL[k].x;
				swq.y = swq.y + w[k] * newDotL[k].y;
			}
			if (k == nPoint) {
				pstar.x = (1 / sw) * swp.x;
				pstar.y = (1 / sw) * swp.y;
				qstar.x = 1 / sw * swq.x;
				qstar.y = 1 / sw * swq.y;
				// Calc miu_s
				miu_s = 0;
				for (k = 0; k < nPoint; k++) {
					if (i == oldDotL[k].x && j == oldDotL[k].y) continue;

					Pi.x = oldDotL[k].x - pstar.x;
					Pi.y = oldDotL[k].y - pstar.y;
					miu_s += w[k] * (Pi.x*Pi.x + Pi.y*Pi.y);
				}

				curV.x -= pstar.x;
				curV.y -= pstar.y;
				curVJ.x = -curV.y, curVJ.y = curV.x;

				for (k = 0; k < nPoint; k++) {
					if (i == oldDotL[k].x && j == oldDotL[k].y) continue;

					Pi.x = oldDotL[k].x - pstar.x;
					Pi.y = oldDotL[k].y - pstar.y;
					PiJ.x = -Pi.y, PiJ.y = Pi.x;
				
					tmpP.x = (Pi.x*curV.x + Pi.y*curV.y) * (newDotL[k].x - qstar.x) -
						(PiJ.x*curV.x + PiJ.y*curV.y) * (newDotL[k].y - qstar.y);
					tmpP.y = -(Pi.x*curVJ.x + Pi.y*curVJ.y) * (newDotL[k].x - qstar.x) +
						(PiJ.x*curVJ.x + PiJ.y*curVJ.y) * (newDotL[k].y - qstar.y);
					tmpP.x *= w[k] / miu_s;
					tmpP.y *= w[k] / miu_s;
					newP.x += tmpP.x;
					newP.y += tmpP.y;
				}
				newP.x += qstar.x;
				newP.y += qstar.y;
			}
			else {
				newP = newDotL[k];
			}

			rDx[j * tarW + i] = newP.x - i;
			rDy[j * tarW + i] = newP.y - j;
		}
	}

	delete[] w;
	if (preScale != 0) {
		for (i = 0; i < nPoint; i++) {
			newDotL[i].x *= ratio;
			newDotL[i].y *= ratio;
		}
	}
}
static int GetNewImg(unsigned char* oriImg, int width, int height, int channels, int gridSize, double* rDx, double* rDy, double transRatio)
{
	int i, j;
	double di, dj;
	double nx, ny;
	int nxi, nyi, nxi1, nyi1;
	double deltaX, deltaY;
	double w, h;
	int ni, nj;
	int tarH = height;
	int tarW = width;
	int stride = width * channels;
	int pos, posa, posb, posc, posd;
	unsigned char* tarImg = (unsigned char*)malloc(sizeof(unsigned char)*channels*width*height);
	for (i = 0; i < tarH; i += gridSize)
		for (j = 0; j < tarW; j += gridSize) {
			ni = i + gridSize, nj = j + gridSize;
			w = h = gridSize;
			if (ni >= tarH) ni = tarH - 1, h = ni - i + 1;
			if (nj >= tarW) nj = tarW - 1, w = nj - j + 1;
			for (di = 0; di < h; di++)
				for (dj = 0; dj < w; dj++) {
					deltaX =
						bilinear_interp(di / h, dj / w, rDx[i * tarW + j], rDx[i * tarW + nj],
							rDx[ni * tarW + j], rDx[ni * tarW + nj]);
					deltaY =
						bilinear_interp(di / h, dj / w, rDy[i * tarW + j], rDy[i * tarW + nj],
							rDy[ni * tarW + j], rDy[ni * tarW + nj]);
					nx = j + dj + deltaX * transRatio;
					ny = i + di + deltaY * transRatio;
					if (nx > width - 1) nx = width - 1;
					if (ny > height - 1) ny = height - 1;
					if (nx < 0) nx = 0;
					if (ny < 0) ny = 0;
					nxi = int(nx);
					nyi = int(ny);
					nxi1 = ceil(nx);
					nyi1 = ceil(ny);
					pos = (int)(i + di) * stride + ((int)(j + dj) *channels);
					posa = nyi * stride + (nxi *channels);
					posb = nyi * stride + (nxi1 *channels);
					posc = nyi1 * stride + (nxi *channels);
					posd = nyi1 * stride + (nxi1 *channels);
					for (int n = 0; n < channels; n++) {
						tarImg[pos + n] = (unsigned char)bilinear_interp(ny - nyi, nx - nxi, oriImg[posa + n], oriImg[posb + n], oriImg[posc + n], oriImg[posd + n]);
					}
				}
		}
	memcpy(oriImg, tarImg, sizeof(unsigned char)*channels*width*height);
	free(tarImg);
	return 0;
};
void F_MLSImageWrapping(unsigned char* oriImg, int width, int height, int channels, int srcPoint[], int dragPoint[], int pointNum, double transRatio, int preScale, int gridSize, int method)
{
	int srcW = width;
	int srcH = height;
	int tarW = width;
	int tarH = height;
	double alpha = 1;
	int len = tarH * tarW;
	std::vector<PointD> oldDotL, newDotL;
	double *rDx = NULL, *rDy = NULL;
	PointD point = { 0 };
	for (int i = 0; i < pointNum; i++)
	{
		point.x = srcPoint[(i << 1)];
		point.y = srcPoint[(i << 1) + 1];
		newDotL.push_back(point);
		point.x = dragPoint[(i << 1)];
		point.y = dragPoint[(i << 1) + 1];
		oldDotL.push_back(point);
	}
	rDx = (double*)malloc(sizeof(double) * len);
	rDy = (double*)malloc(sizeof(double) * len);
	memset(rDx, 0, sizeof(double) * len);
	memset(rDy, 0, sizeof(double) * len);
	if (method != 0)
		calcDelta_Similarity(srcW, srcH, tarW, tarH, alpha, gridSize, pointNum, preScale, rDx, rDy, oldDotL, newDotL);
	else
		calcDelta_rigid(srcW, srcH, tarW, tarH, alpha, gridSize, pointNum, preScale, rDx, rDy, oldDotL, newDotL);
	GetNewImg(oriImg, srcW, srcH, channels, gridSize, rDx, rDy, transRatio);
	if (rDx != NULL)
		free(rDx);
	if (rDy != NULL)
		free(rDy);
};
int F_FaceLift(unsigned char* srcData, int width, int height, int channels, int centerX, int centerY, int rmax, int mx, int my, int strength)
{
	int ret = 0;
	unsigned char* pSrc = srcData;
	int stride = channels * width;
	unsigned char* tempData = (unsigned char*)malloc(sizeof(unsigned char) * height * stride);
	memcpy(tempData, srcData, sizeof(unsigned char) * height * stride);
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			float dis = sqrt((float)(centerX - i) * (centerX - i) + (centerY - j) * (centerY - j));
			if (dis < rmax)
			{
				float xc = (i - centerX) * (i - centerX) + (j - centerY) * (j - centerY);
				float mc = (mx - centerX) * (mx - centerX) + (my - centerY) * (my - centerY);
				float tx = (rmax * rmax - xc);
				float d = tx / (tx + mc * 100.0f / strength);
				d = d * d;
				int px = CLIP3(i - (mx - centerX) * d * (1.0f - dis / rmax), 0, width - 1);
				int py = CLIP3(j - (my - centerY) * d * (1.0f - dis / rmax), 0, height - 1);
				int pos = (px *channels) + py * stride;
				pSrc[0] = tempData[pos];
				pSrc[1] = tempData[pos + 1];
				pSrc[2] = tempData[pos + 2];
			}
			pSrc += channels;
		}
	}
	free(tempData);
	return ret;
};
int F_BigEye(unsigned char* srcData, int width, int height, int channels, int cenX, int cenY, int radius, int intensity)
{
	int ret = 0;
	int nx = CLIP3(cenX - radius, 0, width - 1);
	int ny = CLIP3(cenY - radius, 0, height - 1);
	int nw = CLIP3(cenX + radius, 0, width - 1);
	int nh = CLIP3(cenY + radius, 0, height - 1);
	int D = radius * radius;
	float k0 = intensity / 100.0f;
	int stride = channels * width;
	unsigned char* tempData = (unsigned char*)malloc(sizeof(unsigned char) * height * stride);
	memcpy(tempData, srcData, sizeof(unsigned char) * height * stride);
	for (int j = ny; j < nh; j++)
	{
		for (int i = nx; i < nw; i++)
		{
			float dis = (i - cenX) * (i - cenX) + (j - cenY) * (j - cenY);
			if (dis < D)
			{
				float k = 1.0f - k0 * (1.0f - dis / D);

				//////////////////可以使用双线性插值，效果更好///////////////////////////////////////
				int px = CLIP3((i - cenX) * k + cenX, 0, width - 1);
				int py = CLIP3((j - cenY) * k + cenY, 0, height - 1);
				/////////////////////////////////////////////////////////////////////////////////////
				int pos_new = px * channels + py * stride;
				int pos = i * channels + j * stride;
				srcData[pos] = tempData[pos_new];
				srcData[pos + 1] = tempData[pos_new + 1];
				srcData[pos + 2] = tempData[pos_new + 2];
			}
		}
	}
	free(tempData);
	return ret;
};