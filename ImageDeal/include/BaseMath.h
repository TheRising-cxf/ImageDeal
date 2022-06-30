#pragma once
#include "common.h"
void getWHFromHMatrix(int width, int height, float H[6], int wh[2]);
void F_AffinetransformMetrixCompute(float x1, float y1, float x2, float y2, float x3, float y3, float tx1, float ty1, float tx2, float ty2, float tx3, float ty3, float hMatrix[6]);