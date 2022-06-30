void getWHFromHMatrix(int width, int height, float H[6], int wh[2])
{
	int x0 = (H[0] * 0 + H[1] * 0 + H[2] + 0.5);
	int y0 = (H[3] * 0 + H[4] * 0 + H[5] + 0.5);
	int x1 = (H[0] * (float)(width - 1) + H[1] * (float)(height - 1) + H[2] + 0.5);
	int y1 = (H[3] * (float)(width - 1) + H[4] * (float)(height - 1) + H[5] + 0.5);
	int x2 = (H[0] * (float)(width - 1) + H[1] * 0 + H[2] + 0.5);
	int y2 = (H[3] * (float)(width - 1) + H[4] * 0 + H[5] + 0.5);
	int x3 = (H[0] * 0 + H[1] * (float)(height - 1) + H[2] + 0.5);
	int y3 = (H[3] * 0 + H[4] * (float)(height - 1) + H[5] + 0.5);
	wh[0] = MAX2(x0, MAX2(x1, MAX2(x2, x3))) - MIN2(x0, MIN2(x1, MIN2(x2, x3)));
	wh[1] = MAX2(y0, MAX2(y1, MAX2(y2, y3))) - MIN2(y0, MIN2(y1, MIN2(y2, y3)));
};
void F_AffinetransformMetrixCompute(float x1, float y1, float x2, float y2, float x3, float y3, float tx1, float ty1, float tx2, float ty2, float tx3, float ty3, float hMatrix[6])
{
	//求行列式|A|  
	float detA;
	detA = tx1 * ty2 + tx2 * ty3 + tx3 * ty1 - tx3 * ty2 - tx1 * ty3 - tx2 * ty1;
	// 求伴随矩阵A*  
	float A11, A12, A13, A21, A22, A23, A31, A32, A33;
	A11 = ty2 - ty3;
	A21 = -(ty1 - ty3);
	A31 = ty1 - ty2;
	A12 = -(tx2 - tx3);
	A22 = tx1 - tx3;
	A32 = -(tx1 - tx2);
	A13 = tx2 * ty3 - tx3 * ty2;
	A23 = -(tx1 * ty3 - tx3 * ty1);
	A33 = tx1 * ty2 - tx2 * ty1;
	//求变换矩阵H=A*/|A|  
	//float texMatrix[6]={0};  
	hMatrix[0] = (x1 * A11 + x2 * A21 + x3 * A31) / detA;
	hMatrix[1] = (x1 * A12 + x2 * A22 + x3 * A32) / detA;
	hMatrix[2] = (x1 * A13 + x2 * A23 + x3 * A33) / detA;
	hMatrix[3] = (y1 * A11 + y2 * A21 + y3 * A31) / detA;
	hMatrix[4] = (y1 * A12 + y2 * A22 + y3 * A32) / detA;
	hMatrix[5] = (y1 * A13 + y2 * A23 + y3 * A33) / detA;
};