#include"BaseMath.h"
void GetWHFromHMatrix(int width, int height, float H[6], int wh[2])
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
float Cross2D(Point a, Point b)
{
	return a.x*b.y - b.x*a.y;
}
bool IsRightPoint(Point pt, Line line)
{
	Point v1(line.p2.x - line.p1.x, line.p2.y - line.p1.y);//p1p2
	Point v2(line.p1.x - pt.x, line.p1.y - pt.y);// pp1
	float tmp = Cross2D(v1, v2);
	if (tmp > 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool IsOnLine(Point pt, Line line)
{
	Point v1(line.p2.x - line.p1.x, line.p2.y - line.p1.y);
	Point v2(line.p1.x - pt.x, line.p1.y - pt.y);
	float tmp = Cross2D(v1, v2);
	float minx, miny;
	float maxx, maxy;
	if (line.p1.x < line.p2.x) {
		minx = line.p1.x;
		maxx = line.p2.x;
	}
	else {
		minx = line.p2.x;
		maxx = line.p1.x;
	}

	if (line.p1.y < line.p2.y) {
		miny = line.p1.y;
		maxy = line.p2.y;
	}
	else {
		miny = line.p2.y;
		maxy = line.p1.y;
	}
	if (fabs(tmp) <= 1e-5 && pt.x > minx&&pt.x < maxx&& pt.y>miny && pt.y < maxy)
	{
		return true;
	}
	else
	{
		return false;
	}
};
static float angle3D(Point a, Point b, Point c) {
	Point v1(a - b);
	Point v2(c - b);
	float x1 = v1.x;
	float y1 = v1.y;
	float x2 = v2.x;
	float y2 = v2.y;
	float n = (x1*x2 + y1 * y2);
	float m = sqrtf(x1*x1 + y1 * y1)*sqrtf(x2*x2 + y2 * y2);
	return acos(n / m) * 180 / 3.1415;
}
std::vector<Point>tmpP;
bool CMP(Point p1, Point p2) {
	//p0排首位
	Point p0 = tmpP[0];
	if (p1.x == p0.x && p1.y == p0.y)return true;
	if (p2.x == p0.x && p2.y == p0.y)return false;

	//计算极角（等于0则赋予一个极大值）
	double angle1 = p1.x == p0.x ? MAXNUM : (float)(p1.y - p0.y) / (p1.x - p0.x);
	double angle2 = p2.x == p0.x ? MAXNUM : (float)(p2.y - p0.y) / (p2.x - p0.x);

	//极角排序
	if (angle1 > angle2)return true;
	else if (angle1 == angle2) {
		if (p1.y < p2.y)return true;
		else return false;
	}
	else return false;
}
double cross(Point p1, Point p2, Point p3)
{
	return -((p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y));
}
std::vector<Point> CalcConvexHull(std::vector<Point>p, std::vector<Point>& hullPoints) {
	hullPoints.clear();
	tmpP.assign(p.begin(), p.end());
	std::vector<Point>usePoints;
	for (int i = 1; i < tmpP.size(); i++) {
		if (tmpP[i].x < tmpP[0].x) {
			Point tmp = tmpP[0];
			tmpP[0] = tmpP[i];
			tmpP[i] = tmp;
		}
	}
	std::sort(tmpP.begin() + 1, tmpP.end(), CMP);
	hullPoints.push_back(tmpP[0]);
	hullPoints.push_back(tmpP[1]);
	tmpP.push_back(tmpP[0]);
	for (int i = 2; i < tmpP.size() - 1; i++) {
		if (cross(hullPoints[hullPoints.size() - 1], tmpP[i], tmpP[i + 1]) >= 0) {
			while (hullPoints.size() > 2) {
				if (cross(hullPoints[hullPoints.size() - 2], hullPoints[hullPoints.size() - 1], tmpP[i]) < 0) {
					usePoints.push_back(hullPoints.back());
					hullPoints.pop_back();
					continue;
				}
				break;
			}
			hullPoints.push_back(tmpP[i]);
		}
		else {
			usePoints.push_back(tmpP[i]);
		}
	}
	return usePoints;
}
void DivideHell(std::vector<Point>p, std::vector<Triangle>&hullTriangle) {
	hullTriangle.clear();
	std::vector<Point>hullPoints;
	CalcConvexHull(p, hullPoints);
	std::vector<Point> hpts;
	for (int i = 0; i < p.size(); i++)
		hpts.push_back(p[i]);
	int index = 0;
	while (hpts.size() > 2)
	{
		int tag = index;
		float minangle = 180;  //每次构成相邻的边，优先找角度最小的
		float maxangle = 0;

		for (int i = index; i < hpts.size(); i++)
		{
			float tri_angle = 180.0;

			if (i == 0) {
				tri_angle = angle3D(hpts.back(), hpts[i], hpts[i + 1]);
			}
			else if (i == hpts.size() - 1) {
				tri_angle = angle3D(hpts[i - 1], hpts[i], hpts[0]);
			}
			else {
				tri_angle = angle3D(hpts[i - 1], hpts[i], hpts[i + 1]);
			}
			if (tri_angle < minangle)
			{
				tag = i;
				minangle = tri_angle;
			}
		}
		int tagb = tag - 1;
		int tage = tag + 1;
		if (tag == 0)
			tagb = hpts.size() - 1;
		if (tag == hpts.size() - 1)
			tage = 0;
		hullTriangle.push_back(Triangle(hpts[tagb], hpts[tag], hpts[tage]));
		hpts.erase(hpts.begin() + tag);
	}

}
void GetDelaunay(std::vector<Point> pts, std::vector<Triangle>& hulltins)
{
	hulltins.clear();
	std::vector<Point>hullPoints;
	std::vector<Point> inPoints = CalcConvexHull(pts, hullPoints);
	DivideHell(hullPoints, hulltins);
	int onLine = 0;
	for (int i = 0; i < inPoints.size(); i++)
	{
		std::vector<Triangle> delTin; //保存要删除的三角形 */
		for (int j = 0; j < hulltins.size(); j++)
		{
			if (hulltins[j].isInTriangle(inPoints[i]) == true)
			{
				delTin.push_back(hulltins[j]);
			}
			if (hulltins[j].isOnTriangle(inPoints[i]) == true) {
				delTin.push_back(hulltins[j]);
				onLine++;
			}
		}
		std::vector<Line> borderLines;//寻找离散点的相邻边 */
		std::vector<Triangle> newTri;//寻找离散点的相邻边 */
		if (delTin.size() == 1)
		{

			Line l1 = delTin[0].l1;
			Line l2 = delTin[0].l2;
			Line l3 = delTin[0].l3;
			borderLines.push_back(l1);
			borderLines.push_back(l2);
			borderLines.push_back(l3);

			for (int j = 0; j < borderLines.size(); j++)
			{
				newTri.push_back(Triangle(borderLines[j].p1, borderLines[j].p2, inPoints[i]));
			}
			std::remove(hulltins.begin(), hulltins.end(), delTin[0]);
			hulltins.pop_back();
		}
		else if (delTin.size() == 2)
		{

			Line l[3];
			l[0] = delTin[0].l1;
			l[1] = delTin[0].l2;
			l[2] = delTin[0].l3;
			int index = 0;
			for (int m = 0; m < 3; m++)
				if (delTin[1].containsLine(l[m]) == 0)
				{
					borderLines.push_back(l[m]);
				}
				else
				{
					index = delTin[1].containsLine(l[m]) - 1;
				}
			for (int m = 0; m < 3; m++)
			{
				if (m != index)
					borderLines.push_back(delTin[1].l[m]);
			}


			for (int j = 0; j < borderLines.size(); j++)
			{
				newTri.push_back(Triangle(borderLines[j].p1, borderLines[j].p2, inPoints[i]));
			}
			std::remove(hulltins.begin(), hulltins.end(), delTin[0]);
			hulltins.pop_back();
			std::remove(hulltins.begin(), hulltins.end(), delTin[1]);
			hulltins.pop_back();
		}

		delTin.clear();
		for (int s = 0; s < newTri.size(); s++) {
			for (int j = 0; j < 3; j++) {
				Line line = newTri[s].l[j];
				for (int m = 0; m < hulltins.size(); m++)
				{
					if (hulltins[m].containsLine(line))
					{
						Circle tinCircle = Circle::genTriCircle(hulltins[m]);
						if (tinCircle.isInCircle(Point(inPoints[i])))
						{
							auto it = std::remove(newTri.begin(), newTri.end(), newTri[s]);
							newTri.pop_back();
							int x = hulltins[m].containsLine(line) - 1;
							for (int k = 0; k < 3; k++)
							{
								if (x != k) {
									newTri.push_back(Triangle(hulltins[m].l[k].p1, hulltins[m].l[k].p2, inPoints[i]));
								}
							}
							std::swap(*(std::begin(hulltins) + m), *(std::end(hulltins) - 1));//等同于 swap(demo[1],demo[4])
							hulltins.pop_back();
						}

					}
				}
			}

		}

		for (int m = 0; m < newTri.size(); m++)
		{
			hulltins.push_back(newTri[m]);
		}
	}
};