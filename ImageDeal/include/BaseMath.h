#pragma once
#include "common.h"
#include<algorithm>
#include<iostream>
#include<vector>
#define MAXNUM  100000 
class Point {
public:
	float x, y, z;
	Point() {}
	Point(float x1, float y1, float z1 = 0) { x = x1; y = y1; z = z1; }
	void init(float x1, float y1, float z1 = 0) { x = x1; y = y1; z = z1; }
	bool operator ==(const Point &p)
	{
		return fabs(x - p.x) <= 1e-6&&fabs(y - p.y) <= 1e-6&&fabs(z - p.z) <= 1e-6;
	}
	Point operator +(const Point& p1)
	{
		return Point(p1.x + x, p1.y + y, p1.z + z);
	}
	Point operator -(const Point& p1)
	{
		return Point(x - p1.x, y - p1.y, z - p1.z);
	}
};
struct Line {
	Point p1, p2;
	Line() {}
	Line(Point p1, Point p2) { this->p1 = p1; this->p2 = p2; }
	void init(Point p1, Point p2) { this->p1 = p1; this->p2 = p2; }
	bool operator ==(const Line &l)
	{
		return p1 == l.p1&& p2 == l.p2;
	}
};
bool IsRightPoint(Point pt, Line line);
bool IsOnLine(Point pt, Line line);
struct Triangle {
	Point p1, p2, p3;
	Point p[3];
	Line l1, l2, l3;
	Line l[3];
	Triangle() {}
	Triangle(Point a, Point b, Point c) { init(a, b, c); }
	bool isInTriangle(Point a)
	{
		bool r1 = IsRightPoint(a, l1);
		bool r2 = IsRightPoint(a, l2);
		bool r3 = IsRightPoint(a, l3);
		if (r1 == r2 && r2 == r3)
			return true;
		return false;
	}
	int isOnTriangle(Point a)
	{
		bool r1 = IsOnLine(a, l1);
		bool r2 = IsOnLine(a, l2);
		bool r3 = IsOnLine(a, l3);
		if (r1)
			return 1;
		if (r2)
			return 2;
		if (r3)
			return 3;
		return 0;
	}

	void init(Point a, Point b, Point c) {
		p1 = a;
		p2 = b;
		p3 = c;
		p[0] = a;
		p[1] = b;
		p[2] = c;
		l1 = Line(a, b);
		l2 = Line(b, c);
		l3 = Line(c, a);
		l[0] = l1;
		l[1] = l2;
		l[2] = l3;
	}
	int  containsLine(Line l)
	{
		if ((l.p1 == p1 && l.p2 == p2) || (l.p1 == p2 && l.p2 == p1))
			return 1;
		if ((l.p1 == p2 && l.p2 == p3) || (l.p1 == p3 && l.p2 == p2))
			return 2;
		if ((l.p1 == p1 && l.p2 == p3) || (l.p1 == p3 && l.p2 == p1))
			return 3;
		return 0;
	}
	bool operator ==(const Triangle& tri)
	{
		if (p1 == tri.p1 && p2 == tri.p2&& p3 == tri.p3)
			return true;
		return false;
	}
};
struct Circle {
	double radius;
	Point center;
	Circle() {}
	Circle(Point cent, double r) { center = cent; radius = r; }
	static Circle genTriCircle(Triangle tri) {
		Point p1 = tri.p1;
		Point p2 = tri.p2;
		Point p3 = tri.p3;
		double  x1, x2, x3, y1, y2, y3;
		x1 = p1.x;
		x2 = p2.x;
		x3 = p3.x;
		y1 = p1.y;
		y2 = p2.y;
		y3 = p3.y;
		double a = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
		double b = sqrt((x1 - x3)*(x1 - x3) + (y1 - y3)*(y1 - y3));
		double c = sqrt((x2 - x3)*(x2 - x3) + (y2 - y3)*(y2 - y3));
		double p = (a + b + c) / 2;
		double S = sqrt(p*(p - a)*(p - b)*(p - c));

		double t1 = x1 * x1 + y1 * y1;
		double t2 = x2 * x2 + y2 * y2;
		double t3 = x3 * x3 + y3 * y3;
		double temp = x1 * y2 + x2 * y3 + x3 * y1 - x1 * y3 - x2 * y1 - x3 * y2;
		double x = (t2*y3 + t1 * y2 + t3 * y1 - t2 * y1 - t3 * y2 - t1 * y3) / temp / 2;
		double y = (t3*x2 + t2 * x1 + t1 * x3 - t1 * x2 - t2 * x3 - t3 * x1) / temp / 2;
		double r = a * b*c / (4 * S);
		Point cent(x, y);
		return Circle(cent, r);

	}
	bool isInCircle(Point v)
	{
		double dist = sqrt(pow(center.x - v.x, 2) + pow(center.y - v.y, 2));
		return dist < radius;
	}
};
void DivideHell(std::vector<Point>p, std::vector<Triangle>&hullTriangle);
std::vector<Point> CalcConvexHull(std::vector<Point>p, std::vector<Point>& hullPoints);
void GetDelaunay(std::vector<Point> pts, std::vector<Triangle>& hulltins);
void GetWHFromHMatrix(int width, int height, float H[6], int wh[2]);
void F_AffinetransformMetrixCompute(float x1, float y1, float x2, float y2, float x3, float y3, float tx1, float ty1, float tx2, float ty2, float tx3, float ty3, float hMatrix[6]);