#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include <math.h>
#include <ostream>
#include <cassert>

#define DOUBLE_EPS 1e-12
#define LDOUBLE_EPS 1e-6
#define LHDOUBLE_EPS 1e-5
#define	EQUALZERO(x)	(fabs((x)) < DOUBLE_EPS)
#define	LEQUALZERO(x) (fabs((x)) < LDOUBLE_EPS)
#define PI	3.141592653589793238

/////////////////////////////////////////////////////////////
// Vector2D : 2D vector
/////////////////////////////////////////////////////////////
class   Vector2D
{
public:
	double x, y;

public:
	Vector2D() { x = 0;	y = 0; }
	// constructions
	Vector2D(double xx, double yy) { x = xx;	y = yy; }
	Vector2D(const Vector2D& v) { x = v.x;	y = v.y; }


	// operator
	double	  length() const { return sqrt(x*x + y * y); }
	double	  length2() const { return x * x + y * y; }
	double	  normalize() { double len = length();	if (!EQUALZERO(len)) { x /= len; y /= len; }	return len; }
	double normalizeStrict() { double len = length(); if (len != 0.0) { x /= len; y /= len; } return len; }
	Vector2D& operator=(const Vector2D& v);
	Vector2D& operator+=(const Vector2D& v);
	Vector2D& operator-=(const Vector2D& v);
	Vector2D& operator*=(double u);
	Vector2D& operator/=(double u);
	double& operator[](int idx) {
		assert(idx < 2);
		switch (idx)
		{
		case 0: return x;
		case 1: return y;
		}
	}
	bool operator == (const Vector2D& right) const
	{
		return x == right.x && y == right.y;
	}
	//Vector2D& operator^=(const Vector2D& v);

	bool	Intersect(Vector2D v1, Vector2D v2, Vector2D v3, Vector2D v4);
	bool	Intersect(Vector2D v1, Vector2D v2);

	friend Vector2D operator+(const Vector2D& lv, const Vector2D& rv);
	friend Vector2D operator-(const Vector2D& lv, const Vector2D& rv);
	friend Vector2D operator*(const double u, const Vector2D& rv);
	friend Vector2D operator*(const Vector2D& lv, const double u);
	friend Vector2D operator/(const Vector2D& lv, const double u);
	friend double   operator*(const Vector2D& lv, const Vector2D& rv);
	friend double operator^(const Vector2D& lv, const Vector2D& rv);
	friend std::ostream& operator<< (std::ostream &output, Vector2D& v);

	short	AtWhere(Vector2D v0, Vector2D v1);
	bool	AtRight(Vector2D v0, Vector2D v1);
	bool	AtLeft(Vector2D v0, Vector2D v1);
	bool	OnLine(Vector2D v0, Vector2D v1);
	double	GetArea(Vector2D v);


};

/////////////////////////////////////////////////////////////
// Vector3D : 3D vector
/////////////////////////////////////////////////////////////

class Vector4D;

class   Vector3D
{
public:
	double x, y, z;

	// constructions
	Vector3D() { x = 0;	y = 0;	z = 0; }
	Vector3D(double xx, double yy, double zz) { x = xx;	y = yy;	z = zz; }
	Vector3D(const Vector3D& v) { x = v.x;	y = v.y;	z = v.z; }
	Vector3D(const Vector4D& v);

	// operator
	double	  length() const { return sqrt(x*x + y * y + z * z); }
	double	  length2() const { return x * x + y * y + z * z; }
	double	  normalize() { double len = length();	if (!EQUALZERO(len)) { x /= len; y /= len; z /= len; }	return len; }
	double normalizeStrict() { double len = length(); if (len != 0.0) { x /= len; y /= len; z /= len; } return len; }

	Vector3D& operator=(const Vector3D& v);
	Vector3D& operator=(const Vector4D& v);
	Vector3D& operator+=(const Vector3D& v);
	Vector3D& operator-=(const Vector3D& v);
	Vector3D& operator*=(double u);
	Vector3D& operator/=(double u);
	Vector3D& operator^=(const Vector3D& v);
	double& operator[](int idx) {
		assert(idx < 3);
		switch (idx)
		{
		case 0: return x;
		case 1: return y;
		case 2: return z;
		}
	}

	bool operator < (const Vector3D& right) const
	{
		return x - right.x < -LHDOUBLE_EPS ||
			fabs(x - right.x) <= LHDOUBLE_EPS && y - right.y < -LHDOUBLE_EPS ||
			fabs(x - right.x) <= LHDOUBLE_EPS && fabs(y - right.y) <= LHDOUBLE_EPS && z - right.z < -LHDOUBLE_EPS;
	}

	bool operator == (const Vector3D& right) const
	{
		return x == right.x && y == right.y && z == right.z;
	}

	friend Vector3D operator+(const Vector3D& lv, const Vector3D& rv);
	friend Vector3D operator-(const Vector3D& lv, const Vector3D& rv);
	friend Vector3D operator*(const double u, const Vector3D& rv);
	friend Vector3D operator*(const Vector3D& lv, const double u);
	friend Vector3D operator/(const Vector3D& lv, const double u);
	friend double   operator*(const Vector3D& lv, const Vector3D& rv);
	friend Vector3D operator^(const Vector3D& lv, const Vector3D& rv);
	friend std::ostream& operator<< (std::ostream &output, const Vector3D& v);

};

class   Vector4D
{
public:
	double x, y, z, w;

	// constructions
	Vector4D() { x = 0;	y = 0;	z = 0;	w = 1.0; }
	Vector4D(double xx, double yy, double zz, double ww = 1.0) { x = xx;	y = yy;	z = zz;	w = ww; }
	Vector4D(const Vector4D& v) { x = v.x;	y = v.y;	z = v.z; w = v.w; }
	Vector4D(const Vector3D& v) { x = v.x;	y = v.y;	z = v.z; w = 1.0; }

	double operator* (const Vector4D& right)
	{
		return x * right.x + y * right.y + z * right.z + w * right.w;
	}

	double& operator[](int idx) {
		assert(idx < 4);
		switch (idx)
		{
		case 0: return x;
		case 1: return y;
		case 2: return z;
		case 3: return w;
		}
	}

	bool operator == (const Vector4D& right) const
	{
		return x == right.x && y == right.y && z == right.z && w == right.w;
	}

	Vector4D& operator=(const Vector3D& v)
	{
		x = v.x; y = v.y; z = v.z; w = 1.0;
		return (*this);
	}
};

////twice the area of a face
double Area2(Vector2D &a, Vector2D &b, Vector2D &c);

double SpcDivision(const double &a, const double &b);

Vector3D Rotate(Vector3D point, double angle, Vector3D rAxis);

bool toLeft(const Vector2D &p, const Vector2D &p0, const Vector2D &p1);

#endif