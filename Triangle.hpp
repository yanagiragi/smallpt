#ifndef _TRIANGLE_HPP
#define _TRIANGLE_HPP

#include "Shape.hpp"

class Triangle : public Shape 
{
	public:

		Vec p1, p2, p3;

		Triangle(Vec _p1, Vec _p2, Vec _p3, Vec emi, Vec col, ReflectType mat) : p1(_p1), p2(_p2), p3(_p3), Shape(emi, col, mat){}

		double intersect(const Ray &r) const { // returns distance, 0 if nohit 
			
			//#ifdef MOLLER_TRUMBORE
				Vec p1p2 = p2 - p1;
				Vec p1p3 = p3 - p1;
				Vec pvec = Vec(r.direction) % p1p3; // cross
				double det = p1p2.dot(pvec);
			//#ifdef CULLING
				// if the determinant is negative the triangle is backfacing
				// if the determinant is close to 0, the ray misses the triangle
				if (det < epsilon) return 0;
			//#else
				// // ray and triangle are parallel if det is close to 0
				// if (fabs(det) < epsilon) return 0;
			//#endif
				double invDet = 1 / det;
			
				Vec tvec = r.origin - p1;
				double u = tvec.dot(pvec) * invDet;
				if (u < 0 || u > 1) return 0;
				
				Vec qvec = tvec % p1p2;
				double v = r.direction.dot(qvec) * invDet;
				if (v < 0 || u + v > 1) return 0;
				
				double t = p1p3.dot(qvec) * invDet;
				
				return t;
		}
};

#endif // !TRIANGLE_HPP
