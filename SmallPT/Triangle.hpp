#ifndef _TRIANGLE_HPP
#define _TRIANGLE_HPP

#include "Shape.hpp"
// #include <opencv2/core.hpp>

namespace smallPT {
	class Triangle : public Shape
	{
	public:

		Vec p1, p2, p3, uv1, uv2, uv3, normal, n1, n2, n3;

		// tinyobj::material_t objMaterial;

		bool hasUV;

		Triangle(Vec _p1, Vec _p2, Vec _p3, bool _hasUV, Vec emi, Vec col, ReflectType mat) :
			p1(_p1), p2(_p2), p3(_p3), hasUV(_hasUV), Shape(emi, col, mat)
		{
			Vec p1p2 = (_p1 - _p2);
			Vec p1p3 = (_p1 - _p3);
			normal = (p1p2 % p1p3).norm();
		}

		Triangle(Vec _p1, Vec _p2, Vec _p3, Vec _n1, Vec _n2, Vec _n3, bool _hasUV, Material mat) :
			p1(_p1), p2(_p2), p3(_p3), hasUV(_hasUV), Shape(mat), n1(_n1), n2(_n2), n3(_n3)
		{
			Vec p1p2 = (_p1 - _p2);
			Vec p1p3 = (_p1 - _p3);
			Vec nn = (p1p2 % p1p3).norm();
			normal = (_n1 + _n2 + _n3).norm();

			// printf("%f %f %f, %f %f %f\n", nn.x, nn.y, nn.z, normal.x, normal.y, normal.z);
		}

		void SetUV(Vec _uv1, Vec _uv2, Vec _uv3)
		{
			uv1 = _uv1;
			uv2 = _uv2;
			uv3 = _uv3;
		}

		double Area() const {
			// 1/2 * cross().length
			Vec v1 = p2 - p1;
			Vec v2 = p3 - p2;
			return (v1 % v2).mag() * 0.5;
		}

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

}
#endif // !TRIANGLE_HPP
