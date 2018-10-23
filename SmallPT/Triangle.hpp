#ifndef _TRIANGLE_HPP
#define _TRIANGLE_HPP

#include "Shape.hpp"
#include <opencv2/core.hpp>

class Triangle : public Shape 
{
	public:

		Vec p1, p2, p3, color1, color2, color3, normal;

		Triangle(Vec _p1, Vec _p2, Vec _p3, Vec emi, Vec col, ReflectType mat) : p1(_p1), p2(_p2), p3(_p3), Shape(emi, col, mat)
		{
			Vec p1p2 = (_p1 - _p2);
			Vec p1p3 = (_p1 - _p3);
			normal = (p1p2 % p1p3).norm();
			color1 = color2 = color3 = color;
		}

		void SetColor(Vec uv1, Vec uv2, Vec uv3, cv::Mat tex)
		{
			int u = (uv1.x * (float)tex.cols);
			int v = (uv1.y * (float)tex.rows);
			printf("%f %f %f\n", tex.at<cv::Vec3b>(u,v) [0], tex.at<cv::Vec3b>(u,v) [1], tex.at<cv::Vec3b>(u,v) [2]);
			color1 = Vec(tex.at<cv::Vec3b>(u,v) [0], tex.at<cv::Vec3b>(u,v) [1], tex.at<cv::Vec3b>(u,v) [2]);
			printf("%f %f %f\n", color1.x, color1.y, color1.z);
			color1 = Vec(color1.x / 255.0, color1.y / 255.0, color1.z / 255.0);
			printf("%f %f %f\n", color1.x, color1.y, color1.z);
			
			u = (uv2.x * (float)tex.cols);
			v = (uv2.y * (float)tex.rows);
			color2 = Vec(tex.at<cv::Vec3b>(u,v) [0], tex.at<cv::Vec3b>(u,v) [1], tex.at<cv::Vec3b>(u,v) [2]);
			color2 = Vec(color2.x / 255.0, color2.y / 255.0, color2.z / 255.0);

			u = (uv3.x * (float)tex.cols);
			v = (uv3.y * (float)tex.rows);
			color3 = Vec(tex.at<cv::Vec3b>(u,v) [0], tex.at<cv::Vec3b>(u,v) [1], tex.at<cv::Vec3b>(u,v) [2]);
			color3 = Vec(color3.x / 255.0, color3.y / 255.0, color3.z / 255.0);

			//color = (color1 + color2 + color3);
			// 765 for 3.0 (average term) * 255.0 (color maping [0, 255] to [0, 1] term)
			//color = Vec(color.x / 3.0, color.y / 3.0, color.z / 3.0);
			//color = Vec(0.275, 0.528, 0.615);
			color1 = uv1;
			color2 = uv2;
			color3 = uv3;
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

#endif // !TRIANGLE_HPP
