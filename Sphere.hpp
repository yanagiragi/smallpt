#ifndef SPHERE_HPP
#define SPHERE_HPP

#include "Ray.hpp"
#include "Vec.hpp"

const double epsilon = 1e-4;

enum ReflectType
{
	DIFF, SPEC, REFR
};

class Sphere
{
	public:

		double radius;

		Vec origin;

		Vec emission;

		Vec color;

		ReflectType reflectType;

		Sphere(double rad, Vec ori, Vec emi, Vec col, ReflectType mat) : radius(rad), origin(ori), emission(emi), color(col), reflectType(mat){}

		double intersect_old(const Ray &ray) const
		{
			// Sphere Intersection Equation: (t stands for ray intersect distance)
			// (ray.origin + t * ray.direction - origin) * (ray.origin + t * ray.direction - origin) - radius^2 = 0
			// rewrite to form: at^2 + b^t + c = 0
			double a = ray.direction.dot(ray.direction); // actually is should be 1, since direction of ray is normalized
			double b = 2.0 * ray.direction.dot(ray.origin - origin);
			double c = (ray.origin - origin).dot(ray.origin - origin) - radius * radius;

			if ((b * b - 4.0 * a * c) < 0) // not intersected
				return 0; 

			// solve sphere intersection equation, ans = (-b +- sqrt(b^2 - 4ac)) / 2a
			double t1 = (-1.0 * b + sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
			double t2 = (-1.0 * b - sqrt(b * b - 4.0 * a * c)) / (2.0 * a);

			return t1 > epsilon ? t1 
								: (t2 > epsilon) ? t2 : 0;
		}

		double intersect(const Ray &r) const { // returns distance, 0 if nohit 
			Vec op = origin - r.origin; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0 
			double t, eps = 1e-4, b = op.dot(r.direction), det = b*b - op.dot(op) + radius*radius;
			if (det<0) return 0; else det = sqrt(det);
			return (t = b - det)>eps ? t : ((t = b + det)>eps ? t : 0);
		}
};

#endif // !SPHERE_HPP