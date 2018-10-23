#ifndef _RENDER_HPP
#define _RENDER_HPP

#include "Ray.hpp"
#include "Utils.hpp"

// Xi: Random number seed
Vec radiance(const Ray &ray, int depth, unsigned short *Xi, int includeEmissive = 1) {

	double t;
	int id = 0;
	int MAX_DEPTH = 10;

	if (!intersect(ray, t, id))	// miss, return black color
		return Vec();

	Vec origin;

	bool isTriangleHit = (id >= Spheres.size());

	int hitId = isTriangleHit ? id - Spheres.size() : id;
	Shape &hitShape = isTriangleHit ? (Shape&)(Triangles[hitId]) : (Shape&)(Spheres[hitId]);

	if (isTriangleHit) {
		/*origin = Triangles[hitId].p1 + Triangles[hitId].p2 + Triangles[hitId].p3;
		origin.x = origin.x / 3.0;
		origin.y = origin.y / 3.0;
		origin.z = origin.z / 3.0;*/
	}
	else {
		origin = ((Sphere&)hitShape).origin;
	}

	if (depth > MAX_DEPTH)
		return Vec();

	// get surface point
	Vec hitPoint = ray.origin + ray.direction * t;

	// surface normal
	Vec hitPointNormal;
	
	if(isTriangleHit){
		hitPointNormal = Triangles[hitId].normal;
	}
	else{
		hitPointNormal = (hitPoint - origin).norm();
	}

	// properly oriented surface normal
	// when hit glass surface, ray tracer should determined entering or exiting the surface
	Vec orientedHitPointNormal = hitPointNormal.dot(ray.direction) < 0 ? hitPointNormal : hitPointNormal * -1;

	// object color (BRDF modulator)
	Vec color = hitShape.color;

	// Russian Roulette, stop the recursion randomly based on surface reflectivity
	// p stands for probability, choose maximum component (r,g,b) of the surface color
	double p = color.x > color.y && color.x > color.z ? color.x : color.y > color.z ? color.y : color.z;

	// don't do roulette unless depth > 5
	if (++depth > 5 || !p) { // avoid p = 0
		if (erand48(Xi) < p)
			color = color * (1 / p);
		else
			return hitShape.emission * includeEmissive;
	}

	if (hitShape.reflectType == DIFF)
	{
		// random angle
		double r1 = 2 * pi * erand48(Xi);
		// random distance from hitSphere center
		double r2 = erand48(Xi);
		double randomDistanceFromCenter = sqrt(r2);

		// From Realistic Ray Tracing: 
		// We shouldn't take horizontal and vertical dimensions uniformly to (r, phi)
		// because it does not preserve relative area
		// 簡單來說，射到surface之後，這個射線再度射到Image Plane 的機率分布並不是在水平跟垂直均勻分布的
		// 因此，如果我們分開 randomly sample, 會不符合實際上的機率分布
		// 所以，我們是對 Unit Disk 做 Sample，再用他去反推回去到半球的solid angle 求出新的方向
		// Sampling Unit Disk -> Sampling Unit Hemisphere

		// create orthonormal coordinate frame (w,u,v)
		Vec w = orientedHitPointNormal;
		/*glm::vec3 u = glm::normalize(glm::cross(fabs(w.x) > 0.1 ? glm::vec3(0, 1, 0) : glm::vec3(1, 0, 0), w));
		glm::vec3 v = glm::normalize(glm::cross(u, w));*/
		Vec u = ((fabs(w.x) > 0.1 ? Vec(0, 1) : Vec(1)) % w).norm();
		Vec v = w % u;

		// new random generated ray direction, this ray's origin = hitPoint
		// (cos term, sin term, last term), last term = sqrt(1 - (cos term)^2 + (sin term)^ 2)
		// since cos^2 + sin^2 = 1, last term = sqrt ( 1 - randomDistanceFromCenter ^ 2) = sqrt(1 - r2)
		Vec d = (u * cos(r1) * randomDistanceFromCenter + v * sin(r1) * randomDistanceFromCenter + w * sqrt(1 - r2)).norm();

		//Sampling Lights, loop over any lights
		Vec e;

		for (int i = 0; i < Spheres.size(); ++i) {
			const Sphere &s = Spheres[i];

			if (s.emission.x <= 0 && s.emission.y <= 0 && s.emission.z <= 0)
				continue; // skip non-lights

						  // create random direction towards sphere
			Vec sw = s.origin - hitPoint;
			Vec su = ((fabs(sw.x) > 0.1 ? Vec(0, 1) : Vec(1)) % sw).norm();
			Vec sv = sw % su;

			// Sampling Sphere by Solid Angle
			double cos_a_max = sqrt(1 - s.radius * s.radius / (hitPoint - s.origin).dot(hitPoint - s.origin));
			double esp1 = erand48(Xi), esp2 = erand48(Xi);
			double cos_a = 1 - esp1 + esp1 * cos_a_max;
			double sin_a = sqrt(1 - cos_a * cos_a);
			double phi = 2 * pi * esp2;

			// From Realistic Ray Tracing:
			// [ cos_a ] = [ 1 + esp1(cos_a_max - 1) ]
			// [  phi  ]   [      2 * pi * esp2      ]
			Vec l = su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a;
			l.norm();

			// shoot shadow ray
			if (intersect(Ray(hitPoint, l), t, id) && id == i) {
				double omega = 2 * pi * (1 - cos_a_max); // 1/probability with respect to solid angle

														 // calculating lighting and add to current value
				e = e + color.mult(s.emission * l.dot(orientedHitPointNormal) * omega) * reciprocalPi; // 1/pi for BRDF Energy equals 1
			}
		}

		return hitShape.emission * includeEmissive + e + color.mult(radiance(Ray(hitPoint, d), depth, Xi, 0));

	}
	else if (hitShape.reflectType == SPEC)
	{
		return hitShape.emission + color.mult(radiance(Ray(hitPoint, ray.direction - hitPointNormal * 2 * hitPointNormal.dot(ray.direction)), depth, Xi));
	}
	else { // REFR

		// Glass: Reflect & Refract

		Ray reflectionRay(hitPoint, ray.direction - hitPointNormal * 2 * hitPointNormal.dot(ray.direction));

		bool into = hitPointNormal.dot(orientedHitPointNormal) > 0; // ray should going in?

		double nc = 1; // c is the speed of light in vacuum
		double nt = 1.5; // t is the phase velocity of light in the medium

		// refractiveIndex, glass IOR ~= 1.5
		// Check: https://en.wikipedia.org/wiki/Crown_glass_(optics)
		// Not account for dispersion (to account these, vary index by wavelength)
		double nnt = into ? nc / nt : nt / nc;

		double ddn = ray.direction.dot(orientedHitPointNormal);
		double cos2t; // cos(sida_b)

		// if total internel reflection, reflect
		// happens when angle is too shallow
		// cos2t = cos(sida_b) ^ 2
		cos2t = 1 - nnt * nnt * (1 - ddn * ddn);
		if (cos2t < 0)
		{
			return hitShape.emission + color.mult(radiance(reflectionRay, depth, Xi));
		}

		// other wise, choose refraction or reflection using fresnel
		// Check Slide: page 72
		Vec refractionDir = (ray.direction * nnt - hitPointNormal * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();

		double a = nt - nc;
		double b = nt + nc;
		double R0 = a * a / (b * b); // Reflectance at normal incidence based on IOR
		double c = 1 - (into ? -ddn : refractionDir.dot(hitPointNormal)); // 1 - cos(theta)

		double Re = R0 + (1 - R0) * c * c * c * c * c; // Fresnel Reflectance
		double Tr = 1 - Re;

		double p = 0.25 + 0.5 * Re; // probably of reflecting

		double Rp = Re / p;
		double Tp = Tr / (1 - p);

		// Russuan Roulette
		if (depth > 2)
		{
			if (erand48(Xi) < p)
			{
				return hitShape.emission + radiance(reflectionRay, depth, Xi) * Rp;
			}
			else
			{
				return hitShape.emission + radiance(Ray(hitPoint, refractionDir), depth, Xi) * Tp;
			}
		}
		else {
			return hitShape.emission + radiance(reflectionRay, depth, Xi) * Re + radiance(Ray(hitPoint, refractionDir), depth, Xi) * Tr;
		}
	}
}

#endif // !_RENDER_HPP