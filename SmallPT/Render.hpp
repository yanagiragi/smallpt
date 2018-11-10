#ifndef _RENDER_HPP
#define _RENDER_HPP

#include "Ray.hpp"
#include "Utils.hpp"
#include "Config.hpp"

#include <omp.h>

namespace smallPT
{
	Vec SampleLights(int originalHitId, bool originalisTriangleHit, Vec camRay, Vec hitPoint, Vec orientedHitPointNormal, Vec color, unsigned short *Xi)
	{
		//Sampling Lights, loop over any lights
		Vec e;

		Scene Scene = MainScene;

		double t;
		int id;

		// Sample Spheres
		for (int i = 0; i < Scene.Spheres.size(); ++i) {
			const Sphere &s = Scene.Spheres[i];

			if (s.material.emission.x <= 0 && s.material.emission.y <= 0 && s.material.emission.z <= 0)
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
			if (Scene.intersect(Ray(hitPoint, l), t, id) && id == i) {
				double omega = 2 * pi * (1 - cos_a_max); // 1/probability with respect to solid angle

														 // calculating lighting and add to current value
				e = e + color.mult(s.material.emission * l.dot(orientedHitPointNormal) * omega) * reciprocalPi; // 1/pi for BRDF Energy equals 1
			}
		}

		//printf("%d (%f %f %f)\n", hitId, hitPoint.x, hitPoint.y, hitPoint.z);

		// Sample Triangles	
		for (int i = 0; i < Scene.Triangles.size(); ++i) {
			const Triangle &tri = Scene.Triangles[i];

			if (tri.material.emission.x <= 0 && tri.material.emission.y <= 0 && tri.material.emission.z <= 0)
				continue; // skip non-lights

			double u1 = erand48(Xi);
			double u2 = erand48(Xi);

			// From Pbrt v3 & Minilight
			double su0 = sqrt(u1);
			double esp1 = 1 - su0;
			double esp2 = (1 - u2) * su0;

			Vec hitPointTri = (tri.p3 - tri.p1) * esp1 + (tri.p2 - tri.p1) * esp2 + tri.p1;
			//Vec hitPointTri = tri.p1 * esp1 + tri.p2 * esp2 + tri.p3 * (1 - esp1 - esp2);
			Vec l = hitPointTri - hitPoint;

			// shoot shadow ray
			if (Scene.intersect(Ray(hitPoint, l), t, id) && id == i) {

				/*bool isTriangleHit = (id >= Scene.Spheres.size());
				int hitId = isTriangleHit ? id - Scene.Spheres.size() : id;

				if(originalisTriangleHit && originalHitId == i){
				printf("Occulded ");
				continue;
				}*/

				double area = tri.Area();
				Vec outDirection = hitPoint - hitPointTri;//l.mult(Vec(-1, -1, -1));
														  //Vec outDirection = l.mult(Vec(-1, -1, -1));
				double distance2 = l.dot(l);

				Vec interpolatedNormal = (tri.n3 - tri.n1) * esp1 + (tri.n2 - tri.n1) * esp2 + tri.n1;

				//double cosArea = area * outDirection.dot(interpolatedNormal);
				double cosArea = area * outDirection.dot(tri.normal);

				double omega = cosArea / (distance2 >= 1e-6 ? distance2 : 1e-6); // 1 / probability

				// 看起來應該要 * pi 把 reciprocalPi cancel 掉
				omega *= pi;

				if (cosArea < 0)
					continue;

				// ?: 
				// l is parrellel to orirnetedHitPoointNormal
				// since its should dot with ray from camera
				//printf("%f\n", l.dot(camRay) * omega);
				//printf("%f\n", omega);

				// calculating lighting and add to current value

				// 原始 Sphere 的版本
				//e = e + color.mult(tri.material.emission * l.dot(camRay) * omega) * reciprocalPi; // 1/pi for BRDF Energy equals 1

				// 發現 l.dot(orientedHitPointNormal) = 0 之後嘗試修改的版本
				//e = e + color.mult(tri.material.emission * omega) * reciprocalPi; // 1/pi for BRDF Energy equals 1

				// 顏色 * Emission
				//e = e + color.mult(tri.material.emission) * 100;

				// 感覺 * 100 亮度沒有太大的變化？
				//e = e + color * 100;

				// 牆壁有紅藍顏色了
				e = e + color;

				//e = e + color * omega;
			}
		}

		return e;
	}


	Vec radiance2(const Ray &r, int depth, unsigned short *Xi) {

		Scene Scene = MainScene;

		double t;                               // distance to intersection
		int id = 0;                               // id of intersected object
		if (!Scene.intersect(r, t, id))
			return Vec(); // if miss, return black
		const Sphere &obj = Scene.Spheres[id];        // the hit object

		Vec x = r.origin + r.direction*t;
		Vec n = (x - obj.origin).norm();
		Vec nl = n.dot(r.direction) < 0 ? n : n * -1;
		Vec f = obj.material.color;

		double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl

		if (++depth>5)
			if (erand48(Xi)<p) f = f * (1 / p);
			else
				return obj.material.emission; //R.R.

		if (depth > 10)
			return obj.material.emission;

		if (obj.material.reflectType == DIFF) {                  // Ideal DIFFUSE reflection
			double r1 = 2 * pi *erand48(Xi);
			double r2 = erand48(Xi);
			double r2s = sqrt(r2);

			Vec w = nl, u = ((fabs(w.x)>.1 ? Vec(0, 1) : Vec(1)) % w).norm();
			Vec v = w % u;

			// next camera ray's direction
			Vec d = (u*cos(r1)*r2s + v * sin(r1)*r2s + w * sqrt(1 - r2)).norm();

			//Vec e = SampleLights(id, false, r.direction, x, nl, obj.material.color, Xi);

			// f = obj.c;
			return obj.material.emission + f.mult(radiance2(Ray(x, d), depth, Xi));
		}

		else if (obj.material.reflectType == SPEC) {
			// Ideal SPECULAR reflection
			return obj.material.emission + f.mult(radiance2(Ray(x, r.direction - n * 2 * n.dot(r.direction)), depth, Xi));
		}

		// REFL
		Ray reflRay(x, r.direction - n * 2 * n.dot(r.direction));     // Ideal dielectric REFRACTION
		bool into = n.dot(nl)>0;                // Ray from outside going in?

		double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.direction.dot(nl), cos2t;
		if ((cos2t = 1 - nnt * nnt*(1 - ddn * ddn))<0)    // Total internal reflection
			return obj.material.emission + f.mult(radiance2(reflRay, depth, Xi));

		Vec tdir = (r.direction *nnt - n * ((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm();

		double a = nt - nc, b = nt + nc, R0 = a * a / (b*b), c = 1 - (into ? -ddn : tdir.dot(n));

		double Re = R0 + (1 - R0)*c*c*c*c*c;
		double Tr = 1 - Re;
		double P = .25 + .5*Re;
		double RP = Re / P;
		double TP = Tr / (1 - P);

		return obj.material.emission + f.mult(depth>2 ? (erand48(Xi)<P ?   // Russian roulette
			radiance2(reflRay, depth, Xi)*RP : radiance2(Ray(x, tdir), depth, Xi)*TP) :
			radiance2(reflRay, depth, Xi)*Re + radiance2(Ray(x, tdir), depth, Xi)*Tr);
	}

	Vec radiance(const Ray &ray, int depth, unsigned short *Xi, int includeEmissive = 1) {
		
		// Xi: Random number seed
		Scene Scene = MainScene;

		double t;
		int id = 0;
		int MAX_DEPTH = 10;

		if (!Scene.intersect(ray, t, id))	// miss, return black color
			return Vec();

		Vec origin;

		bool isTriangleHit = (id >= Scene.Spheres.size());

		int hitId = isTriangleHit ? id - Scene.Spheres.size() : id;
		Shape &hitShape = isTriangleHit ? (Shape&)(Scene.Triangles[hitId]) : (Shape&)(Scene.Spheres[hitId]);

		if (isTriangleHit) {
			origin = Scene.Triangles[hitId].p1 + Scene.Triangles[hitId].p2 + Scene.Triangles[hitId].p3;
			origin.x = origin.x / 3.0;
			origin.y = origin.y / 3.0;
			origin.z = origin.z / 3.0;
		}
		else {
			origin = ((Sphere&)hitShape).origin;
		}

		// get surface point
		Vec hitPoint = ray.origin + ray.direction * t;

		// surface normal
		Vec hitPointNormal;

		if (isTriangleHit) {
			hitPointNormal = Scene.Triangles[hitId].normal;
		}
		else {
			hitPointNormal = (hitPoint - origin).norm();
		}

		// properly oriented surface normal
		// when hit glass surface, ray tracer should determined entering or exiting the surface
		Vec orientedHitPointNormal = hitPointNormal.dot(ray.direction) < 0 ? hitPointNormal : hitPointNormal * -1;

		// object color (BRDF modulator)
		Vec color;
		if (isTriangleHit) {
			Triangle tri = (Triangle&)(hitShape);
			float distance1 = (hitPoint - tri.p1).mag();
			float distance2 = (hitPoint - tri.p2).mag();
			float distance3 = (hitPoint - tri.p3).mag();
			float totalDistance = distance1 + distance2 + distance3;

			Vec uv;
			uv.x = (tri.uv1.x * distance1 + tri.uv2.x * distance2 + tri.uv3.x * distance3) / totalDistance;
			uv.y = (tri.uv1.y * distance1 + tri.uv2.y * distance2 + tri.uv3.y * distance3) / totalDistance;

			if (tri.hasUV) {
				int texId = 0;
				/*int w = Scene.Textures[texId].cols;
				int h = Scene.Textures[texId].rows;
				cv::Mat tex = Scene.Textures[texId];
				color = Vec(
				tex.at<cv::Vec3b>(uv.x * w, uv.y * h)[0] / 255.0,
				tex.at<cv::Vec3b>(uv.x * w, uv.y * h)[1] / 255.0,
				tex.at<cv::Vec3b>(uv.x * w, uv.y * h)[2] / 255.0
				);*/
			}
			else {
				color = hitShape.material.color;
			}

			// TODO: wierd uv interpolation?
			// use hitShape.color instead for now.
			color = hitShape.material.color;
		}
		else {
			// no uv support for sphere for now
			color = hitShape.material.color;
		}

		// Russian Roulette, stop the recursion randomly based on surface reflectivity
		// p stands for probability, choose maximum component (r,g,b) of the surface color
		double p = color.x > color.y && color.x > color.z ? color.x : color.y > color.z ? color.y : color.z;

		// don't do roulette unless depth > 5
		if (++depth > 5 || !p) { // avoid p = 0
			if (erand48(Xi) < p)
				color = color * (1 / p);
			else
				return hitShape.material.emission * includeEmissive;
		}

		if (depth > MAX_DEPTH)
			return hitShape.material.emission;

		if (hitShape.material.reflectType == DIFF)
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
			// Sampling Unit Disk . Sampling Unit Hemisphere

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

			Vec e = SampleLights(hitId, isTriangleHit, ray.direction, hitPoint, orientedHitPointNormal, color, Xi);
			//Vec e = Vec(0,0,0);

			//Vec mask = Vec(1, 1, 1);
			/*mask = mask * color;
			mask = mask * d.dot(orientedHitPointNormal);
			mask = mask * 2;*/

			return hitShape.material.emission * includeEmissive + e + color.mult(radiance(Ray(hitPoint, d), depth, Xi, 1));
			//return hitShape.material.emission * includeEmissive + e + color.mult(radiance(Ray(hitPoint, d), depth, Xi, 0));
			//return hitShape.material.emission * includeEmissive + color.mult(radiance(Ray(hitPoint, d), depth, Xi, 0));

		}
		else if (hitShape.material.reflectType == SPEC)
		{
			return hitShape.material.emission + color.mult(radiance(Ray(hitPoint, ray.direction - hitPointNormal * 2 * hitPointNormal.dot(ray.direction)), depth, Xi));
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
				return hitShape.material.emission + color.mult(radiance(reflectionRay, depth, Xi));
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
					return hitShape.material.emission + radiance(reflectionRay, depth, Xi) * Rp;
				}
				else
				{
					return hitShape.material.emission + radiance(Ray(hitPoint, refractionDir), depth, Xi) * Tp;
				}
			}
			else {
				return hitShape.material.emission + radiance(reflectionRay, depth, Xi) * Re + radiance(Ray(hitPoint, refractionDir), depth, Xi) * Tr;
			}
		}
	}

	void UpdateRendering()
	{
		if (currentSpp < limitSpp) {
			if (currentSpp != 0) {
				std::ostringstream outputFilenameStringStream;
				outputFilenameStringStream << SaveImageNamePrefix << currentSpp << "spp";

				std::string outputPpmFilename = outputFilenameStringStream.str() + ".ppm";
				std::string outputPngFilename = outputFilenameStringStream.str() + ".png";

				// Save PPM
				SavePPM(outputPpmFilename, width, height, output);
			}

			++currentSpp;
		}
		else {
			return;
		}

		//const int spp = currentSpp;
		const int spp = limitSpp;
		/*const int width = width;
		const int height = height;
		const int channel = channel;
		Vec* output = output;
		float* pixels = pixels;*/

		Ray Camera(Vec(50.0, 52.0, 295.6), Vec(0.0, -0.042612, -1.0));

		double fov = 0.5135;

		Vec cx = Vec(width * fov / height); // horizontal direction increment, 0.5135f dr
		Vec cy = Vec(cx % Camera.direction).norm() * fov; // vertical direction increment, scalar not works for double

		Vec r; // used for colors of samples

		// set 6 threads
		omp_set_num_threads(6);

		// start ray tracing
		#pragma omp parallel for schedule(dynamic, 1) private(r)		
		for (int y = 0; y < height; ++y) {
			fprintf(stderr, "\r Rendering (%d spp) %5.2f%%", spp, 100.0 * y / (height - 1));

			//unsigned short Xi[3] = { 0, 0, y * y * y};
			unsigned short Xi[3] = { rand(), rand(), rand() };

			for (unsigned short x = 0; x < width; ++x) {
				const int i = (height - y - 1) * width + x;
				const int i2 = i * 2;

				// Store radiance increased after 2x2 subsamples
				Vec increment;

				// for each pixel do 2x2 subsamples
				for (int sy = 0; sy < 2; ++sy) {
					for (int sx = 0; sx < 2; ++sx, r = Vec()) { // clear r				
						for (int s = 0; s < spp; ++s) {
							/*const float r1 = GetRandom(&seeds[i2], &seeds[i2 + 1]) - .5f;
							const float r2 = GetRandom(&seeds[i2], &seeds[i2 + 1]) - .5f;
							const float kcx = (x + r1) * 1.0 / width - .5f;
							const float kcy = (y + r2) * 1.0 / height - .5f;

							Vec dir = cx * kcx + cy * kcy + Camera.direction;*/


							double r1 = 2 * erand48(Xi); // r1 in [0, 2)
							double dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);

							// dx in [-1, 1)
							// Apply tent filter to transform canonical random samples to nonuniform samples (for anti-aliasing)
							// tent filter: https://computergraphics.stackexchange.com/questions/3868/why-use-a-tent-filter-in-path-tracing
							// or : Realistic Ray Tracing, Second Edition: Peter Shirley, R. Keith Morley 						

							double r2 = 2 * erand48(Xi);
							double dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);

							// Here dx and dy determins location of sample within pixel
							// dx and dy are now in [-1,1)
							// (sx * 0.5 + dx) / 2 in [ (-1-sx)/2, (1+sx)/2 )
							// 這邊的意思是 (1+sx)/2 以當前pixel為中心，要sample的距離 = 0.5 + (spp - 1) / 2, 0.5 為 中心pixel的寬度=1 / 2
							// 例如 2Spp, 就是 sample 當前 pixel 中心向左 1.5 向右1.5 的距離 (剛剛好就是自己pixel範圍 + 左右各1/2 pixel寬度的範圍)
							// ((sx * 0.5 + dx)/2 + x) / width 會把這樣的範圍mapping 到 [0,1) 
							// 再扣掉 0.5 mapping 到 [-0.5, 0.5) 做 normalize
							// 最後就可以產生以 (x,y) 這個 pixel 為中心 撒出 落在 x 方向 和 y 方向 增加 [-0.5,0.5) 區間值的 方向

							Vec dir = cx * ((((sx * 0.5 + dx) / 2 + x) / width - 0.5)) + cy * ((((sy * 0.5 + dy) / 2 + y) / height - 0.5)) + Camera.direction;

							// dir may need to normalize, however we force normalize when construct Ray
							// r = r + radiance(Ray(Camera.origin + dir * 140, dir.norm()), 0, Xi) * (1.0 / spp);					
							r = r + radiance(Ray(Camera.origin + dir * 140, dir.norm()), 0, Xi) * (1.0 / spp);
						}
						increment = increment + Vec(clamp(r.x), clamp(r.y), clamp(r.z))* .25; // average radiance, 0.25 for 2 sub spp
					}
				}

				// normalize color if progressive
				/*const float k1 = spp;
				const float k2 = 1.f / (k1 + 1.f);
				if(spp == 1){
				output[i] = increment;
				}
				else{
				output[i] = (output[i] * k1 + increment) * k2;
				} */

				output[i] = increment;

				// flip Y for OpenGL Coordinate
				// Now no flip ?
				const int j = y * width + x;
				pixels[j * channel + 0] = gammaCorrection(output[i].x);
				pixels[j * channel + 1] = gammaCorrection(output[i].y);
				pixels[j * channel + 2] = gammaCorrection(output[i].z);
			}
		}

		fprintf(stderr, "\n");
	}


}

#endif // !_RENDER_HPP