#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008 
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt 
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2 

#include "Ray.hpp"
#include "Vec.hpp"
#include "Sphere.hpp"

const double pi = 3.1415926535;
const double reciprocalPi = 1 / pi;

// Scene Description
Sphere spheres[] = {
	// Radius		Origin							Emission				Color						Material
	Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left
  	Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght
  	Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),//Back
  	Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(),           DIFF),//Frnt
  	Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),//Botm
  	Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top
  	Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, SPEC),//Mirr
  	Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR),//Glas
  	Sphere(1.5, Vec(50,81.6-16.5,81.6),Vec(4,4,4)*100,  Vec(), DIFF),//Lite
};
int numSpheres = 9;

// clamp between 0 and 1
inline double clamp(double x) { return x > 1 ? 1 : x < 0 ? 0 : x; }

// Convert float to int for ppm format, Note that we also apply gamma correction.
//inline int toInt(double x) { return  int(pow(clamp(x), 1.0f / 2.2f) * 255 + 0.5); }
inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }

// get min intersect distance & SphereId (the intersected sphere that is not being occlude from other sphere)
inline bool intersect(const Ray &ray, double &t, int &id)
{
	double infinite = t = 1e20;
	double dis;

	for (int i = 0; i < numSpheres; ++i) {
		dis = spheres[i].intersect(ray);
		if (dis > 0 && dis < t) {
			t = dis;
			id = i;
		}
	}

	return t < infinite;
}

void savePPM(const char* outputFileName, int width, int height, Vec *outputData)
{
	FILE *f = fopen(outputFileName, "w");
	fprintf(f, "P3\n%d %d\n%d\n", width, height, 255); // store value in [0, 255]
	for (int i = 0; i < width * height; ++i) {
		fprintf(f, "%d %d %d ", toInt(outputData[i].x), toInt(outputData[i].y), toInt(outputData[i].z));
	}
	// fclose(f);
	// From smallpt:
	// Instead of fclose()ing the file, I exploit the C++ standard which calls return(0) implicitly, which in turn calls exit(0), which flushes and closes open files. 
}

// Xi: Random number seed
Vec radiance(const Ray &ray, int depth, unsigned short *Xi, int includeEmissive=1) {
	
	double t;
	int id = 0;
	int MAX_DEPTH = 10;
	int numSpheres = 9;

	if (!intersect(ray, t, id))	// miss, return black color
		return Vec();
	
	const Sphere &hitSphere = spheres[id];

	if(depth > MAX_DEPTH)
		return Vec();

	// get surface point
	Vec hitPoint = ray.origin + ray.direction * t;

	// surface normal
	Vec hitPointNormal = (hitPoint - hitSphere.origin).norm();

	// properly oriented surface normal
	// when hit glass surface, ray tracer should determined entering or exiting the surface
	Vec orientedHitPointNormal = hitPointNormal.dot(ray.direction) < 0 ? hitPointNormal : hitPointNormal * -1;

	// object color (BRDF modulator)
	Vec color = hitSphere.color;

	// Russian Roulette, stop the recursion randomly based on surface reflectivity
	// p stands for probability, choose maximum component (r,g,b) of the surface color
	double p = color.x > color.y && color.x > color.z ? color.x : color.y > color.z ? color.y : color.z;
	
	// don't do roulette unless depth > 5
	if (++depth > 5 || !p) { // avoid p = 0
		if (erand48(Xi) < p)
			color = color * (1/p);
		else
			return hitSphere.emission * includeEmissive;
	}

	if (hitSphere.reflectType == DIFF)
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

		for (int i = 0; i < numSpheres; ++i) {
			const Sphere &s = spheres[i];

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

		return hitSphere.emission * includeEmissive + e + color.mult(radiance(Ray(hitPoint, d), depth, Xi, 0));
			
	}
	
	else if (hitSphere.reflectType == SPEC)
	{
		return hitSphere.emission + color.mult(radiance(Ray(hitPoint, ray.direction - hitPointNormal * 2 * hitPointNormal.dot(ray.direction)), depth, Xi));	
	}
	
	// REFR
	else{
		
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
		if(cos2t < 0)
		{
			return hitSphere.emission + color.mult(radiance(reflectionRay, depth, Xi));
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
		if(depth > 2)
		{
			if(erand48(Xi) < p)
			{
				return hitSphere.emission + radiance(reflectionRay, depth, Xi) * Rp;
			}
			else 
			{
				return hitSphere.emission + radiance(Ray(hitPoint, refractionDir), depth, Xi) * Tp;
			}
		}
		else{
				return hitSphere.emission + radiance(reflectionRay, depth, Xi) * Re +  radiance(Ray(hitPoint, refractionDir), depth, Xi) * Tr;
		}
	}

}

int main(int argc, char **argv)
{
	int width = 1024, height = 768; // an image with resolution: 1024 * 768
	int spp = argc == 2 ? atoi(argv[1]) / 4 : 1; // Sample Per Pixel

	Ray Camera(Vec(50.0, 52.0, 295.6), Vec(0.0, -0.042612, -1.0));

	double fov = 0.5135;

	Vec cx = Vec(width * fov / height); // horizontal direction increment, 0.5135f dr
	Vec cy = Vec(cx % Camera.direction).norm() * fov; // vertical direction increment, scalar not works for double

	Vec r; // used for colors of samples
	Vec *output = new Vec[width * height];

	int samps = spp, w = width, h = height;
	
	// start ray tracing
	for (int y = 0; y < height; ++y) {
		fprintf(stderr, "\r Rendering (%d spp) %5.2f%%", spp * 4, 100.0 * y / (height - 1));

		unsigned short Xi[3] = { 0, 0, y * y * y };

		for (unsigned short x = 0; x < width; ++x) {

			// for each pixel do 2x2 subsamples

			// h/0    x      w
			//   -------------	|
			//   -------------	| h - y - 1
			// y .....O-------  |				store as [- - - - - - - O . . . . . . . .]
			//   .............	
			//   .............	
			// 0 .............	

			// . for stored pixel, O for current target pixel, - for not process pixel
			for (int sy = 0, i = (height - y - 1) * width + x; sy < 2; ++sy) {
				for (int sx = 0; sx < 2; ++sx, r = Vec()) { // clear r				
					for (int s = 0; s < spp; ++s) {

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

						Vec dir = cx * ((((sx * 0.5 + dx) / 2 + x) / width - 0.5)) +
							cy * ((((sy * 0.5 + dy) / 2 + y) / height - 0.5)) + Camera.direction;

						// dir may need to normalize, however we force normalize when construct Ray
						r = r + radiance(Ray(Camera.origin + dir * 140, dir.norm()), 0, Xi) * (1.0 / spp);

						// Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) + cy * (((sy + .5 + dy) / 2 + y) / h - .5) + Camera.d;
						// r = r + radiance(Ray(Camera.o + d * 140, d.norm()), 0, Xi) * (1. / samps);
					
					}

					output[i] = output[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z))* .25; // average radiance, 0.25 for 2 spp
				}
			}

		}
	}

	fprintf(stderr, "\n");

	savePPM("output.ppm", width, height, output);

	return 0;
}
