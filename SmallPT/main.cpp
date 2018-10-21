 // smallpt, a Path Tracer by Kevin Beason, 2008 
 // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt 
 // Remove "-fopenmp" for g++ version < 4.2 

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <sstream>

#include "Ray.hpp"
#include "Vec.hpp"
#include "Scene.hpp"
#include "Render.hpp"
#include "Utils.hpp"

int main(int argc, char **argv)
{
	// Setup Scenes
	LoadScene(argv);

	int width = 1024, height = 768; // an image with resolution: 1024 * 768
	int spp = atoi(argv[1]); // Sample Per Pixel

	if (spp <= 0)
		spp = 1;

	Ray Camera(Vec(50.0, 52.0, 295.6), Vec(0.0, -0.042612, -1.0));

	double fov = 0.5135;

	Vec cx = Vec(width * fov / height); // horizontal direction increment, 0.5135f dr
	Vec cy = Vec(cx % Camera.direction).norm() * fov; // vertical direction increment, scalar not works for double

	Vec r; // used for colors of samples
	Vec *output = new Vec[width * height];

	int samps = spp, w = width, h = height;
	
	// set 6 threads
	omp_set_num_threads(6);

	#pragma omp parallel for schedule(dynamic, 1) private(r)
	// start ray tracing
	for (int y = 0; y < height; ++y) {
		fprintf(stderr, "\r Rendering (%d spp) %5.2f%%", spp, 100.0 * y / (height - 1));

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

	SaveResult(argv, width, height, spp, output);

	return 0;
}
