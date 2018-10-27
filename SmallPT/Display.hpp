#ifndef _DISPLAY_HPP
#define _DISPLAY_HPP

#include <GL/freeglut.h>

#include <omp.h>

#include "Ray.hpp"
#include "Vec.hpp"
#include "Utils.hpp"

#include "Config.hpp"

static float GetRandom(unsigned int *seed0, unsigned int *seed1) {
	*seed0 = 36969 * ((*seed0) & 65535) + ((*seed0) >> 16);
	*seed1 = 18000 * ((*seed1) & 65535) + ((*seed1) >> 16);

	unsigned int ires = ((*seed0) << 16) + (*seed1);

	/* Convert to float */
	union {
		float f;
		unsigned int ui;
	} res;
	res.ui = (ires & 0x007fffff) | 0x40000000;

	return (res.f - 2.f) / 2.f;
}

void UpdateRendering()
{
	if(globalConfig::currentSpp < globalConfig::limitSpp){
		if(globalConfig::currentSpp != 0){
			std::ostringstream outputFilenameStringStream;
			outputFilenameStringStream << globalConfig::SaveImageNamePrefix << globalConfig::currentSpp << "spp";

			std::string outputPpmFilename = outputFilenameStringStream.str() + ".ppm";
			std::string outputPngFilename = outputFilenameStringStream.str() + ".png";

			// Save PPM
			SavePPM(outputPpmFilename,globalConfig::width, globalConfig::height, globalConfig::output);
		}
		
        ++ globalConfig::currentSpp;
	}
	else{
		return;
	}

	const int spp = globalConfig::currentSpp;
	const int width = globalConfig::width;
	const int height = globalConfig::height;
	const int channel = globalConfig::channel;
	Vec* output = globalConfig::output;
	float* pixels = globalConfig::pixels;

	Ray Camera(Vec(50.0, 52.0, 295.6), Vec(0.0, -0.042612, -1.0));

	double fov = 0.5135;

	Vec cx = Vec(width * fov / height); // horizontal direction increment, 0.5135f dr
	Vec cy = Vec(cx % Camera.direction).norm() * fov; // vertical direction increment, scalar not works for double

	Vec r; // used for colors of samples
	
	// set 6 threads
	omp_set_num_threads(6);

	#pragma omp parallel for schedule(dynamic, 1) private(r)
	// start ray tracing
	for (int y = 0; y < height; ++y) {
		fprintf(stderr, "\r Rendering (%d spp) %5.2f%%", spp, 100.0 * y / (height - 1));

		//unsigned short Xi[3] = { 0, 0, y * y * y};
		unsigned short Xi[3] = { rand(), rand(), rand()};

		for (unsigned short x = 0; x < width; ++x) {
				const int i =  (height - y - 1) * width + x;
				const int i2 = i * 2;

				// Store radiance increased after 2x2 subsamples
				Vec increment;

				// for each pixel do 2x2 subsamples
				for (int sy = 0; sy < 2; ++sy) {
					for (int sx = 0; sx < 2; ++sx, r = Vec()) { // clear r				
						//for (int s = 0; s < spp; ++s) {                                             
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
							r = r + radiance(Ray(Camera.origin + dir * 140, dir.norm()), 0, Xi);// * (1.0 / spp);					
						//}
							increment = increment + Vec(clamp(r.x), clamp(r.y), clamp(r.z))* .25; // average radiance, 0.25 for 2 sub spp
					}
				}

				// normalize color
				const float k1 = spp;
				const float k2 = 1.f / (k1 + 1.f);
				if(spp == 1){
						output[i] = increment;
				}                                
				else{
						output[i] = (output[i] * k1 + increment) * k2;
				}                                        

				// flip Y for OpenGL Coordinate
				// Now no flip ?
				const int j =  y * width + x;
				pixels[ j * channel + 0] = gammaCorrection(output[i].x);
				pixels[ j * channel + 1] = gammaCorrection(output[i].y);
				pixels[ j * channel + 2] = gammaCorrection(output[i].z);
		}
	}

	fprintf(stderr, "\n");
}

void idleFunc(void) {
    UpdateRendering();        
	glutPostRedisplay();
}

void displayFunc(void) {
	glClear(GL_COLOR_BUFFER_BIT);
        glRasterPos2i(-1, -1);
	glDrawPixels(globalConfig::width, globalConfig::height, GL_RGB, GL_FLOAT, globalConfig::pixels);
	glutSwapBuffers();
}

void reshapeFunc(int newWidth, int newHeight) {
	// no reshape!
	/*globalConfig::width = newWidth;
	globalConfig::height = newHeight;
	glViewport(0, 0, globalConfig::width, globalConfig::height);
	glLoadIdentity();	
	glutPostRedisplay();*/
}

void SetupGlutDisplay(int argc, char **argv)
{
        glutInit(&argc, argv);
        glutInitWindowSize(globalConfig::width, globalConfig::height);
        glutInitWindowPosition(0,0);
        glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);

        glutCreateWindow("Small PT");

        glutReshapeFunc(reshapeFunc);
        glutDisplayFunc(displayFunc);
        glutIdleFunc(idleFunc);

        //glutKeyboardFunc(keyFunc);
        //glutSpecialFunc(specialFunc);
}

#endif