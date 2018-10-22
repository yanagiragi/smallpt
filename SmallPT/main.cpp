 // smallpt, a Path Tracer by Kevin Beason, 2008 
 // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt 
 // Remove "-fopenmp" for g++ version < 4.2 

#include <stdlib.h>
#include <stdio.h>
#include <sstream>

#include "Scene.hpp"
#include "Render.hpp"
#include "Utils.hpp"
#include "Display.hpp"

int main(int argc, char **argv)
{
	// Setup Scenes
	LoadScene(argv);

	// Sample Per Pixel
	limitSpp = atoi(argv[1]) <= 0 ? 1 : atoi(argv[1]); 

 	output = new Vec[width * height];

	pixels = new float [ width * height * channel];
		
	for(int i = 0; i < width * height * channel; ++i)
		pixels[i] = 0;

	//UpdateRendering();

	modelName = argv[2];

	SetupGlutDisplay(argc, argv);
	glutMainLoop();

	//free(pixels);

	return 0;
}
