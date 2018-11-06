 // smallpt, a Path Tracer by Kevin Beason, 2008 
 // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt 
 // Remove "-fopenmp" for g++ version < 4.2 

#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>

#include "Scene.hpp"
#include "Render.hpp"
#include "Utils.hpp"
#include "Display.hpp"

#include "Config.hpp"

void sigintHandler(int arg) 
{
	if (SIGINT == arg)
	{
		printf("\ninterrupted\n");
		smallPT::ReleaseConfig();
		exit(0);
	}
}


int main(int argc, char **argv)
{
	signal(SIGINT, sigintHandler);

	smallPT::InitConfig(
		800, // width
		600, // height
		atoi(argv[1]), // spp
		argv[2] // model loading name
	);

	//	SetupGlutDisplay(argc, argv);
	//	glutMainLoop();

	smallPT::UpdateRendering();
	std::ostringstream outputFilenameStringStream;
	outputFilenameStringStream << smallPT::SaveImageNamePrefix << smallPT::limitSpp << "spp";

	std::string outputPpmFilename = outputFilenameStringStream.str() + ".ppm";
	std::string outputPngFilename = outputFilenameStringStream.str() + ".png";

	// Save PPM
	smallPT::SavePPM(outputPpmFilename, smallPT::width, smallPT::height, smallPT::output);

	smallPT::ReleaseConfig();

	return 0;
}