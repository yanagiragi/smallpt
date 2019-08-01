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
		atoi(argv[2]), // spp
		argv[3] // model loading name
	);

	if (std::string(argv[1]) == "glut") {
		smallPT::SetupGlutDisplay(argc, argv);
		glutMainLoop();
	}
	else {
		// std::string(argv[1]) == "console"
		while (smallPT::currentSpp < smallPT::limitSpp) {
			smallPT::UpdateRendering();	
		}
	}

	smallPT::ReleaseConfig();

	return 0;
}