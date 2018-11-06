#ifndef _CONFIG_HPP
#define _CONFIG_HPP

#include <stdlib.h>
#include <string>

#include "Vec.hpp"
//#include "Scene.hpp"

namespace smallPT
{   
	// GLfloat[3 * width * height] for glut display
	// not exactly same to output
	// short speaking , pixels = flipCoordinate(gammaCorrection(output))
	float *pixels;
	int channel = 3;
    
	// Vec[width * height] for path tracing results
    Vec *output;

	// resolution settings
    int width, height;

	// max spp
    int limitSpp = 4;
    int currentSpp = 0;
    
    Scene MainScene = Scene();
    
    char *modelName;
    std::string SaveImageNamePrefix;

    // static unsigned int *seeds;

	// thresholds & constants
	const double epsilon = 1e-4;
	const double pi = 3.1415926535;
	const double reciprocalPi = 1 / pi;

    void InitConfig(int w, int h, int spp, char *modelName)
    {
        width = w;
        height = h;
        limitSpp = spp > 1 ? spp : 0;
        
        pixels = new float [ width * height * channel];
        for(int i = 0; i < width * height * channel; ++i)
    		pixels[i] = 0;
        
        output = new Vec[width * height];

        /*seeds = new unsigned int[width * height * channel * 2];
        for (int i = 0; i < width * height * channel * 2; ++i) {
            seeds[i] = rand();
            if (seeds[i] < 2)
                seeds[i] = 2;
        }*/

        if (strncmp(modelName, "null", 4) != 0)
		    modelName = modelName;
        else
            modelName = "null";

        // Generate Filename
        std::string modelFileName(modelName);

        if (modelFileName.find_first_of("/") != -1) {
            // 4 for ".obj".length
            modelFileName = modelFileName.substr(modelFileName.find_first_of("/") + 1, modelFileName.length() - modelFileName.find_first_of("/") - 1 - 4);
        }
        else {
            modelFileName = modelFileName.substr(0, modelFileName.length() - 4);
        }
        
        std::ostringstream outputFilenameStringStream;
        outputFilenameStringStream << "images/" << "output_" << modelFileName << "_" ;

        SaveImageNamePrefix = outputFilenameStringStream.str();     
        
        MainScene.LoadScene();
    }

    void ReleaseConfig()
    {
        free(pixels);
        free(output);
        // free(seeds);
    }
}

#endif