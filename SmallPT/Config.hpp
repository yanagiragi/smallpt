#ifndef _CONFIG_HPP
#define _CONFIG_HPP

#include <stdlib.h>
#include <string>

#include "Vec.hpp"
//#include "Scene.hpp"

namespace globalConfig
{   
    float *pixels;
    
    Vec *output;

    int width;
    int height;
    int limitSpp = 4;

    int channel = 3;
    int currentSpp = 0;
    
    Scene MainScene = Scene();
    
    char *modelName;

    std::string SaveImageNamePrefix;

    static unsigned int *seeds;

    bool InitConfig(int width, int height, int spp, char *modelName)
    {
        globalConfig::width = width;
        globalConfig::height = height;
        globalConfig::limitSpp = spp > 1 ? spp : 0;
        
        globalConfig::pixels = new float [ width * height * channel];
        for(int i = 0; i < width * height * channel; ++i)
    		pixels[i] = 0;
        
        globalConfig::output = new Vec[width * height];

        globalConfig::seeds = new unsigned int[width * height * channel * 2];
        for (int i = 0; i < width * height * channel * 2; ++i) {
            globalConfig::seeds[i] = rand();
            if (globalConfig::seeds[i] < 2)
                globalConfig::seeds[i] = 2;
        }

        if (strncmp(modelName, "null", 4) != 0)
		    globalConfig::modelName = modelName;
        else
            globalConfig::modelName = "null";

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
        free(globalConfig::pixels);
        free(globalConfig::output);
        free(globalConfig::seeds);
    }
}

#endif