#ifndef  _UTILS_HPP
#define _UTILS_HPP

#include <math.h> 
#include <string>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "Vec.hpp"

#ifdef _MSC_VER
	#include <random>
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distr(0.0, 1.0);
	inline double erand48(unsigned short *x) {return distr(generator);}
#endif // _MSC_VER

const double pi = 3.1415926535;
const double reciprocalPi = 1 / pi;

// clamp between 0 and 1
inline double clamp(double x) { return x > 1 ? 1 : x < 0 ? 0 : x; }

// Convert float to int for ppm format, Note that we also apply gamma correction.
//inline int toInt(double x) { return  int(pow(clamp(x), 1.0f / 2.2f) * 255 + 0.5); }
inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }

void SavePPM(std::string outputFileName, int width, int height, Vec *outputData)
{
	FILE *f = fopen(outputFileName.c_str(), "w");
	fprintf(f, "P3\n%d %d\n%d\n", width, height, 255); // store value in [0, 255]
	for (int i = 0; i < width * height; ++i) {
		fprintf(f, "%d %d %d ", toInt(outputData[i].x), toInt(outputData[i].y), toInt(outputData[i].z));
	}
	fclose(f);
	// From smallpt:
	// Instead of fclose()ing the file, I exploit the C++ standard which calls return(0) implicitly, which in turn calls exit(0), which flushes and closes open files. 
}

void ConvertPpmToMat(int width, int height, Vec *ppmData, cv::Mat &matData)
{
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			int index = (height - i - 1) * width + j;
			matData.at<cv::Vec3b>(height - i - 1, width - j - 1)[0] = toInt(ppmData[index].x);
			matData.at<cv::Vec3b>(height - i - 1, width - j - 1)[1] = toInt(ppmData[index].y);
			matData.at<cv::Vec3b>(height - i - 1, width - j - 1)[2] = toInt(ppmData[index].z);
		}
	}
}

void SaveResult(char **argv, int width, int height, int spp, Vec* ppmData)
{
	// Generate Filename
	std::string modelFileName(argv[2]);

	if (modelFileName.find_first_of("/") != -1) {
		// 4 for ".obj".length
		modelFileName = modelFileName.substr(modelFileName.find_first_of("/") + 1, modelFileName.length() - modelFileName.find_first_of("/") - 1 - 4);
	}
	else {
		modelFileName = modelFileName.substr(0, modelFileName.length() - 4);
	}
	
	std::ostringstream outputFilenameStringStream;
	outputFilenameStringStream << "images/" << "output_" << modelFileName << "_" << spp << "spp";

	std::string outputPpmFilename = outputFilenameStringStream.str() + ".ppm";
	std::string outputPngFilename = outputFilenameStringStream.str() + ".png";

	// Save PPM
	SavePPM(outputPpmFilename, width, height, ppmData);

	// Save PNG
	cv::Mat Image(height, width, CV_8UC3);
	ConvertPpmToMat(width, height, ppmData, Image);
	cv::imwrite(outputPngFilename, Image);
}

#endif