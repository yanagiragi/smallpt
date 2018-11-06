#ifndef  _UTILS_HPP
#define _UTILS_HPP

#include <math.h> 

//#include <opencv2/core/core.hpp>
//#include <opencv2/highgui/highgui.hpp>
#include <string>


#ifdef _MSC_VER
#include <random>
std::default_random_engine generator;
std::uniform_real_distribution<double> distr(0.0, 1.0);
inline double erand48(unsigned short *x) { return distr(generator); }
#else
#include <opencv2/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#endif // _MSC_VER


#include "Vec.hpp"

namespace smallPT {

	// clamp between 0 and 1
	inline double clamp(double x) { return x > 1 ? 1 : x < 0 ? 0 : x; }

	// Convert float to int for ppm format, Note that we also apply gamma correction.
	//inline int toInt(double x) { return  int(pow(clamp(x), 1.0f / 2.2f) * 255 + 0.5); }
	inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }
	inline double gammaCorrection(double x) { return pow(clamp(x), 1 / 2.2); }

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

	//void ConvertPpmToMat(int width, int height, Vec *ppmData, cv::Mat &matData)
	//{
	//	for (int i = 0; i < height; ++i) {
	//		for (int j = 0; j < width; ++j) {
	//			int index = (height - i - 1) * width + j;
	//			matData.at<cv::Vec3b>(height - i - 1, width - j - 1)[0] = toInt(ppmData[index].x);
	//			matData.at<cv::Vec3b>(height - i - 1, width - j - 1)[1] = toInt(ppmData[index].y);
	//			matData.at<cv::Vec3b>(height - i - 1, width - j - 1)[2] = toInt(ppmData[index].z);
	//		}
	//	}
	//}

	//void SaveResult()
	//{
	//	// Generate Filename
	//	std::ostringstream outputFilenameStringStream;
	//	outputFilenameStringStream << SaveImageNamePrefix << currentSpp << "spp";

	//	std::string outputPpmFilename = outputFilenameStringStream.str() + ".ppm";
	//	std::string outputPngFilename = outputFilenameStringStream.str() + ".png";

	//	// Save PPM
	//	SavePPM(outputPpmFilename, output);

	//	// Save PNG
	//	//cv::Mat Image(height, width, CV_8UC3);
	//	//ConvertPpmToMat(width, height, ppmData, Image);
	//	//cv::imwrite(outputPngFilename, Image);
	//}

	//cv::Mat LoadImage(const std::string &filename)
	//{
	//	cv::Mat image = cv::imread(filename, CV_LOAD_IMAGE_UNCHANGED);
	//	return image;
	//}
}

#endif
