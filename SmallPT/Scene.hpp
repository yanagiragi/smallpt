#ifndef _SCENE_HPP
#define _SCENE_HPP

#include <vector>

#include "Sphere.hpp"
#include "Triangle.hpp"
#include "Model.hpp"
#include "Utils.hpp"

namespace smallPT {

	class Scene
	{
	public:
		// Variables
		std::vector<Sphere> Spheres;
		std::vector<Triangle> Triangles;
		// std::vector<cv::Mat> Textures;

		void LoadScene()
		{
			// always load sphere before loading model
			LoadSphere();

			double x = 25;
			double y = 1;
			double z = 30;

			//obj.translate(Vec(50 + x, y, z));
			Vec bunnyPos = Vec(50, 0, 81.6);

			Material SimpleDiffuseMat(Vec(), Vec(.999, .999, .999), DIFF);
			//char* testModel = "models/bunnyLowWithUV.obj";
			char* testModel = "null";
			//char* testModel  = "models/Plane.obj";
			//bunnyPos = Vec(50, 20, 81.6);

			if (strncmp(testModel, "null", 4) != 0)
				LoadModel(testModel, SimpleDiffuseMat, bunnyPos);

			Vec lightPos = Vec(50, 81.6 - 16.5, 81.6);
			lightPos = Vec(50, 20, 81.6);
			//		lightPos = Vec(50, 10, 5);
			lightPos = Vec(50, 10, 81.6);

			Vec emi = Vec(4, 4, 4) * 100;
			Vec col = Vec(1, 1, 1);
			Material lightMat(emi, col, DIFF);

			char* lightObj = "models/Plane.obj";
			LoadModel(lightObj, lightMat, lightPos);
			// Almost equalilent to
			// Spheres.push_back(Sphere(12, Vec(50, 81.6 - 16.5, 81.6), Vec(4, 4, 4) * 2, Vec(), DIFF)); //6.Lite

			printf(" Sphere\t = %d\n", (int)Spheres.size());
			printf(" Triangle = %d\n", (int)Triangles.size());
		}

		// Defined Functions	
		void LoadSphere()
		{
			Spheres.push_back(Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF)); //0.Left
			Spheres.push_back(Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF)); //1.Rght
			Spheres.push_back(Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF)); //2.Back
			Spheres.push_back(Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), DIFF)); //3.Frnt
			Spheres.push_back(Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF)); //4.Botm
			Spheres.push_back(Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF)); //5.Top

																										 // For Testing
																										 //Spheres.push_back(Sphere(0.05, Vec(50, 81.6, 81.6), Vec(4, 4, 4) * 100000, Vec(), DIFF)); //6.Lite

			Spheres.push_back(Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1)*.999, DIFF)); //7.Mirr

																								// Original
			Spheres.push_back(Sphere(1.5, Vec(50, 81.6 - 16.5, 81.6), Vec(4, 4, 4) * 100, Vec(), DIFF)); //6.Lite

																										 // hard shadow
																										 // Spheres.push_back(Sphere(0.05, Vec(50, 81.6 - 16.5, 81.6), Vec(4, 4, 4) * 100000, Vec(), DIFF)); //6.Lite

																										 // soft shadow
																										 // Spheres.push_back(Sphere(12, Vec(50, 81.6 - 16.5, 81.6), Vec(4, 4, 4) * 2, Vec(), DIFF)); //6.Lite

																										 //Spheres.push_back(Sphere(16.5,Vec(27,16.5,47),			Vec(),			Vec(1,1,1)*.999,	SPEC)); //7.Mirr
																										 //Spheres.push_back(Sphere(16.5,Vec(73,16.5,78),			Vec(),			Vec(1,1,1)*.999,	REFR)); //8.Glas

		}

		void LoadModel(const char *filename, Material mat, Vec origin)
		{
			Model obj = Model(filename);

			obj.translate(origin);

			//cv::Mat image1 = LoadImage("models/TexturesCom_TilesOrnate0043_1_seamless_S.jpg");
			//cv::Mat image1 = LoadImage("models/Untitled.png");
			//Textures.push_back(image1);

			for (int i = 0; i < obj.positions.size(); i += 3) {
				Triangle tri = Triangle(
					obj.positions[i + 0],
					obj.positions[i + 1],
					obj.positions[i + 2],
					obj.normals[i + 0],
					obj.normals[i + 1],
					obj.normals[i + 2],
					obj.hasUV,
					mat);
				if (obj.hasUV) {
					tri.SetUV(obj.texcoords[i + 0], obj.texcoords[i + 1], obj.texcoords[i + 2]);
				}
				Triangles.push_back(tri);
			}
		}

		// get min intersect distance & SphereId (the intersected sphere that is not being occlude from other sphere)
		inline bool intersect(const Ray &ray, double &t, int &id)
		{
			double infinite = t = 1e20;
			double dis;

			for (int i = 0; i < Spheres.size(); ++i) {
				dis = Spheres[i].intersect(ray);
				if (dis > 0 && dis < t) {
					t = dis;
					id = i;
				}
			}

			for (int i = 0; i < Triangles.size(); ++i) {
				dis = Triangles[i].intersect(ray);
				if (dis > 0 && dis < t) {
					t = dis;
					id = Spheres.size() + i;
				}
			}

			return t < infinite;
		}
	};
}

#endif // !_SCENE_HPP