#ifndef _SCENE_HPP
#define _SCENE_HPP

#include <vector>

#include "Sphere.hpp"
#include "Triangle.hpp"
#include "Model.hpp"
#include "Utils.hpp"
#include "Texture.hpp"

namespace smallPT {

	class Scene
	{
	public:
		// Variables
		std::vector<Sphere> Spheres;
		std::vector<Triangle> Triangles;
		
		std::string m_modelName;
		void LoadScene(std::string modelName)
		{
			m_modelName = modelName;
			LoadScene();
		}		

		void LoadScene()
		{
			// always load sphere before loading model
			LoadSphere();

			double x = 25;
			double y = 1;
			double z = 30;

			Vec bunnyPos = Vec(50, -1, 81.6);

			Material SimpleDiffuseMat(Vec(), Vec(.0, .999, .999), DIFF);

			if(m_modelName == "haku")
			{
				Material SkinMat(Vec(), Vec(255.0/255.0, 230.0/255.0, 196.0/255.0), DIFF);
				Material WhiteMat(Vec(), Vec(.999, .999, .999), DIFF);
				Material RedMat(Vec(), Vec(.999, .0, .0), DIFF);
				Material BlackMat(Vec(), Vec(0, 0, 0), DIFF);
				Material HairMat(Vec(), Vec(229.0/255.0, 196.0/255.0, 255.0/255.0), DIFF);
				Material JewelMat(Vec(), Vec(0.999, 0.999, 0.999), REFR);
				Material DressMat(Vec(), Vec(63.0 / 255.0, 61.0 / 255.0, 74.0/255.0), DIFF);
				//Material StockingMat(Vec(), Vec(.067, .067, .067), DIFF);
				Material MouthMat(Vec(), Vec(239.0 / 255.0, 167.0 / 255.0, 176.0/255.0), DIFF);
				
				LoadModel("models/haku/body2.obj", SkinMat, bunnyPos);
				LoadModel("models/haku/white.obj", WhiteMat, bunnyPos);
				LoadModel("models/haku/hair.obj", HairMat, bunnyPos);
				LoadModel("models/haku/earRingUpper.obj", HairMat, bunnyPos);
				LoadModel("models/haku/earRingLower.obj", JewelMat, bunnyPos);
				LoadModel("models/haku/earRingLower.obj", JewelMat, bunnyPos);
				LoadModel("models/haku/eyebrow.obj", JewelMat, bunnyPos);
				LoadModel("models/haku/nicebody.obj", DressMat, bunnyPos);
				LoadModel("models/haku/Gloth.obj", DressMat, bunnyPos);
				LoadModel("models/haku/red.obj", WhiteMat, bunnyPos);
				LoadModel("models/haku/mouth.obj", MouthMat, bunnyPos);
				LoadModel("models/haku/shoe.obj", DressMat, bunnyPos);
			}
			else if (m_modelName != "null")
				LoadModel(m_modelName.c_str(), SimpleDiffuseMat, bunnyPos);
			
			Vec lightPos = Vec(50, 81.6 - 16.5, 81.6);
			// lightPos = Vec(50, 20, 81.6);
			// lightPos = Vec(50, 30, 81.6);
			// lightPos = Vec(50, 10, 5);
			lightPos = Vec(50, 81.3, 81.6);

			Vec emi = Vec(4, 4, 4);
			Vec col = Vec(1, 1, 1);
			Material lightMat(emi, col, DIFF);

			const std::string lightObj("models/HoriPlaneLarge.obj");
			LoadModel(lightObj.c_str(), lightMat, lightPos);
			// Almost equalilent to
			// Spheres.push_back(Sphere(12, Vec(50, 81.6 - 16.5, 81.6), Vec(4, 4, 4) * 2, Vec(), DIFF)); //6.Lite

			printf(" Sphere\t = %d\n", (int)Spheres.size());
			printf(" Triangle = %d\n", (int)Triangles.size());
		}

		// Defined Functions	
		void LoadSphere()
		{
			// Sphere(double rad, Vec ori, Vec emi, Vec col, ReflectType mat)
			Spheres.push_back(Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF)); //0.Left
			Spheres.push_back(Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF)); //1.Rght
			Spheres.push_back(Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF)); //2.Back
			Spheres.push_back(Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), DIFF)); //3.Frnt
			Spheres.push_back(Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF)); //4.Botm
			Spheres.push_back(Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF)); //5.Top

			// For Testing
			// Spheres.push_back(Sphere(12, Vec(50, 81.6 - 16.5, 81.6), Vec(4, 4, 4) * 2, Vec(), DIFF)); //6.Lite

			// Spheres.push_back(Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1)*.999, DIFF)); //7.Mirr

			// Original
			//Spheres.push_back(Sphere(1.5, Vec(50, 81.6 - 16.5, 81.6), Vec(4, 4, 4) * 100, Vec(), DIFF)); //6.Lite

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
