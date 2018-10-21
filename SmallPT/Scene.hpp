#ifndef _SCENE_HPP
#define  _SCENE_HPP

#include <vector>

#include "Vec.hpp"
#include "Shape.hpp"
#include "Sphere.hpp"
#include "Triangle.hpp"
#include "Model.hpp"

// Variables
std::vector<Sphere> Spheres;
std::vector<Triangle> Triangles;

// Defined Functions
void LoadScene(char **argv);
void LoadModel(char *filename);
inline bool intersect(const Ray &ray, double &t, int &id);

void LoadScene(char **argv)
{
	if (strncmp(argv[2], "null", 4) != 0)
		LoadModel(argv[2]);

	Spheres.push_back(Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF)); //0.Left
	Spheres.push_back(Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF)); //1.Rght
	Spheres.push_back(Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF)); //2.Back
	Spheres.push_back(Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), DIFF)); //3.Frnt
	Spheres.push_back(Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF)); //4.Botm
	Spheres.push_back(Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF)); //5.Top
	Spheres.push_back(Sphere(1.5, Vec(50, 81.6 - 16.5, 81.6), Vec(4, 4, 4) * 100, Vec(), DIFF)); //6.Lite
	//Spheres.push_back(Sphere(16.5,Vec(27,16.5,47),			Vec(),			Vec(1,1,1)*.999,	SPEC)); //7.Mirr
	//Spheres.push_back(Sphere(16.5,Vec(73,16.5,78),			Vec(),			Vec(1,1,1)*.999,	REFR)); //8.Glas

	printf(" Sphere\t = %d\n", (int)Spheres.size());
	printf(" Triangle = %d\n", (int)Triangles.size());
}

void LoadModel(char *filename)
{
	Model obj = Model(filename);

	double x = 25;
	double y = 1;
	double z = 30;

	obj.translate(Vec(50 + x, y, z));

	for (int i = 0; i < obj.positions.size(); i += 3) {
		Triangle tri = Triangle(obj.positions[i + 0], obj.positions[i + 1], obj.positions[i + 2], Vec(), Vec(.25, .999, .25), DIFF);
		Triangles.push_back(tri);
	}

	/*Model obj2 = Model(filename);
	obj2.translate(Vec(50 - x, y, z));

	for (int i = 0; i < obj2.positions.size(); i += 3)
	{
		Triangle tri = Triangle(obj2.positions[i + 0], obj2.positions[i + 1], obj2.positions[i + 2], Vec(), Vec(.999, .999, .999), DIFF);
		Triangles.push_back(tri);
	}*/
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


#endif // !_SCENE_HPP