#ifndef MODEL_HPP
#define MODEL_HPP
#define TINYOBJLOADER_IMPLEMENTATION

#include <vector>
#include "Vec.hpp"
#include <tinyobjloader/tiny_obj_loader.h>

class Model
{
	public:

		std::vector<tinyobj::shape_t> shapes;
		std::vector<tinyobj::material_t> materials;
		std::string err;
		std::vector<Vec> positions;

		Model(const char* filename)
		{
			positions = std::vector<Vec>();

			bool result = tinyobj::LoadObj(shapes, materials, err, filename);

			for (size_t i = 0; i < shapes.size(); ++i)
			{
				for (size_t f = 0; f < shapes[i].mesh.indices.size() / 3; f++) {

					int indices0 = shapes[i].mesh.indices[f * 3 + 0];
					int indices1 = shapes[i].mesh.indices[f * 3 + 1];
					int indices2 = shapes[i].mesh.indices[f * 3 + 2];

					Vec pos = Vec(
						shapes[i].mesh.positions[indices0 * 3 + 0],
						shapes[i].mesh.positions[indices0 * 3 + 1],
						shapes[i].mesh.positions[indices0 * 3 + 2]
					);
					positions.push_back(pos);
					
					pos = Vec(
						shapes[i].mesh.positions[indices1 * 3 + 0],
						shapes[i].mesh.positions[indices1 * 3 + 1],
						shapes[i].mesh.positions[indices1 * 3 + 2]
					);
					positions.push_back(pos);
					
					pos = Vec(
						shapes[i].mesh.positions[indices2 * 3 + 0],
						shapes[i].mesh.positions[indices2 * 3 + 1],
						shapes[i].mesh.positions[indices2 * 3 + 2]
					);
					positions.push_back(pos);
				}
			}
		}

		void translate(Vec offset) {
			for(int i = 0; i < positions.size(); ++i)
				positions[i] = positions[i] + offset;
		}

		void scale(double factor){
			for(int i = 0; i < positions.size(); ++i)
				positions[i] = positions[i] * factor;
		}
};

#endif // !SPHERE_HPP