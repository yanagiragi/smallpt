#ifndef MODEL_HPP
#define MODEL_HPP
#define TINYOBJLOADER_IMPLEMENTATION

#include <vector>
#include "Vec.hpp"
#include <tinyobjloader/tiny_obj_loader.h>

namespace smallPT
{
	class Model
	{
		public:
			std::vector<tinyobj::shape_t> shapes;
			std::vector<tinyobj::material_t> materials;
			std::string err;
			std::vector<Vec> positions;
			std::vector<Vec> normals;
			std::vector<Vec> texcoords;
			bool hasUV;

			Model(const char* filename)
			{
				positions = std::vector<Vec>();

				bool result = tinyobj::LoadObj(shapes, materials, err, filename);

				for (size_t i = 0; i < shapes.size(); ++i)
				{
					hasUV = shapes[i].mesh.texcoords.size() > 0;

					//bool hasNormal = shapes[i].mesh.normals.size() > 0;

					for (size_t f = 0; f < shapes[i].mesh.indices.size() / 3; f++) {

						for (int j = 0; j < 3; ++j) {
							int indices = shapes[i].mesh.indices[f * 3 + j];

							Vec pos = Vec(
								shapes[i].mesh.positions[indices * 3 + 0],
								shapes[i].mesh.positions[indices * 3 + 1],
								shapes[i].mesh.positions[indices * 3 + 2]
							);
							positions.push_back(pos);

							Vec normal = Vec(
								shapes[i].mesh.normals[indices * 3 + 0],
								shapes[i].mesh.normals[indices * 3 + 1],
								shapes[i].mesh.normals[indices * 3 + 2]
							);
							// printf("%f %f %f\n", normal.x, normal.y, normal.z);
							normals.push_back(normal);

							if (hasUV) {
								Vec uv = Vec(
									shapes[i].mesh.texcoords[indices * 2 + 0],
									shapes[i].mesh.texcoords[indices * 2 + 1],
									0
								);
								texcoords.push_back(uv);
							}

						}
					}
				}
			}

			void translate(Vec offset) {
				for (int i = 0; i < positions.size(); ++i)
					positions[i] = positions[i] + offset;
			}

			void scale(double factor) {
				for (int i = 0; i < positions.size(); ++i)
					positions[i] = positions[i] * factor;
			}
	};
}
#endif // !SPHERE_HPP
