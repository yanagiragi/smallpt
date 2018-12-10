#ifndef MODEL_HPP
#define MODEL_HPP
#define TINYOBJLOADER_IMPLEMENTATION

#include <vector>
#include <map>
#include <iostream>
#include "Vec.hpp"
#include <tinyobjloader/tiny_obj_loader.h>

#include "Texture.hpp"

namespace smallPT
{
	class Model
	{
		public:
			std::vector<tinyobj::shape_t> shapes;
			
			tinyobj::attrib_t attrib;
			std::vector<tinyobj::material_t> materials;
			std::map<std::string, unsigned char*> textures;
			std::string warn;
			std::string err;		
			
			std::vector<Vec> positions;
			std::vector<Vec> normals;
			std::vector<Vec> texcoords;
			bool hasUV;

			static bool FileExists(const std::string& abs_filename) {
				bool ret;
				FILE* fp = fopen(abs_filename.c_str(), "rb");
				if (fp) {
					ret = true;
					fclose(fp);
				} else {
					ret = false;
				}

				return ret;
			}

			static std::string GetBaseDir(const std::string& filepath) {
				if (filepath.find_last_of("/\\") != std::string::npos)
					return filepath.substr(0, filepath.find_last_of("/\\"));
				return "";
			}

			void LoadMaterial(const char* filename)
			{
				printf("Found %d Materials\n",  materials.size());
				
				// Load diffuse textures
				{
					std::string base_dir = GetBaseDir(std::string(filename));
					for (size_t m = 0; m < materials.size(); m++) {
						
						tinyobj::material_t* mp = &materials[m];

						if (mp->diffuse_texname.length() > 0) {
							// Only load the texture if it is not already loaded
							if (textures.find(mp->diffuse_texname) == textures.end()) {
								int w, h;
								int comp;

								std::string texture_filename = mp->diffuse_texname;
								std::cout << "load texture: " << texture_filename << std::endl;
								if (!FileExists(texture_filename)) {
									// Append base dir.
									texture_filename = base_dir + mp->diffuse_texname;
									if (!FileExists(texture_filename)) {
										std::cerr << "Unable to find file: " << mp->diffuse_texname
													<< std::endl;
										exit(1);
									}
								}

								unsigned char* image = stbi_load(texture_filename.c_str(), &w, &h, &comp, STBI_default);
								if (!image) {
									std::cerr << "Unable to load texture: " << texture_filename
											<< std::endl;
									exit(1);
								}
								std::cout << "Loaded texture: " << texture_filename << ", w = " << w
											<< ", h = " << h << ", comp = " << comp << std::endl;

								textures.insert(std::make_pair(mp->diffuse_texname, image));
								//stbi_image_free(image);								
							}
						}
					}
				}
			}

			void LoadModelOld(const char* filename)
			{
				/*bool result = tinyobj::LoadObj(shapes, materials, err, filename);

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
				}*/
			}

			Model(const char* filename)
			{
				positions = std::vector<Vec>();

				bool result = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filename);
				
				hasUV = attrib.texcoords.size() > 0;

				// Loop over shapes
				for (size_t s = 0; s < shapes.size(); s++) {
					// Loop over faces(polygon)
					size_t index_offset = 0;
					for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
						int fv = shapes[s].mesh.num_face_vertices[f];
						// Loop over vertices in the face.
						for (size_t v = 0; v < fv; v++) {
							// access to vertex
							tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
							tinyobj::real_t vx = attrib.vertices[3*idx.vertex_index+0];
							tinyobj::real_t vy = attrib.vertices[3*idx.vertex_index+1];
							tinyobj::real_t vz = attrib.vertices[3*idx.vertex_index+2];
							tinyobj::real_t nx = attrib.normals[3*idx.normal_index+0];
							tinyobj::real_t ny = attrib.normals[3*idx.normal_index+1];
							tinyobj::real_t nz = attrib.normals[3*idx.normal_index+2];
							
							Vec pos = Vec(vx, vy, vz);
							positions.push_back(pos);

							Vec normal = Vec(nx, ny, nz);
							normals.push_back(normal);

							if (hasUV) {
								tinyobj::real_t tx = attrib.texcoords[2*idx.texcoord_index+0];
								tinyobj::real_t ty = attrib.texcoords[2*idx.texcoord_index+1];
								Vec uv = Vec(tx, ty, 0);
								texcoords.push_back(uv);
							}
							// Optional: vertex colors
							// tinyobj::real_t red = attrib.colors[3*idx.vertex_index+0];
							// tinyobj::real_t green = attrib.colors[3*idx.vertex_index+1];
							// tinyobj::real_t blue = attrib.colors[3*idx.vertex_index+2];
						}
						index_offset += fv;

						// per-face material
						// shapes[s].mesh.material_ids[f];
					}
				}

				// LoadMaterial(filename);
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
