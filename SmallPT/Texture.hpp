#ifndef TEXTURE_HPP
#define TEXTURE_HPP

#include <vector>
#include <iostream>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>


#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

namespace smallPT
{
	class Texture
	{
		public:
			unsigned char* image;
            int width, height;
            int comp;
			
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

			Texture(const char* texture_filename)
			{
                if (!FileExists(texture_filename)) {
                    std::cerr << "Unable to find file: " << texture_filename << std::endl;
                    exit(1);
                }                

                this->image = stbi_load(texture_filename, &width, &height, &comp, STBI_default);
                if (!image) {
                    std::cerr << "Unable to load texture: " << texture_filename << std::endl;
                    exit(1);
                }
			}
            ~Texture()
            {
                printf("Free\n");
                stbi_image_free(image);
            }
	};
}
#endif // !SPHERE_HPP
