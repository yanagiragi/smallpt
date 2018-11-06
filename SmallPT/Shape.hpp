#ifndef SHAPE_HPP
#define SHAPE_HPP

#include "Ray.hpp"
#include "Vec.hpp"
#include "Material.hpp"

namespace smallPT
{
	class Shape
	{
		public:

		const double epsilon = 1e-4;

		Material material;

		Shape(Vec emi, Vec col, ReflectType refl) :
			material(emi, col, refl) {}
		Shape(Material mat) :
			material(mat) {}
	};

}

#endif // !SPHERE_HPP
