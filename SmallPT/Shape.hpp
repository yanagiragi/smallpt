#ifndef SHAPE_HPP
#define SHAPE_HPP

#include "Ray.hpp"
#include "Vec.hpp"
#include "Material.hpp"

const double epsilon = 1e-4;

class Shape
{
	public:

		Material material;

		Shape(Vec emi, Vec col, ReflectType refl) : 
			material(emi, col, refl) {}
};

#endif // !SPHERE_HPP
