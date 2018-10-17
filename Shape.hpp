#ifndef SHAPE_HPP
#define SHAPE_HPP

#include "Ray.hpp"
#include "Vec.hpp"

const double epsilon = 1e-4;

enum ReflectType
{
	DIFF, SPEC, REFR
};

class Shape
{
	public:

		Shape(Vec emi, Vec col, ReflectType refl) : emission(emi), color(col), reflectType(refl) {}

		Vec emission;

		Vec color;

		ReflectType reflectType;

};

#endif // !SPHERE_HPP
