#ifndef _MATERIAL_HPP
#define _MATERIAL_HPP

#include "Vec.hpp"

namespace smallPT
{

	enum ReflectType
	{
		DIFF, SPEC, REFR
	};

	class Material
	{
	public:

		Vec emission;

		Vec color;

		ReflectType reflectType;

		Material(Vec emi, Vec col, ReflectType refl) :
			emission(emi), color(col), reflectType(refl) {}
	};

}

#endif