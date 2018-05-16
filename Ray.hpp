#ifndef RAY_HPP
#define RAY_HPP

#include "Vec.hpp"

class Ray {
	
	public:
		
		Vec origin, direction;

		Ray(Vec ori, Vec dir)
		{
			origin = ori;
			direction = dir;
		}
};

#endif // !RAY_HPP