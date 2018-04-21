

#ifndef 969d7bbbb7bc427f9813bb0bfc7e0abe
#define 969d7bbbb7bc427f9813bb0bfc7e0abe
#include "entity.hpp"

namespace hl::mesh
{
	template<int DIM_WORLD, int DIM_MESH>
	class Mesh
	{
	public:
		static constexpr int dimWorld = DIM_WORLD;
		static constexpr int dimMesh = DIM_MESH;
	
		template<int CODIM>
		using entityIterator_t = Entity<dimWorld, dimMesh, CODIM>*;


	};
}




#endif // 969d7bbbb7bc427f9813bb0bfc7e0abe


