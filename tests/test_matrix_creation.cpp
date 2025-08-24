#include"catch_amalgamated.hpp"
#include"matrixlib/matrix.hpp"

TEST_CASE("Matrix is created with true dimensions", "[creations][copy]")
{
	matrixlib::Matrix original(3,4);

	REQUIRE(original.getRow() == 3);
	REQUIRE(original.getColumn() == 4);

	for(unsigned short i = 0; i < original.getRow(); ++i)
	{
		for(unsigned short j = 0; j < original.getColumn(); ++j)
		{
			original(i,j) = (6*i)/(j+1.0);
		}
	}

	matrixlib::Matrix copy(original);

	REQUIRE(original.getRow() == copy.getRow());
	REQUIRE(original.getColumn() == copy.getColumn());

	for(unsigned short i = 0; i < original.getRow(); ++i)
	{
		for(unsigned short j = 0; j < original.getColumn(); ++j)
		{
			CAPTURE(i,j);
			REQUIRE(copy(i,j) == Catch::Approx(original(i,j)));
		}
	}
} 