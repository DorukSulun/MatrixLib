#include"catch_amalgamated.hpp"
#include"matrixlib/matrix.hpp"

TEST_CASE("Matrix row operations","[row_ops]")
{
	SECTION("swapRows")
	{
		matrixlib::Matrix matrix(2,4);
		for(unsigned short i = 0; i < matrix.getRow(); ++i)
		{
			for(unsigned short j = 0; j < matrix.getColumn(); ++j) matrix(i,j) = (22.0/7.0)*j+i;
		}
		matrixlib::Matrix copy(matrix);
		matrix.swapRows(0,1);

		for(unsigned short i = 0; i < matrix.getColumn(); ++i)
		{	
			CAPTURE(i);
			REQUIRE(matrix(1,i) == Catch::Approx(copy(0,i)));
			REQUIRE(matrix(0,i) == Catch::Approx(copy(1,i)));
		}
	}

	SECTION("scaleRow")
	{
		matrixlib::Matrix matrix(1,4);
		matrix(0,0) = 0.0; matrix(0,1) = 1.0; matrix(0,2) = 2.0; matrix(0,3) = -2.0;

		REQUIRE_THROWS_AS(matrix.scaleRow(10,2.7),std::out_of_range);

		matrix.scaleRow(0,5.0);
		REQUIRE(matrix(0,0) == Catch::Approx(0.0));
		REQUIRE(matrix(0,1) == Catch::Approx(5.0));
		REQUIRE(matrix(0,2) == Catch::Approx(10.0));
		REQUIRE(matrix(0,3) == Catch::Approx(-10.0));
	}

	SECTION("addScaledRow")
	{
		matrixlib::Matrix matrix(2,3);

		REQUIRE_THROWS_AS(matrix.addScaledRow(0,5,1.0),std::out_of_range);
		REQUIRE_THROWS_AS(matrix.addScaledRow(1,1,2.0),std::invalid_argument);

		for(unsigned short i = 0; i < matrix.getColumn(); ++i)
		{
			matrix(0,i) = i*1.0;
			matrix(1,i) = i*2.0;
		}
		matrix.addScaledRow(0,1,2.0);

		for (unsigned short i = 0; i < matrix.getColumn(); ++i)
	    {
	        double expected = 2.0 * i + 2.0 * i; 
	        REQUIRE(matrix(1,i) == Catch::Approx(expected));
	    }
	}
}