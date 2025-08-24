#include"catch_amalgamated.hpp"
#include"matrixlib/matrix.hpp"

TEST_CASE("Arithmetic operations","[arithmetic]")
{
	SECTION("Addition")
	{
		matrixlib::Matrix A(2,2);
		matrixlib::Matrix B(2,2);

		A(0,0) = 1.0; A(0,1) = 2.0;
		A(1,0) = 3.0; A(1,1) = 0.9;

		B(0,0) = 3.0; B(0,1) = 2.0;
		B(1,0) = 1.0; B(1,1) = 3.1;

		matrixlib::Matrix result = A + B;

		for(unsigned short i = 0; i < result.getRow(); ++i)
		{
			for(unsigned short j = 0; j < result.getColumn(); ++j)
			{
				CAPTURE(i,j);
				REQUIRE(result(i,j) == Catch::Approx(A(i,j) + B(i,j)));
			}
		}

		matrixlib::Matrix C(2,3);
		matrixlib::Matrix D(3,2);

		REQUIRE_THROWS_AS(C + D, std::invalid_argument);
	}

	SECTION("Subtraction")
	{
		matrixlib::Matrix A(2,2);
		matrixlib::Matrix B(2,2);

		A(0,0) = 1.0; A(0,1) = 2.0;
		A(1,0) = 3.0; A(1,1) = 0.9;

		B(0,0) = 3.0; B(0,1) = 2.0;
		B(1,0) = 1.0; B(1,1) = 3.1;

		matrixlib::Matrix result = A - B;

		for(unsigned short i = 0; i < result.getRow(); ++i)
		{
			for(unsigned short j = 0; j < result.getColumn(); ++j)
			{
				CAPTURE(i,j);
				REQUIRE(result(i,j) == Catch::Approx(A(i,j) - B(i,j)));
			}
		}

		matrixlib::Matrix C(3,3);

		REQUIRE_THROWS_AS(A - C, std::invalid_argument);
	}

	SECTION("Multiplication")
	{
		matrixlib::Matrix A(2,1);
		matrixlib::Matrix B(1,2);

		A(0,0) = 2.0;		B(0,0) = 5.0; B(0,1) = 10.0;
		A(1,0) = 3.0;

		matrixlib::Matrix result = A*B;		// result 2x2

		for(unsigned short i = 0; i < result.getColumn(); ++i)
		{
			REQUIRE(result(0,i) == Catch::Approx((i+1)*10.0));
			REQUIRE(result(1,i) == Catch::Approx((i+1)*15.0));
		}

		REQUIRE_THROWS_AS(A*result,std::invalid_argument);
	}
}