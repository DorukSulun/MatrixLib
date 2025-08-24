#include"catch_amalgamated.hpp"
#include"matrixlib/matrix.hpp"

TEST_CASE("Special Matrices","[specialMatrices]")
{
	SECTION("Identity")
	{
		matrixlib::Matrix I = matrixlib::Matrix::identity(4);
		const unsigned short n = I.getRow();

		for(unsigned short i = 0; i < n; ++i)
		{
			for(unsigned short j = 0; j < n; ++j)
			{
				if(i == j) REQUIRE(I(i,j) == Catch::Approx(1.0));
				else REQUIRE(I(i,j) == Catch::Approx(0.0));
			}
		}
	}

	SECTION("Zeros")
	{
		matrixlib::Matrix Z = matrixlib::Matrix::zeros(3,9);
		REQUIRE(Z.getRow() == 3);
		REQUIRE(Z.getColumn() == 9);

		for(unsigned short i = 0; i < Z.getRow(); ++i)
		{
			for(unsigned short j = 0; j < Z.getColumn(); ++j) REQUIRE(Z(i,j) == Catch::Approx(0.0));
		}		
	}

	SECTION("Ones")
	{
		matrixlib::Matrix O = matrixlib::Matrix::ones(3,9);
		REQUIRE(O.getRow() == 3);
		REQUIRE(O.getColumn() == 9);

		for(unsigned short i = 0; i < O.getRow(); ++i)
		{
			for(unsigned short j = 0; j < O.getColumn(); ++j)
			{
				REQUIRE(O(i,j) == Catch::Approx(1.0));
			}
		} 
	}

	SECTION("Diagonal")
	{
		matrixlib::Matrix D = matrixlib::Matrix::diagonal(6);
		REQUIRE(D.getRow() == 6);
		REQUIRE(D.getColumn() == 6);
		
		for(unsigned short i = 0; i < D.getRow(); ++i)
		{
			for(unsigned short j = 0; j < D.getColumn(); ++j)
			{
				if(i == j) REQUIRE(D(i,j) == Catch::Approx(1.0));
				else REQUIRE(D(i,j) == Catch::Approx(0.0));
			}
		}
	}

	SECTION("Scalar matrix")
	{
		double val = 2.3571113;
		matrixlib::Matrix S = matrixlib::Matrix::scalar(6,val);

		REQUIRE(S.getRow() == 6);
		REQUIRE(S.getColumn() == 6);

		for(unsigned short i = 0; i < S.getRow(); ++i)
		{
			for(unsigned short j = 0; j < S.getColumn(); ++j)
			{
				if(i == j) REQUIRE(S(i,j) == Catch::Approx(val));
				else REQUIRE(S(i,j) == Catch::Approx(0.0));
			}
		}
	}

	SECTION("Exchange")
	{
		matrixlib::Matrix E = matrixlib::Matrix::exchange(6);
		REQUIRE(E.getRow() == 6);
		REQUIRE(E.getColumn() == 6);

		for(unsigned short i = 0; i < E.getRow(); ++i)
		{
			for(unsigned short j = 0; j < E.getColumn(); ++j)
			{
				if(j == E.getRow()-1-i) REQUIRE(E(i,j) == Catch::Approx(1.0));
				else REQUIRE(E(i,j) == Catch::Approx(0.0));
			}
		}		
	}
}