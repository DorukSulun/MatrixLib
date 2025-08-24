#include"catch_amalgamated.hpp"
#include"matrixlib/matrix.hpp"
#include"matrixlib/utility.hpp"
#include<cmath>

TEST_CASE("Matrix transformations functions","[transformations]")
{
	SECTION("Transpose")
	{
		matrixlib::Matrix matrix(2,3);

		for(unsigned short i = 0; i < matrix.getRow(); ++i)
		{
			for(unsigned short j = 0; j < matrix.getColumn(); ++j) matrix(i,j) = (i == j) ? (3.14*i+j) : (2.71*j+i);
		}

		matrixlib::Matrix matrixT = matrix.transpose();

		REQUIRE(matrixT.getRow() == matrix.getColumn());
		REQUIRE(matrixT.getColumn() == matrix.getRow());

		for(unsigned short i = 0; i < matrix.getRow(); ++i)
		{
			for(unsigned short j = 0; j < matrix.getColumn(); ++j)
			{
				CAPTURE(i,j);
				REQUIRE(matrix(i,j) == Catch::Approx(matrixT(j,i)));
			}
		}
	}

	SECTION("doUpperTriangular")
	{
		SECTION("Throws if matrix is not square.")
		{
			matrixlib::Matrix nonSquare(1,3);
			REQUIRE_THROWS_AS(nonSquare.doUpperTriangular(),std::invalid_argument);
		}

		SECTION("Returns self if already upper triangular")
		{
			matrixlib::Matrix U(4,4);
			const unsigned short n = U.getRow();
			for(unsigned short i = 0; i < n; ++i)
			{
				for(unsigned short j = 0; j < n; ++j) U(i,j) = (i <= j) ? (i*3.24 + j*7.68 + 2.71*3.14) : (0.0);
			}

			matrixlib::Matrix upper = U.doUpperTriangular();
			REQUIRE(upper == U);
		}

		SECTION("Correct upper triangular form without normalization")
		{
			matrixlib::Matrix m(3,3);

			m(0,0) = 2.0; m(0,1) = 3.0; m(0,2) = 5.0;
			m(1,0) = 7.0; m(1,1) = 1.1; m(1,2) = 1.3;
			m(2,0) = 1.7; m(2,1) = 1.9; m(2,2) = 2.3;


			matrixlib::Matrix upper = m.doUpperTriangular(false);

			//std::cout<< upper;

			REQUIRE(matrixlib::isUpperTriangular(upper) == true);
		}

		SECTION("Normalization makes all pivots 1")
		{
			matrixlib::Matrix matrix(4,4);
			constexpr double epsilon = 1e-9;

			const unsigned short n = matrix.getRow();
			for(unsigned short i = 0; i < n; ++i)
			{
				for(unsigned short j = 0; j < n; ++j) matrix(i,j) = (i*2.3 + j*5.7 + 1.1*1.3*1.7*2.3);
			}

			matrixlib::Matrix upper = matrix.doUpperTriangular();
			//std::cout<< upper.doUpperTriangular(true);
			for(unsigned short i = 0; i < n; ++i)
			{
				CAPTURE(i,upper(i,i));
				if(std::abs(upper(i,i)) > epsilon) REQUIRE(std::abs(upper(i,i) - 1.0) < epsilon);
			}
		}

		SECTION("Zero rows handled correctly")
		{
			matrixlib::Matrix m(3, 3);
			m(0,0) = 1; m(0,1) = 2; m(0,2) = 3;
			m(1,0) = 0; m(1,1) = 0; m(1,2) = 0; // zero row
			m(2,0) = 4; m(2,1) = 5; m(2,2) = 6;

			matrixlib::Matrix upper = m.doUpperTriangular(false);

			for (unsigned short j = 0; j < 3; ++j)
			{
				REQUIRE(upper(2, j) == 0.0);
			}
		}
	}

	SECTION("doLowerTriangular")
	{
		SECTION("Throws if matrix is not square.")
		{
			matrixlib::Matrix nonSquare(1,3);
			REQUIRE_THROWS_AS(nonSquare.doLowerTriangular(),std::invalid_argument);
		}

		SECTION("Returns lower triangular form(no normalization)")
		{
			matrixlib::Matrix matrix(4,4);
			const unsigned short n = matrix.getRow();

			for(unsigned short i = 0; i < n; ++i)
			{
				for(unsigned short j = 0; j < n; ++j) matrix(i,j) = i*3.24 + j*2.71 - 22.0/7.0;
			}

			//std::cout<<matrix.doLowerTriangular(false);

			REQUIRE(matrixlib::isLowerTriangular(matrix.doLowerTriangular(false)) == true);
		}

		SECTION("Normalization makes all pivots 1")
		{
			matrixlib::Matrix A(3,3);
			const unsigned short n = A.getRow();

			A(0,0) = 2.357; A(0,1) = 7.113; A(0,2) = 13.17;
			A(1,0) = 19.23; A(1,1) = 29.31; A(1,2) = 37.41;
			A(2,0) = 43.47; A(2,1) = 53.56; A(2,2) = 59.61;


			matrixlib::Matrix normLower = A.doLowerTriangular();

			for(unsigned short i = 0; i < n; ++i)
			{
				REQUIRE(normLower(i,i) == Catch::Approx(1.0).epsilon(1e-6));
			}
		}

		SECTION("Zero rows handled correctly")
		{
			matrixlib::Matrix A(3, 3);
			A(0,0) = 1; A(0,1) = 2; A(0,2) = 3;
			A(1,0) = 0; A(1,1) = 0; A(1,2) = 0; 
			A(2,0) = 4; A(2,1) = 5; A(2,2) = 6;

			matrixlib::Matrix lower = A.doLowerTriangular(false);

			for (unsigned short j = 0; j < 3; ++j) REQUIRE(lower(2, j) == 0.0);
		}
	}

	SECTION("doEchelon")
	{
		SECTION("normalize = true")
		{
			matrixlib::Matrix A(3,4);

			A(0,0) = 1.0; A(0,1) = 2.0; A(0,2) = 5.00; A(0,3) = 3.0;
			A(1,0) = 3.0; A(1,1) = 5.0; A(1,2) = 10.0; A(1,3) = 9.0;
			A(2,0) = 2.0; A(2,1) = 4.0; A(2,2) = 10.0; A(2,3) = 7.0;

			matrixlib::Matrix E(3,4);

			E(0,0) = 1.0; E(0,1) = 2.0; E(0,2) = 5.0; E(0,3) = 3.0;
			E(1,0) = 0.0; E(1,1) = 1.0; E(1,2) = 5.0; E(1,3) = 0.0;
			E(2,0) = 0.0; E(2,1) = 0.0; E(2,2) = 0.0; E(2,3) = 1.0;

			REQUIRE(A.doEchelon() == E);
		}
		SECTION("normalize = false")
		{
			matrixlib::Matrix A(3,4);

			A(0,0) = 1.0; A(0,1) = 2.0; A(0,2) = 5.00; A(0,3) = 3.0;
			A(1,0) = 3.0; A(1,1) = 5.0; A(1,2) = 10.0; A(1,3) = 9.0;
			A(2,0) = 2.0; A(2,1) = 4.0; A(2,2) = 10.0; A(2,3) = 7.0;

			matrixlib::Matrix E(3,4);

			E(0,0) = 1.0; E(0,1) = 2.00; E(0,2) = 5.00; E(0,3) = 3.0;
			E(1,0) = 0.0; E(1,1) = -1.0; E(1,2) = -5.0; E(1,3) = 0.0;
			E(2,0) = 0.0; E(2,1) = 0.00; E(2,2) = 0.00; E(2,3) = 1.0;

			REQUIRE(A.doEchelon(false) == E);
		}

		SECTION("have a zero row")
		{
			matrixlib::Matrix A(3,4);

			A(0,0) = 1.0; A(0,1) = 2.0; A(0,2) = 5.00; A(0,3) = 3.0;
			A(1,0) = 0.0; A(1,1) = 0.0; A(1,2) = 0.00; A(1,3) = 0.0;
			A(2,0) = 3.0; A(2,1) = 5.0; A(2,2) = 10.0; A(2,3) = 9.0;
			

			matrixlib::Matrix E(3,4);

			E(0,0) = 1.0; E(0,1) = 2.0; E(0,2) = 5.0; E(0,3) = 3.0;
			E(1,0) = 0.0; E(1,1) = 1.0; E(1,2) = 5.0; E(1,3) = 0.0;
			E(2,0) = 0.0; E(2,1) = 0.0; E(2,2) = 0.0; E(2,3) = 0.0;

			REQUIRE(A.doEchelon() == E);	
		}

		SECTION("if a pivot is zero")
		{
			matrixlib::Matrix A(3, 3);
			A(0,0) = 0.0; A(0,1) = 1.0; A(0,2) = 2.0;
			A(1,0) = 1.0; A(1,1) = 0.0; A(1,2) = 3.0;
			A(2,0) = 0.0; A(2,1) = 0.0; A(2,2) = 1.0;

			matrixlib::Matrix E(3, 3);
			E(0,0) = 1.0; E(0,1) = 0.0; E(0,2) = 3.0;
			E(1,0) = 0.0; E(1,1) = 1.0; E(1,2) = 2.0;
			E(2,0) = 0.0; E(2,1) = 0.0; E(2,2) = 1.0;

			REQUIRE(A.doEchelon() == E);
		}
	}

	SECTION("toReducedEchelon")
	{
		SECTION("first")
		{
			matrixlib::Matrix m(3,4);

			m(0,0) = 1.0; m(0,1) = 2.0; m(0,2) = -1.0; m(0,3) = -4.0;
			m(1,0) = 2.0; m(1,1) = 3.0; m(1,2) = -1.0; m(1,3) = -11.0;
			m(2,0) = -2.0; m(2,1) = 0.0; m(2,2) = -3.0; m(2,3) = 22.0;

			matrixlib::Matrix rref = m.toReducedEchelon();

			matrixlib::Matrix B(3,4);

			B(0,0) = 1.0; B(0,1) = 0.0; B(0,2) = 0.0; B(0,3) = -8.0;
			B(1,0) = 0.0; B(1,1) = 1.0; B(1,2) = 0.0; B(1,3) = 1.0;
			B(2,0) = 0.0; B(2,1) = 0.0; B(2,2) = 1.0; B(2,3) = -2.0;

			//std::cout<<rref;
			REQUIRE(rref == B);
		}

		SECTION("second")
		{
			matrixlib::Matrix m(4,6);

			m(0,0) = 2.35; m(0,1) = 7.11; m(0,2) = 1.35; m(0,3) = 1.71; m(0,4) = 1.91; m(0,5) = 2.35;
			m(1,0) = 2.95; m(1,1) = 3.11; m(1,2) = 3.75; m(1,3) = 4.11; m(1,4) = 4.31; m(1,5) = 4.75;
			m(2,0) = 5.15; m(2,1) = 5.31; m(2,2) = 5.75; m(2,3) = 5.91; m(2,4) = 6.11; m(2,5) = 6.75;
			m(3,0) = 7.15; m(3,1) = 7.31; m(3,2) = 7.95; m(0,3) = 8.31; m(0,4) = 8.71; m(0,5) = 8.95;

			//std::cout<<m.pivotsCoordinates();
			//std::cout<<m.toReducedEchelon();
		}

		SECTION("third")
		{
			matrixlib::Matrix m(4,6);

			m(0,0) = 2.0; m(0,1) = 7.0; m(0,2) = 1.0; m(0,3) = 1.0; m(0,4) = 1.0; m(0,5) = 2.0;
			m(1,0) = 2.0; m(1,1) = 3.0; m(1,2) = 3.0; m(1,3) = 4.0; m(1,4) = 4.0; m(1,5) = 4.0;
			m(2,0) = 5.0; m(2,1) = 5.0; m(2,2) = 5.0; m(2,3) = 5.0; m(2,4) = 6.0; m(2,5) = 6.0;
			m(3,0) = 7.0; m(3,1) = 7.0; m(3,2) = 7.0; m(0,3) = 8.0; m(0,4) = 8.0; m(0,5) = 8.0;

			//std::cout<<m.pivotsCoordinates();
			//std::cout<<m.toReducedEchelon();
		}
	}
}