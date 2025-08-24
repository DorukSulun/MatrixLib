#include"catch_amalgamated.hpp"
#include"matrixlib/matrix.hpp"
#include"matrixlib/utility.hpp"

TEST_CASE("Check are utility functions","[utility]")
{
	SECTION("isSquare")
	{
		matrixlib::Matrix squareMatrix(4,4);
		REQUIRE(matrixlib::isSquare(squareMatrix) == true);

		matrixlib::Matrix nonSquareMatrix(3,4);
		REQUIRE(matrixlib::isSquare(nonSquareMatrix) == false);

		matrixlib::Matrix singleElement(1,1);
		REQUIRE(matrixlib::isSquare(singleElement) == true);
	}

	SECTION("isZeroRow")
	{
		matrixlib::Matrix matrix(4,4);

		for(unsigned short i = 0; i < matrix.getRow(); ++i) REQUIRE(matrixlib::isZeroRow(matrix,i) == true);

		matrix(0,0) = 22.0/7.0;
		REQUIRE(matrixlib::isZeroRow(matrix,0) == false);

		for(unsigned short i = 1; i < matrix.getRow(); ++i) REQUIRE(matrixlib::isZeroRow(matrix,i) == true);
	}

	SECTION("isSymmetric")
	{
		matrixlib::Matrix nonSquareMatrix(2,3);
		REQUIRE(matrixlib::isSymmetric(nonSquareMatrix) == false);


		matrixlib::Matrix nonSymmetric(2,2);

		nonSymmetric(0,0) = 1.0; nonSymmetric(0,1) = 1.1;
		nonSymmetric(1,0) = 2.0; nonSymmetric(1,1) = 2.1;

		REQUIRE(matrixlib::isSymmetric(nonSymmetric) == false);


		matrixlib::Matrix symmetric(2,2);

		symmetric(0,0) = 1.0; symmetric(0,1) = 1.1;
		symmetric(1,0) = 1.1; symmetric(1,1) = 2.1;

		REQUIRE(matrixlib::isSymmetric(symmetric) == true);
	}

	SECTION("isIdentity")
	{
		matrixlib::Matrix identity(4,4);

		for(unsigned short i = 0; i < identity.getRow(); ++i)
		{
			for(unsigned short j = 0; j < identity.getColumn(); ++j) identity(i,j) = (i == j) ? 1.0 : 0.0;
		}
		REQUIRE(matrixlib::isIdentity(identity) == true);

		matrixlib::Matrix nonIdentity(4,4);
		for(unsigned short i = 0; i < nonIdentity.getRow(); ++i)
		{
			for(unsigned short j = 0; j < nonIdentity.getColumn(); ++j) nonIdentity(i,j) = (i == j) ? 1.0 : 0.0;
		}
		nonIdentity(0,1) = 22.0/7.0;
		REQUIRE(matrixlib::isIdentity(nonIdentity) == false);

		matrixlib::Matrix nonIdentity_(4,4);
		for(unsigned short i = 0; i < nonIdentity_.getRow(); ++i)
		{
			for(unsigned short j = 0; j < nonIdentity_.getColumn(); ++j) nonIdentity_(i,j) = (i == j) ? 1.0 : 0.0;
		}
		nonIdentity_(0,0) = 22.0/7.0;
		REQUIRE(matrixlib::isIdentity(nonIdentity_) == false);	

		matrixlib::Matrix nonSquare(1,2);
		REQUIRE(matrixlib::isIdentity(nonSquare) == false);
	}

	SECTION("isDiagonal")
	{
		matrixlib::Matrix nonSquare(1,3);
		REQUIRE(matrixlib::isDiagonal(nonSquare) == false);

		matrixlib::Matrix diagonal(4,4);
		unsigned short n = diagonal.getRow();
		for(unsigned short i = 0; i < n; ++i)
		{
			for(unsigned short j = 0; j < n; ++j) diagonal(i,j) = (i == j) ? (1.0) : (0.0);
		}
		REQUIRE(matrixlib::isDiagonal(diagonal) == true);

		matrixlib::Matrix nonDiagonal = diagonal;
		nonDiagonal(0,1) = 1.0;
		REQUIRE(matrixlib::isDiagonal(nonDiagonal) == false);

		REQUIRE(matrixlib::isSymmetric(diagonal) == true);
	}

	SECTION("isUpperTriangular")
	{
		matrixlib::Matrix nonSquare(1,3);
		REQUIRE(matrixlib::isUpperTriangular(nonSquare) == false);

		matrixlib::Matrix nonU(4,4);
		unsigned short n = nonU.getRow();
		for(unsigned short i = 0; i < n; ++i)
		{
			for(unsigned short j = 0; j < n; ++j) nonU(i,j) = i*3.14 + j*2.71 + 22.0/7.0;
		}
		REQUIRE(matrixlib::isUpperTriangular(nonU) == false);

		matrixlib::Matrix U(4,4);
		n = U.getRow();
		for(unsigned short i = 0; i < n; ++i)
		{
			for(unsigned short j = 0; j < n; ++j) U(i,j) = (i <= j) ? (1.0) : (0.0);
		}
		REQUIRE(matrixlib::isUpperTriangular(U) == true);

		if(!isDiagonal(U)) REQUIRE(matrixlib::isSymmetric(U) == false);

		REQUIRE(matrixlib::isLowerTriangular(U.transpose()) == true);

		matrixlib::Matrix U_(4,4);
		n = U_.getRow();
		for(unsigned short i = 0; i < n; ++i)
		{
			for(unsigned short j = 0; j < n; ++j) U_(i,j) = (i <= j) ? (i*3.145 + j*2.719 + 22.0/7.0) : (0.0);
		}

		if(matrixlib::isUpperTriangular(U) && matrixlib::isUpperTriangular(U_))
		{
			matrixlib::Matrix mulU = U * U_;
			REQUIRE(matrixlib::isUpperTriangular(mulU) == true);
		} 

		matrixlib::Matrix zero(3,3);
		REQUIRE(matrixlib::isLowerTriangular(zero) == true);
	}

	SECTION("isLowerTriangular")
	{
		matrixlib::Matrix nonSquare(1,3);
		REQUIRE(matrixlib::isLowerTriangular(nonSquare) == false);

		matrixlib::Matrix nonL(4,4);
		unsigned short n = nonL.getRow();
		for(unsigned short i = 0; i < n; ++i)
		{
			for(unsigned short j = 0; j < n; ++j) nonL(i,j) = i*3.14 + j*2.71 + 22.0/7.0;
		}
		REQUIRE(matrixlib::isLowerTriangular(nonL) == false);

		matrixlib::Matrix L(4,4);
		n = L.getRow();
		for(unsigned short i = 0; i < n; ++i)
		{
			for(unsigned short j = 0; j < n; ++j) L(i,j) = (i >= j) ? (1.0) : (0.0);
		}
		REQUIRE(matrixlib::isLowerTriangular(L) == true);

		if(!isDiagonal(L)) REQUIRE(matrixlib::isSymmetric(L) == false);

		REQUIRE(matrixlib::isUpperTriangular(L.transpose()) == true);

		matrixlib::Matrix L_(4,4);
		n = L_.getRow();
		for(unsigned short i = 0; i < n; ++i)
		{
			for(unsigned short j = 0; j < n; ++j) L_(i,j) = (i <= j) ? (i*3.145 + j*2.719 + 22.0/7.0) : (0.0);
		}

		if(matrixlib::isLowerTriangular(L) && matrixlib::isLowerTriangular(L_))
		{
			matrixlib::Matrix mulL = L * L_;
			REQUIRE(matrixlib::isLowerTriangular(mulL) == true);
		} 

		matrixlib::Matrix zero(3,3);
		REQUIRE(matrixlib::isLowerTriangular(zero) == true);
	}

	SECTION("areScalarRows")
	{
		SECTION("The comparison of two zero rows")
		{
			matrixlib::Matrix Z = matrixlib::Matrix::zeros(2,8);
			REQUIRE(matrixlib::areScalarRows(Z,0,1) == true);
		}

		SECTION("Rows not scalar multiple")	
		{
			matrixlib::Matrix A(2,4);

			A(0,0) = 2.3; A(0,1) = 5.7; A(0,2) = 1.1; A(0,3) = 1.3;
			A(1,0) = 1.7; A(1,1) = 1.9; A(1,2) = 2.3; A(1,3) = 2.9;

			REQUIRE(matrixlib::areScalarRows(A,0,1) == false);
		}

		SECTION("Zero row vs scaled row")
		{
			matrixlib::Matrix A(2,4);

			A(0,0) = 0.0; A(0,1) = 0.0; A(0,2) = 0.0; A(0,3) = 0.0;
			A(1,0) = 0.0; A(1,1) = 1.9; A(1,2) = 0.0; A(1,3) = 2.9;

			REQUIRE(matrixlib::areScalarRows(A,0,1) == false);
		}

		SECTION("Same rows are scalar multiples")
		{
			matrixlib::Matrix A(2,8);

			for(unsigned short i = 0; i < A.getColumn(); ++i)
			{
				A(0,i) = ((i+1)/2.3)*3.24;
				A(1,i) = A(0,i)*7.68;
			}

			REQUIRE(matrixlib::areScalarRows(A,0,1) == true);
		}

		SECTION("Floating-point tolerance")
		{
			matrixlib::Matrix M(2, 2);
		    M(0, 0) = 1.0; 		    M(0, 1) = 2.0;
		    M(1, 0) = 2.0 + 1e-14;  M(1, 1) = 4.0 + 1e-14;

        	REQUIRE(areScalarRows(M, 0, 1) == true);
		}

		SECTION("Throws on invalid row index")
		{
			matrixlib::Matrix M = matrixlib::Matrix::identity(3);
	        REQUIRE_THROWS_AS(areScalarRows(M, 0, 3), std::out_of_range);
	        REQUIRE_THROWS_AS(areScalarRows(M, 4, 1), std::out_of_range);
		}
	}

	SECTION("calcTrace")
	{
		matrixlib::Matrix nonSquare(2,3);
		REQUIRE_THROWS_AS(matrixlib::calcTrace(nonSquare), std::invalid_argument);

		matrixlib::Matrix identity(9,9);
		for(unsigned short i = 0; i < identity.getRow(); ++i)
		{
			for(unsigned short j = 0; j < identity.getColumn(); ++j) identity(i,j) = (i == j) ? 1.0 : 0.0;
		}
		REQUIRE(matrixlib::calcTrace(identity) == Catch::Approx(9.0));

		matrixlib::Matrix custom(4,4);
		double result{0.0};
		for(unsigned short i = 0; i < custom.getRow(); ++i)
		{
			for(unsigned short j = 0; j < custom.getColumn(); ++j)
			{
				custom(i,j) = (i == j) ? 2.3*i : 3.14*j;
				if(i == j) result += custom(i,j);
			}
		}
		REQUIRE(matrixlib::calcTrace(custom) == Catch::Approx(result));
	}	
}