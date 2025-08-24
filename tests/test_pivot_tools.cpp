#include"catch_amalgamated.hpp"
#include"matrixlib/matrix.hpp"
#include"matrixlib/utility.hpp"

TEST_CASE("tests for pivot tools","[pivotTools]")
{
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

	SECTION("sortZeroRows")
	{
		matrixlib::Matrix A(4,2);

		A(0,0) = 0.0; A(0,1) = 0.0;
		A(1,0) = 2.5; A(1,1) = 0.3;
		A(2,0) = 0.0; A(2,1) = 0.0;
		A(3,0) = 7.0; A(3,1) = 1.1;

		A.sortZeroRows();

		REQUIRE(matrixlib::isZeroRow(A,0) == false);
		REQUIRE(matrixlib::isZeroRow(A,1) == false);
		REQUIRE(matrixlib::isZeroRow(A,2) == true);
		REQUIRE(matrixlib::isZeroRow(A,3) == true);		

		REQUIRE(A(0,0) == 2.5);
		REQUIRE(A(0,1) == 0.3);
		REQUIRE(A(1,0) == 7.0);
		REQUIRE(A(1,1) == 1.1);

		for (unsigned short j = 0; j < A.getColumn(); ++j)
		{
			REQUIRE(A(2,j) == 0.0);
			REQUIRE(A(3,j) == 0.0);
		}
	}

	SECTION("pivotsCoordinates")
	{
		SECTION("Matrix has only a pivot")
		{
			matrixlib::Matrix m(3,3);

			m(0,0) = 2.0; m(0,1) = 4.0; m(0,2) = 6.0;
			m(1,0) = 0.0; m(1,1) = 0.0; m(1,2) = 0.0;
			m(2,0) = 3.0; m(2,1) = 6.0; m(2,2) = 9.0;

			matrixlib::Matrix pivotsInf = m.pivotsCoordinates();

			REQUIRE(pivotsInf.getRow() == Catch::Approx(1.0));
			REQUIRE(pivotsInf.getColumn() == Catch::Approx(2.0));

			//std::cout<<"matrix has only a pivot:\n"<<pivotsInf;
		}

		SECTION("Matrix has more than one pivot")
		{
			SECTION("first")
			{
				matrixlib::Matrix A(7,3);

				A(0,0) = 1.0; A(0,1) = 0.0; A(0,2) = 0.0;
				A(1,0) = 2.0; A(1,1) = 0.0; A(1,2) = 0.0;
				A(2,0) = 0.0; A(2,1) = 0.0; A(2,2) = 5.0;
				A(3,0) = 0.0; A(3,1) = 1.0; A(3,2) = 0.0;
				A(4,0) = 4.0; A(4,1) = 0.0; A(4,2) = 0.0;
				A(5,0) = 0.0; A(5,1) = 0.0; A(5,2) = 0.0;
				A(6,0) = 2.3; A(6,1) = 5.7; A(6,2) = 1.1;

				matrixlib::Matrix pivots = A.pivotsCoordinates();

				//std::cout<<"matrix has more than one pivot:\n"<<pivots;

				REQUIRE(pivots.getRow() == 3);
				REQUIRE(pivots.getColumn() == 2);

				// first pivot coordinates
				REQUIRE(pivots(0,0) == Catch::Approx(0.0));
				REQUIRE(pivots(0,1) == Catch::Approx(0.0));

				// second pivot
				REQUIRE(pivots(1,0) == Catch::Approx(3.0));
				REQUIRE(pivots(1,1) == Catch::Approx(1.0));

				// third pivot
				REQUIRE(pivots(2,0) == Catch::Approx(2.0));
				REQUIRE(pivots(2,1) == Catch::Approx(2.0));
			}

			SECTION("second")
			{
				matrixlib::Matrix m(3,3);

				m(0,0) = 2.0; m(0,1) = 3.0; m(0,2) = 5.0;
				m(1,0) = 7.0; m(1,1) = 1.0; m(1,2) = 1.0;
				m(2,0) = 1.0; m(2,1) = 1.0; m(2,2) = 2.0;

				matrixlib::Matrix pivots = m.pivotsCoordinates();

				//std::cout<<pivots;

				REQUIRE(pivots(0,0) == 0);
				REQUIRE(pivots(0,1) == 0);

				REQUIRE(pivots(1,0) == 1);
				REQUIRE(pivots(1,1) == 1);

				REQUIRE(pivots(2,0) == 2);
				REQUIRE(pivots(2,1) == 2);
			}

			SECTION("third")
			{
				matrixlib::Matrix m(3,4);

				m(0,0) = 1.0; m(0,1) = 2.0; m(0,2) = 3.0; m(0,3) = 4.0;
				m(1,0) = 0.0; m(1,1) = 2.0; m(1,2) = 5.0; m(1,3) = 6.0;
				m(2,0) = 0.0; m(2,1) = 0.0; m(2,2) = 1.0; m(2,3) = 2.0;

				matrixlib::Matrix pivots = m.pivotsCoordinates();
				//std::cout<<"third:\n"<<pivots;

				REQUIRE(pivots(0,0) == 0);
				REQUIRE(pivots(0,1) == 0);

				REQUIRE(pivots(1,0) == 1);
				REQUIRE(pivots(1,1) == 1);

				REQUIRE(pivots(2,0) == 2);
				REQUIRE(pivots(2,1) == 2);
			}
		}

		SECTION("Matrix has not a pivot")
		{
			matrixlib::Matrix m(3,3);

			m(0,0) = 0.0; m(0,1) = 0.0; m(0,2) = 0.0;
			m(1,0) = 0.0; m(1,1) = 0.0; m(1,2) = 0.0;
			m(2,0) = 0.0; m(2,1) = 0.0; m(2,2) = 0.0;

			matrixlib::Matrix pivotsInf = m.pivotsCoordinates();

			REQUIRE(pivotsInf.getRow() == Catch::Approx(0.0));
			REQUIRE(pivotsInf.getColumn() == Catch::Approx(2.0));

			//std::cout<<"matrix has not a pivot:\n"<<pivotsInf;
		}

		SECTION("Matrix has only a column")
		{
			matrixlib::Matrix m(3,1);

			m(0,0) = 1.0;
			m(1,0) = 2.0;
			m(2,0) = 3.0;

			matrixlib::Matrix pivotsInf = m.pivotsCoordinates();

			//std::cout<<"matrix has only a column:\n"<<pivotsInf;

			REQUIRE(pivotsInf.getRow() == Catch::Approx(1.0));
			REQUIRE(pivotsInf.getColumn() == Catch::Approx(2.0));

			REQUIRE(pivotsInf(0,0) == Catch::Approx(0.0));
			REQUIRE(pivotsInf(0,1) == Catch::Approx(0.0));
		}

		SECTION("Matrix is a null matrix")
		{
			matrixlib::Matrix m(0,0);
			matrixlib::Matrix pivotsInf = m.pivotsCoordinates();

			//std::cout<<"matrix is a null matrix:\n"<<pivotsInf;

			REQUIRE(pivotsInf.getRow() == Catch::Approx(0.0));
			REQUIRE(pivotsInf.getColumn() == Catch::Approx(2.0));
		}

		SECTION("Pivots not on diagonal in square matrix")
		{
			matrixlib::Matrix m(4,4);

			m(0,0) = 0.0; m(0,1) = 1.0; m(0,2) = 2.0; m(0,3) = 3.0;
			m(1,0) = 0.0; m(1,1) = 0.0; m(1,2) = 0.0; m(1,3) = 0.0;
			m(2,0) = 0.0; m(2,1) = 0.0; m(2,2) = 1.0; m(2,3) = 2.0;
			m(3,0) = 0.0; m(3,1) = 0.0; m(3,2) = 0.0; m(3,3) = 0.0;

			matrixlib::Matrix pivots = m.pivotsCoordinates();
			//std::cout<<"pivots not on diagonal in square matrix:\n"<<pivots;

			REQUIRE(pivots.getRow() == 2);
			REQUIRE(pivots.getColumn() == 2);

			REQUIRE(pivots(0,0) == Catch::Approx(0.0));
			REQUIRE(pivots(0,1) == Catch::Approx(1.0));

			REQUIRE(pivots(1,0) == Catch::Approx(2.0));
			REQUIRE(pivots(1,1) == Catch::Approx(2.0));
		}

		SECTION("Pivots not on diagonal in nonSquare matrix")
		{
			matrixlib::Matrix m(3,4);

			m(0,0) = 0.0; m(0,1) = 1.0; m(0,2) = 2.0; m(0,3) = 3.0;
			m(1,0) = 0.0; m(1,1) = 0.0; m(1,2) = 0.0; m(1,3) = 0.0;
			m(2,0) = 0.0; m(2,1) = 0.0; m(2,2) = 1.0; m(2,3) = 2.0;

			matrixlib::Matrix pivots = m.pivotsCoordinates();
			//std::cout<<"pivots not on diagonal in nonSquare matrix:\n"<<pivots;

			REQUIRE(pivots.getRow() == 2);
			REQUIRE(pivots.getColumn() == 2);

			REQUIRE(pivots(0,0) == Catch::Approx(0.0));
			REQUIRE(pivots(0,1) == Catch::Approx(1.0));

			REQUIRE(pivots(1,0) == Catch::Approx(2.0));
			REQUIRE(pivots(1,1) == Catch::Approx(2.0));
		}
	}

	SECTION("normalizePivots")
	{
		SECTION("if pivot < zero")
		{
			matrixlib::Matrix A(5,3);

			A(0,0) = -1.0; A(0,1) = 2.00; A(0,2) = 3.00;
			A(1,0) = 2.00; A(1,1) = -4.0; A(1,2) = -6.0;
			A(2,0) = 0.00; A(2,1) = -2.0; A(2,2) = 0.00;
			A(3,0) = 0.00; A(3,1) = 0.00; A(3,2) = 1.00;
			A(4,0) = 0.00; A(4,1) = 0.00; A(4,2) = -1.0;

			A.normalizePivots();

			matrixlib::Matrix B(5,3);

			B(0,0) = 1.00; B(0,1) = -2.0; B(0,2) = -3.0;
			B(1,0) = 2.00; B(1,1) = -4.0; B(1,2) = -6.0;
			B(2,0) = 0.00; B(2,1) = 1.00; B(2,2) = 0.00;
			B(3,0) = 0.00; B(3,1) = 0.00; B(3,2) = 1.00;
			B(4,0) = 0.00; B(4,1) = 0.00; B(4,2) = -1.0;

			for(unsigned short i = 0; i < 5; ++i)
			{
				for(unsigned short j = 0; j < 3; ++j)
				{
					REQUIRE(A(i,j) == Catch::Approx(B(i,j)));
				}
			}
		}

		SECTION("The case of a nearly zero pivot")
		{
			matrixlib::Matrix A(3, 3);

		    A(0, 0) = 1e-12;    A(0, 1) = 2.0;    A(0, 2) = 3.0;
		    A(1, 0) = 0.000;    A(1, 1) = 4.0;    A(1, 2) = 5.0;
		    A(2, 0) = 0.000;    A(2, 1) = 0.0;    A(2, 2) = 6.0;

		    A.normalizePivots();

		    matrixlib::Matrix B(3, 3);
		    B(0, 0) = 1.0;                 
		    B(0, 1) = 2.0 / 1e-12;
		    B(0, 2) = 3.0 / 1e-12;
		    B(1, 0) = 0.0; B(1, 1) = 1.0;   B(1, 2) = 5.0/4.0;  
		    B(2, 0) = 0.0; B(2, 1) = 0.0;   B(2, 2) = 1.0;  

		    for(unsigned short i = 0; i < 3; ++i)
		    {
		        for(unsigned short j = 0; j < 3; ++j)
		        {
		            REQUIRE(A(i,j) == Catch::Approx(B(i,j)).epsilon(1e-13));
		        }
		    }	
		}

		SECTION("Rows are normalize")
		{
			matrixlib::Matrix A(3, 3);

		    A(0, 0) = 1.0; A(0, 1) = 2.0; A(0, 2) = 3.0;  
		    A(1, 0) = 0.0; A(1, 1) = 1.0; A(1, 2) = 4.0;  
		    A(2, 0) = 0.0; A(2, 1) = 0.0; A(2, 2) = 1.0;  

		    matrixlib::Matrix B = A; 

		    A.normalizePivots();

		    for (unsigned short i = 0; i < 3; ++i)
		    {
		        for (unsigned short j = 0; j < 3; ++j)
		        {
		            REQUIRE(A(i, j) == Catch::Approx(B(i, j)));
		        }
		    }
		}

		SECTION("In case two rows share the same pivot column (which is not supposed to happen)")
		{
			matrixlib::Matrix A(4, 3);

		    A(0, 0) = 2.0; A(0, 1) = 1.0; A(0, 2) = 3.0;  
		    A(1, 0) = 1.0; A(1, 1) = 4.0; A(1, 2) = 6.0;  
		    A(2, 0) = 0.0; A(2, 1) = 0.0; A(2, 2) = 5.0;  
		    A(3, 0) = 0.0; A(3, 1) = 1.0; A(3, 2) = 0.0;  
		    
		    A.normalizePivots();

		    matrixlib::Matrix B(4, 3);
		    B(0, 0) = 1.0;  B(0, 1) = 0.5;  B(0, 2) = 1.5;   
		    B(1, 0) = 0.25;  B(1, 1) = 1.0;  B(1, 2) = 6.0/4.0;   
		    B(2, 0) = 0.0;  B(2, 1) = 0.0;  B(2, 2) = 1.0;   
		    B(3, 0) = 0.0;  B(3, 1) = 1.0;  B(3, 2) = 0.0;   

		    matrixlib::Matrix pivots = A.pivotsCoordinates();
		    //std::cout<<"A pivots:\n"<<pivots;

		    for (unsigned short i = 0; i < 4; ++i)
		    {
		        for (unsigned short j = 0; j < 3; ++j)
		        {
		        	CAPTURE(i,j);
		            REQUIRE(A(i, j) == Catch::Approx(B(i, j)));
		        }
		    }
		}
	}

	SECTION("organizeRows")
	{
		matrixlib::Matrix A(4,4);

		A(0,0) = 2.3; A(0,1) = 2.0; A(0,2) = 3.0; A(0,3) = 5.0;
		A(1,0) = 0.0; A(1,1) = 0.0; A(1,2) = 5.7; A(1,3) = 5.0;
		A(2,0) = 0.0; A(2,1) = 0.0; A(2,2) = 0.0; A(2,3) = 1.1;
		A(3,0) = 0.0; A(3,1) = 1.3; A(3,2) = 0.0; A(3,3) = 0.0;

		A.organizeRows();

		matrixlib::Matrix B(4,4);

		B(0,0) = 2.3; B(0,1) = 2.0; B(0,2) = 3.0; B(0,3) = 5.0;
		B(1,0) = 0.0; B(1,1) = 1.3; B(1,2) = 0.0; B(1,3) = 0.0;
		B(2,0) = 0.0; B(2,1) = 0.0; B(2,2) = 5.7; B(2,3) = 5.0;
		B(3,0) = 0.0; B(3,1) = 0.0; B(3,2) = 0.0; B(3,3) = 1.1;

		REQUIRE(A == B);
	}

	SECTION("applyPivotElimination")
	{
		SECTION("Lower only")
		{
			matrixlib::Matrix A(5, 4);

		    A(0,0) = 2.0; A(0,1) = 4.0; A(0,2) = 6.0; A(0,3) = 8.0;
		    A(1,0) = 0.0; A(1,1) = 1.0; A(1,2) = 2.0; A(1,3) = 3.0;
		    A(2,0) = 0.0; A(2,1) = 0.0; A(2,2) = 0.0; A(2,3) = 0.0;
		    A(3,0) = 1.0; A(3,1) = 2.0; A(3,2) = 3.0; A(3,3) = 4.0;
		    A(4,0) = 0.0; A(4,1) = 0.0; A(4,2) = 1.0; A(4,3) = 2.0;

		    matrixlib::Matrix B(5, 4);

		    B(0,0) = 2.0; B(0,1) = 4.0; B(0,2) = 6.0; B(0,3) = 8.0;
		    B(1,0) = 0.0; B(1,1) = 1.0; B(1,2) = 2.0; B(1,3) = 3.0;
		    B(2,0) = 0.0; B(2,1) = 0.0; B(2,2) = 1.0; B(2,3) = 2.0;
		    B(3,0) = 0.0; B(3,1) = 0.0; B(3,2) = 0.0; B(3,3) = 0.0;
		    B(4,0) = 0.0; B(4,1) = 0.0; B(4,2) = 0.0; B(4,3) = 0.0;

		    matrixlib::Matrix result = A.applyPivotElimination(false,true);

		    REQUIRE(result == B);
		}

		SECTION("Lower + normalize")
		{
			matrixlib::Matrix A(5, 4);

		    A(0,0) = 2.0; A(0,1) = 4.0; A(0,2) = 6.0; A(0,3) = 8.0;
		    A(1,0) = 0.0; A(1,1) = 1.0; A(1,2) = 2.0; A(1,3) = 3.0;
		    A(2,0) = 0.0; A(2,1) = 0.0; A(2,2) = 0.0; A(2,3) = 0.0;
		    A(3,0) = 1.0; A(3,1) = 2.0; A(3,2) = 3.0; A(3,3) = 4.0;
		    A(4,0) = 0.0; A(4,1) = 0.0; A(4,2) = 1.0; A(4,3) = 2.0;

		    matrixlib::Matrix B(5, 4);

		    B(0,0) = 1.0; B(0,1) = 2.0; B(0,2) = 3.0; B(0,3) = 4.0;
		    B(1,0) = 0.0; B(1,1) = 1.0; B(1,2) = 2.0; B(1,3) = 3.0;
		    B(2,0) = 0.0; B(2,1) = 0.0; B(2,2) = 1.0; B(2,3) = 2.0;
		    B(3,0) = 0.0; B(3,1) = 0.0; B(3,2) = 0.0; B(3,3) = 0.0;
		    B(4,0) = 0.0; B(4,1) = 0.0; B(4,2) = 0.0; B(4,3) = 0.0;

		    matrixlib::Matrix result = A.applyPivotElimination(true,true);

		    REQUIRE(result == B);
		}

		SECTION("Upper only")
		{
			matrixlib::Matrix A(5, 4);

		    A(0,0) = 2.0; A(0,1) = 4.0; A(0,2) = 6.0; A(0,3) = 8.0;
		    A(1,0) = 0.0; A(1,1) = 1.0; A(1,2) = 2.0; A(1,3) = 3.0;
		    A(2,0) = 0.0; A(2,1) = 0.0; A(2,2) = 0.0; A(2,3) = 0.0;
		    A(3,0) = 1.0; A(3,1) = 2.0; A(3,2) = 3.0; A(3,3) = 4.0;
		    A(4,0) = 0.0; A(4,1) = 0.0; A(4,2) = 1.0; A(4,3) = 2.0;

		    matrixlib::Matrix B(5, 4);

		    B(0,0) = 2.0; B(0,1) = 0.0; B(0,2) = 0.0; B(0,3) = 0.0;
		    B(1,0) = 0.0; B(1,1) = 1.0; B(1,2) = 0.0; B(1,3) = -1.0;
		    B(2,0) = 0.0; B(2,1) = 0.0; B(2,2) = 1.0; B(2,3) = 2.0;
		    B(3,0) = 0.0; B(3,1) = 0.0; B(3,2) = 0.0; B(3,3) = 0.0;
		    B(4,0) = 0.0; B(4,1) = 0.0; B(4,2) = 0.0; B(4,3) = 0.0;

		    matrixlib::Matrix result = A.applyPivotElimination(false,false,true);

		    REQUIRE(result == B);
		}

		SECTION("Upper + normalize")
		{
			matrixlib::Matrix A(5, 4);

		    A(0,0) = 2.0; A(0,1) = 4.0; A(0,2) = 6.0; A(0,3) = 8.0;
		    A(1,0) = 0.0; A(1,1) = 1.0; A(1,2) = 2.0; A(1,3) = 3.0;
		    A(2,0) = 0.0; A(2,1) = 0.0; A(2,2) = 0.0; A(2,3) = 0.0;
		    A(3,0) = 1.0; A(3,1) = 2.0; A(3,2) = 3.0; A(3,3) = 4.0;
		    A(4,0) = 0.0; A(4,1) = 0.0; A(4,2) = 1.0; A(4,3) = 2.0;

		    matrixlib::Matrix B(5, 4);

		    B(0,0) = 1.0; B(0,1) = 0.0; B(0,2) = 0.0; B(0,3) = 0.0;
		    B(1,0) = 0.0; B(1,1) = 1.0; B(1,2) = 0.0; B(1,3) = -1.0;
		    B(2,0) = 0.0; B(2,1) = 0.0; B(2,2) = 1.0; B(2,3) = 2.0;
		    B(3,0) = 0.0; B(3,1) = 0.0; B(3,2) = 0.0; B(3,3) = 0.0;
		    B(4,0) = 0.0; B(4,1) = 0.0; B(4,2) = 0.0; B(4,3) = 0.0;

		    matrixlib::Matrix result = A.applyPivotElimination(true,false,true);

		    REQUIRE(result == B);
		}

		SECTION("Lower + Upper")
		{
			matrixlib::Matrix A(5, 4);

		    A(0,0) = 2.0; A(0,1) = 4.0; A(0,2) = 6.0; A(0,3) = 8.0;
		    A(1,0) = 0.0; A(1,1) = 1.0; A(1,2) = 2.0; A(1,3) = 3.0;
		    A(2,0) = 0.0; A(2,1) = 0.0; A(2,2) = 0.0; A(2,3) = 0.0;
		    A(3,0) = 1.0; A(3,1) = 2.0; A(3,2) = 3.0; A(3,3) = 4.0;
		    A(4,0) = 0.0; A(4,1) = 0.0; A(4,2) = 1.0; A(4,3) = 2.0;

		    matrixlib::Matrix B(5, 4);

		    B(0,0) = 2.0; B(0,1) = 0.0; B(0,2) = 0.0; B(0,3) = 0.0;
		    B(1,0) = 0.0; B(1,1) = 1.0; B(1,2) = 0.0; B(1,3) = -1.0;
		    B(2,0) = 0.0; B(2,1) = 0.0; B(2,2) = 1.0; B(2,3) = 2.0;
		    B(3,0) = 0.0; B(3,1) = 0.0; B(3,2) = 0.0; B(3,3) = 0.0;
		    B(4,0) = 0.0; B(4,1) = 0.0; B(4,2) = 0.0; B(4,3) = 0.0;

		    matrixlib::Matrix result = A.applyPivotElimination(false,true,true);

		    REQUIRE(result == B);
		}

		SECTION("Lower + Upper + normalize")
		{
			matrixlib::Matrix A(5, 4);

		    A(0,0) = 2.0; A(0,1) = 4.0; A(0,2) = 6.0; A(0,3) = 8.0;
		    A(1,0) = 0.0; A(1,1) = 1.0; A(1,2) = 2.0; A(1,3) = 3.0;
		    A(2,0) = 0.0; A(2,1) = 0.0; A(2,2) = 0.0; A(2,3) = 0.0;
		    A(3,0) = 1.0; A(3,1) = 2.0; A(3,2) = 3.0; A(3,3) = 4.0;
		    A(4,0) = 0.0; A(4,1) = 0.0; A(4,2) = 1.0; A(4,3) = 2.0;

		    matrixlib::Matrix B(5, 4);

		    B(0,0) = 1.0; B(0,1) = 0.0; B(0,2) = 0.0; B(0,3) = 0.0;
		    B(1,0) = 0.0; B(1,1) = 1.0; B(1,2) = 0.0; B(1,3) = -1.0;
		    B(2,0) = 0.0; B(2,1) = 0.0; B(2,2) = 1.0; B(2,3) = 2.0;
		    B(3,0) = 0.0; B(3,1) = 0.0; B(3,2) = 0.0; B(3,3) = 0.0;
		    B(4,0) = 0.0; B(4,1) = 0.0; B(4,2) = 0.0; B(4,3) = 0.0;

		    matrixlib::Matrix result = A.applyPivotElimination(true,true,true);

		    REQUIRE(result == B);
		}

		SECTION("A trial")
		{
			matrixlib::Matrix m(3,3);

			m(0,0) = 2.0; m(0,1) = 3.0; m(0,2) = 5.0;
			m(1,0) = 7.0; m(1,1) = 1.1; m(1,2) = 1.3;
			m(2,0) = 1.7; m(2,1) = 1.9; m(2,2) = 2.3;

			//std::cout<<m.applyPivotElimination(false,true);
			REQUIRE(matrixlib::isUpperTriangular(m.applyPivotElimination(false,true)) == true);
		}
	}

}