#include"utility.hpp"
#include<stdexcept>
#include<cmath>

namespace matrixlib
{
	bool isZeroRow(const Matrix& matrix,unsigned short rowIndex)
	{
		for(unsigned short j = 0; j < matrix.getColumn(); ++j)
		{
			if(matrix(rowIndex,j) != 0) return false; 
		}
		return true;
	}

	bool isZeroColumn(const Matrix& matrix,unsigned short colIndex)
	{
		for(unsigned short i = 0; i < matrix.getRow(); ++i)
		{
			if(matrix(i,colIndex) != 0) return false; 
		}
		return true;
	}

	bool isSymmetric(const Matrix& matrix)
	{
		if(!isSquare(matrix)) return false; 
		Matrix Transpose(matrix);
		Transpose.transpose();
		return matrix == Transpose;
	}

	bool isIdentity(const Matrix& matrix)
	{
		if(!isSquare(matrix)){ return false; }

		constexpr double epsilon = 1e-9;

		for(unsigned short i = 0; i < matrix.getRow(); ++i)
		{
			for(unsigned short j = 0; j < matrix.getColumn(); ++j)
			{
				if(i == j)
				{
					if(std::abs(matrix(i,j) - 1.0) > epsilon) return false; 
				}
				else
				{
					if(std::abs(matrix(i,j)) > epsilon) return false; 
				}
			}
		}
		return true;
	}

	bool isUpperTriangular(const Matrix& matrix)
	{
		if(!isSquare(matrix)) return false;

		constexpr double epsilon = 1e-9;
		
		for(unsigned short i = 0; i < matrix.getRow(); ++i)
		{
			for(unsigned short j = 0; j < matrix.getColumn(); ++j)
			{
				if(i > j)
				{
					if(std::abs(matrix(i,j)) > epsilon) return false;
				}
			}
		}
		return true;
	}

	bool isLowerTriangular(const Matrix& matrix)
	{
		if(!isSquare(matrix)) return false;

		constexpr double epsilon = 1e-9;

		for(unsigned short i = 0; i < matrix.getRow(); ++i)
		{
			for(unsigned short j = 0; j < matrix.getColumn(); ++j)
			{
				if(i < j)
				{
					if(std::abs(matrix(i,j)) > epsilon) return false;
				}
			}
		}
		return true;
	}


	double calcTrace(const Matrix& matrix)
	{
		if(!isSquare(matrix)) throw std::invalid_argument("Matrix must be square to calculate the trace.");

		double result{0.0};

		for(unsigned short i = 0; i < matrix.getRow(); ++i) result += matrix(i,i); 
		
		return result;
	}
}