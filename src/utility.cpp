#include"matrixlib/utility.hpp"
#include<stdexcept>
#include<cmath>

namespace matrixlib
{
	bool isZeroRow(const Matrix& matrix,unsigned short rowIndex)
	{
		constexpr double epsilon = 1e-13;

		for(unsigned short j = 0; j < matrix.getColumn(); ++j)
		{
			if(std::abs(matrix(rowIndex,j)) > epsilon) return false; 
		}
		return true;
	}

	bool isSymmetric(const Matrix& matrix)
	{
		if(!isSquare(matrix)) return false; 

		constexpr double epsilon = 1e-13;
		const unsigned short n = matrix.getRow();

		for(unsigned short i = 0; i < n; ++i)
		{
			for(unsigned short j = 0; j < n; ++j)
			{
				if(std::abs(matrix(i,j) - matrix(j,i)) > epsilon) return false;
			}
		}
		return true;
	}

	bool isIdentity(const Matrix& matrix)
	{
		if(!isSquare(matrix)) return false; 

		constexpr double epsilon = 1e-13;
		const unsigned short n = matrix.getRow();

		for(unsigned short i = 0; i < n; ++i)
		{
			for(unsigned short j = 0; j < n; ++j)
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

	bool isDiagonal(const Matrix& matrix)
	{
		if(!isSquare(matrix)) return false;

		constexpr double epsilon = 1e-13;
		const unsigned short n = matrix.getRow();

		for(unsigned short i = 0; i < n; ++i)
		{
			for(unsigned short j = 0; j < n; ++j)
			{
				if(i != j && std::abs(matrix(i,j)) > epsilon) return false;
			}
		}
		return true;
	}

	bool isUpperTriangular(const Matrix& matrix)
	{
		if(!isSquare(matrix)) return false;

		constexpr double epsilon = 1e-13;
		const unsigned short n = matrix.getRow();

		for(unsigned short i = 0; i < n; ++i)
		{
			for(unsigned short j = i + 1; j < n; ++j)
			{
				if(std::abs(matrix(j,i)) > epsilon) return false;
			}
		}
		return true;
	}

	bool isLowerTriangular(const Matrix& matrix)
	{
		if(!isSquare(matrix)) return false;

		constexpr double epsilon = 1e-13;
		const unsigned short n = matrix.getRow();

		for(unsigned short i = 0; i < n; ++i)
		{
			for(unsigned short j = i + 1; j < n; ++j)
			{
				if(std::abs(matrix(i,j)) > epsilon) return false;
			}
		}
		return true;
	}

	bool areScalarRows(const Matrix& matrix,unsigned short sourceRow,unsigned short targetRow)
	{
		if(sourceRow >= matrix.getRow() || targetRow >= matrix.getRow()) throw std::out_of_range("Row index out of range.");

		constexpr double epsilon = 1e-13;
		double ratio{0.0};
		bool ratioSet = true;

		for(unsigned short i = 0; i < matrix.getColumn(); ++i)
		{
			double a = matrix(sourceRow,i);
			double b = matrix(targetRow,i);

			if(std::abs(a) < epsilon && std::abs(b) < epsilon) continue;
			else if(std::abs(a) < epsilon || std::abs(b) < epsilon) return false;
			
			double currentRatio = b/a;

			if(ratioSet)
			{
				ratio = currentRatio;
				ratioSet = false;
			}
			else if(std::abs(currentRatio - ratio) > epsilon) return false;
			
		}
		return true;
	}


	double calcTrace(const Matrix& matrix)
	{
		if(!isSquare(matrix)) throw std::invalid_argument("Matrix must be square to calculate the trace.");

		double result{0.0};
		const unsigned short n = matrix.getRow();

		for(unsigned short i = 0; i < n; ++i) result += matrix(i,i); 
		
		return result;
	}
}