#include"matrixlib/matrix.hpp"
#include"matrixlib/utility.hpp"
#include<iostream>
#include<stdexcept>    // for std::invalid_argument
#include<cmath>

namespace matrixlib
{
	void Matrix::allocateData()
	{
		if(row == 0 || column == 0) throw std::invalid_argument("Matrix size cannot be zero."); 

		data = new double*[row];
		for(unsigned short i = 0; i < row; ++i) data[i] = new double[column]();	// zero initialize

		#ifdef DEBUG
			std::cout<<"\n\n\t"<<row<<"x"<<column<<" matrix allocated.\n";
		#endif
	}

	void Matrix::deallocateData()
	{
		if(data != nullptr)
		{
			for(unsigned short i = 0; i < row; ++i) delete[] data[i];   // Delete the column array for each row
			delete[] data; 		  										// Delete the array of row pointers
			data = nullptr;
		}
		row = 0;
		column = 0;
	}

	void Matrix::copyDataFrom(const Matrix& other)
	{
		for(unsigned short i = 0; i < row; ++i)
		{
			for(unsigned short j = 0; j < column; ++j) data[i][j] = other.data[i][j];
		}
		#ifdef DEBUG
		std::cout << "Data copied from matrix (" << other.row << "x" << other.column << ").\n";
     	#endif
	}


	Matrix::Matrix(unsigned short row_,unsigned short column_) : row(row_),column(column_){ allocateData(); }
	
	Matrix::~Matrix()
	{
		deallocateData();
		#ifdef DEBUG
		std::cout<<"\n\n\tMatrix destroyed and memory deallocated.\n";
		#endif
	}

	Matrix::Matrix(const Matrix& other) : row(other.row),column(other.column)
	{
		allocateData();
		copyDataFrom(other);
		#ifdef DEBUG
    	std::cout << "Copy constructor called: " << row << "x" << column << " matrix copied.\n";
		#endif
	}


	Matrix& Matrix::operator=(const Matrix& other)
	{
		if(this == &other) return *this;  // Self-assignment check
		deallocateData();
		row = other.row;
		column = other.column;
		allocateData();
		copyDataFrom(other);
		return *this;
	}

	double& Matrix::operator()(unsigned short i, unsigned short j)
	{
		if(i >= row || j >= column) throw std::out_of_range("Matrix access out of bounds."); 
		return data[i][j];
	}

	double operator()(unsigned short i, unsigned short j) const
	{
		if(i >= row || j >= column) throw std::out_of_range("Matrix access out of bounds."); 
		return data[i][j];
	}


	void Matrix::swapRows(unsigned short i,unsigned short j)
	{
		if(i >= row || j >= row) throw std::out_of_range("Row index out of range in swapRows."); 
		else if(i == j) return; 

		double* temp = data[i];
		data[i] = data[j];
		data[j] = temp;
	}

	void Matrix::scaleRow(unsigned short i,double c)
	{
		if(i >= row) throw std::out_of_range("Row index out of range in scaleRow."); 
		else if(c == 0.0) throw std::invalid_argument("Scaling factor can't be zero."); 
		for(unsigned short indexCol = 0; indexCol < column; ++indexCol) data[i][indexCol] *= c;
	}

	void Matrix::addScaledRow(unsigned short sourceRow,unsigned short targetRow,double c)
	{
		if(sourceRow >= row || targetRow >= row) throw std::out_of_range("Row index out of range in addScaledRow."); 
		else if(targetRow == sourceRow) throw std::invalid_argument("Cannot add a scaled row to itself."); 
		else if(c == 0.0) throw std::invalid_argument("Scaling factor cannot be zero."); 
		for(unsigned short indexCol = 0; indexCol < column; ++indexCol) data[targetRow][indexCol] += c * data[sourceRow][indexCol];
	}


	Matrix Matrix::transpose() const
	{
		Matrix result(column,row);
		for(unsigned short i = 0; i < row; ++i)
		{
			for(unsigned short j = 0; j < column; ++j) result.data[j][i] = this->data[i][j];
		}
		return result;
	}

	Matrix Matrix::operator+(const Matrix& other) const
	{
		if(row != other.row || column != other.column) throw std::invalid_argument("Matrix dimension must match for addition."); 
		Matrix result(row,column);
		for(unsigned short i = 0; i < row; ++i)
		{
			for(unsigned short j = 0; j < column; ++j) result.data[i][j] = data[i][j] + other.data[i][j];
		}
		return result;
	}

	Matrix Matrix::operator-(const Matrix& other) const
	{
		if(row != other.row || column != other.column){ throw std::invalid_argument("Matrix dimension must match for subtraction."); }
		Matrix result(row,column);
		for(unsigned short i = 0; i < row; ++i)
		{
			for(unsigned short j = 0; j < column; ++j)
			{
				result.data[i][j] = data[i][j] - other.data[i][j];
			}
		}
		return result;
	}

	Matrix Matrix::operator*(const Matrix& other) const
	{
		if(column != other.row){ throw std::invalid_argument("The matrix multiplication operation must match the column of the first matrix with the row of the second matrix."); }
		Matrix result(row,other.column);

		for(unsigned short i = 0; i < row; ++i)
		{
			for(unsigned short j = 0; j < other.column; ++j)
			{
				for(unsigned short k = 0; k < column; ++k)
				{
					result.data[i][j] += this->data[i][k] * other.data[k][j];
				}
			}
		}
		return result;
	}

	bool Matrix::operator==(const Matrix& other) const
	{
		if(row != other.row || column != other.column) return false; 

		constexpr double epsilon = 1e-9;

		for(unsigned short i = 0; i < row; ++i)
		{
			for(unsigned short j = 0; j < column; ++j)
			{
				if(std::abs((data[i][j]-other.data[i][j]) > epsilon)) return false; 
			}
		}
		return true;
	}


	Matrix Matrix::doUpperTriangular() const
	{
		if(!isSquare(*this)) throw std::invalid_argument("Matrix must be square to convert to upper triangular form."); 

		Matrix result(*this);

		unsigned short locRow{0},locCol{0};
		while(!isUpperTriangular(result))
		{	
			bool pivotFound = true;		
			if(result.data[locRow][locCol] == 0)
			{
				pivotFound = false;
				for(unsigned short i = locRow+1; i < row; ++i)
				{
					if(result.data[i][locCol] != 0)
					{
						result.swapRows(locRow,i);
						pivotFound = true;
						break;
					}
				}
			}

			if(!pivotFound){ ++locCol; continue;}
			
			double pivotValue = result.data[locRow][locCol];
			
			if(pivotValue != 1.0){ result.scaleRow(locRow,1.0/pivotValue); }

			for(unsigned short i = locRow+1; i < row; ++i)
			{
				pivotValue = -result.data[i][locCol];
				result.addScaledRow(locRow,i,pivotValue);
			}
			++locRow;
			++locCol;
		}
		return result;
	}

	Matrix Matrix::doLowerTriangular() const
	{
		if(!isSquare(*this)) throw std::invalid_argument("Matrix must be square to convert to lower triangular form."); 

		Matrix result(*this);

		unsigned short locRow{row-1},locCol{column-1};
		double factor{0.0};
		while(!isLowerTriangular(result))
		{
			for(unsigned short i = 0; i < row; ++i)
			{
				if(result.data[i][i] == 0)
				{
					if(isZeroColumn(result,i)) continue;
					for(unsigned short i_ = 0; i_ < row; ++i_) // for find a pivot
					{
						if(result.data[i_][i] != 0) 
						{
							result.swapRows(i,i_);
							break;
						}
					}
				}
				else
				{
					factor = 1.0/result.data[i][i];
					result.scaleRow(i,factor);					// pivot will 1
				}
			}
			for(int i = locRow; i >= 0; --i)
			{
				for (int j = locCol; j >= 0; --j)
				{
					factor = -result.data[j][i];
					result.addScaledRow(i,j,factor);
				}
			}
		}
		return result;
	}

}