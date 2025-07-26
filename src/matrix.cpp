#include"matrixlib/matrix.hpp"
#include<iostream>
#include<stdexcept>    // for std::invalid_argument

namespace matrixlib
{
	void Matrix::allocateData()
	{
		if(row == 0 || column == 0){throw std::invalid_argument("Matrix size cannot be zero.");}	

		data = new double*[row];
		for(unsigned short i = 0; i < row; ++i)
		{
			data[i] = new double[column]();	// zero initialize
		}

		#ifdef DEBUG
			std::cout<<"\n\n\t"<<row<<"x"<<column<<" matrix allocated.\n";
		#endif
	}

	void Matrix::deallocateData()
	{
		if(data != nullptr)
		{
			for(unsigned short i = 0; i < row; ++i)
			{
				delete[] data[i]; // Delete the column array for each row
			}
			delete[] data; 		  // Delete the array of row pointers
			data = nullptr;
		}
		row = 0;
		column = 0;
	}

	void Matrix::copyDataFrom(const Matrix& other)
	{
		for(unsigned short i = 0; i < row; ++i)
		{
			for(unsigned short j = 0; j < column; ++j)
			{
				data[i][j] = other.data[i][j];
			}
		}
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
		if(this == &other){ return *this; } // Self-assignment check
		deallocateData();
		row = other.row;
		column = other.column;
		allocateData();
		copyDataFrom(other);
	}

}