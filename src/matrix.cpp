#include"matrixlib/matrix.hpp"
#include"matrixlib/utility.hpp"
#include<iostream>
#include<stdexcept>    // for std::invalid_argument
#include<cmath>

namespace matrixlib
{
	void Matrix::allocateData()
	{
		if(row == 0 || column == 0)
		{
			data = nullptr;
			#ifdef DEBUG
			std::cout << "\n\n\tEmpty matrix (" << row << "x" << column << ") allocated.\n";
			#endif
			return;
		}

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


	Matrix::Matrix(unsigned short row_,unsigned short column_,double initial) : row(row_),column(column_)
	{ 
		allocateData(); 
			
		for(unsigned short i = 0; i < row; ++i)
		{
			for(unsigned short j = 0; j < column; ++j)
			{
				data[i][j] = initial;
			}
		}
	}
	
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


    Matrix Matrix::identity(unsigned short n)
	{
		Matrix result(n,n);
		for(unsigned short i = 0; i < n; ++i) result.data[i][i] = 1.0;
		return result;
	}

	Matrix Matrix::zeros(unsigned short i, unsigned short j)
	{
		return Matrix(i,j);
	}

	Matrix Matrix::ones(unsigned short i,unsigned short j)
	{
		return Matrix(i,j,1.0);
	}

	Matrix Matrix::diagonal(unsigned short n)
	{
		Matrix result(n,n);
		for(unsigned short i = 0; i < n; ++i) result.data[i][i] = 1.0;
		return result;
	}

	Matrix Matrix::scalar(unsigned short n,double value)
	{
		Matrix result(n,n);

		for(unsigned short i = 0; i < n; ++i) result(i,i) = value;
		return result;
	}

	Matrix Matrix::exchange(unsigned short n)
	{
		Matrix result(n,n);

		for(unsigned short i = 0; i < result.getRow(); ++i) result.data[i][n-1-i] = 1.0;
		return result;
	}


	Matrix Matrix::pivotsCoordinates() const
	{
		Matrix matrix(*this);
		unsigned short scalarRowsCount = 0;
		bool* isScalar = new bool[row];

		for(unsigned short i = 0; i < row; ++i) isScalar[i] = false;

		for(unsigned short i = 0; i < row; ++i)
		{
			if(isScalar[i]) continue;
			for(unsigned short j = i + 1; j < row; ++j)
			{
				if(areScalarRows(matrix,i,j)) isScalar[j] = true;
			}
		}

		Matrix pivotInf(row,2);
		bool* usedCol = new bool[column];
		for(unsigned short i = 0; i < column; ++i) usedCol[i] = false;

		unsigned short count = 0;
		for(unsigned short i = 0; i < row; ++i)
		{
			if(isScalar[i]) continue;
			for(unsigned short j = 0; j < column; ++j)
			{
				if(std::abs(matrix.data[i][j]) > EPSILON && !usedCol[j])
				{
					pivotInf.data[count][0] = i;
					pivotInf.data[count][1] = j;
					++count;
					usedCol[j] = true;
					break;				
				}
			}
		}

		delete[] isScalar;
		delete[] usedCol;

		Matrix result(count,2);
		for(unsigned short i = 0; i < count; ++i)
		{
			result.data[i][0] = static_cast<unsigned short>(pivotInf.data[i][0]);
    		result.data[i][1] = static_cast<unsigned short>(pivotInf.data[i][1]);
		}
		if(result.data == nullptr) return result;

		for(unsigned short i = 0; i < result.getRow(); ++i)
		{
			for(unsigned short j = 0; j + 1 < result.getRow() - i; ++j)
			{
				if(result.data[j][1] > result.data[j+1][1] || result.data[j][1] == result.data[j+1][1] && result.data[j][0] > result.data[j+1][0]) 
				{ 
					result.swapRows(j,j+1); 
				}
			}
		}

		return result;
	}

	void Matrix::normalizePivots()
	{
		Matrix pivotsInf = this->pivotsCoordinates();

		for(unsigned short i = 0; i < pivotsInf.getRow(); ++i)
		{
			unsigned short row_ = static_cast<unsigned short>(pivotsInf.data[i][0]);
			unsigned short col_ = static_cast<unsigned short>(pivotsInf.data[i][1]);
			double pivot = this->data[row_][col_];

			double factor = 1.0/pivot;

			this->scaleRow(row_,factor);
		}
	}

	void Matrix::organizeRows()
	{
		this->sortZeroRows();
		Matrix pivotsInf = this->pivotsCoordinates();
		Matrix temp(*this);

		for(unsigned short i = 0; i < pivotsInf.getRow(); ++i)
		{
			unsigned short srcRow = static_cast<unsigned short>(pivotsInf.data[i][0]);
			for(unsigned short j = 0; j < column; ++j)
			{
				this->data[i][j] = temp.data[srcRow][j];
			}
		}
	}

	Matrix Matrix::applyPivotElimination(bool normalize,bool lower,bool upper) const
	{
		Matrix matrix(*this);
		bool* isScalar = new bool[row];

		for(unsigned short i = 0; i < row; ++i) isScalar[i] = false;

		for(unsigned short i = 0; i < row; ++i)
		{
			if(isScalar[i]) continue;
			for(unsigned short j = i + 1; j < row; ++j)
			{
				if(areScalarRows(matrix,i,j)) isScalar[j] = true;
			}
		}

		for(unsigned short i = 0; i < row; ++i)
		{
			if(isScalar[i])
			{
				for(unsigned short j = 0; j < column; ++j)
				{
					matrix.data[i][j] = 0.0;
				}
			}
		}

		delete[] isScalar;

		matrix.sortZeroRows();
		matrix.organizeRows();

		Matrix pivotsInf = matrix.pivotsCoordinates();

		if(lower)
		{
			for(unsigned short i = 0; i < pivotsInf.getRow(); ++i)
			{
				unsigned short pivotCol = static_cast<unsigned short>(pivotsInf.data[i][1]);

				if(i == row-1) continue;
				for(unsigned short k = i + 1; k < row; ++k)
				{
					if(std::abs(matrix.data[k][pivotCol]) > EPSILON)
					{
						double factor = -matrix.data[k][pivotCol]/matrix.data[i][pivotCol];

						matrix.addScaledRow(i,k,factor);
					}
				}
			}
		}

		if(upper)
		{
			for(int i = pivotsInf.getRow() - 1; i >= 0; --i)
			{
				unsigned short pivotCol = static_cast<unsigned short>(pivotsInf.data[i][1]);
				
				if(i == 0) continue;
				for(int k = i - 1; k >= 0; --k)
				{
					if(std::abs(matrix.data[k][pivotCol]) > EPSILON)
					{
						double factor = -matrix.data[k][pivotCol]/matrix.data[i][pivotCol];

						matrix.addScaledRow(i,k,factor);
					}
				}
			}
		}

		if(normalize)
		{
			matrix.normalizePivots();
		}

		return matrix;
	}

	void Matrix::sortZeroRows()
	{
		unsigned short writeIndex = 0;

	    for (unsigned short i = 0; i < row; ++i)
	    {
	        if (!isZeroRow(*this, i))  
	        {
	            if (i != writeIndex)
	            {
	                swapRows(i, writeIndex);  
	            }
	            ++writeIndex;
	        }
	    }
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

	double Matrix::operator()(unsigned short i, unsigned short j) const
	{
		if(i >= row || j >= column) throw std::out_of_range("Matrix access out of bounds."); 
		return data[i][j];
	}


	void Matrix::swapRows(unsigned short i,unsigned short j)
	{
		if(i >= row || j >= row) throw std::out_of_range("Row index out of range in swapRows."); 
		if(i == j) return; 

		double* temp = data[i];
		data[i] = data[j];
		data[j] = temp;
	}

	void Matrix::scaleRow(unsigned short i,double c)
	{
		if(i >= row) throw std::out_of_range("Row index out of range in scaleRow."); 
		else if(std::abs(c) < EPSILON) return; 
		for(unsigned short indexCol = 0; indexCol < column; ++indexCol) data[i][indexCol] *= c;
	}

	void Matrix::addScaledRow(unsigned short sourceRow,unsigned short targetRow,double c)
	{
		if(sourceRow >= row || targetRow >= row) throw std::out_of_range("Row index out of range in addScaledRow."); 
		else if(targetRow == sourceRow) throw std::invalid_argument("Cannot add a scaled row to itself."); 
		else if(std::abs(c) < EPSILON) return; 
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

		for(unsigned short i = 0; i < row; ++i)
		{
			for(unsigned short j = 0; j < column; ++j)
			{
				if(std::abs(data[i][j]-other.data[i][j]) > EPSILON) return false; 
			}
		}
		return true;
	}


	Matrix Matrix::doUpperTriangular(bool normalize) const
	{
		if(!isSquare(*this)) throw std::invalid_argument("Matrix must be square to convert to upper triangular form.");
		else if(isUpperTriangular(*this)) return *this;

		Matrix result(*this);
		
		return result.applyPivotElimination(normalize,true);
	}

	Matrix Matrix::doLowerTriangular(bool normalize) const
	{
		if(!isSquare(*this)) throw std::invalid_argument("Matrix must be square to convert to lower triangular form.");
		else if(isLowerTriangular(*this)) return *this;

		Matrix result(*this);
		
		return result.applyPivotElimination(normalize,false,true);
	}

	Matrix Matrix::doEchelon(bool normalize) const
	{
		Matrix result(*this);

		return result.applyPivotElimination(normalize,true);
	}

	Matrix Matrix::toReducedEchelon() const
	{
		Matrix result(*this);
	
		return result.applyPivotElimination(true,true,true);
	}

	
	void Matrix::print(std::ostream& os) const
	{
		for(unsigned short i = 0; i < row; ++i)
		{
			os << "[ ";
			for(unsigned short j = 0; j < column; ++j)
			{
				os << data[i][j] << " ";
			}
			os << " ]\n";
		}
	}

	void Matrix::read(std::istream& is) 
	{
		for(unsigned short i = 0; i < row; ++i)
		{
			for(unsigned short j = 0; j < column; ++j)
			{
				is >> data[i][j];
			}
		}
	}


	std::ostream& operator<<(std::ostream& os,const Matrix& matrix)
	{
		matrix.print(os);
		return os;
	}

	std::istream& operator>>(std::istream& is,Matrix& matrix)
	{
		matrix.read(is);
		return is;
	}
}