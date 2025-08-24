#ifndef MATRIX_HPP
#define MATRIX_HPP

#include<iostream>

namespace matrixlib
{
    inline constexpr double EPSILON = 1e-13;

    class Matrix
    {
        private:
            double** data;
            unsigned short row = 0, column = 0;

            void allocateData();
            void deallocateData();
            void copyDataFrom(const Matrix& other);

        public:
            // Constructor - Destructor
            Matrix(unsigned short row_,unsigned short col_,double initial = 0.0);
            Matrix(const Matrix& other);                                            // Copy constructor
            Matrix& operator=(const Matrix& other);                                 // Copy assignment
            ~Matrix();

            //Special matrices
            static Matrix identity(unsigned short n);
            static Matrix zeros(unsigned short i,unsigned short j);
            static Matrix ones(unsigned short i,unsigned short j);
            static Matrix diagonal(unsigned short n);
            static Matrix scalar(unsigned short n,double value);
            static Matrix exchange(unsigned short n);

            //Pivot tools
            Matrix pivotsCoordinates() const;
            void normalizePivots();
            void organizeRows();
            void sortZeroRows();                                                                         // Organize the zero rows
            Matrix applyPivotElimination(bool normalize = true,bool lower = false,bool upper = false) const;

            // Accessors
            inline unsigned short getRow() const { return row; }
            inline unsigned short getColumn() const { return column; }

            // Element access
            double& operator()(unsigned short i, unsigned short j);             // Write access
            double operator()(unsigned short i, unsigned short j) const;        // Read-only access


            // Elementary Row Operations
            void swapRows(unsigned short i,unsigned short j);                                 // Replace row i with row j
            void scaleRow(unsigned short i,double c);                                         // Let line i be multiplied by scalar c
            void addScaledRow(unsigned short sourceRow,unsigned short targetRow,double c);    // Line source, multiplied by scalar c and added to row target
            

            // Simple Operations
            Matrix transpose() const;

            Matrix operator+(const Matrix& other) const;
            Matrix operator-(const Matrix& other) const;
            Matrix operator*(const Matrix& other) const;

            bool operator==(const Matrix& other) const;
            
            Matrix doEchelon(bool normalize = true) const;
            Matrix toReducedEchelon() const;
            Matrix doUpperTriangular(bool normalize = true) const;
            Matrix doLowerTriangular(bool normalize = true) const;
            Matrix doInverse() const;          

            // Input/Output Function
            void read(std::istream& is = std::cin);
            void print(std::ostream& os = std::cout) const;
    };

    std::ostream& operator<<(std::ostream& os,const Matrix& matrix);
    std::istream& operator>>(std::istream& is,Matrix& matrix);

}

#endif