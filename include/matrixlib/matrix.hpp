#ifndef MATRIX_HPP
#define MATRIX_HPP

#include<iostream>

namespace matrixlib
{

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
            Matrix(unsigned short row_,unsigned short col_);
            Matrix(const Matrix& other);                        // Copy constructor
            Matrix& operator=(const Matrix& other);             // Copy assignment
            ~Matrix();

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
            
            Matrix echelon() const;
            Matrix doUpperTriangular() const;
            Matrix doLowerTriangular() const;
            Matrix inverse() const;             // with Gauss-Jordan Elimination

            // Input/Output Function
            void read(std::istream& is = std::cin);
            void print(std::ostream& os = std::cout) const;
    };

}

#endif