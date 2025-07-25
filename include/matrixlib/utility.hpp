#ifndef UTILITY_HPP
#define UTILITY_HPP

#include"matrix.hpp" // Include Matrix definition

namespace matrixlib
{
    inline bool isSquare(const Matrix& matrix) { return matrix.getRow() == matrix.getColumn(); } // Checks if the matrix is square (rows == columns)
    bool isZeroRow(const Matrix& matrix,unsigned short rowIndex);                                // Checks if the specified row is entirely zero
    bool isSymmetric(const Matrix& matrix);                                                      // Checks if the matrix is symmetric (A == Aᵀ)
    bool isIdentity(const Matrix& matrix);                                                       // Checks if the matrix is an identity matrix
    bool isUpperTriangular(const Matrix& matrix);                                                // Checks if the matrix is upper triangular
    bool isLowerTriangular(const Matrix& matrix);                                                // Checks if the matrix is lower triangular

    double calcTrace(const Matrix& matrix); // Calculates the trace of the matrix (sum of diagonal elements)
    
    unsigned short calcRank(const Matrix& matrix);  // Calculates the rank of the matrix
}


#endif