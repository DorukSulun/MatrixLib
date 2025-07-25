#ifndef MATH_HPP
#define MATH_HPP

#include"utility.hpp"

namespace matrixlib
{
    namespace math
    {
        double calcDeterminant(const Matrix& matrix);

        Matrix inverse() const; // with  A = LU  ->  U^-1 . L^-1 = A^-1
    }
}

#endif