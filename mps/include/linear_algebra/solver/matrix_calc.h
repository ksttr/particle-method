#ifndef __MATRIX_CALC_H__
#define __MATRIX_CALC_H__

#include <cmath>
#include "matrix_def.h"

using namespace std;

namespace linear_algebra
{
    template <typename T>
    matrix<T> cholesky_factorization(const matrix<T> &A);

    template <typename T>
    matrix<T> incomplete_cholesky_factorization(const matrix<T> &A, double precision = 1e-8);

    template <typename T>
    matrix<T> LU(const matrix<T> &A);

    template <typename T>
    matrix<T> cholesky_factorization(const matrix<T> &A)
    {
        matrix<T> L(A.row, A.column, 0);

        L[0][0] = A[0][0];
        for (size_t i = 1; i < L.row; i++)
        {
            for (size_t j = 0; j < i; j++)
            {
                L[i][j] = A[i][j];
                for (size_t k = 0; k < j; k++)
                {
                    L[i][j] -= L[i][k] * L[j][k];
                }
                L[i][j] *= 1 / L[j][j];
            }

            L[i][i] = A[i][i];
            for (size_t k = 0; k < i; k++)
            {
                L[i][i] -= L[i][k] * L[i][k];
            }
            L[i][i] = sqrt(L[i][i]);
        }

        return L;
    }

    template <typename T>
    matrix<T> incomplete_cholesky_factorization(const matrix<T> &A, double precision)
    {
        matrix<T> L(A.row, A.column, 0);

        L[0][0] = A[0][0];
        for (size_t i = 1; i < L.row; i++)
        {
            for (size_t j = 0; j < i; j++)
            {
                if (fabs(A[i][j]) < precision)
                    continue;

                L[i][j] = A[i][j];
                for (size_t k = 0; k < j; k++)
                {
                    L[i][j] -= L[i][k] * L[j][k];
                }
                L[i][j] *= 1 / L[j][j];
            }

            L[i][i] = A[i][i];
            for (size_t k = 0; k < i; k++)
            {
                L[i][i] -= L[i][k] * L[i][k];
            }
            L[i][i] = sqrt(L[i][i]);
        }

        return L;
    }
}

#endif