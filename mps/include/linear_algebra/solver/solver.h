
#ifndef __solver__H__
#define __solver__H__

#include <numeric>
#include "vector_ex.h"
#include "matrix.h"

using namespace std;

namespace linear_algebra
{
    template <typename T>
    vector<T> forward_substitution(const matrix<T> &A, const vector<T> &b);

    template <typename T>
    vector<T> back_substitution(const matrix<T> &A, const vector<T> &b);

    template <typename T>
    class Solver
    {
        matrix<T> A;
        vector<T> b;
        vector<T> x;
        size_t N;

    public:
        Solver(matrix<T> A, vector<T> b);
        ~Solver();

        vector<T> &forward_substitution();
        vector<T> &back_substitution();
        vector<T> &gauss_elimination();
        vector<T> &conjugate_gradient(const vector<T> &initial_guess, const T &precision);
        vector<T> &incomplete_cholesky_factorization_conjugate_gradient(const vector<T> &initial_guess, const T &precision);
    };

    template <typename T>
    Solver<T>::Solver(matrix<T> A, vector<T> b)
        : A(A), b(b)
    {
        x.resize(A.column, 0);
        N = b.size();
    }

    template <typename T>
    Solver<T>::~Solver()
    {
        A.destroy();
        b.clear();
        x.clear();
        b.shrink_to_fit();
        x.shrink_to_fit();
    }

    template <typename T>
    vector<T> &Solver<T>::back_substitution()
    {
        for (size_t row = A.row - 1; row >= 0; row--)
        {
            x[row] = b[row];
            for (size_t i = row + 1; i < A.column; i++)
            {
                x[row] -= A[row][i] * x[i];
            }
            x[row] *= 1 / A[row][row];
        }
        return x;
    }

    template <typename T>
    vector<T> &Solver<T>::forward_substitution()
    {
        for (size_t row = 0; row < A.row; row++)
        {
            x[row] = b[row];
            for (size_t i = 0; i < row; i++)
            {
                x[row] -= A[row][i] * x[i];
            }
            x[row] *= 1 / A[row][row];
        }
        return x;
    }

    template <typename T>
    vector<T> &Solver<T>::conjugate_gradient(const vector<T> &initial_guess, const T &precision)
    {
        x = initial_guess;

        vector<T> p(N);
        auto t = A * x;
        p = b - t;

        vector<T> r(p);
        auto r_square = r * r;

        if (r_square < precision)
            return x;

        for (size_t i = 0; i < N; i++)
        {
            t = A * p;

            auto alpha = r_square / (t * p);

            x += alpha * p;
            r -= alpha * t;

            auto r_square_new = r * r;
            if (r_square_new < precision)
                break;

            auto beta = r_square_new / r_square;

            p = r + beta * p;

            r_square = r_square_new;
        }

        return x;
    }

    template <typename T>
    vector<T> &Solver<T>::incomplete_cholesky_factorization_conjugate_gradient(const vector<T> &initial_guess, const T &precision)
    {
        x = initial_guess;

        matrix<T> L = incomplete_cholesky_factorization(A);
        matrix<T> L_T = L.transpose();

        vector<T> r(N);
        auto k = A * x;
        r = b - (A * x);
        if (r * r < precision)
            return x;

        vector<T> z(N);
        auto w = linear_algebra::forward_substitution(L, r);
        z = linear_algebra::back_substitution(L_T, w);

        vector<T> p(z);

        for (size_t i = 0; i < N; i++)
        {
            auto t = A * p;
            auto dot_r_z = r * z;
            auto alpha = dot_r_z / (t * p);

            x += alpha * p;
            r -= alpha * t;

            if (r * r < precision)
                break;

            w = linear_algebra::forward_substitution(L, r);
            z = linear_algebra::back_substitution(L_T, w);

            auto beta = r * z / dot_r_z;
            p = z + beta * p;
        }

        return x;
    }

    template <typename T>
    vector<T> forward_substitution(const matrix<T> &A, const vector<T> &b)
    {
        vector<T> x(b);

        for (size_t row = 0; row < A.row; row++)
        {
            for (size_t i = 0; i < row; i++)
            {
                x[row] -= A[row][i] * x[i];
            }
            x[row] *= 1 / A[row][row];
        }
        return x;
    }

    template <typename T>
    vector<T> back_substitution(const matrix<T> &A, const vector<T> &b)
    {
        vector<T> x(b);
        size_t N = b.size();

        x[N - 1] = b[N - 1] / A[N - 1][N - 1];
        for (int row = N - 2; row >= 0; row--)
        {
            for (size_t i = row + 1; i < N; i++)
            {
                x[row] -= A[row][i] * x[i];
            }
            x[row] *= 1 / A[row][row];
        }
        return x;
    }
}

#endif