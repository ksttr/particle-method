#ifndef __MATRIX_DEF__H__
#define __MATRIX_DEF__H__

#include <vector>
#include <iostream>

using namespace std;

namespace linear_algebra
{
    template <typename T>
    class matrix
    {
    private:
        vector<vector<T>> A;

    public:
        size_t row;
        size_t column;

        matrix();
        matrix(size_t row, size_t column);
        matrix(size_t row, size_t column, const T &value);
        matrix(const matrix<T> &other);
        matrix(matrix<T> &&other);

        void resize(size_t r, size_t c);
        void clear();
        void shrink_to_fit();
        void destroy();

        matrix<T> &operator=(const matrix<T> &other);
        matrix<T> &operator=(matrix<T> &&other);
        vector<T> &operator[](const size_t i);
        const vector<T> &operator[](const size_t i) const;

        matrix<T> &operator+=(const matrix<T> &other);
        matrix<T> &operator+=(const T &value);
        matrix<T> operator+(const matrix<T> &other) &;
        matrix<T> operator+(const matrix<T> &other) &&;
        matrix<T> operator+(matrix<T> &&other) &;
        matrix<T> operator+(matrix<T> &&other) &&;
        matrix<T> operator+(const T &value) &;
        matrix<T> operator+(const T &value) &&;

        matrix<T> &operator-=(const matrix<T> &other);
        matrix<T> &operator-=(const T &value);
        matrix<T> operator-(const matrix<T> &other) &;
        matrix<T> operator-(const matrix<T> &other) &&;
        matrix<T> operator-(matrix<T> &&other) &;
        matrix<T> operator-(matrix<T> &&other) &&;
        matrix<T> operator-(const T &value) &;
        matrix<T> operator-(const T &value) &&;
        matrix<T> operator-() &;
        matrix<T> operator-() &&;

        matrix<T> &operator*=(const T &value);
        matrix<T> operator*(const T &value) &;
        matrix<T> operator*(const T &value) &&;
        vector<T> operator*(const vector<T> &vec) &;
        vector<T> operator*(const vector<T> &vec) &&;
        vector<T> operator*(vector<T> &&vec) &;
        vector<T> operator*(vector<T> &&vec) &&;
        vector<T> operator*(const vector<T> &vec) const &;
        vector<T> operator*(const vector<T> &&vec) const &;
        matrix<T> operator*(const matrix<T> &other) &;
        matrix<T> operator*(const matrix<T> &other) &&;
        matrix<T> operator*(matrix<T> &&other) &;
        matrix<T> operator*(matrix<T> &&other) &&;
        matrix<T> operator*(const matrix<T> &other) const &;
        matrix<T> operator*(matrix<T> &&other) const &;

        matrix<T> &operator/=(const T &value);
        matrix<T> operator/(const T &value) &;
        matrix<T> operator/(const T &value) &&;

        matrix<T> transpose();

        void print();
        void print() const;

        ~matrix();
    };

    template <typename T>
    matrix<T>::matrix() : row(0), column(0) {}

    template <typename T>
    matrix<T>::matrix(size_t row, size_t column)
        : row(row), column(column)
    {
        A.resize(row);
        for (size_t i = 0; i < row; i++)
        {
            A[i].resize(column);
        }
    }

    template <typename T>
    matrix<T>::matrix(size_t row, size_t column, const T &value)
        : matrix(row, column)
    {
        for (size_t r = 0; r < row; r++)
        {
            for (size_t c = 0; c < column; c++)
            {
                A[r][c] = value;
            }
        }
    }

    template <typename T>
    matrix<T>::matrix(const matrix<T> &other)
        : matrix(other.row, other.column)
    {
        for (size_t r = 0; r < row; r++)
        {
            for (size_t c = 0; c < column; c++)
            {
                A[r][c] = other[r][c];
            }
        }
    }

    template <typename T>
    matrix<T>::matrix(matrix<T> &&other)
        : matrix(other.row, other.column)
    {
        for (size_t r = 0; r < row; r++)
        {
            for (size_t c = 0; c < column; c++)
            {
                A[r][c] = other[r][c];
            }
        }

        other.destroy();
    }

    template <typename T>
    void matrix<T>::resize(size_t r, size_t c)
    {
        row = r;
        column = c;
        A.resize(row);
        for (size_t i = 0; i < row; i++)
        {
            A[i].resize(column);
        }
    }

    template <typename T>
    void matrix<T>::clear()
    {
        for (auto vec = A.rbegin(); vec != A.rend(); ++vec)
        {
            vec->clear();
        }
    }

    template <typename T>
    void matrix<T>::shrink_to_fit()
    {
        for (auto vec = A.rbegin(); vec != A.rend(); ++vec)
        {
            vec->shrink_to_fit();
        }
    }

    template <typename T>
    void matrix<T>::destroy()
    {
        for (auto vec = A.rbegin(); vec != A.rend(); ++vec)
        {
            vec->clear();
            vec->shrink_to_fit();
        }
    }

    template <typename T>
    matrix<T>::~matrix()
    {
        destroy();
    }

    template <typename T>
    matrix<T> &matrix<T>::operator=(const matrix<T> &other)
    {
        if (row != other.row | column != other.column)
            this->resize(other.row, other.column);

        for (size_t r = 0; r < row; r++)
        {
            for (size_t c = 0; c < column; c++)
            {
                A[r][c] = other[r][c];
            }
        }
        return *this;
    }

    template <typename T>
    matrix<T> &matrix<T>::operator=(matrix<T> &&other)
    {
        if (row != other.row | column != other.column)
            this->resize(other.row, other.column);
        for (size_t r = 0; r < row; r++)
        {
            for (size_t c = 0; c < column; c++)
            {
                A[r][c] = other[r][c];
            }
        }
        other.destroy();
        return *this;
    }

    template <typename T>
    vector<T> &matrix<T>::operator[](const size_t i)
    {
        return A[i];
    }

    template <typename T>
    const vector<T> &matrix<T>::operator[](const size_t i) const
    {
        return A[i];
    }

    template <typename T>
    matrix<T> &matrix<T>::operator+=(const matrix<T> &other)
    {
        for (size_t r = 0; r < row; r++)
        {
            for (size_t c = 0; c < column; c++)
            {
                A[r][c] += other[r][c];
            }
        }

        return *this;
    }

    template <typename T>
    matrix<T> &matrix<T>::operator+=(const T &value)
    {

        for (size_t r = 0; r < row; r++)
        {
            for (size_t c = 0; c < column; c++)
            {
                A[r][c] += value;
            }
        }

        return *this;
    }

    template <typename T>
    matrix<T> &matrix<T>::operator-=(const matrix<T> &other)
    {
        for (size_t r = 0; r < row; r++)
        {
            for (size_t c = 0; c < column; c++)
            {
                A[r][c] -= other[r][c];
            }
        }

        return *this;
    }

    template <typename T>
    matrix<T> &matrix<T>::operator-=(const T &value)
    {

        for (size_t r = 0; r < row; r++)
        {
            for (size_t c = 0; c < column; c++)
            {
                A[r][c] -= value;
            }
        }

        return *this;
    }

    template <typename T>
    matrix<T> matrix<T>::operator+(const matrix<T> &other) &
    {
        matrix<T> res(*this);
        res += other;
        return res;
    }

    template <typename T>
    matrix<T> matrix<T>::operator+(const matrix<T> &other) &&
    {
        matrix<T> res(move(*this));
        res += other;
        return res;
    }

    template <typename T>
    matrix<T> matrix<T>::operator+(matrix<T> &&other) &
    {
        matrix<T> res(move(other));
        res += *this;
        return res;
    }

    template <typename T>
    matrix<T> matrix<T>::operator+(matrix<T> &&other) &&
    {
        return *this + move(other);
    }

    template <typename T>
    matrix<T> matrix<T>::operator+(const T &value) &
    {
        matrix<T> res(*this);
        res += value;
        return res;
    }

    template <typename T>
    matrix<T> matrix<T>::operator+(const T &value) &&
    {
        matrix<T> res(move(*this));
        res += value;
        return res;
    }

    template <typename T>
    matrix<T> matrix<T>::operator-(const matrix<T> &other) &
    {
        matrix<T> res(*this);
        res -= other;
        return res;
    }

    template <typename T>
    matrix<T> matrix<T>::operator-(const matrix<T> &other) &&
    {
        matrix<T> res(move(*this));
        res -= other;
        return res;
    }

    template <typename T>
    matrix<T> matrix<T>::operator-(matrix<T> &&other) &
    {
        matrix<T> res(move(other));
        res -= *this;
        return res;
    }

    template <typename T>
    matrix<T> matrix<T>::operator-(matrix<T> &&other) &&
    {
        return *this + move(other);
    }

    template <typename T>
    matrix<T> matrix<T>::operator-(const T &value) &
    {
        matrix<T> res(*this);
        res -= value;
        return res;
    }

    template <typename T>
    matrix<T> matrix<T>::operator-(const T &value) &&
    {
        matrix<T> res(move(*this));
        res -= value;
        return res;
    }

    template <typename T>
    matrix<T> matrix<T>::operator-() &
    {
        matrix<T> res(this->row, this->column, 0);
        res -= *this;
        return res;
    }

    template <typename T>
    matrix<T> matrix<T>::operator-() &&
    {
        matrix<T> res(move(*this));
        for (size_t r = 0; r < res.row; r++)
        {
            for (size_t c = 0; c < res.column; c++)
            {
                res[r][c] = -res[r][c];
            }
        }

        return res;
    }

    template <typename T>
    matrix<T> &matrix<T>::operator*=(const T &value)
    {
        for (size_t r = 0; r < row; r++)
        {
            for (size_t c = 0; c < column; c++)
            {
                A[r][c] *= value;
            }
        }
        return *this;
    }

    template <typename T>
    matrix<T> matrix<T>::operator*(const T &value) &
    {
        matrix<T> res(*this);
        res *= value;
        return res;
    }

    template <typename T>
    matrix<T> matrix<T>::operator*(const T &value) &&
    {
        matrix<T> res(move(*this));
        res *= value;
        return res;
    }

    template <typename T>
    vector<T> matrix<T>::operator*(const vector<T> &vec) &
    {
        vector<T> res(this->row, 0);
        for (size_t r = 0; r < row; r++)
        {
            for (size_t c = 0; c < column; c++)
            {
                res.at(r) += A[r][c] * vec.at(c);
            }
        }

        return res;
    }

    template <typename T>
    vector<T> matrix<T>::operator*(const vector<T> &vec) &&
    {
        vector<T> res(this->row, 0);
        for (size_t r = 0; r < row; r++)
        {
            for (size_t c = 0; c < column; c++)
            {
                res.at(r) += A[r][c] * vec.at(c);
            }
        }

        this->destroy();

        return res;
    }

    template <typename T>
    vector<T> matrix<T>::operator*(const vector<T> &vec) const &
    {
        vector<T> res(this->row, 0);
        for (size_t r = 0; r < row; r++)
        {
            for (size_t c = 0; c < column; c++)
            {
                res.at(r) += A[r][c] * vec.at(c);
            }
        }

        return res;
    }

    template <typename T>
    vector<T> matrix<T>::operator*(const vector<T> &&vec) const &
    {
        vector<T> res(this->row, 0);
        for (size_t r = 0; r < row; r++)
        {
            for (size_t c = 0; c < column; c++)
            {
                res.at(r) += A[r][c] * vec.at(c);
            }
        }

        vec.clear();
        vec.shrink_to_fit();

        return res;
    }
    template <typename T>
    matrix<T> matrix<T>::operator*(const matrix<T> &other) &
    {
        matrix<T> res(row, other.column, 0);
        for (size_t t_r = 0; t_r < row; t_r++)
        {
            for (size_t o_c = 0; o_c < other.column; o_c++)
            {
                for (size_t k = 0; k < column; k++)
                {
                    res[t_r][o_c] += A[t_r][k] * other[k][o_c];
                }
            }
        }
        return res;
    }

    template <typename T>
    matrix<T> matrix<T>::operator*(const matrix<T> &other) &&
    {
        matrix<T> res(row, other.column, 0);
        for (size_t t_r = 0; t_r < row; t_r++)
        {
            for (size_t o_c = 0; o_c < other.column; o_c++)
            {
                for (size_t k = 0; k < column; k++)
                {
                    res[t_r][o_c] += A[t_r][k] * other[k][o_c];
                }
            }
        }
        this->destroy();
        return res;
    }

    template <typename T>
    matrix<T> matrix<T>::operator*(matrix<T> &&other) &
    {
        matrix<T> res(row, other.column, 0);
        for (size_t t_r = 0; t_r < row; t_r++)
        {
            for (size_t o_c = 0; o_c < other.column; o_c++)
            {
                for (size_t k = 0; k < column; k++)
                {
                    res[t_r][o_c] += A[t_r][k] * other[k][o_c];
                }
            }
        }
        other.destroy();
        return res;
    }

    template <typename T>
    matrix<T> matrix<T>::operator*(matrix<T> &&other) &&
    {
        return *this * move(other);
    }

    template <typename T>
    matrix<T> matrix<T>::operator*(const matrix<T> &other) const &
    {
        matrix<T> res(row, other.column, 0);
        for (size_t t_r = 0; t_r < row; t_r++)
        {
            for (size_t o_c = 0; o_c < other.column; o_c++)
            {
                for (size_t k = 0; k < column; k++)
                {
                    res[t_r][o_c] += A[t_r][k] * other[k][o_c];
                }
            }
        }
        return res;
    }

    template <typename T>
    matrix<T> matrix<T>::operator*(matrix<T> &&other) const &
    {
        matrix<T> res(row, other.column, 0);
        for (size_t t_r = 0; t_r < row; t_r++)
        {
            for (size_t o_c = 0; o_c < other.column; o_c++)
            {
                for (size_t k = 0; k < column; k++)
                {
                    res[t_r][o_c] += A[t_r][k] * other[k][o_c];
                }
            }
        }
        other.destroy();
        return res;
    }

    template <typename T>
    matrix<T> operator*(const T &value, const matrix<T> &mat)
    {
        matrix<T> res(mat);
        res *= value;
        return res;
    }

    template <typename T>
    matrix<T> operator*(const T &value, const matrix<T> &&mat)
    {
        matrix<T> res(move(mat));
        res *= value;
        return res;
    }

    template <typename T>
    matrix<T> &matrix<T>::operator/=(const T &value)
    {
        for (size_t r = 0; r < row; r++)
        {
            for (size_t c = 0; c < column; c++)
            {
                A[r][c] /= value;
            }
        }
        return *this;
    }

    template <typename T>
    matrix<T> matrix<T>::operator/(const T &value) &
    {
        matrix<T> res(*this);
        res /= value;
        return res;
    }

    template <typename T>
    matrix<T> matrix<T>::operator/(const T &value) &&
    {
        matrix<T> res(move(*this));
        res /= value;
        return res;
    }

    template <typename T>
    matrix<T> matrix<T>::transpose()
    {
        matrix<T> res(column, row);
        for (size_t r = 0; r < res.row; r++)
        {
            for (size_t c = 0; c < res.column; c++)
            {
                res[r][c] = A[c][r];
            }
        }
        return res;
    }

    template <typename T>
    void matrix<T>::print()
    {
        for (size_t i = 0; i < row; i++)
        {
            for (size_t j = 0; j < column; j++)
            {
                cout << A[i][j];
            }
            cout << endl;
        }
    }

    template <typename T>
    void matrix<T>::print() const
    {
        for (size_t i = 0; i < row; i++)
        {
            for (size_t j = 0; j < column; j++)
            {
                cout << A[i][j];
            }
            cout << endl;
        }
    }
}

#endif