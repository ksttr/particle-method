#ifndef __VECTOR_EX_HPP__
#define __VECTOR_EX_HPP__

#include <numeric>
#include <vector>

using namespace std;

namespace linear_algebra
{
    template <typename T>
    vector<T> operator+=(vector<T> &left, const vector<T> &right)
    {
        for (size_t i = 0; i < left.size(); i++)
        {
            left[i] += right[i];
        }
        return left;
    }

    template <typename T>
    vector<T> operator+=(vector<T> &left, const T &right)
    {
        for (size_t i = 0; i < left.size(); i++)
        {
            left[i] += right;
        }
        return left;
    }

    template <typename T>
    vector<T> operator+(const vector<T> &left, const vector<T> &right)
    {
        vector<T> res(left);
        for (size_t i = 0; i < res.size(); i++)
        {
            res[i] += right[i];
        }
        return res;
    }

    template <typename T>
    vector<T> operator+(vector<T> &&left, const vector<T> &right)
    {
        vector<T> res(move(left));
        for (size_t i = 0; i < res.size(); i++)
        {
            res[i] += right[i];
        }
        return res;
    }

    template <typename T>
    vector<T> operator+(const vector<T> &left, vector<T> &&right)
    {
        vector<T> res(move(right));
        for (size_t i = 0; i < res.size(); i++)
        {
            res[i] += left[i];
        }
        return res;
    }

    template <typename T>
    vector<T> operator+(vector<T> &&left, vector<T> &&right)
    {
        return move(left) + right;
    }

    template <typename T>
    vector<T> operator+(const vector<T> &left, const T &value)
    {
        vector<T> res(left);
        for (size_t i = 0; i < res.size(); i++)
        {
            res[i] += value;
        }
        return res;
    }

    template <typename T>
    vector<T> operator+(vector<T> &&left, const T &value)
    {
        vector<T> res(move(left));
        for (size_t i = 0; i < res.size(); i++)
        {
            res[i] += value;
        }
        return res;
    }

    template <typename T>
    vector<T> operator-=(vector<T> &left, const vector<T> &right)
    {
        for (size_t i = 0; i < left.size(); i++)
        {
            left[i] -= right[i];
        }
        return left;
    }

    template <typename T>
    vector<T> operator-=(vector<T> &left, const T &right)
    {
        for (size_t i = 0; i < left.size(); i++)
        {
            left[i] -= right;
        }
        return left;
    }

    template <typename T>
    vector<T> operator-(const vector<T> &me)
    {
        auto res = vector<T>(me);
        for (size_t i = 0; i < res.size(); i++)
        {
            res[i] = -me[i];
        }
        return res;
    }

    template <typename T>
    vector<T> operator-(vector<T> &&me)
    {
        auto res = vector<T>(move(me));
        for (size_t i = 0; i < res.size(); i++)
        {
            res[i] = -res[i];
        }
        return res;
    }

    template <typename T>
    vector<T> operator-(const vector<T> &left, const vector<T> &right)
    {
        vector<T> res(left);
        res -= right;
        return res;
    }

    template <typename T>
    vector<T> operator-(vector<T> &&left, const vector<T> &right)
    {
        vector<T> res(move(left));
        res -= right;
        return res;
    }

    template <typename T>
    vector<T> operator-(const vector<T> &left, vector<T> &&right)
    {
        vector<T> res(move(right));
        res -= left;
        return -res;
    }

    template <typename T>
    vector<T> operator-(vector<T> &&left, vector<T> &&right)
    {
        return move(left) - right;
    }

    template <typename T>
    vector<T> operator-(const vector<T> &left, const T &value)
    {
        vector<T> res(left);
        for (size_t i = 0; i < res.size(); i++)
        {
            res[i] -= value;
        }
        return res;
    }

    template <typename T>
    vector<T> operator-(const vector<T> &&left, const T &value)
    {
        vector<T> res(move(left));
        for (size_t i = 0; i < res.size(); i++)
        {
            res[i] -= value;
        }
        return res;
    }

    template <typename T>
    T operator*(const vector<T> &left, const vector<T> &right)
    {
        T res = 0;
        res = inner_product(left.begin(), left.end(), right.begin(), static_cast<T>(0));
        return res;
    }

    template <typename T>
    vector<T> operator*(const T &left, const vector<T> &right)
    {
        vector<T> res(right);
        for (size_t i = 0; i < res.size(); i++)
        {
            res[i] *= left;
        }
        return res;
    }

    template <typename T>
    vector<T> operator*(const T &left, const vector<T> &&right)
    {
        vector<T> res(move(right));
        for (size_t i = 0; i < res.size(); i++)
        {
            res[i] *= left;
        }
        return res;
    }

    template <typename T>
    vector<T> operator*(const vector<T> &left, const T &right)
    {
        vector<T> res(left);
        for (size_t i = 0; i < res.size(); i++)
        {
            res[i] *= right;
        }
        return res;
    }

    template <typename T>
    vector<T> operator*(const vector<T> &&left, const T &right)
    {
        vector<T> res(move(left));
        for (size_t i = 0; i < res.size(); i++)
        {
            res[i] *= right;
        }
        return res;
    }

    template <typename T>
    bool operator>(const vector<T> &left, const vector<T> &right)
    {
        for (size_t i = 0; i < left.size(); i++)
        {
            if (left[i] > right[i])
            {
            }
            else
            {
                return false;
            }
        }

        return true;
    }

    template <typename T>
    bool operator<(const vector<T> &left, const vector<T> &right)
    {
        for (size_t i = 0; i < left.size(); i++)
        {
            if (left[i] < right[i])
            {
            }
            else
            {
                return false;
            }
        }

        return true;
    }

    template <typename T>
    bool operator>=(const vector<T> &left, const vector<T> &right)
    {
        for (size_t i = 0; i < left.size(); i++)
        {
            if (left[i] >= right[i])
            {
            }
            else
            {
                return false;
            }
        }

        return true;
    }

    template <typename T>
    bool operator<=(const vector<T> &left, const vector<T> &right)
    {
        for (size_t i = 0; i < left.size(); i++)
        {
            if (left[i] <= right[i])
            {
            }
            else
            {
                return false;
            }
        }

        return true;
    }

}

#endif