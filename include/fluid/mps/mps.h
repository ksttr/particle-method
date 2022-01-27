#ifndef __MPS_H__
#define __MPS_H__

#include <fluid/mps/particle.h>
#include <linear_algebra/solver/matrix.h>

using namespace linear_algebra;

namespace particle_method
{
    template <typename T, size_t D>
    class mps
    {
    private:
        matrix<T> pressure_field;

    public:
    };

}

#endif