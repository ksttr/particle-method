#include <linear_algebra/solver/matrix.h>
#include <fluid/mps/particle.h>
#include <fluid/mps/mps.h>

using namespace linear_algebra;
using namespace particle_method;

int main()
{
    matrix<double> A(2, 2, 0);
    A[0][0] = 1;
    A.print();

    Particle<double, 2> p(0);

    return 0;
}