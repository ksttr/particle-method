#include <linear_algebra/solver/matrix.h>
#include "particle.h"
#include "mps2d.h"
#include "initial_condition.h"
#include "filesio.h"

#include <string>

using namespace std;
using namespace linear_algebra;
using namespace particle_method;

int main()
{
    const double VISCOSITY = 0.000001;
    const double FLUID_DENSITY = 1000.0;
    const double GRAVITATIONAL_ACCELERATION = -9.8;
    const double COMPRESSIBILITY = 5e-10;
    const int END_OF_TIME = 2000;
    const double dt = 0.001;
    const double Lx = 0.6;
    const double Ly = 0.6;
    const double lx = 0.20;
    const double ly = 0.45;
    const size_t Nx = 50;
    const size_t Ny = 50;

    mps2d<double> mps_obj(VISCOSITY, FLUID_DENSITY, GRAVITATIONAL_ACCELERATION, COMPRESSIBILITY, dt, END_OF_TIME);
    InitialCondition2D<double> initial_condition(Nx, Ny);
    auto particles = initial_condition.dam_break(Lx, Ly, lx, ly);
    auto data_dir = initial_condition.create_data_directories_for_dam_break(Lx, Ly, lx, ly, dt);
    mps_obj.set_initial_condition(particles);
    mps_obj.time_evolution(data_dir, 10);

    cout << mps_obj.get_initial_particle_distance() << endl;

    return 0;
}