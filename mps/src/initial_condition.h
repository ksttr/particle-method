#ifndef __INITIAL_CONDITION_H__
#define __INITIAL_CONDITION_H__

#include <iostream>
#include <vector>
#include <filesystem>
#include "particle.h"
#include "filesio.h"
#include <linear_algebra/solver/matrix.h>

using namespace std;
using namespace linear_algebra;
namespace fs = filesystem;

namespace particle_method
{
    template <typename T>
    class InitialCondition2D
    {
    private:
        const size_t Nx;
        const size_t Ny;
        size_t N;
        vector<Particle<T, 2>> particles;

    public:
        InitialCondition2D(const size_t Nx, const size_t Ny);
        vector<Particle<T, 2>> &dam_break(const T Lx, const T Ly, const T lx, const T ly);
        fs::path create_data_directories_for_dam_break(const T &Lx, const T &Ly, const T &lx, const T &ly, const T &dt, const string &data_root = "./data");
    };

    template <typename T>
    InitialCondition2D<T>::InitialCondition2D(const size_t Nx, const size_t Ny)
        : Nx(Nx), Ny(Ny)
    {
    }

    template <typename T>
    vector<Particle<T, 2>> &InitialCondition2D<T>::dam_break(const T Lx, const T Ly, const T lx, const T ly)
    {
        T dx = Lx / Nx;
        T dy = Ly / Ny;

        for (size_t ix = 0; ix < Nx + 8; ix++)
        {
            for (size_t iy = 0; iy < Ny + 8; iy++)
            {
                vector<T> position(2, 0);
                vector<T> velocity(2, 0);
                position[0] = (static_cast<T>(ix) - 4) * dx;
                position[1] = (static_cast<T>(iy) - 4) * dy;

                if (ix <= 1 || ix >= Nx + 6 || iy <= 1)
                {
                    Particle<T, 2> particle(position, velocity);
                    particle.set_type(particleType::Dummy);
                    particles.push_back(particle);
                }
                else if (ix <= 3 || ix >= Nx + 4 || iy <= 3)
                {
                    Particle<T, 2> particle(position, velocity);
                    particle.set_type(particleType::Wall);
                    particles.push_back(particle);
                }
                else if (position[0] < lx && position[1] < ly)
                {
                    Particle<T, 2> particle(position, velocity);
                    particle.set_type(particleType::Fluid);
                    particles.push_back(particle);
                }
            }
        }

        return particles;
    }

    template <typename T>
    fs::path InitialCondition2D<T>::create_data_directories_for_dam_break(const T &Lx, const T &Ly, const T &lx, const T &ly, const T &dt, const string &data_root)
    {
        auto parent_directory = parent_directory_name(Lx, Ly, lx, ly, dt, data_root);

        vector<string> sub_dir_array = {};

        create_directories(parent_directory, sub_dir_array);

        return parent_directory;
    }
}

#endif