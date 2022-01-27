#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include <vector>

using namespace std;

namespace particle_method
{
    template <typename T, size_t D>
    class Particle
    {
    private:
        bool dummy;
        bool wall;
        bool surface;

        T initial_density;

    public:
        T density;
        vector<T> position;
        vector<T> velocity;
        T pressure = 0;

        Particle();
        Particle(T density);
        Particle(T density, vector<T> position, vector<T> velocity);
        ~Particle();

        bool is_dummy();
        bool is_wall();
        bool is_surface(T beta = 0.97);
    };

    template <typename T, size_t D>
    Particle<T, D>::Particle()
        : position(D, 0), velocity(D, 0) {}

    template <typename T, size_t D>
    Particle<T, D>::Particle(T density)
        : initial_density(density), density(density), position(D, 0), velocity(D, 0) {}

    template <typename T, size_t D>
    Particle<T, D>::Particle(T density, vector<T> position, vector<T> velocity)
        : initial_density(density), density(density), position(position), velocity(velocity) {}

    template <typename T, size_t D>
    Particle<T, D>::~Particle()
    {
        position.clear();
        position.shrink_to_fit();
        velocity.clear();
        velocity.shrink_to_fit();
    }

    template <typename T, size_t D>
    bool Particle<T, D>::is_dummy()
    {
        return dummy;
    }

    template <typename T, size_t D>
    bool Particle<T, D>::is_wall()
    {
        return wall;
    }

    template <typename T, size_t D>
    bool Particle<T, D>::is_surface(T beta)
    {
        if (density < beta * initial_density)
            return true;
        else
            return false;
    }
}

#endif