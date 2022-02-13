#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include <vector>

using namespace std;

namespace particle_method
{
    enum class particleType
    {
        Fluid,
        Wall,
        Dummy,
        Ghost // for inflow/outflow boundary condition
    };

    template <typename T, int D>
    class Particle
    {
    private:
        particleType type;

        size_t bucket_number;

        T initial_density;

    public:
        T density;
        vector<T> position;
        vector<T> velocity;
        T pressure = 0;

        Particle();
        Particle(T density);
        Particle(vector<T> position, vector<T> velocity);
        Particle(T density, vector<T> position, vector<T> velocity);
        ~Particle();

        void set_initial_density(const T &particle_density);
        T get_initial_density();

        void set_type(particleType particle_type);
        particleType get_type();

        bool is_dummy();
        bool is_wall();
        bool is_fluid();
        bool is_inner_particle(T beta = 0.97);
        bool is_surface(T beta = 0.97);
    };

    template <typename T, int D>
    Particle<T, D>::Particle()
        : position(D, 0), velocity(D, 0) {}

    template <typename T, int D>
    Particle<T, D>::Particle(vector<T> position, vector<T> velocity)
        : position(position), velocity(velocity) {}

    template <typename T, int D>
    Particle<T, D>::Particle(T density)
        : initial_density(density), density(density), position(D, 0), velocity(D, 0) {}

    template <typename T, int D>
    Particle<T, D>::Particle(T density, vector<T> position, vector<T> velocity)
        : initial_density(density), density(density), position(position), velocity(velocity) {}

    template <typename T, int D>
    Particle<T, D>::~Particle()
    {
        position.clear();
        position.shrink_to_fit();
        velocity.clear();
        velocity.shrink_to_fit();
    }

    template <typename T, int D>
    void Particle<T, D>::set_initial_density(const T &particle_density)
    {
        initial_density = particle_density;
    }

    template <typename T, int D>
    T Particle<T, D>::get_initial_density()
    {
        return initial_density;
    }

    template <typename T, int D>
    void Particle<T, D>::set_type(particleType particle_type)
    {
        type = particle_type;
    }

    template <typename T, int D>
    particleType Particle<T, D>::get_type()
    {
        return type;
    }

    template <typename T, int D>
    bool Particle<T, D>::is_dummy()
    {
        if (type == particleType::Dummy)
            return true;
        else
            return false;
    }

    template <typename T, int D>
    bool Particle<T, D>::is_wall()
    {
        if (type == particleType::Wall)
            return true;
        else
            return false;
    }

    template <typename T, int D>
    bool Particle<T, D>::is_fluid()
    {
        if (type == particleType::Fluid)
            return true;
        else
            return false;
    }

    template <typename T, int D>
    bool Particle<T, D>::is_surface(T beta)
    {
        if (type == particleType::Dummy || type == particleType::Ghost)
            return false;

        if (density < beta * initial_density)
            return true;
        else
            return false;
    }

    template <typename T, int D>
    bool Particle<T, D>::is_inner_particle(T beta)
    {
        if (type == particleType::Dummy || type == particleType::Ghost)
            return false;

        if (density < beta * initial_density)
            return false;
        else
            return true;
    }
}

#endif