#ifndef __MPS2D_H__
#define __MPS2D_H__

#include <cmath>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <cstdlib>

#include <linear_algebra/solver/matrix.h>
#include <linear_algebra/solver/vector_ex.h>
#include <linear_algebra/solver/solver.h>

#include "particle.h"
#include "filesio.h"

using namespace linear_algebra;
using namespace std;
namespace fs = filesystem;

namespace particle_method
{
    template <typename T>
    class mps2d
    {
    private:
        int time_step_count = 0;

        int D = 2;           // dimension
        size_t N;            // number of particles
        T l_0;               // particle distance at the initial time.
        T r_e_for_laplacian; // radius of influence for Laplacian model.
        T r_e_for_gradient;  // radius of influence for gradient model.
        T r_e_for_density;   // radius of influence for particle density.
        T collision_radius;
        T restitution_coefficient;

        matrix<T> weight;
        vector<Particle<T, 2>> particles;
        vector<T> pressure;

        T n0_for_density;
        T n0_for_gradient;
        T n0_for_laplacian;
        T lambda;

    public:
        const T VISCOSITY;
        const T FLUID_DENSITY;
        const vector<T> GRAVITY;
        const T COMPRESSIBILITY;

        const T dt;
        const size_t END_OF_TIME;

        mps2d(const T &VISCOSITY, const T &FLUID_DENSITY, const T &GRAVITATIONAL_ACCELERATION, const T &COMPRESSIBILITY, const T &dt, const size_t &END_OF_TIME);
        ~mps2d();
        size_t get_particle_number();
        void set_initial_particle_density(const T &n0);
        T get_initial_particle_distance();
        void set_initial_condition(const vector<Particle<T, 2>> &init_particles);
        T get_particle_min_distance();

        T calc_standard_density(const T r_e);
        T calc_lambda(const T r_e);

        T weight_value(T influence_distance, const Particle<T, 2> &particle1, const Particle<T, 2> &particle2);
        matrix<T> weight_function(T influence_distance);
        vector<T> weighted_average(const matrix<T> &f, const T &n0, const matrix<T> &weight_func);
        vector<T> calc_minimum_value(const vector<T> &func);

        matrix<T> calc_laplacian_model(const matrix<T> &func);
        matrix<T> calc_gradient_model(const vector<T> &func);

        void calc_particle_density(const matrix<T> &weight_func);
        void calc_pressure(const matrix<T> &weight_func);

        matrix<T> calc_collision_velocity(const T &collision_radius, const T &restitution_coefficient);

        void time_integral_first_step();
        void time_integral_second_step();
        void time_integral();
        void time_evolution(const fs::path &data_dir, const size_t print_count);
        bool check_finite();

        fs::path create_directory_for(const string &child_dir, const fs::path &parent_dir);
        void print_out(const fs::path &data_directory);
    };

    template <typename T>
    mps2d<T>::mps2d(const T &VISCOSITY, const T &FLUID_DENSITY, const T &GRAVITATIONAL_ACCELERATION, const T &COMPRESSIBILITY, const T &dt, const size_t &END_OF_TIME)
        : weight(1, 1, 0), VISCOSITY(VISCOSITY), FLUID_DENSITY(FLUID_DENSITY), GRAVITY{0, GRAVITATIONAL_ACCELERATION}, COMPRESSIBILITY(COMPRESSIBILITY), dt(dt), END_OF_TIME(END_OF_TIME) {}

    template <typename T>
    mps2d<T>::~mps2d()
    {
        weight.destroy();
        pressure.clear();
        pressure.shrink_to_fit();
        particles.clear();
        particles.shrink_to_fit();
    }

    template <typename T>
    size_t mps2d<T>::get_particle_number()
    {
        return N;
    }

    template <typename T>
    void mps2d<T>::set_initial_particle_density(const T &n0)
    {
        for (size_t i = 0; i < N; i++)
            particles[i].set_initial_density(n0);
    }

    template <typename T>
    T mps2d<T>::get_initial_particle_distance()
    {
        return l_0;
    }

    template <typename T>
    T mps2d<T>::get_particle_min_distance()
    {
        vector<T> diff_particles = particles[1].position - particles[0].position;
        T min_distance = sqrt(diff_particles * diff_particles);
        for (size_t i = 1; i < N - 1; i++)
        {
            diff_particles = particles[i + 1].position - particles[i].position;
            if (min_distance * min_distance > diff_particles * diff_particles)
                min_distance = sqrt(diff_particles * diff_particles);
        }

        return min_distance;
    }

    template <typename T>
    void mps2d<T>::set_initial_condition(const vector<Particle<T, 2>> &init_particles)
    {
        particles = init_particles;
        N = init_particles.size();

        pressure.resize(N);
        for (size_t i = 0; i < N; i++)
        {
            pressure[i] = init_particles[i].pressure;
        }

        l_0 = get_particle_min_distance();
        r_e_for_density = 2.1 * l_0;
        r_e_for_gradient = 2.1 * l_0;
        r_e_for_laplacian = 3.1 * l_0;
        collision_radius = 0.5 * l_0;
        restitution_coefficient = 0.2;
        weight = weight_function(r_e_for_density);

        n0_for_density = calc_standard_density(r_e_for_density);
        n0_for_gradient = calc_standard_density(r_e_for_gradient);
        n0_for_laplacian = calc_standard_density(r_e_for_laplacian);
        lambda = calc_lambda(r_e_for_laplacian);
        set_initial_particle_density(n0_for_density);
        calc_particle_density(weight);
        calc_pressure(weight);
    }

    template <typename T>
    void mps2d<T>::calc_particle_density(const matrix<T> &weight_func)
    {
        for (size_t i = 0; i < N; i++)
        {
            particles[i].density = 0;
            for (size_t j = 0; j < N; j++)
            {
                if (j == i)
                    continue;
                particles[i].density += weight_func[i][j];
            }
        }
    }

    template <typename T>
    T mps2d<T>::calc_standard_density(const T r_e)
    {
        T res_density = 0;
        size_t num = ceil(r_e / l_0);
        size_t particle_num = (2 * num + 1) * (2 * num + 1);
        size_t central_num = (2 * num + 1) * num + num;
        vector<Particle<T, 2>> init_particles(particle_num);
        for (size_t i = 0; i < particle_num; i++)
        {
            auto ix = i / (2 * num + 1);
            auto iy = i % (2 * num + 1);
            init_particles[i].position[0] = (static_cast<T>(ix) - static_cast<T>(num)) * l_0;
            init_particles[i].position[1] = (static_cast<T>(iy) - static_cast<T>(num)) * l_0;
        }

        for (size_t i = 0; i < particle_num; i++)
        {
            if (i == central_num)
                continue;

            res_density += weight_value(r_e, init_particles[central_num], init_particles[i]);
        }

        return res_density;
    }

    template <typename T>
    T mps2d<T>::calc_lambda(const T r_e)
    {
        T res_lambda = 0;
        T n0 = 0;
        size_t num = ceil(r_e / l_0);
        size_t particle_num = (2 * num + 1) * (2 * num + 1);
        size_t central_num = (2 * num + 1) * num + num;
        vector<Particle<T, 2>> init_particles(particle_num);
        for (size_t i = 0; i < particle_num; i++)
        {
            auto ix = i / (2 * num + 1);
            auto iy = i % (2 * num + 1);
            init_particles[i].position[0] = (static_cast<T>(ix) - static_cast<T>(num)) * l_0;
            init_particles[i].position[1] = (static_cast<T>(iy) - static_cast<T>(num)) * l_0;
        }

        for (size_t i = 0; i < particle_num; i++)
        {
            if (i == central_num)
                continue;

            auto distance = init_particles[central_num].position - init_particles[i].position;
            auto distance_square = distance * distance;
            auto w = weight_value(r_e, init_particles[central_num], init_particles[i]);
            n0 += w;
            res_lambda += distance_square * w;
        }
        res_lambda *= 1 / n0;
        return res_lambda;
    }

    template <typename T>
    vector<T> mps2d<T>::calc_minimum_value(const vector<T> &func)
    {
        vector<T> min_value(N, 0);

        for (size_t i = 0; i < N; i++)
        {
            if (particles[i].is_dummy())
                continue;

            min_value[i] = func[i];
            for (size_t j = 0; j < N; j++)
            {
                if (j == i)
                    continue;

                if (particles[j].is_dummy())
                    continue;

                auto diff_particles = particles[i].position - particles[j].position;
                auto distance = sqrt(diff_particles * diff_particles);

                if (distance >= r_e_for_gradient)
                    continue;

                if (min_value[i] > func[j])
                    min_value[i] = func[j];
            }
        }

        return min_value;
    }

    template <typename T>
    T mps2d<T>::weight_value(T influence_distance, const Particle<T, 2> &particle1, const Particle<T, 2> &particle2)
    {
        auto diff_particles = particle1.position - particle2.position;
        auto distance = sqrt(diff_particles * diff_particles);

        if (distance >= influence_distance)
            return 0.0;
        else
            return influence_distance / distance - 1;
    }

    template <typename T>
    matrix<T> mps2d<T>::weight_function(T influence_distance)
    {
        matrix<T> weight_func(N, N, 0);
        for (size_t i = 0; i < N; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                if (i == j)
                    continue;

                weight_func[i][j] = weight_value(influence_distance, particles[i], particles[j]);
            }
        }

        return weight_func;
    }

    template <typename T>
    vector<T> mps2d<T>::weighted_average(const matrix<T> &f, const T &n0, const matrix<T> &weight_func)
    {
        vector<T> average(N, 0);
        for (size_t i = 0; i < average.size(); i++)
        {
            for (size_t j = 0; j < average.size(); j++)
            {
                if (j == i)
                    continue;

                average[i] += f[i][j] * weight_func[i][j];
            }
            average[i] *= 1 / n0;
        }
        return average;
    }

    template <typename T>
    matrix<T> mps2d<T>::calc_collision_velocity(const T &collision_radius, const T &restitution_coefficient)
    {
        matrix<T> collision_velocity(N, 2, 0);
        for (size_t i = 0; i < N; i++)
        {
            if (!particles[i].is_fluid())
                continue;

            collision_velocity[i] = particles[i].velocity;

            for (size_t j = 0; j < N; j++)
            {
                if (i == j)
                    continue;

                auto relative_velocity = particles[j].velocity - particles[i].velocity;
                auto relative_position = particles[j].position - particles[i].position;
                auto distance = sqrt(relative_position * relative_position);

                if (distance < collision_radius)
                {
                    auto relative_velocity_collision_direction = -relative_velocity * relative_position / distance;
                    collision_velocity[i] -= (1 + restitution_coefficient) / 2 * relative_velocity_collision_direction / distance * relative_position;
                }
            }
        }
        return collision_velocity;
    }

    template <typename T>
    matrix<T> mps2d<T>::calc_gradient_model(const vector<T> &func)
    {
        matrix<T> res(2, N, 0);
        matrix<T> grad_x_func(N, N, 0);
        matrix<T> grad_y_func(N, N, 0);
        weight = weight_function(r_e_for_gradient);
        auto min_func = calc_minimum_value(func);
        for (size_t i = 0; i < N; i++)
        {
            if (!particles[i].is_fluid())
                continue;

            for (size_t j = 0; j < N; j++)
            {
                if (j == i)
                    continue;

                if (particles[j].is_dummy())
                    continue;

                auto diff = particles[j].position - particles[i].position;
                auto distance_square = diff * diff;
                grad_x_func[i][j] = static_cast<T>(D) * (func[j] - min_func[i]) * (diff[0] / distance_square);
                grad_y_func[i][j] = static_cast<T>(D) * (func[j] - min_func[i]) * (diff[1] / distance_square);
            }

            res[0] = weighted_average(grad_x_func, n0_for_gradient, weight);
            res[1] = weighted_average(grad_y_func, n0_for_gradient, weight);
        }

        return res.transpose();
    }

    template <typename T>
    matrix<T> mps2d<T>::calc_laplacian_model(const matrix<T> &func)
    {
        matrix<T> res(2, N, 0);
        matrix<T> Delta_func_x(N, N, 0);
        matrix<T> Delta_func_y(N, N, 0);
        weight = weight_function(r_e_for_laplacian);
        T coefficient = 2 * static_cast<T>(D) / lambda;
        for (size_t i = 0; i < N; i++)
        {
            if (!particles[i].is_fluid())
                continue;

            for (size_t j = 0; j < N; j++)
            {
                if (j == i)
                    continue;

                Delta_func_x[i][j] = coefficient * (func[j][0] - func[i][0]);
                Delta_func_y[i][j] = coefficient * (func[j][1] - func[i][1]);
            }
        }
        res[0] = weighted_average(Delta_func_x, n0_for_laplacian, weight);
        res[1] = weighted_average(Delta_func_y, n0_for_laplacian, weight);

        return res.transpose();
    }

    template <typename T>
    void mps2d<T>::calc_pressure(const matrix<T> &weight_func)
    {
        matrix<T> A(N, N, 0);
        vector<T> b(N, 0);

        vector<T> x_0(pressure);
        T precision = 1.0e-6;

        T coefficient = 2 * static_cast<T>(D) / (n0_for_laplacian * lambda);
        T stability_coefficient = 0.2;

        for (size_t i = 0; i < N; i++)
        {
            if (particles[i].is_inner_particle())
                b[i] = stability_coefficient / (dt * dt) * ((particles[i].density - n0_for_density) / n0_for_density);

            for (size_t j = 0; j < N; j++)
            {
                if (i != j)
                {
                    if (particles[i].is_dummy() | particles[j].is_dummy() | particles[i].is_surface())
                    {
                        A[i][j] = 0.0;
                        A[i][i] += coefficient * weight_func[i][j] / FLUID_DENSITY;
                    }
                    else if (particles[i].is_inner_particle())
                    {
                        A[i][j] = -coefficient * weight_func[i][j] / FLUID_DENSITY;
                        A[i][i] += coefficient * weight_func[i][j] / FLUID_DENSITY;
                    }
                }
            }
            A[i][i] += COMPRESSIBILITY / (dt * dt);
        }

        auto solver = Solver<T>(A, b);
        pressure = solver.incomplete_cholesky_factorization_conjugate_gradient(x_0, precision);

        for (size_t i = 0; i < N; i++)
        {
            if (pressure[i] < 0.0)
                pressure[i] = 0.0;

            particles[i].pressure = pressure[i];
        }
    }

    template <typename T>
    void mps2d<T>::time_integral_first_step()
    {
        matrix<T> velocity_field(N, 2, 0);
        for (size_t i = 0; i < N; i++)
        {

            velocity_field[i] = particles[i].velocity;
        }

        auto nu_laplace_u = VISCOSITY * calc_laplacian_model(velocity_field);

        for (size_t i = 0; i < N; i++)
        {
            if (particles[i].is_fluid())
            {
                particles[i].velocity += (nu_laplace_u[i] + GRAVITY) * dt;
                particles[i].position += particles[i].velocity * dt;
            }
        }

        auto collision_velocity = calc_collision_velocity(collision_radius, restitution_coefficient);

        for (size_t i = 0; i < N; i++)
        {
            if (particles[i].is_fluid())
            {
                particles[i].position += (collision_velocity[i] - particles[i].velocity) * dt;
                particles[i].velocity = collision_velocity[i];
            }
        }
    }

    template <typename T>
    void mps2d<T>::time_integral_second_step()
    {
        weight = weight_function(r_e_for_density);
        calc_particle_density(weight);
        calc_pressure(weight);
        auto grad_pressure_per_density = calc_gradient_model(pressure) / FLUID_DENSITY;

        for (size_t i = 0; i < N; i++)
        {
            if (particles[i].is_fluid())
            {
                particles[i].velocity -= grad_pressure_per_density[i] * dt;
                particles[i].position -= grad_pressure_per_density[i] * dt * dt;
            }
        }
    }

    template <typename T>
    void mps2d<T>::time_integral()
    {
        time_integral_first_step();
        time_integral_second_step();
    }

    template <typename T>
    bool mps2d<T>::check_finite()
    {
        for (size_t i = 0; i < N; i++)
        {
            if (isfinite(particles[i].pressure))
            {
                continue;
            }
            else
            {
                return false;
            }
        }

        return true;
    }

    template <typename T>
    void mps2d<T>::time_evolution(const fs::path &data_dir, const size_t print_count)
    {
        print_out(data_dir);

        while (time_step_count < END_OF_TIME)
        {
            time_integral();
            ++time_step_count;

            if (time_step_count / print_count * print_count == time_step_count)
                print_out(data_dir);

            if (!check_finite())
            {
                cerr << "ERROR: Simulation is broken (simulation count: " << time_step_count << " )." << endl;
                exit(1);
            }
        }
    }

    template <typename T>
    fs::path mps2d<T>::create_directory_for(const string &child_dir, const fs::path &parent_dir)
    {
        fs::path child = child_dir;
        auto filepath = parent_dir / child;

        if (time_step_count == 0)
        {
            vector<string> sub_dir_array = {};
            create_directories(filepath, sub_dir_array);
        }

        fs::path filename = child_dir + string_format("%06d", time_step_count) + ".dat";
        filepath /= filename;

        return filepath;
    }

    template <typename T>
    void mps2d<T>::print_out(const fs::path &data_directory)
    {
        cout << "StepCount: " << time_step_count << " "
             << "ParticleNumber: " << N << endl;

        auto particle_type_path = create_directory_for("type", data_directory);
        auto pressure_path = create_directory_for("pressure", data_directory);
        auto position_path = create_directory_for("position", data_directory);
        auto velocity_path = create_directory_for("velocity", data_directory);

        ofstream out_stream_particle_type(particle_type_path, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
        for (size_t n = 0; n < N; n++)
        {
            auto type_value = static_cast<int>(particles[n].get_type());
            out_stream_particle_type.write(reinterpret_cast<char *>(&type_value), sizeof(int));
        }

        ofstream out_stream_pressure(pressure_path, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
        for (size_t n = 0; n < N; n++)
        {
            out_stream_pressure.write(reinterpret_cast<char *>(&particles[n].pressure), sizeof(decltype(particles[n].pressure)));
        }

        ofstream out_stream_position(position_path, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
        for (size_t n = 0; n < N; n++)
        {
            out_stream_position.write(reinterpret_cast<char *>(&particles[n].position[0]), particles[n].position.size() * sizeof(T));
        }

        ofstream out_stream_velocity(velocity_path, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
        for (size_t n = 0; n < N; n++)
        {
            out_stream_velocity.write(reinterpret_cast<char *>(&particles[n].velocity[0]), particles[n].velocity.size() * sizeof(T));
        }
    }
}

#endif