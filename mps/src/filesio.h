#ifndef __FILESIO_H__
#define __FILESIO_H__

#include <string>
#include <cstdio>
#include <vector>
#include <cassert>
#include <filesystem>
#include <regex>

namespace fs = std::filesystem;

namespace particle_method
{
    template <typename... Args>
    std::string string_format(const std::string &fmt, Args... args);

    template <typename T>
    fs::path parent_directory_name(const double &Lx, const double &Ly, const double &lx, const double &ly, const double &dt, std::string data_root = "./data/");

    void create_directories(fs::path parent_directory, const std::vector<std::string> &sub_directories);

    template <typename... Args>
    std::string string_format(const std::string &fmt, Args... args)
    {
        size_t len = std::snprintf(nullptr, 0, fmt.c_str(), args...);
        std::vector<char> buf(len + 1);
        std::snprintf(&buf[0], len + 1, fmt.c_str(), args...);
        return std::string(&buf[0], &buf[0] + len);
    }

    template <typename T>
    fs::path parent_directory_name(const T &Lx, const T &Ly, const T &lx, const T &ly, const T &dt, std::string data_root)
    {
        fs::path parent_directory = data_root;

        std::regex re("\\.");
        auto Lx_string = std::regex_replace(std::to_string(Lx), re, "");
        auto Ly_string = std::regex_replace(std::to_string(Ly), re, "");
        auto lx_string = std::regex_replace(std::to_string(lx), re, "");
        auto ly_string = std::regex_replace(std::to_string(ly), re, "");
        auto dt_string = std::regex_replace(std::to_string(dt), re, "");

        fs::path dirname = "Lx" + Lx_string + "_Ly" + Ly_string + "_lx" + lx_string + "_ly" + ly_string + "_dt" + dt_string;

        parent_directory /= dirname;

        return parent_directory;
    }
}

#endif