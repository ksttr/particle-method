#include <iostream>
#include <string>

#include "filesio.h"

namespace fs = std::filesystem;

namespace particle_method
{

    void create_directories(fs::path parent_directory, const std::vector<std::string> &sub_directories)
    {
        if (sub_directories.empty())
        {
            if (!fs::is_directory(parent_directory))
                fs::create_directories(parent_directory);
        }
        else
        {

            for (const auto &element : sub_directories)
            {
                fs::path current_path = parent_directory;
                current_path += element;

                if (!fs::is_directory(current_path))
                    fs::create_directories(current_path);
            }
        }
    }
}