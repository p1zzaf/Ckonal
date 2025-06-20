#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include "vprmodelio.h"
// #include "../core/konal_constants.hpp" 

int readvpr(std::vector<Ckonal::real_t> &model_data,
    std::vector<Ckonal::real_t> &coordsmin, std::vector<Ckonal::uint_t> &coordsnum, std::vector<Ckonal::real_t> &coordsstep,
    std::string input_vprfile){
    
        // 打开文件
    std::ifstream file(input_vprfile);
    if (!file.is_open()) {
        std::cerr << "Failed to open " << input_vprfile << "!!" << std::endl;
        return -1;
    }

    // read data
    std::string str1line;
    Ckonal::real_t vmin, vmax;
    Ckonal::real_t x_min, y_min, z_min, x_step, y_step, z_step;
    Ckonal::uint_t nx, ny, nz;
    Ckonal::real_t var1;
    
    // read first line
    std::getline(file, str1line);
    std::istringstream iss(str1line);
    iss >> vmin >> vmax;

    // read second line
    std::getline(file, str1line);
    std::string vprfile;
    iss.str(str1line);
    iss >> vprfile;

    // read third line
    std::getline(file, str1line);
    iss.str(str1line);
    iss >> x_min >> x_step >> nx >> y_min >> y_step >> ny >> z_min >> z_step >> nz;
    coordsmin.push_back(x_min);
    coordsmin.push_back(y_min);
    coordsmin.push_back(z_min);
    coordsstep.push_back(x_step);
    coordsstep.push_back(y_step);
    coordsstep.push_back(z_step);
    coordsnum.push_back(nx);
    coordsnum.push_back(ny);
    coordsnum.push_back(nz);

    for(Ckonal::uint_t ix=0; ix<nx; ix++){
        for(Ckonal::uint_t iy=0; iy<ny; iy++){
            std::getline(file, str1line);
            iss.str(str1line);
            for(Ckonal::uint_t iz=0; iz<nz; iz++){
                iss >> var1;
                model_data.push_back(var1);
            }
        }
        // skip blank line
        std::getline(file, str1line);
    }
    
    std::cout << "**readvpr** read (" << nx << ", " << ny << ", " << nz << ") data from " << input_vprfile << " file." << std::endl;

    return 0;
}

