#ifndef VPRMODELIO_H
#define VPRMODELIO_H

#include <vector>
#include "core/konal_constants.hpp"

// function: read vpr format model to a 1D vetcor
int readvpr(std::vector<Ckonal::real_t> &model_data,
    std::vector<Ckonal::real_t> &coordsmin, std::vector<Ckonal::uint_t> &coordsnum, std::vector<Ckonal::real_t> &coordsstep,
    std::string input_vprfile);

//

#endif