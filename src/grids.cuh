//
// Created by leo on 6/27/23.
//

#include "data.cuh"
#include "utils.cuh"

#ifndef CUHPC_GRIDS_CUH
#define CUHPC_GRIDS_CUH

class grids : public data {

public:
    double dx, dy, dz;
    double h, dv;
    grids(Idim3 _NodeSize, Idim3 _GridSize, Idim3 _id, int _gid, bool pinned, DOMAIN_INFO _domain_info, bool pts);
    void as_point_data(DOMAIN_INFO _domain_info);
    void as_cell_data(DOMAIN_INFO _domain_info);
};


#endif //CUHPC_GRIDS_CUH
