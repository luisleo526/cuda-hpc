//
// Created by leo on 6/27/23.
//

#include "grids.cuh"
#include "utils.cuh"

#ifndef CUHPC_CCD_CUH
#define CUHPC_CCD_CUH

class ccd {

public:
    const int NUM_OF_DIMS = 3;
    const int NUM_OF_TYPES = 3;

    int gid;
    int *Nodesize;
    double *mesh_size;
    double ****A, ****B, ****AA, ****BB, ***S, ***SS;

    ccd(Idim3 _Nodesize, int _gid, grids *_grid);

};


#endif //CUHPC_CCD_CUH
