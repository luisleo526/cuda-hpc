//
// Created by leo on 6/27/23.
//

#include "data.cuh"
#include "grids.cuh"
#include "ccd.cuh"
#include "utils.cuh"
#include "constants.cuh"

#ifndef CUHPC_VARIABLE_CUH
#define CUHPC_VARIABLE_CUH

class variable{

public:
    int derivative_level;
    int size;
    Idim3 NodeSize, GridSize, id;
    int gid, dim, spatial_dim;
    data **f, **df, **ddf, **flux;
    grids *grid;
    ccd *ccd_solver;
    variable(Idim3 _NodeSize, Idim3 _GridSize, Idim3 _id, int _gid, int _dim, bool pinned, int _derivative_level, grids *_grid);
    int map_index(int i, int j, int k, int d);
    void LinkNeighbor();
    void LinkNeighbor(variable *_r, variable *_l);
    void LinkNeighbor(variable *_r, variable *_l, variable *_u, variable *_d);
    void LinkNeighbor(variable *_r, variable *_l, variable *_u, variable *_d, variable *_f, variable *_b);

    void get_derivative_SEC(int level, int d, int direction);
    void assign_CCD_source(int d, int direction, int type, int I, int J, int K);
    void get_derivative_CCD(int d, int direction);
    void get_derivative_UCCD(int d, int direction, variable *vel);
    void get_derivative_USEC(int d, int direction, variable *vel);

};


#endif //CUHPC_VARIABLE_CUH
