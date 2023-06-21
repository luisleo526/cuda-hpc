//
// Created by 溫晧良 on 2023/6/21.
//

#include <cuda.h>
#include <assert.h>

#ifndef CUHPC_DATA_CUH
#define CUHPC_DATA_CUH


class data {

public:

    dim3 NodeSize, GridSize, id;
    int gid, dim, spatial_dim;
    double **host, **device;

    data *neir, *neil, *neiu, *neid, *neif, *neib;
    int **bc_type;

    data(dim3 _NodeSize, dim3 _GridSize, dim3 _id, int _gid, int _dim);

    void LinkNeighbor(data *_r, data *_l, data *_u, data *_d, data *_f, data *_b);

    void AssignBC(int d, int bc_r, int bc_l);

    void AssignBC(int d, int bc_r, int bc_l, int bc_f, int bc_b);

    void AssignBC(int d, int bc_r, int bc_l, int bc_f, int bc_b, int bc_u, int bc_d);

    void set(double val, int i, int j, int k, int d);

    double get(int i, int j, int k, int d);

    void apply_bc_xr(int dim);

    void apply_bc_xl(int d);

    void apply_bc_yr(int d);

    void apply_bc_yl(int d);

    void apply_bc_zr(int d);

    void apply_bc_zl(int d);

};


#endif //CUHPC_DATA_CUH
