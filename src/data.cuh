//
// Created by 溫晧良 on 2023/6/21.
//

#include <cuda.h>
#include <assert.h>
#include <iostream>
#include "utils.cuh"
#include "constants.cuh"

#ifndef CUHPC_DATA_CUH
#define CUHPC_DATA_CUH

class data {

private:

    void XL_NEUMANN(int d);
    void XL_DIRICHLET(int d);
    void XL_PERIODIC(int d);
    void XL_NO_SLIP(int d);
    void XL_FREE_SLIP(int d);
    void XR_NEUMANN(int d);
    void XR_DIRICHLET(int d);
    void XR_PERIODIC(int d);
    void XR_NO_SLIP(int d);
    void XR_FREE_SLIP(int d);

public:

    Idim3 NodeSize, GridSize, id;
    int gid, dim, spatial_dim;
    int size;
    double *host, *device;

    data *devptr;

    data **neighbor;
    BC::INFO **bc_type;

    data(Idim3 _NodeSize, Idim3 _GridSize, Idim3 _id, int _gid, int _dim, bool pinned);

    void LinkNeighbor();
    void LinkNeighbor(data *_r, data *_l);
    void LinkNeighbor(data *_r, data *_l, data *_u, data *_d);
    void LinkNeighbor(data *_r, data *_l, data *_u, data *_d, data *_f, data *_b);

    void set(double val, int i, int j, int k, int d);
    void set(double val, int i, int j, int k);

    int map_index(int i, int j, int k, int d);

    double get(int i, int j, int k, int d);
    double get(int i, int j, int k);

    void assign_bc(int d, BC::INFO xr, BC::INFO xl);
    void assign_bc(int d, BC::INFO xr, BC::INFO xl, BC::INFO yr, BC::INFO yl);
    void assign_bc(int d, BC::INFO xr, BC::INFO xl, BC::INFO yr, BC::INFO yl, BC::INFO zr, BC::INFO zl);

    void apply_bc_x(int d);
    void apply_bc_y(int d);
    void apply_bc_z(int d);

    void ToDevice();
    void ToHost();

};


#endif //CUHPC_DATA_CUH
