//
// Created by 溫晧良 on 2023/6/21.
//

#include <cuda.h>
#include <assert.h>
#include <iostream>

#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)

#ifndef CUHPC_DATA_CUH
#define CUHPC_DATA_CUH

namespace BC {
    int const CELL_CENTER_NEUMANN = 0;
    int const CELL_CENTER_DIRICHLET = 1;
    int const CELL_CENTER_NO_SLIP = 2;
    int const CELL_CENTER_FREE_SLIP = 3;
    int const CELL_CENTER_PERIODIC = 4;

    int const CELL_FACE_X_NEUMANN = 10;
    int const CELL_FACE_X_DIRICHLET = 11;
    int const CELL_FACE_X_NO_SLIP = 12;
    int const CELL_FACE_X_FREE_SLIP = 13;
    int const CELL_FACE_X_PERIODIC = 14;

    int const CELL_FACE_Y_NEUMANN = 20;
    int const CELL_FACE_Y_DIRICHLET = 21;
    int const CELL_FACE_Y_NO_SLIP = 22;
    int const CELL_FACE_Y_FREE_SLIP = 23;
    int const CELL_FACE_Y_PERIODIC = 24;

    int const CELL_FACE_Z_NEUMANN = 30;
    int const CELL_FACE_Z_DIRICHLET = 31;
    int const CELL_FACE_Z_NO_SLIP = 32;
    int const CELL_FACE_Z_FREE_SLIP = 33;
    int const CELL_FACE_Z_PERIODIC = 34;

    struct INFO {
        int type = CELL_CENTER_NEUMANN;
        double value = 0.0;

        INFO() : type(CELL_CENTER_NEUMANN), value(0.0) {};

        INFO(int t, double v) : type(t), value(v) {};
        INFO(int t) : type(t), value(0.0) {};
    };
}

struct Idim3 {
    int x, y, z;

    Idim3() : x(1), y(1), z(1) {};

    Idim3(int i, int j, int k) : x(i), y(j), z(k) {};
};

class data {

private:
    int const XL = 0;
    int const XR = 1;
    int const YL = 2;
    int const YR = 3;
    int const ZL = 4;
    int const ZR = 5;

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

    int map_index(int i, int j, int k, int d);

    double get(int i, int j, int k, int d);

    void assign_bc(int d, BC::INFO xr, BC::INFO xl);

    void assign_bc(int d, BC::INFO xr, BC::INFO xl, BC::INFO yr, BC::INFO yl);

    void assign_bc(int d, BC::INFO xr, BC::INFO xl, BC::INFO yr, BC::INFO yl, BC::INFO zr, BC::INFO zl);

    void apply_bc_x(int d);

    void apply_bc_y(int d);

    void apply_bc_z(int d);

    void ToDevice();

};


#endif //CUHPC_DATA_CUH
