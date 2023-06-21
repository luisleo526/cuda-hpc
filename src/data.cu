//
// Created by 溫晧良 on 2023/6/21.
//

#include "data.cuh"

data::data(dim3 _NodeSize, dim3 _GridSize, dim3 _id, int _gid, int _dim) {

    NodeSize = dim3(_NodeSize.x, _NodeSize.y, _NodeSize.z);
    GridSize = dim3(_GridSize.x, _GridSize.y, _GridSize.z);
    id = dim3(_id.x, _id.y, _id.z);

    gid = _gid;
    dim = _dim;

    spatial_dim = 0;
    if (NodeSize.x > 1) spatial_dim++;
    if (NodeSize.y > 1) spatial_dim++;
    if (NodeSize.z > 1) spatial_dim++;

    host = new double *[dim];
    device = new double *[dim];
    bc_type = new int *[dim];

    cudaSetDevice(gid);
    for (int i = 0; i < dim; i++) {
        host[i] = new double[NodeSize.x * NodeSize.y * NodeSize.z];
        cudaMalloc((void **) &device[i], NodeSize.x * NodeSize.y * NodeSize.z * sizeof(double));
    }
    cudaDeviceReset();

}

void data::set(double val, int i, int j, int k, int d) {
    host[d][i + j * NodeSize.x + k * NodeSize.x * NodeSize.y] = val;
}

double data::get(int i, int j, int k, int d) {


}

void data::LinkNeighbor(data *_r, data *_l, data *_u, data *_d, data *_f, data *_b) {
    neir = _r;
    neil = _l;
    neiu = _u;
    neid = _d;
    neif = _f;
    neib = _b;
}

void data::AssignBC(int d, int bc_r, int bc_l) {
    assert(spatial_dim == 1);
    bc_type[d] = new int[6]{bc_r, bc_l, -1, -1, -1, -1};
}

void data::AssignBC(int d, int bc_r, int bc_l, int bc_f, int bc_b) {
    assert(spatial_dim == 2);
    bc_type[d] = new int[6]{bc_r, bc_l, bc_f, bc_b, -1, -1};
}

void data::AssignBC(int d, int bc_r, int bc_l, int bc_f, int bc_b, int bc_u, int bc_d) {
    assert(spatial_dim == 3);
    bc_type[d] = new int[6]{bc_r, bc_l, bc_f, bc_b, bc_u, bc_d};
}
