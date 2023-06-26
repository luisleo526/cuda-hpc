//
// Created by 溫晧良 on 2023/6/21.
//

#include "data.cuh"


data::data(Idim3 _NodeSize, Idim3 _GridSize, Idim3 _id, int _gid, int _dim, bool pinned) {

    NodeSize = Idim3(_NodeSize.x, _NodeSize.y, _NodeSize.z);
    Idim3 oNodeSize = Idim3(_NodeSize.x, _NodeSize.y, _NodeSize.z);
    GridSize = Idim3(_GridSize.x, _GridSize.y, _GridSize.z);
    id = Idim3(_id.x, _id.y, _id.z);

    gid = _gid;
    dim = _dim;

    spatial_dim = 0;
    if (NodeSize.x > 1) {
        spatial_dim++;
        oNodeSize.x += 6;
    }
    if (NodeSize.y > 1) {
        spatial_dim++;
        oNodeSize.y += 6;
    }
    if (NodeSize.z > 1) {
        spatial_dim++;
        oNodeSize.z += 6;
    }

    size = oNodeSize.x * oNodeSize.y * oNodeSize.z * sizeof(double);
    bc_type = (BC::INFO **) malloc(dim * sizeof(BC::INFO *));
    neighbor = (data **) malloc(6 * sizeof(data *));

    for (int i = 0; i < dim; i++) bc_type[i] = (BC::INFO *) malloc(6 * sizeof(BC::INFO));

    if (!pinned) {
        host = (double *) malloc(dim * size);
    } else {
        cudaMallocHost((void **) &host, dim * size);
        cudaCheckErrors("cudaMallocHost failed");
    }

    cudaSetDevice(gid);
    cudaCheckErrors("cudaSetDevice failed");
    cudaMalloc((void **) &device, dim * size);
    cudaCheckErrors("cudaMalloc failed");

}

int data::map_index(int i, int j, int k, int d) {
/*               |                         |
 *    -3, -2, -1 | 0, 1, 2, ...,  N-2, N-1 | N, N+1, N+2
 *               |                         |
 */

    assert((i > -4) && (i < NodeSize.x + 3));
    assert((j > -4 && j < NodeSize.y + 3 && spatial_dim > 1) || (j == 0 && spatial_dim < 2));
    assert((k > -4 && j < NodeSize.z + 3 && spatial_dim > 2) || (k == 0 && spatial_dim < 3));
    assert(d < dim);

    int I = i + 3;
    int J = spatial_dim > 1 ? j + 3 : j;
    int K = spatial_dim > 2 ? k + 3 : k;

//    std::cout << "(" << I << "," << J << "," << K << ")" << std::endl;
    return I + J * (NodeSize.x + 6) + K * (NodeSize.x + 6) * (NodeSize.y + 6) + d * size / sizeof(double);
}

void data::set(double val, int i, int j, int k, int d) {
    host[map_index(i, j, k, d)] = val;
}

double data::get(int i, int j, int k, int d) {
    return host[map_index(i, j, k, d)];
}

void data::LinkNeighbor() {
    neighbor[XL] = this;
    neighbor[XR] = this;
    neighbor[YL] = this;
    neighbor[YR] = this;
    neighbor[ZL] = this;
    neighbor[ZR] = this;
}

void data::LinkNeighbor(data *_r, data *_l) {
    neighbor[XL] = _l;
    neighbor[XR] = _r;
    neighbor[YL] = this;
    neighbor[YR] = this;
    neighbor[ZL] = this;
    neighbor[ZR] = this;
}

void data::LinkNeighbor(data *_r, data *_l, data *_u, data *_d) {
    neighbor[XL] = _l;
    neighbor[XR] = _r;
    neighbor[YL] = _d;
    neighbor[YR] = _u;
    neighbor[ZL] = this;
    neighbor[ZR] = this;
}

void data::LinkNeighbor(data *_r, data *_l, data *_u, data *_d, data *_f, data *_b) {
    neighbor[XL] = _l;
    neighbor[XR] = _r;
    neighbor[YL] = _d;
    neighbor[YR] = _u;
    neighbor[ZL] = _b;
    neighbor[ZR] = _f;
}

void data::assign_bc(int d, BC::INFO xr, BC::INFO xl) {
    BC::INFO buffer;
    bc_type[d][XL] = xl;
    bc_type[d][XR] = xr;
    bc_type[d][YL] = buffer;
    bc_type[d][YR] = buffer;
    bc_type[d][ZL] = buffer;
    bc_type[d][ZR] = buffer;
}

void data::assign_bc(int d, BC::INFO xr, BC::INFO xl, BC::INFO yr, BC::INFO yl) {
    BC::INFO buffer;
    bc_type[d][XL] = xl;
    bc_type[d][XR] = xr;
    bc_type[d][YL] = yl;
    bc_type[d][YR] = yr;
    bc_type[d][ZL] = buffer;
    bc_type[d][ZR] = buffer;
}

void data::assign_bc(int d, BC::INFO xr, BC::INFO xl, BC::INFO yr, BC::INFO yl, BC::INFO zr, BC::INFO zl) {
    bc_type[d][XL] = xl;
    bc_type[d][XR] = xr;
    bc_type[d][YL] = yl;
    bc_type[d][YR] = yr;
    bc_type[d][ZL] = zl;
    bc_type[d][ZR] = zr;
}

void data::apply_bc_x(int d) {

    if (id.x == 0) {
        switch (bc_type[d][XL].type % 10) {

            case BC::CELL_CENTER_NEUMANN:
                XL_NEUMANN(d);
                break;

            case BC::CELL_CENTER_DIRICHLET:
                XL_DIRICHLET(d);
                break;

            case BC::CELL_CENTER_PERIODIC:
                XL_PERIODIC(d);
                break;

            case BC::CELL_CENTER_NO_SLIP:
                if (bc_type[d][XL].type == BC::CELL_FACE_X_NO_SLIP) {
                    for (int k = 0; k < NodeSize.z; k++) {
                        for (int j = 0; j < NodeSize.y; j++) {
                            set(0.0, -1, j, k, d);
                            set(-get(0, j, k, d), -2, j, k, d);
                            set(-get(1, j, k, d), -3, j, k, d);
                        }
                    }
                } else {
                    XL_NO_SLIP(d);
                }
                break;

            case BC::CELL_CENTER_FREE_SLIP:
                if (bc_type[d][XL].type == BC::CELL_FACE_X_FREE_SLIP) {
                    for (int k = 0; k < NodeSize.z; k++) {
                        for (int j = 0; j < NodeSize.y; j++) {
                            set(0.0, -1, j, k, d);
                            set(get(0, j, k, d), -2, j, k, d);
                            set(get(1, j, k, d), -3, j, k, d);
                        }
                    }
                } else {
                    XL_FREE_SLIP(d);
                }

                break;
        }
    } else {
        XL_PERIODIC(d);
    }

    if (id.x == GridSize.x - 1) {
        switch (bc_type[d][XR].type % 10) {

            case BC::CELL_CENTER_NEUMANN:
                XR_NEUMANN(d);
                break;

            case BC::CELL_CENTER_DIRICHLET:
                XR_DIRICHLET(d);
                break;

            case BC::CELL_CENTER_PERIODIC:
                XR_PERIODIC(d);
                break;

            case BC::CELL_CENTER_NO_SLIP:
                if (bc_type[d][XR].type == BC::CELL_FACE_X_NO_SLIP) {
                    for (int k = 0; k < NodeSize.z; k++) {
                        for (int j = 0; j < NodeSize.y; j++) {
                            set(0.0, NodeSize.x - 1, j, k, d);
                            set(-get(NodeSize.x - 2, j, k, d), NodeSize.x, j, k, d);
                            set(-get(NodeSize.x - 3, j, k, d), NodeSize.x + 1, j, k, d);
                            set(-get(NodeSize.x - 4, j, k, d), NodeSize.x + 2, j, k, d);
                        }
                    }
                } else {
                    XR_NO_SLIP(d);
                }
                break;

            case BC::CELL_CENTER_FREE_SLIP:
                if (bc_type[d][XR].type == BC::CELL_FACE_X_FREE_SLIP) {
                    for (int k = 0; k < NodeSize.z; k++) {
                        for (int j = 0; j < NodeSize.y; j++) {
                            set(0.0, NodeSize.x - 1, j, k, d);
                            set(get(NodeSize.x - 2, j, k, d), NodeSize.x, j, k, d);
                            set(get(NodeSize.x - 3, j, k, d), NodeSize.x + 1, j, k, d);
                            set(get(NodeSize.x - 4, j, k, d), NodeSize.x + 2, j, k, d);
                        }
                    }
                } else {
                    XR_FREE_SLIP(d);
                }

                break;
        }
    } else {
        XR_PERIODIC(d);
    }

}

void data::XL_NEUMANN(int d) {
    for (int k = 0; k < NodeSize.z; k++) {
        for (int j = 0; j < NodeSize.y; j++) {
            set(get(0, j, k, d) - bc_type[d][XL].value, -1, j, k, d);
            set(get(-1, j, k, d) - bc_type[d][XL].value, -2, j, k, d);
            set(get(-2, j, k, d) - bc_type[d][XL].value, -3, j, k, d);
        }
    }
}

void data::XL_DIRICHLET(int d) {
    for (int k = 0; k < NodeSize.z; k++) {
        for (int j = 0; j < NodeSize.y; j++) {
            set(bc_type[d][XL].value, -3, j, k, d);
            set(bc_type[d][XL].value, -2, j, k, d);
            set(bc_type[d][XL].value, -1, j, k, d);
        }
    }
}

void data::XL_PERIODIC(int d) {
    for (int k = 0; k < NodeSize.z; k++) {
        for (int j = 0; j < NodeSize.y; j++) {
            set(neighbor[XL]->get(NodeSize.x - 1, j, k, d), -1, j, k, d);
            set(neighbor[XL]->get(NodeSize.x - 2, j, k, d), -2, j, k, d);
            set(neighbor[XL]->get(NodeSize.x - 3, j, k, d), -3, j, k, d);
        }
    }
}

void data::XL_NO_SLIP(int d) {
    for (int k = 0; k < NodeSize.z; k++) {
        for (int j = 0; j < NodeSize.y; j++) {
            set(-get(0, j, k, d), -1, j, k, d);
            set(-get(1, j, k, d), -2, j, k, d);
            set(-get(2, j, k, d), -3, j, k, d);
        }
    }
}

void data::XL_FREE_SLIP(int d) {
    for (int k = 0; k < NodeSize.z; k++) {
        for (int j = 0; j < NodeSize.y; j++) {
            set(get(0, j, k, d), -1, j, k, d);
            set(get(1, j, k, d), -2, j, k, d);
            set(get(2, j, k, d), -3, j, k, d);
        }
    }
}

void data::XR_NEUMANN(int d) {
    for (int k = 0; k < NodeSize.z; k++) {
        for (int j = 0; j < NodeSize.y; j++) {
            set(get(NodeSize.x - 1, j, k, d) + bc_type[d][XR].value, NodeSize.x, j, k, d);
            set(get(NodeSize.x, j, k, d) + bc_type[d][XR].value, NodeSize.x + 1, j, k, d);
            set(get(NodeSize.x + 1, j, k, d) + bc_type[d][XR].value, NodeSize.x + 2, j, k, d);
        }
    }
}

void data::XR_DIRICHLET(int d) {
    for (int k = 0; k < NodeSize.z; k++) {
        for (int j = 0; j < NodeSize.y; j++) {
            set(bc_type[d][XR].value, NodeSize.x, j, k, d);
            set(bc_type[d][XR].value, NodeSize.x + 1, j, k, d);
            set(bc_type[d][XR].value, NodeSize.x + 2, j, k, d);
        }
    }
}

void data::XR_PERIODIC(int d) {
    for (int k = 0; k < NodeSize.z; k++) {
        for (int j = 0; j < NodeSize.y; j++) {
            set(neighbor[XR]->get(0, j, k, d), NodeSize.x, j, k, d);
            set(neighbor[XR]->get(1, j, k, d), NodeSize.x + 1, j, k, d);
            set(neighbor[XR]->get(2, j, k, d), NodeSize.x + 2, j, k, d);
        }
    }
}

void data::XR_NO_SLIP(int d) {
    for (int k = 0; k < NodeSize.z; k++) {
        for (int j = 0; j < NodeSize.y; j++) {
            set(-get(NodeSize.x - 1, j, k, d), NodeSize.x, j, k, d);
            set(-get(NodeSize.x - 2, j, k, d), NodeSize.x + 1, j, k, d);
            set(-get(NodeSize.x - 3, j, k, d), NodeSize.x + 2, j, k, d);
        }
    }
}

void data::XR_FREE_SLIP(int d) {
    for (int k = 0; k < NodeSize.z; k++) {
        for (int j = 0; j < NodeSize.y; j++) {
            set(get(NodeSize.x - 1, j, k, d), NodeSize.x, j, k, d);
            set(get(NodeSize.x - 2, j, k, d), NodeSize.x + 1, j, k, d);
            set(get(NodeSize.x - 3, j, k, d), NodeSize.x + 2, j, k, d);
        }
    }
}

void data::ToDevice() {

    cudaSetDevice(gid);
    cudaCheckErrors("cudaSetDevice fail");
    cudaMemcpy(device, host, dim * size, cudaMemcpyHostToDevice);
    cudaCheckErrors("cudaMemcpy fail");

}


