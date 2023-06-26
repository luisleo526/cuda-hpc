//
// Created by leo on 6/27/23.
//

#include "variable.cuh"

variable::variable(Idim3 _NodeSize, Idim3 _GridSize, Idim3 _id, int _gid, int _dim, bool pinned,
                   int _derivative_level, grids *_grid) {

    grid = _grid;
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

    if (spatial_dim == 1) {
        X = 0;
        Y = -1;
        Z = -1;
        XX = 0;
        YY = -1;
        ZZ = -1;
        XY = -1;
        XZ = -1;
        YZ = -1;
        CURV = -1;
    } else if (spatial_dim == 2) {
        X = 0;
        Y = 1;
        XX = 0;
        YY = 1;
        XY = 2;
        CURV = 3;
        Z = -1;
        ZZ = -1;
        XZ = -1;
        YZ = -1;
    } else {
        X = 0;
        Y = 1;
        Z = 2;
        XX = 0;
        YY = 1;
        ZZ = 2;
        XY = 3;
        XZ = 4;
        YZ = 5;
        CURV = 6;
    }

    derivative_level = _derivative_level;

    f = (data **) malloc(dim * sizeof(data *));
    for (int i = 0; i < dim; i++) f[i] = new data(NodeSize, GridSize, id, gid, 1, pinned);

    if (derivative_level > 0) {
        df = (data **) malloc(dim * sizeof(data *));
        for (int i = 0; i < dim; i++) {
            df[i] = new data(NodeSize, GridSize, id, gid, spatial_dim, pinned);
        }
    }

    if (derivative_level > 1) {
        ddf = (data **) malloc(dim * sizeof(data *));
        for (int i = 0; i < dim; i++) {
            ddf[i] = new data(NodeSize, GridSize, id, gid, 1 + 3 * (spatial_dim - 1), pinned);
        }
    }

}

int variable::map_index(int i, int j, int k, int d) {

    int I = i + 3;
    int J = spatial_dim > 1 ? j + 3 : j;
    int K = spatial_dim > 2 ? k + 3 : k;

    return I + J * (NodeSize.x + 6) + K * (NodeSize.x + 6) * (NodeSize.y + 6) + d * size / sizeof(double);
}

void variable::LinkNeighbor() {

    for (int d = 0; d < dim; d++) {
        f[d]->LinkNeighbor();
    }

    if (derivative_level > 0) {
        for (int d = 0; d < dim; d++) {
            df[d]->LinkNeighbor();
        }
    }

    if (derivative_level > 1) {
        for (int d = 0; d < dim; d++) {
            ddf[d]->LinkNeighbor();
        }
    }
}

void variable::LinkNeighbor(variable *_r, variable *_l) {

    for (int d = 0; d < dim; d++) {
        f[d]->LinkNeighbor(_r->f[d], _l->f[d]);
    }

    if (derivative_level > 0) {
        for (int d = 0; d < dim; d++) {
            df[d]->LinkNeighbor(_r->df[d], _l->df[d]);
        }
    }

    if (derivative_level > 1) {
        for (int d = 0; d < dim; d++) {
            ddf[d]->LinkNeighbor(_r->ddf[d], _l->ddf[d]);
        }
    }
}

void variable::LinkNeighbor(variable *_r, variable *_l, variable *_u, variable *_d) {

    for (int d = 0; d < dim; d++) {
        f[d]->LinkNeighbor(_r->f[d], _l->f[d], _u->f[d], _d->f[d]);
    }

    if (derivative_level > 0) {
        for (int d = 0; d < dim; d++) {
            df[d]->LinkNeighbor(_r->df[d], _l->df[d], _u->df[d], _d->df[d]);
        }
    }

    if (derivative_level > 1) {
        for (int d = 0; d < dim; d++) {
            ddf[d]->LinkNeighbor(_r->ddf[d], _l->ddf[d], _u->ddf[d], _d->ddf[d]);
        }
    }

}

void variable::LinkNeighbor(variable *_r, variable *_l, variable *_u, variable *_d, variable *_f, variable *_b) {

    for (int d = 0; d < dim; d++) {
        f[d]->LinkNeighbor(_r->f[d], _l->f[d], _u->f[d], _d->f[d], _f->f[d], _b->f[d]);
    }

    if (derivative_level > 0) {
        for (int d = 0; d < dim; d++) {
            df[d]->LinkNeighbor(_r->df[d], _l->df[d], _u->df[d], _d->df[d], _f->df[d], _b->df[d]);
        }
    }

    if (derivative_level > 1) {
        for (int d = 0; d < dim; d++) {
            ddf[d]->LinkNeighbor(_r->ddf[d], _l->ddf[d], _u->ddf[d], _d->ddf[d], _f->ddf[d], _b->ddf[d]);
        }
    }

}

void variable::get_derivative_SEC(int level) {

    assert(level <= derivative_level && level > 0);

    if (level == 1) {
        for (int d = 0; d < dim; d++) {
            for (int k = 0; k < NodeSize.z; k++) {
                for (int j = 0; j < NodeSize.y; j++) {
                    for (int i = 0; i < NodeSize.x; i++) {
                        double val = 0.5 * (f[d]->get(i + 1, j, k, 0) - f[d]->get(i - 1, j, k, 0)) / grid->dx;
                        df[d]->set(val, i, j, k, X);
                    }
                }
            }
        }
        if (spatial_dim > 1) {
            for (int d = 0; d < dim; d++) {
                for (int k = 0; k < NodeSize.z; k++) {
                    for (int j = 0; j < NodeSize.y; j++) {
                        for (int i = 0; i < NodeSize.x; i++) {
                            double val = 0.5 * (f[d]->get(i, j + 1, k, 0) - f[d]->get(i, j - 1, k, 0)) / grid->dy;
                            df[d]->set(val, i, j, k, Y);
                        }
                    }
                }
            }
        }
        if (spatial_dim > 2) {
            for (int d = 0; d < dim; d++) {
                for (int k = 0; k < NodeSize.z; k++) {
                    for (int j = 0; j < NodeSize.y; j++) {
                        for (int i = 0; i < NodeSize.x; i++) {
                            double val = 0.5 * (f[d]->get(i, j, k + 1, 0) - f[d]->get(i, j, k - 1, 0)) / grid->dz;
                            df[d]->set(val, i, j, k, Z);
                        }
                    }
                }
            }
        }
    } else if (level == 2) {
        for (int d = 0; d < dim; d++) {
            for (int k = 0; k < NodeSize.z; k++) {
                for (int j = 0; j < NodeSize.y; j++) {
                    for (int i = 0; i < NodeSize.x; i++) {
                        double val = (f[d]->get(i + 1, j, k, 0) - 2.0 * f[d]->get(i, j, k, 0) +
                                      f[d]->get(i - 1, j, k, 0)) / grid->dx / grid->dx;
                        ddf[d]->set(val, i, j, k, XX);
                    }
                }
            }
        }
        if (spatial_dim > 1) {
            for (int d = 0; d < dim; d++) {
                for (int k = 0; k < NodeSize.z; k++) {
                    for (int j = 0; j < NodeSize.y; j++) {
                        for (int i = 0; i < NodeSize.x; i++) {
                            double val = (f[d]->get(i, j + 1, k, 0) - 2.0 * f[d]->get(i, j, k, 0) +
                                          f[d]->get(i, j - 1, k, 0)) / grid->dy / grid->dy;
                            ddf[d]->set(val, i, j, k, YY);
                        }
                    }
                }
            }
        }
        if (spatial_dim > 2) {
            for (int d = 0; d < dim; d++) {
                for (int k = 0; k < NodeSize.z; k++) {
                    for (int j = 0; j < NodeSize.y; j++) {
                        for (int i = 0; i < NodeSize.x; i++) {
                            double val = (f[d]->get(i, j, k + 1, 0) - 2.0 * f[d]->get(i, j, k, 0) +
                                          f[d]->get(i, j, k - 1, 0)) / grid->dz / grid->dz;
                            ddf[d]->set(val, i, j, k, ZZ);
                        }
                    }
                }
            }
        }
    }
}
