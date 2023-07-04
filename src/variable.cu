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

    derivative_level = _derivative_level;

    /*
     *   f = {f1, f2, f3 ..., fn }
     *   df = { {df1/dx, df1/dy, df1/dz}, {df2/dx, df2/dy, df2/dz}, ..., {dfn/dx, dfn/dy, dfn/dz} }
     */

    f = (data **) malloc(dim * sizeof(data *));
    flux = (data **) malloc(dim * sizeof(data *));
    for (int i = 0; i < dim; i++) {
        f[i] = new data(NodeSize, GridSize, id, gid, 1, pinned);
        flux[i] = new data(NodeSize, GridSize, id, gid, 1, pinned);
    }

    if (derivative_level > 0) {
        df = (data **) malloc(dim * sizeof(data *));
        for (int i = 0; i < dim; i++) {
            df[i] = new data(NodeSize, GridSize, id, gid, 3, pinned);
        }
    }

    if (derivative_level > 1) {
        ddf = (data **) malloc(dim * sizeof(data *));
        for (int i = 0; i < dim; i++) {
            ddf[i] = new data(NodeSize, GridSize, id, gid, 7, pinned);
        }
    }

    ccd_solver = new ccd(NodeSize, gid, grid);

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

void variable::get_derivative_SEC(int level, int d, int direction) {

    assert(level <= derivative_level && level > 0);
    assert(d < dim);
    assert(direction < spatial_dim);

    if (level == 1) {

        switch (direction) {
            case DIM::X:
                for (int k = 0; k < NodeSize.z; k++) {
                    for (int j = 0; j < NodeSize.y; j++) {
                        for (int i = 0; i < NodeSize.x; i++) {
                            double val = 0.5 * (f[d]->get(i + 1, j, k) - f[d]->get(i - 1, j, k)) / grid->dx;
                            df[d]->set(val, i, j, k, DIM::X);
                        }
                    }
                }
                break;

            case DIM::Y:
                for (int k = 0; k < NodeSize.z; k++) {
                    for (int j = 0; j < NodeSize.y; j++) {
                        for (int i = 0; i < NodeSize.x; i++) {
                            double val = 0.5 * (f[d]->get(i, j + 1, k) - f[d]->get(i, j - 1, k)) / grid->dy;
                            df[d]->set(val, i, j, k, DIM::Y);
                        }
                    }
                }
                break;

            case DIM::Z:
                for (int k = 0; k < NodeSize.z; k++) {
                    for (int j = 0; j < NodeSize.y; j++) {
                        for (int i = 0; i < NodeSize.x; i++) {
                            double val = 0.5 * (f[d]->get(i, j, k + 1) - f[d]->get(i, j, k - 1)) / grid->dz;
                            df[d]->set(val, i, j, k, DIM::Z);
                        }
                    }
                }
                break;
        }

    } else if (level == 2) {

        switch (direction) {

            case DIM::X:
                for (int k = 0; k < NodeSize.z; k++) {
                    for (int j = 0; j < NodeSize.y; j++) {
                        for (int i = 0; i < NodeSize.x; i++) {
                            double val = (f[d]->get(i + 1, j, k) - 2.0 * f[d]->get(i, j, k) +
                                          f[d]->get(i - 1, j, k)) / grid->dx / grid->dx;
                            ddf[d]->set(val, i, j, k, DIM::XX);
                        }
                    }
                }
                break;

            case DIM::Y:
                for (int k = 0; k < NodeSize.z; k++) {
                    for (int j = 0; j < NodeSize.y; j++) {
                        for (int i = 0; i < NodeSize.x; i++) {
                            double val = (f[d]->get(i, j + 1, k) - 2.0 * f[d]->get(i, j, k) +
                                          f[d]->get(i, j - 1, k)) / grid->dy / grid->dy;
                            ddf[d]->set(val, i, j, k, DIM::YY);
                        }
                    }
                }
                break;

            case DIM::Z:
                for (int k = 0; k < NodeSize.z; k++) {
                    for (int j = 0; j < NodeSize.y; j++) {
                        for (int i = 0; i < NodeSize.x; i++) {
                            double val = (f[d]->get(i, j, k + 1) - 2.0 * f[d]->get(i, j, k) +
                                          f[d]->get(i, j, k - 1)) / grid->dz / grid->dz;
                            ddf[d]->set(val, i, j, k, DIM::ZZ);
                        }
                    }
                }
                break;
        }

    }

}

void variable::get_derivative_CCD(int d, int direction) {

    assert(d < dim);
    assert(direction < spatial_dim);
    assert(derivative_level > 0);

    switch (direction) {
        case DIM::X:
            for (int k = 0; k < NodeSize.z; k++) {
                for (int j = 0; j < NodeSize.y; j++) {

                    assign_CCD_source(d, DIM::X, CCD::CENTER, INT_MIN, j, k);

                    twin_bks(ccd_solver->A[CCD::CENTER][DIM::X], ccd_solver->B[CCD::CENTER][DIM::X],
                             ccd_solver->AA[CCD::CENTER][DIM::X], ccd_solver->BB[CCD::CENTER][DIM::X],
                             ccd_solver->S[CCD::CENTER][DIM::X], ccd_solver->SS[CCD::CENTER][DIM::X],
                             ccd_solver->Nodesize[DIM::X]);

                    for (int i = 0; i < NodeSize.x + 6; i++) {
                        df[d]->set(ccd_solver->S[CCD::CENTER][DIM::X][i], i - 3, j, k, DIM::X);
                        if (derivative_level > 1) {
                            ddf[d]->set(ccd_solver->SS[CCD::CENTER][DIM::X][i], i - 3, j, k, DIM::XX);
                        }
                    }
                }
            }
            break;

        case DIM::Y:
            for (int k = 0; k < NodeSize.z; k++) {
                for (int i = 0; i < NodeSize.x; i++) {

                    assign_CCD_source(d, DIM::Y, CCD::CENTER, i, INT_MIN, k);

                    twin_bks(ccd_solver->A[CCD::CENTER][DIM::Y], ccd_solver->B[CCD::CENTER][DIM::Y],
                             ccd_solver->AA[CCD::CENTER][DIM::Y], ccd_solver->BB[CCD::CENTER][DIM::Y],
                             ccd_solver->S[CCD::CENTER][DIM::Y], ccd_solver->SS[CCD::CENTER][DIM::Y],
                             ccd_solver->Nodesize[DIM::Y]);

                    for (int j = 0; j < NodeSize.y + 6; j++) {
                        df[d]->set(ccd_solver->S[CCD::CENTER][DIM::Y][j], i, j - 3, k, DIM::Y);
                        if (derivative_level > 1) {
                            ddf[d]->set(ccd_solver->SS[CCD::CENTER][DIM::Y][j], i, j - 3, k, DIM::YY);
                        }
                    }
                }
            }
            break;

        case DIM::Z:
            for (int j = 0; j < NodeSize.y; j++) {
                for (int i = 0; i < NodeSize.x; i++) {

                    assign_CCD_source(d, DIM::Z, CCD::CENTER, i, j, INT_MIN);

                    twin_bks(ccd_solver->A[CCD::CENTER][DIM::Z], ccd_solver->B[CCD::CENTER][DIM::Z],
                             ccd_solver->AA[CCD::CENTER][DIM::Z], ccd_solver->BB[CCD::CENTER][DIM::Z],
                             ccd_solver->S[CCD::CENTER][DIM::Z], ccd_solver->SS[CCD::CENTER][DIM::Z],
                             ccd_solver->Nodesize[DIM::Z]);

                    for (int k = 0; k < NodeSize.z + 6; k++) {
                        df[d]->set(ccd_solver->S[CCD::CENTER][DIM::Z][k], i, j, k - 3, DIM::Z);
                        if (derivative_level > 1) {
                            ddf[d]->set(ccd_solver->SS[CCD::CENTER][DIM::Z][k], i, j, k - 3, DIM::ZZ);
                        }
                    }

                }
            }
    }

}

void variable::assign_CCD_source(int d, int direction, int type, int I, int J, int K) {

    double bda1 = -3.5;
    double bdb1 = 4.0;
    double bdc1 = -0.5;
    double bda2 = 9.0;
    double bdb2 = -12.0;
    double bdc2 = 3.0;

    double s, ss;

    switch (direction) {
        case DIM::X:
            for (int i = -2; i < NodeSize.x + 2; i++) {

                s = 0.0;
                ss = 0.0;
                for (int dig = 0; dig < 3; dig++) {
                    s += CCD::PARAM1[type][dig] * f[d]->get(i + dig - 1, J, K, 0);
                    ss += CCD::PARAM2[type][dig] * f[d]->get(i + dig - 1, J, K, 0);
                }
                s = s / grid->dx;
                ss = ss / grid->dx / grid->dx;

                ccd_solver->S[type][DIM::X][i + 3] = s;
                ccd_solver->SS[type][DIM::X][i + 3] = ss;
            }


            s = (bda1 * f[d]->get(-3, J, K, 0) + bdb1 * f[d]->get(-2, J, K, 0) +
                 bdc1 * f[d]->get(-1, J, K, 0)) / grid->dx;
            ss = (bda2 * f[d]->get(-3, J, K, 0) + bdb2 * f[d]->get(-2, J, K, 0) +
                  bdc2 * f[d]->get(-1, J, K, 0)) / grid->dx / grid->dx;

            ccd_solver->S[type][DIM::X][0] = s;
            ccd_solver->SS[type][DIM::X][0] = ss;

            s = -(bda1 * f[d]->get(NodeSize.x + 2, J, K, 0) + bdb1 * f[d]->get(NodeSize.x + 1, J, K, 0) +
                  bdc1 * f[d]->get(NodeSize.x, J, K, 0)) / grid->dx;
            ss = (bda2 * f[d]->get(NodeSize.x + 2, J, K, 0) + bdb2 * f[d]->get(NodeSize.x + 1, J, K, 0) +
                  bdc2 * f[d]->get(NodeSize.x, J, K, 0)) / grid->dx / grid->dx;

            ccd_solver->S[type][DIM::X][NodeSize.x + 5] = s;
            ccd_solver->SS[type][DIM::X][NodeSize.x + 5] = ss;

            break;

        case DIM::Y:
            for (int j = 0; j < NodeSize.y; j++) {

                s = 0.0;
                ss = 0.0;
                for (int dig = 0; dig < 3; dig++) {
                    s += CCD::PARAM1[type][dig] * f[d]->get(I, j + dig - 1, K, 0);
                    ss += CCD::PARAM2[type][dig] * f[d]->get(I, j + dig - 1, K, 0);
                }
                s = s / grid->dx;
                ss = ss / grid->dx / grid->dx;

                ccd_solver->S[type][DIM::Y][j + 3] = s;
                ccd_solver->SS[type][DIM::Y][j + 3] = ss;
            }


            s = (bda1 * f[d]->get(I, -3, K, 0) + bdb1 * f[d]->get(I, -2, K, 0) +
                 bdc1 * f[d]->get(I, -1, K, 0)) / grid->dy;
            ss = (bda2 * f[d]->get(I, -3, K, 0) + bdb2 * f[d]->get(I, -2, K, 0) +
                  bdc2 * f[d]->get(I, -1, K, 0)) / grid->dy / grid->dy;

            ccd_solver->S[type][DIM::Y][0] = s;
            ccd_solver->SS[type][DIM::Y][0] = ss;

            s = -(bda1 * f[d]->get(I, NodeSize.y + 2, K, 0) + bdb1 * f[d]->get(I, NodeSize.y + 1, K, 0) +
                  bdc1 * f[d]->get(I, NodeSize.y, K, 0)) / grid->dy;
            ss = (bda2 * f[d]->get(I, NodeSize.y + 2, K, 0) + bdb2 * f[d]->get(I, NodeSize.y + 1, K, 0) +
                  bdc2 * f[d]->get(I, NodeSize.y, K, 0)) / grid->dy / grid->dy;

            ccd_solver->S[type][DIM::Y][NodeSize.y + 5] = s;
            ccd_solver->SS[type][DIM::Y][NodeSize.y + 5] = ss;

            break;

        case DIM::Z:
            for (int k = 0; k < NodeSize.z; k++) {

                s = 0.0;
                ss = 0.0;
                for (int dig = 0; dig < 3; dig++) {
                    s += CCD::PARAM1[type][dig] * f[d]->get(I, J, k + dig - 1, 0);
                    ss += CCD::PARAM2[type][dig] * f[d]->get(I, J, k + dig - 1, 0);
                }
                s = s / grid->dx;
                ss = ss / grid->dx / grid->dx;

                ccd_solver->S[type][DIM::Z][k + 3] = s;
                ccd_solver->SS[type][DIM::Z][k + 3] = ss;
            }

            s = (bda1 * f[d]->get(I, J, -3, 0) + bdb1 * f[d]->get(I, J, -2, 0) +
                 bdc1 * f[d]->get(I, J, -1, 0)) / grid->dz;
            ss = (bda2 * f[d]->get(I, J, -3, 0) + bdb2 * f[d]->get(I, J, -2, 0) +
                  bdc2 * f[d]->get(I, J, -1, 0)) / grid->dz / grid->dz;

            ccd_solver->S[type][DIM::Z][0] = s;
            ccd_solver->SS[type][DIM::Z][0] = ss;

            s = -(bda1 * f[d]->get(I, J, NodeSize.z + 2, 0) + bdb1 * f[d]->get(I, J, NodeSize.z + 1, 0) +
                  bdc1 * f[d]->get(I, J, NodeSize.z, 0)) / grid->dz;

            ss = (bda2 * f[d]->get(I, J, NodeSize.z + 2, 0) + bdb2 * f[d]->get(I, J, NodeSize.z + 1, 0) +
                  bdc2 * f[d]->get(I, J, NodeSize.z, 0)) / grid->dz / grid->dz;

            ccd_solver->S[type][DIM::Z][NodeSize.z + 5] = s;
            ccd_solver->SS[type][DIM::Z][NodeSize.z + 5] = ss;

            break;

    }
}

void variable::get_derivative_UCCD(int d, int direction, variable *vel) {

    assert(d < dim);
    assert(direction < spatial_dim && direction < vel->spatial_dim);
    assert(derivative_level > 0);

    switch (direction) {
        case DIM::X:

            for (int k = 0; k < NodeSize.z; k++) {
                for (int j = 0; j < NodeSize.y; j++) {

                    assign_CCD_source(d, DIM::X, CCD::UPWIND, INT_MIN, j, k);
                    assign_CCD_source(d, DIM::X, CCD::DOWNWIND, INT_MIN, j, k);

                    twin_bks(ccd_solver->A[CCD::UPWIND][DIM::X], ccd_solver->B[CCD::UPWIND][DIM::X],
                             ccd_solver->AA[CCD::UPWIND][DIM::X], ccd_solver->BB[CCD::UPWIND][DIM::X],
                             ccd_solver->S[CCD::UPWIND][DIM::X], ccd_solver->SS[CCD::UPWIND][DIM::X],
                             ccd_solver->Nodesize[DIM::X]);

                    twin_bks(ccd_solver->A[CCD::DOWNWIND][DIM::X], ccd_solver->B[CCD::DOWNWIND][DIM::X],
                             ccd_solver->AA[CCD::DOWNWIND][DIM::X], ccd_solver->BB[CCD::DOWNWIND][DIM::X],
                             ccd_solver->S[CCD::DOWNWIND][DIM::X], ccd_solver->SS[CCD::DOWNWIND][DIM::X],
                             ccd_solver->Nodesize[DIM::X]);

                    for (int i = 0; i < NodeSize.x + 6; i++) {

                        if (vel->f[VAR::U]->get(i - 3, j, k) > 0.0) {
                            df[d]->set(ccd_solver->S[CCD::UPWIND][DIM::X][i], i - 3, j, k, DIM::X);
                        } else {
                            df[d]->set(ccd_solver->S[CCD::DOWNWIND][DIM::X][i], i - 3, j, k, DIM::X);
                        }

                        if (derivative_level > 1) {
                            auto val = (ccd_solver->SS[CCD::UPWIND][DIM::X][i] +
                                        ccd_solver->SS[CCD::DOWNWIND][DIM::X][i]) / 2.0;
                            ddf[d]->set(val, i - 3, j, k, DIM::XX);
                        }
                    }
                }
            }
            break;

        case DIM::Y:
            for (int k = 0; k < NodeSize.z; k++) {
                for (int i = 0; i < NodeSize.x; i++) {
                    assign_CCD_source(d, DIM::Y, CCD::UPWIND, i, INT_MIN, k);
                    assign_CCD_source(d, DIM::Y, CCD::DOWNWIND, i, INT_MIN, k);

                    twin_bks(ccd_solver->A[CCD::UPWIND][DIM::Y], ccd_solver->B[CCD::UPWIND][DIM::Y],
                             ccd_solver->AA[CCD::UPWIND][DIM::Y], ccd_solver->BB[CCD::UPWIND][DIM::Y],
                             ccd_solver->S[CCD::UPWIND][DIM::Y], ccd_solver->SS[CCD::UPWIND][DIM::Y],
                             ccd_solver->Nodesize[DIM::Y]);

                    twin_bks(ccd_solver->A[CCD::DOWNWIND][DIM::Y], ccd_solver->B[CCD::DOWNWIND][DIM::Y],
                             ccd_solver->AA[CCD::DOWNWIND][DIM::Y], ccd_solver->BB[CCD::DOWNWIND][DIM::Y],
                             ccd_solver->S[CCD::DOWNWIND][DIM::Y], ccd_solver->SS[CCD::DOWNWIND][DIM::Y],
                             ccd_solver->Nodesize[DIM::Y]);

                    for (int j = 0; j < NodeSize.y + 6; j++) {

                        if (vel->f[VAR::V]->get(i, j - 3, k) > 0.0) {
                            df[d]->set(ccd_solver->S[CCD::UPWIND][DIM::Y][j], i, j - 3, k, DIM::Y);
                        } else {
                            df[d]->set(ccd_solver->S[CCD::DOWNWIND][DIM::Y][j], i, j - 3, k, DIM::Y);
                        }

                        if (derivative_level > 1) {
                            auto val = (ccd_solver->SS[CCD::UPWIND][DIM::Y][j] +
                                        ccd_solver->SS[CCD::DOWNWIND][DIM::Y][j]) / 2.0;
                            ddf[d]->set(val, i, j - 3, k, DIM::YY);
                        }

                    }
                }
            }
            break;

        case DIM::Z:

            for (int j = 0; j < NodeSize.y; j++) {
                for (int i = 0; i < NodeSize.x; i++) {

                    assign_CCD_source(d, DIM::Z, CCD::UPWIND, i, j, INT_MIN);
                    assign_CCD_source(d, DIM::Z, CCD::DOWNWIND, i, j, INT_MIN);

                    twin_bks(ccd_solver->A[CCD::UPWIND][DIM::Z], ccd_solver->B[CCD::UPWIND][DIM::Z],
                             ccd_solver->AA[CCD::UPWIND][DIM::Z], ccd_solver->BB[CCD::UPWIND][DIM::Z],
                             ccd_solver->S[CCD::UPWIND][DIM::Z], ccd_solver->SS[CCD::UPWIND][DIM::Z],
                             ccd_solver->Nodesize[DIM::Z]);

                    twin_bks(ccd_solver->A[CCD::DOWNWIND][DIM::Z], ccd_solver->B[CCD::DOWNWIND][DIM::Z],
                             ccd_solver->AA[CCD::DOWNWIND][DIM::Z], ccd_solver->BB[CCD::DOWNWIND][DIM::Z],
                             ccd_solver->S[CCD::DOWNWIND][DIM::Z], ccd_solver->SS[CCD::DOWNWIND][DIM::Z],
                             ccd_solver->Nodesize[DIM::Z]);

                    for (int k = 0; k < NodeSize.z + 6; k++) {

                        if (vel->f[VAR::W]->get(i, j, k - 3) > 0.0) {
                            df[d]->set(ccd_solver->S[CCD::UPWIND][DIM::Z][k], i, j, k - 3, DIM::Z);
                        } else {
                            df[d]->set(ccd_solver->S[CCD::DOWNWIND][DIM::Z][k], i, j, k - 3, DIM::Z);
                        }

                        if (derivative_level > 1) {
                            auto val = (ccd_solver->SS[CCD::UPWIND][DIM::Z][k] +
                                        ccd_solver->SS[CCD::DOWNWIND][DIM::Z][k]) / 2.0;
                            ddf[d]->set(val, i, j, k - 3, DIM::ZZ);
                        }

                    }

                }
            }
            break;

    }

}

void variable::get_derivative_USEC(int d, int direction, variable *vel) {

    assert(d < dim);
    assert(direction < spatial_dim && direction < vel->spatial_dim);
    assert(derivative_level > 0);

    switch (direction) {
        case DIM::X:
            for (int k = 0; k < NodeSize.z; k++) {
                for (int j = 0; j < NodeSize.y; j++) {
                    for (int i = 0; i < NodeSize.x; i++) {
                        if (vel->f[VAR::U]->get(i, j, k) > 0.0) {
                            double val =
                                    (3.0 * f[d]->get(i, j, k) - 4.0 * f[d]->get(i - 1, j, k) + f[d]->get(i - 2, j, k)) /
                                    2.0 / grid->dx;
                            df[d]->set(val, i, j, k, DIM::X);
                        } else {
                            double val =
                                    (-3.0 * f[d]->get(i, j, k) + 4.0 * f[d]->get(i + 1, j, k) -
                                     f[d]->get(i + 2, j, k)) / 2.0 / grid->dx;
                            df[d]->set(val, i, j, k, DIM::X);
                        }
                    }
                }
            }
            break;
        case DIM::Y:
            for (int k = 0; k < NodeSize.z; k++) {
                for (int j = 0; j < NodeSize.y; j++) {
                    for (int i = 0; i < NodeSize.x; i++) {
                        if (vel->f[VAR::V]->get(i, j, k) > 0.0) {
                            double val =
                                    (3.0 * f[d]->get(i, j, k) - 4.0 * f[d]->get(i, j - 1, k) + f[d]->get(i, j - 2, k)) /
                                    2.0 / grid->dy;
                            df[d]->set(val, i, j, k, DIM::Y);
                        } else {
                            double val =
                                    (-3.0 * f[d]->get(i, j, k) + 4.0 * f[d]->get(i, j + 1, k) -
                                     f[d]->get(i, j + 2, k)) / 2.0 / grid->dy;
                            df[d]->set(val, i, j, k, DIM::Y);
                        }
                    }
                }
            }
            break;
        case DIM::Z:
            for (int k = 0; k < NodeSize.z; k++) {
                for (int j = 0; j < NodeSize.y; j++) {
                    for (int i = 0; i < NodeSize.x; i++) {
                        if (vel->f[VAR::W]->get(i, j, k) > 0.0) {
                            double val =
                                    (3.0 * f[d]->get(i, j, k) - 4.0 * f[d]->get(i, j, k - 1) + f[d]->get(i, j, k - 2)) /
                                    2.0 / grid->dz;
                            df[d]->set(val, i, j, k, DIM::Z);
                        } else {
                            double val =
                                    (-3.0 * f[d]->get(i, j, k) + 4.0 * f[d]->get(i, j, k + 1) -
                                     f[d]->get(i, j, k + 2)) / 2.0 / grid->dz;
                            df[d]->set(val, i, j, k, DIM::Z);
                        }
                    }
                }
            }
            break;
    }

}
