//
// Created by leo on 6/27/23.
//

#include "ccd.cuh"

ccd::ccd(Idim3 _Nodesize, int _gid, grids *_grid) {

    Idim3 Node(_Nodesize.x, _Nodesize.y, _Nodesize.z);
    gid = _gid;

    mesh_size = (double *) malloc(NUM_OF_DIMS * sizeof(double));
    mesh_size[DIM::X] = _grid->dx;
    mesh_size[DIM::Y] = _grid->dy;
    mesh_size[DIM::Z] = _grid->dz;

    if (Node.x > 1) Node.x += 6;
    if (Node.y > 1) Node.y += 6;
    if (Node.z > 1) Node.z += 6;

    Nodesize = (int *) malloc(NUM_OF_DIMS * sizeof(int));
    Nodesize[DIM::X] = Node.x;
    Nodesize[DIM::Y] = Node.y;
    Nodesize[DIM::Z] = Node.z;

    A = (double ****) malloc(NUM_OF_TYPES * sizeof(double ***));
    B = (double ****) malloc(NUM_OF_TYPES * sizeof(double ***));
    AA = (double ****) malloc(NUM_OF_TYPES * sizeof(double ***));
    BB = (double ****) malloc(NUM_OF_TYPES * sizeof(double ***));
    S = (double ***) malloc(NUM_OF_TYPES * sizeof(double **));
    SS = (double ***) malloc(NUM_OF_TYPES * sizeof(double **));

    for (int type = 0; type < NUM_OF_TYPES; type++) {

        S[type] = (double **) malloc(NUM_OF_DIMS * sizeof(double *));
        SS[type] = (double **) malloc(NUM_OF_DIMS * sizeof(double *));
        for (int d = 0; d < NUM_OF_DIMS; d++) {
            S[type][d] = (double *) malloc(Nodesize[d] * sizeof(double));
            SS[type][d] = (double *) malloc(Nodesize[d] * sizeof(double));
        }

        A[type] = (double ***) malloc(NUM_OF_DIMS * sizeof(double **));
        B[type] = (double ***) malloc(NUM_OF_DIMS * sizeof(double **));
        AA[type] = (double ***) malloc(NUM_OF_DIMS * sizeof(double **));
        BB[type] = (double ***) malloc(NUM_OF_DIMS * sizeof(double **));
        for (int d = 0; d < NUM_OF_DIMS; d++) {
            A[type][d] = (double **) malloc(3 * sizeof(double *));
            B[type][d] = (double **) malloc(3 * sizeof(double *));
            AA[type][d] = (double **) malloc(3 * sizeof(double *));
            BB[type][d] = (double **) malloc(3 * sizeof(double *));
            for (int diag = 0; diag < 3; diag++) {
                A[type][d][diag] = (double *) malloc(Nodesize[d] * sizeof(double));
                B[type][d][diag] = (double *) malloc(Nodesize[d] * sizeof(double));
                AA[type][d][diag] = (double *) malloc(Nodesize[d] * sizeof(double));
                BB[type][d][diag] = (double *) malloc(Nodesize[d] * sizeof(double));
            }
        }
    }

    double a1 = 0.875;
    double b1 = 0.1251282341599089;
    double b2 = -0.2487176584009104;
    double b3 = 0.0001282341599089;

    for (int d = 0; d < NUM_OF_DIMS; d++) {
        for (int i = 0; i < Nodesize[d]; i++) {
            A[CCD::CENTER][d][0][i] = 7.0 / 16.0;
            A[CCD::CENTER][d][1][i] = 1.0;
            A[CCD::CENTER][d][2][i] = 7.0 / 16.0;

            B[CCD::CENTER][d][0][i] = mesh_size[d] / 16.0;
            B[CCD::CENTER][d][1][i] = 0.0;
            B[CCD::CENTER][d][2][i] = -mesh_size[d] / 16.0;

            AA[CCD::CENTER][d][0][i] = -9.0 / 8.0 / mesh_size[d];
            AA[CCD::CENTER][d][1][i] = 0.0;
            AA[CCD::CENTER][d][2][i] = 9.0 / 8.0 / mesh_size[d];

            BB[CCD::CENTER][d][0][i] = -1.0 / 8.0;
            BB[CCD::CENTER][d][1][i] = 1.0;
            BB[CCD::CENTER][d][2][i] = -1.0 / 8.0;
            // -----------------------------------------------------
            A[CCD::UPWIND][d][0][i] = a1;
            A[CCD::UPWIND][d][1][i] = 1.0;
            A[CCD::UPWIND][d][2][i] = 0.0;

            B[CCD::UPWIND][d][0][i] = b1 * mesh_size[d];
            B[CCD::UPWIND][d][1][i] = b2 * mesh_size[d];
            B[CCD::UPWIND][d][2][i] = b3 * mesh_size[d];

            AA[CCD::UPWIND][d][0][i] = -9.0 / 8.0 / mesh_size[d];
            AA[CCD::UPWIND][d][1][i] = 0.0;
            AA[CCD::UPWIND][d][2][i] = 9.0 / 8.0 / mesh_size[d];

            BB[CCD::UPWIND][d][0][i] = -1.0 / 8.0;
            BB[CCD::UPWIND][d][1][i] = 1.0;
            BB[CCD::UPWIND][d][2][i] = -1.0 / 8.0;
            // -----------------------------------------------------
            A[CCD::DOWNWIND][d][0][i] = 0.0;
            A[CCD::DOWNWIND][d][1][i] = 1.0;
            A[CCD::DOWNWIND][d][2][i] = a1;

            B[CCD::DOWNWIND][d][0][i] = -b3 * mesh_size[d];
            B[CCD::DOWNWIND][d][1][i] = -b2 * mesh_size[d];
            B[CCD::DOWNWIND][d][2][i] = -b1 * mesh_size[d];

            AA[CCD::DOWNWIND][d][0][i] = -9.0 / 8.0 / mesh_size[d];
            AA[CCD::DOWNWIND][d][1][i] = 0.0;
            AA[CCD::DOWNWIND][d][2][i] = 9.0 / 8.0 / mesh_size[d];

            BB[CCD::DOWNWIND][d][0][i] = -1.0 / 8.0;
            BB[CCD::DOWNWIND][d][1][i] = 1.0;
            BB[CCD::DOWNWIND][d][2][i] = -1.0 / 8.0;
        }
    }

    double alpha_1 = 2.0;
    double beta_1 = -1.0;
    double alpha_2 = 5.0;
    double beta_2 = -6.0;

    for (int type = 0; type < NUM_OF_TYPES; type++) {
        for (int d = 0; d < NUM_OF_DIMS; d++) {
            A[type][d][0][0] = 0.0;
            A[type][d][1][0] = 1.0;
            A[type][d][2][0] = alpha_1;

            B[type][d][0][0] = 0.0;
            B[type][d][1][0] = 0.0;
            B[type][d][2][0] = beta_1 * mesh_size[d];

            AA[type][d][0][0] = 0.0;
            AA[type][d][1][0] = 0.0;
            AA[type][d][2][0] = beta_2 / mesh_size[d];

            BB[type][d][0][0] = 0.0;
            BB[type][d][1][0] = 1.0;
            BB[type][d][2][0] = alpha_2;

            A[type][d][0][Nodesize[d] - 1] = alpha_1;
            A[type][d][1][Nodesize[d] - 1] = 1.0;
            A[type][d][2][Nodesize[d] - 1] = 0.0;

            B[type][d][0][Nodesize[d] - 1] = -beta_1 * mesh_size[d];
            B[type][d][1][Nodesize[d] - 1] = 0.0;
            B[type][d][2][Nodesize[d] - 1] = 0.0;

            AA[type][d][0][Nodesize[d] - 1] = -beta_2 / mesh_size[d];
            AA[type][d][1][Nodesize[d] - 1] = 0.0;
            AA[type][d][2][Nodesize[d] - 1] = 0.0;

            BB[type][d][0][Nodesize[d] - 1] = alpha_2;
            BB[type][d][1][Nodesize[d] - 1] = 1.0;
            BB[type][d][2][Nodesize[d] - 1] = 0.0;
        }
    }

    for (int type = 0; type < NUM_OF_TYPES; type++) {
        for (int d = 0; d < NUM_OF_DIMS; d++) {
            twin_dec(A[type][d], B[type][d], AA[type][d], BB[type][d], Nodesize[d]);
        }
    }

}
