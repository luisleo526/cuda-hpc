//
// Created by leo on 6/27/23.
//

#ifndef CUHPC_CONSTANTS_CUH
#define CUHPC_CONSTANTS_CUH

namespace CUDA {
    const bool PINNED_MEM = true;
    const bool UNPINNED_MEM = false;
}

namespace VAR {
    const int U = 0;
    const int V = 0;
    const int W = 0;
    const int SCALAR = 0;
    const int WO_DERI = 0;
    const int W_FIRST_DERI = 1;
    const int W_SECOND_DERI = 2;
}

namespace CCD {
    const int UPWIND = 0;
    const int DOWNWIND = 1;
    const int CENTER = 2;
    const double CENTER_PARAM[3][3] = {-15.0 / 16.0, 0.0, 15.0 / 16.0};
}

namespace DIM {
    const int X = 0;
    const int Y = 1;
    const int Z = 2;
    const int XX = 0;
    const int YY = 1;
    const int ZZ = 2;
    const int XY = 3;
    const int YZ = 4;
    const int XZ = 5;
    const int CURV = 6;
}

namespace GRID {
    int const XL = 0;
    int const XR = 1;
    int const YL = 2;
    int const YR = 3;
    int const ZL = 4;
    int const ZR = 5;
}

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
}

#endif //CUHPC_CONSTANTS_CUH
