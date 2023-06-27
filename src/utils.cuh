//
// Created by leo on 6/27/23.
//

#include "constants.cuh"

#ifndef CUHPC_UTILS_CUH
#define CUHPC_UTILS_CUH

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

namespace BC {
    struct INFO {
        int type = CELL_CENTER_NEUMANN;
        double value = 0.0;

        INFO() : type(CELL_CENTER_NEUMANN), value(0.0) {};

        INFO(int t, double v) : type(t), value(v) {};

        INFO(int t) : type(t), value(0.0) {};
    };
}

struct DOMAIN_INFO {
    double xl, yl, zl;
    double xr, yr, zr;

    DOMAIN_INFO() : xl(0.0), yl(0.0), zl(0.0), xr(0.0), yr(0.0), zr(0.0) {}

    DOMAIN_INFO(double xl, double xr) : xl(xl), xr(xr), yl(0.0), yr(0.0), zl(0.0), zr(0.0) {}

    DOMAIN_INFO(double xl, double xr, double yl, double yr) : xl(xl), xr(xr), yl(yl), yr(yr), zl(0.0), zr(0.0) {}

    DOMAIN_INFO(double xl, double xr, double yl, double yr, double zl, double zr) : xl(xl), xr(xr), yl(yl), yr(yr),
                                                                                    zl(zl), zr(zr) {}
};

void twin_dec(double **a, double **b, double **aa, double **bb, int n);

void twin_bks(double **a, double **b, double **aa, double **bb, double *s, double *ss, int n);

#endif //CUHPC_UTILS_CUH
