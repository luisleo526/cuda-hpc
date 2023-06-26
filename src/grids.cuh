//
// Created by leo on 6/27/23.
//

#include "data.cuh"

#ifndef CUHPC_GRIDS_CUH
#define CUHPC_GRIDS_CUH

struct DOMAIN_INFO {
    double xl, yl, zl;
    double xr, yr, zr;

    DOMAIN_INFO() : xl(0.0), yl(0.0), zl(0.0), xr(0.0), yr(0.0), zr(0.0) {}

    DOMAIN_INFO(double xl, double xr) : xl(xl), xr(xr), yl(0.0), yr(0.0), zl(0.0), zr(0.0) {}

    DOMAIN_INFO(double xl, double xr, double yl, double yr) : xl(xl), xr(xr), yl(yl), yr(yr), zl(0.0), zr(0.0) {}

    DOMAIN_INFO(double xl, double xr, double yl, double yr, double zl, double zr) : xl(xl), xr(xr), yl(yl), yr(yr),
                                                                                    zl(zl), zr(zr) {}
};

class grids : public data {

public:
    double dx, dy, dz;
    grids(Idim3 _NodeSize, Idim3 _GridSize, Idim3 _id, int _gid, bool pinned, DOMAIN_INFO _domain_info);

};


#endif //CUHPC_GRIDS_CUH
