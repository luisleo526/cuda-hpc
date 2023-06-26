//
// Created by leo on 6/27/23.
//

#include "grids.cuh"

grids::grids(Idim3 _NodeSize, Idim3 _GridSize, Idim3 _id, int _gid, bool pinned, DOMAIN_INFO _domain_info) :
        data(_NodeSize, _GridSize, _id, _gid, 3, pinned) {

    dx = (_domain_info.xr - _domain_info.xl) / GridSize.x / NodeSize.x;
    dy = (_domain_info.yr - _domain_info.yl) / GridSize.y / NodeSize.y;
    dz = (_domain_info.zr - _domain_info.zl) / GridSize.z / NodeSize.z;

    double xs = _domain_info.xl + dx * NodeSize.x * id.x;
    double ys = spatial_dim > 1 ? _domain_info.yl + dy * NodeSize.y * id.y : 0.0;
    double zs = spatial_dim > 2 ? _domain_info.zl + dz * NodeSize.z * id.z : 0.0;

    if (spatial_dim == 1) {
        for (int i = 0; i < dim; i++) {
            assign_bc(i, BC::INFO(BC::CELL_CENTER_NEUMANN), BC::INFO(BC::CELL_CENTER_NEUMANN));
        }
    } else if (spatial_dim == 2) {
        for (int i = 0; i < dim; i++) {
            assign_bc(i,
                      BC::INFO(BC::CELL_CENTER_NEUMANN), BC::INFO(BC::CELL_CENTER_NEUMANN),
                      BC::INFO(BC::CELL_CENTER_NEUMANN), BC::INFO(BC::CELL_CENTER_NEUMANN));
        }
    } else {
        for (int i = 0; i < dim; i++) {
            assign_bc(i,
                      BC::INFO(BC::CELL_CENTER_NEUMANN), BC::INFO(BC::CELL_CENTER_NEUMANN),
                      BC::INFO(BC::CELL_CENTER_NEUMANN), BC::INFO(BC::CELL_CENTER_NEUMANN),
                      BC::INFO(BC::CELL_CENTER_NEUMANN), BC::INFO(BC::CELL_CENTER_NEUMANN));
        }
    }

    for (int k = 0; k < NodeSize.z; k++) {
        for (int j = 0; j < NodeSize.y; j++) {
            for (int i = 0; i < NodeSize.x; i++) {
                set(xs + (i + 0.5) * dx, i, j, k, 0);
                set(ys + (j + 0.5) * dy, i, j, k, 1);
                set(zs + (k + 0.5) * dz, i, j, k, 2);
            }
        }
    }
}
