//
// Created by luis on 2023/7/4.
//

#include "runge_kutta_solver.cuh"

runge_kutta_solver::runge_kutta_solver(Idim3 _NodeSize, Idim3 _GridSize, Idim3 id, int gid, bool pinned) {

    NodeSize.x = _NodeSize.x;
    NodeSize.y = _NodeSize.y;
    NodeSize.z = _NodeSize.z;

    spatial_dim = 0;
    if (NodeSize.x > 1) spatial_dim++;
    if (NodeSize.y > 1) spatial_dim++;
    if (NodeSize.z > 1) spatial_dim++;

    solution = new data(_NodeSize, _GridSize, id, gid, 3, pinned);
    src = new data(_NodeSize, _GridSize, id, gid, 3, pinned);

}

void runge_kutta_solver::solve_euler(double dt, variable *phi, int d, variable *vel,
                                     void (*rhs)(variable *, variable *, int)) {

    (*rhs)(phi, vel, d);
    for (int k = 0; k < phi->NodeSize.z; k++) {
        for (int j = 0; j < phi->NodeSize.y; j++) {
            for (int i = 0; i < phi->NodeSize.x; i++) {
                double phi_new = phi->f[d]->get(i, j, k) + dt * phi->flux[d]->get(i, j, k);
                phi->f[d]->set(phi_new, i, j, k);
            }
        }
    }

}

void runge_kutta_solver::solve_tvdrk3(double dt, variable *phi, int d, variable *vel,
                                      void (*rhs)(variable *, variable *, int)) {

    // 1st stage
    (*rhs)(phi, vel, d);
    for (int k = 0; k < phi->NodeSize.z; k++) {
        for (int j = 0; j < phi->NodeSize.y; j++) {
            for (int i = 0; i < phi->NodeSize.x; i++) {
                src->set(phi->flux[d]->get(i, j, k), i, j, k, VAR::STAGE_1);
                double phi_new = phi->f[d]->get(i, j, k) + dt * src->get(i, j, k, VAR::STAGE_1);
                phi->f[d]->set(phi_new, i, j, k);
            }
        }
    }
    phi->f[d]->apply_bc_x(0);

    // 2nd stage
    (*rhs)(phi, vel, d);
    for (int k = 0; k < phi->NodeSize.z; k++) {
        for (int j = 0; j < phi->NodeSize.y; j++) {
            for (int i = 0; i < phi->NodeSize.x; i++) {
                src->set(phi->flux[d]->get(i, j, k), i, j, k, VAR::STAGE_2);
                double phi_new = phi->f[d]->get(i, j, k) + dt * (src->get(i, j, k, VAR::STAGE_2) -
                                                                 3.0 * src->get(i, j, k, VAR::STAGE_1)) / 4.0;
                phi->f[d]->set(phi_new, i, j, k);
            }
        }
    }
    phi->f[d]->apply_bc_x(0);

    // 3rd stage
    (*rhs)(phi, vel, d);
    for (int k = 0; k < phi->NodeSize.z; k++) {
        for (int j = 0; j < phi->NodeSize.y; j++) {
            for (int i = 0; i < phi->NodeSize.x; i++) {
                src->set(phi->flux[d]->get(i, j, k), i, j, k, VAR::STAGE_3);
                double phi_new = phi->f[d]->get(i, j, k) + dt * (8.0 * src->get(i, j, k, VAR::STAGE_3) -
                                                                 src->get(i, j, k, VAR::STAGE_2) -
                                                                 src->get(i, j, k, VAR::STAGE_1)) / 12.0;
                phi->f[d]->set(phi_new, i, j, k);
            }
        }
    }
    phi->f[d]->apply_bc_x(0);

}
