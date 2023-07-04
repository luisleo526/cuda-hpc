//
// Created by luis on 2023/7/4.
//

#include "variable.cuh"
#include "flux.cuh"

#ifndef CUHPC_RUNGE_KUTTA_SOLVER_CUH
#define CUHPC_RUNGE_KUTTA_SOLVER_CUH

class runge_kutta_solver {
public:
    int spatial_dim;
    Idim3 NodeSize;
    data *solution, *src;

    runge_kutta_solver(Idim3 _NodeSize, Idim3 _GridSize, Idim3 id, int gid, bool pinned);

    void solve_euler(double dt, variable *phi, int d, variable *vel, void(*rhs)(variable *, variable *, int));
    void solve_tvdrk3(double dt, variable *phi, int d, variable *vel, void(*rhs)(variable *, variable *, int));
};


#endif //CUHPC_RUNGE_KUTTA_SOLVER_CUH
