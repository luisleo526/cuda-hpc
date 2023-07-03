//
// Created by luis on 2023/7/4.
//

#include "variable.cuh"
#include "flux.cuh"

#ifndef CUHPC_RUNGE_KUTTA_SOLVER_CUH
#define CUHPC_RUNGE_KUTTA_SOLVER_CUH

class runge_kutta_solver {
public:
    data *solution, *src;
};


#endif //CUHPC_RUNGE_KUTTA_SOLVER_CUH
